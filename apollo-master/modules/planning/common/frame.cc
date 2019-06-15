/******************************************************************************
 * Copyright 2017 The Apollo Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/

/**
 * @file frame.cc
 **/
#include "modules/planning/common/frame.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <list>
#include <string>
#include <utility>

#include "modules/routing/proto/routing.pb.h"

#include "modules/common/adapters/adapter_manager.h"
#include "modules/common/configs/vehicle_config_helper.h"
#include "modules/common/log.h"
#include "modules/common/math/vec2d.h"
#include "modules/common/vehicle_state/vehicle_state_provider.h"
#include "modules/map/hdmap/hdmap_util.h"
#include "modules/planning/common/planning_gflags.h"
#include "modules/planning/reference_line/reference_line_provider.h"

namespace apollo {
namespace planning {

using apollo::common::ErrorCode;
using apollo::common::Status;
using apollo::common::VehicleStateProvider;
using apollo::common::adapter::AdapterManager;
using apollo::common::math::Box2d;
using apollo::common::math::Vec2d;
using apollo::prediction::PredictionObstacles;

FrameHistory::FrameHistory()
    : IndexedQueue<uint32_t, Frame>(FLAGS_max_history_frame_num) {}

Frame::Frame(uint32_t sequence_num,
             const common::TrajectoryPoint &planning_start_point,
             const double start_time, const common::VehicleState &vehicle_state)
    : sequence_num_(sequence_num),
      planning_start_point_(planning_start_point),
      start_time_(start_time),
      vehicle_state_(vehicle_state) {
  if (FLAGS_enable_lag_prediction) {
    lag_predictor_.reset(
        new LagPrediction(FLAGS_lag_prediction_min_appear_num,
                          FLAGS_lag_prediction_max_disappear_num));
  }
}

const common::TrajectoryPoint &Frame::PlanningStartPoint() const {
  return planning_start_point_;
}

const common::VehicleState &Frame::vehicle_state() const {
  return vehicle_state_;
}

bool Frame::Rerouting() {
  auto *adapter_manager = AdapterManager::instance();
  if (adapter_manager->GetRoutingResponse()->Empty()) {
    AERROR << "No previous routing available";
    return false;
  }
  auto request = adapter_manager->GetRoutingResponse()
                     ->GetLatestObserved()
                     .routing_request();
  request.clear_header();
  AdapterManager::FillRoutingRequestHeader("planning", &request);
  auto point = common::util::MakePointENU(
      vehicle_state_.x(), vehicle_state_.y(), vehicle_state_.z());
  double s = 0.0;
  double l = 0.0;
  hdmap::LaneInfoConstPtr lane;
  if (hdmap_->GetNearestLaneWithHeading(point, 5.0, vehicle_state_.heading(),
                                        M_PI / 3.0, &lane, &s, &l) != 0) {
    AERROR << "Failed to find nearest lane from map at position: "
           << point.DebugString() << ", heading:" << vehicle_state_.heading();
    return false;
  }
  request.clear_waypoint();
  auto *start_point = request.add_waypoint();
  start_point->set_id(lane->id().id());
  start_point->set_s(s);
  start_point->mutable_pose()->CopyFrom(point);
  for (const auto &waypoint :
       ReferenceLineProvider::instance()->FutureRouteWaypoints()) {
    request.add_waypoint()->CopyFrom(waypoint);
  }
  if (request.waypoint_size() <= 1) {
    AERROR << "Failed to find future waypoints";
    return false;
  }
  AdapterManager::PublishRoutingRequest(request);
  return true;
}

std::list<ReferenceLineInfo> &Frame::reference_line_info() {
  return reference_line_info_;
}

bool Frame::CreateReferenceLineInfo() {
  std::list<ReferenceLine> reference_lines;
  std::list<hdmap::RouteSegments> segments;
  if (!ReferenceLineProvider::instance()->GetReferenceLines(&reference_lines,
                                                            &segments)) {
    AERROR << "Failed to create reference line";
    return false;
  }
  DCHECK_EQ(reference_lines.size(), segments.size());

  auto forword_limit =
      ReferenceLineProvider::LookForwardDistance(vehicle_state_);

  for (auto &ref_line : reference_lines) {
    if (!ref_line.Shrink(Vec2d(vehicle_state_.x(), vehicle_state_.y()),
                         FLAGS_look_backward_distance, forword_limit)) {
      AERROR << "Fail to shrink reference line.";
      return false;
    }
  }
  for (auto &seg : segments) {
    if (!seg.Shrink(Vec2d(vehicle_state_.x(), vehicle_state_.y()),
                    FLAGS_look_backward_distance, forword_limit)) {
      AERROR << "Fail to shrink routing segments.";
      return false;
    }
  }

  reference_line_info_.clear();
  bool near_destination = false;
  auto ref_line_iter = reference_lines.begin();
  auto segments_iter = segments.begin();
  while (ref_line_iter != reference_lines.end()) {
    if (segments_iter->StopForDestination()) {
      near_destination = true;
    }
    reference_line_info_.emplace_back(vehicle_state_, planning_start_point_,
                                      *ref_line_iter, *segments_iter);
    ++ref_line_iter;
    ++segments_iter;
  }
  if (near_destination && CreateDestinationObstacle() < 0) {
    AERROR << "Failed to create the destination obstacle";
    return false;
  }
  if (FLAGS_enable_change_lane_decider &&
      !change_lane_decider_.Apply(&reference_line_info_)) {
    AERROR << "Failed to apply change lane decider";
    return false;
  }

  if (reference_line_info_.size() == 2) {
    common::math::Vec2d xy_point(vehicle_state_.x(), vehicle_state_.y());
    common::SLPoint first_sl;
    if (!reference_line_info_.front().reference_line().XYToSL(xy_point,
                                                              &first_sl)) {
      return false;
    }
    common::SLPoint second_sl;
    if (!reference_line_info_.back().reference_line().XYToSL(xy_point,
                                                             &second_sl)) {
      return false;
    }
    const double offset = first_sl.l() - second_sl.l();
    reference_line_info_.front().SetOffsetToOtherReferenceLine(offset);
    reference_line_info_.back().SetOffsetToOtherReferenceLine(-offset);
  }

  // delay the time-consumping reference_line_info init() step to planner.
  return true;
}

const Obstacle *Frame::AddStaticVirtualObstacle(const std::string &id,
                                                const Box2d &box) {
  const auto *object = obstacles_.Find(id);
  if (object) {
    AWARN << "obstacle " << id << " already exist.";
    return object;
  }
  auto *ptr =
      obstacles_.Add(id, *Obstacle::CreateStaticVirtualObstacles(id, box));
  if (!ptr) {
    AERROR << "Failed to create virtual obstacle " << id;
  }
  return ptr;
}

int Frame::CreateDestinationObstacle() {
  const auto &routing =
      AdapterManager::GetRoutingResponse()->GetLatestObserved();
  if (routing.routing_request().waypoint_size() < 2) {
    ADEBUG << "routing_request has no end";
    return -1;
  }
  const auto &routing_end = *routing.routing_request().waypoint().rbegin();
  const auto lane = hdmap_->GetLaneById(hdmap::MakeMapId(routing_end.id()));
  if (!lane) {
    AERROR << "Failed to find lane for destination : "
           << routing_end.ShortDebugString();
    return -2;
  }

  double dest_lane_s =
      std::max(0.0, routing_end.s() - FLAGS_virtual_stop_wall_length -
                        FLAGS_stop_distance_destination);
  auto dest_point = lane->GetSmoothPoint(dest_lane_s);
  double left_width = 0.0;
  double right_width = 0.0;
  lane->GetWidth(dest_lane_s, &left_width, &right_width);
  // check if destination point is in planning range
  Box2d destination_box{{dest_point.x(), dest_point.y()},
                        lane->Heading(dest_lane_s),
                        FLAGS_virtual_stop_wall_length,
                        left_width + right_width};
  // add destination's projection to each reference line info
  AddStaticVirtualObstacle(FLAGS_destination_obstacle_id, destination_box);
  return 0;
}

Status Frame::Init() {
  hdmap_ = hdmap::HDMapUtil::BaseMapPtr();
  vehicle_state_ = common::VehicleStateProvider::instance()->vehicle_state();
  const auto &point = common::util::MakePointENU(
      vehicle_state_.x(), vehicle_state_.y(), vehicle_state_.z());
  if (std::isnan(point.x()) || std::isnan(point.y())) {
    AERROR << "init point is not set";
    return Status(ErrorCode::PLANNING_ERROR, "init point is not set");
  }
  ADEBUG << "Enabled align prediction time ? : " << std::boolalpha
         << FLAGS_align_prediction_time;

  // prediction
  if (FLAGS_enable_prediction && AdapterManager::GetPrediction() &&
      !AdapterManager::GetPrediction()->Empty()) {
    if (FLAGS_enable_lag_prediction && lag_predictor_) {
      lag_predictor_->GetLaggedPrediction(&prediction_);
    } else {
      prediction_.CopyFrom(
          AdapterManager::GetPrediction()->GetLatestObserved());
    }
    if (FLAGS_align_prediction_time) {
      AlignPredictionTime(vehicle_state_.timestamp(), &prediction_);
    }
    for (auto &ptr : Obstacle::CreateObstacles(prediction_)) {
      AddObstacle(*ptr);
    }
  }
  const auto *collision_obstacle = FindCollisionObstacle();
  if (collision_obstacle) {
    AERROR << "Found collision with obstacle: " << collision_obstacle->Id();
    return Status(ErrorCode::PLANNING_ERROR,
                  "Collision found with " + collision_obstacle->Id());
  }
  if (!CreateReferenceLineInfo()) {
    AERROR << "Failed to init reference line info";
    return Status(ErrorCode::PLANNING_ERROR,
                  "failed to init reference line info");
  }

  return Status::OK();
}

const Obstacle *Frame::FindCollisionObstacle() const {
  if (obstacles_.Items().empty()) {
    return nullptr;
  }
  const auto &param =
      common::VehicleConfigHelper::instance()->GetConfig().vehicle_param();
  Vec2d position(vehicle_state_.x(), vehicle_state_.y());
  Vec2d vec_to_center(
      (param.front_edge_to_center() - param.back_edge_to_center()) / 2.0,
      (param.left_edge_to_center() - param.right_edge_to_center()) / 2.0);
  Vec2d center(position + vec_to_center.rotate(vehicle_state_.heading()));
  Box2d adc_box(center, vehicle_state_.heading(), param.length(),
                param.width());
  const double adc_half_diagnal = adc_box.diagonal() / 2.0;
  for (const auto &obstacle : obstacles_.Items()) {
    if (obstacle->IsVirtual()) {
      continue;
    }

    double center_dist =
        adc_box.center().DistanceTo(obstacle->PerceptionBoundingBox().center());
    if (center_dist > obstacle->PerceptionBoundingBox().diagonal() / 2.0 +
                          adc_half_diagnal + FLAGS_max_collision_distance) {
      ADEBUG << "Obstacle : " << obstacle->Id() << " is too far to collide";
      continue;
    }
    if (obstacle->PerceptionPolygon().DistanceTo(adc_box) <
        FLAGS_max_collision_distance) {
      AERROR << "Found collision with obstacle " << obstacle->Id();
      return obstacle;
    }
  }
  return nullptr;
}

uint32_t Frame::SequenceNum() const { return sequence_num_; }

std::string Frame::DebugString() const {
  return "Frame: " + std::to_string(sequence_num_);
}

void Frame::RecordInputDebug(planning_internal::Debug *debug) {
  if (!debug) {
    ADEBUG << "Skip record input into debug";
    return;
  }
  auto *planning_data = debug->mutable_planning_data();
  auto *adc_position = planning_data->mutable_adc_position();
  const auto &localization =
      AdapterManager::GetLocalization()->GetLatestObserved();
  adc_position->CopyFrom(localization);

  const auto &chassis = AdapterManager::GetChassis()->GetLatestObserved();
  auto debug_chassis = planning_data->mutable_chassis();
  debug_chassis->CopyFrom(chassis);

  auto debug_routing = planning_data->mutable_routing();
  debug_routing->CopyFrom(
      AdapterManager::GetRoutingResponse()->GetLatestObserved());

  planning_data->mutable_prediction_header()->CopyFrom(prediction_.header());
}

void Frame::AlignPredictionTime(const double planning_start_time,
                                PredictionObstacles *prediction_obstacles) {
  if (!prediction_obstacles || !prediction_obstacles->has_header() ||
      !prediction_obstacles->header().has_timestamp_sec()) {
    return;
  }
  double prediction_header_time =
      prediction_obstacles->header().timestamp_sec();
  for (auto &obstacle : *prediction_obstacles->mutable_prediction_obstacle()) {
    for (auto &trajectory : *obstacle.mutable_trajectory()) {
      for (auto &point : *trajectory.mutable_trajectory_point()) {
        point.set_relative_time(prediction_header_time + point.relative_time() -
                                planning_start_time);
      }
      if (!trajectory.trajectory_point().empty() &&
          trajectory.trajectory_point().begin()->relative_time() < 0) {
        auto it = trajectory.trajectory_point().begin();
        while (it != trajectory.trajectory_point().end() &&
               it->relative_time() < 0) {
          ++it;
        }
        trajectory.mutable_trajectory_point()->erase(
            trajectory.trajectory_point().begin(), it);
      }
    }
  }
}

Obstacle *Frame::Find(const std::string &id) { return obstacles_.Find(id); }

void Frame::AddObstacle(const Obstacle &obstacle) {
  obstacles_.Add(obstacle.Id(), obstacle);
}

const ReferenceLineInfo *Frame::FindDriveReferenceLineInfo() {
  double min_cost = std::numeric_limits<double>::infinity();
  for (const auto &reference_line_info : reference_line_info_) {
    if (reference_line_info.ReachedDestination() ||
        (reference_line_info.IsDrivable() &&
         reference_line_info.Cost() < min_cost)) {
      drive_reference_line_info_ = &reference_line_info;
      min_cost = reference_line_info.Cost();
    }
  }
  return drive_reference_line_info_;
}

const ReferenceLineInfo *Frame::DriveReferenceLineInfo() const {
  return drive_reference_line_info_;
}

const std::vector<const Obstacle *> Frame::obstacles() const {
  return obstacles_.Items();
}

}  // namespace planning
}  // namespace apollo
