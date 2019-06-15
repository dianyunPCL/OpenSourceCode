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
 * @file
 **/

#include "modules/planning/tasks/traffic_decider/signal_light.h"

#include <limits>
#include <vector>

#include "modules/common/proto/pnc_point.pb.h"
#include "modules/common/util/map_util.h"
#include "modules/planning/proto/planning_internal.pb.h"

#include "modules/common/adapters/adapter_manager.h"
#include "modules/common/vehicle_state/vehicle_state_provider.h"
#include "modules/planning/common/frame.h"
#include "modules/planning/common/planning_gflags.h"

namespace apollo {
namespace planning {

using apollo::common::adapter::AdapterManager;
using apollo::perception::TrafficLight;
using apollo::perception::TrafficLightDetection;
using apollo::common::util::WithinBound;

SignalLight::SignalLight(const RuleConfig& config) : TrafficRule(config) {}

bool SignalLight::ApplyRule(Frame* frame,
                            ReferenceLineInfo* const reference_line_info) {
  if (!FindValidSignalLight(reference_line_info)) {
    return true;
  }
  ReadSignals();
  MakeDecisions(frame, reference_line_info);
  return true;
}

void SignalLight::ReadSignals() {
  detected_signals_.clear();
  if (AdapterManager::GetTrafficLightDetection()->Empty()) {
    return;
  }
  if (AdapterManager::GetTrafficLightDetection()->GetDelaySec() >
      FLAGS_signal_expire_time_sec) {
    ADEBUG << "traffic signals msg is expired: "
           << AdapterManager::GetTrafficLightDetection()->GetDelaySec();
    return;
  }
  const TrafficLightDetection& detection =
      AdapterManager::GetTrafficLightDetection()->GetLatestObserved();
  for (int j = 0; j < detection.traffic_light_size(); j++) {
    const TrafficLight& signal = detection.traffic_light(j);
    detected_signals_[signal.id()] = &signal;
  }
}

bool SignalLight::FindValidSignalLight(
    ReferenceLineInfo* const reference_line_info) {
  const std::vector<hdmap::PathOverlap>& signal_lights =
      reference_line_info->reference_line().map_path().signal_overlaps();
  if (signal_lights.size() <= 0) {
    ADEBUG << "No signal lights from reference line.";
    return false;
  }
  signal_lights_from_path_.clear();
  for (const hdmap::PathOverlap& signal_light : signal_lights) {
    if (signal_light.start_s + FLAGS_stop_max_distance_buffer >
        reference_line_info->AdcSlBoundary().end_s()) {
      signal_lights_from_path_.push_back(signal_light);
    }
  }
  return signal_lights_from_path_.size() > 0;
}

void SignalLight::MakeDecisions(Frame* frame,
                                ReferenceLineInfo* const reference_line_info) {
  planning_internal::SignalLightDebug* signal_light_debug =
      reference_line_info->mutable_debug()
          ->mutable_planning_data()
          ->mutable_signal_light();
  signal_light_debug->set_adc_front_s(
      reference_line_info->AdcSlBoundary().end_s());
  signal_light_debug->set_adc_speed(
      common::VehicleStateProvider::instance()->linear_velocity());

  bool has_stop = false;
  for (auto& signal_light : signal_lights_from_path_) {
    const TrafficLight signal = GetSignal(signal_light.object_id);
    double stop_deceleration =
        GetStopDeceleration(reference_line_info, &signal_light);

    planning_internal::SignalLightDebug::SignalDebug* signal_debug =
        signal_light_debug->add_signal();
    signal_debug->set_adc_stop_deacceleration(stop_deceleration);
    signal_debug->set_color(signal.color());
    signal_debug->set_light_id(signal_light.object_id);
    signal_debug->set_light_stop_s(signal_light.start_s);

    if ((signal.color() == TrafficLight::RED &&
         stop_deceleration < FLAGS_stop_max_deceleration) ||
        (signal.color() == TrafficLight::UNKNOWN &&
         stop_deceleration < FLAGS_stop_max_deceleration) ||
        (signal.color() == TrafficLight::YELLOW &&
         stop_deceleration < FLAGS_max_deacceleration_for_yellow_light_stop)) {
      if (FLAGS_right_turn_creep_forward &&
          reference_line_info->IsRightTurnPath()) {
        SetCreepForwardSignalDecision(reference_line_info, &signal_light);
      }
      if (CreateStopObstacle(frame, reference_line_info, &signal_light)) {
        has_stop = true;
        signal_debug->set_is_stop_wall_created(true);
      }
    }
    if (has_stop) {
      reference_line_info->SetJunctionRightOfWay(signal_light.start_s,
                                                 false);  // not protected
    } else {
      reference_line_info->SetJunctionRightOfWay(signal_light.start_s, true);
      // is protected
    }
  }
}

void SignalLight::SetCreepForwardSignalDecision(
    const ReferenceLineInfo* reference_line_info,
    hdmap::PathOverlap* const signal_light) const {
  CHECK_NOTNULL(signal_light);

  constexpr double kMaxCreepSpeed = 1.0;
  if (reference_line_info->AdcPlanningPoint().v() > kMaxCreepSpeed) {
    ADEBUG << "Do not creep forward due to large speed.";
    return;
  }

  constexpr double kCreepBuff = 3.0;
  const auto& path_decision = reference_line_info->path_decision();
  for (const auto& path_obstacle : path_decision.path_obstacles().Items()) {
    const auto& st_boundary = path_obstacle->st_boundary();
    const double stop_s =
        signal_light->start_s - FLAGS_stop_distance_traffic_light;
    if (reference_line_info->AdcSlBoundary().end_s() + st_boundary.min_s() <
        stop_s + kCreepBuff) {
      AERROR << "Do not creep forward because obstacles are close.";
      return;
    }
  }
  signal_light->start_s = reference_line_info->AdcSlBoundary().end_s() +
                          FLAGS_stop_distance_traffic_light + kCreepBuff;
  ADEBUG << "Creep forward s = " << signal_light->start_s;
}

TrafficLight SignalLight::GetSignal(const std::string& signal_id) {
  const auto* result =
      apollo::common::util::FindPtrOrNull(detected_signals_, signal_id);
  if (result == nullptr) {
    TrafficLight traffic_light;
    traffic_light.set_id(signal_id);
    traffic_light.set_color(TrafficLight::UNKNOWN);
    traffic_light.set_confidence(0.0);
    traffic_light.set_tracking_time(0.0);
    return traffic_light;
  }
  return *result;
}

double SignalLight::GetStopDeceleration(
    ReferenceLineInfo* const reference_line_info,
    const hdmap::PathOverlap* signal_light) {
  double adc_speed =
      common::VehicleStateProvider::instance()->linear_velocity();
  if (adc_speed < FLAGS_stop_min_speed) {
    return 0.0;
  }
  double stop_distance = 0;
  double adc_front_s = reference_line_info->AdcSlBoundary().end_s();
  double stop_line_s = signal_light->start_s;

  if (stop_line_s > adc_front_s) {
    stop_distance = stop_line_s - adc_front_s;
  } else {
    stop_distance = stop_line_s + FLAGS_stop_max_distance_buffer - adc_front_s;
  }
  if (stop_distance < 1e-5) {
    return std::numeric_limits<double>::max();
  }
  return (adc_speed * adc_speed) / (2 * stop_distance);
}

bool SignalLight::CreateStopObstacle(
    Frame* frame, ReferenceLineInfo* const reference_line_info,
    const hdmap::PathOverlap* signal_light) {
  const auto& reference_line = reference_line_info->reference_line();
  const double stop_s =
      signal_light->start_s - FLAGS_stop_distance_traffic_light;
  const double box_center_s =
      signal_light->start_s + FLAGS_virtual_stop_wall_length / 2.0;
  if (!WithinBound(0.0, reference_line.Length(), stop_s) ||
      !WithinBound(0.0, reference_line.Length(), box_center_s)) {
    ADEBUG << "signal " << signal_light->object_id
           << " is not on reference line";
    return false;
  }
  double heading = reference_line.GetReferencePoint(stop_s).heading();
  double left_lane_width = 0.0;
  double right_lane_width = 0.0;
  reference_line.GetLaneWidth(signal_light->start_s, &left_lane_width,
                              &right_lane_width);

  auto box_center = reference_line.GetReferencePoint(box_center_s);
  common::math::Box2d stop_box{box_center, heading,
                               FLAGS_virtual_stop_wall_length,
                               left_lane_width + right_lane_width};

  PathObstacle* stop_wall =
      reference_line_info->AddObstacle(frame->AddStaticVirtualObstacle(
          FLAGS_signal_light_virtual_object_id_prefix + signal_light->object_id,
          stop_box));
  auto* path_decision = reference_line_info->path_decision();
  auto stop_point = reference_line.GetReferencePoint(stop_s);
  ObjectDecisionType stop;
  auto stop_decision = stop.mutable_stop();
  stop_decision->set_reason_code(StopReasonCode::STOP_REASON_SIGNAL);
  stop_decision->set_distance_s(-FLAGS_stop_distance_traffic_light);
  stop_decision->set_stop_heading(heading);
  stop_decision->mutable_stop_point()->set_x(stop_point.x());
  stop_decision->mutable_stop_point()->set_y(stop_point.y());
  stop_decision->mutable_stop_point()->set_z(0.0);

  if (!path_decision->MergeWithMainStop(stop.stop(), stop_wall->Id(),
                                        reference_line_info->reference_line(),
                                        reference_line_info->AdcSlBoundary())) {
    ADEBUG << "signal " << signal_light->object_id
           << " is not the cloest stop.";
    return false;
  }
  path_decision->AddLongitudinalDecision(
      RuleConfig::RuleId_Name(config_.rule_id()), stop_wall->Id(), stop);
  return true;
}

}  // namespace planning
}  // namespace apollo
