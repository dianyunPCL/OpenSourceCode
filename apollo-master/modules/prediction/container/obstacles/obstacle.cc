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

#include "modules/prediction/container/obstacles/obstacle.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <unordered_set>
#include <utility>

#include "modules/common/log.h"
#include "modules/common/math/math_utils.h"
#include "modules/common/util/map_util.h"
#include "modules/prediction/common/prediction_gflags.h"
#include "modules/prediction/common/prediction_map.h"
#include "modules/prediction/common/road_graph.h"
#include "modules/prediction/network/rnn_model/rnn_model.h"

namespace apollo {
namespace prediction {

using ::apollo::common::ErrorCode;
using ::apollo::common::Point3D;
using ::apollo::common::math::KalmanFilter;
using ::apollo::common::util::FindOrDie;
using ::apollo::common::util::FindOrNull;
using ::apollo::hdmap::LaneInfo;
using ::apollo::perception::PerceptionObstacle;

namespace {

double Damp(const double x, const double sigma) {
  return 1 / (1 + exp(1 / (std::fabs(x) + sigma)));
}

}  // namespace

Obstacle::Obstacle() {
  double heading_filter_param = FLAGS_heading_filter_param;
  if (FLAGS_heading_filter_param < 0.0 || FLAGS_heading_filter_param > 1.0) {
    heading_filter_param = 0.98;
  }
  heading_filter_ = common::DigitalFilter{{1.0, 1.0 - heading_filter_param},
                                          {heading_filter_param}};
}

PerceptionObstacle::Type Obstacle::type() const { return type_; }

int Obstacle::id() const { return id_; }

double Obstacle::timestamp() const {
  if (feature_history_.size() > 0) {
    return feature_history_.front().timestamp();
  } else {
    return 0.0;
  }
}

const Feature& Obstacle::feature(size_t i) const {
  CHECK(i < feature_history_.size());
  return feature_history_[i];
}

Feature* Obstacle::mutable_feature(size_t i) {
  CHECK(i < feature_history_.size());
  return &feature_history_[i];
}

const Feature& Obstacle::latest_feature() const {
  CHECK_GT(feature_history_.size(), 0);
  return feature_history_.front();
}

Feature* Obstacle::mutable_latest_feature() {
  CHECK_GT(feature_history_.size(), 0);
  return &(feature_history_.front());
}

size_t Obstacle::history_size() const { return feature_history_.size(); }

const KalmanFilter<double, 4, 2, 0>& Obstacle::kf_lane_tracker(
    const std::string& lane_id) {
  return FindOrDie(kf_lane_trackers_, lane_id);
}

const KalmanFilter<double, 6, 2, 0>& Obstacle::kf_motion_tracker() const {
  return kf_motion_tracker_;
}

const KalmanFilter<double, 2, 2, 4>& Obstacle::kf_pedestrian_tracker() const {
  return kf_pedestrian_tracker_;
}

bool Obstacle::IsOnLane() {
  if (feature_history_.size() > 0) {
    if (feature_history_.front().has_lane() &&
        (feature_history_.front().lane().current_lane_feature_size() > 0 ||
         feature_history_.front().lane().nearby_lane_feature_size() > 0)) {
      ADEBUG << "Obstacle [" << id_ << "] is on lane.";
      return true;
    }
  }
  ADEBUG << "Obstacle [" << id_ << "] is not on lane.";
  return false;
}

bool Obstacle::IsNearJunction() {
  if (feature_history_.size() <= 0) {
    return false;
  }
  double pos_x = latest_feature().position().x();
  double pos_y = latest_feature().position().y();
  return PredictionMap::instance()->NearJunction({pos_x, pos_y},
                                                 FLAGS_junction_search_radius);
}

void Obstacle::Insert(const PerceptionObstacle& perception_obstacle,
                      const double timestamp) {
  if (feature_history_.size() > 0 &&
      timestamp <= feature_history_.front().timestamp()) {
    AERROR << "Obstacle [" << id_ << "] received an older frame ["
           << std::setprecision(20) << timestamp
           << "] than the most recent timestamp [ "
           << feature_history_.front().timestamp() << "].";
    return;
  }

  Feature feature;
  if (SetId(perception_obstacle, &feature) == ErrorCode::PREDICTION_ERROR) {
    return;
  }
  if (SetType(perception_obstacle) == ErrorCode::PREDICTION_ERROR) {
    return;
  }

  // Set obstacle observation for KF tracking
  SetStatus(perception_obstacle, timestamp, &feature);

  // Update KF
  if (!kf_motion_tracker_.IsInitialized()) {
    InitKFMotionTracker(feature);
  }
  UpdateKFMotionTracker(feature);
  if (type_ == PerceptionObstacle::PEDESTRIAN) {
    if (!kf_pedestrian_tracker_.IsInitialized()) {
      InitKFPedestrianTracker(feature);
    }
    UpdateKFPedestrianTracker(feature);
  }

  // Update obstacle status based on KF if enabled
  if (FLAGS_enable_kf_tracking) {
    UpdateStatus(&feature);
  }

  // Set obstacle lane features
  SetCurrentLanes(&feature);
  SetNearbyLanes(&feature);
  SetLaneGraphFeature(&feature);
  UpdateKFLaneTrackers(&feature);

  // Insert obstacle feature to history
  InsertFeatureToHistory(feature);

  // Set obstacle motion status
  SetMotionStatus();

  // Trim historical features
  Trim();
}

void Obstacle::SetStatus(const PerceptionObstacle& perception_obstacle,
                         const double timestamp, Feature* feature) {
  SetTimestamp(perception_obstacle, timestamp, feature);
  SetPosition(perception_obstacle, feature);
  SetVelocity(perception_obstacle, feature);
  SetAcceleration(feature);
  SetTheta(perception_obstacle, feature);
  SetLengthWidthHeight(perception_obstacle, feature);
}

void Obstacle::UpdateStatus(Feature* feature) {
  // Update motion belief
  if (!kf_motion_tracker_.IsInitialized()) {
    ADEBUG << "Obstacle [" << id_ << "] has not initialized motion tracker.";
    return;
  }
  auto state = kf_motion_tracker_.GetStateEstimate();

  feature->mutable_t_position()->set_x(state(0, 0));
  feature->mutable_t_position()->set_y(state(1, 0));
  feature->mutable_t_position()->set_z(feature->position().z());

  double velocity_x = state(2, 0);
  double velocity_y = state(3, 0);
  double speed = std::hypot(velocity_x, velocity_y);
  double velocity_heading = std::atan2(velocity_y, velocity_x);
  if (FLAGS_enable_adjust_velocity_heading) {
    UpdateVelocity(feature->theta(), &velocity_x, &velocity_y,
                   &velocity_heading, &speed);
  }
  feature->mutable_velocity()->set_x(velocity_x);
  feature->mutable_velocity()->set_y(velocity_y);
  feature->mutable_velocity()->set_z(velocity_heading);
  feature->set_speed(speed);
  feature->set_velocity_heading(std::atan2(state(3, 0), state(2, 0)));

  double acc_x = common::math::Clamp(state(4, 0), FLAGS_min_acc, FLAGS_max_acc);
  double acc_y = common::math::Clamp(state(5, 0), FLAGS_min_acc, FLAGS_max_acc);
  double acc =
      acc_x * std::cos(velocity_heading) + acc_y * std::sin(velocity_heading);
  feature->mutable_acceleration()->set_x(acc_x);
  feature->mutable_acceleration()->set_y(acc_y);
  feature->mutable_acceleration()->set_z(feature->acceleration().z());
  feature->set_acc(acc);

  ADEBUG << "Obstacle [" << id_ << "] has tracked position [" << std::fixed
         << std::setprecision(6) << feature->t_position().x() << ", "
         << std::fixed << std::setprecision(6) << feature->t_position().y()
         << ", " << std::fixed << std::setprecision(6)
         << feature->t_position().z() << "]";
  ADEBUG << "Obstacle [" << id_ << "] has tracked velocity [" << std::fixed
         << std::setprecision(6) << feature->velocity().x() << ", "
         << std::fixed << std::setprecision(6) << feature->velocity().y()
         << ", " << std::fixed << std::setprecision(6)
         << feature->velocity().z() << "]";
  ADEBUG << "Obstacle [" << id_ << "] has tracked acceleration [" << std::fixed
         << std::setprecision(6) << feature->acceleration().x() << ", "
         << std::fixed << std::setprecision(6) << feature->acceleration().y()
         << ", " << std::fixed << std::setprecision(6)
         << feature->acceleration().z() << "]";
  ADEBUG << "Obstacle [" << id_ << "] has tracked velocity heading ["
         << std::fixed << std::setprecision(6) << feature->velocity_heading()
         << "].";

  // Update pedestrian motion belief
  if (type_ == PerceptionObstacle::PEDESTRIAN) {
    if (!kf_pedestrian_tracker_.IsInitialized()) {
      ADEBUG << "Obstacle [" << id_
             << "] has not initialized pedestrian tracker.";
      return;
    }
    feature->mutable_t_position()->set_x(
        kf_pedestrian_tracker_.GetStateEstimate()(0, 0));
    feature->mutable_t_position()->set_y(
        kf_pedestrian_tracker_.GetStateEstimate()(1, 0));
    ADEBUG << "Obstacle [" << id_ << "] has tracked pedestrian position ["
           << std::setprecision(6) << feature->t_position().x() << std::fixed
           << ", " << std::setprecision(6) << feature->t_position().y()
           << std::fixed << ", " << std::setprecision(6)
           << feature->t_position().z() << std::fixed << "]";
  }
}

ErrorCode Obstacle::SetId(const PerceptionObstacle& perception_obstacle,
                          Feature* feature) {
  if (!perception_obstacle.has_id()) {
    AERROR << "Obstacle has no ID.";
    return ErrorCode::PREDICTION_ERROR;
  }

  int id = perception_obstacle.id();
  if (id_ < 0) {
    id_ = id;
    ADEBUG << "Obstacle has id [" << id_ << "].";
  } else {
    if (id_ != id) {
      AERROR << "Obstacle [" << id_ << "] has a mismatched ID [" << id
             << "] from perception obstacle.";
      return ErrorCode::PREDICTION_ERROR;
    } else {
    }
  }
  return ErrorCode::OK;
}

ErrorCode Obstacle::SetType(const PerceptionObstacle& perception_obstacle) {
  if (perception_obstacle.has_type()) {
    type_ = perception_obstacle.type();
    ADEBUG << "Obstacle [" << id_ << "] has type [" << type_ << "].";
  } else {
    AERROR << "Obstacle [" << id_ << "] has no type.";
    return ErrorCode::PREDICTION_ERROR;
  }
  return ErrorCode::OK;
}

void Obstacle::SetTimestamp(const PerceptionObstacle& perception_obstacle,
                            const double timestamp, Feature* feature) {
  double ts = timestamp;
  if (perception_obstacle.has_timestamp() &&
      perception_obstacle.timestamp() > 0.0) {
    ts = perception_obstacle.timestamp();
  }
  feature->set_timestamp(ts);

  ADEBUG << "Obstacle [" << id_ << "] has timestamp [" << std::fixed
         << std::setprecision(6) << ts << "].";
}

void Obstacle::SetPosition(const PerceptionObstacle& perception_obstacle,
                           Feature* feature) {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  if (perception_obstacle.has_position()) {
    if (perception_obstacle.position().has_x()) {
      x = perception_obstacle.position().x();
    }
    if (perception_obstacle.position().has_y()) {
      y = perception_obstacle.position().y();
    }
    if (perception_obstacle.position().has_z()) {
      z = perception_obstacle.position().z();
    }
  }

  feature->mutable_position()->set_x(x);
  feature->mutable_position()->set_y(y);
  feature->mutable_position()->set_z(z);

  ADEBUG << "Obstacle [" << id_ << "] has position [" << std::fixed
         << std::setprecision(6) << x << ", " << std::fixed
         << std::setprecision(6) << y << ", " << std::fixed
         << std::setprecision(6) << z << "].";
}

void Obstacle::SetVelocity(const PerceptionObstacle& perception_obstacle,
                           Feature* feature) {
  double velocity_x = 0.0;
  double velocity_y = 0.0;
  double velocity_z = 0.0;

  if (perception_obstacle.has_velocity()) {
    if (perception_obstacle.velocity().has_x()) {
      velocity_x = perception_obstacle.velocity().x();
    }
    if (perception_obstacle.velocity().has_y()) {
      velocity_y = perception_obstacle.velocity().y();
    }
    if (perception_obstacle.velocity().has_z()) {
      velocity_z = perception_obstacle.velocity().z();
    }
  }

  double speed = std::hypot(velocity_x, velocity_y);
  double velocity_heading = perception_obstacle.theta();

  if (FLAGS_enable_adjust_velocity_heading && history_size() > 0) {
    double diff_x =
        feature->position().x() - feature_history_.front().position().x();
    double diff_y =
        feature->position().y() - feature_history_.front().position().y();
    double prev_obstacle_size = std::max(feature_history_.front().length(),
                                         feature_history_.front().width());
    double obstacle_size =
        std::max(perception_obstacle.length(), perception_obstacle.width());
    double size_diff = std::abs(obstacle_size - prev_obstacle_size);
    double shift_thred =
        std::max(obstacle_size * FLAGS_valid_position_diff_rate_threshold,
                 FLAGS_valid_position_diff_threshold);
    double size_diff_thred =
        FLAGS_split_rate * std::min(obstacle_size, prev_obstacle_size);
    if (std::fabs(diff_x) > shift_thred && std::fabs(diff_y) > shift_thred &&
        size_diff < size_diff_thred) {
      double shift_heading = std::atan2(diff_y, diff_x);
      double angle_diff = ::apollo::common::math::NormalizeAngle(
          shift_heading - velocity_heading);
      if (std::fabs(angle_diff) > FLAGS_max_lane_angle_diff) {
        velocity_heading = shift_heading;
      }
    }
    double filtered_heading = heading_filter_.Filter(velocity_heading);
    if (type_ == PerceptionObstacle::BICYCLE ||
        type_ == PerceptionObstacle::PEDESTRIAN) {
      velocity_heading = filtered_heading;
    }
    velocity_x = speed * std::cos(velocity_heading);
    velocity_y = speed * std::sin(velocity_heading);
  }

  feature->mutable_velocity()->set_x(velocity_x);
  feature->mutable_velocity()->set_y(velocity_y);
  feature->mutable_velocity()->set_z(velocity_z);
  feature->set_velocity_heading(velocity_heading);
  feature->set_speed(speed);

  ADEBUG << "Obstacle [" << id_ << "] has velocity [" << std::fixed
         << std::setprecision(6) << velocity_x << ", " << std::fixed
         << std::setprecision(6) << velocity_y << ", " << std::fixed
         << std::setprecision(6) << velocity_z << "]";
  ADEBUG << "Obstacle [" << id_ << "] has velocity heading [" << std::fixed
         << std::setprecision(6) << velocity_heading << "] ";
  ADEBUG << "Obstacle [" << id_ << "] has speed [" << std::fixed
         << std::setprecision(6) << speed << "].";
}

void Obstacle::UpdateVelocity(const double theta, double* velocity_x,
                              double* velocity_y, double* velocity_heading,
                              double* speed) {
  *speed = std::hypot(*velocity_x, *velocity_y);
  double angle_diff =
      ::apollo::common::math::NormalizeAngle(*velocity_heading - theta);
  if (std::fabs(angle_diff) <= FLAGS_max_lane_angle_diff) {
    *velocity_heading = theta;
    *velocity_x = *speed * std::cos(*velocity_heading);
    *velocity_y = *speed * std::sin(*velocity_heading);
  }
}

void Obstacle::SetAcceleration(Feature* feature) {
  double acc_x = 0.0;
  double acc_y = 0.0;
  double acc_z = 0.0;
  double acc = 0.0;

  if (feature_history_.size() > 0) {
    double curr_ts = feature->timestamp();
    double prev_ts = feature_history_.front().timestamp();

    const Point3D& curr_velocity = feature->velocity();
    const Point3D& prev_velocity = feature_history_.front().velocity();

    if (curr_ts > prev_ts) {
      /*
       * A damp function is to punish acc calculation for low speed
       * and reward it for high speed
       */
      double damping_x = Damp(curr_velocity.x(), 0.001);
      double damping_y = Damp(curr_velocity.y(), 0.001);
      double damping_z = Damp(curr_velocity.z(), 0.001);

      acc_x = (curr_velocity.x() - prev_velocity.x()) / (curr_ts - prev_ts);
      acc_y = (curr_velocity.y() - prev_velocity.y()) / (curr_ts - prev_ts);
      acc_z = (curr_velocity.z() - prev_velocity.z()) / (curr_ts - prev_ts);

      acc_x *= damping_x;
      acc_y *= damping_y;
      acc_z *= damping_z;

      acc_x =
          common::math::Clamp(acc_x * damping_x, FLAGS_min_acc, FLAGS_max_acc);
      acc_y =
          common::math::Clamp(acc_y * damping_y, FLAGS_min_acc, FLAGS_max_acc);
      acc_z =
          common::math::Clamp(acc_z * damping_z, FLAGS_min_acc, FLAGS_max_acc);

      double heading = feature->velocity_heading();
      acc = acc_x * std::cos(heading) + acc_y * std::sin(heading);
    }
  }

  feature->mutable_acceleration()->set_x(acc_x);
  feature->mutable_acceleration()->set_y(acc_y);
  feature->mutable_acceleration()->set_z(acc_z);
  feature->set_acc(acc);

  ADEBUG << "Obstacle [" << id_ << "] has acceleration [" << std::fixed
         << std::setprecision(6) << acc_x << ", " << std::fixed
         << std::setprecision(6) << acc_y << ", " << std::fixed
         << std::setprecision(6) << acc_z << "]";
  ADEBUG << "Obstacle [" << id_ << "] has acceleration value [" << std::fixed
         << std::setprecision(6) << acc << "].";
}

void Obstacle::SetTheta(const PerceptionObstacle& perception_obstacle,
                        Feature* feature) {
  double theta = 0.0;
  if (perception_obstacle.has_theta()) {
    theta = perception_obstacle.theta();
  }
  feature->set_theta(theta);

  ADEBUG << "Obstacle [" << id_ << "] has theta [" << std::fixed
         << std::setprecision(6) << theta << "].";
}

void Obstacle::SetLengthWidthHeight(
    const PerceptionObstacle& perception_obstacle, Feature* feature) {
  double length = 0.0;
  double width = 0.0;
  double height = 0.0;

  if (perception_obstacle.has_length()) {
    length = perception_obstacle.length();
  }
  if (perception_obstacle.has_width()) {
    width = perception_obstacle.width();
  }
  if (perception_obstacle.has_height()) {
    height = perception_obstacle.height();
  }

  feature->set_length(length);
  feature->set_width(width);
  feature->set_height(height);

  ADEBUG << "Obstacle [" << id_ << "] has dimension [" << std::fixed
         << std::setprecision(6) << length << ", " << std::fixed
         << std::setprecision(6) << width << ", " << std::fixed
         << std::setprecision(6) << height << "].";
}

void Obstacle::InitKFMotionTracker(const Feature& feature) {
  // Set transition matrix F
  // constant acceleration dynamic model
  Eigen::Matrix<double, 6, 6> F;
  F.setIdentity();
  kf_motion_tracker_.SetTransitionMatrix(F);

  // Set observation matrix H
  Eigen::Matrix<double, 2, 6> H;
  H.setIdentity();
  kf_motion_tracker_.SetObservationMatrix(H);

  // Set covariance of transition noise matrix Q
  // make the noise this order:
  // noise(x/y) < noise(vx/vy) < noise(ax/ay)
  Eigen::Matrix<double, 6, 6> Q;
  Q.setIdentity();
  Q *= FLAGS_q_var;
  kf_motion_tracker_.SetTransitionNoise(Q);

  // Set observation noise matrix R
  Eigen::Matrix<double, 2, 2> R;
  R.setIdentity();
  R *= FLAGS_r_var;
  kf_motion_tracker_.SetObservationNoise(R);

  // Set current state covariance matrix P
  // make the covariance this order:
  // cov(x/y) < cov(vx/vy) < cov(ax/ay)
  Eigen::Matrix<double, 6, 6> P;
  P.setIdentity();
  P *= FLAGS_p_var;

  // Set initial state
  Eigen::Matrix<double, 6, 1> x;
  x(0, 0) = feature.position().x();
  x(1, 0) = feature.position().y();
  x(2, 0) = feature.velocity().x();
  x(3, 0) = feature.velocity().y();
  x(4, 0) = feature.acceleration().x();
  x(5, 0) = feature.acceleration().y();

  kf_motion_tracker_.SetStateEstimate(x, P);
}

void Obstacle::UpdateKFMotionTracker(const Feature& feature) {
  double delta_ts = 0.0;
  if (feature_history_.size() > 0) {
    delta_ts = feature.timestamp() - feature_history_.front().timestamp();
  }
  if (delta_ts > FLAGS_double_precision) {
    // Set tansition matrix and predict
    auto F = kf_motion_tracker_.GetTransitionMatrix();
    F(0, 2) = delta_ts;
    F(0, 4) = 0.5 * delta_ts * delta_ts;
    F(1, 3) = delta_ts;
    F(1, 5) = 0.5 * delta_ts * delta_ts;
    F(2, 4) = delta_ts;
    F(3, 5) = delta_ts;
    kf_motion_tracker_.SetTransitionMatrix(F);
    kf_motion_tracker_.Predict();

    // Set observation and correct
    Eigen::Matrix<double, 2, 1> z;
    z(0, 0) = feature.position().x();
    z(1, 0) = feature.position().y();
    kf_motion_tracker_.Correct(z);
  }
}

void Obstacle::InitKFLaneTracker(const std::string& lane_id,
                                 const double beta) {
  KalmanFilter<double, 4, 2, 0> kf;

  // transition matrix: update delta_t at each processing step
  Eigen::Matrix<double, 4, 4> F;
  F.setIdentity();
  F(1, 1) = beta;
  kf.SetTransitionMatrix(F);

  // observation matrix
  Eigen::Matrix<double, 2, 4> H;
  H.setZero();
  H(0, 0) = 1.0;
  H(1, 1) = 1.0;
  kf.SetObservationMatrix(H);

  // Set covariance of transition noise matrix Q
  Eigen::Matrix<double, 4, 4> Q;
  Q.setIdentity();
  Q *= FLAGS_q_var;
  kf.SetTransitionNoise(Q);

  // Set observation noise matrix R
  Eigen::Matrix<double, 2, 2> R;
  R.setIdentity();
  R *= FLAGS_r_var;
  kf.SetObservationNoise(R);

  // Set current state covariance matrix P
  Eigen::Matrix<double, 4, 4> P;
  P.setIdentity();
  P *= FLAGS_p_var;
  kf.SetStateCovariance(P);

  kf_lane_trackers_.emplace(lane_id, std::move(kf));
}

void Obstacle::UpdateKFLaneTrackers(Feature* feature) {
  if (!feature->has_lane()) {
    return;
  }
  std::unordered_set<std::string> lane_ids;
  for (auto& lane_feature : feature->lane().current_lane_feature()) {
    if (!lane_feature.lane_id().empty()) {
      lane_ids.insert(lane_feature.lane_id());
    }
  }
  for (auto& lane_feature : feature->lane().nearby_lane_feature()) {
    if (!lane_feature.lane_id().empty()) {
      lane_ids.insert(lane_feature.lane_id());
    }
  }
  if (lane_ids.empty()) {
    return;
  }

  auto iter = kf_lane_trackers_.begin();
  while (iter != kf_lane_trackers_.end()) {
    if (lane_ids.find(iter->first) == lane_ids.end()) {
      iter = kf_lane_trackers_.erase(iter);
    } else {
      ++iter;
    }
  }

  double ts = feature->timestamp();
  double speed = feature->speed();
  double acc = feature->acc();
  for (auto& lane_feature : feature->lane().current_lane_feature()) {
    std::string lane_id = lane_feature.lane_id();
    if (lane_id.empty()) {
      continue;
    }
    double s = lane_feature.lane_s();
    double l = lane_feature.lane_l();
    UpdateKFLaneTracker(lane_id, s, l, speed, acc, ts, FLAGS_go_approach_rate);
  }
  for (auto& nearby_lane_feature : feature->lane().nearby_lane_feature()) {
    std::string lane_id = nearby_lane_feature.lane_id();
    if (lane_id.empty()) {
      continue;
    }
    double s = nearby_lane_feature.lane_s();
    double l = nearby_lane_feature.lane_l();
    UpdateKFLaneTracker(lane_id, s, l, speed, acc, ts,
                        FLAGS_cutin_approach_rate);
  }

  if (FLAGS_enable_kf_tracking && id_ >= 0) {
    UpdateLaneBelief(feature);
  }
}

void Obstacle::UpdateKFLaneTracker(const std::string& lane_id,
                                   const double lane_s, const double lane_l,
                                   const double lane_speed,
                                   const double lane_acc,
                                   const double timestamp, const double beta) {
  auto* kf_ptr = FindOrNull(kf_lane_trackers_, lane_id);
  if (kf_ptr != nullptr) {
    double delta_ts = 0.0;
    if (feature_history_.size() > 0) {
      delta_ts = timestamp - feature_history_.front().timestamp();
    }
    if (delta_ts > FLAGS_double_precision) {
      auto F = kf_ptr->GetTransitionMatrix();
      F(0, 2) = delta_ts;
      F(0, 3) = 0.5 * delta_ts * delta_ts;
      F(2, 3) = delta_ts;
      kf_ptr->SetTransitionMatrix(F);
      kf_ptr->Predict();

      Eigen::Matrix<double, 2, 1> z;
      z(0, 0) = lane_s;
      z(1, 0) = lane_l;
      kf_ptr->Correct(z);
    }
  } else {
    InitKFLaneTracker(lane_id, beta);
    auto& kf = FindOrDie(kf_lane_trackers_, lane_id);
    Eigen::Matrix<double, 4, 1> state;
    state.setZero();
    state(0, 0) = lane_s;
    state(1, 0) = lane_l;
    state(2, 0) = lane_speed;
    state(3, 0) = lane_acc;

    auto P = kf.GetStateCovariance();
    kf.SetStateEstimate(state, P);
  }
}

void Obstacle::UpdateLaneBelief(Feature* feature) {
  if ((!feature->has_lane()) || (!feature->lane().has_lane_feature())) {
    return;
  }
  const std::string& lane_id = feature->lane().lane_feature().lane_id();
  if (lane_id.empty()) {
    return;
  }

  auto* kf_ptr = FindOrNull(kf_lane_trackers_, lane_id);
  if (kf_ptr == nullptr) {
    return;
  }

  double lane_speed = kf_ptr->GetStateEstimate()(2, 0);
  double lane_acc = common::math::Clamp(kf_ptr->GetStateEstimate()(3, 0),
                                        FLAGS_min_acc, FLAGS_max_acc);

  ADEBUG << "Obstacle [" << id_ << "] has tracked lane speed [" << std::fixed
         << std::setprecision(6) << lane_speed << "]";
  ADEBUG << "Obstacle [" << id_ << "] has tracked lane acceleration ["
         << std::fixed << std::setprecision(6) << lane_acc << "]";
}

void Obstacle::InitKFPedestrianTracker(const Feature& feature) {
  // Set transition matrix F
  Eigen::Matrix<double, 2, 2> F;
  F.setIdentity();
  kf_pedestrian_tracker_.SetTransitionMatrix(F);

  // Set observation matrix H
  Eigen::Matrix<double, 2, 2> H;
  H.setIdentity();
  kf_pedestrian_tracker_.SetObservationMatrix(H);

  // Set control matrix
  Eigen::Matrix<double, 2, 4> B;
  B.setZero();
  kf_pedestrian_tracker_.SetControlMatrix(B);

  // Set covariance of transition noise matrix Q
  Eigen::Matrix<double, 2, 2> Q;
  Q.setIdentity();
  Q *= FLAGS_q_var;
  kf_pedestrian_tracker_.SetTransitionNoise(Q);

  // Set observation noise matrix R
  Eigen::Matrix<double, 2, 2> R;
  R.setIdentity();
  R *= FLAGS_r_var;
  kf_pedestrian_tracker_.SetObservationNoise(R);

  // Set current state covariance matrix P
  Eigen::Matrix<double, 2, 2> P;
  P.setIdentity();
  P *= FLAGS_p_var;

  // Set initial state
  Eigen::Matrix<double, 2, 1> x;
  x.setZero();
  x(0, 0) = feature.position().x();
  x(1, 0) = feature.position().y();

  kf_pedestrian_tracker_.SetStateEstimate(x, P);
}

void Obstacle::UpdateKFPedestrianTracker(const Feature& feature) {
  double delta_ts = 0.0;
  if (!feature_history_.empty()) {
    delta_ts = feature.timestamp() - feature_history_.front().timestamp();
  }
  if (delta_ts > std::numeric_limits<double>::epsilon()) {
    Eigen::Matrix<double, 2, 4> B = kf_pedestrian_tracker_.GetControlMatrix();
    B(0, 0) = delta_ts;
    B(0, 2) = 0.5 * delta_ts * delta_ts;
    B(1, 1) = delta_ts;
    B(1, 3) = 0.5 * delta_ts * delta_ts;
    kf_pedestrian_tracker_.SetControlMatrix(B);

    // Set control vector
    Eigen::Matrix<double, 4, 1> u;
    u(0, 0) = feature.velocity().x();
    u(1, 0) = feature.velocity().y();
    if (FLAGS_enable_pedestrian_acc) {
      u(2, 0) = feature.acceleration().x();
      u(3, 0) = feature.acceleration().y();
    }

    kf_pedestrian_tracker_.Predict(u);

    // Set observation vector
    Eigen::Matrix<double, 2, 1> z;
    z(0, 0) = feature.position().x();
    z(1, 0) = feature.position().y();
    kf_pedestrian_tracker_.Correct(z);
  }
}

void Obstacle::SetCurrentLanes(Feature* feature) {
  PredictionMap* map = PredictionMap::instance();

  Eigen::Vector2d point(feature->position().x(), feature->position().y());
  double heading = feature->velocity_heading();
  std::vector<std::shared_ptr<const LaneInfo>> current_lanes;
  map->OnLane(current_lanes_, point, heading, FLAGS_lane_search_radius, true,
              &current_lanes);
  current_lanes_ = current_lanes;
  if (current_lanes_.empty()) {
    ADEBUG << "Obstacle [" << id_ << "] has no current lanes.";
    kf_lane_trackers_.clear();
    return;
  }
  Lane lane;
  if (feature->has_lane()) {
    lane = feature->lane();
  }
  double min_heading_diff = std::numeric_limits<double>::infinity();
  for (std::shared_ptr<const LaneInfo> current_lane : current_lanes) {
    if (current_lane == nullptr) {
      continue;
    }

    int turn_type = map->LaneTurnType(current_lane->id().id());
    std::string lane_id = current_lane->id().id();
    double s = 0.0;
    double l = 0.0;
    map->GetProjection(point, current_lane, &s, &l);
    if (s < 0.0) {
      continue;
    }

    common::math::Vec2d vec_point(point[0], point[1]);
    double distance = 0.0;
    common::PointENU nearest_point =
        current_lane->GetNearestPoint(vec_point, &distance);
    double nearest_point_heading =
        map->PathHeading(current_lane, nearest_point);
    double angle_diff = common::math::AngleDiff(heading, nearest_point_heading);
    double left = 0.0;
    double right = 0.0;
    current_lane->GetWidth(s, &left, &right);
    LaneFeature* lane_feature = lane.add_current_lane_feature();
    lane_feature->set_lane_turn_type(turn_type);
    lane_feature->set_lane_id(lane_id);
    lane_feature->set_lane_s(s);
    lane_feature->set_lane_l(l);
    lane_feature->set_angle_diff(angle_diff);
    lane_feature->set_dist_to_left_boundary(left - l);
    lane_feature->set_dist_to_right_boundary(right + l);
    if (std::fabs(angle_diff) < min_heading_diff) {
      lane.mutable_lane_feature()->CopyFrom(*lane_feature);
      min_heading_diff = std::fabs(angle_diff);
    }
    ADEBUG << "Obstacle [" << id_ << "] has current lanes ["
           << lane_feature->ShortDebugString() << "].";
  }

  if (lane.has_lane_feature()) {
    ADEBUG << "Obstacle [" << id_ << "] has one current lane ["
           << lane.lane_feature().ShortDebugString() << "].";
  }

  feature->mutable_lane()->CopyFrom(lane);
}

void Obstacle::SetNearbyLanes(Feature* feature) {
  PredictionMap* map = PredictionMap::instance();

  Eigen::Vector2d point(feature->position().x(), feature->position().y());
  double theta = feature->velocity_heading();
  std::vector<std::shared_ptr<const LaneInfo>> nearby_lanes;
  map->NearbyLanesByCurrentLanes(point, theta, FLAGS_lane_search_radius,
                                 current_lanes_, &nearby_lanes);
  if (nearby_lanes.empty()) {
    ADEBUG << "Obstacle [" << id_ << "] has no nearby lanes.";
    return;
  }

  for (std::shared_ptr<const LaneInfo> nearby_lane : nearby_lanes) {
    if (nearby_lane == nullptr) {
      continue;
    }

    // Ignore bike and sidewalk lanes for vehicles
    if (type_ == PerceptionObstacle::VEHICLE &&
        nearby_lane->lane().has_type() &&
        (nearby_lane->lane().type() == ::apollo::hdmap::Lane::BIKING ||
         nearby_lane->lane().type() == ::apollo::hdmap::Lane::SIDEWALK)) {
      ADEBUG << "Obstacle [" << id_ << "] ignores disqualified lanes.";
      continue;
    }

    double s = -1.0;
    double l = 0.0;
    map->GetProjection(point, nearby_lane, &s, &l);
    if (s < 0.0 || s >= nearby_lane->total_length()) {
      continue;
    }
    int turn_type = map->LaneTurnType(nearby_lane->id().id());
    double heading = feature->velocity_heading();
    double angle_diff = 0.0;
    hdmap::MapPathPoint nearest_point;
    if (!map->ProjectionFromLane(nearby_lane, s, &nearest_point)) {
      angle_diff = common::math::AngleDiff(nearest_point.heading(), heading);
    }

    double left = 0.0;
    double right = 0.0;
    nearby_lane->GetWidth(s, &left, &right);

    LaneFeature* lane_feature =
        feature->mutable_lane()->add_nearby_lane_feature();

    lane_feature->set_lane_turn_type(turn_type);
    lane_feature->set_lane_id(nearby_lane->id().id());
    lane_feature->set_lane_s(s);
    lane_feature->set_lane_l(l);
    lane_feature->set_angle_diff(angle_diff);
    lane_feature->set_dist_to_left_boundary(left - l);
    lane_feature->set_dist_to_right_boundary(right + l);
    ADEBUG << "Obstacle [" << id_ << "] has nearby lanes ["
           << lane_feature->ShortDebugString() << "]";
  }
}

void Obstacle::SetLaneGraphFeature(Feature* feature) {
  PredictionMap* map = PredictionMap::instance();
  double speed = feature->speed();
  double acc = feature->acc();
  double road_graph_distance =
      speed * FLAGS_prediction_duration +
      0.5 * acc * FLAGS_prediction_duration * FLAGS_prediction_duration +
      FLAGS_min_prediction_length;
  int curr_lane_count = 0;
  for (auto& lane : feature->lane().current_lane_feature()) {
    std::shared_ptr<const LaneInfo> lane_info = map->LaneById(lane.lane_id());
    RoadGraph road_graph(lane.lane_s(), road_graph_distance, lane_info);
    LaneGraph lane_graph;
    road_graph.BuildLaneGraph(&lane_graph);
    if (lane_graph.lane_sequence_size() > 0) {
      ++curr_lane_count;
    }
    int seq_id =
        feature->mutable_lane()->mutable_lane_graph()->lane_sequence_size();
    for (const auto& lane_seq : lane_graph.lane_sequence()) {
      feature->mutable_lane()
          ->mutable_lane_graph()
          ->add_lane_sequence()
          ->CopyFrom(lane_seq);
      feature->mutable_lane()
          ->mutable_lane_graph()
          ->mutable_lane_sequence(seq_id)
          ->set_lane_sequence_id(seq_id);
      ++seq_id;
      ADEBUG << "Obstacle [" << id_ << "] set a lane sequence ["
             << lane_seq.ShortDebugString() << "].";
    }
    if (curr_lane_count >= FLAGS_max_num_current_lane) {
      break;
    }
  }

  int nearby_lane_count = 0;
  for (auto& lane : feature->lane().nearby_lane_feature()) {
    std::shared_ptr<const LaneInfo> lane_info = map->LaneById(lane.lane_id());
    RoadGraph road_graph(lane.lane_s(), road_graph_distance, lane_info);
    LaneGraph lane_graph;
    road_graph.BuildLaneGraph(&lane_graph);
    if (lane_graph.lane_sequence_size() > 0) {
      ++nearby_lane_count;
    }
    int seq_id =
        feature->mutable_lane()->mutable_lane_graph()->lane_sequence_size();
    for (const auto& lane_seq : lane_graph.lane_sequence()) {
      feature->mutable_lane()
          ->mutable_lane_graph()
          ->add_lane_sequence()
          ->CopyFrom(lane_seq);
      feature->mutable_lane()
          ->mutable_lane_graph()
          ->mutable_lane_sequence(seq_id)
          ->set_lane_sequence_id(seq_id);
      ADEBUG << "Obstacle [" << id_ << "] set a lane sequence ["
             << lane_seq.ShortDebugString() << "].";
    }
    if (nearby_lane_count >= FLAGS_max_num_nearby_lane) {
      break;
    }
  }

  if (feature->has_lane() && feature->lane().has_lane_graph()) {
    SetLanePoints(feature);
  }

  ADEBUG << "Obstacle [" << id_ << "] set lane graph features.";
}

void Obstacle::SetLanePoints(Feature* feature) {
  if (feature == nullptr || !feature->velocity_heading()) {
    AERROR << "Null feature or no velocity heading.";
    return;
  }
  PredictionMap* map = PredictionMap::instance();

  LaneGraph* lane_graph = feature->mutable_lane()->mutable_lane_graph();
  double heading = feature->velocity_heading();
  double x = feature->position().x();
  double y = feature->position().y();
  Eigen::Vector2d position(x, y);
  for (int i = 0; i < lane_graph->lane_sequence_size(); ++i) {
    LaneSequence* lane_sequence = lane_graph->mutable_lane_sequence(i);
    if (lane_sequence->lane_segment_size() <= 0) {
      continue;
    }
    int lane_index = 0;
    LaneSegment* lane_segment = lane_sequence->mutable_lane_segment(lane_index);
    double start_s = lane_sequence->lane_segment(lane_index).start_s();
    double total_s = 0.0;
    double lane_seg_s = start_s;
    while (lane_index < lane_sequence->lane_segment_size()) {
      if (lane_seg_s > lane_segment->end_s()) {
        start_s = lane_seg_s - lane_segment->end_s();
        lane_seg_s = start_s;
        ++lane_index;
        if (lane_index < lane_sequence->lane_segment_size()) {
          lane_segment = lane_sequence->mutable_lane_segment(lane_index);
        } else {
          lane_segment = nullptr;
        }
      } else {
        std::string lane_id = lane_segment->lane_id();
        std::shared_ptr<const LaneInfo> lane_info = map->LaneById(lane_id);
        if (lane_info == nullptr) {
          break;
        }
        LanePoint lane_point;
        Eigen::Vector2d lane_point_pos =
            map->PositionOnLane(lane_info, lane_seg_s);
        double lane_point_heading = map->HeadingOnLane(lane_info, lane_seg_s);
        double lane_point_width = map->LaneTotalWidth(lane_info, lane_seg_s);
        double lane_point_angle_diff =
            common::math::AngleDiff(lane_point_heading, heading);
        lane_point.mutable_position()->set_x(lane_point_pos[0]);
        lane_point.mutable_position()->set_y(lane_point_pos[1]);
        lane_point.set_heading(lane_point_heading);
        lane_point.set_width(lane_point_width);
        double lane_s = -1.0;
        double lane_l = 0.0;
        map->GetProjection(position, lane_info, &lane_s, &lane_l);
        lane_point.set_relative_s(total_s);
        lane_point.set_relative_l(0.0 - lane_l);
        lane_point.set_angle_diff(lane_point_angle_diff);
        lane_segment->set_lane_turn_type(map->LaneTurnType(lane_id));
        lane_segment->add_lane_point()->CopyFrom(lane_point);
        total_s += FLAGS_target_lane_gap;
        lane_seg_s += FLAGS_target_lane_gap;
      }
    }
  }
  ADEBUG << "Obstacle [" << id_ << "] has lane segments and points.";
}

void Obstacle::SetMotionStatus() {
  int history_size = static_cast<int>(feature_history_.size());
  if (history_size < 2) {
    ADEBUG << "Obstacle [" << id_ << "] has no history and "
           << "is considered moving.";
    if (history_size > 0) {
      feature_history_.front().set_is_still(false);
    }
    return;
  }

  double start_x = 0.0;
  double start_y = 0.0;
  double avg_drift_x = 0.0;
  double avg_drift_y = 0.0;
  int len = std::min(history_size, FLAGS_still_obstacle_history_length);
  CHECK_GT(len, 1);

  auto feature_riter = feature_history_.rbegin();
  start_x = feature_riter->position().x();
  start_y = feature_riter->position().y();
  ++feature_riter;
  while (feature_riter != feature_history_.rend()) {
    avg_drift_x += (feature_riter->position().x() - start_x) / (len - 1);
    avg_drift_y += (feature_riter->position().y() - start_y) / (len - 1);
    ++feature_riter;
  }

  double delta_ts = feature_history_.front().timestamp() -
                    feature_history_.back().timestamp();
  double std = FLAGS_still_obstacle_position_std;
  double speed_sensibility =
      std::sqrt(2 * history_size) * 4 * std / ((history_size + 1) * delta_ts);
  double speed = feature_history_.front().speed();
  double speed_threshold = FLAGS_still_obstacle_speed_threshold;
  if (type_ == PerceptionObstacle::PEDESTRIAN ||
      type_ == PerceptionObstacle::BICYCLE) {
    speed_threshold = FLAGS_still_pedestrian_speed_threshold;
  }
  if (speed < speed_threshold) {
    ADEBUG << "Obstacle [" << id_ << "] has a small speed [" << speed
           << "] and is considered stationary.";
    feature_history_.front().set_is_still(true);
  } else if (speed_sensibility < speed_threshold) {
    ADEBUG << "Obstacle [" << id_ << "]"
           << "] considered moving [sensibility = " << speed_sensibility << "]";
    feature_history_.front().set_is_still(false);
  } else {
    double distance = std::hypot(avg_drift_x, avg_drift_y);
    double distance_std = std::sqrt(2.0 / len) * std;
    if (distance > 2.0 * distance_std) {
      ADEBUG << "Obstacle [" << id_ << "] is moving.";
      feature_history_.front().set_is_still(false);
    } else {
      ADEBUG << "Obstacle [" << id_ << "] is stationary.";
      feature_history_.front().set_is_still(true);
    }
  }
}

void Obstacle::InsertFeatureToHistory(const Feature& feature) {
  feature_history_.emplace_front(feature);
  ADEBUG << "Obstacle [" << id_ << "] inserted a frame into the history.";
}

void Obstacle::Trim() {
  if (feature_history_.size() < 2) {
    return;
  }
  int count = 0;
  const double latest_ts = feature_history_.front().timestamp();
  while (!feature_history_.empty() &&
         latest_ts - feature_history_.back().timestamp() >=
             FLAGS_max_history_time) {
    feature_history_.pop_back();
    ++count;
  }
  if (count > 0) {
    ADEBUG << "Obstacle [" << id_ << "] trimmed " << count
           << " historical features";
  }
}

void Obstacle::SetRNNStates(const std::vector<Eigen::MatrixXf>& rnn_states) {
  rnn_states_ = rnn_states;
}

void Obstacle::GetRNNStates(std::vector<Eigen::MatrixXf>* rnn_states) {
  rnn_states->clear();
  rnn_states->insert(rnn_states->end(), rnn_states_.begin(), rnn_states_.end());
}

void Obstacle::InitRNNStates() {
  if (network::RnnModel::instance()->IsOk()) {
    network::RnnModel::instance()->ResetState();
    network::RnnModel::instance()->State(&rnn_states_);
    rnn_enabled_ = true;
    ADEBUG << "Success to initialize rnn model.";
  } else {
    AWARN << "Fail to initialize rnn model.";
    rnn_enabled_ = false;
  }
}

bool Obstacle::RNNEnabled() const { return rnn_enabled_; }

}  // namespace prediction
}  // namespace apollo
