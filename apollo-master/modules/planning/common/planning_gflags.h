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

#ifndef MODULES_PLANNING_COMMON_PLANNING_GFLAGS_H_
#define MODULES_PLANNING_COMMON_PLANNING_GFLAGS_H_

#include "gflags/gflags.h"

DECLARE_bool(planning_test_mode);
DECLARE_double(test_duration);

DECLARE_string(planning_config_file);
DECLARE_string(planning_adapter_config_filename);
DECLARE_int32(planning_loop_rate);
DECLARE_string(rtk_trajectory_filename);
DECLARE_uint64(rtk_trajectory_forward);
DECLARE_double(rtk_trajectory_resolution);
DECLARE_double(look_backward_distance);
DECLARE_double(look_forward_short_distance);
DECLARE_double(look_forward_long_distance);
DECLARE_double(look_forward_time_sec);
DECLARE_bool(enable_reference_line_stitching);
DECLARE_double(look_forward_extend_distance);
DECLARE_double(reference_line_stitch_overlap_distance);
DECLARE_double(reference_line_lateral_buffer);
DECLARE_double(prepare_rerouting_time);
DECLARE_double(rerouting_cooldown_time);

DECLARE_bool(enable_smooth_reference_line);
DECLARE_bool(enable_spiral_reference_line);
DECLARE_double(spiral_smoother_max_deviation);
DECLARE_int32(spiral_smoother_num_iteration);
DECLARE_double(spiral_smoother_piecewise_length);
DECLARE_double(spiral_reference_line_resolution);

DECLARE_bool(prioritize_change_lane);
DECLARE_bool(reckless_change_lane);
DECLARE_double(change_lane_fail_freeze_time);
DECLARE_double(change_lane_success_freeze_time);
DECLARE_double(change_lane_min_length);
DECLARE_bool(enable_change_lane_decider);
DECLARE_double(change_lane_speed_relax_percentage);
DECLARE_bool(enable_side_vehicle_st_boundary);

DECLARE_double(max_collision_distance);
DECLARE_bool(publish_estop);
DECLARE_bool(enable_trajectory_stitcher);

DECLARE_int32(max_history_frame_num);

// parameters for trajectory stitching and reinit planning starting point.
DECLARE_double(replan_lateral_distance_threshold);
DECLARE_double(replan_longitudinal_distance_threshold);
DECLARE_bool(estimate_current_vehicle_state);

// parameter for reference line
DECLARE_bool(enable_reference_line_provider_thread);
DECLARE_double(default_reference_line_width);
DECLARE_double(smoothed_reference_line_max_diff);

// parameters for trajectory planning
DECLARE_double(planning_upper_speed_limit);
DECLARE_double(trajectory_time_length);
DECLARE_double(trajectory_time_min_interval);
DECLARE_double(trajectory_time_max_interval);
DECLARE_double(trajectory_time_high_density_period);

// parameters for trajectory sanity check
DECLARE_bool(enable_trajectory_check);
DECLARE_double(speed_lower_bound);
DECLARE_double(speed_upper_bound);

DECLARE_double(longitudinal_acceleration_lower_bound);
DECLARE_double(longitudinal_acceleration_upper_bound);

DECLARE_double(lateral_jerk_bound);

DECLARE_double(longitudinal_jerk_lower_bound);
DECLARE_double(longitudinal_jerk_upper_bound);

DECLARE_double(dl_bound);
DECLARE_double(kappa_bound);
DECLARE_double(dkappa_bound);

// STBoundary
DECLARE_double(st_max_s);
DECLARE_double(st_max_t);

// Decision Part
DECLARE_double(static_obstacle_speed_threshold);
DECLARE_bool(enable_nudge_decision);
DECLARE_bool(enable_nudge_slowdown);
DECLARE_bool(try_history_decision);
DECLARE_double(static_decision_nudge_l_buffer);
DECLARE_double(lateral_ignore_buffer);
DECLARE_double(min_stop_distance_obstacle);
DECLARE_double(max_stop_distance_obstacle);
DECLARE_double(stop_distance_destination);
DECLARE_double(stop_distance_traffic_light);
DECLARE_double(destination_check_distance);
DECLARE_double(nudge_distance_obstacle);
DECLARE_double(follow_min_distance);
DECLARE_double(yield_min_distance);
DECLARE_double(follow_time_buffer);
DECLARE_double(follow_min_time_sec);
DECLARE_double(within_lane_bound);

DECLARE_string(destination_obstacle_id);
DECLARE_double(virtual_stop_wall_length);
DECLARE_double(virtual_stop_wall_height);
DECLARE_string(reference_line_end_obstacle_id);

DECLARE_double(prediction_total_time);
DECLARE_bool(align_prediction_time);
DECLARE_bool(enable_lag_prediction);
DECLARE_int32(lag_prediction_min_appear_num);
DECLARE_double(lag_prediction_max_disappear_num);
DECLARE_int32(trajectory_point_num_for_debug);
DECLARE_double(lag_prediction_protection_distance);
DECLARE_double(perception_confidence_threshold);

DECLARE_bool(enable_record_debug);
DECLARE_bool(enable_prediction);
DECLARE_bool(enable_traffic_light);

DECLARE_double(turn_signal_distance);
DECLARE_bool(right_turn_creep_forward);

// speed decider
DECLARE_double(low_speed_obstacle_threshold);
DECLARE_double(decelerating_obstacle_threshold);

// QpSt optimizer
DECLARE_bool(enable_slowdown_profile_generator);
DECLARE_double(slowdown_profile_deceleration);
DECLARE_bool(enable_follow_accel_constraint);

// traffic decision
/// common
DECLARE_double(stop_max_distance_buffer);
DECLARE_double(stop_min_speed);
DECLARE_double(stop_max_deceleration);
DECLARE_double(signal_expire_time_sec);
DECLARE_double(max_valid_stop_distance);

/// Clear Zone
DECLARE_string(clear_zone_virtual_object_id_prefix);
/// triffic light
DECLARE_string(signal_light_virtual_object_id_prefix);
DECLARE_double(max_deacceleration_for_yellow_light_stop);
/// crosswalk
DECLARE_bool(enable_crosswalk);
DECLARE_string(crosswalk_virtual_object_id_prefix);
DECLARE_double(crosswalk_expand_distance);
DECLARE_double(crosswalk_strick_l_distance);
DECLARE_double(crosswalk_loose_l_distance);
/// stop_sign
DECLARE_bool(enable_stop_sign);
DECLARE_string(stop_sign_virtual_object_id_prefix);
DECLARE_double(stop_duration_for_stop_sign);

DECLARE_bool(enable_sqp_solver);

/// thread pool
DECLARE_int32(num_thread_planning_thread_pool);
DECLARE_bool(enable_multi_thread_in_dp_poly_path);
DECLARE_bool(enable_multi_thread_in_dp_st_graph);

#endif  // MODULES_PLANNING_COMMON_PLANNING_GFLAGS_H
