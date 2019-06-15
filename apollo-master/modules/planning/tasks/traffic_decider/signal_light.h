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

#ifndef MODULES_PLANNING_TASKS_TRAFFIC_DECIDER_SIGNAL_LIGHT_H_
#define MODULES_PLANNING_TASKS_TRAFFIC_DECIDER_SIGNAL_LIGHT_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "modules/perception/proto/traffic_light_detection.pb.h"

#include "modules/planning/tasks/traffic_decider/traffic_rule.h"

namespace apollo {
namespace planning {

class SignalLight : public TrafficRule {
 public:
  explicit SignalLight(const RuleConfig& config);

  virtual ~SignalLight() = default;

  bool ApplyRule(Frame* frame, ReferenceLineInfo* const reference_line_info);

 private:
  void ReadSignals();
  bool FindValidSignalLight(ReferenceLineInfo* const reference_line_info);
  apollo::perception::TrafficLight GetSignal(const std::string& signal_id);
  void MakeDecisions(Frame* frame,
                     ReferenceLineInfo* const reference_line_info);
  double GetStopDeceleration(ReferenceLineInfo* const reference_line_info,
                             const hdmap::PathOverlap* signal_light);
  bool CreateStopObstacle(Frame* frame,
                          ReferenceLineInfo* const reference_line_info,
                          const hdmap::PathOverlap* signal_light);
  void SetCreepForwardSignalDecision(
      const ReferenceLineInfo* reference_line_info,
      hdmap::PathOverlap* const signal_light) const;

  std::vector<hdmap::PathOverlap> signal_lights_from_path_;
  std::unordered_map<std::string, const apollo::perception::TrafficLight*>
      detected_signals_;
};

}  // namespace planning
}  // namespace apollo

#endif  // MODULES_PLANNING_TASKS_TRAFFIC_DECIDER_SIGNAL_LIGHT_H_
