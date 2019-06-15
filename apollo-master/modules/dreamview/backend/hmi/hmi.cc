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

#include "modules/dreamview/backend/hmi/hmi.h"

#include <cstdlib>
#include <thread>
#include <vector>

#include "gflags/gflags.h"
#include "modules/common/adapters/adapter_manager.h"
#include "modules/common/kv_db/kv_db.h"
#include "modules/common/util/http_client.h"
#include "modules/common/util/json_util.h"
#include "modules/common/util/map_util.h"
#include "modules/common/util/string_tokenizer.h"
#include "modules/common/util/string_util.h"
#include "modules/common/util/util.h"
#include "modules/control/proto/pad_msg.pb.h"
#include "modules/data/proto/static_info.pb.h"
#include "modules/dreamview/backend/common/dreamview_gflags.h"
#include "modules/dreamview/backend/hmi/vehicle_manager.h"
#include "modules/monitor/proto/system_status.pb.h"

DEFINE_string(global_flagfile, "modules/common/data/global_flagfile.txt",
              "Global flagfile shared by all modules.");

DEFINE_string(map_data_path, "/apollo/modules/map/data", "Path to map data.");

DEFINE_string(vehicle_data_path, "/apollo/modules/calibration/data",
              "Path to vehicle data.");

DEFINE_string(ota_service_url, "http://180.76.145.202:5000/query",
              "OTA service url. [Attention! It's still in experiment.]");
DEFINE_string(ota_vehicle_info_file, "modules/tools/ota/vehicle_info.pb.txt",
              "Vehicle info to request OTA.");

namespace apollo {
namespace dreamview {
namespace {

using apollo::canbus::Chassis;
using apollo::common::adapter::AdapterManager;
using apollo::common::util::FindOrNull;
using apollo::common::util::GetProtoFromASCIIFile;
using apollo::common::util::JsonUtil;
using apollo::common::util::StringTokenizer;
using apollo::control::DrivingAction;
using apollo::data::VehicleInfo;
using google::protobuf::Map;
using Json = WebSocketHandler::Json;

// Convert a string to be title-like. E.g.: "hello_world" -> "Hello World".
std::string TitleCase(const std::string &origin,
                      const std::string &delimiter = "_") {
  std::vector<std::string> parts = StringTokenizer::Split(origin, delimiter);
  for (auto &part : parts) {
    if (!part.empty()) {
      // Upper case the first char.
      part[0] = toupper(part[0]);
    }
  }

  return apollo::common::util::PrintIter(parts);
}

// List subdirs and return a dict of {subdir_title: subdir_path}.
Map<std::string, std::string> ListDirAsDict(const std::string &dir) {
  Map<std::string, std::string> result;
  const auto subdirs = apollo::common::util::ListSubDirectories(dir);
  for (const auto &subdir : subdirs) {
    const auto subdir_title = TitleCase(subdir);
    const auto subdir_path = apollo::common::util::StrCat(dir, "/", subdir);
    result.insert({subdir_title, subdir_path});
  }
  return result;
}

// Send PadMessage to change driving mode to target mode.
// Retry for several times to try to guarantee the result.
bool GuaranteeDrivingMode(const Chassis::DrivingMode target_mode,
                          const bool reset_first) {
  if (reset_first) {
    if (!GuaranteeDrivingMode(Chassis::COMPLETE_MANUAL, false)) {
      return false;
    }
  }

  control::PadMessage pad;
  switch (target_mode) {
    case Chassis::COMPLETE_MANUAL:
      pad.set_action(DrivingAction::RESET);
      break;
    case Chassis::COMPLETE_AUTO_DRIVE:
      pad.set_action(DrivingAction::START);
      break;
    default:
      AFATAL << "Unknown action to change driving mode to " << target_mode;
  }

  constexpr int kMaxTries = 3;
  constexpr auto kTryInterval = std::chrono::milliseconds(500);
  auto* chassis = CHECK_NOTNULL(AdapterManager::GetChassis());
  for (int i = 0; i < kMaxTries; ++i) {
    // Send driving action periodically until entering target driving mode.
    AdapterManager::FillPadHeader("HMI", &pad);
    AdapterManager::PublishPad(pad);

    std::this_thread::sleep_for(kTryInterval);
    chassis->Observe();
    if (chassis->Empty()) {
      AERROR << "No Chassis message received!";
    } else if (chassis->GetLatestObserved().driving_mode() == target_mode) {
      return true;
    }
  }
  AERROR << "Failed to change driving mode to " << target_mode;
  return false;
}

}  // namespace

HMI::HMI(WebSocketHandler *websocket, MapService *map_service)
    : websocket_(websocket), map_service_(map_service) {
  CHECK(common::util::GetProtoFromFile(FLAGS_hmi_config_filename, &config_))
      << "Unable to parse HMI config file " << FLAGS_hmi_config_filename;
  // If the module path doesn't exist, remove it from list.
  auto *modules = config_.mutable_modules();
  for (auto iter = modules->begin(); iter != modules->end();) {
    const auto &conf = iter->second;
    if (conf.has_path() && !common::util::PathExists(conf.path())) {
      iter = modules->erase(iter);
    } else {
      ++iter;
    }
  }

  // If the default mode is unavailable, select the first one.
  const auto &modes = config_.modes();
  if (!ContainsKey(modes, status_.current_mode())) {
    CHECK(!modes.empty());
    status_.set_current_mode(modes.begin()->first);
  }

  // Get available maps and vehicles by listing data directory.
  *config_.mutable_available_maps() = ListDirAsDict(FLAGS_map_data_path);
  *config_.mutable_available_vehicles() =
      ListDirAsDict(FLAGS_vehicle_data_path);
  ADEBUG << "Loaded HMI config: " << config_.DebugString();

  // Register websocket message handlers.
  if (websocket_) {
    RegisterMessageHandlers();
  }
}

void HMI::RegisterMessageHandlers() {
  // Send current config and status to new HMI client.
  websocket_->RegisterConnectionReadyHandler(
      [this](WebSocketHandler::Connection *conn) {
        websocket_->SendData(
            conn, JsonUtil::ProtoToTypedJson("HMIConfig", config_).dump());
        websocket_->SendData(
            conn, JsonUtil::ProtoToTypedJson("HMIStatus", status_).dump());
      });

  // HMI client asks for executing module command.
  websocket_->RegisterMessageHandler(
      "ExecuteModuleCommand",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {module: "module_name", command: "command_name"}.
        // If module_name is "all", then run the command on all modules.
        std::string module;
        std::string command;
        if (JsonUtil::GetStringFromJson(json, "module", &module) &&
            JsonUtil::GetStringFromJson(json, "command", &command)) {
          RunComponentCommand(config_.modules(), module, command);
        } else {
          AERROR << "Truncated module command.";
        }
      });

  // HMI client asks for executing tool command.
  websocket_->RegisterMessageHandler(
      "ExecuteToolCommand",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {tool: "tool_name", command: "command_name"}.
        std::string tool;
        std::string command;
        if (JsonUtil::GetStringFromJson(json, "tool", &tool) &&
            JsonUtil::GetStringFromJson(json, "command", &command)) {
          RunComponentCommand(config_.tools(), tool, command);
        } else {
          AERROR << "Truncated tool command.";
        }
      });

  // HMI client asks for executing mode command.
  websocket_->RegisterMessageHandler(
      "ExecuteModeCommand",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {command: "command_name"}.
        // Supported commands are: "start", "stop".
        std::string command;
        if (JsonUtil::GetStringFromJson(json, "command", &command)) {
          RunModeCommand(command);
        } else {
          AERROR << "Truncated mode command.";
        }
      });

  // HMI client asks for changing driving mode.
  websocket_->RegisterMessageHandler(
      "ChangeDrivingMode",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {new_mode: "DrivingModeName"}.
        // DrivingModeName should be one of canbus::Chassis::DrivingMode.
        // For now it is either COMPLETE_MANUAL or COMPLETE_AUTO_DRIVE.
        std::string new_mode;
        if (JsonUtil::GetStringFromJson(json, "new_mode", &new_mode)) {
          ChangeDrivingModeTo(new_mode);
        } else {
          AERROR << "Truncated ChangeDrivingMode request.";
        }
      });

  // HMI client asks for changing map.
  websocket_->RegisterMessageHandler(
      "ChangeMap",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {new_map: "MapName"}.
        // MapName should be a key of config_.available_maps.
        std::string new_map;
        if (JsonUtil::GetStringFromJson(json, "new_map", &new_map)) {
          ChangeMapTo(new_map);
        } else {
          AERROR << "Truncated ChangeMap request.";
        }
      });

  // HMI client asks for changing vehicle.
  websocket_->RegisterMessageHandler(
      "ChangeVehicle",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {new_vehicle: "VehicleName"}.
        // VehicleName should be a key of config_.available_vehicles.
        std::string new_vehicle;
        if (JsonUtil::GetStringFromJson(json, "new_vehicle", &new_vehicle)) {
          ChangeVehicleTo(new_vehicle);
        } else {
          AERROR << "Truncated ChangeVehicle request.";
        }
      });

  // HMI client asks for changing mode.
  websocket_->RegisterMessageHandler(
      "ChangeMode",
      [this](const Json &json, WebSocketHandler::Connection *conn) {
        // json should contain {new_mode: "ModeName"}.
        // ModeName should be a key of config_.modes.
        std::string new_mode;
        if (JsonUtil::GetStringFromJson(json, "new_mode", &new_mode)) {
          ChangeModeTo(new_mode);
        } else {
          AERROR << "Truncated ChangeMode request.";
        }
      });

  // Received new system status, broadcast to clients.
  AdapterManager::AddSystemStatusCallback(
      [this](const monitor::SystemStatus &system_status) {
        *status_.mutable_system_status() = system_status;
        BroadcastHMIStatus();
      });
}

void HMI::BroadcastHMIStatus() const {
  // In unit tests, we may leave websocket_ as NULL and skip broadcasting.
  if (websocket_) {
    websocket_->BroadcastData(
        JsonUtil::ProtoToTypedJson("HMIStatus", status_).dump());
  }
}

int HMI::RunComponentCommand(const Map<std::string, Component> &components,
                             const std::string &component_name,
                             const std::string &command_name) {
  const auto *component = FindOrNull(components, component_name);
  if (component == nullptr) {
    AERROR << "Cannot find component " << component_name;
    return -1;
  }
  const auto *cmd = FindOrNull(component->supported_commands(), command_name);
  if (cmd == nullptr) {
    AERROR << "Cannot find command " << component_name << "." << command_name;
    return -1;
  }
  AINFO << "Execute system command: " << *cmd;
  const int ret = std::system(cmd->c_str());

  AERROR_IF(ret != 0) << "Command returns " << ret << ": " << *cmd;
  return ret;
}

void HMI::RunModeCommand(const std::string &command_name) {
  RunModeCommand(status_.current_mode(), command_name);
}

void HMI::RunModeCommand(const std::string &mode,
                         const std::string &command_name) {
  const Mode &mode_conf = config_.modes().at(mode);
  if (command_name == "start" || command_name == "stop") {
    // Run the command on all live modules.
    for (const auto &module : mode_conf.live_modules()) {
      RunComponentCommand(config_.modules(), module, command_name);
    }
  }
}

void HMI::ChangeDrivingModeTo(const std::string &new_mode) {
  Chassis::DrivingMode mode;
  if (!Chassis::DrivingMode_Parse(new_mode, &mode)) {
    AERROR << "Unknown driving mode " << new_mode;
    return;
  }
  const bool reset_first = (mode != Chassis::COMPLETE_MANUAL);
  GuaranteeDrivingMode(mode, reset_first);
}

void HMI::ChangeMapTo(const std::string &map_name) {
  if (status_.current_map() == map_name) {
    return;
  }
  const auto *map_dir = FindOrNull(config_.available_maps(), map_name);
  if (map_dir == nullptr) {
    AERROR << "Unknown map " << map_name;
    return;
  }
  status_.set_current_map(map_name);
  apollo::common::KVDB::Put("apollo:dreamview:map", map_name);

  FLAGS_map_dir = *map_dir;
  // Append new map_dir flag to global flagfile.
  std::ofstream fout(FLAGS_global_flagfile, std::ios_base::app);
  CHECK(fout) << "Fail to open " << FLAGS_global_flagfile;
  fout << "\n--map_dir=" << *map_dir << std::endl;
  // Also reload simulation map.
  CHECK(map_service_->ReloadMap(true)) << "Failed to load map from "
                                       << *map_dir;
  RunModeCommand("stop");
  BroadcastHMIStatus();
}

void HMI::ChangeVehicleTo(const std::string &vehicle_name) {
  if (status_.current_vehicle() == vehicle_name) {
    return;
  }
  const auto *vehicle = FindOrNull(config_.available_vehicles(), vehicle_name);
  if (vehicle == nullptr) {
    AERROR << "Unknown vehicle " << vehicle_name;
    return;
  }
  status_.set_current_vehicle(vehicle_name);
  apollo::common::KVDB::Put("apollo:dreamview:vehicle", vehicle_name);

  CHECK(VehicleManager::instance()->UseVehicle(*vehicle));
  RunModeCommand("stop");
  // Check available updates for current vehicle.
  // CheckOTAUpdates();
  BroadcastHMIStatus();
}

void HMI::ChangeModeTo(const std::string &mode_name) {
  if (status_.current_mode() == mode_name) {
    return;
  }
  if (!ContainsKey(config_.modes(), mode_name)) {
    AERROR << "Unknown mode " << mode_name;
    return;
  }
  const std::string previous_mode = status_.current_mode();
  status_.set_current_mode(mode_name);
  apollo::common::KVDB::Put("apollo:dreamview:mode", mode_name);

  RunModeCommand(previous_mode, "stop");
  BroadcastHMIStatus();
}

void HMI::CheckOTAUpdates() {
  VehicleInfo vehicle_info;
  if (!GetProtoFromASCIIFile(FLAGS_ota_vehicle_info_file, &vehicle_info)) {
    return;
  }

  Json ota_request;
  ota_request["car_type"] = apollo::common::util::StrCat(
      VehicleInfo::Brand_Name(vehicle_info.brand()),
      ".", VehicleInfo::Model_Name(vehicle_info.model()));
  ota_request["vin"] = vehicle_info.license().vin();
  ota_request["tag"] = std::getenv("DOCKER_IMG");

  Json ota_response;
  const auto status = apollo::common::util::HttpClient::Post(
      FLAGS_ota_service_url, ota_request, &ota_response);
  if (status.ok()) {
    CHECK(JsonUtil::GetStringFromJson(ota_response, "tag",
                                      status_.mutable_ota_update()));
    AINFO << "Found available OTA update: " << status_.ota_update();
  }
}

}  // namespace dreamview
}  // namespace apollo
