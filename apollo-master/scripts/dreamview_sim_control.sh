#!/usr/bin/env bash

###############################################################################
# Copyright 2017 The Apollo Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###############################################################################


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "${DIR}/.."

source "${DIR}/apollo_base.sh"

# run function from apollo_base.sh
# run command_name module_name

function start() {
    echo "Start roscore..."
    ROSCORELOG="${APOLLO_ROOT_DIR}/data/log/roscore.out"
    nohup roscore </dev/null >"${ROSCORELOG}" 2>&1 &
    if [ "$HOSTNAME" == "in_release_docker" ]; then
        supervisord -c /apollo/modules/tools/supervisord/release.conf >& /tmp/supervisord.start.log
        echo "Started supervisord with release conf"
    else
        supervisord -c /apollo/modules/tools/supervisord/dev.conf >& /tmp/supervisord.start.log
        echo "Started supervisord with dev conf"
    fi
    supervisorctl start monitor > /dev/null
    supervisorctl start sim_control > /dev/null
    echo "Dreamview is running at http://localhost:8888"
}

function stop() {
    supervisorctl stop sim_control > /dev/null 2>&1 &
    supervisorctl stop monitor > /dev/null 2>&1 &
    pkill -f roscore
}

case $1 in
  start)
    start
    ;;
  stop)
    stop
    ;;
  *)
    start
    ;;
esac
