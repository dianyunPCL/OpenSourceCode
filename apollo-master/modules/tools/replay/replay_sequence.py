#!/usr/bin/env python

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
"""
This program can replay a planning output pb file via ros
"""
import os.path
import sys
import argparse
import rospy
from std_msgs.msg import String
from google.protobuf import text_format

import common.proto_utils as proto_utils


def generate_message(filename, pb_type):
    f_handle = file(filename, 'r')
    message = pb_type()
    text_format.Merge(f_handle.read(), message)
    f_handle.close()
    return message


def seq_publisher(seq_num, period):
    """publisher"""
    rospy.init_node('replay_node', anonymous=True)
    messages = {}
    for topic, msg_type in proto_utils.topic_pb_dict.iteritems():
        messages[topic] = {}
        filename = str(seq_num) + "_" + topic + ".pb.txt"
        print "trying to load pb file:", filename
        messages[topic]["publisher"] = rospy.Publisher(
            topic, msg_type, queue_size=1)
        pb_msg = msg_type()
        if not proto_utils.get_pb_from_file(filename, pb_msg):
            print topic, " pb is none"
        messages[topic]["value"] = pb_msg

    rate = rospy.Rate(int(1.0 / period))  # 10hz
    while not rospy.is_shutdown():
        for topic, module_features in messages:
            if module_features["value"] is not None:
                module_features["publisher"].publish(module_features["value"])
        rate.sleep()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="replay a set of pb files with the same sequence number")
    parser.add_argument(
        "seq",
        action="store",
        type=int,
        default=-1,
        help="set sequence number to replay")
    parser.add_argument(
        "--period",
        action="store",
        type=float,
        default=1,
        help="set the topic publish time duration")
    args = parser.parse_args()
    try:
        seq_publisher(args.seq, args.period)

    except rospy.ROSInterruptException:
        print "failed to replay message"
