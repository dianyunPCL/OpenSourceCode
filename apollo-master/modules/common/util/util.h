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
 * @brief Some util functions.
 */

#ifndef MODULES_COMMON_UTIL_UTIL_H_
#define MODULES_COMMON_UTIL_UTIL_H_

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "google/protobuf/util/message_differencer.h"

#include "modules/common/math/vec2d.h"
#include "modules/common/proto/geometry.pb.h"
#include "modules/common/proto/pnc_point.pb.h"
#include "modules/perception/proto/perception_obstacle.pb.h"

/**
 * @namespace apollo::common::util
 * @brief apollo::common::util
 */
namespace apollo {
namespace common {
namespace util {

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename ProtoA, typename ProtoB>
bool IsProtoEqual(const ProtoA& a, const ProtoB& b) {
  return google::protobuf::util::MessageDifferencer::Equals(a, b);
}

struct PairHash {
  template <typename T, typename U>
  size_t operator()(const std::pair<T, U>& pair) const {
    return std::hash<T>()(pair.first) ^ std::hash<U>()(pair.second);
  }
};

template <typename T>
bool WithinBound(T start, T end, T value) {
  return value >= start && value <= end;
}

/**
 * @brief create a SL point
 * @param s the s value
 * @param l the l value
 * @return a SLPoint instance
 */
SLPoint MakeSLPoint(const double s, const double l);

template <typename T>
common::math::Vec2d MakeVec2d(const T& t) {
  return common::math::Vec2d(t.x(), t.y());
}

PointENU MakePointENU(const double x, const double y, const double z);

PointENU MakePointENU(const math::Vec2d& xy);

apollo::perception::Point MakePerceptionPoint(const double x, const double y,
                                              const double z);

SpeedPoint MakeSpeedPoint(const double s, const double t, const double v,
                          const double a, const double da);

PathPoint MakePathPoint(const double x, const double y, const double z,
                        const double theta, const double kappa,
                        const double dkappa, const double ddkappa);

/**
 * uniformly slice a segment [start, end] to num + 1 pieces
 * the result sliced will contain the n + 1 points that slices the provided
 * segment. `start` and `end` will be the first and last element in `sliced`.
 */
void uniform_slice(double start, double end, uint32_t num,
                   std::vector<double>* sliced);

template <typename Container>
typename Container::value_type MaxElement(const Container& elements) {
  return *std::max_element(elements.begin(), elements.end());
}

template <typename Container>
typename Container::value_type MinElement(const Container& elements) {
  return *std::min_element(elements.begin(), elements.end());
}

template <typename T>
std::unordered_set<T> Intersection(const std::unordered_set<T>& s1,
                                   const std::unordered_set<T>& s2) {
  if (s1.size() < s2.size()) {
    std::unordered_set<T> result;
    for (const auto& v : s1) {
      if (s2.count(v) > 0) {
        result.insert(v);
      }
    }
    return result;
  } else {
    return intersection(s2, s1);
  }
}

/**
 * calculate the distance beteween Point u and Point v, which are all have
 * member function x() and y() in XY dimension.
 * @param u one point that has member function x() and y().
 * @param b one point that has member function x() and y().
 * @return sqrt((u.x-v.x)^2 + (u.y-v.y)^2), i.e., the Euclid distance on XY
 * dimension.
 */
template <typename U, typename V>
double DistanceXY(const U& u, const V& v) {
  return std::hypot(u.x() - v.x(), u.y() - v.y());
}

/**
 * Check if two points u and v are the same point on XY dimension.
 * @param u one point that has member function x() and y().
 * @param v one point that has member function x() and y().
 * @return sqrt((u.x-v.x)^2 + (u.y-v.y)^2) < epsilon, i.e., the Euclid distance
 * on XY dimension.
 */
template <typename U, typename V>
bool SamePointXY(const U& u, const V& v) {
  constexpr double kMathEpsilonSqr = 1e-8 * 1e-8;
  return (u.x() - v.x()) * (u.x() - v.x()) < kMathEpsilonSqr &&
         (u.y() - v.y()) * (u.y() - v.y()) < kMathEpsilonSqr;
}

}  // namespace util
}  // namespace common
}  // namespace apollo

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, std::pair<A, B>& p) {
  return os << "first: " << p.first << ", second: " << p.second;
}

#endif  // MODULES_COMMON_UTIL_UTIL_H_
