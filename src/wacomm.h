#ifndef WACOMM_H
#define WACOMM_H


#include <boost/multi_array.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

typedef boost::multi_array_ref<float, 1> float1d;
typedef boost::multi_array_ref<float, 2> float2d;
typedef boost::multi_array_ref<float, 3> float3d;
typedef boost::multi_array_ref<double, 1> double1d;
typedef boost::multi_array_ref<double, 2> double2d;
typedef boost::multi_array_ref<double, 3> double3d;

template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
}

#endif 
