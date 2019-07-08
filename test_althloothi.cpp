#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/detail/tolerance_manip.hpp>
// #include <iostream>
// #include <iomanip> // for std::setprecision
#include <cmath>
#include "althloothi.hpp"

const double pi{std::acos(-1)};

using quaternion = std::array<double, 4>;

quaternion operator- (const quaternion& rot) {
    return quaternion{-rot[0], -rot[1], -rot[2], -rot[3]};
}

// template <typename T, std::size_t N>
// std::ostream& operator<< (std::ostream &ios, const std::array<T, N> &arr) {
//     ios << "{";
//     for (int i = 0; i < N - 1; ++i)
//         ios << arr[i] << ", ";
//     ios << arr[N-1] <<  "}";
//
//     return ios;
// }
//
// template <typename T>
// std::ostream& operator<< (std::ostream &ios, const std::vector<T> &vec) {
//     int n = vec.size();
//     ios << "{";
//     for (int i = 0; i < n - 1; ++i)
//         ios << vec[i] << ", ";
//     ios << vec[n-1] <<  "}";
//
//     return ios;
// }

#define CHECK_CLOSE_COLLECTION(a_vec, b_vec, tolerance) { \
    assert(a_vec.size() == b_vec.size()); \
    for (std::size_t i = 0; i < a_vec.size(); ++i) { \
        BOOST_TEST(a_vec[i] == b_vec[i], tolerance); \
    } \
}

BOOST_AUTO_TEST_CASE(test_align)
{

    {
        // spheres
        std::array<std::vector<double>, 3> Sh1 = {{{0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 1, 0}}};
        std::array<std::vector<double>, 3> Sh0 = {{{0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 1, 0}}};
        auto qvec = align<quaternion>(Sh1, Sh0);

        // no rotation
        CHECK_CLOSE_COLLECTION(qvec, quaternion({1, 0, 0, 0}), boost::test_tools::tolerance(1e-12));
    }

    {
        // ellipsoids
        std::array<std::vector<double>, 3> Sh1 = {{{0, 0, 0, 2}, {0, 1, 0, 0}, {0, 0, 1, 0}}};
        std::array<std::vector<double>, 3> Sh0 = {{{0, 0, 0, 1}, {0, 2, 0, 0}, {0, 0, 1, 0}}};
        auto qvec = align<quaternion>(Sh1, Sh0);

        // rotation around z-axis
        auto rot = quaternion{std::cos(pi/4), 0, 0, std::sin(pi/4)};
        if (qvec[0] * rot[0] < 0) rot = -rot;
        CHECK_CLOSE_COLLECTION(qvec, rot, boost::test_tools::tolerance(1e-12));
    }

    {
        // ellipsoids
        std::array<std::vector<double>, 3> Sh1 = {{{0, 0, 0, 2}, {0, 1, 0, 0}, {0, 0, 1, 0}}};
        std::array<std::vector<double>, 3> Sh0 = {{{0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 2, 0}}};
        auto qvec = align<quaternion>(Sh1, Sh0);

        // rotation around y-axis
        auto rot = quaternion({std::cos(-pi/4), 0, std::sin(-pi/4), 0});
        if (qvec[0] * rot[0] < 0) rot = -rot;
        CHECK_CLOSE_COLLECTION(qvec, rot, boost::test_tools::tolerance(1e-12));
    }

    {
        std::array<std::vector<double>, 3> Sh1 = {{{0, 0, 0, 1}, {0, 2, 0, 0}, {0, 0, 4, 0}}};
        std::array<std::vector<double>, 3> Sh0 = {{{0, 0, 4, 0}, {0, 0, 0, 1}, {0, 2, 0, 0}}};
        auto qvec = align<quaternion>(Sh1, Sh0);

        auto rot = quaternion{std::cos(pi/3), std::sin(pi/3)/std::sqrt(3), std::sin(pi/3)/std::sqrt(3), std::sin(pi/3)/std::sqrt(3)};
        if (qvec[0] * rot[0] < 0) rot = -rot;
        CHECK_CLOSE_COLLECTION(qvec, rot, boost::test_tools::tolerance(1e-12));
    }

    {
        // some wavy shape
        std::array<std::vector<double>, 3> Sh1 = {{{0.,0.,0.,1.,0.5,0.}, {0.,2.,0.,0.,1.,1.}, {0.,0.,3.,0.,1.,0.}}};
        std::array<std::vector<double>, 3> Sh0 = {{{0.,0.,3.,0.,1.,0.}, {0.,0.,0.,1.,0.5,0.}, {0.,2.,0.,0.,1.,1.}}};
        auto qvec = align<quaternion>(Sh1, Sh0);

        auto rot = quaternion{std::cos(pi/3), std::sin(pi/3)/std::sqrt(3), std::sin(pi/3)/std::sqrt(3), std::sin(pi/3)/std::sqrt(3)};
        if (qvec[0] * rot[0] < 0) rot = -rot;
        CHECK_CLOSE_COLLECTION(qvec, rot, boost::test_tools::tolerance(1e-12));
    }
}
