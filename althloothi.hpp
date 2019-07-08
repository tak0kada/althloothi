#include <vector>
#include <array>
#include <cassert>
#include <numeric>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
// #include <x86intrin.h>


using quaternion = std::array<double, 4>;
using matrix = std::array<std::vector<double>, 3>;

template <typename matrix>
std::vector<int> dim(const matrix&);

template <>
std::vector<int> dim(const std::array<std::vector<double>, 3> m) {
    return {3, static_cast<int>(m[0].size())};
}


// Althloothi, Salah, Mohammad H. Mahoor, and Richard M. Voyles.
// "A robust method for rotation estimation using spherical harmonics representation."
// IEEE Transactions on Image Processing 22.6 (2013): 2306-2316.
//
// align two objects represented by 3xN matrix, whose i-th column is a vector for sphere harmonics coefficients
template <typename quaternion, typename matrix>
quaternion align(const matrix& from, const matrix& to) {
#ifdef DEBUG
    assert(dim(from)[0] == 3);
    assert(dim(to)[0] == 3);
    assert(dim(from)[1] == dim(to)[1]);
#endif
    const int N{dim(from)[1]};

    // p2311 (18)
    // Covariance matrix of obj and target
    Eigen::MatrixXd Cov(3,3);
    for (int i = 0; i < 3; ++i) {
        // calculate mean of coeffcient of object 'from'
        const double mu_from = std::accumulate(from[i].begin(), from[i].end(), 0);

        for (int j = 0; j < 3; ++j) {
            // calculate mean of coeffcient of object 'to'
            const double mu_to = std::accumulate(to[j].begin(), to[j].end(), 0);

            std::vector<double> c(N);
            for (int n = 0; n < N; ++n) {
                c[n] = (from[i][n] - mu_from) * (to[j][n] - mu_to) / N;
            }
            Cov(i, j) = std::accumulate(c.begin(), c.end(), 0);
        }
    }

    //p.2311 (19)
    Eigen::MatrixXd Q(4, 4);
    Q(0,0) =  Cov(0,0) + Cov(1,1) + Cov(2,2);
    Q(0,1) =  Cov(1,2) - Cov(2,1);
    Q(0,2) =  Cov(2,0) - Cov(0,2);
    Q(0,3) =  Cov(0,1) - Cov(1,0);
    Q(1,0) =  Q(0,1);
    Q(1,1) =  Cov(0,0) - Cov(1,1) - Cov(2,2);
    Q(1,2) =  Cov(0,1) + Cov(1,0);
    Q(1,3) =  Cov(2,0) + Cov(0,2);
    Q(2,0) =  Q(0,2);
    Q(2,1) =  Q(1,2);
    Q(2,2) = -Cov(0,0) + Cov(1,1) - Cov(2,2);
    Q(2,3) =  Cov(1,2) + Cov(2,1);
    Q(3,0) =  Q(0,3);
    Q(3,1) =  Q(1,3);
    Q(3,2) =  Q(2,3);
    Q(3,3) = -Cov(0,0) - Cov(1,1) + Cov(2,2);

    // Eigen decomposition
    Eigen::EigenSolver<Eigen::MatrixXd> solver{Q};
    const auto val = solver.eigenvalues();
    const auto vec = solver.eigenvectors();

    // find the largest eigen value and correspoinding eigen vector
    // quarternion to rotate obj is the eigenvector corresponds to the largest eigenvalue
    std::size_t idx{0};
    // type of val(idx) is std::complex
    // we don't need imaginary part because we know the element is real value
    double max = val(idx).real();
    for (int i = 0; i < 4; ++i) {
        if (val(i).real() > max) {
            max = val(idx).real();
            idx = i;
        }
    }

    quaternion ret{
        // Eigen matrix is column-wise by default
        vec.col(idx)(0).real(), vec.col(idx)(1).real(), vec.col(idx)(2).real(), vec.col(idx)(3).real()};
    return ret;
}
