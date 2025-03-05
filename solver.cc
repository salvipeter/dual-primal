#include <Eigen/Dense>

#include "solver.hh"

namespace QEFSolver {

using namespace Eigen;
using namespace Geometry;

Vector3D solve(const VectorVector &Lhs, const DoubleVector &Rhs, double tau) {
  MatrixXd A(Lhs.size(), 3);
  VectorXd b(Rhs.size());
  for (size_t i = 0; i < Lhs.size(); ++i) {
    for (size_t j = 0; j < 3; ++j)
      A(i, j) = Lhs[i][j];
    b(i) = Rhs[i];
  }

  JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
  if (tau > 0)
    svd.setThreshold(1 / tau);
  return Vector3D(Vector3d(svd.solve(b)).data());

  // Manual version - not needed

  // VectorXd values = svd.singularValues();
  // auto s0 = values(0);
  // size_t m = values.size();
  // for (size_t i = 0; i < m; ++i)
  //   if (values(i) * tau < s0)
  //     values(i) = 0;
  //   else
  //     values(i) = 1.0 / values(i);
  // MatrixXd pseudo_inverse =
  //   svd.matrixV() * values.asDiagonal() * svd.matrixU().transpose();

  // Vector3d x = pseudo_inverse * b;
  // return Vector3D(x.data());
}

}
