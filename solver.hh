#pragma once

#include <geometry.hh>

namespace QEFSolver {

// Finds minimum norm solution of Lhs x = Rhs by singular value decomposition.
// All singular values s_i less than s_0/tau are regarded as zeroes.
// The default setting of tau = 0 uses Eigen's default behavior.
Geometry::Vector3D solve(const Geometry::VectorVector &Lhs,
                         const Geometry::DoubleVector &Rhs,
                         double tau = 0);

}
