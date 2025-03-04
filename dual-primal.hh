#pragma once

// Based on:
//   Y. Ohtake, A.G. Belyaev:
//     Dual/primal mesh optimization for polygonized implicit surfaces.
//       Proceedings of the 7th ACM symposium on Solid modeling and applications,
//       pp. 171-178, 2002.
//     https://doi.org/10.1145/566282.56630

#include <functional>
#include <geometry.hh>          // https://github.com/salvipeter/libgeom/

struct DualPrimal {
  using ValueGradient = std::pair<double, Geometry::Vector3D>;
  using Surface = std::function<ValueGradient(const Geometry::Point3D &)>;

  Surface fdf;                  // implicit surace (value + gradient)
  Geometry::TriMesh primal;     // the original Marching Cubes mesh

  double eps = 1e-3;            // precision - see "Optimizing dual mesh" on p. 173
  double tau = 1e3;             // singular value threshold - see below Eq.(1) on p. 174
  double c = 2.0;               // curvature dependency - see Eq. (4) on p. 174
  size_t m = 3, n = 3;          // iterations - see the end of p. 175

  void optimize();

private:
  struct DualMesh {
    Geometry::PointVector verts;               // at projected primal mass centers
    std::vector<std::vector<size_t>> faces;    // correspond to primal vertices
    std::vector<std::vector<size_t>> adjacent; // (at most 3) adjacent (dual) vertices
  } dual;
  Geometry::Point3D project(const Geometry::Point3D &result, double edge_length) const;
  void createDual();
  void optimizeDual();
  void optimizePrimal();
  void resamplePrimal();
  void subdividePrimal();
};
