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
  struct Mesh {
    Geometry::PointVector verts;
    std::vector<std::vector<size_t>> faces;
    void writeOBJ(std::string filename);
  };

  Surface fdf;                  // implicit surace (value + gradient)
  Mesh primal;                  // the original Marching Cubes mesh

  double eps = 1e-4;            // precision - see "Optimizing dual mesh" on p. 173
  size_t iter_halve = 20;       // after how many iterations should lambda be halved, cf. p. 173
  size_t iter_trials = 5;       // maximum number of projection trials (not in the paper)
  double tau = 1e3;             // singular value threshold - see below Eq.(1) on p. 174
  double c = 2.0;               // curvature dependency - see Eq. (4) on p. 174
  size_t m = 3, n = 3;          // iterations - see the end of p. 175

  void optimize();

private:
  Mesh dual;
  std::vector<std::vector<size_t>> ff_map; // face-face (dual vertex-vertex) adjacency
  std::vector<std::vector<size_t>> vf_map; // vertex-face (dual face-vertex) adjacency

  Geometry::Point3D project(const Geometry::Point3D &result, double edge_length) const;
  void createDual();
  void optimizeDual();
  void optimizePrimal();
  void resamplePrimal();
  void subdividePrimal();
};
