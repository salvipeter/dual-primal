#include <algorithm>
#include <cmath>

#include "dual-primal.hh"
#include "marching.hh"          // https://github.com/salvipeter/marching/
#include "mc.h"                 // https://github.com/agostonsipos/implici_marching_cubes/

using namespace Geometry;

auto star(const Point3D &p) {
  auto x = p[0] * p[0], y = p[1] * p[1], z = p[2] * p[2];
  auto v = 400 * (x * y + y * z + x * z) - std::pow(1 - x - y - z, 3);
  auto dx = 800 * p[0] * (y + z) + 6 * p[0] * std::pow(1 - x - y - z, 2);
  auto dy = 800 * p[1] * (x + z) + 6 * p[1] * std::pow(1 - x - y - z, 2);
  auto dz = 800 * p[2] * (x + y) + 6 * p[2] * std::pow(1 - x - y - z, 2);
  return std::make_pair(v, Vector3D(dx, dy, dz));
}

auto hunt(const Point3D &p) {
  auto x = p[0] * p[0], y = p[1] * p[1], z = p[2] * p[2];
  auto v = 4 * std::pow(x + y + z - 13, 3) + 27 * std::pow(3 * x + y - 4 * z - 12, 2);
  auto dx = 24 * p[0] * std::pow(x + y + z - 13, 2) + 324 * p[0] * (3 * x + y - 4 * z - 12);
  auto dy = 24 * p[1] * std::pow(x + y + z - 13, 2) + 108 * p[1] * (3 * x + y - 4 * z - 12);
  auto dz = 24 * p[2] * std::pow(x + y + z - 13, 2) - 432 * p[2] * (3 * x + y - 4 * z - 12);
  return std::make_pair(v, Vector3D(dx, dy, dz));
}

auto spheres(const Point3D &p) {
  auto x = p[0] * p[0], y = p[1] * p[1], z = p[2] * p[2];
  auto v = x + y + z - 4.0/3.0 * std::abs(p[0]) - 5.0/9.0;
  auto dx = 2 * p[0] - 4.0/3.0 * (x > 0 ? 1 : -1);
  auto dy = 2 * p[1];
  auto dz = 2 * p[2];
  return std::make_pair(v, Vector3D(dx, dy, dz));
}

double chebyshev(size_t n, double x) {
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  return 2 * x * chebyshev(n - 1, x) - chebyshev(n - 2, x);
}

double chebyshevd(size_t n, double x) {
  if (n < 2)
    return n;
  return 2 * chebyshev(n - 1, x) + 2 * x * chebyshevd(n - 1, x) - chebyshevd(n - 2, x);
}

// Banchoff-Chmutov, n is even
auto chmutov(size_t n) {
  return [=](const Point3D &p) {
    auto v = chebyshev(n, p[0]) + chebyshev(n, p[1]) + chebyshev(n, p[2]);
    auto dx = chebyshevd(n, p[0]);
    auto dy = chebyshevd(n, p[1]);
    auto dz = chebyshevd(n, p[2]);
    return std::make_pair(v, Vector3D(dx, dy, dz));
  };
}

int main() {
  DualPrimal dp;
  dp.fdf = spheres;
  // dp.fdf = hunt;
  // dp.fdf = chmutov(4);
  auto mesh = MarchingCubes::isosurface([&](const Point3D &p) { return dp.fdf(p).first; },
                                        {0, 0, 0}, 4.5, 4, 4);
  // auto mesh = IMC::marching_cubes([&](const Point3D &p) { return dp.fdf(p).first; }, 0,
  //                                 { { { -4, -4, -4 }, { 4, 4, 4 } } },
  //                                 { { 20, 20, 20 } });
  dp.primal.verts = mesh.points();
  std::transform(mesh.triangles().begin(), mesh.triangles().end(),
                 std::back_inserter(dp.primal.faces),
                 [](const TriMesh::Triangle &t) {
                   return std::vector<size_t>(t.begin(), t.end());
                 });
  dp.optimize();
  dp.primal.writeOBJ("/tmp/dp.obj");
}
