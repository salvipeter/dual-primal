#include <algorithm>
#include <cmath>
#include <map>

#include "dual-primal.hh"

using namespace Geometry;

void DualPrimal::optimize() {
  for (size_t i = 0; i <= n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      createDual();
      optimizeDual();
      resamplePrimal();
    }
    createDual();
    optimizeDual();
    optimizePrimal();
    if (i < n)
      subdividePrimal();
  }
}

Point3D DualPrimal::project(const Point3D &p, double edge_length) const {
  auto [v, d] = fdf(p);
  if (std::abs(v) / d.norm() < eps)
    return p;
  auto lambda = edge_length / 2;
  // TODO: "If the number of iterations is too large, we set lambda to half its previous value
  //        in order to catch thin components of the implicit surface." ???
  Point3D q = p, r;
  while (true) {                // v, d are the value/gradient at q
    d *= -1;
    d.normalize();
    r = q + d * lambda;
    auto [v1, d1] = fdf(r);
    if (v * v1 < 0)
      break;
    if ((r - p).norm() > edge_length)
      return p;
    q = r;
    v = v1;
    d = d1;
  }
  // Bisection between q and r
  double lo = 0, hi = 1;
  Point3D result;
  while (true) {
    double mid = (lo + hi) / 2;
    result = q + (r - q) * mid;
    auto [v1, d1] = fdf(result);
    if (std::abs(v1) / d1.norm() < eps)
      break;
    if (v * v1 < 0)
      hi = mid;
    else
      lo = mid;
  }
  return result;
}

void DualPrimal::createDual() {
  dual.verts.clear();
  dual.faces.clear();
  dual.faces.resize(primal.points().size());
  dual.adjacent.resize(primal.triangles().size());

  std::map<std::pair<size_t, size_t>, size_t> edgeToFace;
  auto addEdge = [&](size_t i, size_t j, size_t face) {
    if (i > j)
      std::swap(i, j);
    if (edgeToFace.contains({i, j})) {
      auto other = edgeToFace[{i, j}];
      dual.adjacent[face].push_back(other);
      dual.adjacent[other].push_back(face);
    } else
      edgeToFace[{i, j}] = face;
  };

  size_t count = 0;
  for (const auto &tri : primal.triangles()) {
    Point3D center(0, 0, 0);
    for (size_t i = 0; i < 3; ++i) {
      center += primal[tri[i]];
      dual.faces[tri[i]].push_back(count);
      addEdge(tri[i], tri[(i+1)%3], count);
    }
    dual.verts.push_back(center / 3);
    count++;
  }
}

void DualPrimal::optimizeDual() {
  size_t count = 0;
  for (const auto &tri : primal.triangles()) {
    double e = 0;
    for (size_t i = 0; i < 3; ++i)
      e += (primal[tri[(i+1)%3]] - primal[tri[i]]).norm();
    e /= 3;
    dual.verts[count] = project(dual.verts[count], e);
    count++;
  }
}

void DualPrimal::optimizePrimal() {

}

void DualPrimal::resamplePrimal() {
  size_t n_verts = primal.points().size();
  for (size_t v = 0; v < n_verts; ++v) {
    double w = 0;
    Point3D p(0, 0, 0);
    for (auto i : dual.faces[v]) {
      auto pi = dual.verts[i];
      auto mi = fdf(pi).second.normalize();
      double ki = 0.0;
      for (auto j : dual.adjacent[i]) {
        auto pj = dual.verts[j];
        auto mj = fdf(pj).second.normalize();
        ki += std::acos(std::clamp(mi * mj, -1.0, 1.0)) / (pi - pj).norm();
      }
      auto wi = 1 + c * ki;
      p += pi * wi;
      w += wi;
    }
    p /= w;
    primal[v] = p;
  }
}

void DualPrimal::subdividePrimal() {

}
