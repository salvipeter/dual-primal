#include <algorithm>
#include <fstream>
#include <cmath>
#include <map>
#include <set>

#include <Eigen/Dense>

#include "dual-primal.hh"

using namespace Eigen;

using Real = DualPrimal::Real;
using Point3D = DualPrimal::Point3D;
using Vector3D = DualPrimal::Vector3D;
using IndexPair = std::pair<size_t, size_t>;
using DoubleVector = std::vector<double>;
using VectorVector = std::vector<Vector3D>;

void DualPrimal::Mesh::writeOBJ(std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &p : verts)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (const auto &face : faces) {
    f << 'f';
    for (auto i : face)
      f << ' ' << i + 1;
    f << std::endl;
  }
}

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
  Point3D q = p, r;
  size_t iter = 0;
  while (true) {                // v, d are the value/gradient at q
    if (++iter % iter_halve == 0) {
      if (++iter == iter_halve * iter_trials)
        return p;
      // Start over with smaller lambda
      auto vd = fdf(p);
      v = vd.first; d = vd.second;
      q = p;
      lambda /= 2;
    }
    d *= -v;
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
  dual.faces.resize(primal.verts.size());
  ff_map.clear();
  ff_map.resize(primal.faces.size());
  vf_map.clear();
  vf_map.resize(primal.verts.size());

  std::map<IndexPair, size_t> edgeToFace;
  auto addEdge = [&](size_t i, size_t j, size_t face) {
    if (i > j)
      std::swap(i, j);
    if (edgeToFace.contains({i, j})) {
      auto other = edgeToFace[{i, j}];
      ff_map[face].push_back(other);
      ff_map[other].push_back(face);
    } else
      edgeToFace[{i, j}] = face;
  };

  for (size_t ti = 0; ti < primal.faces.size(); ++ti) {
    const auto &tri = primal.faces[ti];
    Point3D center(0, 0, 0);
    for (size_t i = 0; i < 3; ++i) {
      vf_map[tri[i]].push_back(ti);
      center += primal.verts[tri[i]];
      dual.faces[tri[i]].push_back(ti);
      addEdge(tri[i], tri[(i+1)%3], ti);
    }
    dual.verts.push_back(center / 3);
  }
}

void DualPrimal::optimizeDual() {
  for (size_t ti = 0; ti < primal.faces.size(); ++ti) {
    const auto &tri = primal.faces[ti];
    double e = 0;
    for (size_t i = 0; i < 3; ++i)
      e += (primal.verts[tri[(i+1)%3]] - primal.verts[tri[i]]).norm();
    e /= 3;
    dual.verts[ti] = project(dual.verts[ti], e);
  }
}

static Vector3D solve(const VectorVector &Lhs, const DoubleVector &Rhs, double tau) {
  Matrix<Real, Dynamic, Dynamic> A(Lhs.size(), 3);
  Vector<Real, Dynamic> b(Rhs.size());
  for (size_t i = 0; i < Lhs.size(); ++i) {
    for (size_t j = 0; j < 3; ++j)
      A(i, j) = Lhs[i][j];
    b(i) = Rhs[i];
  }

  JacobiSVD svd(A, ComputeThinU | ComputeThinV);
  if (tau > 0)
    svd.setThreshold(1 / tau);
  return svd.solve(b);
}

void DualPrimal::optimizePrimal() {
  for (size_t i = 0; i < primal.verts.size(); ++i) {
    auto &x = primal.verts[i];
    VectorVector Lhs(3, Vector3D(0, 0, 0));
    DoubleVector Rhs(3, 0);
    for (auto j : vf_map[i]) {
      const auto &Pj = dual.verts[j];
      auto m = fdf(Pj).second.normalized();
      for (size_t k = 0; k < 3; ++k) {
        Lhs[k] += m * m[k];
        Rhs[k] += m.dot(Pj - x) * m[k];
      }
    }
    x += solve(Lhs, Rhs, tau);
  }
}

void DualPrimal::resamplePrimal() {
  size_t n_verts = primal.verts.size();
  for (size_t v = 0; v < n_verts; ++v) {
    double w = 0;
    Point3D p(0, 0, 0);
    for (auto i : dual.faces[v]) {
      auto pi = dual.verts[i];
      auto mi = fdf(pi).second.normalized();
      double ki = 0.0;
      for (auto j : ff_map[i]) {
        auto pj = dual.verts[j];
        auto mj = fdf(pj).second.normalized();
        ki += std::acos(std::clamp<Real>(mi.dot(mj), -1, 1)) / (pi - pj).norm();
      }
      auto wi = 1 + c * ki;
      p += pi * wi;
      w += wi;
    }
    p /= w;
    primal.verts[v] = p;
  }
}

void DualPrimal::subdividePrimal() {
  std::set<size_t> subdivide;
  std::map<size_t, std::pair<IndexPair, size_t>> halve;
  std::function<void(size_t)> sub = [&](size_t i) {
    halve.erase(i);
    subdivide.insert(i);
    for (auto adj : ff_map[i])
      if (!subdivide.contains(adj)) {
        if (halve.contains(adj))
          sub(adj);
        else {
          const auto &tri = primal.faces[adj];
          const auto &other = primal.faces[i];
          size_t ci = -1;
          for (size_t j = 0; j < 3; ++j)
            if (std::find(other.begin(), other.end(), tri[j]) == other.end()) {
              ci = j;
              break;
            }
          auto a = tri[(ci+1)%3], b = tri[(ci+2)%3];
          if (a > b)
            std::swap(a, b);
          halve[adj] = std::make_pair(std::make_pair(a, b), ci);
        }
      }
  };

  // Check for subdivision
  for (size_t ti = 0; ti < primal.faces.size(); ++ti) {
    const auto &tri = primal.faces[ti];
    const auto &a = primal.verts[tri[0]], &b = primal.verts[tri[1]], &c = primal.verts[tri[2]];
    auto ab = (a + b) / 2, ac = (a + c) / 2, bc = (b + c) / 2;
    PointVector centroids = {
      (a + ab + ac) / 3, (b + ab + bc) / 3, (c + ac + bc) / 3, (ab + ac + bc) / 3
    };
    auto normal = (b - a).cross(c - a);
    auto area = normal.norm() / 2;
    normal.normalize();
    double err = 0;
    for (size_t i = 0; i < 4; ++i)
      err += 1 - std::abs(normal.dot(fdf(centroids[i]).second.normalized()));
    err *= area;
    if (err > eps)
      sub(ti);
  }

  // Subdivide
  std::map<IndexPair, size_t> midpoints;
  for (auto ti : subdivide) {
    auto &tri = primal.faces[ti];
    auto t = tri;
    for (size_t i = 0; i < 3; ++i) {
      auto a = t[i], b = t[(i+1)%3];
      if (a > b)
        std::swap(a, b);
      if (!midpoints.contains({a,b})) {
        midpoints[{a,b}] = primal.verts.size();
        primal.verts.push_back((primal.verts[a] + primal.verts[b]) / 2);
      }
      tri[i] = midpoints[{a,b}];
    }
    auto t2 = tri;
    primal.faces.push_back({ t[0], t2[0], t2[2] });
    primal.faces.push_back({ t2[0], t[1], t2[1] });
    primal.faces.push_back({ t2[2], t2[1], t[2] });
  }

  // Halve
  for (const auto &[ti, abc] : halve) {
    auto &tri = primal.faces[ti];
    auto [ab, ci] = abc;
    auto m = midpoints[ab];
    auto a = tri[(ci+1)%3], c = tri[ci];
    tri[(ci+1)%3] = m;
    primal.faces.push_back({ m, c, a });
  }
}
