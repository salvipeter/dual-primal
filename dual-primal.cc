#include <algorithm>
#include <fstream>
#include <cmath>
#include <map>
#include <set>

#include "dual-primal.hh"

using namespace Geometry;

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
  dual.faces.resize(primal.verts.size());
  adjacent.resize(primal.faces.size());

  std::map<std::pair<size_t, size_t>, size_t> edgeToFace;
  auto addEdge = [&](size_t i, size_t j, size_t face) {
    if (i > j)
      std::swap(i, j);
    if (edgeToFace.contains({i, j})) {
      auto other = edgeToFace[{i, j}];
      adjacent[face].push_back(other);
      adjacent[other].push_back(face);
    } else
      edgeToFace[{i, j}] = face;
  };

  for (size_t ti = 0; ti < primal.faces.size(); ++ti) {
    const auto &tri = primal.faces[ti];
    Point3D center(0, 0, 0);
    for (size_t i = 0; i < 3; ++i) {
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

void DualPrimal::optimizePrimal() {

}

void DualPrimal::resamplePrimal() {
  size_t n_verts = primal.verts.size();
  for (size_t v = 0; v < n_verts; ++v) {
    double w = 0;
    Point3D p(0, 0, 0);
    for (auto i : dual.faces[v]) {
      auto pi = dual.verts[i];
      auto mi = fdf(pi).second.normalize();
      double ki = 0.0;
      for (auto j : adjacent[i]) {
        auto pj = dual.verts[j];
        auto mj = fdf(pj).second.normalize();
        ki += std::acos(std::clamp(mi * mj, -1.0, 1.0)) / (pi - pj).norm();
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
  std::map<size_t, size_t> halve;
  std::function<void(size_t)> sub = [&](size_t i) {
    subdivide.insert(i);
    for (auto adj : adjacent[i])
      if (!subdivide.contains(adj)) {
        if (halve.contains(adj)) {
          halve.erase(adj);
          sub(adj);
        } else
          halve[adj] = i;
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
    VectorVector normals = {
      ((ab - a) ^ (ac - a)).normalize(),
      ((ab - b) ^ (bc - b)).normalize(),
      ((ac - c) ^ (bc - c)).normalize(),
      ((ac - ab) ^ (bc - ab)).normalize()
    };
    auto area = ((b - a) ^ (c - a)).norm() / 2;
    double err = 0;
    for (size_t i = 0; i < 4; ++i)
      err += 1 - std::abs(normals[i] * fdf(centroids[i]).second.normalize());
    err *= area;
    if (err > eps)
      sub(ti);
  }

  // Subdivide
  std::map<std::pair<size_t, size_t>, size_t> midpoints;
  for (auto ti : subdivide) {
    auto &tri = primal.faces[ti];
    auto mi = primal.verts.size();
    for (size_t i = 0; i < 3; ++i) {
      auto a = tri[i], b = tri[(i+1)%3];
      if (a > b)
        std::swap(a, b);
      midpoints[{a,b}] = mi + i;
      primal.verts.push_back((primal.verts[a] + primal.verts[b]) / 2);
    }
    primal.faces.push_back({ tri[0], mi, mi + 2 });
    primal.faces.push_back({ mi, tri[1], mi + 1 });
    primal.faces.push_back({ mi + 2, mi + 1, tri[2] });
    tri[0] = mi;
    tri[1] = mi + 1;
    tri[2] = mi + 2;
  }

  // Halve
  for (auto [ti1, ti2] : halve) {
    auto &tri = primal.faces[ti1];
    const auto &other = primal.faces[ti2];
    size_t ci = -1;
    for (size_t i = 0; i < 3; ++i)
      if (std::find(other.begin(), other.end(), tri[i]) == other.end()) {
        ci = i;
        break;
      }
    auto a = tri[(ci+1)%3], b = tri[(ci+2)%3], c = tri[ci];
    size_t m;
    if (a > b)
      m = midpoints[{b,a}];
    else
      m = midpoints[{a,b}];
    tri[(ci+1)%3] = m;
    primal.faces.push_back({ m, c, a });
  }
}
