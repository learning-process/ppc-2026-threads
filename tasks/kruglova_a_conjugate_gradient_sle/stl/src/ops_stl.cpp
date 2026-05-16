#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <thread>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

namespace {

template <typename Func>
void LaunchParallel(int total, int num_threads, Func &&func) {
  if (num_threads <= 1) {
    func(0, total, 0);
    return;
  }
  std::vector<std::thread> workers;
  int chunk = total / num_threads;
  for (int i = 0; i < num_threads; ++i) {
    int start = i * chunk;
    int end = (i == num_threads - 1) ? total : (i + 1) * chunk;
    workers.emplace_back(func, start, end, i);
  }
  for (auto &w : workers) {
    w.join();
  }
}

static void MatVec(int n, int num_threads, const std::vector<double> &a, const std::vector<double> &p,
                   std::vector<double> &ap) {
  LaunchParallel(n, num_threads, [&](int start, int end, int) {
    for (int i = start; i < end; ++i) {
      double sum = 0.0;
      for (int j = 0; j < n; ++j) {
        sum += a[(i * n) + j] * p[j];
      }
      ap[i] = sum;
    }
  });
}

static double Dot(int n, int num_threads, const std::vector<double> &v1, const std::vector<double> &v2,
                  std::vector<double> &buffer) {
  LaunchParallel(n, num_threads, [&](int start, int end, int tid) {
    double local = 0.0;
    for (int i = start; i < end; ++i) {
      local += v1[i] * v2[i];
    }
    buffer[tid] = local;
  });
  return std::accumulate(buffer.begin(), buffer.begin() + num_threads, 0.0);
}

static double UpdateXR(int n, int num_threads, double alpha, const std::vector<double> &p,
                       const std::vector<double> &ap, std::vector<double> &x, std::vector<double> &r,
                       std::vector<double> &buffer) {
  LaunchParallel(n, num_threads, [&](int start, int end, int tid) {
    double local_rs = 0.0;
    for (int i = start; i < end; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * ap[i];
      local_rs += r[i] * r[i];
    }
    buffer[tid] = local_rs;
  });
  return std::accumulate(buffer.begin(), buffer.begin() + num_threads, 0.0);
}

static void UpdateP(int n, int num_threads, double beta, const std::vector<double> &r, std::vector<double> &p) {
  LaunchParallel(n, num_threads, [&](int start, int end, int) {
    for (int i = start; i < end; ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
  });
}

}  // namespace

KruglovaAConjGradSleSTL::KruglovaAConjGradSleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.size > 0 && in.A.size() == static_cast<size_t>(in.size) * in.size &&
         in.b.size() == static_cast<size_t>(in.size);
}

bool KruglovaAConjGradSleSTL::PreProcessingImpl() {
  GetOutput().assign(GetInput().size, 0.0);
  return true;
}

bool KruglovaAConjGradSleSTL::RunImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().b;
  const int n = GetInput().size;
  if (n <= 0) {
    return true;
  }

  int num_threads = (n >= 250) ? static_cast<int>(std::thread::hardware_concurrency()) : 1;
  if (num_threads < 1) {
    num_threads = 1;
  }

  auto &x = GetOutput();
  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);
  std::vector<double> partial_buffer(num_threads);

  double rsold = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
  const double tolerance = 1e-8;

  for (int iter = 0; iter < n; ++iter) {
    MatVec(n, num_threads, a, p, ap);

    double p_ap = Dot(n, num_threads, p, ap, partial_buffer);
    if (std::abs(p_ap) < 1e-16) {
      break;
    }

    const double alpha = rsold / p_ap;
    double rsnew = UpdateXR(n, num_threads, alpha, p, ap, x, r, partial_buffer);

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    UpdateP(n, num_threads, (rsnew / rsold), r, p);
    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
