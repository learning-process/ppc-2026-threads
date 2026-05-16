#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

namespace kruglova_a_conjugate_gradient_sle {

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
  auto &x = GetOutput();

  if (n <= 0) {
    return true;
  }

  int num_threads = (n >= 250) ? static_cast<int>(std::thread::hardware_concurrency()) : 1;
  if (num_threads < 1) {
    num_threads = 1;
  }

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  std::vector<double> partial_buffer(num_threads);

  auto launch_parallel = [&](int total, auto &&func) {
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
  };

  double rsold = 0.0;
  for (int i = 0; i < n; ++i) {
    rsold += r[i] * r[i];
  }

  const double tolerance = 1e-8;
  const int max_iter = n;

  for (int iter = 0; iter < max_iter; ++iter) {
    // 1. ap = A * p
    if (num_threads > 1) {
      launch_parallel(n, [&](int start, int end, int /*tid*/) {
        for (int i = start; i < end; ++i) {
          double sum = 0.0;
          for (int j = 0; j < n; ++j) {
            sum += a[i * n + j] * p[j];
          }
          ap[i] = sum;
        }
      });
    } else {
      for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
          sum += a[i * n + j] * p[j];
        }
        ap[i] = sum;
      }
    }

    // 2. p_ap = p * ap
    double p_ap = 0.0;
    if (num_threads > 1) {
      launch_parallel(n, [&](int start, int end, int tid) {
        double local = 0.0;
        for (int i = start; i < end; ++i) {
          local += p[i] * ap[i];
        }
        partial_buffer[tid] = local;
      });
      for (double val : partial_buffer) {
        p_ap += val;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        p_ap += p[i] * ap[i];
      }
    }

    if (std::abs(p_ap) < 1e-16) {
      break;
    }
    const double alpha = rsold / p_ap;

    // 3. Обновление x, r и вычисление rsnew
    double rsnew = 0.0;
    if (num_threads > 1) {
      launch_parallel(n, [&](int start, int end, int tid) {
        double local = 0.0;
        for (int i = start; i < end; ++i) {
          x[i] += alpha * p[i];
          r[i] -= alpha * ap[i];
          local += r[i] * r[i];
        }
        partial_buffer[tid] = local;
      });
      for (double val : partial_buffer) {
        rsnew += val;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
        rsnew += r[i] * r[i];
      }
    }

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    // 4. p = r + beta * p
    const double beta = rsnew / rsold;
    if (num_threads > 1) {
      launch_parallel(n, [&](int start, int end, int /*tid*/) {
        for (int i = start; i < end; ++i) {
          p[i] = r[i] + beta * p[i];
        }
      });
    } else {
      for (int i = 0; i < n; ++i) {
        p[i] = r[i] + beta * p[i];
      }
    }

    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
