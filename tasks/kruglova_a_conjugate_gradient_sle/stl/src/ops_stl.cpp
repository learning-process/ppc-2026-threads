#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <numeric>
#include <thread>
#include <vector>

namespace kruglova_a_conjugate_gradient_sle {

namespace {

template <typename Func>
void launch_parallel_tasks(int total, int num_threads, Func &&func) {
  if (num_threads <= 1) {
    func(0, total, 0);
    return;
  }

  std::vector<std::thread> workers;
  int chunk = total / num_threads;
  workers.reserve(num_threads);

  for (int i = 0; i < num_threads; ++i) {
    int start = i * chunk;
    int end = (i == num_threads - 1) ? total : (i + 1) * chunk;
    workers.emplace_back(func, start, end, i);
  }
  for (auto &w : workers) {
    w.join();
  }
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

  double rsold = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
  const double tolerance = 1e-8;

  for (int iter = 0; iter < n; ++iter) {
    launch_parallel_tasks(n, num_threads, [&](int start, int end, int) {
      for (int i = start; i < end; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
          sum += a[i * n + j] * p[j];
        }
        ap[i] = sum;
      }
    });

    double p_ap = 0.0;
    launch_parallel_tasks(n, num_threads, [&](int start, int end, int tid) {
      double local = 0.0;
      for (int i = start; i < end; ++i) {
        local += p[i] * ap[i];
      }
      partial_buffer[tid] = local;
    });
    p_ap = std::accumulate(partial_buffer.begin(), partial_buffer.begin() + num_threads, 0.0);

    if (std::abs(p_ap) < 1e-16) {
      break;
    }

    const double alpha = rsold / p_ap;

    double rsnew = 0.0;
    launch_parallel_tasks(n, num_threads, [&](int start, int end, int tid) {
      double local_rs = 0.0;
      for (int i = start; i < end; ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * ap[i];
        local_rs += r[i] * r[i];
      }
      partial_buffer[tid] = local_rs;
    });
    rsnew = std::accumulate(partial_buffer.begin(), partial_buffer.begin() + num_threads, 0.0);

    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    const double beta = rsnew / rsold;
    launch_parallel_tasks(n, num_threads, [&](int start, int end, int) {
      for (int i = start; i < end; ++i) {
        p[i] = r[i] + (beta * p[i]);
      }
    });

    rsold = rsnew;
  }
  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
