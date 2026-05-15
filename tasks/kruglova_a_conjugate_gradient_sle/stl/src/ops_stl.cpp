#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

namespace kruglova_a_conjugate_gradient_sle {

static double ParallelDotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  std::size_t n = v1.size();
  if (n == 0) {
    return 0.0;
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }
  std::size_t chunk_size = (n + num_threads - 1) / num_threads;

  std::vector<std::future<double>> futures;
  for (unsigned int i = 0; i < num_threads; ++i) {
    std::size_t start = i * chunk_size;
    std::size_t end = std::min(start + chunk_size, n);
    if (start >= n) {
      break;
    }

    futures.push_back(std::async(std::launch::async, [&v1, &v2, start, end]() {
      double sum = 0.0;
      for (std::size_t j = start; j < end; ++j) {
        sum += v1[j] * v2[j];
      }
      return sum;
    }));
  }

  double total_sum = 0.0;
  for (auto &f : futures) {
    total_sum += f.get();
  }
  return total_sum;
}

KruglovaAConjGradSleSTL::KruglovaAConjGradSleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KruglovaAConjGradSleSTL::ValidationImpl() {
  const auto &in = GetInput();
  return in.size > 0 && in.A.size() == static_cast<std::size_t>(in.size) * in.size &&
         in.b.size() == static_cast<std::size_t>(in.size);
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

  // Инициализация
  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  double rsold = ParallelDotProduct(r, r);
  const double tolerance = 1e-9;  // Ужесточили tolerance
  const int max_iter = n * 10;    // Увеличили максимальное число итераций

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    // 1. Matrix-Vector Multiply: ap = A * p
    std::vector<std::future<void>> futures;
    std::size_t chunk = (static_cast<std::size_t>(n) + num_threads - 1) / num_threads;

    for (unsigned int t = 0; t < num_threads; ++t) {
      std::size_t start = t * chunk;
      std::size_t end = std::min(start + chunk, static_cast<std::size_t>(n));
      if (start >= static_cast<std::size_t>(n)) {
        break;
      }

      futures.push_back(std::async(std::launch::async, [&a, &p, &ap, n, start, end]() {
        for (std::size_t i = start; i < end; ++i) {
          double sum = 0.0;
          std::size_t row_off = i * static_cast<std::size_t>(n);
          for (int j = 0; j < n; ++j) {
            sum += a[row_off + static_cast<std::size_t>(j)] * p[static_cast<std::size_t>(j)];
          }
          ap[i] = sum;
        }
      }));
    }
    for (auto &f : futures) {
      f.get();
    }
    futures.clear();

    double p_ap = ParallelDotProduct(p, ap);
    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    double alpha = rsold / p_ap;

    // 2. Update x and r: x = x + alpha*p, r = r - alpha*ap
    for (unsigned int t = 0; t < num_threads; ++t) {
      std::size_t start = t * chunk;
      std::size_t end = std::min(start + chunk, static_cast<std::size_t>(n));
      if (start >= static_cast<std::size_t>(n)) {
        break;
      }

      futures.push_back(std::async(std::launch::async, [&x, &r, &p, &ap, alpha, start, end]() {
        for (std::size_t i = start; i < end; ++i) {
          x[i] += alpha * p[i];
          r[i] -= alpha * ap[i];
        }
      }));
    }
    for (auto &f : futures) {
      f.get();
    }
    futures.clear();

    double rsnew = ParallelDotProduct(r, r);
    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    double beta = rsnew / rsold;

    // 3. Update p: p = r + beta * p
    for (unsigned int t = 0; t < num_threads; ++t) {
      std::size_t start = t * chunk;
      std::size_t end = std::min(start + chunk, static_cast<std::size_t>(n));
      if (start >= static_cast<std::size_t>(n)) {
        break;
      }

      futures.push_back(std::async(std::launch::async, [&p, &r, beta, start, end]() {
        for (std::size_t i = start; i < end; ++i) {
          p[i] = r[i] + beta * p[i];
        }
      }));
    }
    for (auto &f : futures) {
      f.get();
    }

    rsold = rsnew;
  }

  return true;
}

bool KruglovaAConjGradSleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kruglova_a_conjugate_gradient_sle
