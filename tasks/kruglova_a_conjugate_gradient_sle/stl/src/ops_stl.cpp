#include <algorithm>
#include <cmath>
#include <future>
#include <numeric>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"
#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

namespace kruglova_a_conjugate_gradient_sle {

double ParallelDotProduct(const std::vector<double> &v1, const std::vector<double> &v2) {
  size_t n = v1.size();
  if (n == 0) {
    return 0.0;
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }

  size_t chunk_size = (n + num_threads - 1) / num_threads;
  std::vector<std::future<double>> futures;
  futures.reserve(num_threads);

  for (unsigned int t = 0; t < num_threads; ++t) {
    size_t start = t * chunk_size;
    size_t end = std::min(start + chunk_size, n);
    if (start >= n) {
      break;
    }

    futures.emplace_back(std::async(std::launch::async, [&v1, &v2, start, end]() -> double {
      double sum = 0.0;
      for (size_t i = start; i < end; ++i) {
        sum += v1[i] * v2[i];
      }
      return sum;
    }));
  }

  double total = 0.0;
  for (auto &f : futures) {
    total += f.get();
  }
  return total;
}

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
  const auto &A = GetInput().A;
  const auto &b = GetInput().b;
  const int n = GetInput().size;
  auto &x = GetOutput();

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> Ap(n, 0.0);

  double rsold = ParallelDotProduct(r, r);
  const double tolerance = 1e-6;
  const int max_iter = n * 6;

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 4;
  }
  const size_t chunk = (static_cast<size_t>(n) + num_threads - 1) / num_threads;

  for (int iter = 0; iter < max_iter; ++iter) {
    // A * p → Ap
    {
      std::vector<std::future<void>> futures;
      for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(start + chunk, static_cast<size_t>(n));
        if (start >= static_cast<size_t>(n)) {
          break;
        }

        futures.emplace_back(std::async(std::launch::async, [&A, &p, &Ap, n, start, end]() {
          for (size_t i = start; i < end; ++i) {
            double sum = 0.0;
            size_t row = i * static_cast<size_t>(n);
            for (int j = 0; j < n; ++j) {
              sum += A[row + j] * p[j];
            }
            Ap[i] = sum;
          }
        }));
      }
      for (auto &f : futures) {
        f.get();
      }
    }

    double pAp = ParallelDotProduct(p, Ap);
    if (std::abs(pAp) < 1e-14) {
      break;
    }

    double alpha = rsold / pAp;

    {
      std::vector<std::future<void>> futures;
      for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(start + chunk, static_cast<size_t>(n));
        if (start >= static_cast<size_t>(n)) {
          break;
        }

        futures.emplace_back(std::async(std::launch::async, [&x, &r, &p, &Ap, alpha, start, end]() {
          for (size_t i = start; i < end; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
          }
        }));
      }
      for (auto &f : futures) {
        f.get();
      }
    }

    double rsnew = ParallelDotProduct(r, r);
    if (std::sqrt(rsnew) < tolerance) {
      break;
    }

    double beta = rsnew / rsold;

    {
      std::vector<std::future<void>> futures;
      for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(start + chunk, static_cast<size_t>(n));
        if (start >= static_cast<size_t>(n)) {
          break;
        }

        futures.emplace_back(std::async(std::launch::async, [&p, &r, beta, start, end]() {
          for (size_t i = start; i < end; ++i) {
            p[i] = r[i] + beta * p[i];
          }
        }));
      }
      for (auto &f : futures) {
        f.get();
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
