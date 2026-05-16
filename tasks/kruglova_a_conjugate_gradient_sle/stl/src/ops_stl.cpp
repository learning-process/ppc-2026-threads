#include "kruglova_a_conjugate_gradient_sle/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "kruglova_a_conjugate_gradient_sle/common/include/common.hpp"

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

  std::vector<double> r = b;
  std::vector<double> p = r;
  std::vector<double> ap(n, 0.0);

  // Определение количества доступных вычислительных ядер
  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads <= 0) {
    num_threads = 2;
  }
  if (n < num_threads * 10) {
    num_threads = 1;  // Для маленьких размеров векторов считаем последовательно
  }

  // Универсальный планировщик параллельного распределения задач по потокам
  auto launch_parallel = [&](int total_elements, int threads_count, auto &&func) {
    std::vector<std::thread> workers;
    workers.reserve(threads_count);
    int chunk = total_elements / threads_count;
    int rem = total_elements % threads_count;
    int start = 0;
    for (int i = 0; i < threads_count; ++i) {
      int end = start + chunk + (i < rem ? 1 : 0);
      workers.emplace_back(func, start, end, i);
      start = end;
    }
    for (auto &w : workers) {
      w.join();
    }
  };

  // Вычисление стартового значения rsold
  double rsold = 0.0;
  if (num_threads > 1) {
    std::vector<double> partial_rs(num_threads, 0.0);
    launch_parallel(n, num_threads, [&](int start, int end, int tid) {
      double local_sum = 0.0;
      for (int i = start; i < end; ++i) {
        local_sum += r[i] * r[i];
      }
      partial_rs[tid] = local_sum;
    });
    for (double val : partial_rs) {
      rsold += val;
    }
  } else {
    for (int i = 0; i < n; ++i) {
      rsold += r[i] * r[i];
    }
  }

  const double tolerance = 1e-8;
  const int max_iter = n * 2;

  for (int iter = 0; iter < max_iter; ++iter) {
    // 1. Матрично-векторное умножение: ap = A * p
    if (num_threads > 1) {
      launch_parallel(n, num_threads, [&](int start, int end, int tid) {
        for (int i = start; i < end; ++i) {
          double sum = 0.0;
          size_t row_offset = static_cast<size_t>(i) * n;
          for (int j = 0; j < n; ++j) {
            sum += a[row_offset + j] * p[j];
          }
          ap[i] = sum;
        }
      });
    } else {
      for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        size_t row_offset = static_cast<size_t>(i) * n;
        for (int j = 0; j < n; ++j) {
          sum += a[row_offset + j] * p[j];
        }
        ap[i] = sum;
      }
    }

    // 2. Вычисление скалярного произведения: p_ap = p * ap
    double p_ap = 0.0;
    if (num_threads > 1) {
      std::vector<double> partial_p_ap(num_threads, 0.0);
      launch_parallel(n, num_threads, [&](int start, int end, int tid) {
        double local_sum = 0.0;
        for (int i = start; i < end; ++i) {
          local_sum += p[i] * ap[i];
        }
        partial_p_ap[tid] = local_sum;
      });
      for (double val : partial_p_ap) {
        p_ap += val;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        p_ap += p[i] * ap[i];
      }
    }

    if (std::abs(p_ap) < 1e-15) {
      break;
    }

    const double alpha = rsold / p_ap;

    // 3. Параллельное обновление векторов x, r и вычисление rsnew
    double rsnew = 0.0;
    if (num_threads > 1) {
      std::vector<double> partial_rs_new(num_threads, 0.0);
      launch_parallel(n, num_threads, [&](int start, int end, int tid) {
        double local_sum = 0.0;
        for (int i = start; i < end; ++i) {
          x[i] += alpha * p[i];
          r[i] -= alpha * ap[i];
          local_sum += r[i] * r[i];
        }
        partial_rs_new[tid] = local_sum;
      });
      for (double val : partial_rs_new) {
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

    // 4. Расчет коэффициента beta и обновление вектора направлений p
    const double beta = rsnew / rsold;
    if (num_threads > 1) {
      launch_parallel(n, num_threads, [&](int start, int end, int tid) {
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
