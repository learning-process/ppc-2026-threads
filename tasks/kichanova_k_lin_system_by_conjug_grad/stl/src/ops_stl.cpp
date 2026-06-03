#include "kichanova_k_lin_system_by_conjug_grad/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "kichanova_k_lin_system_by_conjug_grad/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kichanova_k_lin_system_by_conjug_grad {

namespace {

class ThreadPool {
 public:
  static ThreadPool &Instance() {
    static ThreadPool pool(ppc::util::GetNumThreads());
    return pool;
  }

  void ParallelFor(int n, const std::function<void(int, int)> &func) {
    if (n <= 0) {
      return;
    }

    int num_threads = static_cast<int>(threads_.size());
    if (num_threads == 0) {
      func(0, n);
      return;
    }

    std::vector<std::future<void>> futures;
    futures.reserve(num_threads);

    int chunk_size = n / num_threads;
    int remainder = n % num_threads;

    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      int start = thread_idx * chunk_size;
      int end = start + chunk_size;
      if (thread_idx == num_threads - 1) {
        end += remainder;
      }

      futures.push_back(Enqueue([func, start, end]() {
        if (start < end) {
          func(start, end);
        }
      }));
    }

    for (auto &future : futures) {
      future.get();
    }
  }

  double ParallelReduce(int n, const std::function<double(int, int)> &func) {
    if (n <= 0) {
      return 0.0;
    }

    int num_threads = static_cast<int>(threads_.size());
    if (num_threads == 0) {
      return func(0, n);
    }

    std::vector<std::future<double>> futures;
    futures.reserve(num_threads);

    int chunk_size = n / num_threads;
    int remainder = n % num_threads;

    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      int start = thread_idx * chunk_size;
      int end = start + chunk_size;
      if (thread_idx == num_threads - 1) {
        end += remainder;
      }

      futures.push_back(Enqueue([func, start, end]() {
        if (start >= end) {
          return 0.0;
        }
        return func(start, end);
      }));
    }

    double result = 0.0;
    for (auto &future : futures) {
      result += future.get();
    }
    return result;
  }

 private:
  explicit ThreadPool(size_t num_threads) {
    if (num_threads == 0) {
      num_threads = 1;
    }
    for (size_t i = 0; i < num_threads; ++i) {
      threads_.emplace_back([this]() { Worker(); });
    }
  }

  ~ThreadPool() {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      stop_ = true;
    }
    condition_.notify_all();
    for (auto &thread : threads_) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  template <typename F>
  auto Enqueue(F &&f) -> std::future<decltype(f())> {
    auto task = std::make_shared<std::packaged_task<decltype(f())()>>(std::forward<F>(f));
    std::future<decltype(f())> result = task->get_future();
    {
      std::unique_lock<std::mutex> lock(mutex_);
      if (stop_) {
        throw std::runtime_error("Enqueue on stopped ThreadPool");
      }
      tasks_.emplace([task]() { (*task)(); });
    }
    condition_.notify_one();
    return result;
  }

  void Worker() {
    while (true) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        condition_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
        if (stop_ && tasks_.empty()) {
          return;
        }
        if (tasks_.empty()) {
          continue;
        }
        task = std::move(tasks_.front());
        tasks_.pop();
      }
      if (task) {
        task();
      }
    }
  }

  std::vector<std::thread> threads_;
  std::queue<std::function<void()>> tasks_;
  std::mutex mutex_;
  std::condition_variable condition_;
  bool stop_ = false;
};

double ComputeDotProduct(const std::vector<double> &a, const std::vector<double> &b, int n) {
  if (n <= 0) {
    return 0.0;
  }
  return ThreadPool::Instance().ParallelReduce(n, [&](int start, int end) {
    double sum = 0.0;
    for (int i = start; i < end; ++i) {
      sum += a[i] * b[i];
    }
    return sum;
  });
}

void ComputeMatrixVectorProduct(const std::vector<double> &a, const std::vector<double> &v, std::vector<double> &result,
                                int n) {
  if (n <= 0) {
    return;
  }
  const auto stride = static_cast<size_t>(n);
  ThreadPool::Instance().ParallelFor(n, [&](int start, int end) {
    for (int i = start; i < end; ++i) {
      double sum = 0.0;
      const double *a_row = &a[static_cast<size_t>(i) * stride];
      for (int j = 0; j < n; ++j) {
        sum += a_row[j] * v[j];
      }
      result[i] = sum;
    }
  });
}

void UpdateSolution(std::vector<double> &x, const std::vector<double> &p, double alpha, int n) {
  if (n <= 0) {
    return;
  }
  ThreadPool::Instance().ParallelFor(n, [&](int start, int end) {
    for (int i = start; i < end; ++i) {
      x[i] += alpha * p[i];
    }
  });
}

void UpdateResidual(std::vector<double> &r, const std::vector<double> &ap, double alpha, int n) {
  if (n <= 0) {
    return;
  }
  ThreadPool::Instance().ParallelFor(n, [&](int start, int end) {
    for (int i = start; i < end; ++i) {
      r[i] -= alpha * ap[i];
    }
  });
}

void UpdateSearchDirection(std::vector<double> &p, const std::vector<double> &r, double beta, int n) {
  if (n <= 0) {
    return;
  }
  ThreadPool::Instance().ParallelFor(n, [&](int start, int end) {
    for (int i = start; i < end; ++i) {
      p[i] = r[i] + (beta * p[i]);
    }
  });
}

}  // namespace

KichanovaKLinSystemByConjugGradSTL::KichanovaKLinSystemByConjugGradSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KichanovaKLinSystemByConjugGradSTL::ValidationImpl() {
  const InType &input_data = GetInput();
  if (input_data.n <= 0) {
    return false;
  }
  if (input_data.A.size() != static_cast<size_t>(input_data.n) * static_cast<size_t>(input_data.n)) {
    return false;
  }
  if (input_data.b.size() != static_cast<size_t>(input_data.n)) {
    return false;
  }
  return true;
}

bool KichanovaKLinSystemByConjugGradSTL::PreProcessingImpl() {
  GetOutput().assign(GetInput().n, 0.0);
  return true;
}

bool KichanovaKLinSystemByConjugGradSTL::RunImpl() {
  const InType &input_data = GetInput();
  OutType &x = GetOutput();

  int n = input_data.n;
  if (n <= 0) {
    return false;
  }

  const std::vector<double> &a = input_data.A;
  const std::vector<double> &b = input_data.b;
  double epsilon = input_data.epsilon;

  std::vector<double> r(n);
  std::vector<double> p(n);
  std::vector<double> ap(n);

  for (int i = 0; i < n; i++) {
    r[i] = b[i];
    p[i] = r[i];
  }

  double rr_old = ComputeDotProduct(r, r, n);
  double residual_norm = std::sqrt(rr_old);
  if (residual_norm < epsilon) {
    return true;
  }

  int max_iter = n * 1000;
  for (int iter = 0; iter < max_iter; iter++) {
    ComputeMatrixVectorProduct(a, p, ap, n);

    double p_ap = ComputeDotProduct(p, ap, n);
    if (std::abs(p_ap) < 1e-30) {
      break;
    }

    double alpha = rr_old / p_ap;
    UpdateSolution(x, p, alpha, n);
    UpdateResidual(r, ap, alpha, n);

    double rr_new = ComputeDotProduct(r, r, n);
    residual_norm = std::sqrt(rr_new);
    if (residual_norm < epsilon) {
      break;
    }

    double beta = rr_new / rr_old;
    UpdateSearchDirection(p, r, beta, n);

    rr_old = rr_new;
  }

  return true;
}

bool KichanovaKLinSystemByConjugGradSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kichanova_k_lin_system_by_conjug_grad
