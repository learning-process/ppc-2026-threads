#include "timofeev_n_radix_batcher_sort/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "timofeev_n_radix_batcher_sort/common/include/common.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

TimofeevNRadixBatcherALL::TimofeevNRadixBatcherALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool TimofeevNRadixBatcherALL::ValidationImpl() {
  return GetInput().size() >= 2;
}

bool TimofeevNRadixBatcherALL::PreProcessingImpl() {
  return true;
}

void TimofeevNRadixBatcherALL::CompExch(int &a, int &b, int digit) {
  int b_r = b % (digit * 10) / digit;
  int a_r = a % (digit * 10) / digit;
  if (b_r < a_r) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherALL::BubbleSort(std::vector<int> &arr, int digit, int left, int right) {
  for (int i = left; i <= right; i++) {
    for (int j = 0; j + 1 < right - left; j++) {
      CompExch(arr[left + j], arr[left + j + 1], digit);
    }
  }
}

void TimofeevNRadixBatcherALL::ComparR(int &a, int &b) {
  if (a > b) {
    std::swap(a, b);
  }
}

void TimofeevNRadixBatcherALL::OddEvenMerge(std::vector<int> &arr, int lft, int n) {
  if (n <= 1) {
    return;
  }

  int otstup = n / 2;
  for (int i = 0; i < otstup; i += 1) {
    if (arr[lft + i] > arr[lft + otstup + i]) {
      std::swap(arr[lft + i], arr[lft + otstup + i]);
    }
  }

  for (otstup = n / 4; otstup > 0; otstup /= 2) {
    int h = otstup * 2;
    for (int start = otstup; start + otstup < n; start += h) {
      for (int i = 0; i < otstup; i += 1) {
        ComparR(arr[lft + start + i], arr[lft + start + i + otstup]);
      }
    }
  }
}

void TimofeevNRadixBatcherALL::BubbleSortAux(size_t &num_threads, int &max_x, std::vector<int> &r_in, size_t &piece,
                                             std::vector<std::thread> &threads) {
  for (size_t t_s = 0; t_s < num_threads; ++t_s) {
    threads.emplace_back([&, t_s]() {
      for (int k = 1; k <= max_x; k *= 10) {
        BubbleSort(r_in, k, static_cast<int>(piece * t_s), static_cast<int>((piece * t_s) + piece));
      }
    });
  }
}

void TimofeevNRadixBatcherALL::OddMergeAux(std::vector<int> &r_in, size_t &piece, std::vector<std::thread> &threads,
                                           size_t &n_n) {
  for (size_t c_p = piece * 2; c_p <= n_n; c_p *= 2) {
    threads.clear();

    for (size_t i = 0; i < n_n; i += c_p) {
      threads.emplace_back([&, i]() { OddEvenMerge(r_in, static_cast<int>(i), static_cast<int>(c_p)); });
    }

    for (auto &th : threads) {
      th.join();
    }
  }
}

void TimofeevNRadixBatcherALL::ProcessLocalArray(std::vector<int> &local_arr, size_t num_threads) {
  int n = static_cast<int>(local_arr.size());
  int max_x = *std::ranges::max_element(local_arr.begin(), local_arr.end());

  std::vector<std::thread> threads;
  size_t piece = n / num_threads;
  threads.reserve(num_threads);

  BubbleSortAux(num_threads, max_x, local_arr, piece, threads);

  for (auto &th : threads) {
    th.join();
  }

  size_t nnn = n;
  OddMergeAux(local_arr, piece, threads, nnn);
}

void TimofeevNRadixBatcherALL::ProcessLocalArrayWOSort(std::vector<int> &local_arr, size_t num_threads,
                                                       size_t &elements_per_process) {
  int n = static_cast<int>(local_arr.size());

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  size_t nnn = n;
  OddMergeAux(local_arr, elements_per_process, threads, nnn);
}

void TimofeevNRadixBatcherALL::PrepAux(int &n, int &m, std::vector<int> &in, int &max_x, size_t &num_threads,
                                       size_t &n_thr, size_t &number_procs) {
  n = static_cast<int>(in.size());
  m = n;

  // дополняем
  while (n % 2 == 0) {
    n /= 2;
  }
  if (n > 1) {
    n = static_cast<int>(in.size());
    int p = 1;
    while (p < n) {
      p *= 2;
    }
    n = p;
  } else {
    n = m;
  }

  max_x = *std::ranges::max_element(in);
  if (n != m) {
    in.resize(n, max_x);
  }

  n_thr = std::thread::hardware_concurrency();
  while ((n_thr & (n_thr - 1)) != 0 && n_thr > 1) {  // и то, и то - степени двойки, поэтому достаточно вот такого цикла
    n_thr--;
  }
  if (std::cmp_greater(n_thr, static_cast<size_t>(n))) {
    n_thr = n;
  }
  num_threads = 1;
  while (num_threads * 2 <= n_thr / number_procs &&
         n / num_threads >= 4) {  // хотя бы не более, чем число n_thr / число пр-ов, и <условие с const>
    num_threads *= 2;
  }
}

void TimofeevNRadixBatcherALL::HandleZero(std::vector<int> &global_array, size_t &total_elements,
                                          size_t &total_elements_primal, int &num_processes,
                                          std::vector<int> &local_array, int &maxxx, size_t &num_threads_per_process,
                                          size_t &elements_per_process) {
  global_array.assign(GetInput().begin(), GetInput().end());

  size_t numnum = std::thread::hardware_concurrency();

  int ttttt = static_cast<int>(total_elements);
  int mmmmm = static_cast<int>(total_elements_primal);
  size_t nnnnn = num_processes;

  PrepAux(ttttt, mmmmm, global_array, maxxx, num_threads_per_process, numnum, nnnnn);

  total_elements = ttttt;
  total_elements_primal = mmmmm;
  num_processes = static_cast<int>(nnnnn);

  elements_per_process = total_elements / num_processes;

  local_array.resize(elements_per_process);

  for (int i = 0; i < num_processes; ++i) {
    if (i == 0) {
      std::copy(global_array.begin(), global_array.begin() + static_cast<long>(elements_per_process),
                local_array.begin());
    } else {
      MPI_Send(&elements_per_process, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Send(global_array.data() + static_cast<long>(i * elements_per_process),
               static_cast<int>(elements_per_process), MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
}

bool TimofeevNRadixBatcherALL::RunImpl() {
  int world_size = 0;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  size_t total_elements = GetInput().size();
  size_t total_elements_primal = total_elements;
  size_t num_threads_per_process = 1;
  int num_processes = world_size;
  while ((num_processes & (num_processes - 1)) != 0 && num_processes > 1) {
    num_processes--;
  }
  std::vector<int> global_array;
  int maxxx = 0;
  size_t elements_per_process = 0;
  std::vector<int> local_array;

  if (world_rank == 0) {
    HandleZero(global_array, total_elements, total_elements_primal, num_processes, local_array, maxxx,
               num_threads_per_process, elements_per_process);
  } else if (world_rank < num_processes) {
    MPI_Recv(&elements_per_process, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    local_array.resize(elements_per_process);
    MPI_Recv(local_array.data(), static_cast<int>(elements_per_process), MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  if (world_rank < num_processes) {
    ProcessLocalArray(local_array, num_threads_per_process);
  }

  if (world_rank == 0) {
    std::ranges::copy(local_array, global_array.begin());

    for (int i = 1; i < num_processes; ++i) {
      MPI_Recv(global_array.data() + static_cast<int>(i * elements_per_process), static_cast<int>(elements_per_process),
               MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    ProcessLocalArrayWOSort(global_array, num_threads_per_process, elements_per_process);

    if (total_elements != total_elements_primal) {
      global_array.resize(total_elements_primal);
    }

    for (int i = 1; i < world_size; ++i) {
      MPI_Send(&total_elements_primal, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Send(global_array.data(), static_cast<int>(total_elements_primal), MPI_INT, i, 0, MPI_COMM_WORLD);
    }

  } else {
    if (world_rank < num_processes) {
      MPI_Send(local_array.data(), static_cast<int>(elements_per_process), MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Recv(&total_elements_primal, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    global_array.resize(total_elements_primal);
    MPI_Recv(global_array.data(), static_cast<int>(total_elements_primal), MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
  }

  GetOutput() = global_array;
  return true;
}

bool TimofeevNRadixBatcherALL::PostProcessingImpl() {
  return true;
}

}  // namespace timofeev_n_radix_batcher_sort_threads
