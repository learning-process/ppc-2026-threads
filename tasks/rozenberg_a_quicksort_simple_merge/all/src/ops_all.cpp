#include "rozenberg_a_quicksort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <stack>
#include <utility>
#include <vector>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"

namespace rozenberg_a_quicksort_simple_merge {

RozenbergAQuicksortSimpleMergeALL::RozenbergAQuicksortSimpleMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  InType empty;
  GetInput().swap(empty);

  for (const auto &elem : in) {
    GetInput().push_back(elem);
  }

  GetOutput().clear();
}

bool RozenbergAQuicksortSimpleMergeALL::ValidationImpl() {
  return (!(GetInput().empty())) && (GetOutput().empty());
}

bool RozenbergAQuicksortSimpleMergeALL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

std::pair<int, int> RozenbergAQuicksortSimpleMergeALL::Partition(InType &data, int left, int right) {
  const int pivot = data[left + ((right - left) / 2)];
  int i = left;
  int j = right;

  while (i <= j) {
    while (data[i] < pivot) {
      i++;
    }
    while (data[j] > pivot) {
      j--;
    }

    if (i <= j) {
      std::swap(data[i], data[j]);
      i++;
      j--;
    }
  }
  return {i, j};
}

void RozenbergAQuicksortSimpleMergeALL::PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right,
                                                      int i, int j) {
  if (j - left > right - i) {
    if (left < j) {
      stack.emplace(left, j);
    }
    if (i < right) {
      stack.emplace(i, right);
    }
  } else {
    if (i < right) {
      stack.emplace(i, right);
    }
    if (left < j) {
      stack.emplace(left, j);
    }
  }
}

void RozenbergAQuicksortSimpleMergeALL::Quicksort(InType &data, int low, int high) {
  if (low >= high) {
    return;
  }

  std::stack<std::pair<int, int>> stack;

  stack.emplace(low, high);

  while (!stack.empty()) {
    const auto [left, right] = stack.top();
    stack.pop();

    if (left < right) {
      const auto [i, j] = Partition(data, left, right);
      PushSubarrays(stack, left, right, i, j);
    }
  }
}

InType RozenbergAQuicksortSimpleMergeALL::Merge(const InType &v1, const InType &v2) {
  InType res;
  res.reserve(v1.size() + v2.size());
  auto it1 = v1.begin();
  auto it2 = v2.begin();
  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 <= *it2) {
      res.push_back(*it1++);
    } else {
      res.push_back(*it2++);
    }
  }
  res.insert(res.end(), it1, v1.end());
  res.insert(res.end(), it2, v2.end());
  return res;
}

void RozenbergAQuicksortSimpleMergeALL::ThreadMerge(InType &data, int left, int mid, int right) {
  std::vector<int> temp(right - left + 1);
  int i = left;
  int j = mid + 1;
  int k = 0;

  while (i <= mid && j <= right) {
    if (data[i] <= data[j]) {
      temp[k++] = data[i++];
    } else {
      temp[k++] = data[j++];
    }
  }

  while (i <= mid) {
    temp[k++] = data[i++];
  }
  while (j <= right) {
    temp[k++] = data[j++];
  }

  for (int idx = 0; idx < k; ++idx) {
    data[left + idx] = temp[idx];
  }
}

void RozenbergAQuicksortSimpleMergeALL::ThreadQuicksort(InType &local_data) {
  int num_threads = omp_get_max_threads();
  int local_n = static_cast<int>(local_data.size());

  if (local_n > num_threads) {
    std::vector<int> thr_borders(num_threads + 1);
    int thr_chunk = local_n / num_threads;

    for (int i = 0; i < num_threads; i++) {
      thr_borders[i] = i * thr_chunk;
    }
    thr_borders[num_threads] = local_n;

#pragma omp parallel for default(none) shared(local_data, thr_borders, num_threads) num_threads(num_threads)
    for (int i = 0; i < num_threads; i++) {
      Quicksort(local_data, thr_borders[i], thr_borders[i + 1] - 1);
    }

    for (int i = 1; i < num_threads; i++) {
      ThreadMerge(local_data, 0, thr_borders[i] - 1, thr_borders[i + 1] - 1);
    }
  } else {
    Quicksort(local_data, 0, local_n - 1);
  }
}

bool RozenbergAQuicksortSimpleMergeALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n = 0;
  if (rank == 0) {
    n = static_cast<int>(GetInput().size());
  }
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> send_counts(size);
  std::vector<int> displs(size, 0);
  int chunk = n / size;
  int rem = n % size;
  for (int i = 0; i < size; i++) {
    send_counts[i] = chunk + (i < rem ? 1 : 0);
    if (i > 0) {
      displs[i] = displs[i - 1] + send_counts[i - 1];
    }
  }

  std::vector<int> local_data(send_counts[rank]);
  MPI_Scatterv(rank == 0 ? GetInput().data() : nullptr, send_counts.data(), displs.data(), MPI_INT, local_data.data(),
               send_counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

  ThreadQuicksort(local_data);

  if (rank != 0) {
    MPI_Send(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, 0, 0, MPI_COMM_WORLD);
  } else {
    InType total_data = local_data;
    for (int i = 1; i < size; ++i) {
      InType recv_buf(send_counts[i]);
      MPI_Recv(recv_buf.data(), send_counts[i], MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      auto merged = Merge(total_data, recv_buf);
      total_data.swap(merged);
    }
    GetOutput() = total_data;
  }

  return true;
}

bool RozenbergAQuicksortSimpleMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace rozenberg_a_quicksort_simple_merge
