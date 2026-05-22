#include "vasiliev_m_shell_sort_batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "util/include/util.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

VasilievMShellSortBatcherMergeALL::VasilievMShellSortBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool VasilievMShellSortBatcherMergeALL::ValidationImpl() {
  return !GetInput().empty();
}

bool VasilievMShellSortBatcherMergeALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool VasilievMShellSortBatcherMergeALL::RunImpl() {
  int rank = 0;
  int process_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  auto &vec = GetInput();
  int n = static_cast<int>(vec.size());

  std::vector<int> counts(process_count);
  std::vector<int> displs(process_count);

  if (rank == 0) {
    CalcCountsAndDispls(n, process_count, counts, displs);
  }

  MPI_Bcast(counts.data(), process_count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(displs.data(), process_count, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<ValType> local_data(static_cast<size_t>(counts[rank]));
  MPI_Scatterv(vec.data(), counts.data(), displs.data(), MPI_INT, local_data.data(), counts[rank], MPI_INT, 0,
               MPI_COMM_WORLD);

  const int threads = std::max(1, ppc::util::GetNumThreads());
  std::vector<size_t> bounds = ChunkBoundaries(static_cast<size_t>(counts[rank]), threads);
  size_t chunk_count = bounds.size() - 1;

  ShellSort(local_data, bounds);

  std::vector<ValType> buffer(static_cast<size_t>(counts[rank]));
  for (size_t size = 1; size < chunk_count; size *= 2) {
    CycleMerge(local_data, buffer, bounds, size);
    local_data.swap(buffer);
  }

  if (rank == 0) {
    GetOutput().resize(static_cast<size_t>(n));
  }

  MPI_Gatherv(local_data.data(), counts[rank], MPI_INT, rank == 0 ? GetOutput().data() : nullptr, counts.data(),
              displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    auto root_bounds = BoundsFromCounts(counts);
    size_t root_chunk_count = root_bounds.size() - 1;
    std::vector<ValType> root_buffer(static_cast<size_t>(n));
    for (size_t size = 1; size < root_chunk_count; size *= 2) {
      CycleMerge(GetOutput(), root_buffer, root_bounds, size);
      GetOutput().swap(root_buffer);
    }
  }

  if (rank != 0) {
    GetOutput().resize(static_cast<size_t>(n));
  }
  MPI_Bcast(GetOutput().data(), n, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool VasilievMShellSortBatcherMergeALL::PostProcessingImpl() {
  return true;
}

void VasilievMShellSortBatcherMergeALL::CalcCountsAndDispls(int n, int process_count, std::vector<int> &counts,
                                                            std::vector<int> &displs) {
  int chunk = n / process_count;
  int remainder = n % process_count;

  for (int i = 0; i < process_count; i++) {
    counts[i] = chunk + (i < remainder ? 1 : 0);
  }

  displs[0] = 0;
  for (int i = 1; i < process_count; i++) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
}

std::vector<size_t> VasilievMShellSortBatcherMergeALL::BoundsFromCounts(const std::vector<int> &counts) {
  std::vector<size_t> bounds;
  bounds.reserve(counts.size() + 1);
  bounds.push_back(0);
  for (int c : counts) {
    bounds.push_back(bounds.back() + static_cast<size_t>(c));
  }
  return bounds;
}

std::vector<size_t> VasilievMShellSortBatcherMergeALL::ChunkBoundaries(size_t vec_size, int threads) {
  size_t chunks = static_cast<size_t>(std::max(1, std::min(threads, static_cast<int>(vec_size))));
  std::vector<size_t> bounds;
  bounds.reserve(chunks + 1);

  size_t chunk_size = vec_size / chunks;
  size_t remainder = vec_size % chunks;
  bounds.push_back(0);

  for (size_t i = 0; i < chunks; i++) {
    bounds.push_back(bounds.back() + chunk_size + (i < remainder ? 1 : 0));
  }
  return bounds;
}

void VasilievMShellSortBatcherMergeALL::ShellSort(std::vector<ValType> &vec, std::vector<size_t> &bounds) {
  size_t chunk_count = bounds.size() - 1;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, chunk_count), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t chunk = range.begin(); chunk < range.end(); chunk++) {
      size_t first = bounds[chunk];
      size_t last = bounds[chunk + 1];
      size_t n = last - first;

      for (size_t gap = n / 2; gap > 0; gap /= 2) {
        for (size_t i = first + gap; i < last; i++) {
          ValType tmp = vec[i];
          size_t j = i;
          while (j >= first + gap && vec[j - gap] > tmp) {
            vec[j] = vec[j - gap];
            j -= gap;
          }
          vec[j] = tmp;
        }
      }
    }
  });
}

void VasilievMShellSortBatcherMergeALL::CycleMerge(std::vector<ValType> &vec, std::vector<ValType> &buffer,
                                                   std::vector<size_t> &bounds, size_t size) {
  const size_t chunk_count = bounds.size() - 1;
  const size_t merge_count = (chunk_count + (2 * size) - 1) / (2 * size);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, merge_count), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t idx = range.begin(); idx < range.end(); idx++) {
      const size_t l = idx * 2 * size;
      const size_t mid = std::min(l + size, chunk_count);
      const size_t r = std::min(l + (2 * size), chunk_count);

      const size_t start = bounds[l];
      const size_t middle = bounds[mid];
      const size_t end = bounds[r];

      if (mid == r) {
        std::copy(vec.begin() + static_cast<std::ptrdiff_t>(start), vec.begin() + static_cast<std::ptrdiff_t>(end),
                  buffer.begin() + static_cast<std::ptrdiff_t>(start));
      } else {
        std::vector<ValType> l_vect(vec.begin() + static_cast<std::ptrdiff_t>(start),
                                    vec.begin() + static_cast<std::ptrdiff_t>(middle));
        std::vector<ValType> r_vect(vec.begin() + static_cast<std::ptrdiff_t>(middle),
                                    vec.begin() + static_cast<std::ptrdiff_t>(end));
        std::vector<ValType> merged = BatcherMerge(l_vect, r_vect);
        for (size_t i = 0; i < merged.size(); i++) {
          buffer[start + i] = merged[i];
        }
      }
    }
  });
}

std::vector<ValType> VasilievMShellSortBatcherMergeALL::BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r) {
  std::vector<ValType> even_l;
  std::vector<ValType> odd_l;
  std::vector<ValType> even_r;
  std::vector<ValType> odd_r;

  SplitEvenOdd(l, even_l, odd_l);
  SplitEvenOdd(r, even_r, odd_r);

  std::vector<ValType> even = Merge(even_l, even_r);
  std::vector<ValType> odd = Merge(odd_l, odd_r);

  std::vector<ValType> res;
  res.reserve(l.size() + r.size());

  for (size_t i = 0; i < even.size() || i < odd.size(); i++) {
    if (i < even.size()) {
      res.push_back(even[i]);
    }
    if (i < odd.size()) {
      res.push_back(odd[i]);
    }
  }

  for (size_t i = 1; i + 1 < res.size(); i += 2) {
    if (res[i] > res[i + 1]) {
      std::swap(res[i], res[i + 1]);
    }
  }

  return res;
}

void VasilievMShellSortBatcherMergeALL::SplitEvenOdd(std::vector<ValType> &vec, std::vector<ValType> &even,
                                                     std::vector<ValType> &odd) {
  even.reserve(even.size() + (vec.size() / 2) + 1);
  odd.reserve(odd.size() + (vec.size() / 2));

  for (size_t i = 0; i < vec.size(); i += 2) {
    even.push_back(vec[i]);
    if (i + 1 < vec.size()) {
      odd.push_back(vec[i + 1]);
    }
  }
}

std::vector<ValType> VasilievMShellSortBatcherMergeALL::Merge(std::vector<ValType> &a, std::vector<ValType> &b) {
  std::vector<ValType> merged;
  merged.reserve(a.size() + b.size());
  size_t i = 0;
  size_t j = 0;

  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      merged.push_back(a[i++]);
    } else {
      merged.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    merged.push_back(a[i++]);
  }
  while (j < b.size()) {
    merged.push_back(b[j++]);
  }

  return merged;
}

}  // namespace vasiliev_m_shell_sort_batcher_merge
