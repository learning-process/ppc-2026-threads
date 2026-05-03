#include "vasiliev_m_shell_sort_batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "vasiliev_m_shell_sort_batcher_merge/common/include/common.hpp"

namespace vasiliev_m_shell_sort_batcher_merge {

VasilievMShellSortBatcherMergeSTL::VasilievMShellSortBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool VasilievMShellSortBatcherMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool VasilievMShellSortBatcherMergeSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool VasilievMShellSortBatcherMergeSTL::RunImpl() {
  auto &vec = GetInput();
  const size_t n = vec.size();

  if (vec.empty()) {
    return false;
  }

  int threads = std::max(1, ppc::util::GetNumThreads());

  std::vector<size_t> bounds = ChunkBoundaries(n, threads);
  size_t chunk_count = bounds.size() - 1;

  ShellSort(vec, bounds, threads);

  std::vector<ValType> buffer(n);
  for (size_t size = 1; size < chunk_count; size *= 2) {
    CycleMerge(vec, buffer, bounds, size, threads);
    vec.swap(buffer);
  }

  GetOutput() = vec;
  return true;
}

bool VasilievMShellSortBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

std::vector<size_t> VasilievMShellSortBatcherMergeSTL::ChunkBoundaries(size_t vec_size, int threads) {
  size_t chunks = std::max(1, std::min(threads, static_cast<int>(vec_size)));

  std::vector<size_t> bounds;
  bounds.reserve(chunks + 1);

  size_t chunk_size = vec_size / chunks;
  size_t remainder = vec_size % chunks;

  bounds.push_back(0);

  for (size_t i = 0; i < chunks; i++) {
    if (i < remainder) {
      bounds.push_back(bounds.back() + chunk_size + 1);
    } else {
      bounds.push_back(bounds.back() + chunk_size);
    }
  }
  return bounds;
}

void VasilievMShellSortBatcherMergeSTL::ShellSort(std::vector<ValType> &vec, std::vector<size_t> &bounds,
                                                  size_t num_threads) {
  const size_t chunk_count = bounds.size() - 1;
  const size_t worker_count = std::min(chunk_count, num_threads);
  std::vector<std::thread> threads;
  threads.reserve(worker_count);

  const size_t base = chunk_count / worker_count;
  const size_t rem = chunk_count % worker_count;
  size_t current = 0;
  for (size_t thread_id = 0; thread_id < worker_count; thread_id++) {
    const size_t count = base + (thread_id < rem ? 1 : 0);
    const size_t begin = current;
    const size_t end = current + count;
    current = end;

    threads.emplace_back([&, begin, end]() {
      for (size_t chunk = begin; chunk < end; chunk++) {
        ShellSortChunk(vec, bounds[chunk], bounds[chunk + 1]);
      }
    });
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
}

void VasilievMShellSortBatcherMergeSTL::CycleMerge(std::vector<ValType> &vec, std::vector<ValType> &buffer,
                                                   std::vector<size_t> &bounds, size_t size, size_t num_threads) {
  const size_t chunk_count = bounds.size() - 1;
  const size_t merge_count = (chunk_count + (2 * size) - 1) / (2 * size);
  const size_t worker_count = std::min(merge_count, num_threads);
  std::vector<std::thread> threads;
  threads.reserve(worker_count);

  const size_t base = merge_count / worker_count;
  const size_t rem = merge_count % worker_count;
  size_t current = 0;
  for (size_t thread_id = 0; thread_id < worker_count; thread_id++) {
    const size_t count = base + (thread_id < rem ? 1 : 0);
    const size_t begin = current;
    const size_t end = current + count;
    current = end;

    threads.emplace_back([&, begin, end]() {
      for (size_t idx = begin; idx < end; idx++) {
        const size_t l = idx * 2 * size;
        const size_t mid = std::min(l + size, chunk_count);
        const size_t r = std::min(l + (2 * size), chunk_count);

        const size_t start = bounds[l];
        const size_t middle = bounds[mid];
        const size_t end_pos = bounds[r];
        MergeChunks(vec, buffer, start, middle, end_pos);
      }
    });
  }

  for (auto &th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
}

std::vector<ValType> VasilievMShellSortBatcherMergeSTL::BatcherMerge(std::vector<ValType> &l, std::vector<ValType> &r) {
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

void VasilievMShellSortBatcherMergeSTL::SplitEvenOdd(std::vector<ValType> &vec, std::vector<ValType> &even,
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

std::vector<ValType> VasilievMShellSortBatcherMergeSTL::Merge(std::vector<ValType> &a, std::vector<ValType> &b) {
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

void VasilievMShellSortBatcherMergeSTL::ShellSortChunk(std::vector<ValType> &vec, size_t first, size_t last) {
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

void VasilievMShellSortBatcherMergeSTL::MergeChunks(std::vector<ValType> &vec, std::vector<ValType> &buffer,
                                                    size_t start, size_t middle, size_t end_pos) {
  if (middle == end_pos) {
    std::copy(vec.begin() + static_cast<std::ptrdiff_t>(start), vec.begin() + static_cast<std::ptrdiff_t>(end_pos),
              buffer.begin() + static_cast<std::ptrdiff_t>(start));
    return;
  }
  std::vector<ValType> l_vect(vec.begin() + static_cast<std::ptrdiff_t>(start),
                              vec.begin() + static_cast<std::ptrdiff_t>(middle));
  std::vector<ValType> r_vect(vec.begin() + static_cast<std::ptrdiff_t>(middle),
                              vec.begin() + static_cast<std::ptrdiff_t>(end_pos));

  std::vector<ValType> merged = BatcherMerge(l_vect, r_vect);
  for (size_t i = 0; i < merged.size(); i++) {
    buffer[start + i] = merged[i];
  }
}

}  // namespace vasiliev_m_shell_sort_batcher_merge
