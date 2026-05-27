#include "spichek_d_radix_sort_for_integers_with_simple_merging/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "util/include/util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

SpichekDRadixSortALL::SpichekDRadixSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SpichekDRadixSortALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return !GetInput().empty();
  }
  return true;
}

bool SpichekDRadixSortALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetOutput() = GetInput();
  }
  return true;
}

bool SpichekDRadixSortALL::RunImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0 && !GetOutput().empty()) {
    int num_threads = ppc::util::GetNumThreads();
    int n = static_cast<int>(GetOutput().size());

    if (num_threads <= 1 || n < num_threads) {
      RadixSort(GetOutput());
    } else {
      std::vector<std::vector<int>> local_data = SplitData(GetOutput(), num_threads);

      tbb::parallel_for(0, num_threads, [&](int i) { RadixSort(local_data[i]); });

      GetOutput() = MergeData(local_data);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool SpichekDRadixSortALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return std::ranges::is_sorted(GetOutput());
  }
  return true;
}

std::vector<std::vector<int>> SpichekDRadixSortALL::SplitData(const std::vector<int> &data, int num_parts) {
  int n = static_cast<int>(data.size());
  std::vector<std::vector<int>> parts(num_parts);
  int base_size = n / num_parts;
  int remainder = n % num_parts;
  int offset = 0;

  for (int i = 0; i < num_parts; ++i) {
    int current_size = base_size + (i < remainder ? 1 : 0);
    parts[i] = std::vector<int>(data.begin() + offset, data.begin() + offset + current_size);
    offset += current_size;
  }
  return parts;
}

std::vector<int> SpichekDRadixSortALL::MergeData(std::vector<std::vector<int>> &parts) {
  std::vector<int> result = std::move(parts[0]);
  int num_parts = static_cast<int>(parts.size());
  for (int i = 1; i < num_parts; ++i) {
    std::vector<int> temp;
    temp.reserve(result.size() + parts[i].size());
    std::ranges::merge(result, parts[i], std::back_inserter(temp));
    result = std::move(temp);
  }
  return result;
}

void SpichekDRadixSortALL::RadixSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  int min_val = *std::ranges::min_element(data);
  if (min_val < 0) {
    for (auto &x : data) {
      x -= min_val;
    }
  }

  int max_val = *std::ranges::max_element(data);
  for (int shift = 0; (max_val >> shift) > 0; shift += 8) {
    std::vector<int> output(data.size());
    std::vector<int> count(256, 0);

    for (int x : data) {
      count[(x >> shift) & 255]++;
    }
    for (int i = 1; i < 256; i++) {
      count[i] += count[i - 1];
    }
    for (int i = static_cast<int>(data.size()) - 1; i >= 0; i--) {
      output[count[(data[i] >> shift) & 255] - 1] = data[i];
      count[(data[i] >> shift) & 255]--;
    }
    data = std::move(output);
  }

  if (min_val < 0) {
    for (auto &x : data) {
      x += min_val;
    }
  }
}

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
