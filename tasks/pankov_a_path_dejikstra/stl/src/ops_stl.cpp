#include "pankov_a_path_dejikstra/stl/include/ops_stl.hpp"

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "pankov_a_path_dejikstra/common/include/common.hpp"
#include "util/include/util.hpp"

namespace pankov_a_path_dejikstra {

PankovAPathDejikstraSTL::PankovAPathDejikstraSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool PankovAPathDejikstraSTL::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool PankovAPathDejikstraSTL::PreProcessingImpl() {
  GetOutput() = 2 * GetInput();
  return GetOutput() > 0;
}

bool PankovAPathDejikstraSTL::RunImpl() {
  for (InType i = 0; i < GetInput(); i++) {
    for (InType j = 0; j < GetInput(); j++) {
      for (InType k = 0; k < GetInput(); k++) {
        std::vector<InType> tmp(i + j + k, 1);
        GetOutput() += std::accumulate(tmp.begin(), tmp.end(), 0);
        GetOutput() -= i + j + k;
      }
    }
  }

  const int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads(num_threads);
  GetOutput() *= num_threads;

  std::atomic<int> counter(0);
  for (int i = 0; i < num_threads; i++) {
    threads[i] = std::thread([&]() { counter++; });
    threads[i].join();
  }

  GetOutput() /= counter;
  return GetOutput() > 0;
}

bool PankovAPathDejikstraSTL::PostProcessingImpl() {
  GetOutput() -= GetInput();
  return GetOutput() > 0;
}

}  // namespace pankov_a_path_dejikstra
