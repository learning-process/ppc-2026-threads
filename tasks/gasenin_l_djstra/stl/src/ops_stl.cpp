#include "gasenin_l_djstra/stl/include/ops_stl.hpp"

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "gasenin_l_djstra/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gasenin_l_djstra {

GaseninLDjstraSTL::GaseninLDjstraSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool GaseninLDjstraSTL::ValidationImpl() {
  return (GetInput() > 0) && (GetOutput() == 0);
}

bool GaseninLDjstraSTL::PreProcessingImpl() {
  GetOutput() = 2 * GetInput();
  return GetOutput() > 0;
}

bool GaseninLDjstraSTL::RunImpl() {
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

bool GaseninLDjstraSTL::PostProcessingImpl() {
  GetOutput() -= GetInput();
  return GetOutput() > 0;
}

}  // namespace gasenin_l_djstra
