#include "romanova_v_linear_histogram_stretch/stl/include/ops_stl.hpp"

#include <atomic>
#include <numeric>
#include <thread>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

RomanovaVLinHistogramStretchSTL::RomanovaVLinHistogramStretchSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool RomanovaVLinHistogramStretchSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchSTL::RunImpl() {
  return true;
}

bool RomanovaVLinHistogramStretchSTL::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
