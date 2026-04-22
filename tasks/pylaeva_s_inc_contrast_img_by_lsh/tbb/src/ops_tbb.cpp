#include "pylaeva_s_inc_contrast_img_by_lsh/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

PylaevaSIncContrastImgByLshTBB::PylaevaSIncContrastImgByLshTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool PylaevaSIncContrastImgByLshTBB::ValidationImpl() {
  return !(GetInput().empty());
}

bool PylaevaSIncContrastImgByLshTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool PylaevaSIncContrastImgByLshTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  size_t size = input.size();

  if (size == 0) {
    return false;
  }

  using MinMax = std::pair<uint8_t, uint8_t>;

  auto [min_pixel, max_pixel] = tbb::parallel_reduce(tbb::blocked_range<size_t>(1, size), MinMax{input[0], input[0]},
                                                     [&](const auto &range, MinMax init) -> MinMax {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      init.first = std::min(init.first, input[i]);
      init.second = std::max(init.second, input[i]);
    }
    return init;
  }, [](const MinMax &a, const MinMax &b) -> MinMax {
    return {std::min(a.first, b.first), std::max(a.second, b.second)};
  });

  if (min_pixel == max_pixel) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const auto &r) {
      const uint8_t *src = input.data() + r.begin();
      uint8_t *dst = output.data() + r.begin();
      std::copy(src, src + r.size(), dst);
    }, tbb::simple_partitioner());
    return true;
  }

  const float scale = 255.0F / static_cast<float>(max_pixel - min_pixel);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const auto &r) {
    for (auto i = r.begin(); i != r.end(); ++i) {
      output[i] = static_cast<uint8_t>(std::round((input[i] - min_pixel) * scale));
    }
  }, tbb::static_partitioner());

  return true;
}

bool PylaevaSIncContrastImgByLshTBB::PostProcessingImpl() {
  return true;
}

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
