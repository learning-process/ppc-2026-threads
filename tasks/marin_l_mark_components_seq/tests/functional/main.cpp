#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "marin_l_mark_components_seq/common/include/common.hpp"
#include "marin_l_mark_components_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace marin_l_mark_components_seq {

namespace {

using Image = std::vector<std::vector<int>>;
using Labels = std::vector<std::vector<int>>;

Image MakeImage(int height, int width, int fill_value = 0) {
  return {static_cast<size_t>(height), std::vector<int>(static_cast<size_t>(width), fill_value)};
}

bool HaveSameDimensions(const Labels &first, const Labels &second) {
  if (first.size() != second.size() || first.empty()) {
    return false;
  }
  return first[0].size() == second[0].size();
}

Labels NormalizeLabels(int height, int width, const Labels &source) {
  Labels normalized = source;

  std::vector<int> label_mapping(100000, 0);
  int next_component_id = 1;

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < width; ++col_idx) {
      const int label_value = normalized[row_idx][col_idx];
      if (label_value == 0) {
        continue;
      }

      if (static_cast<size_t>(label_value) >= label_mapping.size()) {
        label_mapping.resize(static_cast<size_t>(label_value) + 1000ULL, 0);
      }

      if (label_mapping[static_cast<size_t>(label_value)] == 0) {
        label_mapping[static_cast<size_t>(label_value)] = next_component_id++;
      }

      normalized[row_idx][col_idx] = label_mapping[static_cast<size_t>(label_value)];
    }
  }
  return normalized;
}

bool SameComponentStructure(const Labels &first, const Labels &second) {
  if (!HaveSameDimensions(first, second)) {
    return false;
  }

  const int height = static_cast<int>(first.size());
  const int width = static_cast<int>(first[0].size());

  return NormalizeLabels(height, width, first) == NormalizeLabels(height, width, second);
}

bool IsInsideBounds(int height, int width, int row, int col) {
  return row >= 0 && row < height && col >= 0 && col < width;
}

void FloodFillComponent(const Image &binary_image, Labels &result_labels, int height, int width, int start_row,
                        int start_col, int component_label, std::vector<std::pair<int, int>> &stack) {
  stack.emplace_back(start_row, start_col);
  result_labels[start_row][start_col] = component_label;

  while (!stack.empty()) {
    const auto [current_row, current_col] = stack.back();
    stack.pop_back();

    const std::array<std::pair<int, int>, 4> directions{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

    for (const auto &dir : directions) {
      const int next_row = current_row + dir.first;
      const int next_col = current_col + dir.second;

      if (IsInsideBounds(height, width, next_row, next_col) && binary_image[next_row][next_col] == 1 &&
          result_labels[next_row][next_col] == 0) {
        result_labels[next_row][next_col] = component_label;
        stack.emplace_back(next_row, next_col);
      }
    }
  }
}

Labels ComputeReferenceLabels(const Image &binary) {
  const int height = static_cast<int>(binary.size());
  const int width = height != 0 ? static_cast<int>(binary[0].size()) : 0;

  Labels result(static_cast<size_t>(height), std::vector<int>(static_cast<size_t>(width), 0));

  int current_label = 0;

  std::vector<std::pair<int, int>> stack;
  stack.reserve(static_cast<size_t>(height * width / 4));

  for (int row_idx = 0; row_idx < height; ++row_idx) {
    for (int col_idx = 0; col_idx < width; ++col_idx) {
      if (binary[row_idx][col_idx] == 1 && result[row_idx][col_idx] == 0) {
        ++current_label;
        FloodFillComponent(binary, result, height, width, row_idx, col_idx, current_label, stack);
      }
    }
  }
  return result;
}

void FillSingleBlob(Image &img, int height, int width) {
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      img[row][col] = 1;
    }
  }
}

void FillTwoBlocks(Image &img, int height, int width) {
  for (int row = 1; row < (height / 2); ++row) {
    for (int col = 1; col < (width / 2); ++col) {
      img[row][col] = 1;
    }
  }

  for (int row = (height / 2) + 1; row < height - 1; ++row) {
    for (int col = (width / 2) + 1; col < width - 1; ++col) {
      img[row][col] = 1;
    }
  }
}

void FillChecker(Image &img, int height, int width) {
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      img[row][col] = (row + col) % 2;
    }
  }
}

void FillDiagonal(Image &img, int height, int width) {
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      img[row][col] = (row % 3 == col % 3) ? 1 : 0;
    }
  }
}

void FillRandom(Image &img, int height, int width, std::mt19937 &gen, std::uniform_real_distribution<double> &dist) {
  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      img[row][col] = (dist(gen) < 0.25) ? 1 : 0;
    }
  }
}

}  // namespace

class MarinLRunFuncTestComponents : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "x" + std::to_string(std::get<1>(test_param)) + "_" +
           std::get<2>(test_param);
  }

 protected:
  void SetUp() override {
    const auto params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const int width = std::get<0>(params);
    const int height = std::get<1>(params);
    const std::string &scenario_name = std::get<2>(params);

    input_data_.binary = MakeImage(height, width);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> probability_dist(0.0, 1.0);

    if (scenario_name == "SingleBlob") {
      FillSingleBlob(input_data_.binary, height, width);
    } else if (scenario_name == "TwoBlocks") {
      FillTwoBlocks(input_data_.binary, height, width);
    } else if (scenario_name == "Checker") {
      FillChecker(input_data_.binary, height, width);
    } else if (scenario_name == "Diagonal") {
      FillDiagonal(input_data_.binary, height, width);
    } else if (scenario_name == "Random") {
      FillRandom(input_data_.binary, height, width, gen, probability_dist);
    }

    expected_output_.labels = ComputeReferenceLabels(input_data_.binary);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return SameComponentStructure(expected_output_.labels, output_data.labels);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  OutType expected_output_{};
};

namespace {

TEST_P(MarinLRunFuncTestComponents, MarkComponentsSeq) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParams{std::make_tuple(5, 5, "SingleBlob"), std::make_tuple(8, 10, "TwoBlocks"),
                                          std::make_tuple(7, 7, "Checker"),    std::make_tuple(12, 12, "Diagonal"),
                                          std::make_tuple(15, 10, "Random"),   std::make_tuple(4, 6, "TwoBlocks")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MarinLMarkComponentsSEQ, InType>(kTestParams, PPC_SETTINGS_marin_l_mark_components_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = MarinLRunFuncTestComponents::PrintFuncTestName<MarinLRunFuncTestComponents>;

INSTANTIATE_TEST_SUITE_P(ComponentLabelingTests, MarinLRunFuncTestComponents, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace marin_l_mark_components_seq
