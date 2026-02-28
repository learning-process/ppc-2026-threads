#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"

#include "leonova_a_radix_merge_sort/seq/include/ops_seq.hpp"

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace leonova_a_radix_merge_sort {

namespace {

// Вспомогательные функции для уменьшения cognitive complexity

std::string GetSizeCategory(size_t size) {
  if (size == 1) {
    return "_single";
  }
  if (size <= 32) {
    return "_small";
  }
  if (size <= 100) {
    return "_medium";
  }
  return "_large";
}

void AnalyzeArrayProperties(const std::vector<int64_t>& input, 
                           bool& has_negative, bool& all_same, 
                           bool& is_sorted, bool& is_reverse) {
  has_negative = false;
  all_same = true;
  is_sorted = true;
  is_reverse = true;
  
  if (input.empty()) {
    return;
  }
  
  for (size_t i = 0; i < input.size(); ++i) {
    if (input[i] < 0) {
      has_negative = true;
    }
    
    if (i > 0) {
      if (input[i] != input[i - 1]) {
        all_same = false;
      }
      if (input[i] < input[i - 1]) {
        is_sorted = false;
      }
      if (input[i] > input[i - 1]) {
        is_reverse = false;
      }
    }
  }
}

std::string GetPropertyString(bool has_negative, bool all_same, 
                              bool is_sorted, bool is_reverse, size_t size) {
  std::string result;
  
  if (all_same && size > 1) {
    result += "_allsame";
  }
  if (is_sorted && size > 1) {
    result += "_sorted";
  }
  if (is_reverse && size > 1) {
    result += "_reverse";
  }
  if (has_negative) {
    result += "_withneg";
  }
  
  return result;
}

std::string GetThresholdString(size_t size) {
  if (size == 32) {
    return "_threshold";
  }
  if (size == 33) {
    return "_abovethreshold";
  }
  return "";
}

}  // namespace

class LeonovaARadixMergeSortRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    const auto &input = std::get<0>(test_param);

    static int test_counter = 0;
    std::string name = "size_" + std::to_string(input.size());
    name += "_test" + std::to_string(++test_counter);
    
    // Добавляем категорию размера
    name += GetSizeCategory(input.size());
    
    // Анализируем свойства массива - инициализируем каждую переменную отдельно
    bool has_negative = false;
    bool all_same = true;
    bool is_sorted = true;
    bool is_reverse = true;
    
    AnalyzeArrayProperties(input, has_negative, all_same, is_sorted, is_reverse);
    
    // Добавляем свойства
    name += GetPropertyString(has_negative, all_same, is_sorted, is_reverse, input.size());
    
    // Добавляем информацию о пороге
    name += GetThresholdString(input.size());

    return name;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    return output_data == expected_output_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(LeonovaARadixMergeSortRunFuncTests, RadixMergeSort) {
  ExecuteTest(GetParam());
}

// Тестовые данные для целых чисел
const std::array<TestType, 24> kTestParam = {
    // 1 элемент
    std::make_tuple(std::vector<int64_t>{42}, 
                    std::vector<int64_t>{42}),
    
    // 2 элемента
    std::make_tuple(std::vector<int64_t>{3, 2}, 
                    std::vector<int64_t>{2, 3}),
    
    // 5 элементов, случайные
    std::make_tuple(std::vector<int64_t>{5, 1, 4, 2, 3}, 
                    std::vector<int64_t>{1, 2, 3, 4, 5}),
    
    // С отрицательными
    std::make_tuple(std::vector<int64_t>{-5, 0, 5, -10, 10}, 
                    std::vector<int64_t>{-10, -5, 0, 5, 10}),
    
    // Дубликаты
    std::make_tuple(std::vector<int64_t>{2, 2, 1, 1, 3, 3}, 
                    std::vector<int64_t>{1, 1, 2, 2, 3, 3}),
    
    // Уже отсортированный
    std::make_tuple(std::vector<int64_t>{1, 2, 3, 4, 5}, 
                    std::vector<int64_t>{1, 2, 3, 4, 5}),
    
    // Обратно отсортированный
    std::make_tuple(std::vector<int64_t>{5, 4, 3, 2, 1}, 
                    std::vector<int64_t>{1, 2, 3, 4, 5}),
    
    // Крайние значения int64_t
    std::make_tuple(std::vector<int64_t>{INT64_MAX, INT64_MIN, 0, -1, 1}, 
                    std::vector<int64_t>{INT64_MIN, -1, 0, 1, INT64_MAX}),
    
    // Большие числа
    std::make_tuple(std::vector<int64_t>{1000000, -1000000, 500000, -500000}, 
                    std::vector<int64_t>{-1000000, -500000, 500000, 1000000}),
    
    // Ровно 32 элемента (порог RadixSort)
    std::make_tuple(
        []() {
            std::vector<int64_t> v(32);
            for (int i = 0; i < 32; ++i) {
              v[i] = 31 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(32);
            for (int i = 0; i < 32; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // 33 элемента (чуть выше порога)
    std::make_tuple(
        []() {
            std::vector<int64_t> v(33);
            for (int i = 0; i < 33; ++i) {
              v[i] = 32 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(33);
            for (int i = 0; i < 33; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // 50 элементов
    std::make_tuple(
        []() {
            std::vector<int64_t> v(50);
            for (int i = 0; i < 50; ++i) {
              v[i] = 49 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(50);
            for (int i = 0; i < 50; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // 100 элементов
    std::make_tuple(
        []() {
            std::vector<int64_t> v(100);
            for (int i = 0; i < 100; ++i) {
              v[i] = 99 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(100);
            for (int i = 0; i < 100; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // Все отрицательные
    std::make_tuple(std::vector<int64_t>{-5, -2, -8, -1, -3}, 
                    std::vector<int64_t>{-8, -5, -3, -2, -1}),
    
    // Все положительные
    std::make_tuple(std::vector<int64_t>{5, 2, 8, 1, 3}, 
                    std::vector<int64_t>{1, 2, 3, 5, 8}),
    
    // Чередование положительных и отрицательных
    std::make_tuple(std::vector<int64_t>{1, -1, 2, -2, 3, -3, 4, -4}, 
                    std::vector<int64_t>{-4, -3, -2, -1, 1, 2, 3, 4}),
    
    // Все одинаковые
    std::make_tuple(std::vector<int64_t>(25, 7), 
                    std::vector<int64_t>(25, 7)),
    
    // Степени двойки
    std::make_tuple(std::vector<int64_t>{1, 2, 4, 8, 16, 32, 64, 128}, 
                    std::vector<int64_t>{1, 2, 4, 8, 16, 32, 64, 128}),
    
    // Степени двойки в обратном порядке
    std::make_tuple(std::vector<int64_t>{128, 64, 32, 16, 8, 4, 2, 1}, 
                    std::vector<int64_t>{1, 2, 4, 8, 16, 32, 64, 128}),
    
    // 64 элемента (степень двойки)
    std::make_tuple(
        []() {
            std::vector<int64_t> v(64);
            for (int i = 0; i < 64; ++i) {
              v[i] = 63 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(64);
            for (int i = 0; i < 64; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // 127 элементов (простое число)
    std::make_tuple(
        []() {
            std::vector<int64_t> v(127);
            for (int i = 0; i < 127; ++i) {
              v[i] = 126 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(127);
            for (int i = 0; i < 127; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // 256 элементов
    std::make_tuple(
        []() {
            std::vector<int64_t> v(256);
            for (int i = 0; i < 256; ++i) {
              v[i] = 255 - i;
            }
            return v;
        }(),
        []() {
            std::vector<int64_t> v(256);
            for (int i = 0; i < 256; ++i) {
              v[i] = i;
            }
            return v;
        }()
    ),
    
    // Случайный набор
    std::make_tuple(std::vector<int64_t>{-10, 5, -3, 8, 0, -7, 2, 9, -1, 4}, 
                    std::vector<int64_t>{-10, -7, -3, -1, 0, 2, 4, 5, 8, 9}),
    
    // Граничные значения
    std::make_tuple(std::vector<int64_t>{INT64_MAX, INT64_MIN, 0, INT64_MAX - 1, INT64_MIN + 1}, 
                    std::vector<int64_t>{INT64_MIN, INT64_MIN + 1, 0, INT64_MAX - 1, INT64_MAX})
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<LeonovaARadixMergeSortSEQ, InType>(kTestParam, PPC_SETTINGS_leonova_a_radix_merge_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = LeonovaARadixMergeSortRunFuncTests::PrintFuncTestName<LeonovaARadixMergeSortRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixMergeSortTests, LeonovaARadixMergeSortRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace leonova_a_radix_merge_sort