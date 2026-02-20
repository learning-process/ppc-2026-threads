/*
#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// #include "barkalova_m_mult_matrix_ccs/all/include/ops_all.hpp"
#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
// #include "barkalova_m_mult_matrix_ccs/omp/include/ops_omp.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
// #include "barkalova_m_mult_matrix_ccs/stl/include/ops_stl.hpp"
// #include "barkalova_m_mult_matrix_ccs/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMMultMatrixCcsFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    int width = -1;
    int height = -1;
    int channels = -1;
    std::vector<uint8_t> img;
    // Read image
    {
      std::string abs_path = ppc::util::GetAbsoluteTaskPath(std::string(PPC_ID_barkalova_m_mult_matrix_ccs), "pic.jpg");
      auto *data = stbi_load(abs_path.c_str(), &width, &height, &channels, 0);
      if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + std::string(stbi_failure_reason()));
      }
      img = std::vector<uint8_t>(data, data + (static_cast<ptrdiff_t>(width * height * channels)));
      stbi_image_free(data);
      if (std::cmp_not_equal(width, height)) {
        throw std::runtime_error("width != height: ");
      }
    }

    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = width - height + std::min(std::accumulate(img.begin(), img.end(), 0), channels);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (input_data_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {

TEST_P(BarkalovaMMultMatrixCcsFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BarkalobaMMultMatrixCcsSEQ, InType>(kTestParam, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    BarkalovaMMultMatrixCcsFuncTestsThreads::PrintFuncTestName<BarkalovaMMultMatrixCcsFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, BarkalovaMMultMatrixCcsFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs
*/


/*
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

// Тест с известными значениями (фиксированные матрицы)
class BarkalovaMMultMatrixCcsFixedTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "_" +
           std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_case = std::get<0>(params);
    
    // Создаем матрицы с известными значениями в зависимости от тестового случая
    CCSMatrix A, B;
    
    switch(test_case) {
      case 1: {
        // Тест 1: Простое умножение (2x2) * (2x2)
        // A = [1+2i, 0    ; 0,    3+4i]
        // B = [5+6i, 0    ; 0,    7+8i]
        // Результат: [(1+2i)(5+6i), 0; 0, (3+4i)(7+8i)]
        
        A = CCSMatrix(2, 2);
        A.values = {Complex(1.0, 2.0), Complex(3.0, 4.0)};
        A.row_indices = {0, 1};
        A.col_ptrs = {0, 1, 2};
        A.nnz = 2;
        
        B = CCSMatrix(2, 2);
        B.values = {Complex(5.0, 6.0), Complex(7.0, 8.0)};
        B.row_indices = {0, 1};
        B.col_ptrs = {0, 1, 2};
        B.nnz = 2;
        
        // Ожидаемый результат
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {
          Complex(5.0*1.0 - 6.0*2.0, 1.0*6.0 + 2.0*5.0),  // (1+2i)*(5+6i) = -7+16i
          Complex(7.0*3.0 - 8.0*4.0, 3.0*8.0 + 4.0*7.0)   // (3+4i)*(7+8i) = -11+52i
        };
        expected_result_.row_indices = {0, 1};
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.nnz = 2;
        break;
      }
      
      case 2: {
        // Тест 2: Умножение с разными размерностями (2x3) * (3x2)
        // A = [1+1i, 2+2i, 0    ; 0,    3+3i, 4+4i]
        // B = [5+5i, 0    ; 6+6i, 7+7i ; 0,    8+8i]
        
        A = CCSMatrix(2, 3);
        A.values = {Complex(1.0, 1.0), Complex(2.0, 2.0), Complex(3.0, 3.0), Complex(4.0, 4.0)};
        A.row_indices = {0, 0, 1, 1};
        A.col_ptrs = {0, 2, 3, 4};
        A.nnz = 4;
        
        B = CCSMatrix(3, 2);
        B.values = {Complex(5.0, 5.0), Complex(6.0, 6.0), Complex(7.0, 7.0), Complex(8.0, 8.0)};
        B.row_indices = {0, 1, 1, 2};
        B.col_ptrs = {0, 1, 4};
        B.nnz = 4;
        
        // Ожидаемый результат: C[0][0] = (1+i)*(5+i5) + (2+2i)*(6+6i) = (10i) + (24i) = 34i
        // C[0][1] = (2+2i)*(7+7i) = 28i
        // C[1][0] = (3+3i)*(6+6i) = 36i
        // C[1][1] = (3+3i)*(7+7i) + (4+4i)*(8+8i) = 42i + 64i = 106i
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {
          Complex(0.0, 34.0),   // C[0][0]
          Complex(0.0, 28.0),   // C[0][1]
          Complex(0.0, 36.0),   // C[1][0]
          Complex(0.0, 106.0)   // C[1][1]
        };
        expected_result_.row_indices = {0, 0, 1, 1};
        expected_result_.col_ptrs = {0, 2, 4};
        expected_result_.nnz = 4;
        break;
      }
      
      case 3: {
        // Тест 3: Умножение с нулевым результатом
        // A = [1+2i, 0    ; 3+4i, 0]
        // B = [0,    5+6i; 0,    7+8i]
        
        A = CCSMatrix(2, 2);
        A.values = {Complex(1.0, 2.0), Complex(3.0, 4.0)};
        A.row_indices = {0, 1};
        A.col_ptrs = {0, 2, 2};
        A.nnz = 2;
        
        B = CCSMatrix(2, 2);
        B.values = {Complex(5.0, 6.0), Complex(7.0, 8.0)};
        B.row_indices = {1, 0};
        B.col_ptrs = {0, 0, 2};
        B.nnz = 2;
        
        // Результат должен быть нулевой матрицей (нет пересечений ненулевых элементов)
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.col_ptrs = {0, 0, 0};
        expected_result_.nnz = 0;
        break;
      }
      
      case 4: {
        // Тест 4: Умножение с частичным перекрытием
        // A = [1+1i, 2+2i, 0    ; 0,    3+3i, 4+4i]
        // B = [5+5i, 0    ; 6+6i, 0    ; 7+7i, 8+8i]
        
        A = CCSMatrix(2, 3);
        A.values = {Complex(1.0, 1.0), Complex(2.0, 2.0), Complex(3.0, 3.0), Complex(4.0, 4.0)};
        A.row_indices = {0, 0, 1, 1};
        A.col_ptrs = {0, 2, 3, 4};
        A.nnz = 4;
        
        B = CCSMatrix(3, 2);
        B.values = {Complex(5.0, 5.0), Complex(6.0, 6.0), Complex(7.0, 7.0), Complex(8.0, 8.0)};
        B.row_indices = {0, 1, 2, 2};
        B.col_ptrs = {0, 2, 4};
        B.nnz = 4;
        
        // Ожидаемый результат с частичным перекрытием
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {
          Complex(0.0, 34.0),   // C[0][0]: (1+i)*(5+5i) + (2+2i)*(6+6i)
          Complex(0.0, 28.0),   // C[0][1]: (2+2i)*(7+7i)
          Complex(0.0, 36.0),   // C[1][0]: (3+3i)*(6+6i)
          Complex(0.0, 70.0)    // C[1][1]: (3+3i)*(7+7i) + (4+4i)*(8+8i)
        };
        expected_result_.row_indices = {0, 0, 1, 1};
        expected_result_.col_ptrs = {0, 2, 4};
        expected_result_.nnz = 4;
        break;
      }
      
      default:
        throw std::runtime_error("Unknown test case");
    }
    
    input_data_ = std::make_pair(A, B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем размерность
    if (output_data.rows != expected_result_.rows || output_data.cols != expected_result_.cols) {
      return false;
    }
    
    // Проверяем количество ненулевых элементов
    if (output_data.nnz != expected_result_.nnz) {
      return false;
    }
    
    // Если результат нулевой, проверяем пустые векторы
    if (expected_result_.nnz == 0) {
      return output_data.values.empty() && output_data.row_indices.empty();
    }
    
    // Проверяем структуру col_ptrs
    if (output_data.col_ptrs.size() != expected_result_.col_ptrs.size()) {
      return false;
    }
    
    for (size_t i = 0; i < output_data.col_ptrs.size(); ++i) {
      if (output_data.col_ptrs[i] != expected_result_.col_ptrs[i]) {
        return false;
      }
    }
    
    // Проверяем индексы строк
    if (output_data.row_indices.size() != expected_result_.row_indices.size()) {
      return false;
    }
    
    for (size_t i = 0; i < output_data.row_indices.size(); ++i) {
      if (output_data.row_indices[i] != expected_result_.row_indices[i]) {
        return false;
      }
    }
    
    // Проверяем значения с учетом погрешности
    const double eps = 1e-10;
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i].real() - expected_result_.values[i].real()) > eps ||
          std::abs(output_data.values[i].imag() - expected_result_.values[i].imag()) > eps) {
        return false;
      }
    }
    
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CCSMatrix expected_result_;
};

namespace {

// Параметризованный тест с известными значениями
TEST_P(BarkalovaMMultMatrixCcsFixedTest, FixedMatrixMultiplication) {
  ExecuteTest(GetParam());
}

// Тестовые параметры: (test_case, placeholder1, placeholder2)
const std::array<TestType, 4> kFixedTestParams = {
    std::make_tuple(1, 0, 0),  // диагональные матрицы
    std::make_tuple(2, 0, 0),  // прямоугольные матрицы
    std::make_tuple(3, 0, 0),  // нулевой результат
    std::make_tuple(4, 0, 0)   // частичное перекрытие
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BarkalovaMMultMatrixCcsSEQ, InType>(kFixedTestParams, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFixedTestName = BarkalovaMMultMatrixCcsFixedTest::PrintFuncTestName<BarkalovaMMultMatrixCcsFixedTest>;

INSTANTIATE_TEST_SUITE_P(FixedMatrixTests, BarkalovaMMultMatrixCcsFixedTest, kGtestValues, kFixedTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs
*/


#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMMultMatrixCcsFixedTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + 
           std::to_string(std::get<1>(test_param)) + "_" +
           std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_case = std::get<0>(params);
    
    CCSMatrix A, B;
    
    switch(test_case) {
      case 1: {
        // Тест 1: Простейшие диагональные матрицы 2x2
        A = CCSMatrix(2, 2);
        A.values = {Complex(1.0, 0.0), Complex(1.0, 0.0)};
        A.row_indices = {0, 1};
        A.col_ptrs = {0, 1, 2};
        A.nnz = 2;
        
        B = CCSMatrix(2, 2);
        B.values = {Complex(2.0, 0.0), Complex(2.0, 0.0)};
        B.row_indices = {0, 1};
        B.col_ptrs = {0, 1, 2};
        B.nnz = 2;
        
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {Complex(2.0, 0.0), Complex(2.0, 0.0)};
        expected_result_.row_indices = {0, 1};
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.nnz = 2;
        break;
      }
      
      case 2: {
        // Тест 2: Матрицы 2x2 с разными значениями
        // A = [1, 2; 3, 4] в CCS формате
        A = CCSMatrix(2, 2);
        A.values = {Complex(1.0, 0.0), Complex(3.0, 0.0),    // столбец 0: строки 0,1
                    Complex(2.0, 0.0), Complex(4.0, 0.0)};   // столбец 1: строки 0,1
        A.row_indices = {0, 1, 0, 1};
        A.col_ptrs = {0, 2, 4};
        A.nnz = 4;
        
        // B = [5, 6; 7, 8] в CCS формате
        B = CCSMatrix(2, 2);
        B.values = {Complex(5.0, 0.0), Complex(7.0, 0.0),    // столбец 0: строки 0,1
                    Complex(6.0, 0.0), Complex(8.0, 0.0)};   // столбец 1: строки 0,1
        B.row_indices = {0, 1, 0, 1};
        B.col_ptrs = {0, 2, 4};
        B.nnz = 4;
        
        // Результат: C = A * B
        // C[0][0] = 1*5 + 2*7 = 5 + 14 = 19
        // C[0][1] = 1*6 + 2*8 = 6 + 16 = 22
        // C[1][0] = 3*5 + 4*7 = 15 + 28 = 43
        // C[1][1] = 3*6 + 4*8 = 18 + 32 = 50
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {Complex(19.0, 0.0), Complex(43.0, 0.0),    // столбец 0
                                   Complex(22.0, 0.0), Complex(50.0, 0.0)};   // столбец 1
        expected_result_.row_indices = {0, 1, 0, 1};
        expected_result_.col_ptrs = {0, 2, 4};
        expected_result_.nnz = 4;
        break;
      }
      
      case 3: {
        // Тест 3: Умножение с нулевым результатом
        // A = [1, 0; 2, 0] - ненулевые только в первом столбце
        A = CCSMatrix(2, 2);
        A.values = {Complex(1.0, 0.0), Complex(2.0, 0.0)};
        A.row_indices = {0, 1};
        A.col_ptrs = {0, 2, 2};  // первый столбец: 2 элемента, второй: 0
        A.nnz = 2;
        
        // B = [0, 3; 0, 4] - ненулевые только во втором столбце
        B = CCSMatrix(2, 2);
        B.values = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
        B.row_indices = {0, 1};  // оба элемента во втором столбце
        B.col_ptrs = {0, 0, 2};  // первый столбец: 0, второй: 2
        B.nnz = 2;
        
        // Результат должен быть нулевой матрицей
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.col_ptrs = {0, 0, 0};
        expected_result_.values.clear();
        expected_result_.row_indices.clear();
        expected_result_.nnz = 0;
        break;
      }
      
      case 4: {
        // Тест 4: Прямоугольные матрицы 2x3 * 3x2
        // A: 2x3
        // столбец 0: [0]=1, [1]=4
        // столбец 1: [0]=2, [1]=5
        // столбец 2: [0]=3, [1]=6
        A = CCSMatrix(2, 3);
        A.values = {Complex(1.0, 0.0), Complex(4.0, 0.0),    // столбец 0
                    Complex(2.0, 0.0), Complex(5.0, 0.0),    // столбец 1
                    Complex(3.0, 0.0), Complex(6.0, 0.0)};   // столбец 2
        A.row_indices = {0, 1, 0, 1, 0, 1};
        A.col_ptrs = {0, 2, 4, 6};
        A.nnz = 6;
        
        // B: 3x2
        // столбец 0: [0]=7, [1]=9, [2]=11
        // столбец 1: [0]=8, [1]=10, [2]=12
        B = CCSMatrix(3, 2);
        B.values = {Complex(7.0, 0.0), Complex(9.0, 0.0), Complex(11.0, 0.0),   // столбец 0
                    Complex(8.0, 0.0), Complex(10.0, 0.0), Complex(12.0, 0.0)}; // столбец 1
        B.row_indices = {0, 1, 2, 0, 1, 2};
        B.col_ptrs = {0, 3, 6};
        B.nnz = 6;
        
        // Результат: C = A * B (2x2)
        // C[0][0] = 1*7 + 2*9 + 3*11 = 7 + 18 + 33 = 58
        // C[0][1] = 1*8 + 2*10 + 3*12 = 8 + 20 + 36 = 64
        // C[1][0] = 4*7 + 5*9 + 6*11 = 28 + 45 + 66 = 139
        // C[1][1] = 4*8 + 5*10 + 6*12 = 32 + 50 + 72 = 154
        expected_result_ = CCSMatrix(2, 2);
        expected_result_.values = {Complex(58.0, 0.0), Complex(139.0, 0.0),   // столбец 0
                                   Complex(64.0, 0.0), Complex(154.0, 0.0)};  // столбец 1
        expected_result_.row_indices = {0, 1, 0, 1};
        expected_result_.col_ptrs = {0, 2, 4};
        expected_result_.nnz = 4;
        break;
      }
      
      default:
        throw std::runtime_error("Unknown test case");
    }
    
    input_data_ = std::make_pair(A, B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Проверяем размерность
    if (output_data.rows != expected_result_.rows || output_data.cols != expected_result_.cols) {
      std::cout << "Dimension mismatch: got (" << output_data.rows << "x" << output_data.cols 
                << "), expected (" << expected_result_.rows << "x" << expected_result_.cols << ")" << std::endl;
      return false;
    }
    
    // Проверяем количество ненулевых элементов
    if (output_data.nnz != expected_result_.nnz) {
      std::cout << "NNZ mismatch: got " << output_data.nnz << ", expected " << expected_result_.nnz << std::endl;
      return false;
    }
    
    // Если результат нулевой, проверяем пустые векторы
    if (expected_result_.nnz == 0) {
      return output_data.values.empty() && output_data.row_indices.empty();
    }
    
    // Проверяем структуру col_ptrs
    if (output_data.col_ptrs.size() != expected_result_.col_ptrs.size()) {
      std::cout << "col_ptrs size mismatch" << std::endl;
      return false;
    }
    
    for (size_t i = 0; i < output_data.col_ptrs.size(); ++i) {
      if (output_data.col_ptrs[i] != expected_result_.col_ptrs[i]) {
        std::cout << "col_ptrs[" << i << "] mismatch: got " << output_data.col_ptrs[i] 
                  << ", expected " << expected_result_.col_ptrs[i] << std::endl;
        return false;
      }
    }
    
    // Проверяем индексы строк
    if (output_data.row_indices.size() != expected_result_.row_indices.size()) {
      std::cout << "row_indices size mismatch" << std::endl;
      return false;
    }
    
    for (size_t i = 0; i < output_data.row_indices.size(); ++i) {
      if (output_data.row_indices[i] != expected_result_.row_indices[i]) {
        std::cout << "row_indices[" << i << "] mismatch: got " << output_data.row_indices[i] 
                  << ", expected " << expected_result_.row_indices[i] << std::endl;
        return false;
      }
    }
    
    // Проверяем значения с учетом погрешности
    const double eps = 1e-10;
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i].real() - expected_result_.values[i].real()) > eps ||
          std::abs(output_data.values[i].imag() - expected_result_.values[i].imag()) > eps) {
        std::cout << "Value[" << i << "] mismatch: got (" << output_data.values[i].real() 
                  << "," << output_data.values[i].imag() << "), expected (" 
                  << expected_result_.values[i].real() << "," << expected_result_.values[i].imag() << ")" << std::endl;
        return false;
      }
    }
    
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CCSMatrix expected_result_;
};

namespace {

TEST_P(BarkalovaMMultMatrixCcsFixedTest, FixedMatrixMultiplication) {
  ExecuteTest(GetParam());
}

// Тестовые параметры: (test_case, placeholder1, placeholder2)
const std::array<TestType, 4> kFixedTestParams = {
    std::make_tuple(1, 0, 0),  // простейшие диагональные
    std::make_tuple(2, 0, 0),  // плотные 2x2
    std::make_tuple(3, 0, 0),  // нулевой результат
    std::make_tuple(4, 0, 0)   // прямоугольные
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BarkalovaMMultMatrixCcsSEQ, InType>(kFixedTestParams, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFixedTestName = BarkalovaMMultMatrixCcsFixedTest::PrintFuncTestName<BarkalovaMMultMatrixCcsFixedTest>;

INSTANTIATE_TEST_SUITE_P(FixedMatrixTests, BarkalovaMMultMatrixCcsFixedTest, kGtestValues, kFixedTestName);

}  // namespace

}  // namespace barkalova_m_mult_matrix_ccs