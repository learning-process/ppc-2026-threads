#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "sabutay_sparse_complex_ccs_multfix/common/include/common.hpp"
#include "sabutay_sparse_complex_ccs_multfix/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_multfix {
namespace {

using Z = std::complex<double>;
constexpr double kCmpTol = 1e-12;

struct Dense {
  int m{0};
  int n{0};
  std::vector<Z> v;

  static Dense New(int rows, int cols) {
    return {.m = rows, .n = cols, .v = std::vector<Z>(static_cast<std::size_t>(rows * cols), Z{})};
  }

  Z &At(int row, int col) {
    return v[((static_cast<std::size_t>(row) * static_cast<std::size_t>(n)) + static_cast<std::size_t>(col))];
  }
  [[nodiscard]] const Z &At(int row, int col) const {
    return v[((static_cast<std::size_t>(row) * static_cast<std::size_t>(n)) + static_cast<std::size_t>(col))];
  }
};

Dense Matmul(const Dense &a_mat, const Dense &b_mat) {
  const int m = a_mat.m;
  const int k = a_mat.n;
  const int p = b_mat.n;
  Dense c_res = Dense::New(m, p);
  for (int row = 0; row < m; ++row) {
    for (int tcol = 0; tcol < p; ++tcol) {
      Z acc{};
      for (int inner = 0; inner < k; ++inner) {
        acc += a_mat.At(row, inner) * b_mat.At(inner, tcol);
      }
      c_res.At(row, tcol) = acc;
    }
  }
  return c_res;
}

CCS DenseToCcs(const Dense &dense, double drop) {
  CCS g;
  g.row_count = dense.m;
  g.col_count = dense.n;
  g.col_start.assign(static_cast<std::size_t>(g.col_count) + 1U, 0);
  g.row_index.clear();
  g.nz.clear();
  for (int jcol = 0; jcol < dense.n; ++jcol) {
    for (int irow = 0; irow < dense.m; ++irow) {
      if (std::abs(dense.At(irow, jcol)) > drop) {
        g.row_index.push_back(irow);
        g.nz.push_back(dense.At(irow, jcol));
      }
    }
    g.col_start[static_cast<std::size_t>(jcol) + 1U] = static_cast<int>(g.nz.size());
  }
  return g;
}

bool DenseClose(const Z &a, const Z &b) {
  return (std::abs(a.real() - b.real()) < kCmpTol) && (std::abs(a.imag() - b.imag()) < kCmpTol);
}

bool BuildUsedMaskFromSparseColumn(const CCS &got, const Dense &ref, int jcol, int col_begin, int col_end,
                                   std::vector<int> *used) {
  for (int pidx = col_begin; pidx < col_end; ++pidx) {
    const int irow = got.row_index[static_cast<std::size_t>(pidx)];
    if (irow < 0 || irow >= ref.m) {
      return false;
    }
    if (std::abs(ref.At(irow, jcol)) <= kCmpTol) {
      return false;
    }
    if (used->at(static_cast<std::size_t>(irow)) != 0) {
      return false;
    }
    (*used)[static_cast<std::size_t>(irow)] = 1;
    if (!DenseClose(got.nz[static_cast<std::size_t>(pidx)], ref.At(irow, jcol))) {
      return false;
    }
  }
  return true;
}

bool DenseColumnMatchesUsedMask(const Dense &ref, int jcol, const std::vector<int> &used) {
  for (int irow = 0; irow < ref.m; ++irow) {
    const bool has_ref = std::abs(ref.At(irow, jcol)) > kCmpTol;
    if (has_ref) {
      if (used[static_cast<std::size_t>(irow)] == 0) {
        return false;
      }
    } else {
      if (used[static_cast<std::size_t>(irow)] != 0) {
        return false;
      }
    }
  }
  return true;
}

bool ColumnNonzerosMatchDenseView(const CCS &got, const Dense &ref, int jcol) {
  const int b0 = got.col_start[static_cast<std::size_t>(jcol)];
  const int b1 = got.col_start[static_cast<std::size_t>(jcol) + 1U];
  std::vector<int> used(static_cast<std::size_t>(ref.m), 0);
  if (!BuildUsedMaskFromSparseColumn(got, ref, jcol, b0, b1, &used)) {
    return false;
  }
  return DenseColumnMatchesUsedMask(ref, jcol, used);
}

bool CcsEqualsDenseView(const CCS &got, const Dense &ref) {
  if (got.row_count != ref.m || got.col_count != ref.n) {
    return false;
  }
  for (int jcol = 0; jcol < ref.n; ++jcol) {
    if (!ColumnNonzerosMatchDenseView(got, ref, jcol)) {
      return false;
    }
  }
  return true;
}

}  // namespace

class SabutayARunFuncTestsThreadsFIX : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static auto PrintTestParam(const TestType &id) -> std::string {
    return std::to_string(id);
  }

 protected:
  void SetUp() override {
    const auto case_id = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const double drop = 1e-15;
    if (case_id == 0) {
      // 3x3 and 3x2 pattern; nontrivial off-diagonals, purely real factors.
      Dense a = Dense::New(3, 3);
      a.At(0, 0) = Z(1, 0);
      a.At(0, 1) = Z(2, 0);
      a.At(1, 2) = Z(-1, 0);
      a.At(2, 0) = Z(0.5, 0);
      Dense b = Dense::New(3, 2);
      b.At(0, 0) = Z(3, 0);
      b.At(1, 1) = Z(1, 0);
      b.At(2, 0) = Z(0, 0);
      b.At(2, 1) = Z(2, 0);
      a_in_ = DenseToCcs(a, drop);
      b_in_ = DenseToCcs(b, drop);
      expected_dense_ = Matmul(a, b);
    } else if (case_id == 1) {
      // Complex phases; 2x2 times 2x2.
      Dense a = Dense::New(2, 2);
      a.At(0, 0) = Z(0, 1);
      a.At(1, 0) = Z(2, -0.5);
      a.At(0, 1) = Z(1, 1);
      Dense b = Dense::New(2, 2);
      b.At(0, 0) = Z(0.5, 0.25);
      b.At(1, 0) = Z(1, 0);
      b.At(0, 1) = Z(0, -1);
      b.At(1, 1) = Z(3, 0);
      a_in_ = DenseToCcs(a, drop);
      b_in_ = DenseToCcs(b, drop);
      expected_dense_ = Matmul(a, b);
    } else {
      // Single nonzeros on a chain; 4x1 * 1x3 -> 4x3 rank-one.
      Dense a = Dense::New(4, 1);
      a.At(1, 0) = Z(1, 1);
      Dense b = Dense::New(1, 3);
      b.At(0, 0) = Z(1, 0);
      b.At(0, 1) = Z(0, 2);
      b.At(0, 2) = Z(-0.5, 0.25);
      a_in_ = DenseToCcs(a, drop);
      b_in_ = DenseToCcs(b, drop);
      expected_dense_ = Matmul(a, b);
    }
  }

  bool CheckTestOutputData(OutType &output) override {
    return CcsEqualsDenseView(output, expected_dense_);
  }

  InType GetTestInputData() override {
    return std::make_tuple(a_in_, b_in_);
  }

 private:
  CCS a_in_;
  CCS b_in_;
  Dense expected_dense_;
};

namespace {

TEST_P(SabutayARunFuncTestsThreadsFIX, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam{0, 1, 2};

const auto kTestTasksList = ppc::util::AddFuncTask<SabutaySparseComplexCcsMultFixSEQ, InType>(
    kTestParam, PPC_SETTINGS_sabutay_sparse_complex_ccs_multfix);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SabutayARunFuncTestsThreadsFIX::PrintFuncTestName<SabutayARunFuncTestsThreadsFIX>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, SabutayARunFuncTestsThreadsFIX, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabutay_sparse_complex_ccs_multfix
