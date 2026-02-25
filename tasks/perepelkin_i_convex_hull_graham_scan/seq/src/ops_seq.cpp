#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "util/include/util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

PerepelkinIConvexHullGrahamScanSEQ::PerepelkinIConvexHullGrahamScanSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool PerepelkinIConvexHullGrahamScanSEQ::ValidationImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanSEQ::PreProcessingImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanSEQ::RunImpl() {
  return true;
}

bool PerepelkinIConvexHullGrahamScanSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace perepelkin_i_convex_hull_graham_scan
