#pragma once

#include <vector>

#include "shkenev_i_constr_hull_for_binary_image_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace shkenev_i_constr_hull_for_binary_image_seq {

class ShkenevIConstrHullSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit ShkenevIConstrHullSeq(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void BinarizeImage(uint8_t threshold = 128);
  void ExtractConnectedRegions();
  static std::vector<ImagePoint> ComputeConvexEnvelope(const std::vector<ImagePoint> &region);
  static long long ComputeCrossProduct(const ImagePoint &p1, const ImagePoint &p2, const ImagePoint &p3);

  BinaryImageData working_image_;
};

}  // namespace shkenev_i_constr_hull_for_binary_image_seq
