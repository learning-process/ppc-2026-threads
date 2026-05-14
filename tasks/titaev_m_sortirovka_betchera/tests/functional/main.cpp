#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"

namespace titaev_m_sortirovka_betchera {

TEST(titaev_m_sortirovka_betchera_omp, Test_Sorting_OMP) {
  const int size = 100;
  std::vector<double> in(size);
  std::vector<double> out(size);
  for (int i = 0; i < size; i++) {
    in[i] = size - i;
  }

  auto taskData = std::make_shared<ppc::task::TaskData>();
  taskData->inputs.push_back(reinterpret_cast<uint8_t *>(in.data()));
  taskData->inputs_count.push_back(in.size());
  taskData->outputs.push_back(reinterpret_cast<uint8_t *>(out.data()));
  taskData->outputs_count.push_back(out.size());

  auto task = std::make_shared<TitaevSortirovkaBetcheraOMP>(taskData);
  ASSERT_TRUE(task->ValidationImpl());
  task->PreProcessingImpl();
  task->RunImpl();
  task->PostProcessingImpl();
  ASSERT_TRUE(std::is_sorted(out.begin(), out.end()));
}

}  // namespace titaev_m_sortirovka_betchera
