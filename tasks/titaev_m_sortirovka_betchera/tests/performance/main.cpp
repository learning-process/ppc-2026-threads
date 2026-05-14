#include <gtest/gtest.h>

#include <vector>

#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"
#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {
// Вспомогательная функция для генерации задач производительности
std::shared_ptr<ppc::task::Task<InType, OutType>> CreateOmpTask(const InType &in) {
  auto taskData = std::make_shared<ppc::task::TaskData>();
  // В тестах производительности данные обычно копируются во внутренние векторы
  // Мы создаем TaskData, который будет инициализирован внутри PipelineRun/TaskRun
  return std::make_shared<TitaevSortirovkaBetcheraOMP>(taskData);
}
}  // namespace

TEST(titaev_m_sortirovka_betchera_omp, test_pipeline_run) {
  const int count = 200000;
  InType in(count, 1.0);
  OutType out(count, 0.0);

  auto taskData = std::make_shared<ppc::task::TaskData>();
  taskData->inputs.push_back(reinterpret_cast<uint8_t *>(in.data()));
  taskData->inputs_count.push_back(in.size());
  taskData->outputs.push_back(reinterpret_cast<uint8_t *>(out.data()));
  taskData->outputs_count.push_back(out.size());

  auto task = std::make_shared<TitaevSortirovkaBetcheraOMP>(taskData);
  auto perfAttr = std::make_shared<ppc::core::PerfAttr>();
  perfAttr->num_running = 10;
  const auto t0 = std::chrono::high_resolution_clock::now();
  perfAttr->current_timer = [&] {
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(now - t0).count();
  };

  auto perfResults = std::make_shared<ppc::core::PerfResults>();
  auto perfAnalyzer = std::make_shared<ppc::core::Perf>(task);
  perfAnalyzer->PipelineRun(perfAttr, perfResults);
  ppc::core::Perf::PrintPerfStatistic(perfResults);
  ASSERT_TRUE(std::is_sorted(out.begin(), out.end()));
}

}  // namespace titaev_m_sortirovka_betchera
