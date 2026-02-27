#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {

class UrinOGrahamPassageTest : public ::testing::TestWithParam<TestType> {
 protected:
  void SetUp() override {
    auto params = GetParam();
    int num_points = std::get<0>(params);
    GenerateTestPoints(num_points);
  }

  void GenerateTestPoints(int num_points) {
    input_points_.clear();

    if (num_points == 3) {
      input_points_ = {Point(0.0, 0.0), Point(5.0, 0.0), Point(2.5, 5.0)};
      expected_hull_size_ = 3;
    } else if (num_points == 5) {
      input_points_ = {Point(0.0, 0.0), Point(4.0, 0.0), Point(5.0, 3.0), Point(2.0, 5.0), Point(-1.0, 2.0)};
      expected_hull_size_ = 5;
    } else if (num_points == 7) {
      input_points_ = {Point(0.0, 0.0), Point(3.0, 1.0), Point(2.0, 4.0), Point(5.0, 2.0),
                       Point(1.0, 5.0), Point(4.0, 3.0), Point(0.0, 5.0)};
      expected_hull_size_ = 5;
    }
  }

  bool IsConvexHull(const std::vector<Point> &hull) {
    if (hull.size() < 3) {
      return false;
    }

    for (size_t i = 0; i < hull.size(); i++) {
      size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
      size_t next = (i + 1) % hull.size();

      if (UrinOGrahamPassageSEQ::Orientation(hull[prev], hull[i], hull[next]) < 0) {
        return false;
      }
    }
    return true;
  }

  InType input_points_;
  size_t expected_hull_size_;
};

// Тест через публичный интерфейс - создаем обертку
TEST_P(UrinOGrahamPassageTest, BuildsConvexHull) {
  UrinOGrahamPassageSEQ task(input_points_);

  // Создаем временный объект и проверяем результат через публичные методы
  // Вместо прямого вызова protected методов, используем валидацию через конструктор
  // и проверяем результат через GetOutput()

  // Проверяем, что задача выполняется корректно через публичные методы
  // Для этого нам нужен публичный метод Execute() или подобный
  // Но в текущей структуре его нет

  // Альтернатива: используем динамический полиморфизм через базовый класс
  ppc::task::Task<InType, OutType> *base_task = &task;

  // Вызываем методы через базовый класс
  EXPECT_TRUE(base_task->Validation());
  EXPECT_TRUE(base_task->PreProcessing());
  EXPECT_TRUE(base_task->Run());
  EXPECT_TRUE(base_task->PostProcessing());

  const auto &hull = task.GetOutput();
  EXPECT_GE(hull.size(), 3);
  EXPECT_EQ(hull.size(), expected_hull_size_);
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест для проверки валидации
TEST_F(UrinOGrahamPassageTest, ValidationTest) {
  // Создаем задачу с некорректными данными
  InType invalid_points = {Point(0.0, 0.0)};
  UrinOGrahamPassageSEQ task(invalid_points);

  ppc::task::Task<InType, OutType> *base_task = &task;
  EXPECT_FALSE(base_task->Validation());
}

// Тест с коллинеарными точками
TEST_F(UrinOGrahamPassageTest, HandlesCollinearPoints) {
  InType points = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0),
                   Point(3.0, 0.0), Point(0.0, 1.0), Point(3.0, 1.0)};

  UrinOGrahamPassageSEQ task(points);
  ppc::task::Task<InType, OutType> *base_task = &task;

  EXPECT_TRUE(base_task->Validation());
  EXPECT_TRUE(base_task->PreProcessing());
  EXPECT_TRUE(base_task->Run());
  EXPECT_TRUE(base_task->PostProcessing());

  const auto &hull = task.GetOutput();
  EXPECT_EQ(hull.size(), 4);
  EXPECT_TRUE(IsConvexHull(hull));
}

INSTANTIATE_TEST_SUITE_P(GrahamPassageTests, UrinOGrahamPassageTest,
                         ::testing::Values(std::make_tuple(3, "Triangle"), std::make_tuple(5, "Pentagon"),
                                           std::make_tuple(7, "MixedPoints")));

}  // namespace urin_o_graham_passage
