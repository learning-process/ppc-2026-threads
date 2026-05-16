  try {
    static_cast<void>(task.PostProcessing());
  } catch (const std::runtime_error &) {
    return true;
  }

  return false;
}

OutType RunTaskTwiceBeforePostProcessing(const InType &input) {
  DoubleSortEvenOddBatcherSTL task(input);
  if (!task.Validation() || !task.PreProcessing() || !task.Run() || !task.Run() || !task.PostProcessing()) {
    throw std::runtime_error("Repeated run pipeline failed");
  }

  return task.GetOutput();
}

OutType RunWithInputMutationAfterPreProcessing(const InType &input, ValueType first_value, ValueType second_value) {
  DoubleSortEvenOddBatcherSTL task(input);
  if (!task.Validation() || !task.PreProcessing()) {
    throw std::runtime_error("Preprocessing pipeline failed");
  }

  auto &input_ref = task.GetInput();
  input_ref[0] = first_value;
  input_ref[1] = second_value;

  if (!task.Run() || !task.PostProcessing()) {
    throw std::runtime_error("Run pipeline failed");
  }

  return task.GetOutput();
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsEmptyInput) {
  CheckMatchesStdSort({});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsSingleElement) {
  CheckMatchesStdSort({42.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsAlreadySortedInput) {
  CheckMatchesStdSort({-7.0, -2.0, -0.0, 0.0, 1.5, 3.0, 4.0, 9.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsReverseSortedInput) {
  CheckMatchesStdSort({9.0, 4.0, 3.0, 1.5, 0.0, -0.0, -2.0, -7.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsOddSizedInput) {
  CheckMatchesStdSort({3.0, -1.0, 2.0, 0.0, 5.0});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsDenseDuplicateInput) {
  CheckMatchesStdSort(MakeDenseDuplicateInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, SortsAllEqualInput) {
  CheckMatchesStdSort(InType(257, 3.5));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForPrimeSizedRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(997, 20260320));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForLargeRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(1024, 20260321));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForExtremesAndSignedZeros) {
  CheckMatchesStdSort(MakeExtremeInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForAlternatingMagnitudeInput) {
  CheckMatchesStdSort(MakeAlternatingMagnitudeInput());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortWhenInputSizeIsLessThanThreadCount) {
  const StlThreadCountGuard guard(8);
  CheckMatchesStdSort({7.0, -4.0, 2.5});
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForOddNumberOfBlocks) {
  const StlThreadCountGuard guard(5);
  CheckMatchesStdSort(GenerateRandomInput(23, 20260322));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForParallelMergeWithOddCarryBlock) {
  const StlThreadCountGuard guard(7);
  CheckMatchesStdSort(GenerateRandomInput(17, 20260325));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, MatchesStdSortForSmallRandomInput) {
  CheckMatchesStdSort(GenerateRandomInput(23, 20260323));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, ValidationRejectsPreparedOutput) {
  EXPECT_TRUE(ValidationRejectsPreparedOutputImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, ThrowsIfPreProcessingBeforeValidation) {
  EXPECT_TRUE(ThrowsIfPreProcessingBeforeValidationImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, ThrowsIfRunBeforePreProcessing) {
  EXPECT_TRUE(ThrowsIfRunBeforePreProcessingImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, ThrowsIfPostProcessingBeforeRun) {
  EXPECT_TRUE(ThrowsIfPostProcessingBeforeRunImpl());
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, AllowsRepeatedRunBeforePostProcessing) {
  const InType input{9.0, -1.0, 5.0, 3.0, -7.0, 11.0, 0.0, 2.0};
  const auto expected_output = RunTaskPipeline(input);
  EXPECT_EQ(RunTaskTwiceBeforePostProcessing(input), expected_output);
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, UsesInputSnapshotFromPreProcessing) {
  EXPECT_EQ(RunWithInputMutationAfterPreProcessing({5.0, 4.0, 3.0, 2.0, 1.0}, -100.0, 200.0),
            (OutType{1.0, 2.0, 3.0, 4.0, 5.0}));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, KeepsElementMultiplicity) {
  const auto input = MakeDenseDuplicateInput();
  const auto output = RunTaskPipeline(input);

  EXPECT_TRUE(std::ranges::is_sorted(output));
  EXPECT_TRUE(std::ranges::is_permutation(output, input));
}

TEST_P(GusevDoubleSortEvenOddBatcherStlEnabled, KeepsOutputEmptyAfterRunningOnEmptyInput) {
  const auto output = RunTaskPipeline({});
  EXPECT_TRUE(output.empty());
}

std::string PrintStlFunctionalParamName(const ::testing::TestParamInfo<int> &info) {
  static_cast<void>(info);
  return "enabled";
}

INSTANTIATE_TEST_SUITE_P(gusev_d_double_sort_even_odd_batcher_stl_enabled, GusevDoubleSortEvenOddBatcherStlEnabled,
                         ::testing::Values(0), PrintStlFunctionalParamName);

}  // namespace
