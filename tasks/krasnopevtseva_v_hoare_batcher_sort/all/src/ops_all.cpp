#include "krasnopevtseva_v_hoare_batcher_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <queue>
#include <stack>
#include <utility>
#include <vector>
}

bool KrasnopevtsevaVHoareBatcherSortALL::ValidationImpl() {
  return !GetInput().empty();
  const auto &input = GetInput();
  return !input.empty();
}

bool KrasnopevtsevaVHoareBatcherSortALL::PreProcessingImpl() {
                                                      std::vector<int> &sizes, int par_if_greater) {
                                                        int pack = static_cast<int>(pointers.size());
                                                        for (int step = 1; pack > 1; step *= 2, pack /= 2) {
#pragma omp parallel for default(none) shared(pointers, sizes, pack, step, thread_input_size, \
                                                  par_if_greater) if ((thread_input_size / step) > par_if_greater)
                                                          for (int off = 0; off < pack / 2; ++off) {
                                                            auto idx1 = static_cast<std::size_t>(2 * step) *
                                                                        static_cast<std::size_t>(off);
                                                            auto idx2 = idx1 + static_cast<std::size_t>(step);
                                                            BatcherMergeBlocksStep(pointers[idx1], sizes[idx1],
                                                                                   pointers[idx2], sizes[idx2]);
                                                            bool do_parallel =
                                                                (thread_input_size / step) > par_if_greater;

                                                            if (do_parallel) {
#pragma omp parallel for default(none) shared(pointers, sizes, pack, step)
                                                              for (int off = 0; off < pack / 2; ++off) {
                                                                auto idx1 = static_cast<std::size_t>(2 * step) *
                                                                            static_cast<std::size_t>(off);
                                                                auto idx2 = idx1 + static_cast<std::size_t>(step);
                                                                BatcherMergeBlocksStep(pointers[idx1], sizes[idx1],
                                                                                       pointers[idx2], sizes[idx2]);
                                                              }
                                                            } else {
                                                              for (int off = 0; off < pack / 2; ++off) {
                                                                auto idx1 = static_cast<std::size_t>(2 * step) *
                                                                            static_cast<std::size_t>(off);
                                                                auto idx2 = idx1 + static_cast<std::size_t>(step);
                                                                BatcherMergeBlocksStep(pointers[idx1], sizes[idx1],
                                                                                       pointers[idx2], sizes[idx2]);
                                                              }
                                                            }

                                                            if ((pack / 2) - 1 == 0) {
                                                              BatcherMergeBlocksStep(pointers[0],
                                                                                     sizes[sizes.size() - 1],
                                                                                     pointers[pointers.size() - 1],
                                                                                     sizes[sizes.size() - 1]);
                                                            }
                                                          }

                                                          void
                                                          KrasnopevtsevaVHoareBatcherSortALL::ParallelSortChunksOpenMP(
                                                              std::vector<int> & res, int n, int numthreads) {
                                                            void KrasnopevtsevaVHoareBatcherSortALL::ParallelSortChunks(
                                                                std::vector<int> & arr, int n, int numthreads) {
                                                              if (n <= 0) {
                                                                return;
                                                              }

                                                              for (int i = 0; i < numthreads; ++i) {
                                                                std::ptrdiff_t offset =
                                                                    static_cast<std::ptrdiff_t>(i) *
                                                                    static_cast<std::ptrdiff_t>(thread_input_size);
                                                                pointers[i] = res.data() + offset;
                                                                pointers[i] = arr.data() + offset;
                                                                sizes[i] = thread_input_size;
                                                              }
                                                              sizes.back() += thread_input_remainder_size;

#pragma omp parallel for default(none) shared(res, pointers, sizes, numthreads)
#pragma omp parallel for default(none) shared(arr, pointers, sizes, numthreads)
                                                              for (int i = 0; i < numthreads; ++i) {
                                                                int left = static_cast<int>(pointers[i] - res.data());
                                                                int left = static_cast<int>(pointers[i] - arr.data());
                                                                int right = left + sizes[i] - 1;
                                                                if (left < right) {
                                                                  QuickSort(res, left, right);
                                                                  QuickSort(arr, left, right);
                                                                }
                                                              }

                                                              if (n < 1000) {
                                                                QuickSort(data, 0, n - 1);
                                                              } else {
                                                                ParallelSortChunksOpenMP(data, n, numthreads);
                                                                ParallelSortChunks(data, n, numthreads);
                                                              }
                                                            }

                                                            bool KrasnopevtsevaVHoareBatcherSortALL::RunImpl() {
                                                              int rank = 0;
                                                              int proc_size = 1;

                                                              MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                                                              MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

                                                              const auto &input = GetInput();
                                                              int n = static_cast<int>(input.size());
                                                              std::vector<int> result = input;

                                                              if (proc_size == 1 || n <= 10000) {
                                                                std::vector<int> result = input;
                                                                SortLocalData(result);
                                                                if (rank == 0) {
                                                                  GetOutput() = std::move(result);
                                                                }
                                                                if (result.size() <= 1) {
                                                                  GetOutput() = result;
                                                                  return true;
                                                                }

                                                                int local_n = n / proc_size;
                                                                int remainder = n % proc_size;
                                                                std::vector<int> send_counts(proc_size);
                                                                std::vector<int> displs(proc_size);

                                                                for (int i = 0; i < proc_size; ++i) {
                                                                  send_counts[i] = local_n + (i < remainder ? 1 : 0);
                                                                  displs[i] =
                                                                      (i == 0) ? 0 : displs[i - 1] + send_counts[i - 1];
                                                                }

                                                                int my_n = send_counts[rank];
                                                                std::vector<int> local_data(my_n);

                                                                if (rank == 0) {
                                                                  std::vector<int> temp_input = input;
                                                                  MPI_Scatterv(temp_input.data(), send_counts.data(),
                                                                               displs.data(), MPI_INT,
                                                                               local_data.data(), my_n, MPI_INT, 0,
                                                                               MPI_COMM_WORLD);
                                                                } else {
                                                                  MPI_Scatterv(nullptr, send_counts.data(),
                                                                               displs.data(), MPI_INT,
                                                                               local_data.data(), my_n, MPI_INT, 0,
                                                                               MPI_COMM_WORLD);
                                                                }

                                                                if (my_n > 0) {
                                                                  SortLocalData(local_data);
                                                                }

                                                                if (rank == 0) {
                                                                  std::vector<int> global_data(n);
                                                                  MPI_Gatherv(local_data.data(), my_n, MPI_INT,
                                                                              global_data.data(), send_counts.data(),
                                                                              displs.data(), MPI_INT, 0,
                                                                              MPI_COMM_WORLD);

                                                                  std::ranges::sort(global_data);
                                                                  GetOutput() = std::move(global_data);
                                                                } else {
                                                                  MPI_Gatherv(local_data.data(), my_n, MPI_INT, nullptr,
                                                                              nullptr, nullptr, MPI_INT, 0,
                                                                              MPI_COMM_WORLD);
                                                                }

                                                                SortLocalData(result);
                                                                GetOutput() = std::move(result);
                                                                return true;
                                                              }

                                                            }  // namespace krasnopevtseva_v_hoare_batcher_sort
