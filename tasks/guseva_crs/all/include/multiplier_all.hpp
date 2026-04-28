#pragma once
#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/multiplier.hpp"

namespace guseva_crs {

class MultiplierAll : public Multiplier {
  static void PerformCalculation(std::size_t ind3, std::size_t ind4, const CRS &a, const CRS &bt, double &sum,
                                 const std::vector<int> &temp) {
    for (std::size_t k = ind3; k < ind4; k++) {
      std::size_t bcol = bt.cols[k];
      int aind = temp[bcol];

      if (aind != -1) {
        sum += a.values[aind] * bt.values[k];
      }
    }
  }

  static void ComputeLocalRow(std::size_t global_i, std::size_t n, const CRS &a, const CRS &bt,
                              std::vector<std::size_t> &columns, std::vector<double> &values, std::size_t &row_nnz) {
    std::vector<int> temp(n, -1);

    std::size_t ind1 = a.row_ptrs[global_i];
    std::size_t ind2 = a.row_ptrs[global_i + 1];
    for (std::size_t j = ind1; j < ind2; j++) {
      std::size_t col = a.cols[j];
      temp[col] = static_cast<int>(j);
    }

    for (std::size_t j = 0; j < n; j++) {
      double sum = 0;
      std::size_t ind3 = bt.row_ptrs[j];
      std::size_t ind4 = bt.row_ptrs[j + 1];

      PerformCalculation(ind3, ind4, a, bt, sum, temp);

      if (std::fabs(sum) > kZERO) {
        columns.push_back(j);
        values.push_back(sum);
        row_nnz++;
      }
    }
  }

  static void ComputeLocalResults(std::size_t start_row, std::size_t local_nrows, std::size_t n, const CRS &a,
                                  const CRS &bt, std::vector<std::vector<std::size_t>> &local_columns,
                                  std::vector<std::vector<double>> &local_values,
                                  std::vector<std::size_t> &local_row_index) {
#pragma omp parallel for default(none) \
    shared(n, a, bt, local_columns, local_values, local_row_index, start_row, local_nrows)
    for (std::size_t local_i = 0; local_i < local_nrows; local_i++) {
      std::size_t global_i = start_row + local_i;
      ComputeLocalRow(global_i, n, a, bt, local_columns[local_i], local_values[local_i], local_row_index[local_i]);
    }
  }

  static void FlattenLocalData(const std::vector<std::vector<std::size_t>> &local_columns,
                               const std::vector<std::vector<double>> &local_values,
                               std::vector<std::size_t> &flat_columns, std::vector<double> &flat_values,
                               std::vector<int> &row_sizes) {
    for (std::size_t i = 0; i < local_columns.size(); i++) {
      row_sizes[i] = static_cast<int>(local_columns[i].size());
      flat_columns.insert(flat_columns.end(), local_columns[i].begin(), local_columns[i].end());
      flat_values.insert(flat_values.end(), local_values[i].begin(), local_values[i].end());
    }
  }

  struct ProcessData {
    std::size_t start_row{};
    std::size_t local_nrows{};
    std::vector<int> row_sizes;
    std::vector<std::size_t> flat_columns;
    std::vector<double> flat_values;
  };

  static ProcessData ReceiveProcessData(int source, std::size_t p_start_row, std::size_t p_local_nrows) {
    ProcessData data;
    data.start_row = p_start_row;
    data.local_nrows = p_local_nrows;
    data.row_sizes.resize(p_local_nrows);

    MPI_Recv(data.row_sizes.data(), static_cast<int>(p_local_nrows), MPI_INT, source, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    int total_nz = 0;
    for (std::size_t i = 0; i < p_local_nrows; i++) {
      total_nz += data.row_sizes[i];
    }

    if (total_nz > 0) {
      data.flat_columns.resize(total_nz);
      data.flat_values.resize(total_nz);
      MPI_Recv(data.flat_columns.data(), total_nz, MPI_UNSIGNED_LONG, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(data.flat_values.data(), total_nz, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return data;
  }

  static void AssembleResultFromProcess(const ProcessData &data, std::vector<std::vector<std::size_t>> &columns,
                                        std::vector<std::vector<double>> &values) {
    std::size_t offset = 0;
    for (std::size_t local_i = 0; local_i < data.local_nrows; local_i++) {
      std::size_t global_row = data.start_row + local_i;
      int row_size = data.row_sizes[local_i];

      if (row_size > 0) {
        columns[global_row].resize(row_size);
        values[global_row].resize(row_size);

        for (int j = 0; j < row_size; j++) {
          columns[global_row][j] = data.flat_columns[offset + j];
          values[global_row][j] = data.flat_values[offset + j];
        }
        offset += static_cast<std::size_t>(row_size);
      }
    }
  }

  static CRS BuildFinalMatrix(std::size_t n, std::vector<std::vector<std::size_t>> &columns,
                              std::vector<std::vector<double>> &values) {
    CRS result;
    result.row_ptrs.resize(n + 1, 0);

    std::size_t nz = 0;
    for (std::size_t i = 0; i < n; i++) {
      result.row_ptrs[i] = nz;
      nz += columns[i].size();
    }
    result.row_ptrs[n] = nz;

    result.cols.reserve(nz);
    result.values.reserve(nz);
    for (std::size_t i = 0; i < n; i++) {
      result.cols.insert(result.cols.end(), columns[i].begin(), columns[i].end());
      result.values.insert(result.values.end(), values[i].begin(), values[i].end());
    }

    result.nz = nz;
    result.ncols = n;
    result.nrows = n;

    return result;
  }

  static void SendLocalData(int dest, const std::vector<int> &row_sizes, const std::vector<std::size_t> &flat_columns,
                            const std::vector<double> &flat_values) {
    std::vector<int> row_sizes_copy = row_sizes;
    MPI_Send(row_sizes_copy.data(), static_cast<int>(row_sizes_copy.size()), MPI_INT, dest, 0, MPI_COMM_WORLD);

    if (!flat_columns.empty()) {
      std::vector<std::size_t> flat_columns_copy = flat_columns;
      std::vector<double> flat_values_copy = flat_values;
      MPI_Send(flat_columns_copy.data(), static_cast<int>(flat_columns_copy.size()), MPI_UNSIGNED_LONG, dest, 1,
               MPI_COMM_WORLD);
      MPI_Send(flat_values_copy.data(), static_cast<int>(flat_values_copy.size()), MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
    }
  }

 public:
  [[nodiscard]] CRS Multiply(const CRS &a, const CRS &b) const override {
    int rank = -1;
    int num_procs = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if ((rank < 0) || (num_procs < 0)) {
      throw std::runtime_error("MPI rank or world size is incorrect");
    }

    std::size_t n = a.nrows;

    if (num_procs == 1) {
      return MultiplySerial(a, b);
    }

    std::size_t rows_per_proc = n / static_cast<std::size_t>(num_procs);
    std::size_t remainder = n % static_cast<std::size_t>(num_procs);
    std::size_t start_row =
        (static_cast<std::size_t>(rank) * rows_per_proc) + std::min(static_cast<std::size_t>(rank), remainder);
    std::size_t local_nrows = rows_per_proc + (std::cmp_less(rank, remainder) ? 1 : 0);

    auto bt = this->Transpose(b);

    std::vector<std::vector<std::size_t>> local_columns(local_nrows);
    std::vector<std::vector<double>> local_values(local_nrows);
    std::vector<std::size_t> local_row_index(local_nrows, 0);

    ComputeLocalResults(start_row, local_nrows, n, a, bt, local_columns, local_values, local_row_index);

    std::vector<std::size_t> flat_columns;
    std::vector<double> flat_values;
    std::vector<int> row_sizes(local_nrows);
    FlattenLocalData(local_columns, local_values, flat_columns, flat_values, row_sizes);

    CRS result;

    if (rank == 0) {
      std::vector<std::vector<std::size_t>> columns(n);
      std::vector<std::vector<double>> values(n);

      for (int pp = 0; pp < num_procs; pp++) {
        std::size_t p_start_row =
            (static_cast<std::size_t>(pp) * rows_per_proc) + std::min(static_cast<std::size_t>(pp), remainder);
        std::size_t p_local_nrows = rows_per_proc + (std::cmp_less(pp, remainder) ? 1 : 0);

        if (pp == 0) {
          ProcessData self_data;
          self_data.start_row = p_start_row;
          self_data.local_nrows = p_local_nrows;
          self_data.row_sizes = row_sizes;
          self_data.flat_columns = flat_columns;
          self_data.flat_values = flat_values;
          AssembleResultFromProcess(self_data, columns, values);
        } else {
          ProcessData received_data = ReceiveProcessData(pp, p_start_row, p_local_nrows);
          AssembleResultFromProcess(received_data, columns, values);
        }
      }

      result = BuildFinalMatrix(n, columns, values);
    } else {
      SendLocalData(0, row_sizes, flat_columns, flat_values);
    }

    BroadcastResult(result, rank);

    return result;
  }

  [[nodiscard]] CRS MultiplySerial(const CRS &a, const CRS &b) const {
    std::size_t n = a.nrows;
    auto bt = this->Transpose(b);

    std::vector<std::vector<std::size_t>> columns(n);
    std::vector<std::vector<double>> values(n);
    std::vector<std::size_t> row_index(n + 1, 0);

#pragma omp parallel for default(none) shared(n, a, bt, columns, values, row_index)
    for (std::size_t i = 0; i < n; i++) {
      std::vector<int> temp(n, -1);

      std::size_t ind1 = a.row_ptrs[i];
      std::size_t ind2 = a.row_ptrs[i + 1];
      for (std::size_t j = ind1; j < ind2; j++) {
        std::size_t col = a.cols[j];
        temp[col] = static_cast<int>(j);
      }

      for (std::size_t j = 0; j < n; j++) {
        double sum = 0;
        std::size_t ind3 = bt.row_ptrs[j];
        std::size_t ind4 = bt.row_ptrs[j + 1];

        PerformCalculation(ind3, ind4, a, bt, sum, temp);

        if (std::fabs(sum) > kZERO) {
          columns[i].push_back(j);
          values[i].push_back(sum);
          row_index[i]++;
        }
      }
    }

    std::size_t nz = 0;
    for (std::size_t i = 0; i < n; i++) {
      std::size_t tmp = row_index[i];
      row_index[i] = nz;
      nz += tmp;
    }
    row_index[n] = nz;

    CRS result;
    result.cols.reserve(nz);
    result.values.reserve(nz);
    for (std::size_t i = 0; i < n; i++) {
      result.cols.insert(result.cols.end(), columns[i].begin(), columns[i].end());
      result.values.insert(result.values.end(), values[i].begin(), values[i].end());
    }
    result.row_ptrs = row_index;
    result.nz = nz;
    result.ncols = n;
    result.nrows = n;

    return result;
  }

  static void BroadcastResult(CRS &result, int rank) {
    std::size_t nrows = result.nrows;
    std::size_t ncols = result.ncols;
    std::size_t nz = result.nz;

    MPI_Bcast(&nrows, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncols, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank != 0) {
      result.nrows = nrows;
      result.ncols = ncols;
      result.nz = nz;
      if (nrows > 0) {
        result.row_ptrs.resize(nrows + 1);
      }
      if (nz > 0) {
        result.cols.resize(nz);
        result.values.resize(nz);
      }
    }

    if (nrows > 0) {
      MPI_Bcast(result.row_ptrs.data(), static_cast<int>(nrows + 1), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    }

    if (nz > 0) {
      MPI_Bcast(result.cols.data(), static_cast<int>(nz), MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(result.values.data(), static_cast<int>(nz), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }
};

}  // namespace guseva_crs
