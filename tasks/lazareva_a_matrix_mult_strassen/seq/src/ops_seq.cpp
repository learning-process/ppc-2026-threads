#include "lazareva_a_matrix_mult_strassen/seq/include/ops_seq.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "lazareva_a_matrix_mult_strassen/common/include/common.hpp"
#include "task/include/task.hpp"

namespace lazareva_a_matrix_mult_strassen {

LazarevaATestTaskSEQ::LazarevaATestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool LazarevaATestTaskSEQ::ValidationImpl() {
  const int n = GetInput().n;
  if (n <= 0) {
    return false;
  }
  const auto expected = static_cast<size_t>(n) * static_cast<size_t>(n);
  return std::cmp_equal(GetInput().a.size(), expected) && std::cmp_equal(GetInput().b.size(), expected);
}

bool LazarevaATestTaskSEQ::PreProcessingImpl() {
  n_ = GetInput().n;
  padded_n_ = NextPowerOfTwo(n_);
  a_ = PadMatrix(GetInput().a, n_, padded_n_);
  b_ = PadMatrix(GetInput().b, n_, padded_n_);
  const auto padded_size = static_cast<size_t>(padded_n_) * static_cast<size_t>(padded_n_);
  result_.assign(padded_size, 0.0);
  return true;
}

bool LazarevaATestTaskSEQ::RunImpl() {
  result_ = Strassen(a_, b_, padded_n_);
  return true;
}

bool LazarevaATestTaskSEQ::PostProcessingImpl() {
  GetOutput() = UnpadMatrix(result_, padded_n_, n_);
  return true;
}

int LazarevaATestTaskSEQ::NextPowerOfTwo(int n) {
  if (n <= 0) {
    return 1;
  }
  int p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

std::vector<double> LazarevaATestTaskSEQ::PadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size, 0.0);
  for (int i = 0; i < old_n; ++i) {
    for (int j = 0; j < old_n; ++j) {
      const auto dst = (static_cast<ptrdiff_t>(i) * new_n) + j;
      const auto src = (static_cast<ptrdiff_t>(i) * old_n) + j;
      result[static_cast<size_t>(dst)] = m[static_cast<size_t>(src)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::UnpadMatrix(const std::vector<double> &m, int old_n, int new_n) {
  const auto new_size = static_cast<size_t>(new_n) * static_cast<size_t>(new_n);
  std::vector<double> result(new_size);
  for (int i = 0; i < new_n; ++i) {
    for (int j = 0; j < new_n; ++j) {
      const auto dst = (static_cast<ptrdiff_t>(i) * new_n) + j;
      const auto src = (static_cast<ptrdiff_t>(i) * old_n) + j;
      result[static_cast<size_t>(dst)] = m[static_cast<size_t>(src)];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::Add(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::Sub(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

void LazarevaATestTaskSEQ::Split(const std::vector<double> &parent, int n, std::vector<double> &a11,
                                 std::vector<double> &a12, std::vector<double> &a21, std::vector<double> &a22) {
  const int half = n / 2;
  const auto half_size = static_cast<size_t>(half) * static_cast<size_t>(half);
  a11.resize(half_size);
  a12.resize(half_size);
  a21.resize(half_size);
  a22.resize(half_size);

  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < half; ++j) {
      const auto idx = static_cast<size_t>((static_cast<ptrdiff_t>(i) * half) + j);
      a11[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j)];
      a12[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j + half)];
      a21[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j)];
      a22[idx] = parent[static_cast<size_t>((static_cast<ptrdiff_t>(i + half) * n) + j + half)];
    }
  }
}

std::vector<double> LazarevaATestTaskSEQ::Merge(const std::vector<double> &c11, const std::vector<double> &c12,
                                                const std::vector<double> &c21, const std::vector<double> &c22, int n) {
  const int full = n * 2;
  const auto full_size = static_cast<size_t>(full) * static_cast<size_t>(full);
  std::vector<double> result(full_size);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const auto src = static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j);
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j)] = c11[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i) * full) + j + n)] = c12[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j)] = c21[src];
      result[static_cast<size_t>((static_cast<ptrdiff_t>(i + n) * full) + j + n)] = c22[src];
    }
  }
  return result;
}

std::vector<double> LazarevaATestTaskSEQ::NaiveMult(const std::vector<double> &a, const std::vector<double> &b, int n) {
  const auto size = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> c(size, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      const double aik = a[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + k)];
      for (int j = 0; j < n; ++j) {
        c[static_cast<size_t>((static_cast<ptrdiff_t>(i) * n) + j)] +=
            aik * b[static_cast<size_t>((static_cast<ptrdiff_t>(k) * n) + j)];
      }
    }
  }
  return c;
}

std::vector<double> LazarevaATestTaskSEQ::Strassen(const std::vector<double> &root_a, const std::vector<double> &root_b,
                                                   int root_n) {
  std::vector<std::vector<double>> results;
  std::vector<StrassenNode> nodes;
  std::vector<int> call_stack;

  results.emplace_back();
  {
    StrassenNode root;
    root.a = root_a;
    root.b = root_b;
    root.n = root_n;
    root.result_slot = 0;
    root.expanded = false;
    nodes.push_back(std::move(root));
  }
  call_stack.push_back(0);

  while (!call_stack.empty()) {
    const int node_idx = call_stack.back();
    call_stack.pop_back();

    const int cur_n = nodes[static_cast<size_t>(node_idx)].n;
    const int cur_slot = nodes[static_cast<size_t>(node_idx)].result_slot;

    if (cur_n <= 64) {
      results[static_cast<size_t>(cur_slot)] =
          NaiveMult(nodes[static_cast<size_t>(node_idx)].a, nodes[static_cast<size_t>(node_idx)].b, cur_n);
      nodes[static_cast<size_t>(node_idx)].a = {};
      nodes[static_cast<size_t>(node_idx)].b = {};
      continue;
    }

    if (!nodes[static_cast<size_t>(node_idx)].expanded) {
      const int half = cur_n / 2;

      std::vector<double> a11;
      std::vector<double> a12;
      std::vector<double> a21;
      std::vector<double> a22;
      std::vector<double> b11;
      std::vector<double> b12;
      std::vector<double> b21;
      std::vector<double> b22;

      Split(nodes[static_cast<size_t>(node_idx)].a, cur_n, a11, a12, a21, a22);
      Split(nodes[static_cast<size_t>(node_idx)].b, cur_n, b11, b12, b21, b22);

      nodes[static_cast<size_t>(node_idx)].a = {};
      nodes[static_cast<size_t>(node_idx)].b = {};
      nodes[static_cast<size_t>(node_idx)].expanded = true;

      std::vector<std::pair<std::vector<double>, std::vector<double>>> args;
      args.reserve(7);
      args.emplace_back(Add(a11, a22, half), Add(b11, b22, half));
      args.emplace_back(Add(a21, a22, half), std::vector<double>(b11));
      args.emplace_back(std::vector<double>(a11), Sub(b12, b22, half));
      args.emplace_back(std::vector<double>(a22), Sub(b21, b11, half));
      args.emplace_back(Add(a11, a12, half), std::vector<double>(b22));
      args.emplace_back(Sub(a21, a11, half), Add(b11, b12, half));
      args.emplace_back(Sub(a12, a22, half), Add(b21, b22, half));

      const int base_slot = static_cast<int>(results.size());
      for (int k = 0; k < 7; ++k) {
        nodes[static_cast<size_t>(node_idx)].child_slots[k] = base_slot + k;
        results.emplace_back();
      }

      call_stack.push_back(node_idx);

      for (int k = 6; k >= 0; --k) {
        StrassenNode child;
        child.a = std::move(args[static_cast<size_t>(k)].first);
        child.b = std::move(args[static_cast<size_t>(k)].second);
        child.n = half;
        child.result_slot = base_slot + k;
        child.expanded = false;
        const int child_idx = static_cast<int>(nodes.size());
        nodes.push_back(std::move(child));
        call_stack.push_back(child_idx);
      }

    } else {
      const int half = cur_n / 2;
      const int *cs = nodes[static_cast<size_t>(node_idx)].child_slots;

      const auto &m1 = results[static_cast<size_t>(cs[0])];
      const auto &m2 = results[static_cast<size_t>(cs[1])];
      const auto &m3 = results[static_cast<size_t>(cs[2])];
      const auto &m4 = results[static_cast<size_t>(cs[3])];
      const auto &m5 = results[static_cast<size_t>(cs[4])];
      const auto &m6 = results[static_cast<size_t>(cs[5])];
      const auto &m7 = results[static_cast<size_t>(cs[6])];

      auto c11 = Add(Sub(Add(m1, m4, half), m5, half), m7, half);
      auto c12 = Add(m3, m5, half);
      auto c21 = Add(m2, m4, half);
      auto c22 = Add(Sub(Add(m1, m3, half), m2, half), m6, half);

      results[static_cast<size_t>(cur_slot)] = Merge(c11, c12, c21, c22, half);
    }
  }

  return results[0];
}

}  // namespace lazareva_a_matrix_mult_strassen
