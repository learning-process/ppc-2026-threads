#include "mityaeva_radix/seq/include/sorter_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace mityaeva_radix {
void SorterSeq::CountingSortAsc(std::vector<double> &source, std::vector<double> &destination, int byte) {
  auto *mas = reinterpret_cast<unsigned char *>(source.data());
  std::vector<int> counter(256, 0);
  auto size = source.size();
  for (std::size_t i = 0; i < size; i++) {
    counter[mas[(8 * i) + byte]]++;
  }
  int j = 0;
  for (; j < 256; j++) {
    if (counter[j] != 0) {
      break;
    }
  }
  int tem = counter[j];
  counter[j] = 0;
  j++;
  for (; j < 256; j++) {
    int b = counter[j];
    counter[j] = tem;
    tem += b;
  }
  for (std::size_t i = 0; i < size; i++) {
    destination[counter[mas[(8 * i) + byte]]] = source[i];
    counter[mas[(8 * i) + byte]]++;
  }
}

void SorterSeq::CountingSortDesc(std::vector<double> &source, std::vector<double> &destination, int byte) {
  auto *mas = reinterpret_cast<unsigned char *>(source.data());
  std::vector<int> count(256, 0);
  auto size = source.size();
  for (std::size_t i = 0; i < size; i++) {
    count[mas[(8 * i) + byte]]++;
  }
  int sum = 0;
  std::vector<int> pos(256, 0);
  for (int i = 255; i >= 0; i--) {
    sum += count[i];
    pos[i] = sum - count[i];
  }
  for (std::size_t i = 0; i < size; i++) {
    int byte_val = mas[(8 * i) + byte];
    destination[pos[byte_val]] = source[i];
    pos[byte_val]++;
  }
}

void SorterSeq::LSDSortDouble(std::vector<double> &inp) {
  if (inp.size() <= 1) {
    return;
  }
  auto count_negative = std::ranges::count_if(inp, [](auto x) { return x < 0; });
  std::vector<double> negative;
  negative.reserve(count_negative);
  std::vector<double> positive;
  positive.reserve(inp.size() - count_negative);
  for (auto x : inp) {
    if (x < 0) {
      negative.push_back(x);
    } else {
      positive.push_back(x);
    }
  }
  std::vector<double> out_n(negative.size());
  std::vector<double> out_p(positive.size());
  for (int i = 0; i < 8; i++) {
    if (i % 2 == 0) {
      CountingSortDesc(negative, out_n, i);
      CountingSortAsc(positive, out_p, i);
    } else {
      CountingSortDesc(out_n, negative, i);
      CountingSortAsc(out_p, positive, i);
    }
  }
  inp.clear();
  inp.insert(inp.begin(), negative.begin(), negative.end());
  inp.insert(inp.end(), positive.begin(), positive.end());
};

}  // namespace mityaeva_radix
