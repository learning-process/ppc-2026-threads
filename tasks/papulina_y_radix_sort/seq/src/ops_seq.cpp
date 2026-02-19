#include "papulina_y_radix_sort/seq/include/ops_seq.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"

namespace papulina_y_radix_sort {

PapulinaYRadixSortSEQ::PapulinaYRadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}
bool PapulinaYRadixSortSEQ::ValidationImpl() {
  return true;
}

bool PapulinaYRadixSortSEQ::PreProcessingImpl() {
  return true;
}

bool PapulinaYRadixSortSEQ::RunImpl() {
  double *result = GetInput().data();
  size_t size = GetInput().size();

  RadixSort(result, static_cast<int>(size));

  GetOutput() = std::vector<double>(size);
  for (size_t i = 0; i < size; i++) {
    GetOutput()[i] = result[i];
  }

  return true;
}

bool PapulinaYRadixSortSEQ::PostProcessingImpl() {
  return true;
}
void PapulinaYRadixSortSEQ::SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size) {
  auto *byte_view = reinterpret_cast<unsigned char *>(bytes);  // просматриваем как массив байтов
  std::array<int, 256> counter = {0};

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    *(counter.data() + index) += 1;
  }
  int tmp = 0;
  int j = 0;
  for (; j < 256; j++) {
    if (*(counter.data() + j) != 0) {
      tmp = *(counter.data() + j);
      *(counter.data() + j) = 0;
      j++;
      break;
    }
  }
  for (; j < 256; j++) {
    int a = *(counter.data() + j);
    *(counter.data() + j) = tmp;
    tmp += a;
  }
  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    out[*(counter.data() + index)] = bytes[i];
    *(counter.data() + index) += 1;
  }
}
uint64_t PapulinaYRadixSortSEQ::InBytes(double d) {
  uint64_t bits = 0;
  memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}
double PapulinaYRadixSortSEQ::FromBytes(uint64_t bits) {
  double d = NAN;
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  memcpy(&d, &bits, sizeof(double));
  return d;
}
void PapulinaYRadixSortSEQ::RadixSort(double *arr, int size) {
  std::vector<uint64_t> bytes(size);
  std::vector<uint64_t> out(size);

  for (int i = 0; i < size; i++) {
    bytes[i] = InBytes(arr[i]);
  }

  SortByByte(bytes.data(), out.data(), 0, size);
  SortByByte(out.data(), bytes.data(), 1, size);
  SortByByte(bytes.data(), out.data(), 2, size);
  SortByByte(out.data(), bytes.data(), 3, size);
  SortByByte(bytes.data(), out.data(), 4, size);
  SortByByte(out.data(), bytes.data(), 5, size);
  SortByByte(bytes.data(), out.data(), 6, size);
  SortByByte(out.data(), bytes.data(), 7, size);

  for (int i = 0; i < size; i++) {
    arr[i] = FromBytes(bytes[i]);
  }
}
}  // namespace papulina_y_radix_sort
