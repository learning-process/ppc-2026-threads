#include "egashin_k_radix_simple_merge/seq/include/ops_seq.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include "egashin_k_radix_simple_merge/common/include/common.hpp"

namespace egashin_k_radix_simple_merge {

namespace {

constexpr uint64_t kSignBit = 0x8000000000000000ULL;

}  // namespace

EgashinKRadixSimpleMergeSEQ::EgashinKRadixSimpleMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool EgashinKRadixSimpleMergeSEQ::ValidationImpl() {
  return true;
}

bool EgashinKRadixSimpleMergeSEQ::PreProcessingImpl() {
  result_ = GetInput();
  return true;
}

bool EgashinKRadixSimpleMergeSEQ::RunImpl() {
  RadixSort(result_);
  return true;
}

bool EgashinKRadixSimpleMergeSEQ::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

void EgashinKRadixSimpleMergeSEQ::CountingPass(const std::vector<uint64_t> &source, std::vector<uint64_t> &destination,
                                               int byte_index) {
  std::array<size_t, 256> count{};

  for (uint64_t value : source) {
    const auto byte = static_cast<uint8_t>((value >> (byte_index * 8)) & 0xFFU);
    count.at(byte)++;
  }

  std::array<size_t, 256> position{};
  for (size_t i = 1; i < count.size(); ++i) {
    position.at(i) = position.at(i - 1) + count.at(i - 1);
  }

  for (uint64_t value : source) {
    const auto byte = static_cast<uint8_t>((value >> (byte_index * 8)) & 0xFFU);
    destination[position.at(byte)++] = value;
  }
}

void EgashinKRadixSimpleMergeSEQ::RadixSort(std::vector<double> &data) {
  if (data.size() < 2) {
    return;
  }

  std::vector<uint64_t> keys(data.size());
  std::vector<uint64_t> buffer(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    uint64_t bits = 0;
    std::memcpy(&bits, &data[i], sizeof(data[i]));
    keys[i] = ((bits & kSignBit) != 0U) ? ~bits : (bits ^ kSignBit);
  }

  auto *source = &keys;
  auto *destination = &buffer;
  for (int byte_index = 0; byte_index < 8; ++byte_index) {
    CountingPass(*source, *destination, byte_index);
    std::swap(source, destination);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    const uint64_t bits = (((*source)[i] & kSignBit) != 0U) ? ((*source)[i] ^ kSignBit) : ~(*source)[i];
    std::memcpy(&data[i], &bits, sizeof(data[i]));
  }
}

}  // namespace egashin_k_radix_simple_merge
