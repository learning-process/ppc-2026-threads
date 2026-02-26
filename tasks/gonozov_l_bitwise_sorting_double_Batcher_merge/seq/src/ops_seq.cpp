#include "gonozov_l_bitwise_sorting_double_Batcher_merge/seq/include/ops_seq.hpp"

#include <chrono>
#include <thread>

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>


#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_Batcher_merge {

GonozovLBitSortBatcherMergeSEQ::GonozovLBitSortBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSEQ::ValidationImpl() {
  return !GetInput().empty(); // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeSEQ::PreProcessingImpl() {
  return true;
}

namespace {
/// double -> uint64_t
uint64_t double_to_sortable_int(double d) {
    uint64_t bits;
    std::memcpy(&bits, &d, sizeof(double));

    if ((bits >> 63) != 0) { // Отрицательное число
        return ~bits;  // Инвертируем ВСЕ биты
    }
    else { // Положительное число или ноль
        return bits ^ 0x8000000000000000ULL; // Инвертируем ТОЛЬКО знаковый бит
    }
}

// uint64_t -> double
double sortable_int_to_double(uint64_t bits) {
    if ((bits >> 63) != 0) { // Если старший бит установлен (было положительное)
        bits = bits ^ 0x8000000000000000ULL; // Возвращаем знаковый бит
    }
    else { // Если старший бит не установлен (было отрицательное)
        bits = ~bits; // Инвертируем все биты обратно
    }

    double result;
    std::memcpy(&result, &bits, sizeof(double));
    return result;
}

void radix_sort_double(std::vector<double>& data) {
    if (data.empty()) return;
    
    // Преобразуем в сортируемые целые числа
    std::vector<uint64_t> keys(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        keys[i] = double_to_sortable_int(data[i]);
    }
    
    const int RADIX = 256; // 8 бит за проход
    std::vector<uint64_t> temp_keys(data.size());
    
    // 8 проходов для 64-битных чисел (8 байт)
    for (int pass = 0; pass < 8; ++pass) {
        std::vector<size_t> count(RADIX, 0);
        int shift = pass * 8;
        
        // Подсчет
        for (uint64_t key : keys) {
            uint8_t byte = (key >> shift) & 0xFF;
            count[byte]++;
        }
        
        // Накопление
        for (int i = 1; i < RADIX; ++i) {
            count[i] += count[i - 1];
        }
        
        // Распределение
        for (int i = keys.size() - 1; i >= 0; --i) {
            uint8_t byte = (keys[i] >> shift) & 0xFF;
            temp_keys[--count[byte]] = keys[i];
        }
        
        keys.swap(temp_keys);
    }
    
    // Преобразуем обратно
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = sortable_int_to_double(keys[i]);
    }
}




void batcher_odd_even_merge(std::vector<double>& arr, int low, int high) {
    if (high - low <= 1) return;
    
    int mid = (low + high) / 2;
    batcher_odd_even_merge(arr, low, mid);
    batcher_odd_even_merge(arr, mid, high);
    
    // Сравниваем и меняем элементы из двух половин
    for (int i = low; i < mid; ++i) {
        int j = i + mid - low; // Индекс во второй половине
        if (j < high && arr[i] > arr[j]) {
            std::swap(arr[i], arr[j]);
        }
    }
}

void hybrid_sort_double(std::vector<double>& data) {
    radix_sort_double(data);
    batcher_odd_even_merge(data, 0, static_cast<int>(data.size()));
}

} // namespace

bool GonozovLBitSortBatcherMergeSEQ::RunImpl() {
  std::vector<double> array = GetInput();
  hybrid_sort_double(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_Batcher_merge
