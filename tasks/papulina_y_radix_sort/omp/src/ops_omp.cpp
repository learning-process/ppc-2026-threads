#include "papulina_y_radix_sort/omp/include/ops_omp.hpp"


#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring> 
#include <omp.h>
#include <utility>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"

namespace papulina_y_radix_sort {

PapulinaYRadixSortOMP::PapulinaYRadixSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}
bool PapulinaYRadixSortOMP::ValidationImpl() {
  return true;
}

bool PapulinaYRadixSortOMP::PreProcessingImpl() {
  return true;
}
bool PapulinaYRadixSortOMP::RunImpl() {
  double * result = GetInput().data();
  int size = static_cast<int>(GetInput().size());
  int threads_count = std::min(omp_get_max_threads(), std::max(1, size / 1000));

  std::vector<std::span<double>> chunks;
  std::vector<int> chunks_offsets;

  int chunk_size = size/threads_count;
  int reminder = size%threads_count;

  int offset = 0;
  for(int i=0; i<threads_count; i++){
    int real_chunk_size = chunk_size + (i < reminder ? 1 : 0);
    if ( real_chunk_size > 0){
      chunks_offsets.push_back(offset);
      chunks.emplace_back(result + offset, real_chunk_size);
      offset+=real_chunk_size;
    }
  }
  threads_count = static_cast<int>(chunks.size()); // тк возможно chunks.size() < threads_count(каким-то потокам ничего не распределилось из данных)
  #pragma omp parallel for default(none) shared(result, chunks, threads_count) num_threads(threads_count)
  for(int i=0; i<threads_count; i++){
    RadixSort(chunks[i].data(), static_cast<int>(chunks[i].size()));
  }
  MergeChunks(chunks,chunks_offsets,result);

  GetOutput() = std::vector<double>(size);
  for (int i = 0; i < size; i++) {
    GetOutput()[i] = result[i];
    //std::cout << result [i] << " ";
  }
  //std::cout << std::endl;

  return true;
}
void PapulinaYRadixSortOMP::Merge(std::span<double> & res, const std::span<double> & left, const std::span<double> & right){
  int i=0;
  int j=0;
  int k=0;
  while(static_cast<size_t>(i)<left.size() && static_cast<size_t>(j)<right.size()){
    if(left[i] <= right[j]){
      res[k++] = left[i++];
    }
    else{
      res[k++] = right[j++];
    }
  }
  while(static_cast<size_t>(i)<left.size() ){
    res[k++] = left[i++];
  }
  while(static_cast<size_t>(j)<right.size() ){
    res[k++] = right[j++];
  }  

}
void PapulinaYRadixSortOMP::MergeChunks(std::vector<std::span<double>> chunks, std::vector<int> chunks_offsets, double * result){
  if(chunks.size() <= 1) {return;}
  int n = static_cast<int>(GetInput().size());
  int chunks_count = static_cast<int>(chunks.size());
  std::vector<double> buffer(n);

  for(int step = 1; step < chunks_count; step*=2){
    #pragma omp parallel for default(none) shared(chunks_count, step, chunks_offsets, n, result) num_threads((chunks_count + 1)/2)
    for(int i=0; i < chunks_count; i+=2*step){
      if(i+step<chunks_count){
        int left_idx = i;
        int right_idx = i+step;

        int left_start = chunks_offsets[left_idx];
        
        int right_end = (i+2*step<chunks_count) ? chunks_offsets[i+2*step] : n;
        
        std::span<double> left(result + left_start, chunks_offsets[right_idx] - left_start);
        std::span<double> right(result + chunks_offsets[right_idx],right_end - chunks_offsets[right_idx]);
        std::span<double> res(result + left_start, right_end - left_start);
        Merge(res,left,right);
      }
    }
    std::vector<int> new_offsets;
    for(int i=0; i<chunks_count; i+=2*step){
      new_offsets.push_back(chunks_offsets[i]);
    }
    new_offsets.push_back(n);
    chunks_offsets = std::move(new_offsets);

    std::vector<std::span<double>> new_chunks;
    new_chunks.reserve(chunks_offsets.size() - 1); 
    for(int i=0; i<static_cast<int>(chunks_offsets.size())-1; i++){
        new_chunks.emplace_back(result+chunks_offsets[i], chunks_offsets[i+1]-chunks_offsets[i]);
    }
    chunks = std::move(new_chunks);
    chunks_count = static_cast<int>(chunks.size());
  }
}
void PapulinaYRadixSortOMP::RadixSort(double *arr, int size) {
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
bool PapulinaYRadixSortOMP::PostProcessingImpl() {
  return true;
}
void PapulinaYRadixSortOMP::SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size) {
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
uint64_t PapulinaYRadixSortOMP::InBytes(double d) {
  uint64_t bits = 0;
  memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}
double PapulinaYRadixSortOMP::FromBytes(uint64_t bits) {
  double d = NAN;
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  memcpy(&d, &bits, sizeof(double));
  return d;
}
}  // namespace papulina_y_radix_sort
