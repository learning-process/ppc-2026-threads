#include "orehov_n_Jarvis_pass_seq/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "orehov_n_Jarvis_pass_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace orehov_n_Jarvis_pass_seq {

OrehovNJarvisPassSEQ::OrehovNJarvisPassSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassSEQ::ValidationImpl() {
  return (!GetInput().empty());
}

bool OrehovNJarvisPassSEQ::PreProcessingImpl() {
  std::set<Point> tmp(GetInput().begin(), GetInput().end());
  input.assign(tmp.begin(), tmp.end());
  return true;
}

bool OrehovNJarvisPassSEQ::RunImpl() {
  if (input.size() == 1 || input.size() == 2){
    res = input;
    return true;
  }

  Point current = FindFirstElem();
  res.push_back(current);

  while(true){
    Point next = current == input[0] ? input[1] : input[0];
    for (int i = 0; i < input.size(); i++){
      if (current == input[i] || next == input[i]){
        continue;
      }
      int orient = CheckLeft(current, next, input[i]); 
      if (orient > 0){
        next = input[i];
      }
      if (orient == 0){
        if (distance(current, next) < distance(current, input[i])) next = input[i];
      }
    }
    if (next == res[0]) break;

    current = next;
    res.push_back(next);
  }

  return true;
}

double OrehovNJarvisPassSEQ::CheckLeft(Point A, Point B, Point C) const {
  return (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
}

Point OrehovNJarvisPassSEQ::FindFirstElem() const {
  Point current =  input[0];
  for (auto f: input){
    if (f.x < current.x || (f.y < current.y && f.x == current.x)){
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassSEQ::distance(Point A, Point B) const{
  return std::sqrt(pow(A.y - B.y, 2) + pow(A.x - B.x, 2));
}

bool OrehovNJarvisPassSEQ::PostProcessingImpl() {
  GetOutput() = res;
  return true;
}

}  // namespace orehov_n_Jarvis_pass_seq
