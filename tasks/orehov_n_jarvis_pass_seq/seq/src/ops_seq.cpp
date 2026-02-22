#include "orehov_n_jarvis_pass_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "orehov_n_jarvis_pass_seq/common/include/common.hpp"

namespace orehov_n_jarvis_pass_seq {

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
  input_.assign(tmp.begin(), tmp.end());
  return true;
}

bool OrehovNJarvisPassSEQ::RunImpl() {
  if (input_.size() == 1 || input_.size() == 2){
    res_ = input_;
    return true;
  }

  Point current = FindFirstElem();
  res_.push_back(current);

  while(true){
    Point next = FindNext(current);
    if (next == res_[0]) {break;}

    current = next;
    res_.push_back(next);
  }

  return true;
}

Point OrehovNJarvisPassSEQ::FindNext(Point current) const {
  Point next = current == input_[0] ? input_[1] : input_[0];
  for (auto p : input_){
    if (current == p || next == p){
      continue;
    }
    double orient = CheckLeft(current, next, p); 
    if (orient > 0){
      next = p;
    }
    if (orient == 0){
      if (Distance(current, next) < Distance(current, p)) {next = p;}
    }
  }
  return next;
}

double OrehovNJarvisPassSEQ::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassSEQ::FindFirstElem() const {
  Point current =  input_[0];
  for (auto f: input_){
    if (f.x < current.x || (f.y < current.y && f.x == current.x)){
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassSEQ::Distance(Point a, Point b){
  return std::sqrt(pow(a.y - b.y, 2) + pow(a.x - b.x, 2));
}

bool OrehovNJarvisPassSEQ::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace orehov_n_jarvis_pass_seq
