namespace {

void ReferenceMultiply(const std::vector<double>& a,
                       const std::vector<double>& b,
                       std::vector<double>& c,
                       size_t n) {
  // Заменяем std::fill на цикл
  for (size_t i = 0; i < c.size(); ++i) {
    c[i] = 0.0;
  }
  
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      const double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
}

}  // namespace