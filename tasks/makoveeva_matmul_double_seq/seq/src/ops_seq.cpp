bool MatmulDoubleSeqTask::RunImpl() {
  if (n_ <= 0) {
    return false;
  }

  // Очищаем C_ перед вычислениями - заменяем std::fill на цикл
  for (size_t i = 0; i < C_.size(); ++i) {
    C_[i] = 0.0;
  }

  const int n_int = static_cast<int>(n_);
  const int block_size = CalculateBlockSize(n_int);
  const int num_blocks = CalculateNumBlocks(n_int, block_size);

  for (int ib = 0; ib < num_blocks; ++ib) {
    for (int jb = 0; jb < num_blocks; ++jb) {
      for (int kb = 0; kb < num_blocks; ++kb) {
        const int i_start = ib * block_size;
        const int i_end = std::min(i_start + block_size, n_int);
        const int j_start = jb * block_size;
        const int j_end = std::min(j_start + block_size, n_int);
        const int k_start = kb * block_size;
        const int k_end = std::min(k_start + block_size, n_int);

        ProcessBlock(A_, B_, C_, n_int, i_start, i_end, j_start, j_end, k_start, k_end);
      }
    }
  }

  GetOutput() = C_;
  return true;
}
