#pragma once

#include "abstract_matrix.h"

template<class T>
class MutableMatrix : public AbstractMatrix<T> {
 public:
  MutableMatrix(size_t n, size_t m) : AbstractMatrix<T>(n, m) {}
  virtual ~MutableMatrix() = default;
  virtual inline T& At(size_t index1, size_t index2) = 0;
  virtual MutableMatrix<T>& operator=(const AbstractMatrix<T>& matrix) = 0;
};
