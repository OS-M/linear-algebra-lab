#pragma once

#include "abstract_matrix.h"

template<class T>
class DiagonalBoxMatrix : public AbstractMatrix<T> {
 public:
  explicit DiagonalBoxMatrix(size_t n) : AbstractMatrix<T>(n, n) {}
  virtual ~DiagonalBoxMatrix() = default;

  inline T At(size_t index1, size_t index2) const override {
    if (index1 == index2) {
      return this->n_;
    }
    if (index1 == 0 || index2 == 0
        || index1 == this->n_ - 1 || index2 == this->n_ - 1) {
      return 1;
    }
    return 0;
  }
};
