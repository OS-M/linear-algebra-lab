#pragma once

#include "abstract_matrix.h"
#include "matrix.h"

template<class T>
class TransposedMatrix : public AbstractMatrix<T> {
 public:
  explicit TransposedMatrix(const AbstractMatrix<T>& matrix) :
      AbstractMatrix<T>(matrix.Rows(), matrix.Cols()), matrix_{matrix} {}

  inline T At(size_t index1, size_t index2) const override {
    return matrix_.At(index2, index1);
  }
  size_t Rows() const override {
    return matrix_.Cols();
  }
  size_t Cols() const override {
    return matrix_.Rows();
  }

 private:
  const AbstractMatrix<T>& matrix_;
};
