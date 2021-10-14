#pragma once

#include "abstract_matrix.h"
#include "matrix.h"

template<class T>
class MutableMatrix : public AbstractMatrix<T> {
 public:
  MutableMatrix(size_t n, size_t m) : AbstractMatrix<T>(n, m) {}
  virtual ~MutableMatrix() = default;
  inline typename AbstractMatrix<T>::DataAccessor operator[](
      size_t index) noexcept {
    return AbstractMatrix<T>::DataAccessor::DataAccessor(
        static_cast<AbstractMatrix<T>*>(this), index);
  }
  virtual inline T& At(size_t index1, size_t index2) = 0;
};
