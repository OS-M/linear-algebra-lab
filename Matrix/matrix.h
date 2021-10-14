#pragma once

#include <iostream>
#include <iomanip>
#include "abstract_matrix.h"

const double kEps = 1e-6;

template<typename T>
class Matrix : public AbstractMatrix<T> {
 public:
  explicit Matrix(size_t size) : AbstractMatrix<T>(size, size),
                                 data_size_{size * size} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = T();
    }
  }
  explicit Matrix(size_t n, size_t m) : AbstractMatrix<T>(n, m),
                                        data_size_{n * m} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = T();
    }
  }
  Matrix(size_t n, size_t m, T default_) : AbstractMatrix<T>(n, m),
                                           data_size_{n * m} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = default_;
    }
  }
  Matrix(const Matrix<T>& matrix) : AbstractMatrix<T>(matrix.n_, matrix.m_) {
    *this = matrix;
  }
  Matrix(Matrix<T>&& matrix) noexcept: AbstractMatrix<T>(matrix.n_, matrix.m_) {
    *this = std::move(matrix);
  }
  virtual ~Matrix() {
    Reset();
  }
  Matrix& operator=(const Matrix<T>& other) {
    if (other.data_ != data_) {
      Reset();
      this->n_ = other.n_;
      this->m_ = other.m_;
      data_size_ = other.data_size_;
      data_ = new T[data_size_];
      for (int i = 0; i < data_size_; i++) {
        data_[i] = other.data_[i];
      }
    }
    return *this;
  }
  Matrix& operator=(Matrix<T>&& other) noexcept {
    std::swap(this->n_, other.n_);
    std::swap(this->m_, other.m_);
    std::swap(data_size_, other.data_size_);
    std::swap(data_, other.data_);
    return *this;
  }
  bool operator==(const Matrix& other) const {
    if (data_size_ != other.data_size_) {
      return false;
    }
    if (data_ == other.data_) {
      return true;
    }
    for (size_t i = 0; i < data_size_; i++) {
      if (std::abs(data_[i] - other.data_[i]) > kEps) {
        return false;
      }
    }
    return true;
  }
  Matrix operator*(const Matrix& other) const {
    if (this->Size().second != other.Size().first) {
      throw std::runtime_error(
          "Bad matrix sizes " + std::to_string(this->Size().second) + ' '
              + std::to_string(other.Size().first));
    }
    Matrix<T> result(this->n_, other.m_);
    for (int i = 0; i < this->n_; i++) {
      for (int k = 0; k < this->m_; k++) {
        for (int j = 0; j < other.m_; j++) {
          result[i][j] += this->At(i, k) * other.At(k, j);
        }
      }
    }
    return result;
  }
  void Randomize(int max = 1000) override {
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        this->At(i, j) = rand() % max;
      }
    }
  }
  Matrix<T> Transposed() const {
    Matrix<T> res(this->m_, this->n_);
    for (int i = 0; i < this->n_; i++) {
      for (int j = 0; j < this->m_; j++) {
        res[j][i] = this->At(i, j);
      }
    }
    return res;
  }

  template<class U>
  friend std::ostream& operator<<(std::ostream& stream,
                                  const Matrix<U>& matrix);
  inline const T& At(size_t index1, size_t index2) const override {
    return data_[this->m_ * index1 + index2];
  }
  inline T& At(size_t index1, size_t index2) override {
    return data_[this->m_ * index1 + index2];
  }

 private:
  void Reset() {
    delete[] data_;
    data_size_ = 0;
  }
  size_t data_size_{0};
  T* data_{nullptr};
};
