#pragma once

#include <iostream>
#include <iomanip>
#include "mutable_matrix.h"

const double kEps = 1e-6;

template<typename T>
class Matrix : public MutableMatrix<T> {
 public:
  explicit Matrix(size_t size) : Matrix(size, size) {}
  explicit Matrix(size_t n, size_t m) : Matrix(n, m, T()) {}
  explicit Matrix(const std::initializer_list<std::initializer_list<T>>& list) :
      MutableMatrix<T>(list.size(), list.begin()->size()),
      data_size_{this->n_ * this->m_},
      data_{new T[data_size_]} {
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        this->At(i, j) = *((list.begin() + i)->begin() + j);
      }
    }
  }
  Matrix(size_t n, size_t m, T default_) : MutableMatrix<T>(n, m),
                                           data_size_{n * m},
                                           data_{new T[data_size_]} {
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = default_;
    }
  }
  static Matrix<T> FromAbstract(const AbstractMatrix<T>& matrix) {
    Matrix<T> res(matrix.Rows(), matrix.Cols());
    res = matrix;
    return res;
  }
  Matrix(const Matrix<T>& matrix) :
      MutableMatrix<T>(matrix.n_, matrix.m_) {
    *this = matrix;
  }
  Matrix(Matrix<T>&& matrix) noexcept: MutableMatrix<T>(matrix.n_, matrix.m_) {
    *this = std::move(matrix);
  }
  virtual ~Matrix() {
    delete[] data_;
  }
  static Matrix<T> Ones(size_t n) {
    Matrix<T> res(n);
    for (int i = 0; i < n; i++) {
      res.At(i, i) = 1;
    }
    return res;
  }
  Matrix operator+(const AbstractMatrix<T>& other) const {
    Matrix res(*this);
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        res.At(i, j) += other.At(i, j);
      }
    }
    return res;
  }
  Matrix operator-() const {
    Matrix res(*this);
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        res.At(i, j) *= -1;
      }
    }
    return res;
  }
  Matrix operator-(const AbstractMatrix<T>& other) const {
    Matrix res(*this);
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        res.At(i, j) -= other.At(i, j);
      }
    }
    return res;
  }
  Matrix& operator=(const AbstractMatrix<T>& other) override {
    if (!this->data_ || *this != other) {
      this->n_ = other.Rows();
      this->m_ = other.Cols();
      data_size_ = this->n_ * this->m_;
      delete[] data_;
      data_ = new T[data_size_];
      for (int i = 0; i < this->Rows(); i++) {
        for (int j = 0; j < this->Cols(); j++) {
          this->At(i, j) = other.At(i, j);
        }
      }
    }
    return *this;
  }
  Matrix& operator=(const Matrix<T>& other) {
    *this = static_cast<const AbstractMatrix<T>&>(other);
    return *this;
  }
  Matrix& operator=(Matrix<T>&& other) noexcept {
    std::swap(this->n_, other.n_);
    std::swap(this->m_, other.m_);
    std::swap(data_size_, other.data_size_);
    std::swap(data_, other.data_);
    return *this;
  }
  Matrix operator*(const Matrix& other) const {
    return static_cast<const AbstractMatrix<T>&>(*this) * other;
  }
  void Randomize(int max = 1000) override {
    for (int i = 0; i < this->Rows(); i++) {
      for (int j = 0; j < this->Cols(); j++) {
        this->At(i, j) = rand() % max;
      }
    }
  }
  Matrix operator*(T value) const {
    for (int i = 0; i < data_size_; i++) {
      data_[i] *= value;
    }
    return *this;
  }
  Matrix operator/(T value) const {
    for (int i = 0; i < data_size_; i++) {
      data_[i] /= value;
    }
    return *this;
  }

  inline T At(size_t index1, size_t index2) const override {
#ifndef NDEBUG
    if (index1 >= this->Rows() || index2 >= this->Cols()) {
      throw std::out_of_range(
          "Indexes " + std::to_string(index1) + " " + std::to_string(index2)
              + " out of size " + std::to_string(this->Rows()) + " "
              + std::to_string(this->Cols()));
    }
#endif
    return data_[this->m_ * index1 + index2];
  }

  inline T& At(size_t index1, size_t index2) override {
#ifndef NDEBUG
    if (index1 >= this->Rows() || index2 >= this->Cols()) {
      throw std::out_of_range(
          "Indexes " + std::to_string(index1) + " " + std::to_string(index2)
              + " out of size " + std::to_string(this->Rows()) + " "
              + std::to_string(this->Cols()));
    }
#endif
    return data_[this->m_ * index1 + index2];
  }

 private:
  size_t data_size_{0};
  T* data_{nullptr};
};

template<class T>
Matrix<T> operator*(const AbstractMatrix<T>& lhs,
                    const AbstractMatrix<T>& rhs) {
  if (lhs.Cols() != rhs.Rows()) {
    throw std::runtime_error(
        "Bad matrix sizes " + std::to_string(lhs.Cols()) + ' '
            + std::to_string(rhs.Rows()));
  }
  Matrix<T> result(lhs.Rows(), rhs.Cols());
  for (int i = 0; i < lhs.Rows(); i++) {
    for (int k = 0; k < lhs.Cols(); k++) {
      for (int j = 0; j < rhs.Cols(); j++) {
        result.At(i, j) += lhs.At(i, k) * rhs.At(k, j);
      }
    }
  }
  return result;
}

template<class T>
bool operator==(const AbstractMatrix<T>& lhs,
                const AbstractMatrix<T>& rhs) {
  if (lhs.Rows() != rhs.Rows() || lhs.Cols() != rhs.Cols()) {
    return false;
  }
  for (int i = 0; i < lhs.Rows(); i++) {
    for (int j = 0; j < lhs.Cols(); j++) {
      if (std::abs(lhs.At(i, j) - rhs.At(i, j)) > kEps) {
        return false;
      }
    }
  }
  return true;
}
