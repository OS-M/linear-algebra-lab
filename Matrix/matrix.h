#pragma once

#include <iostream>
#include <iomanip>

const double kEps = 1e-6;

template<typename T>
class Matrix {
 public:
  struct DataAccessor {
    inline T& operator[](size_t index2) noexcept {
      return matrix_->data_[matrix_->m_ * index1_ + index2];
    }
    inline const T& operator[](size_t index2) const noexcept {
      return matrix_->data_[matrix_->m_ * index1_ + index2];
    }
    operator T() const {
      if (matrix_->Size().second == 1) {
        return (*this)[0];
      } else if (matrix_->Size().first == 1) {
        return matrix_->data_[index1_];
      }
      throw std::runtime_error("Cannot cast");
    }
   private:
    friend class Matrix<T>;

    Matrix<T>* matrix_;
    size_t index1_{0};
    DataAccessor(Matrix<T>* matrix, size_t index1) : matrix_{matrix},
                                                     index1_{index1} {}
  };

  explicit Matrix(size_t size) : n_{size}, m_{size}, data_size_{size * size} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = T();
    }
  }
  explicit Matrix(size_t n, size_t m) : n_{n}, m_{m}, data_size_{n * m} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = T();
    }
  }
  Matrix(size_t n, size_t m, T default_) : n_{n},
                                           m_{m},
                                           data_size_{n * m} {
    data_ = new T[data_size_];
    for (size_t i = 0; i < data_size_; i++) {
      data_[i] = default_;
    }
  }
  Matrix(const Matrix<T>& matrix) {
    *this = matrix;
  }
  Matrix(Matrix<T>&& matrix) noexcept {
    *this = std::move(matrix);
  }
  ~Matrix() {
    Reset();
  }
  inline DataAccessor operator[](size_t index) noexcept {
    return DataAccessor(this, index);
  }
  inline const DataAccessor operator[](size_t index) const noexcept {
    return DataAccessor(const_cast<Matrix<T>*>(this), index);
  }
  bool Square() const {
    return n_ == m_;
  }
  Matrix& operator=(const Matrix<T>& other) {
    if (other.data_ != data_) {
      Reset();
      n_ = other.n_;
      m_ = other.m_;
      data_size_ = other.data_size_;
      data_ = new T[data_size_];
      for (int i = 0; i < data_size_; i++) {
        data_[i] = other.data_[i];
      }
    }
    return *this;
  }
  Matrix& operator=(Matrix<T>&& other) noexcept {
    std::swap(n_, other.n_);
    std::swap(m_, other.m_);
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
    Matrix<T> result(n_, other.m_);
    for (int i = 0; i < n_; i++) {
      for (int k = 0; k < m_; k++) {
        for (int j = 0; j < other.m_; j++) {
          result[i][j] += this->At(i, k) * other.At(k, j);
        }
      }
    }
    return result;
  }
  std::string ToWolframString() const {
    std::stringstream res;
    res << "{";
    for (int i = 0; i < this->Size().first; i++) {
      res << "{";
      for (int j = 0; j < this->Size().second; j++) {
        res << std::fixed << std::setprecision(4) << this->At(i, j);
        if (j + 1 != this->Size().second) {
          res << ",";
        }
      }
      res << "}";
      if (i + 1 != this->Size().first) {
        res << ",";
      }
    }
    res << "}";
    return res.str();
  }
  void Randomize(int max = 1000) {
    for (int i = 0; i < this->Size().first; i++) {
      for (int j = 0; j < this->Size().second; j++) {
        this->At(i, j) = rand() % max;
      }
    }
  }
  Matrix<T> Transposed() const {
    Matrix<T> res(m_, n_);
    for (int i = 0; i < n_; i++) {
      for (int j = 0; j < m_; j++) {
        res[j][i] = this->At(i, j);
      }
    }
    return res;
  }

  template<class U>
  friend std::ostream& operator<<(std::ostream& stream,
                                  const Matrix<U>& matrix);
  inline std::pair<size_t, size_t> Size() const noexcept {
    return {n_, m_};
  }
  inline const T& At(size_t index1, size_t index2) const {
    // ASSERT(index1 < matrix_size_) << "Index1 " << index1 << " out ob bounds [0; " << matrix_size_ << ")";
    // ASSERT(index2 < matrix_size_) << "Index2 " << index2 << " out ob bounds [0; " << matrix_size_ << ")";
    return data_[m_ * index1 + index2];
  }
  inline T& At(size_t index1, size_t index2) {
    // ASSERT(index1 < matrix_size_) << "Index1 " << index1 << " out ob bounds [0; " << matrix_size_ << ")";
    // ASSERT(index2 < matrix_size_) << "Index2 " << index2 << " out ob bounds [0; " << matrix_size_ << ")";
    return data_[m_ * index1 + index2];
  }

 private:
  void Reset() {
    delete[] data_;
    data_size_ = n_ = m_ = 0;
  }
  size_t n_{0};
  size_t m_{0};
  size_t data_size_{0};
  T* data_{nullptr};
};

template<class U>
std::ostream& operator<<(std::ostream& stream, const Matrix<U>& matrix) {
  size_t maxlen = 0;
  for (size_t i = 0; i < matrix.Size().first; i++) {
    for (size_t j = 0; j < matrix.Size().second; j++) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(5) << matrix[i][j];
      maxlen = std::max(maxlen, ss.str().length());
    }
  }
  stream << "[";
  for (size_t i = 0; i < matrix.Size().first; i++) {
    if (i != 0) {
      std::cout << ' ';
    }
    for (size_t j = 0; j < matrix.Size().second; j++) {
      stream << std::fixed << std::setprecision(5) << std::setw(maxlen)
             << matrix[i][j];
      if (i + 1 < matrix.Size().first || j + 1 < matrix.Size().second) {
        stream << ", ";
      }
    }
    if (i + 1 < matrix.Size().first) {
      stream << '\n';
    }
  }
  stream << "]";
  return stream;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& stream, const std::pair<T, U>& pair) {
  stream << "(" << pair.first << ", " << pair.second << ")";
  return stream;
}
