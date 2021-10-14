#pragma once

#include <cstdio>
#include <valarray>
#include <sstream>
#include <iomanip>

template<class T>
class AbstractMatrix {
 protected:
  class DataAccessor {
   public:
    DataAccessor(AbstractMatrix<T>* matrix, size_t index1) :
        matrix_{matrix},
        index1_{index1} {}
    inline T& operator[](size_t index2) noexcept {
      return matrix_->At(index1_, index2);
    }
    inline const T& operator[](size_t index2) const noexcept {
      return matrix_->At(index1_, index2);
    }
    operator T() const {
      if (matrix_->Size().second == 1) {
        return (*this)[0];
      } else if (matrix_->Size().first == 1) {
        return matrix_->At(index1_, 0);
      }
      throw std::runtime_error("Cannot cast");
    }

   private:
    AbstractMatrix<T>* matrix_;
    size_t index1_{0};
  };

 public:
  AbstractMatrix(size_t n, size_t m) : n_{n}, m_{m} {};
  virtual ~AbstractMatrix() = default;
  inline DataAccessor operator[](size_t index) noexcept {
    return DataAccessor(this, index);
  }
  inline const DataAccessor operator[](size_t index) const noexcept {
    return DataAccessor(const_cast<AbstractMatrix<T>*>(this), index);
  }
  bool IsOneDimensional() const {
    return n_ == 1 || m_ == 1;
  }
  bool IsSquare() const {
    return n_ == m_;
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
  virtual void Randomize(int max = 1000) = 0;

  template<class U>
  friend std::ostream& operator<<(std::ostream& stream,
                                  const AbstractMatrix<U>& matrix);
  inline std::pair<size_t, size_t> Size() const noexcept {
    return {n_, m_};
  }
  size_t Rows() const {
    return n_;
  }
  size_t Cols() const {
    return m_;
  }

  virtual inline const T& At(size_t index1, size_t index2) const = 0;
  virtual inline T& At(size_t index1, size_t index2) = 0;

 protected:
  size_t n_{0};
  size_t m_{0};
};

template<class U>
std::ostream& operator<<(std::ostream& stream,
                         const AbstractMatrix<U>& matrix) {
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
      stream << ' ';
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
