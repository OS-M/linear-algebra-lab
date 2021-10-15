#pragma once

#include <cstdio>
#include <valarray>
#include <sstream>
#include <iomanip>

template<class T>
class TransposedMatrix;

template<class T>
class AbstractMatrix {
 public:
  AbstractMatrix(size_t n, size_t m) : n_{n}, m_{m} {};
  virtual ~AbstractMatrix() = default;
  bool IsOneDimensional() const {
    return n_ == 1 || m_ == 1;
  }
  bool IsSquare() const {
    return n_ == m_;
  }
  virtual void Randomize(int max = 1000) {};
  virtual size_t Rows() const {
    return n_;
  }
  virtual size_t Cols() const {
    return m_;
  }
  virtual inline T At(size_t index1, size_t index2) const = 0;
  TransposedMatrix<T> Transposed() const {
    return TransposedMatrix<T>(*this);
  }
  std::string ToWolframString() const {
    std::stringstream res;
    res << "{";
    for (int i = 0; i < this->Size().first; i++) {
      res << "{";
      for (int j = 0; j < this->Cols(); j++) {
        res << std::fixed << std::setprecision(4) << this->At(i, j);
        if (j + 1 != this->Cols()) {
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

 protected:
  size_t n_{0};
  size_t m_{0};
};

template<class U>
std::ostream& operator<<(std::ostream& stream,
                         const AbstractMatrix<U>& matrix) {
  size_t maxlen = 0;
  for (size_t i = 0; i < matrix.Rows(); i++) {
    for (size_t j = 0; j < matrix.Cols(); j++) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(5) << matrix.At(i, j);
      maxlen = std::max(maxlen, ss.str().length());
    }
  }
  stream << "[";
  for (size_t i = 0; i < matrix.Rows(); i++) {
    if (i != 0) {
      stream << ' ';
    }
    for (size_t j = 0; j < matrix.Cols(); j++) {
      stream << std::fixed << std::setprecision(5) << std::setw(maxlen)
             << matrix.At(i, j);
      if (i + 1 < matrix.Rows() || j + 1 < matrix.Cols()) {
        stream << ", ";
      }
    }
    if (i + 1 < matrix.Rows()) {
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
