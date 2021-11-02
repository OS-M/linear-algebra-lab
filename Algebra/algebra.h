#pragma once

#include <numeric>
#include "Matrix/abstract_matrix.h"
#include "Matrix/matrix.h"

namespace algebra {

template<class T>
void GetLu(const AbstractMatrix<T>& a,
           MutableMatrix<T>& l,
           MutableMatrix<T>& u);

template<class T>
void GetLdlt(const AbstractMatrix<T>& a,
             MutableMatrix<T>& l,
             MutableMatrix<T>& d);

template<class T>
Matrix<T> SolveLxb(const AbstractMatrix<T>& l,
                   const AbstractMatrix<T>& b);

template<class T>
Matrix<T> SolveUxb(const AbstractMatrix<T>& u,
                   const AbstractMatrix<T>& b);

template<class T>
Matrix<T> LuSolve(const AbstractMatrix<T>& l,
                  const AbstractMatrix<T>& u,
                  const AbstractMatrix<T>& b);

template<class T>
Matrix<T> LdltSolve(const AbstractMatrix<T>& l,
                    const AbstractMatrix<T>& d,
                    const AbstractMatrix<T>& b);

template<class T>
Matrix<T> GaussSeidelSolve(const AbstractMatrix<T>& a,
                           const AbstractMatrix<T>& b,
                           T eps,
                           int check_metric_every = 10,
                           int iteration_limit = 10000,
                           int* iterations = nullptr);

template<class T>
Matrix<T> GaussSolve(const AbstractMatrix<T>& a,
                     const AbstractMatrix<T>& b,
                     T eps = 1e-10);

template<class T>
T Norm2D(const AbstractMatrix<T>& a) {
  T res = 0;
  for (int i = 0; i < a.Rows(); i++) {
    for (int j = 0; j < a.Cols(); j++) {
      res += a.At(i, j) * a.At(i, j);
    }
  }
  return std::sqrt(res);
}

}

namespace algebra {

template<class T>
void GetLu(const AbstractMatrix<T>& a,
           MutableMatrix<T>& l,
           MutableMatrix<T>& u) {
  if (!a.IsSquare()) {
    throw std::runtime_error("A is not square");
  }

  u = a;
  auto n = a.Rows();
  for (int k = 0; k < n - 1; k++) {
    for (int i = k + 1; i < n; i++) {
      u.At(i, k) /= u.At(k, k);
      for (int j = k + 1; j < n; j++) {
        u.At(i, j) -= u.At(i, k) * u.At(k, j);
      }
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      l.At(i, j) = u.At(i, j);
      if (i == j) {
        l.At(i, j) = 1;
      } else {
        u.At(i, j) = 0;
      }
    }
  }
}

template<class T>
void GetLdlt(const AbstractMatrix<T>& a,
             MutableMatrix<T>& l,
             MutableMatrix<T>& d) {
  int n = a.Rows();
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T sum = a.At(j, i);
      for (int k = 0; k < i; k++) {
        sum -= l.At(i, k) * d.At(k, 0) * l.At(j, k);
      }
      if (i == j) {
        // if (sum <= 0) {
        // throw std::runtime_error("A is not positive defined");
        // }
        d.At(i, 0) = sum;
        l.At(i, i) = 1;
      } else {
        l.At(j, i) = sum / d.At(i, 0);
      }
    }
  }
}

template<class T>
Matrix<T> SolveLxb(const AbstractMatrix<T>& l,
                   const AbstractMatrix<T>& b) {
  Matrix<T> x(b.Rows(), 1);
  for (int i = 0; i < x.Rows(); i++) {
    T sum = b.At(i, 0);
    for (int j = 0; j < i; j++) {
      sum -= x.At(j, 0) * l.At(i, j);
    }
    x.At(i, 0) = sum / l.At(i, i);
  }
  return x;
}

template<class T>
Matrix<T> SolveUxb(const AbstractMatrix<T>& u,
                   const AbstractMatrix<T>& b) {
  Matrix<T> x(b.Rows(), 1);
  for (int i = x.Rows() - 1; i >= 0; i--) {
    T sum = b.At(i, 0);
    for (int j = x.Rows() - 1; j > i; j--) {
      sum -= x.At(j, 0) * u.At(i, j);
    }
    x.At(i, 0) = sum / u.At(i, i);
  }
  return x;
}

template<class T>
Matrix<T> LuSolve(const AbstractMatrix<T>& l,
                  const AbstractMatrix<T>& u,
                  const AbstractMatrix<T>& b) {
  return SolveUxb(u, SolveLxb(l, b));
}

template<class T>
Matrix<T> LdltSolve(const AbstractMatrix<T>& l,
                    const AbstractMatrix<T>& d,
                    const AbstractMatrix<T>& l_t,
                    const AbstractMatrix<T>& b) {
  auto n = d.Rows();
  auto y = SolveLxb(l, b);
  for (int i = 0; i < n; i++) {
    y.At(i, 0) /= d.At(i, 0);
  }
  return SolveUxb(l_t, y);
}

template<class T>
Matrix<T> GaussSeidelSolve(const AbstractMatrix<T>& a,
                           const AbstractMatrix<T>& b,
                           T eps,
                           int check_metric_every,
                           int iteration_limit,
                           int* iterations) {
  Matrix<T> x(b.Rows(), 1);
  int counter = 0;
  while (counter++ < iteration_limit) {
    for (int i = 0; i < b.Rows(); i++) {
      T now = 0;
      for (int j = 0; j < b.Rows(); j++) {
        if (i != j) {
          now += x.At(j, 0) * a.At(i, j);
        }
      }
      x.At(i, 0) = (b.At(i, 0) - now) / a.At(i, i);
    }
    if (counter % check_metric_every == 0 && Norm2D(a * x - b) < eps) {
      break;
    }
  }
  if (iterations) {
    *iterations = counter;
  }
  return x;
}

template<class T>
Matrix<T> GaussSolve(const AbstractMatrix<T>& a_,
                     const AbstractMatrix<T>& b_,
                     T eps) {
  auto n = a_.Rows();
  auto a = Matrix<T>::FromAbstract(a_);
  auto b = Matrix<T>::FromAbstract(b_);
  std::vector<size_t> reindex(n);
  std::iota(reindex.begin(), reindex.end(), 0);

  for (int i = 0; i < n; i++) {
    size_t index_of_max = i;
    for (int k = i + 1; k < n; k++) {
      if (a.At(reindex[index_of_max], i) < a.At(reindex[k], i)) {
        index_of_max = k;
      }
    }
    std::swap(reindex[i], reindex[index_of_max]);
    if (std::abs(a.At(reindex[i], i)) < eps) {
      continue;
    }
    for (int j = i + 1; j < n; j++) {
      if (std::abs(a.At(reindex[j], i)) < eps) {
        continue;
      }
      auto m = a.At(reindex[j], i) / a.At(reindex[i], i);
      for (int k = i + 1; k < n; k++) {
        a.At(reindex[j], k) -= m * a.At(reindex[i], k);
      }
      a.At(reindex[j], i) = 0;
      b.At(reindex[j], 0) -= m * b.At(reindex[i], 0);
    }
  }

  auto swapped_a = a;
  auto swapped_b = b;
  for (int i = 0; i < n; i++) {
    swapped_b.At(i, 0) = b.At(reindex[i], 0);
    for (int j = 0; j < n; j++) {
      swapped_a.At(i, j) = a.At(reindex[i], j);
    }
  }
  return SolveUxb(swapped_a, swapped_b);
}

}
