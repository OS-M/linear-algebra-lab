#pragma once

#include "Matrix/matrix.h"

namespace algebra {

template<class T>
void GetLu(const Matrix<T>& a, Matrix<T>& l, Matrix<T>& u);

template<class T>
void GetLdlt(const Matrix<T>& a, Matrix<T>& l, Matrix<T>& d);

template<class T>
Matrix<T> SolveLxb(const Matrix<T>& l,
                   const Matrix<T>& b);

template<class T>
Matrix<T> SolveUxb(const Matrix<T>& u,
                   const Matrix<T>& b);

template<class T>
Matrix<T> LuSolve(const Matrix<T>& l,
                  const Matrix<T>& u,
                  const Matrix<T>& b);

template<class T>
Matrix<T> LdltSolve(const Matrix<T>& l,
                    const Matrix<T>& d,
                    const Matrix<T>& b);

}

namespace algebra {

template<class T>
void GetLu(const Matrix<T>& a, Matrix<T>& l, Matrix<T>& u) {
  if (!a.Square()) {
    throw std::runtime_error("A is not square");
  }
  auto n = a.Size().first;
  u = a;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      l[j][i] = u[j][i] / u[i][i];
    }
  }

  for (int k = 1; k < n; k++) {
    for (int i = k - 1; i < n; i++) {
      for (int j = i; j < n; j++) {
        l[j][i] = u[j][i] / u[i][i];
      }
    }

    for (int i = k; i < n; i++) {
      for (int j = k - 1; j < n; j++) {
        u[i][j] = u[i][j] - l[i][k - 1] * u[k - 1][j];
      }
    }
  }
}

template<class T>
void GetLdlt(const Matrix<T>& a, Matrix<T>& l, Matrix<T>& d) {
  int n = a.Size().first;
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      T sum = a[j][i];
      for (int k = 0; k < i; k++) {
        sum -= l[i][k] * d[k] * l[j][k];
      }
      if (i == j) {
        // if (sum <= 0) {
        // throw std::runtime_error("A is not positive defined");
        // }
        d[i][0] = sum;
        l[i][i] = 1;
      } else {
        l[j][i] = sum / d[i];
      }
    }
  }
}

template<class T>
Matrix<T> SolveLxb(const Matrix<T>& l,
                   const Matrix<T>& b) {
  Matrix<T> x(b.Size().first, 1);
  for (int i = 0; i < x.Size().first; i++) {
    T sum = b[i][0];
    for (int j = 0; j < i; j++) {
      sum -= x[j][0] * l[i][j];
    }
    x[i][0] = sum / l[i][i];
  }
  return x;
}

template<class T>
Matrix<T> SolveUxb(const Matrix<T>& u,
                   const Matrix<T>& b) {
  Matrix<T> x(b.Size().first, 1);
  for (int i = x.Size().first - 1; i >= 0; i--) {
    T sum = b[i][0];
    for (int j = x.Size().first - 1; j > i; j--) {
      sum -= x[j][0] * u[i][j];
    }
    x[i][0] = sum / u[i][i];
  }
  return x;
}

template<class T>
Matrix<T> LuSolve(const Matrix<T>& l,
                  const Matrix<T>& u,
                  const Matrix<T>& b) {
  return SolveUxb(u, SolveLxb(l, b));
}

template<class T>
Matrix<T> LdltSolve(const Matrix<T>& l,
                    const Matrix<T>& d,
                    const Matrix<T>& b) {
  auto n = d.Size().first;
  auto y = SolveLxb(l, b);
  for (int i = 0; i < n; i++) {
    y[i][0] /= d[i];
  }
  return SolveUxb(l.Transposed(), y);
}

}
