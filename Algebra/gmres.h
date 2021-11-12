#pragma once

#include "Matrix/matrix.h"

namespace algebra {

template<class T>
std::pair<T, T> GetGivens(T rho, T sigma);

template<class T>
std::pair<Matrix<T>, int> Gmres(const AbstractMatrix<T>& a,
                                const Matrix<T>& b,
                                const std::function<T(
                                    const AbstractMatrix<T>&)>& norm,
                                int max_iter,
                                T threshold = 1e-10) {
  int n = a.Rows();
  Matrix<T> x0(n, 1);
  // x0.At(0, 0) = norm(b);
  x0 = b / norm(b);

  auto r0 = b - a * x0;
  T norm_r0 = norm(r0);
  T norm_b = norm(b);
  T error = norm_r0 / norm_b;

  Matrix<T> cos_v(max_iter, 1);
  Matrix<T> sin_v(max_iter, 1);
  Matrix<T> beta(max_iter, 1);
  beta.At(0, 0) = norm_r0;
  Matrix<T> y(max_iter, 1);
  Matrix<T> h(max_iter, max_iter + 1);
  Matrix<T> v(n, max_iter);
  Matrix<T> g_n(0);
  Matrix<T> delta_x(0);
  Matrix<T> v_iter(n, 1);

  auto v_1 = r0 / norm_r0;
  for (int i = 0; i < n; ++i) {
    v.At(i, 0) = v_1.At(i, 0);
  }

  int iter = 0;
  while (error > threshold && iter + 1 < max_iter) {
    for (int i = 0; i < n; ++i) {
      v_iter.At(i, 0) = v.At(i, iter);
    }

    // Arnoldi
    auto w = a * v_iter;
    for (int j = 0; j < iter + 1; ++j) {
      for (int i = 0; i < n; ++i) {
        h.At(j, iter) += w.At(i, 0) * v.At(i, j);
      }
      for (int k = 0; k < n; ++k) {
        w.At(k, 0) -= h.At(j, iter) * v.At(k, j);
      }
    }
    h.At(iter + 1, iter) = norm(w);

    for (int i = 0; i < n; ++i) {
      v.At(i, iter + 1) = w.At(i, 0) / h.At(iter + 1, iter);
    }

    for (int i = 0; i < iter; ++i) {
      auto temp = cos_v.At(i, 0) * h.At(i, iter)
          + sin_v.At(i, 0) * h.At(i + 1, iter);
      h.At(i + 1, iter) = -sin_v.At(i, 0) * h.At(i, iter)
          + cos_v.At(i, 0) * h.At(i + 1, iter);
      h.At(i, iter) = temp;
    }

    auto[cos, sin] = GetGivens(h.At(iter, iter), h.At(iter + 1, iter));
    cos_v.At(iter, 0) = cos;
    sin_v.At(iter, 0) = sin;

    h.At(iter, iter) = cos_v.At(iter, 0) * h.At(iter, iter)
        + sin_v.At(iter, 0) * h.At(iter + 1, iter);
    h.At(iter + 1, iter) = 0;

    beta.At(iter + 1, 0) = -sin_v.At(iter, 0) * beta.At(iter, 0);
    beta.At(iter, 0) = cos_v.At(iter, 0) * beta.At(iter, 0);

    error = std::fabs(beta.At(iter + 1, 0)) / norm_b;
    iter += 1;
  }

  if (error > threshold) {
    std::cerr << "No\n";
  }

  g_n = Matrix<T>(iter + 1, 1);
  for (int i = 0; i < std::min(g_n.Rows(), beta.Rows()); i++) {
    g_n.At(i, 0) = beta.At(i, 0);
  }
  for (int i = iter - 1; i >= 0; i--) {
    y.At(i, 0) = g_n.At(i, 0) / h.At(i, i);
    if (std::fabs(y.At(i, 0)) < threshold) {
      continue;
    }
    for (int k = i - 1; k >= 0; k--) {
      g_n.At(k, 0) -= h.At(k, i) * y.At(i, 0);
    }
  }
  delta_x = Matrix<T>(n, 1);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < iter; ++j) {
      delta_x.At(i, 0) += v.At(i, j) * y.At(j, 0);
    }
  }
  return {x0 + delta_x, iter};
}

template<class T>
std::pair<T, T> GetGivens(T rho, T sigma) {
  if (std::fabs(rho) < 1e-10) {
    return {0, 1};
  } else {
    T d = sqrt(rho * rho + sigma * sigma);
    return {std::fabs(rho) / d, std::fabs(rho) / d * sigma / rho};
  }
}

}
