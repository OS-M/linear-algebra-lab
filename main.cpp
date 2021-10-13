#include <iostream>

#include "Matrix/matrix.h"
#include "TimeMeasurer/time_measurer.h"

void TestMatrixMult() {
  for (int k = 0; k < 1000; k++) {
    Matrix<double> a(100);
    Matrix<double> b(100);
    a.Randomize();
    b.Randomize();
    auto c = a * b;
  }
}

template<class T>
void LU(const Matrix<T>& a, Matrix<T>& l, Matrix<T>& u) {
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
Matrix<T> Solve(const Matrix<T>& l, const Matrix<T>& u, const Matrix<T>& b) {
  // Ly = b;
  Matrix<T> y(b.Size().first, 1);
  for (int i = 0; i < y.Size().first; i++) {
    T sum = b[i][0];
    for (int j = 0; j < i; j++) {
      sum -= y[j][0] * l[i][j];
    }
    y[i][0] = sum;
  }
  // Ux = y;
  Matrix<T> x(b.Size().first, 1);
  for (int i = x.Size().first - 1; i >= 0; i--) {
    T sum = y[i][0];
    for (int j = x.Size().first - 1; j > i; j--) {
      sum -= x[j][0] * u[i][j];
    }
    x[i][0] = sum / u[i][i];
  }
  return x;
}

int main() {
  srand(time(0));
  TimeMeasurer time_measurer;
  int n = 3;
  Matrix<double> a(n);
  a[0][0] = a[1][1] = 2;
  a[1][0] = 3;
  a.Randomize();
  Matrix<double> l(n);
  Matrix<double> u(n);
  LU(a, l, u);
  // std::cout << l << '\n' << u << '\n' << l * u << '\n' << a << '\n';
  Matrix<double> b(n, 1);
  b[0][0] = 4;
  b[0][1] = 5;
  b.Randomize();
  std::cout << a.ToWolframString() << '\n' << b.ToWolframString() << '\n';
  auto solve = Solve(l, u, b);
  std::cout << solve.ToWolframString() << '\n';
  std::cout << a * solve;

  // TestMatrixMult();
  // Matrix<double> a(2, 3);
  // a.Randomize(10);
  // Matrix<double> b(3, 2);
  // b.Randomize(10);
  // std::cout << a.ToWolframString() << '\n' << b.ToWolframString() << '\n' <<
  //           (a * b).ToWolframString() << '\n';

  // std::cout << time_measurer.GetDurationString();
  return 0;
}
