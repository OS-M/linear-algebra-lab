#include <iostream>

#include "Matrix/matrix.h"
#include "Algebra/algebra.h"
#include "TimeMeasurer/time_measurer.h"
#include "Matrix/diagonal_box_matrix.h"
#include "Matrix/transposed_matrix.h"

void TestMatrixMult() {
  for (int k = 0; k < 1000; k++) {
    Matrix<double> a(100);
    Matrix<double> b(100);
    a.Randomize();
    b.Randomize();
    auto c = a * b;
  }
}

void TestSolvers(int n, int test_count, int seed) {
  Matrix<double> a(n);
  a.Randomize();
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      a.At(i, j) = a.At(j, i);
    }
  }
  Matrix<double> b(n, 1);

  Matrix<double> lu_l(n);
  Matrix<double> lu_u(n);
  {
    TimeMeasurer time_measurer;
    algebra::GetLu(a, lu_l, lu_u);
    std::cout << "Lu prepare: " << time_measurer.GetDurationString() << '\n';
  }
  if (lu_l * lu_u != a) {
    exit(228);
  }

  Matrix<double> ldlt_l(n);
  Matrix<double> ldlt_d(n, 1);
  {
    TimeMeasurer time_measurer;
    algebra::GetLdlt(a, ldlt_l, ldlt_d);
    std::cout << "Ldlt prepare: " << time_measurer.GetDurationString() << '\n';
  }

  srand(seed);
  {
    TimeMeasurer time_measurer;
    for (int test = 0; test < test_count; test++) {
      b.Randomize();
      auto lu_x = algebra::LuSolve(lu_l, lu_u, b);
    }
    std::cout << "Lu: " << time_measurer.GetDurationString() << '\n';
  }
  srand(seed);
  {
    TimeMeasurer time_measurer;
    for (int test = 0; test < test_count; test++) {
      b.Randomize();
      auto ldlt_x = algebra::LdltSolve(ldlt_l, ldlt_d, b);
    }
    std::cout << "Ldlt: " << time_measurer.GetDurationString() << '\n';
  }
}

int main() {
  TestSolvers(1000, 5000, time(0));

  int n = 3;
  Matrix<double> b(n, 1);
  for (int i = 0; i < n; i++) {
    b.At(i, 0) = 1;
  }
  auto a = DiagonalBoxMatrix<double>(n);
  auto solve = algebra::GaussSeidelSolve(a, b, 1e-10);
  std::cout << solve << '\n';
  std::cout << Matrix<double>::FromAbstract(a) * solve << '\n';
  return 0;
}
