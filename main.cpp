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
  Matrix<double> ldlt_d_full(n, n);
  for (int i = 0; i < n; i++) {
    ldlt_d_full.At(i, i) = ldlt_d.At(i, 0);
  }
  if (ldlt_l * ldlt_d_full * ldlt_l.Transposed() != a) {
    exit(229);
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
  srand(228);
  // TestSolvers(100, 5000, time(0));

  // {
  //   Matrix<double> a{{0.593086, -0.633130, 0.193700},
  //                    {-0.633130, 0.559382, -0.108334},
  //                    {0.193700, -0.108334 ,-0.800050}};
  //   Matrix<double> l(3);
  //   Matrix<double> d(3, 1);
  //   algebra::GetLdlt(a, l, d);
  //   Matrix<double> ldlt_d_full(3);
  //   for (int i = 0; i < 3; i++) {
  //     ldlt_d_full.At(i, i) = d.At(i, 0);
  //   }
  //   std::cout << l << ldlt_d_full << l * ldlt_d_full * l.Transposed();
  // }

  // {
  //   Matrix<double> a{{2, -1, 0},
  //                    {-1, 1, 4},
  //                    {1, 2, 3}};
  //   Matrix<double> b{{0}, {13}, {14}};
  //   std::cout << a << b;
  //   auto solve = algebra::GaussSolve(a, b);
  //   std::cout << solve << a * solve;
  // }

  {
    for (int i = 0; i < 100; i++) {
      Matrix<double> a(3);
      a.Randomize(100);
      Matrix<double> b(3, 1);
      b.Randomize(100);
      // std::cout << a << b;
      auto solve = algebra::GaussSolve(a, b);
      // std::cout << solve << a * solve;
      if (a * solve != b) {
        std::cout << b << a * solve << "===========\n";
      }
    }
  }


  return 0;
  int n = 3;
  Matrix<double> b(n, 1);
  for (int i = 0; i < n; i++) {
    b.At(i, 0) = 1;
  }
  auto a = DiagonalBoxMatrix<double>(n);
  auto solve = algebra::GaussSeidelSolve(a, b, 1e-10);
  std::cout << solve;
  std::cout << Matrix<double>::FromAbstract(a) * solve;
  return 0;
}
