#include <iostream>

#include "Matrix/matrix.h"
#include "Algebra/algebra.h"
#include "TimeMeasurer/time_measurer.h"
#include "Matrix/diagonal_box_matrix.h"

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
      a[i][j] = a[j][i];
    }
  }
  Matrix<double> b(n, 1);

  Matrix<double> lu_l(n);
  Matrix<double> lu_u(n);
  algebra::GetLu(a, lu_l, lu_u);
  std::cout << "Ready lu\n";

  Matrix<double> ldlt_l(n);
  Matrix<double> ldlt_d(n, 1);
  algebra::GetLdlt(a, ldlt_l, ldlt_d);
  std::cout << "Ready ldlt\n";

  srand(seed);
  {
    TimeMeasurer time_measurer;
    for (int test = 0; test < test_count; test++) {
      b.Randomize();
      auto lu_x = algebra::LuSolve(lu_l, lu_u, b);
      // if (a * lu_x != b) {
      //   std::cerr << a.ToWolframString() << '\n' << b.ToWolframString() << '\n';
      //   std::cerr << "lu\n";
      //   return;
      // }
    }
    std::cout << "Lu: " << time_measurer.GetDurationString() << '\n';
  }
  srand(seed);
  {
    TimeMeasurer time_measurer;
    for (int test = 0; test < test_count; test++) {
      b.Randomize();
      auto ldlt_x = algebra::LdltSolve(ldlt_l, ldlt_d, b);
      // if (a * ldlt_x != b) {
      //   std::cerr << a.ToWolframString() << '\n' << b.ToWolframString() << '\n';
      //   std::cerr << "ldlt\n";
      //   return;
      // }
    }
    std::cout << "Ldlt: " << time_measurer.GetDurationString() << '\n';
  }
}

int main() {
  srand(time(0));
  TimeMeasurer time_measurer;

  std::cout << DiagonalBoxMatrix<double>(6);
  TestSolvers(200, 1000, time(0));

  // int n = 3;
  // Matrix<double> a(n);
  // a.Randomize();
  // for (int i = 0; i < n; i++) {
  //   for (int j = i + 1; j < n; j++) {
  //     a[i][j] = a[j][i];
  //   }
  // }
  // Matrix<double> b(n, 1);
  // b.Randomize();
  // std::cout << "A = \n" << a << '\n';
  // std::cout << "b = \n" << b << '\n';

  // Matrix<double> l(n);
  // Matrix<double> u(n);
  // GetLu(a, l, u);
  // auto solve = LuSolve(l, u, b);
  // std::cout << solve.ToWolframString() << '\n';
  // std::cout << a * solve;

  // Matrix<double> l(n);
  // Matrix<double> d(n, 1);
  // Matrix<double> d_m(n);
  // algebra::GetLdlt(a, l, d);
  // for (int i = 0; i < n; i++) {
  //   d_m[i][i] = d[i];
  // }
  // // std::cout << l << '\n' << d_m << '\n' << l * d_m * l.Transposed() << '\n';
  // auto solve = algebra::LdltSolve(l, d, b);
  // std::cout << a * solve << '\n';

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
