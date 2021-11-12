#include <iostream>
#include <fstream>

#include "Matrix/matrix.h"
#include "Algebra/algebra.h"
#include "TimeMeasurer/time_measurer.h"
#include "Matrix/diagonal_box_matrix.h"
#include "Matrix/transposed_matrix.h"
#include "Algebra/gmres.h"

void TestMatrixMult() {
  for (int k = 0; k < 1000; k++) {
    Matrix<double> a(100);
    Matrix<double> b(100);
    a.Randomize();
    b.Randomize();
    auto c = a * b;
  }
}

std::vector<double> TestTask1Solvers(int n, int test_count, int seed) {
  std::vector<double> durations;

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
    durations.push_back(time_measurer.GetDuration());
  }
  if (lu_l * lu_u != a) {
    exit(228);
  }

  Matrix<double> ldlt_l(n);
  Matrix<double> ldlt_d(n, 1);
  Matrix<double> ldlt_l_transposed(n, 1);
  {
    TimeMeasurer time_measurer;
    algebra::GetLdlt(a, ldlt_l, ldlt_d);
    ldlt_l_transposed = ldlt_l;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        std::swap(ldlt_l_transposed.At(i, j), ldlt_l_transposed.At(j, i));
      }
    }
    std::cout << "Ldlt prepare: " << time_measurer.GetDurationString() << '\n';
    durations.push_back(time_measurer.GetDuration());
  }
  Matrix<double> ldlt_d_full(n, n);
  for (int i = 0; i < n; i++) {
    ldlt_d_full.At(i, i) = ldlt_d.At(i, 0);
  }
  if (ldlt_l * ldlt_d_full * ldlt_l_transposed != a) {
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
    durations.push_back(time_measurer.GetDuration());
  }
  srand(seed);
  {
    TimeMeasurer time_measurer;
    for (int test = 0; test < test_count; test++) {
      b.Randomize();
      auto ldlt_x = algebra::LdltSolve(ldlt_l, ldlt_d, ldlt_l_transposed, b);
    }
    std::cout << "Ldlt: " << time_measurer.GetDurationString() << '\n';
    durations.push_back(time_measurer.GetDuration());
  }
  return durations;
}

void Task1(int solve_count = 1000) {
  // std::vector<int> sizes{4000, 4400, 4800};
  std::vector<int> sizes{40, 80, 120, 160, 200};
  std::vector<double> lu_prepare;
  std::vector<double> ldlt_prepare;
  std::vector<double> lu_solve;
  std::vector<double> ldlt_solve;
  for (auto size: sizes) {
    auto times = TestTask1Solvers(size, solve_count, 228);
    lu_prepare.push_back(times[0]);
    ldlt_prepare.push_back(times[1]);
    lu_solve.push_back(times[2]);
    ldlt_solve.push_back(times[3]);
  }
  std::ofstream out("../task1.txt");
  for (auto size: sizes) {
    out << size << ' ';
  }
  out << '\n';
  out << "LU Prepare/LDLt Prepare/LU Solve/LDLt solve\n";
  for (auto time: lu_prepare) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: ldlt_prepare) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: lu_solve) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: ldlt_solve) {
    out << time << ' ';
  }
  out << '\n';
}

void Task2() {
  std::vector<int> sizes{500, 1000, 1500, 2000, 2500, 3000, 3500, 4000};
  // std::vector<int> sizes{10};
  std::vector<double> seidel_time;
  std::vector<double> gmres_time;
  std::vector<int> seidel_iters;
  std::vector<int> gmres_iters;
  for (auto size: sizes) {
    DiagonalBoxMatrix<double> a(size);
    Matrix<double> b(size, 1, 1.);
    {
      seidel_iters.push_back(0);
      TimeMeasurer time_measurer;
      auto ans = algebra::GaussSeidelSolve(a, b, 1e-10, 1, 1e3,
                                           &seidel_iters.back());
      seidel_time.push_back(time_measurer.GetDuration());
      // std::cerr << ans << '\n';
    }
    {
      TimeMeasurer time_measurer;
      auto ans = algebra::Gmres<double>(a, b, algebra::EuclideanNorm<double>,
                                        size + 1, 1e-10);
      gmres_time.push_back(time_measurer.GetDuration());
      gmres_iters.push_back(ans.second);
      // std::cerr << ans.first << '\n';
    }
  }
  std::ofstream out("../task2.txt");
  for (auto size: sizes) {
    out << size << ' ';
  }
  out << '\n';
  out << "Gauss Seidel/Gauss Seidel/GMRES/GMRES\n";
  for (auto time: seidel_iters) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: seidel_time) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: gmres_iters) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: gmres_time) {
    out << time << ' ';
  }
  out << '\n';
  for (int i = 0; i < sizes.size(); i++) {
    std::cout << seidel_iters[i] << ' ' << seidel_time[i] << " : " <<
              gmres_iters[i] << ' ' << gmres_time[i] << '\n';
  }
}

void Task3(double min, double max, int seed, double eps) {
  int size = 1000;
  std::vector<double> sparsity{0., 0.5, 0.9, 0.95, 0.99};
  // std::vector<double> sparsity{0.5};
  std::vector<double> gauss_times;
  std::vector<double> gmres_times;
  for (auto sp: sparsity) {
    Matrix<double> a(size);
    algebra::FillRandomNonDegenerate(a, sp, min, max, seed, eps);
    Matrix<double> b(size, 1, 1);
    {
      TimeMeasurer time_measurer;
      auto ans = algebra::GaussSolve(a, b, eps).first;
      gauss_times.push_back(time_measurer.GetDuration());
      if (a * ans != b) {
        std::cerr << "no\n";
      }
    }
    {
      TimeMeasurer time_measurer;
      auto ans = algebra::Gmres<double>(a, b,
                                        algebra::EuclideanNorm<double>,
                                        size + 1, eps).first;
      gmres_times.push_back(time_measurer.GetDuration());
      if (a * ans != b) {
        std::cerr << "no2\n";
      }
    }
    std::cerr << sp << '\n';
  }
  for (int i = 0; i < sparsity.size(); i++) {
    std::cout << gauss_times[i] << " : " << gmres_times[i] << '\n';
  }
  std::ofstream out("../task3.txt");
  for (auto sp: sparsity) {
    out << sp << ' ';
  }
  out << '\n';
  out << "Gauss/GMRES\n";
  for (auto time: gauss_times) {
    out << time << ' ';
  }
  out << '\n';
  for (auto time: gmres_times) {
    out << time << ' ';
  }
  out << '\n';
}

int main() {
  // srand(228);
  srand(time(0));

  // Task1(2000);
  // Task2();
  Task3(-10, 10, 228, 1e-10);

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

  // {
  //   Matrix<double> a{{1, 2, 0},
  //                    {2, 4, 0},
  //                    {3, 9, 0}};
  //   std::cout << a.ToWolframString() << '\n';
  //   std::cout << algebra::GaussSolve(a, Matrix<double>(3, 1)).second;
  // }

  // {
  //   for (int i = 0; i < 100; i++) {
  //     int n = 100;
  //     Matrix<double> a(n);
  //     algebra::FillRandomNonDegenerate(a, 0.5, 1., 10., 228);
  //     // std::cerr << a;
  //     // std::cout << a << a.ToWolframString() << '\n';
  //     Matrix<double> b(n, 1);
  //     b.Randomize();
  //     // std::cout << b << '\n';
  //     // auto ans = algebra::GaussSolve(a, b);
  //     auto ans = algebra::Gmres<double>(a, b, algebra::EuclideanNorm<double>, n + 10);
  //     // std::cout << algebra::EuclideanNorm(a * ans.first - b) << '\n';
  //     if (a * ans.first != b) {
  //       std::cerr << b << a * ans.first;
  //     }
  //     std::cout << ans.second << '\n';
  //     // break;
  //   }
  // }

  // {
  // int n = 10;
  // Matrix<double> a(n);
  // algebra::FillRandomNonDegenerate(a, 0.8, 1., 10., 228);
  // std::cout << a;
  // std::cout << algebra::GaussSolve(a, Matrix<double>(n, 1)).second;
  // }

  // {
  //   for (int i = 0; i < 100; i++) {
  //     Matrix<double> a(3);
  //     a.Randomize(100);
  //     Matrix<double> b(3, 1);
  //     b.Randomize(100);
  //     // std::cout << a << b;
  //     auto [solve, rank] = algebra::GaussSolve(a, b);
  //     // std::cout << solve << a * solve;
  //     if (a * solve != b) {
  //       std::cout << b << a * solve << "===========\n";
  //     }
  //   }
  // }

  // {
  //   int n = 4000;
  //   Matrix<double> b(n, 1);
  //   for (int i = 0; i < n; i++) {
  //     b.At(i, 0) = 1;
  //   }
  //   int iters = 0;
  //   auto a = DiagonalBoxMatrix<double>(n);
  //   auto solve = algebra::GaussSeidelSolve(a, b, 1e-10, 1, 10000, &iters);
  //   std::cout << solve;
  //   std::cout << Matrix<double>::FromAbstract(a) * solve;
  //   std::cout << iters;
  //   return 0;
  // }
}
