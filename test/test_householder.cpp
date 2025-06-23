#ifndef USING_IMPORT_STD_MOD
#  include <print>
#  include <vector>
#else
import std;
#endif

import math;

void test_householder()
{
  jf::matrix::dmat f1(3, 3);
  f1.set_value({2, 5, 1, -4, 5, 0, -1, 2, 1});
  auto [Q, R] = f1.houseHolder();

  std::print("value of Q\n{}\nvalue of R\n{}\n\n", Q, R);
}
void test_eigenvalues()
{
  jf::matrix::dmat f2(3, 3);
  f2.set_value({2, 0, 1, -1, 4, -1, -1, 2, 0});
  auto eigs1 = f2.eigenvalues();

  std::print("the eigenvalues of f2 are\n");
  for (const auto& v : eigs1)
    std::print("{} | ", v);
  std::print("\n");

  jf::matrix::dmat f3(4, 4);
  f3.set_value({6, 2, 1, 3, 2, 5, 2, 1, 1, 2, 4, 2, 3, 1, 2, 7});
  auto eigs2 = f3.eigenvalues();

  std::print("the eigenvalues of f3 are\n");
  for (const auto& v : eigs2)
    std::print("{} | ", v);

  jf::matrix::dmat f4(2, 2);
  f4.set_value({-6, 3, 4, 5});
  auto eigs3 = f4.eigenvalues();

  std::print("\nthe eigenvalues of f4 are\n");
  for (const auto& v : eigs3)
    std::print("{} | ", v);
  std::print("\nfim\n");
}

void test_eigenvectors()
{
  jf::matrix::dmat mat1(3, 3);
  mat1.set_value({2, 0, 1, -1, 4, -1, -1, 2, 0});
  auto eigs1 = mat1.eigenvalues();

  std::print("mat1 = \n{}\nthe eigenvalues of mat1 are\n", mat1);
  for (const auto& v : eigs1)
    std::print("{} | ", v);
  std::print("\n");

  auto eigenV1 = mat1.eigenvectors(eigs1);

  std::print(
      "the eigenvectors of mat1 are (one eigenvector in each column)\n{}\n",
      eigenV1);

  jf::matrix::dmat mat2(4, 4);
  mat2.set_value({6, 2, 1, 3, 2, 5, 2, 1, 1, 2, 4, 2, 3, 1, 2, 7});
  auto eigs2 = mat2.eigenvalues();

  std::print("mat2 = \n{}\nthe eigenvalues of mat2 are\n", mat2);
  for (const auto& v : eigs2)
    std::print("{} | ", v);
  std::print("\n");

  auto eigenV2 = mat2.eigenvectors(eigs2);

  std::print(
      "the eigenvectors of mat2 are (one eigenvector in each column)\n{}\n",
      eigenV2);

  auto idlamb = mat2.identity();
  for (std::size_t i {}; i < idlamb.rows(); ++i) {
    idlamb(i, i) = eigs2[i];
  }
  auto AV = mat2 * eigenV2;
  auto VL = eigenV2 * idlamb;

  std::print("Test AV = VL\nAV\n{}\nVL\n{}\n", AV, VL);

  jf::matrix::dmat mat3(2, 2);
  mat3.set_value({-6, 3, 4, 5});
  auto eigs3 = mat3.eigenvalues();

  std::print("mat3 = \n{}\nthe eigenvalues of mat3 are\n", mat3);
  for (const auto& v : eigs3)
    std::print("{} | ", v);
  std::print("\n");

  auto eigenV3 = mat3.eigenvectors(eigs3[0]);

  std::print("the eigenvector for elgenvalue [{}] of mat3 is\n{}\n",
             eigs3[0],
             eigenV3);
}

int main()
{
  // test_householder();
  // test_eigenvalues();
  test_eigenvectors();
  return 0;
}
