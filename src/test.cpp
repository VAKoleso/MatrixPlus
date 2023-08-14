#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(test_move_constructor, test_move_constructor1) {
  S21Matrix expected(1, 2);
  expected(0, 0) = 5;
  expected(0, 1) = 7;

  S21Matrix matrixA(1, 2);
  matrixA(0, 0) = 5;
  matrixA(0, 1) = 7;

  S21Matrix matrixB(std::move(matrixA));

  ASSERT_EQ(matrixB.getcols(), expected.getcols());

  for (int i = 0; i < matrixB.getrows(); i++) {
    for (int j = 0; j < matrixB.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixB(i, j), EPS);
    }
  }
}

TEST(test_sum, test_sum1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 10.2;
  expected(0, 1) = 20.4;
  expected(1, 0) = 30.6;
  expected(1, 1) = 40.8;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 5.1;
  matrixA(0, 1) = 10.2;
  matrixA(1, 0) = 15.3;
  matrixA(1, 1) = 20.4;

  S21Matrix matrixB(matrixA);

  matrixA.sum_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_sum, test_sum2) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 0;
  expected(0, 1) = 0;
  expected(1, 0) = 0;
  expected(1, 1) = 0;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = -1;
  matrixA(0, 1) = -2;
  matrixA(1, 0) = -3;
  matrixA(1, 1) = -4;

  S21Matrix matrixB(2, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;

  matrixA.sum_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_sum, test_sum3) {
  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 5;
  matrixA(1, 0) = 5;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;
  matrixB(2, 0) = 5;
  matrixB(2, 1) = 6;

  ASSERT_THROW(matrixA.sum_matrix(matrixB), std::invalid_argument);
}

TEST(test_sum, test_sum4) {
  S21Matrix expected(3, 2);
  expected(0, 0) = 5;
  expected(0, 1) = 5;
  expected(1, 0) = 5;
  expected(1, 1) = 5;
  expected(2, 0) = 5;
  expected(2, 1) = 5;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 0;
  matrixA(0, 1) = 1;
  matrixA(1, 0) = 2;
  matrixA(1, 1) = 3;
  matrixA(2, 0) = 4;
  matrixA(2, 1) = 5;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 5;
  matrixB(0, 1) = 4;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 2;
  matrixB(2, 0) = 1;
  matrixB(2, 1) = 0;

  matrixA.sum_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_sub, test_sub1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 0;
  expected(0, 1) = 0;
  expected(1, 0) = 0;
  expected(1, 1) = 0;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 10;
  matrixA(0, 1) = 20;
  matrixA(1, 0) = 30;
  matrixA(1, 1) = 40;

  S21Matrix matrixB(matrixA);
  matrixA.sub_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_sub, test_sub2) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 3;
  expected(0, 1) = 1;
  expected(1, 0) = -1;
  expected(1, 1) = -3;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 5;
  matrixA(0, 1) = 4;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 2;

  S21Matrix matrixB(2, 2);
  matrixB(0, 0) = 2;
  matrixB(0, 1) = 3;
  matrixB(1, 0) = 4;
  matrixB(1, 1) = 5;

  matrixA.sub_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_sub, test_sub3) {
  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1;
  matrixA(1, 0) = 2;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 3;
  matrixB(0, 1) = 4;
  matrixB(1, 0) = 5;
  matrixB(1, 1) = 6;
  matrixB(2, 0) = 7;
  matrixB(2, 1) = 8;

  ASSERT_THROW(matrixA.sub_matrix(matrixB), std::invalid_argument);
}

TEST(test_sub, test_sub4) {
  S21Matrix expected(3, 2);
  expected(0, 0) = -9;
  expected(0, 1) = -8;
  expected(1, 0) = -7;
  expected(1, 1) = -26;
  expected(2, 0) = -25;
  expected(2, 1) = -24;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;
  matrixA(2, 0) = 5;
  matrixA(2, 1) = 6;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 10;
  matrixB(0, 1) = 10;
  matrixB(1, 0) = 10;
  matrixB(1, 1) = 30;
  matrixB(2, 0) = 30;
  matrixB(2, 1) = 30;

  matrixA.sub_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_mul_num, test_mul_num1) {
  S21Matrix expected(3, 2);
  expected(0, 0) = 4;
  expected(0, 1) = 8;
  expected(1, 0) = 12;
  expected(1, 1) = 16;
  expected(2, 0) = 20;
  expected(2, 1) = 24;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;
  matrixA(2, 0) = 5;
  matrixA(2, 1) = 6;

  matrixA.mul_number(4);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_mul_num, test_mul_num2) {
  S21Matrix expected(3, 4);
  expected(0, 0) = -1;
  expected(0, 1) = -2;
  expected(0, 2) = -3;
  expected(0, 3) = -4;

  expected(1, 0) = -5;
  expected(1, 1) = -6;
  expected(1, 2) = -7;
  expected(1, 3) = -8;

  expected(2, 0) = -9;
  expected(2, 1) = -10;
  expected(2, 2) = -11;
  expected(2, 3) = -12;

  S21Matrix matrixA(3, 4);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(0, 2) = 3;
  matrixA(0, 3) = 4;

  matrixA(1, 0) = 5;
  matrixA(1, 1) = 6;
  matrixA(1, 2) = 7;
  matrixA(1, 3) = 8;

  matrixA(2, 0) = 9;
  matrixA(2, 1) = 10;
  matrixA(2, 2) = 11;
  matrixA(2, 3) = 12;

  matrixA.mul_number(-1);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_mul, test_mul1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 7;
  expected(0, 1) = 10;
  expected(1, 0) = 15;
  expected(1, 1) = 22;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;

  S21Matrix matrixB(2, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;

  matrixA.mul_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_mul, test_mul2) {
  S21Matrix expected(1, 3);
  expected(0, 0) = 29;
  expected(0, 1) = 42.2;
  expected(0, 2) = 5;

  S21Matrix matrixA(1, 2);
  matrixA(0, 0) = 5;
  matrixA(0, 1) = 1;

  S21Matrix matrixB(2, 3);
  matrixB(0, 0) = 4;
  matrixB(0, 1) = 8;
  matrixB(0, 2) = -2;
  matrixB(1, 0) = 9;
  matrixB(1, 1) = 2.2;
  matrixB(1, 2) = 15;

  matrixA.mul_matrix(matrixB);

  ASSERT_EQ(matrixA.getcols(), expected.getcols());

  for (int i = 0; i < matrixA.getrows(); i++) {
    for (int j = 0; j < matrixA.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), matrixA(i, j), EPS);
    }
  }
}

TEST(test_mul, test_mul3) {
  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1.1;
  matrixA(1, 0) = 3;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;
  matrixB(2, 0) = 4;
  matrixB(2, 1) = 4;

  ASSERT_THROW(matrixA.mul_matrix(matrixB), std::invalid_argument);
}

TEST(test_eq, test_eq1) {
  bool expected = true;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1.1;
  matrixA(0, 1) = 2.2;
  matrixA(1, 0) = 3.3;
  matrixA(1, 1) = 4.4;
  matrixA(2, 0) = 5.5;
  matrixA(2, 1) = 6.6;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 1.1;
  matrixB(0, 1) = 2.2;
  matrixB(1, 0) = 3.3;
  matrixB(1, 1) = 4.4;
  matrixB(2, 0) = 5.5;
  matrixB(2, 1) = 6.6;

  ASSERT_EQ(expected, matrixA.eq_matrix(matrixB));
}

TEST(test_eq, test_eq2) {
  bool expected = false;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 0;
  matrixA(0, 1) = 1;
  matrixA(1, 0) = 2;
  matrixA(1, 1) = 3;
  matrixA(2, 0) = 4;
  matrixA(2, 1) = 5;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;
  matrixB(2, 0) = 5;
  matrixB(2, 1) = 6;

  ASSERT_EQ(expected, matrixA.eq_matrix(matrixB));
}

TEST(test_eq, test_eq3) {
  bool expected = false;

  S21Matrix matrixA(3, 2);
  S21Matrix matrixB(3, 3);

  ASSERT_EQ(expected, matrixA.eq_matrix(matrixB));
}

TEST(test_transpose, test_transpose1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 1;
  expected(0, 1) = 3;
  expected(1, 0) = 2;
  expected(1, 1) = 4;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;

  S21Matrix actual = matrixA.transpose();

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_determinant, test_determinant1) {
  double expected = -40;

  S21Matrix matrixA(3, 3);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(0, 2) = 3;

  matrixA(1, 0) = 0;
  matrixA(1, 1) = 4;
  matrixA(1, 2) = 2;

  matrixA(2, 0) = 5;
  matrixA(2, 1) = 2;
  matrixA(2, 2) = 1;

  ASSERT_EQ(expected, matrixA.determinant());
}

TEST(test_determinant, test_determinant2) {
  double expected = -7900;

  S21Matrix matrixA(5, 5);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(0, 2) = 3;
  matrixA(0, 3) = 4;
  matrixA(0, 4) = 3;

  matrixA(1, 0) = 0;
  matrixA(1, 1) = 4;
  matrixA(1, 2) = 2;
  matrixA(1, 3) = 3;
  matrixA(1, 4) = 13;

  matrixA(2, 0) = 5;
  matrixA(2, 1) = 2;
  matrixA(2, 2) = 1;
  matrixA(2, 3) = 3;
  matrixA(2, 4) = 3;

  matrixA(3, 0) = 5;
  matrixA(3, 1) = 12;
  matrixA(3, 2) = 1;
  matrixA(3, 3) = 3;
  matrixA(3, 4) = 3;

  matrixA(4, 0) = 5;
  matrixA(4, 1) = 2;
  matrixA(4, 2) = 1;
  matrixA(4, 3) = 8;
  matrixA(4, 4) = 3;

  ASSERT_EQ(matrixA.determinant(), expected);
}

TEST(test_determinant, test_determinant3) {
  double expected = 1;

  S21Matrix matrixA(1, 1);
  matrixA(0, 0) = 1;

  ASSERT_EQ(matrixA.determinant(), expected);
}

TEST(test_determinant, test_determinant4) {
  S21Matrix matrixA(1, 3);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 1;
  matrixA(0, 2) = 1;

  ASSERT_THROW(matrixA.determinant(), std::invalid_argument);
}

TEST(test_calc_complements, test_calc_complements1) {
  S21Matrix expected;
  expected(0, 0) = 0;
  expected(0, 1) = 10;
  expected(0, 2) = -20;

  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;

  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  S21Matrix matrixA(3, 3);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(0, 2) = 3;

  matrixA(1, 0) = 0;
  matrixA(1, 1) = 4;
  matrixA(1, 2) = 2;

  matrixA(2, 0) = 5;
  matrixA(2, 1) = 2;
  matrixA(2, 2) = 1;

  S21Matrix actual = matrixA.calc_complements();

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_calc_complements, test_calc_complements2) {
  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 0;
  matrixA(1, 1) = 4;
  matrixA(2, 0) = 5;
  matrixA(2, 1) = 2;

  ASSERT_THROW(matrixA.calc_complements(), std::invalid_argument);
}

TEST(test_inverse_matrix, test_inverse_matrix1) {
  S21Matrix expected(3, 3);
  expected(0, 0) = 1;
  expected(0, 1) = -1;
  expected(0, 2) = 1;

  expected(1, 0) = -38;
  expected(1, 1) = 41;
  expected(1, 2) = -34;

  expected(2, 0) = 27;
  expected(2, 1) = -29;
  expected(2, 2) = 24;

  S21Matrix matrixA(3, 3);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(0, 2) = 7;

  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;
  matrixA(1, 2) = 4;

  matrixA(2, 0) = 5;
  matrixA(2, 1) = -2;
  matrixA(2, 2) = -3;

  S21Matrix actual = matrixA.inverse_matrix();

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_inverse_matrix, test_inverse_matrix2) {
  S21Matrix matrixA(3, 3);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(0, 2) = 3;

  matrixA(1, 0) = 1;
  matrixA(1, 1) = 2;
  matrixA(1, 2) = 3;

  matrixA(2, 0) = 1;
  matrixA(2, 1) = 2;
  matrixA(2, 2) = 3;

  ASSERT_THROW(matrixA.inverse_matrix(), std::invalid_argument);
}

TEST(test_operators, test_mul_operator_1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 34;
  expected(0, 1) = 25;
  expected(1, 0) = 30;
  expected(1, 1) = 39;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix matrixB(matrixA);
  S21Matrix actual = matrixA * matrixB;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_mul_operator2) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(1, 0) = 12;
  expected(1, 1) = 6;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix actual = matrixA * 2;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_mul_operator3) {
  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix matrixB(1, 1);
  matrixB(0, 0) = 2;

  ASSERT_THROW(matrixA * matrixB, std::invalid_argument);
}

TEST(test_operators, test_mul_operator4) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(1, 0) = 12;
  expected(1, 1) = 6;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix actual = 2 * matrixA;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_sum_operator1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 4;
  expected(0, 1) = 10;
  expected(1, 0) = 4;
  expected(1, 1) = 20;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 2;
  matrixA(1, 1) = 10;

  S21Matrix matrixB(matrixA);
  S21Matrix actual = matrixA + matrixB;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_sum_operator2) {
  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix matrixB(1, 2);
  matrixB(0, 0) = 2;
  matrixB(0, 1) = 5;

  ASSERT_THROW(matrixA + matrixB, std::invalid_argument);
}

TEST(test_operators, test_sub_operator1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 0;
  expected(0, 1) = 0;
  expected(1, 0) = 0;
  expected(1, 1) = 0;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 20;
  matrixA(0, 1) = 50;
  matrixA(1, 0) = 20;
  matrixA(1, 1) = 60;

  S21Matrix matrixB(matrixA);
  S21Matrix actual = matrixA - matrixB;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_sub_operator2) {
  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 5;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 3;

  S21Matrix matrixB(1, 2);
  matrixB(0, 0) = 2;
  matrixB(0, 1) = 5;

  ASSERT_THROW(matrixA - matrixB, std::invalid_argument);
}

TEST(test_operators, test_eq1) {
  bool expected = true;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;
  matrixA(2, 0) = 4;
  matrixA(2, 1) = 4;

  S21Matrix matrixB(3, 2);
  matrixB(0, 0) = 1;
  matrixB(0, 1) = 2;
  matrixB(1, 0) = 3;
  matrixB(1, 1) = 4;
  matrixB(2, 0) = 4;
  matrixB(2, 1) = 4;

  bool actual = matrixA == matrixB;

  ASSERT_EQ(expected, actual);
}

TEST(test_operators, test_assignment1) {
  S21Matrix expected(3, 2);
  expected(0, 0) = 1;
  expected(0, 1) = 2;
  expected(1, 0) = 3;
  expected(1, 1) = 4;
  expected(2, 0) = 4;
  expected(2, 1) = 4;

  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 3;
  matrixA(1, 1) = 4;
  matrixA(2, 0) = 4;
  matrixA(2, 1) = 4;

  S21Matrix actual = matrixA;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_addition_assignment1) {
  S21Matrix expected(3, 2);
  expected(0, 0) = 2;
  expected(0, 1) = 4;
  expected(1, 0) = 6;
  expected(1, 1) = 8;
  expected(2, 0) = 8;
  expected(2, 1) = 8;

  S21Matrix actual(3, 2);
  actual(0, 0) = 1;
  actual(0, 1) = 2;
  actual(1, 0) = 3;
  actual(1, 1) = 4;
  actual(2, 0) = 4;
  actual(2, 1) = 4;

  actual += actual;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_addition_assignment2) {
  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 4;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 8;
  matrixA(2, 0) = 8;
  matrixA(2, 1) = 8;

  S21Matrix matrixB(1, 1);
  matrixB(0, 0) = 1;

  ASSERT_THROW(matrixA += matrixB, std::invalid_argument);
}

TEST(test_operators, test_substraction_assignment2) {
  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 4;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 8;
  matrixA(2, 0) = 8;
  matrixA(2, 1) = 8;

  S21Matrix matrixB(1, 1);
  matrixB(0, 0) = 1;

  ASSERT_THROW(matrixA -= matrixB, std::invalid_argument);
}

TEST(test_operators, test_multiplication_assignment1) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 9;
  expected(0, 1) = 10;
  expected(1, 0) = 20;
  expected(1, 1) = 24;

  S21Matrix actual(2, 2);
  actual(0, 0) = 1;
  actual(0, 1) = 2;
  actual(1, 0) = 4;
  actual(1, 1) = 4;

  actual *= actual;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_multiplication_assignment2) {
  S21Matrix expected(2, 2);
  expected(0, 0) = 2;
  expected(0, 1) = 4;
  expected(1, 0) = 8;
  expected(1, 1) = 8;

  S21Matrix actual(2, 2);
  actual(0, 0) = 1;
  actual(0, 1) = 2;
  actual(1, 0) = 4;
  actual(1, 1) = 4;

  actual *= 2;

  ASSERT_EQ(actual.getcols(), expected.getcols());

  for (int i = 0; i < actual.getrows(); i++) {
    for (int j = 0; j < actual.getcols(); j++) {
      ASSERT_NEAR(expected(i, j), actual(i, j), EPS);
    }
  }
}

TEST(test_operators, test_multiplication_assignment3) {
  S21Matrix matrixA(3, 2);
  matrixA(0, 0) = 2;
  matrixA(0, 1) = 4;
  matrixA(1, 0) = 6;
  matrixA(1, 1) = 8;
  matrixA(2, 0) = 8;
  matrixA(2, 1) = 8;

  S21Matrix matrixB(1, 1);
  matrixB(0, 0) = 1;

  ASSERT_THROW(matrixA *= matrixB, std::invalid_argument);
}

TEST(test_operators, test_parathesis1) {
  double expected = 2;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 4;
  matrixA(1, 1) = 4;

  double actual = matrixA(0, 1);

  ASSERT_EQ(actual, expected);
}

TEST(test_operators, test_parathesis2) {
  double expected = 5;

  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 4;
  matrixA(1, 1) = 4;

  matrixA(0, 1) = 5;
  double actual = matrixA(0, 1);

  ASSERT_EQ(actual, expected);
}

TEST(test_operators, test_parathesis3) {
  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 4;
  matrixA(1, 1) = 4;

  ASSERT_THROW(matrixA(0, 5) = 2, std::out_of_range);
}

TEST(test_operators, test_parathesis4) {
  S21Matrix matrixA(2, 2);
  matrixA(0, 0) = 1;
  matrixA(0, 1) = 2;
  matrixA(1, 0) = 4;
  matrixA(1, 1) = 4;

  ASSERT_THROW(matrixA(0, 5), std::out_of_range);
}

TEST(test_getters, test_getters1) {
  double expected = 2;

  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1;
  matrixA(1, 0) = 4;

  double actual = matrixA.getrows();

  ASSERT_EQ(actual, expected);
}

TEST(test_getters, test_getters2) {
  double expected = 1;

  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1;
  matrixA(1, 0) = 4;

  double actual = matrixA.getcols();

  ASSERT_EQ(actual, expected);
}

TEST(test_setters, test_setters1) {
  double expected = 3;

  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1;
  matrixA(1, 0) = 4;

  matrixA.setrows(3);
  double actual = matrixA.getrows();
  matrixA.setrows(2);

  ASSERT_EQ(actual, expected);
}

TEST(test_setters, test_setters2) {
  double expected = 3;

  S21Matrix matrixA(2, 1);
  matrixA(0, 0) = 1;
  matrixA(1, 0) = 4;

  matrixA.setcols(3);
  double actual = matrixA.getcols();
  matrixA.setcols(1);

  ASSERT_EQ(actual, expected);
}

TEST(test_constructor, test_constructor_throw) {
  ASSERT_THROW(S21Matrix matrixA(-1, 2), std::out_of_range);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
