#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <cmath>
#include <iostream>

#define FALSE 0
#define TRUE 1

const double EPS = 1e-7;

class S21Matrix {
 private:
  // Attributes
  int _rows, _cols;  // Rows and columns
  double** _matrix;  // Pointer to the memory where the matrix is allocated
  void minor_matrix(const S21Matrix* other, int row, int column);

  void init_matrix(int rows, int cols);
  void copy_matrix(double** matrix);

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  void sum_matrix(const S21Matrix& other);
  bool eq_matrix(const S21Matrix& other);
  void sub_matrix(const S21Matrix& other);
  void mul_number(const double num);
  void mul_matrix(const S21Matrix& other);
  S21Matrix transpose();
  S21Matrix calc_complements();
  double determinant();
  S21Matrix inverse_matrix();

  int getrows();
  int getcols();
  void setrows(int rows);
  void setcols(int cols);

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double& one);
  friend S21Matrix operator*(const double& one, const S21Matrix& two);

  S21Matrix operator=(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other);
  S21Matrix operator=(const S21Matrix& other);
  S21Matrix operator+=(const S21Matrix& other);
  S21Matrix operator-=(const S21Matrix& other);
  S21Matrix operator*=(const S21Matrix& other);
  S21Matrix operator*=(const double& one);
  double& operator()(int row, int col) const;
};

#endif  // SRC_S21_MATRIX_H_
