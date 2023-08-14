#include "s21_matrix_oop.h"

// Иницилизация матрицы
void S21Matrix::init_matrix(int rows, int cols) {
  _rows = rows;
  _cols = cols;
  _matrix = new double*[rows];
  for (int i = 0; i < rows; i++) {
    _matrix[i] = new double[cols]();
  }
}

// Копирование матрицы
void S21Matrix::copy_matrix(double** matrix) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = matrix[i][j];
    }
  }
}

// Записать в _rows
void S21Matrix::setrows(int rows) { _rows = rows; }

// Записать в _cols
void S21Matrix::setcols(int cols) { _cols = cols; }

// Считать _rows
int S21Matrix::getrows() { return _rows; }

// Считать _cols
int S21Matrix::getcols() { return _cols; }

// Базовый конструктор, инициализирующий матрицу некоторой заранее заданной
// размерностью
S21Matrix::S21Matrix() : _rows(3), _cols(3) { init_matrix(_rows, _cols); }

// Параметризированный конструктор с количеством строк и столбцов
S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1) {
    throw std::out_of_range("Number of rows or columns is less than zero");
  }
  _rows = rows;
  _cols = cols;
  init_matrix(rows, cols);
}

// Деструктор
S21Matrix::~S21Matrix() {
  if (_matrix) {
    for (int i = 0; i < _rows; i++) {
      delete[] _matrix[i];
    }
    delete[] _matrix;
  }
}

// Конструктор копирования
S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other._rows, other._cols) {
  _rows = other._rows;
  _cols = other._cols;
  copy_matrix(other._matrix);
}

// Конструктор переноса
S21Matrix::S21Matrix(S21Matrix&& other)
    : _rows(other._rows), _cols(other._cols), _matrix(other._matrix) {
  _rows = other._rows;
  _cols = other._cols;
  copy_matrix(other._matrix);
  other._rows = 0;
  other._cols = 0;
  other._matrix = nullptr;
}

// Сложение двух матриц
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->sum_matrix(other);
  return *temp;
}

// Вычитание одной матрицы из другой
S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->sub_matrix(other);
  return *temp;
}

// Умножение матриц
S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->mul_matrix(other);
  return *temp;
}

// Умножение матрицы на число
S21Matrix S21Matrix::operator*(const double& one) {
  S21Matrix* temp(this);
  temp->mul_number(one);
  return *temp;
}

// Умножение числа на матрицу
S21Matrix operator*(const double& one, const S21Matrix& two) {
  S21Matrix res(two);
  res.mul_number(one);
  return res;
}

// Присвоение сложения (sum_matrix)
S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->sum_matrix(other);
  return *temp;
}

// Присвоение разности (sub_matrix)
S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->sub_matrix(other);
  return *temp;
}

// Присвоение умножения (mul_matrix)
S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("sizes of both matrixes are not equal");
  }
  S21Matrix* temp(this);
  temp->mul_matrix(other);
  return *temp;
}

// Присвоение умножения (mul_number)
S21Matrix S21Matrix::operator*=(const double& other) {
  S21Matrix* temp(this);
  temp->mul_number(other);
  return *temp;
}

// Присвоение матрице значений другой матрицы
S21Matrix S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    this->~S21Matrix();
    init_matrix(other._rows, other._cols);
    copy_matrix(other._matrix);
  }
  return *this;
}

// Проверка на равенство матриц (eq_matrix)
bool S21Matrix::operator==(const S21Matrix& other) { return eq_matrix(other); }

// Индексация по элементам матрицы (строка, колонка)
double& S21Matrix::operator()(int row, int col) const {
  if (_rows < row || _cols < col || col < 0 || row < 0) {
    throw std::out_of_range("Number of rows or columns is less than zero");
  }
  return _matrix[row][col];
}

// Создает новую транспонированную матрицу из текущей и возвращает ее
S21Matrix S21Matrix::transpose() {
  S21Matrix res(_cols, _rows);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      res._matrix[j][i] = _matrix[i][j];
    }
  }
  return res;
}

// Минор матрицы
void S21Matrix::minor_matrix(const S21Matrix* other, int rows, int cols) {
  int x = 0, y = 0;
  for (int i = 0; i < other->_rows; i++) {
    if (i != rows) {
      for (int j = 0; j < other->_cols; j++) {
        if (j != cols) {
          _matrix[x][y] = other->_matrix[i][j];
          y++;
        }
      }
      x++;
      y = 0;
    }
  }
}

// Вычисляет матрицу алгебраических дополнений текущей матрицы и возвращает ее
S21Matrix S21Matrix::calc_complements() {
  if (_rows != _cols || _rows < 2) {
    throw std::invalid_argument("Matrix is not square");
  }
  S21Matrix res(_rows, _cols);
  S21Matrix tmp(_rows - 1, _cols - 1);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      tmp.minor_matrix(this, i, j);
      res._matrix[i][j] = pow(-1, i + j) * tmp.determinant();
    }
  }
  return res;
}

// Вычисляет и возвращает обратную матрицу
S21Matrix S21Matrix::inverse_matrix() {
  double deter = determinant();
  if (fabs(deter) < EPS) {
    throw std::invalid_argument("Determinant is zero");
  }

  S21Matrix tmp = calc_complements();
  S21Matrix res = tmp.transpose();
  res.mul_number(1.0 / deter);
  return res;
}

// Вычисляет и возвращает определитель текущей матрицы
double S21Matrix::determinant() {
  if (_rows != _cols) {
    throw std::invalid_argument("Matrix is not square");
  }

  double res = 0;
  if (_rows == 1) {
    res = _matrix[0][0];
  } else if (_rows == 2) {
    res = _matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1];
  } else {
    S21Matrix tmp(_rows - 1, _cols - 1);
    int sign = 1;
    for (int i = 0; i < _rows; i++) {
      tmp.minor_matrix(this, i, 0);
      res += sign * _matrix[i][0] * tmp.determinant();
      sign *= -1;
    }
  }
  return res;
}

// Проверяет матрицы на равенство между собой
bool S21Matrix::eq_matrix(const S21Matrix& other) {
  bool tmp = TRUE;
  if (_rows == other._rows && _cols == other._cols) {
    for (int i = 0; i < _rows && tmp == TRUE; i++) {
      for (int j = 0; j < _cols && tmp == TRUE; j++) {
        if (fabs(_matrix[i][j] - other._matrix[i][j]) > EPS) {
          tmp = FALSE;
        }
      }
    }
  } else {
    tmp = FALSE;
  }
  return tmp;
}

// Умножает текущую матрицу на вторую
void S21Matrix::mul_matrix(const S21Matrix& other) {
  if (other._rows != _cols)
    throw std::invalid_argument(
        "Number of columns in first matrix are not equal to number of "
        "rows in "
        "second matrix");
  S21Matrix res(_rows, other._cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int k = 0; k < other._rows; k++) {
        res._matrix[i][j] += _matrix[i][k] * other._matrix[k][j];
      }
    }
  }
  *this = res;
}

// Умножает текущую матрицу на число
void S21Matrix::mul_number(const double num) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] * num;
    }
  }
}

// Прибавляет вторую матрицы к текущей
void S21Matrix::sum_matrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("Sizes of both Matrixes are not equal");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] + other._matrix[i][j];
    }
  }
}

// Вычитает из текущей матрицы другую
void S21Matrix::sub_matrix(const S21Matrix& other) {
  if (_rows != other._rows || _cols != other._cols) {
    throw std::invalid_argument("Sizes of both Matrixes are not equal");
  }
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] - other._matrix[i][j];
    }
  }
}
