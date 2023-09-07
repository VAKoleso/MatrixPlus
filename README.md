# MatrixPlus

>

Реализация библиотеки s21_matrix_oop.h


## Содержание

0. [Введение](#введение)
1. [Реализация библиотеки ](#реализация-библиотеки)
2. [Информация](#информация)
3. [Операции над матрицами](#операции-над-матрицами)


![Matrix](misc/images/Matrix.png)



## Введение

В данном проекте нам предстояло еще реализовать библиотеку для работы с матрицами, используя объектно-ориентированный подход. Объектно-ориентированный подход позволяет реализовать библиотеку для работы с матрицами в виде отдельного класса, над объектами которого определены операции, представленные как методами, так и стандартными операторами +, -, * и т.д.


## Реализация библиотеки s21_matrix_oop.h

- Программа разработана на языке C++ стандарта C++17 с использованием компилятора gcc
- Матрицы реализованы в виде класса `S21Matrix`
- Доступ к приватным полям `rows` и `columns` реализован через accessor и mutator
- Решение оформлено как статическая библиотека
- Реализованы операции и перегруженые операторы, описанные [здесь](#операции-над-матрицами)
- Подготовлено покрытие unit-тестами функций библиотеки c помощью библиотеки GTest
- Предусмотрен Makefile для сборки библиотеки и тестов (с целями all, clean, test, s21_matrix_oop.a)

## Информация

### Пример класса матрицы на C++

```cpp
class S21Matrix {
    private:
        // Attributes
        int _rows, _cols;         // Rows and columns
        double **_matrix;         // Pointer to the memory where the matrix is allocated

    public:
        S21Matrix();              // Default constructor
        ~S21Matrix();             // Destructor

        void sum_matrix(const S21Matrix& other); 
        // Other methods..
}
```

### Напоминание основных положений о матрице

Матрица - прямоугольная таблица чисел, расположенных в m стрках и n столбцах

```
    1 2 3
A = 4 5 6
    7 8 9
```

```
     1  2  3  4
В =  5  6  7  8
     9 10 11 12
```

Получить нужный элемент можно при помощи индексов, так
A[1,1] = 1, где первый индекс - номер строки, второй - номер столбца.

Порядок матрицы — это число ее строк или столбцов. \
Главная диагональ квадратной матрицы — это диагональ, идущая из левого верхнего в правый нижний угол. \
Прямоугольная матрица (В) — это матрица, у которой число строк не равно числу столбцов. \
Квадратная матрица (А) — это матрица у которой число строк равно числу столбцов.

### Операции над матрицами

Ниже представлено краткое описание операций над матрицами, которые реализованы в разрабатываемой библиотеке.

| Операция    | Описание   | Исключительные ситуации |
| ----------- | ----------- | ----------- |
| `bool eq_matrix(const S21Matrix& other)` | Проверяет матрицы на равенство между собой |  |
| `void sum_matrix(const S21Matrix& other)` | Прибавляет вторую матрицы к текущей | различная размерность матриц |
| `void sub_matrix(const S21Matrix& other)` | Вычитает из текущей матрицы другую | различная размерность матриц |
| `void mul_number(const double num)` | Умножает текущую матрицу на число |  |
| `void mul_matrix(const S21Matrix& other)` | Умножает текущую матрицу на вторую | число столбцов первой матрицы не равно числу строк второй матрицы |
| `S21Matrix transpose()` | Создает новую транспонированную матрицу из текущей и возвращает ее |  |
| `S21Matrix calc_complements()` | Вычисляет матрицу алгебраических дополнений текущей матрицы и возвращает ее | матрица не является квадратной |
| `double determinant()` | Вычисляет и возвращает определитель текущей матрицы | матрица не является квадратной |
| `S21Matrix inverse_matrix()` | Вычисляет и возвращает обратную матрицу | определитель матрицы равен 0 |

Помимо реализации данных операций, также реализованы конструкторы и деструкторы:

| Метод    | Описание   |
| ----------- | ----------- |
| `S21Matrix()` | Базовый конструктор, инициализирующий матрицу некоторой заранее заданной размерностью |  
| `S21Matrix(int rows, int cols)` | Параметризированный конструктор с количеством строк и столбцов | 
| `S21Matrix(const S21Matrix& other)` | Конструктор копирования |
| `S21Matrix(S21Matrix&& other)` | Конструктор переноса |
| `~S21Matrix()` | Деструктор |

А также перегружены следующие операторы, частично соответствующие операциям выше:

| Оператор    | Описание   | Исключительные ситуации |
| ----------- | ----------- | ----------- |
| `+`      | Сложение двух матриц  | различная размерность матриц |
| `-`   | Вычитание одной матрицы из другой | различная размерность матриц |
| `*`  | Умножение матриц и умножение матрицы на число | число столбцов первой матрицы не равно числу строк второй матрицы |
| `==`  | Проверка на равенство матриц (`eq_matrix`) | |
| `=`  | Присвоение матрице значений другой матрицы | |
| `+=`  | Присвоение сложения (`sum_matrix`)   | различная размерность матриц |
| `-=`  | Присвоение разности (`sub_matrix`) | различная размерность матриц |
| `*=`  | Присвоение умножения (`mul_matrix`/`mul_number`) | число столбцов первой матрицы не равно числу строк второй матрицы |
| `(int i, int j)`  | Индексация по элементам матрицы (строка, колонка) | индекс за пределами матрицы |



