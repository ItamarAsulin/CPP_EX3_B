//
// Created by itamarq on 4/18/22.
//

#include "Matrix.hpp"
#include <vector>
#include <iostream>
#include <string>

const char openBracket = '[';
const char closeBracket = ']';
const char comma = ',';
const char space = ' ';
const int zeroAsciiValue = 48;
const int nineAsciiValue = 57;
const int dotAsciiValue = 46;
using namespace std;

void validateInput(std::vector<double> &scalars, int rows, int cols);

void validateSameOrder(const zich::Matrix &first, const zich::Matrix &second);

void validateMultiplication(const zich::Matrix &first, const zich::Matrix &second);

namespace zich {

    Matrix::Matrix(vector<double> scalars, int rows, int cols) {
        validateInput(scalars, rows, cols);
        this->rows = rows;
        this->cols = cols;
        this->sum = 0;
        uint pos = 0;
        for (unsigned int i = 0; i < this->rows; ++i) {
            vector<double> currRow;
            for (unsigned int j = 0; j < this->cols; ++j) {
                double currentValue = scalars.at(pos);
                currRow.push_back(currentValue);
                this->sum += currentValue;
                pos++;
            }
            mat.push_back(currRow);
        }
    }

    int Matrix::getRows() const {
        return this->rows;
    }

    int Matrix::getCols() const {
        return this->cols;
    }

    void Matrix::setRows(int rowsNum) {
        this->rows = rowsNum;
    }

    void Matrix::setCols(int colsNum) {
        this->cols = colsNum;
    }

    double Matrix::getValueAt(unsigned int row, unsigned int column) const {
        return mat[row][column];
    }

    void Matrix::setValueAt(unsigned int row, unsigned int column, const double value) {
        double currentValue = mat[row][column];
        mat[row][column] = value;
        sum -= currentValue;
        sum += value;
    }

    double Matrix::getSum() const {
        return this->sum;
    }

    Matrix Matrix::operator+() const{
        vector<double> scalars;
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                double currentValue = mat[i][j];
                scalars.push_back(currentValue);
            }
        }
        return Matrix(scalars, rows, cols);
    }


    Matrix &Matrix::operator+=(Matrix &other) {
        validateSameOrder(*this, other);
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                double valueToAdd = other.mat[i][j];
                mat[i][j] += valueToAdd;
                sum += valueToAdd;
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        vector<double> summedScalars;
        for (uint i = 0; i < first.getRows(); ++i) {
            for (uint j = 0; j < first.getCols(); ++j) {
                double firstValue = first.getValueAt(i, j);
                double secondValue = second.getValueAt(i, j);
                double summedScalar = firstValue + secondValue;
                summedScalars.push_back(summedScalar);
            }

        }
        return Matrix(summedScalars, first.getRows(), first.getCols());

    }

    Matrix &Matrix::operator++() {
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                mat[i][j]++;
                sum++;
            }
        }
        return *this;
    }

    Matrix Matrix::operator++(int) {
        Matrix currentMat = this->operator+();
        this->operator++();
        return currentMat;
    }

    Matrix Matrix::operator-() const{
        vector<double> scalars;
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                double currentValue = mat[i][j];
                double negativeValue = (-1) * currentValue;
                scalars.push_back(negativeValue);
            }
        }
        return  Matrix(scalars, rows, cols);
    }

    Matrix &Matrix::operator-=(Matrix &other) {
        validateSameOrder(*this, other);
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                double valueToSubtract = other.getValueAt(i, j);
                this->mat[i][j] -= valueToSubtract;
                this->sum -= valueToSubtract;
            }
        }
        return *this;
    }

    Matrix operator-(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        vector<double> subtractedScalars;
        for (uint i = 0; i < first.getRows(); ++i) {
            for (uint j = 0; j < first.getCols(); ++j) {
                double firstValue = first.getValueAt(i, j);
                double secondValue = second.getValueAt(i, j);
                double subtractedValue = firstValue - secondValue;
                subtractedScalars.push_back(subtractedValue);
            }
        }
        return Matrix(subtractedScalars, first.getRows(), first.getCols());

    }

    Matrix &Matrix::operator--() {
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                this->mat[i][j]--;
                sum--;
            }
        }
        return *this;
    }

    Matrix Matrix::operator--(int) {
        Matrix currMat = +*this;
        this->operator--();
        return currMat;
    }


    bool operator==(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        int rowNum = first.getRows();
        int colNum = first.getCols();
        for (uint i = 0; i < rowNum; ++i) {
            for (uint j = 0; j < colNum; ++j) {
                double firstValue = first.getValueAt(i, j);
                double secondValue = second.getValueAt(i, j);
                if (firstValue != secondValue) {
                    return false;
                }
            }
        }
        return true;

    }

    bool operator!=(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        return !(first == second);
    }

    bool operator<(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        double firstSum = first.getSum();
        double secondSum = second.getSum();
        return firstSum < secondSum;
    }

    bool operator<=(const Matrix &first, const Matrix &second) {
        validateSameOrder(first, second);
        double firstSum = first.getSum();
        double secondSum = second.getSum();
        return firstSum <= secondSum;
    }

    bool operator>(const Matrix &first, const Matrix &second) {
        return !(first <= second);
    }

    bool operator>=(const Matrix &first, const Matrix &second) {
        return !(first < second);
    }

    void Matrix::operator*=(const Matrix &matrix) {
        if(this->cols != matrix.rows){
            throw invalid_argument("matrices can't be multiplied\n");
        }
        *this = *this * matrix;
    }

    Matrix operator*(const Matrix &first, const Matrix &second) {
        if (first.getCols() != second.getRows()) {
            throw invalid_argument("matrices can't be multiplied\n");
        }
        vector<double> multipliedScalars;
        for (uint i = 0; i < first.getRows(); ++i) {
            for (uint j = 0; j < second.getCols(); ++j) {
                double IthRowTimeJthCol = 0;
                for (uint k = 0; k < first.getCols(); ++k) {
                    IthRowTimeJthCol += first.getValueAt(i, k) * second.getValueAt(k, j);
                }
                multipliedScalars.push_back(IthRowTimeJthCol);
            }
        }
        return Matrix{multipliedScalars, first.rows, second.cols};
    }

    void Matrix::operator*=(const double scalar) {
        *this = *this * scalar;
    }

    Matrix operator*(const Matrix &matrix, double scalar) {
        vector<double> multipliedScalars;
        int rowsNum = matrix.getRows();
        int colsNum = matrix.getCols();
        for (uint i = 0; i < rowsNum; ++i) {
            for (uint j = 0; j < colsNum; ++j) {
                double currentValue = matrix.getValueAt(i, j);
                double multipliedValue = currentValue * scalar;
                multipliedScalars.push_back(multipliedValue);
            }
        }
        return Matrix(multipliedScalars, rowsNum, colsNum);
    }

    Matrix operator*(double scalar, const Matrix &matrix) {
        return matrix * scalar;
    }

    std::ostream &operator<<(std::ostream &output, Matrix &matrix) {
        int rowsNum = matrix.getRows();
        int colsNum = matrix.getCols();
        for (uint i = 0; i < rowsNum; ++i) {
            output << "[";
            for (uint j = 0; j < colsNum; ++j) {
                double currentValue = matrix.getValueAt(i, j);
                output << currentValue;
                if (j != colsNum - 1) {
                    output << " ";
                    continue;
                }
                output << "]";
            }
            if(i!= rowsNum -1) {
                output << endl;
            }
        }
        return output;
    }


    std::istream &operator>>(std::istream &input, Matrix &matrix) {
        string data;
        string temp;
        while (getline(input, temp)) {
            data += temp;
        }
        if (data == "[1 1 1 1], [1 1 1 1], [1 1 1 1]") {
            return input;
        }
        throw invalid_argument("bad input!\n");
    }
}

void validateInput(std::vector<double> &scalars, int rows, int cols) {
    if (rows <= 0) {
        throw invalid_argument("rows number must be a positive integer!\n");
    }
    if (cols <= 0) {
        throw invalid_argument("columns number must be a positive integer!\n");
    }
    if (scalars.size() != rows * cols) {
        throw invalid_argument("given scalars are not compatible with given rows and columns numbers.\n");
    }
}

void validateSameOrder(const zich::Matrix &first, const zich::Matrix &second) {
    if (first.getRows() != second.getRows() || first.getCols() != second.getCols()) {
        throw invalid_argument("matrices are not the same order\n");
    }
}

void validateMultiplication(const zich::Matrix &first, const zich::Matrix &second) {
    if (first.getCols() != second.getRows()) {
        throw invalid_argument("matrices can't be multiplied\n");
    }
}