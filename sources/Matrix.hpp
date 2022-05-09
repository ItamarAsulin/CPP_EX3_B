//
// Created by itamarq on 4/18/22.
//

#ifndef MATRIX_CALCULATOR_B_MATRIX_HPP
#define MATRIX_CALCULATOR_B_MATRIX_HPP
#include <iostream>
#include <vector>
using namespace std;
namespace zich {
    class Matrix {
    private:
        int rows;
        int cols;
        double sum;
        vector<vector<double>> mat;

        void free();
    public:
        Matrix();
        Matrix(std::vector<double> scalars, int rows, int cols);

         int  getRows() const;
         int  getCols() const;
        void setRows(int rows);
        void setCols(int cols);
        double  getValueAt(unsigned int row, unsigned int column)const;
        void setValueAt(unsigned int row, unsigned int column, const double value);
        double  getSum()const;

        Matrix operator+() const;
        Matrix& operator+=(Matrix& other);
        friend Matrix operator+(const Matrix& first, const Matrix &second);
        Matrix& operator++();
        Matrix operator++(int);
        Matrix operator-() const;
        Matrix& operator-=(Matrix& other);
        friend Matrix operator-(const Matrix& first, const Matrix &second);
        Matrix& operator--();
        Matrix operator--(int);

        friend bool operator==(const Matrix& first,const  Matrix &second);
        friend bool operator!=(const Matrix& first, const Matrix &second);
        friend bool operator<(const Matrix& first, const Matrix &second);
        friend bool operator<=(const Matrix& first, const Matrix &second);
        friend bool operator>(const Matrix& first, const Matrix &second);
        friend bool operator>=(const Matrix& first, const Matrix &second);

        friend Matrix operator*(const Matrix& first, const Matrix &second);
        void operator*=(double scalar);
        void operator*=(const Matrix &matrix);
        friend Matrix operator*(const Matrix& matrix, double scalar);
        friend Matrix operator*( double scalar, const Matrix& matrix);


        friend std::ostream& operator<< (std::ostream& output, Matrix& matrix);
        friend std::istream& operator>> (std::istream& input, Matrix& matrix);


    };

}


#endif //MATRIX_CALCULATOR_B_MATRIX_HPP
