#include "Matrices.h"
#include <iostream>
Matrices::Matrices(std::vector<std::vector<float>>) {};

std::vector<std::vector<float>> matrixTranspose(std::vector<std::vector<float>> matrix){
    int row = matrix.size();
    int column = matrix[0].size();
    std::vector<std::vector<float>> result(row, std::vector<float>(column));
    for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            result[i][j]=matrix[j][i];
        }
    }
    return result;
}
std::vector <std::vector<float>> matrixInverse(std::vector<std::vector<float>> matrix){
    int row = matrix.size();
    int column = matrix[0].size();
    std::vector<std::vector<float>> result(row, std::vector<float>(column));
    //compute inverse
    return result;

}

std::vector<std::vector<float>> matrixMultiplication(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2) {
    int row1 = matrix1.size();
    int column1 = matrix1[0].size();
    int row2 = matrix2.size();
    int column2 = matrix2[0].size();
    if (column2 == 1) {
        std::vector<float> result;
        result.resize(row1);
    } else{
        std::vector<std::vector<float>> result(row1, std::vector<float>(column2));
    }

    if (column1 != row2) {
        std::cout << "cannot multiply matrices";
    } else {
        for (i = 0; i < row1; ++i) {
            for (j = 0; j < column2; ++j) {
                for (k = 0; k < column1; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }
    return result;

}

std::vector<std::vector<float>> matrixMultiplicationScalar(std::vector<std::vector<float>> matrix, float scalar) {
    int row = matrix.size();
    int column = matrix[0].size();
    std::vector<std::vector<float>> result(row, std::vector<float>(column));
    for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            result[i][j]=matrix[i][j]*scalar;
        }
    }
    return result;
}
