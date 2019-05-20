#include "Matrices.h"
#include <iostream>

Matrices::Matrices(std::vector<std::vector<float>>) {};

std::vector<std::vector<float>> matrixTranspose(std::vector<std::vector<float>> matrix){
    int column = matrix.size();
    int row = matrix[0].size();
    std::vector<std::vector<float>> result(column, std::vector<float>(row));
    for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            result[i][j]=matrix[j][i];
        }
    }
    return result;
}

std::vector <std::vector<float>> matrixInverse(std::vector<std::vector<float>> matrix){
    int column = matrix.size();
    int row = matrix[0].size();
    std::vector<std::vector<float>> result(column, std::vector<float>(row));
    //compute inverse
    return result;

}

std::vector<std::vector<float>> matrixMultiplication(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2) {
    int column1 = matrix1.size();
    int row1 = matrix1[0].size();
    int column2 = matrix2.size();
    int row2 = matrix2[0].size();

    std::vector<std::vector<float>> result(column1, std::vector<float>(row2));

    if (row1 != column2) {
        std::cout << "cannot multiply matrices";
    } else {
        for (int i = 0; i < row1; i++) {
            for (int j = 0; j < column2; j++) {
                for (int k = 0; k < row2; k++) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }
    return result;
}

std::vector<std::vector<float>> matrixMultiplicationScalar(std::vector<std::vector<float>> matrix, float scalar) {
    int column = matrix.size();
    int row = matrix[0].size();
    std::vector<std::vector<float>> result(row, std::vector<float>(column));
    for(int i=0;i<row;i++){
        for(int j=0;j<column;j++){
            result[i][j]=matrix[i][j]*scalar;
        }
    }
    return result;
}

std::vector<float> matrixMultiplicationVector(std::vector<std::vector<float>> matrix, std::vector<float> vec) {
    int column = matrix.size();
    int row = matrix[0].size();
    int rowvec = vec.size();

    std::vector<float> result;
    result.resize(rowvec);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < rowvec; j++) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

std::vector<float> vectorMultiplicationScalar(std::vector<float> vec, float scalar){
    std::vector<float> result;
    result.resize(vec.size());

    for (int i = 0; i<vec.size(); i++){
        result[i]=vec[i]*scalar;
    }
    return result;
}

std::vector<std::vector<float>> matrixSubtraction(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2) {
    int column1 = matrix1.size();
    int row1 = matrix1[0].size();
    int column2 = matrix2.size();
    int row2 = matrix2[0].size();
    std::vector<std::vector<float>> result(row1, std::vector<float>(column1));

    for(int i=0; i<row1; i++) {
        for(int j=0; j<column1; j++) {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
    return result;
}

std::vector<float> vectorSubtraction(std::vector<float> vec1, std::vector<float> vec2) {
    int ii, size = vec1.size();
    for (ii=0; ii<size; ii++) {
        vec1[ii] = vec1[ii] - vec2[ii];
    }
    return vec1;
}