#pragma once

#include <vector>

class Matrices {
public:
    Matrices(std::vector<std::vector<float>>);
    std::vector<std::vector<float>> matrixTranspose(std::vector<std::vector<float>> matrix);
    std::vector<std::vector<float>> matrixInverse(std::vector<std::vector<float>> matrix);
    std::vector<std::vector<float>> matrixMultiplication(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2);
    std::vector<std::vector<float>> matrixMultiplicationScalar(std::vector<std::vector<float>> matrix, float scalar);
    std::vector<float> matrixMultiplicationVector(std::vector<std::vecotor<float>> matrix, std::vector<float> vec);
    std::vector<float> vectorMultiplicationScalar(std::vector<float> vec, float scalar);
    std::vector<vector<float>> matrixSubtraction(std::vector<vector<float>> matrix1, std::vector<vector<float>> matri2);
};