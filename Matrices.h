#pragma once

#include <vector>

class Matrices {
public:
    Matrices(std::vector<std::vector<float>>);
    std::vector<std::vector<float>> matrixTranspose(std::vector<std::vector<float>> matrix);
    std::vector<std::vector<float>> matrixInverse(std::vector<std::vector<float>> matrix);
    std::vector<std::vector<float>> matrixMultiplication(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2);
    std::vector<std::vector<float>> matrixMultiplicationScalar(std::vector<std::vector<float>> matrix, float scalar);
};