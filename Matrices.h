#pragma once

#include <vector>

class Matrices {
public:
    Matrices(vector <vector<float>>);
    vector <vector<float>> matrixTranspose(vector <vector<float>> matrix);
    vector <vector<float>> matrixInverse(vector <vector<float>> matrix);
}