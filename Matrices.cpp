#include "Matrices.h"

Matrices::Matrices(vector < float > < vector > ) {};

vector <vector<float>> matrixTranspose(vector <vector<float>> matrix){
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
vector <vector<float>> matrixInverse(vector <vector<float>> matrix){
    std::vector<std::vector<float>> result(row, std::vector<float>(column));
    //compute inverse
    return result;

}
