#include <iostream>
#include <vector>
using namespace std;


class COO{
public:
  int _M, _N;
  int _value_size;
  vector<int> _value;
  vector<int> _row;
  vector<int> _column;
  void initSparseMatrix_Random();
  void initSparseMatrix_fromNormal(vector<vector<int>> A); 
  vector<vector<int>> COOtoNormal();
};
