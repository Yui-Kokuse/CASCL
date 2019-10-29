#include "COO.h"
#include <random>
using namespace std;

#define rep(i,a,b) for(int i=a; i<b; i++)
#define REP(i,n) rep(i,0,n)
#define ROW_SIZE 10
#define COLUMN_SIZE 4


void COO::initSparseMatrix_Random()
{
  _M = ROW_SIZE; _N = COLUMN_SIZE; _value_size = 10;
  _value.resize(_value_size);
  _row.resize(_value_size);
  _column.resize(_value_size);

  // rand
  mt19937 mt(5);
  uniform_int_distribution<> rand01(1, 3);
  uniform_int_distribution<> randM(0, _M-1);
  uniform_int_distribution<> randN(0, _N-1);

  REP(i, _value_size){
    _value[i] = rand01(mt);
    _row[i] = randM(mt);
    _column[i] = randN(mt);
    
    REP(j, i){
      // it takes O(_value_size) time ! too slow !
      // It prevents duplication of (row, column) pair
      if( _row[j]==_row[i] && _column[j]==_column[i] ){ // duprication
        i--;
        break;
      }
    }
    
  }
}

void COO::initSparseMatrix_fromNormal(vector<vector<int>> A)
{
  // require A as _M × _N matrix

  int iter = 0;
  _M = A.size(); _N = A[0].size();
  _value_size = _M * _N; // as max
  _value.resize(_value_size);
  _row.resize(_value_size);
  _column.resize(_value_size);


  int sparse_iter = 0;
  REP(i, _M){
    REP(j, _N){
      if(A[i][j] != 0){
        _value[sparse_iter] = A[i][j];
        _row[sparse_iter] = i;
        _column[sparse_iter] = j;
        sparse_iter++;
      }
    }
  }

  _value_size = sparse_iter;
  _value.resize(_value_size);
  _row.resize(_value_size);
  _column.resize(_value_size);
}


vector<vector<int>> COO::COOtoNormal()
{
  vector<vector<int>> out;
  out.resize(_M); REP(i, _M){ out[i].resize(_N); }

  REP(i, _value_size){
    out.at(_row[i]).at(_column[i]) = _value[i];
  }

  return out;
}
