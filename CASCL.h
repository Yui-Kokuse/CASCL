// Parity-Check aided Pruning
// Encoding : x = wG ( G is Sparse Matrix(type:COO) )
// Decoding SCL (with CRC)

// author : Yusuke Oki

#include <iostream>
#include <cstdint>
#include <vector>
#include <stack>
#include <cmath>
#include "COO.h"
using namespace std;


class CA_SCL
{

public:
  // Channel Setting
  int _m; // log(N)
  int _list_size;
  int _list_max;
  int _block_length; // power(2, _m)
  int _info_length;
  int _effective_info_length;  

  // CRC setting
  int _crc_length;
  vector<int> _crc_gen;
  vector<int> crc_encode(vector<int> info);
  int crc_decode(vector<int> y);

  // Polar Code
  double _design_ep; 
  COO _Gen; // Generate Mattix
  vector<int> _frozen; // if i is frozen bit's index, _frozen[i] = 0
  
  void init();
  vector<int> encode(vector<int> info_bits);
  vector<int> decode_scl_llr(vector<double> llr);

private:
  std::vector<std::vector<double *>> _arrayPointer_LLR;
  std::vector<double> _pathMetric_LLR;

  stack<int> _inactivePathIndices;
  vector<int> _activePath;
  vector<vector<double *>> _arrayPointer_P;
  vector<vector<int *>> _arrayPointer_C;
  vector<int *> _arrayPointer_Info;
  vector<vector<int>> _pathIndexToArrayIndex;
  vector<stack<int>> _inactiveArrayIndices;
  vector<vector<int>> _arrayReferenceCount;

  vector<int> _channel_order;
  vector<int> _bit_rev_order;

  // function for SCL Decoding
  void initializeDataStructures();
  int assignInitialPath();
  int clonePath(int l);
  void killPath(int l);
  double * getArrayPointer_LLR(int lam, int l);
  int * getArrayPointer_C(int lam, int l);
  void recursivelyCalcLLR(int lam, int phi);
  void recursivelyUpdateC(int lam, int phi);
  void continuePaths_FrozenBit(int phi);
  void continuePaths_UnfrozenBit(int phi);
  int findMostProbablePath();

  void create_bit_rev_order(); 
  void initialize_frozen_bits();  // function for init frozen bit

  // function for Encoding
  void initGenerateMatrix();
};
