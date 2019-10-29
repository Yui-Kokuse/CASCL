#include <iostream>
#include <fstream>
#include <cmath>  // log
#include <vector> // vector
#include <algorithm> // sort
#include "CASCL.h"
#define filename "GenerateMatrix_2048.txt"


using namespace std;
#define rep(i,a,b) for(int i=a;i<b;i++)
#define REP(i,n) rep(i,0,n)


void CA_SCL::init()
{
  // CRC (CRC-16)
  _crc_length = 16;
  _crc_gen.resize(_crc_length+1);
  _crc_gen = {1,1,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,0, 1 }; 

  // PolarCode
  _m = 11;
  _block_length = 2048;
  _info_length = 1040; // Rate = 1/2
  _effective_info_length = _info_length - _crc_length; 
  _design_ep = 0.32;
  _list_size = 1;
  _list_max = 1; 

  _frozen.resize(_block_length);
  create_bit_rev_order();
  initialize_frozen_bits();
  initGenerateMatrix();
}


// CRC function
vector<int> CA_SCL::crc_encode(vector<int> info)
{
  int K = info.size();
  int N = K + _crc_length;
  vector<int> encoded_word(N);
  vector<int> rem(N); // remainder
  REP(i, K){ rem[i] = info[i]; }

  REP(i, K){
    if(rem[i] == 1){
      int j_gen = 0;
      rep(j_rem,i,N-K+1+i){
        rem[j_rem] = (rem[j_rem] + _crc_gen[j_gen])%2;
        j_gen++;
      }
    }
  }

  REP(i,K){ encoded_word[i] = info[i]; }
  rep(i,K,N){ encoded_word[i] = rem[i]; }

  return encoded_word; 
}


int CA_SCL::crc_decode(vector<int> received_word)
{
  // if CRC confirms received_word true, it returns '1'
  // else, it returns '0' 

  int N = received_word.size();
  int K = N - _crc_length;

  vector<int> rem(N); // remainder
  REP(i, N){ rem[i] = received_word[i]; } 
  int judge = 1;

  // division
  REP(i, K){
    if(rem[i] == 1){
      int j_gen = 0;
      rep(j_rem,i,N-K+1+i){
        rem[j_rem] = (rem[j_rem] + _crc_gen[j_gen])%2;
        j_gen++;
      }
    }
  }

  REP(i, N){
    if(rem[i] == 1){
      judge = 0;
      break;
    }
  }

  return judge;
}

void CA_SCL::initialize_frozen_bits()
{
  vector<double> channel_vec(_block_length);

  // init channel capacity
  REP(i, _block_length) {
    channel_vec.at(i) = _design_ep;
  }

  // calculate polar channel capacity
  REP(k, _m) {
    int increment = 1 << k;
    REP(j, increment){
      for (int i = 0; i < _block_length ; i += 2 * increment) {
        double c1 = channel_vec.at(i + j);
        double c2 = channel_vec.at(i + j + increment);
        channel_vec.at(i + j) = c1 + c2 - c1*c2;
        channel_vec.at(i + j + increment) = c1*c2;
      }
    }
  }

  // argsort : channel order
  _channel_order.resize(_block_length);
  std::size_t n_t(0);
  std::generate(std::begin(_channel_order), std::end(_channel_order), [&]{ return n_t++; });
  std::sort(  std::begin(_channel_order),
                std::end(_channel_order),
                [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] < channel_vec[_bit_rev_order.at(i2)]; } );

  REP(i, _info_length){
    _frozen.at(_channel_order.at(i)) = 1; // unfrozen
  }
  // else frozen

}


void CA_SCL::create_bit_rev_order() {
  _bit_rev_order.resize(_block_length);
  REP(i, _block_length) {
    int to_be_reversed = i;
    _bit_rev_order.at(i) = (uint16_t) ((to_be_reversed & 1) << (_m - 1));
    for (int j = (int)(_m-1);j; --j) {
      to_be_reversed >>= 1;
      _bit_rev_order.at(i) += (to_be_reversed & 1) << (j - 1);
    }
  }
}


void CA_SCL::initGenerateMatrix()
{
  vector<vector<int>> A;
  ifstream ifs;
  ifs.open(filename);

  // file -> A
  A.resize(_block_length);
  REP(i, _block_length){
    A[i].resize(_block_length);
    REP(j, _block_length){
      ifs >> A[i][j];
    }
  }

  _Gen.initSparseMatrix_fromNormal(A);
}


vector<int> CA_SCL::encode(vector<int> info_bits)
{
  /*
    input : info_bits( size = _effective_info_length )
  */
  vector<int> info_bits_crc(_info_length);
  vector<int> info_bits_padded(_block_length);
  vector<int> coded_bits(_block_length);

  // CRC encoding
  info_bits_crc = crc_encode(info_bits);

  // pad information (frozen and unfrozen)
  int iter_info = 0;
  int tmp = 0;
  REP(i, _block_length){
    if(_frozen[i] == 0){
      info_bits_padded[i] = 0; // pad frozen bit
    }
    else{
      info_bits_padded[i] = info_bits_crc[iter_info]; // pad unfrozen bit
      iter_info++;
    }
  }

  // x = uG ( it takes O(NlogN) )
  REP(i, _Gen._value_size){
    coded_bits.at(_Gen._column[i]) ^= _Gen._value[i] * info_bits_padded.at(_Gen._row[i]);
  }

  return coded_bits;
}


// SCL - Decoding

void CA_SCL::initializeDataStructures()
{
  while(_inactivePathIndices.size()){
    _inactivePathIndices.pop();
  };
  _activePath.resize(_list_size);

  _pathMetric_LLR.resize(_list_size);
  _arrayPointer_LLR.resize(_m+1);
  REP(i, _m+1){ _arrayPointer_LLR[i].resize(_list_size); }

  _arrayPointer_C.resize(_m+1);
  REP(i, _m+1){ _arrayPointer_C[i].resize(_list_size); }

  _arrayPointer_Info.resize(_list_size);

  _pathIndexToArrayIndex.resize(_m+1);
  REP(i, _m+1){ _pathIndexToArrayIndex[i].resize(_list_size); }

  _inactiveArrayIndices.resize(_m+1);
  REP(i, _m+1){
    while(_inactiveArrayIndices[i].size()){ _inactiveArrayIndices[i].pop(); };
  }

  _arrayReferenceCount.resize(_m+1);
  REP(i, _m+1){ _arrayReferenceCount[i].resize(_list_size); }

  REP(s, _list_size){
    _arrayPointer_Info[s] = new int[_block_length]();
    REP(lam, _m+1){ 
      _arrayPointer_LLR[lam][s] = new double[(int)(pow(2, _m-lam))]();
      _arrayPointer_C[lam][s] = new int[2 * (int)(pow(2, _m-lam))]();
      _arrayReferenceCount[lam][s] = 0;
      _inactiveArrayIndices[lam].push(s);
    }
  }

  REP(l, _list_size){
    _activePath[l] = 0;
    _inactivePathIndices.push(l);
    _pathMetric_LLR[l] = 0;
  }
}


int CA_SCL::assignInitialPath()
{
  int l = _inactivePathIndices.top();
  _inactivePathIndices.pop();
  _activePath[l] = 1;

  REP(lam, _m+1){
    int s = _inactiveArrayIndices[lam].top();
    _inactiveArrayIndices[lam].pop();
    _pathIndexToArrayIndex[lam][l] = s;
    _arrayReferenceCount[lam][s] = 1;
  }

  return l;
}


int CA_SCL::clonePath(int l)
{
  int lp = _inactivePathIndices.top(); _inactivePathIndices.pop();
  _activePath[lp] = 1;
  
  _pathMetric_LLR[lp] = _pathMetric_LLR[l];

  REP(lam, _m+1){
    int s = _pathIndexToArrayIndex[lam][l];
    _pathIndexToArrayIndex[lam][lp] = s;
    _arrayReferenceCount[lam][s]++;
  }

  return lp;
}


void CA_SCL::killPath(int l)
{
  _activePath[l] = 0;
  _inactivePathIndices.push(l);
  _pathMetric_LLR[l] = 0;

  REP(lam, _m+1){
    int s = _pathIndexToArrayIndex[lam][l];
    _arrayReferenceCount[lam][s]--;
    if(_arrayReferenceCount[lam][s] == 0){
      _inactiveArrayIndices[lam].push(s);
    }
  }
}


double * CA_SCL::getArrayPointer_LLR(int lam, int l)
{
  int s = _pathIndexToArrayIndex[lam][l];
  int sp;

  if(_arrayReferenceCount[lam][s] == 1){ sp = s; }
  else{
    sp = _inactiveArrayIndices[lam].top();  _inactiveArrayIndices[lam].pop();

    // copy
    std::copy(_arrayPointer_C[lam][s], _arrayPointer_C[lam][s] +  (int)(pow(2, _m-lam+1)) ,  _arrayPointer_C[lam][sp]);
    std::copy(_arrayPointer_LLR[lam][s], _arrayPointer_LLR[lam][s] + (int)(pow(2, _m-lam)) ,  _arrayPointer_LLR[lam][sp]);

    _arrayReferenceCount[lam][s]--;
    _arrayReferenceCount[lam][sp] = 1;
    _pathIndexToArrayIndex[lam][l] = sp;
  }

  return _arrayPointer_LLR[lam][sp];
}


int * CA_SCL::getArrayPointer_C(int lam, int l)
{
  int s = _pathIndexToArrayIndex[lam][l];
  int sp;

  if(_arrayReferenceCount[lam][s] == 1){ sp = s; }
  else{
    sp = _inactiveArrayIndices[lam].top();  _inactiveArrayIndices[lam].pop();

    // copy
    std::copy(_arrayPointer_LLR[lam][s], _arrayPointer_LLR[lam][s] +  ( (int)pow(2, _m-lam) ),  _arrayPointer_LLR[lam][sp]);
    std::copy(_arrayPointer_C[lam][s], _arrayPointer_C[lam][s] +  ( (int)pow(2, _m-lam+1)),  _arrayPointer_C[lam][sp]);


    _arrayReferenceCount[lam][s]--;
    _arrayReferenceCount[lam][sp] = 1;
    _pathIndexToArrayIndex[lam][l] = sp;
  }

  return _arrayPointer_C[lam][sp];
}


void CA_SCL::recursivelyCalcLLR(int lam, int phi)
{
  if(lam==0){ return; }
  int psi = phi/2;

  if(phi%2 == 0){ recursivelyCalcLLR(lam-1, psi); }

  REP(l, _list_size){
    if(_activePath[l]==0){ continue; }
    double * llr_lam = getArrayPointer_LLR(lam, l);
    double * llr_lam_1 = getArrayPointer_LLR(lam-1, l);
    int * c_lam = getArrayPointer_C(lam, l);

    REP( beta, (int)( pow(2, _m-lam) ) ){
      if(phi%2 == 0){
        if( 40 > max(abs(llr_lam_1[2*beta]), abs(llr_lam_1[2*beta+1])) ){
          llr_lam[beta] = log( (exp(llr_lam_1[2*beta] + llr_lam_1[2*beta+1]) + 1) /
                               (exp(llr_lam_1[2*beta]) + exp(llr_lam_1[2*beta+1])) );
        }
        else{
          llr_lam[beta] = (double) ( (llr_lam_1[2*beta]<0) ? -1 : (llr_lam_1[2*beta]>0) ) *
                                   ( (llr_lam_1[2*beta+1]<0) ? -1 : (llr_lam_1[2*beta+1]>0) ) *
                                   min( abs(llr_lam_1[2*beta]), abs(llr_lam_1[2*beta+1]) );
        }
      }
      else{
        int up = c_lam[2*beta];
        llr_lam[beta] = (1-2*up) * llr_lam_1[2*beta] + llr_lam_1[2*beta+1];
      }
    }
  }
}


void CA_SCL::recursivelyUpdateC(int lam, int phi)
{
  int psi = phi/2;
  REP(l, _list_size){
    if(_activePath[l] == 0){ continue; }

    int *c_lam = getArrayPointer_C(lam, l);
    int *c_lam_1 = getArrayPointer_C(lam-1, l);

    REP( beta, (int)( pow(2, _m-lam) ) ){
      c_lam_1[2*(2*beta) + (psi%2)] = (int)( ( c_lam[2*beta] + c_lam[2*beta+1] )%2 );
      c_lam_1[2*(2*beta+1) + (psi%2)] = c_lam[2 * beta + 1];
    }
  }

  if(psi%2 == 1){ recursivelyUpdateC( (int)(lam-1), psi); }
}



void CA_SCL::continuePaths_FrozenBit(int phi)
{

  REP(l, _list_size){
    if(_activePath[l] == 0){ continue; }
    int *c_m = getArrayPointer_C(_m, l);
    c_m[phi%2] = 0;
    double *llr_p = getArrayPointer_LLR(_m, l);
    _pathMetric_LLR[l] += log( 1 + exp(-llr_p[0]) );

    _arrayPointer_Info[l][phi] = 0;
  }
}


void CA_SCL::continuePaths_UnfrozenBit(int phi)
{
  vector<double> ProbForks(2*_list_size);
  vector<double> prob;
  vector<int> contForks(2* _list_size);

  int i = 0;
  REP(l, _list_size){
    if(_activePath[l] == 0){
      ProbForks[2*l+1] = NAN;
      ProbForks[2*l]   = NAN;
    }
    else{
      double *llr_p = getArrayPointer_LLR(_m, l);
      ProbForks[2*l]   = -( _pathMetric_LLR[l] + log(1 + exp(-llr_p[0])) );
      ProbForks[2*l+1] = -( _pathMetric_LLR[l] + log(1 + exp( llr_p[0])) );

      prob.push_back( ProbForks[2*l] );
      prob.push_back( ProbForks[2*l+1] );
      i++;
    }


  }

  // i = _list_size - 1
  int rho = _list_size;
  if((2*i) < _list_size){ rho = 2*i; }


  // START - implemention of [2] Algorithm13 line 14 -
  REP(l, 2*_list_size){ contForks[l] = 0; }
  sort( prob.begin(), prob.end(), greater<double>() );

  double threshold = prob[rho-1];
  int num_paths_continued = 0;

  REP(l, 2*_list_size){
    if(ProbForks[l] > threshold){
      contForks[l] = 1;
      num_paths_continued++;
    }
    if(num_paths_continued == rho){ break; }
  }

  if( num_paths_continued < rho ){
    REP(l, 2*_list_size){
      if(ProbForks[l] == threshold){
        contForks[l] = 1;
        num_paths_continued++;
      }
      if(num_paths_continued == rho){ break; }
    }
  }
  // END - implemention of [2] Algorithm13 line 14 -



  // kill-off non-continuing paths
  REP(l, _list_size){
    if(_activePath[l] == 0){ continue; }
    if( contForks[2*l]==0 && contForks[2*l+1] == 0){ killPath(l); }
  }

 

  // continue reveant paths, and duplicate if necessary
  REP(l, _list_size){
    if( contForks[2*l]==0 && contForks[2*l+1]==0 ){ continue; }

    int * c_m = getArrayPointer_C(_m, l);

    if( contForks[2*l]==1 && contForks[2*l+1]==1 ){

      c_m[phi%2] = 0;
      int lp = clonePath(l);
      c_m = getArrayPointer_C(_m, lp);
      c_m[phi%2] = 1;
      copy(_arrayPointer_Info[l], _arrayPointer_Info[l] + phi,  _arrayPointer_Info[lp]);
      _arrayPointer_Info[l][phi] = 0;
      _arrayPointer_Info[lp][phi] = 1;

      double *llr_p = getArrayPointer_LLR(_m, l);
      _pathMetric_LLR[l] += log( 1 + exp(-llr_p[0]) );
      llr_p = getArrayPointer_LLR(_m, lp);
      _pathMetric_LLR[lp] += log( 1 + exp( llr_p[0]) );
    }
    else{
    	if(contForks[2*l] == 1){
    		c_m[phi%2] = 0;
      	_arrayPointer_Info[l][phi] = 0;

      	double *llr_p = getArrayPointer_LLR(_m, l);
      	_pathMetric_LLR[l] += log( 1 + exp(-llr_p[0]) );
    	}
    	else{
    		c_m[phi%2] = 1;
      	_arrayPointer_Info[l][phi] = 1;

      	double *llr_p = getArrayPointer_LLR(_m, l);
      	_pathMetric_LLR[l] += log( 1 + exp( llr_p[0]) );
    	}
    }


  }
}

int CA_SCL::findMostProbablePath()
{
  // CRC-Check
  vector<int> estimated(_info_length);
  REP(l, _list_size){
    int * C0 = _arrayPointer_Info[l];
    // pointer to vector
    int tmp = 0;
    REP(beta, _block_length){
      if(_frozen[beta] != 0){ // unfrozen
        estimated[tmp] = C0[beta];
        tmp += 1;
      } 
    }

    if(crc_decode(estimated) == 1){ return l; }
  }

  // else pathmetric check
  int lp = 0;
  double p_llr = numeric_limits<double>::max();
  REP(l, _list_size){
    if(_activePath[l] == 0){ continue; }

    if(_pathMetric_LLR[l] < p_llr){
      p_llr = _pathMetric_LLR[l];
      lp = l;
    }
  }

  return lp;
}

vector<int> CA_SCL::decode_scl_llr(vector<double> llr)
{

  _list_size = _list_max;
  initializeDataStructures();
  
  int l = assignInitialPath();
  double *llr_0 = getArrayPointer_LLR(0,l);

  REP(beta, _block_length){
    llr_0[beta] = llr[beta]; // init LLR
  }

  REP(phi, _block_length){ // main loop

    recursivelyCalcLLR(_m, phi);

    if(_frozen[phi] == 0){ continuePaths_FrozenBit(phi); }
    else{ continuePaths_UnfrozenBit(phi); }

    if(phi%2 == 1){ recursivelyUpdateC(_m, phi); }

  } // main loop end


  l = findMostProbablePath();
  int * C0 = _arrayPointer_Info[l];
  vector<int> estimated_info_bits(_effective_info_length);

  // pointer to vector
  int tmp = 0;
  REP(beta, _block_length){
    if(_frozen[beta] != 0){ // unfrozen
      estimated_info_bits[tmp] = C0[beta];
      tmp += 1;
      if(tmp == _effective_info_length -1){ break; }   
    } 
  }

  // delete data
  REP(s, _list_size){
    delete[] _arrayPointer_Info[s];
    REP(lam, _m+1){
      delete[] _arrayPointer_LLR[lam][s];
      delete[] _arrayPointer_C[lam][s];
    }
  }

  return estimated_info_bits;
}
