#include <iostream>
#include <vector>
#include <random>
#include <chrono> // measuring time
#include <iomanip>  // std::setprecision
#include "CASCL.h"
#include <omp.h> // parallel

using namespace std;
#define rep(i,a,b) for(int i=a;i<b;i++)
#define REP(i,n) rep(i,0,n)


vector<double> err_simulate(CA_SCL scl, double SNR, int Nsim)
{
  // rand settings
  std::mt19937 mt(15);
  std::uniform_real_distribution<double> rand01(0.0, 1.0); // uniform distribution [0, 1]

  int N_now = Nsim;

  int k = scl._effective_info_length;
  int n = scl._block_length;
  double Rate = (double)k / (double)n;

  double ber_avg, fer_avg;
  double err_count = 0.0;
  double fer_count = 0.0;
  vector<double> res(2);
  

  REP(i,Nsim){
    vector<int> info_bits(k, 0);
    vector<int> info_bits_padded(n);
    vector<int> trans_word(n, 0); // transmit word
    vector<double> received_signal(n, 0); // received word
    vector<int> estimated_info_bits;

    // make information bit 
    REP(j, k){
      if(rand01(mt) < 0.0 ){ info_bits[j] = 1; }
    }

    // encode
    trans_word =  scl.encode(info_bits);

    // plus noise : AWGN
    double variance = pow(10.0f, - SNR/10) / (2*Rate);
    std::normal_distribution<> gauss_noise(0.0, sqrt(variance) );
    REP(j, n){
      if(trans_word[j] == 0){ received_signal[j] = 1.0 + gauss_noise(mt); }
      else{ received_signal[j] = - 1.0 + gauss_noise(mt); }
    }
    // init LLR
    vector<double> llr(n);
    REP(j, n){ llr[j] =  2.0  * received_signal[j] / variance; }
    
    // Decode
    estimated_info_bits = scl.decode_scl_llr(llr);

    // Error Check & Calculate bit-error-rate
    REP(j,k){
        if(estimated_info_bits[j] != info_bits[j]){ err_count += 1.0; }
    }

    // FER check
    REP(j,k){
      if(estimated_info_bits[j] != info_bits[j]){ 
        fer_count += 1.0;
        break;  
      }
    }

    if(fer_count > 100){
      N_now = i+1;
      break;
    }
  } // simulate loop end


  // calculate average code error rate
  ber_avg = err_count / (double)N_now / (double)k; 
  fer_avg = fer_count / (double)N_now; 
  res[0] = fer_avg; res[1] = ber_avg;  

  return res;
}

int main()
{
  CA_SCL plscl;
  plscl.init();

  double SNR = 1.0;
  double Nsim = 4000000;
  vector<double> res(2);

  cout << "# Nsim =" << Nsim << endl;
  cout << "# block length = " << plscl._block_length << endl;
  cout << "# info length (without CRC)= " << plscl._effective_info_length << endl;
  cout << "# info length (with CRC)= " << plscl._info_length << endl;
  cout << "# list size = " << plscl._list_max << endl;
  cout << endl;

  chrono::system_clock::time_point  start, end;

  start = std::chrono::system_clock::now(); // clock start
  while(SNR < 3.1){
    res = err_simulate(plscl, SNR, Nsim);
    cout  << scientific << setprecision(2) << SNR  << " ";
    cout << scientific << setprecision(3)  << res[0] << "   " << res[1] << endl; // 
    SNR = SNR + 0.25;
  }
  end = std::chrono::system_clock::now(); // clock end

  cout << endl;
  double elapsed = chrono::duration_cast<std::chrono::seconds>(end-start).count(); //taling time [sec]
  cout << "# it takes " << elapsed << " [sec]" << endl;

  return 0;
}
