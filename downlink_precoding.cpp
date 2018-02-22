#include <itpp/itcomm.h>
#include <iostream>
#include <complex>
#include <cstdlib>

using namespace itpp;
using std::cout;
using std::endl;

namespace itpp {
    template<class Num_T> Num_T* begin(Vec<Num_T> &v) {
        return &v[0];
    }

    template<class Num_T> Num_T* end(Vec<Num_T> &v) {
        return &v[v.size()];
    }
};

int main(int argc, char *argv[]){
int iter=10;	  // No. of iterations
int C=3;		  // # of cells in a network
int U=2;	      // # of active users in each cell (worst case)
int N=4;          // # of maximum transmitting antennas at each base station
int P=9; 	      // total length of past samples of y/h_bcu used
it_file ff;
Real_Timer tt; 
tt.tic();
if (argc > 1) {
        iter = atoi(argv[1]);
        cout << "No. of iterations = " << iter << "\n";
}
if (argc > 2) {
        C = atoi(argv[2]);
        cout << "C = " << C << "\n";
}
if (argc > 3) {
        U = atoi(argv[3]);
        cout << "U = " << U << "\n";
}
if (argc > 4) {
        N = atoi(argv[4]);
        cout << "N = " << N << "\n";
}
if (argc > 5) {
        P = atoi(argv[5]);
        cout << "P = " << P << "\n";
}
int l=3;        // # of past samples considered
vec final_sum_rate="0.0"; // stores the total sum
vec avg_sum_rate="0.0";   // average of data rates for U users
ivec Nt_vals = linspace_fixed_step(20, N, 10);
int Nt;	
int no_iterations_in_antenna;
//no_iterations_in_antenna=((N-Nt_vals[0])/10)+1;
//vec Data_Rate[no_iterations_in_antenna];
int count=0;

QPSK qpsk;                     //The QPSK modulator class   
AWGN_Channel awgn_channel;     //The AWGN channel class
/* Generation of H_bc/h_bcu  which is coming from cell to base station */
RNG_randomize();
double Ec=1.0;
double Eb=Ec/2.0;
vec EbN0dB,N0,EbN0;

EbN0dB=linspace(0.0,20,10);
EbN0=inv_dB(EbN0dB);
N0=Eb/EbN0;   //This will also work

cout<<"No. of Cells Selected ="<<C<<endl;
cout<<"No. of Users Selected ="<<U<<endl;
cout<<"No. of antennas at each base station ="<<N<<endl;


















return 0;
}
