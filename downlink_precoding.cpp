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
	int P=9; 	      // total length of past samples of y/h_bcu used (Samples for operating)
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
	Nt=N; // look here loop terminating guy

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

	/* QPSK Symbol generation */
	int Number_of_bits;
	bvec transmitted_bits;
	cvec transmitted_symbols;
	Number_of_bits=2*U;  // Focusing on one cell only & others will follow
	transmitted_bits = randb(Number_of_bits);
	transmitted_symbols = qpsk.modulate_bits(transmitted_bits);

    // Corr matrix of the channel coefficients
    cmat R_bcu[C][C][U];
    cmat Rsqr_bcu[C][C][U];  // Rsqr_bcu represents the square root of the corresponding complex matrix

    for(int i=0;i<C;i++) {
    	for(int j=0;j<C;j++){
    	    for(int k=0;k<U;k++){
    	        cmat unsymmetricR_bcu = randn_c(Nt,Nt);  //dummy cmat to construct correlation mat
                R_bcu[i][j][k]= unsymmetricR_bcu*hermitian_transpose(unsymmetricR_bcu);
                Rsqr_bcu[i][j][k]= sqrtm(R_bcu[i][j][k]); 
			}
		}
	}
	
    // Declaring channel coeff.
    cvec v_bcu[C][C][U];
    cvec h_bcu[C][C][U];
    for(int i=0;i<C;i++) {
    	for(int j=0;j<C;j++) {
    		for(int k=0;k<U;k++){
    			v_bcu[i][j][k]= randn_c(Nt);
                h_bcu[i][j][k]= Rsqr_bcu[i][j][k] * v_bcu[i][j][k] ; 
                //cout<<"hbcu"<<h_bcu[i][j][k]<<endl;
    	    }
        }
	}
	/* alpha will tell about the dopler spread */
	vec alpha_bcu[C][C];  // alpha with range 0.005 to 0.05 uniformly
    for(int i=0;i<C;i++) {
		for(int j=0;j<C;j++){
			vec fDTs= 0.1 + (0.3-0.1)*randu(U);
		    alpha_bcu[i][j]=besselj(0,2*3.14*fDTs);
		    //cout<<alpha_bcu[i][j][k]<<endl;
		}
	}
	cout<<"Channel coefficients declared"<<endl;
	
	/* Creating P samples of h_bcu using AR(1) evolution */
    cvec h_pbcu[P][C][C][U];
    for(int i=0;i<C;i++) { 									// cell iterator
	    for(int j=0;j<C;j++){								// base station iterator
	        for(int k=0;k<U;k++){							// user number iterator
	            h_pbcu[0][i][j][k]=h_bcu[i][j][k];
	            //cout<<h_pbcu[0][i][j][k]<<endl;
                for(int p=1;p<P;p++){ 						// time instant iterator
    	            h_pbcu[p][i][j][k]=(alpha_bcu[i][j][k]*h_pbcu[p-1][i][j][k])+(sqrtm(1-sqr(alpha_bcu[i][j][k])*R_bcu[i][j][k]))*randn_c(Nt);
                    //cout<<"h_pbcu"<<h_pbcu[p][i][j][k]<<endl;
                }
			}
		}
	}

	/* Combining all channel coefficients of all users for a given cell to a given base station */
    cmat H_pbc[P][C][C];
	for(int p=0;p<P;p++){
    	for(int i=0;i<C;i++){
    	    for(int j=0;j<C;j++){
    	        H_pbc[p][i][j]=zeros_c(Nt,U);
                for(int k=0; k<U;k++){
    	            H_pbc[p][i][j].set_col(k,h_pbcu[p][i][j][k]);
                    //cout<<H_pbc[p][i][j]<<endl;
                }
			}
		}
	}

	/* To find the sum rate of all users in the center cell-> b=0 & c=0  */
    double pl_dB; // path loss
    double d=1.35;  // in Km
    pl_dB=12.81+ 3.76* log10(d); // PL = 128.1 + 37.6 log_{10}(d)
    double pl;
    pl=1/(inv_dB(pl_dB));
	cout<<"path loss"<<pl<<endl; // pl = 0.0466729
    /* Multiply pl with H_pbc[0][c] c=1:C-1 */
    for(int p=0;p<P;p++){
	    for(int i=1;i<C;i++){
	        H_pbc[p][0][i]=pl*H_pbc[p][0][i];
	        //cout<<H_pbc[p][0][i]<<endl;
        }
    }
    
    cvec sum_hpbcu[P][U];
	for(int p=0;p<P;p++){
		for(int k=0;k<U;k++){
			sum_hpbcu[p][k]=zeros_c(Nt);
            //sum_hbcu[k].clear();
            for(int i=0;i<C;i++){
	            sum_hpbcu[p][k]=sum_hpbcu[p][k]+H_pbc[p][0][i].get_col(k);
                //cout<<H_pbc[p][0][i].get_col(k)<<endl;
                //cout<<sum_hpbcu[p][k]<<endl;
			}
		}
	}

    // y at pth instant of all users in cell 0
    cvec y_tilde[P][U];
    cvec n[P][U];
    for(int i=0;i<U;i++){
        for(int p=0;p<P;p++){
            n[p][i]=sqrt(N0(0)/2)*randn_c(Nt); //noise
            y_tilde[p][i]=sum_hpbcu[p][i]+ n[p][i];
            //cout<<"y_tilde["<<p<<"]["<<i<<"]"<<y_tilde[p][i]<<endl;
        }
    }

	vec a_bb[U];
    vec v=linspace(1,l,l);
    //cout<<v<<endl;

    for(int i=0;i<U;i++){
		//cout<<"alpha_bcu a value"<<alpha_bcu[0][0][i]<<endl;
		a_bb[i]=pow(alpha_bcu[0][0][i], v);
		cout<<"a vec"<<a_bb[i]<<endl;
    }

	//cout<<a_bb[0]<<endl;
	mat A_u[U];
	for(int k=0;k<U;k++){
		    A_u[k].set_size(l,l);
		    A_u[k]=zeros(l,l);
		    A_u[k].set_row(0,a_bb[k]);
		    //cout<<dfg[k]<<endl;
		    for(int i=0;i<l-1;i++){
		            //vec g=eye(l).get_row(i+1);
		            //A_u[k].set_row(i,g);
		            A_u[k].set_row(i+1,eye(l).get_row(i));
		            //cout<<A_u[k]<<endl;
		    }
	}

	cvec hsbu_n_n[U];
	cvec hsbu_n_n_1[U];
	cmat psbu_n_n[U];
	cmat psbu_n_n_1[U];
	
	/* Initialisation of kalman filter */
	for(int i=0;i<U;i++){
	    for(int j=l-1;j>=0;j--){
            hsbu_n_n[i]=concat(hsbu_n_n[i],y_tilde[j][i]);
	    }
	}
	for(int i=0;i<U;i++){
	    psbu_n_n[i]=outer_product(hsbu_n_n[i],conj(hsbu_n_n[i]));
	}
	mat Qv=eye(Nt*l); //Qv is yet to define properly right now unit variance
	cmat K_gain[U][P-l];
	mat I=eye(l);
	mat c=I.get_row(0);
	mat c_kron_col=kron(c,eye(Nt));
	mat c_kron=transpose(c_kron_col);
	cvec h_hat[U];
	int Debug1=0;
	
	/* Implementing kalman filter */
	for(int p=0;p<U;p++){						// for a selected user

	    cvec hsu_n_n_1;
	    cmat psu_n_n_1;
	    cvec hsu_n_n=hsbu_n_n[p];
	    cmat psu_n_n=psbu_n_n[p];
	    //cout<<psu_n_n;
	    cmat k_gain_p;
	    cmat Q_w=0*randn_c(Nt*l,Nt*l);    // Q_ set to zero matrix
	    //hsu_n_n_1=c_kron*hsu_n_n;
	    int Debug1=0;
	    for(int i=0;i<P-l;i++){					// for the selected instance
		    if (Debug1){
		    	cout << " y_tilde index = "<<i<<","<<p<< endl;
		    }
            hsu_n_n_1=kron(A_u[p],eye(Nt))*hsu_n_n;
            psu_n_n_1=kron(A_u[p],eye(Nt))*psu_n_n*hermitian_transpose(kron(A_u[p],eye(Nt)))+Q_w;
            k_gain_p=psu_n_n_1* hermitian_transpose(c_kron)* inv(c_kron*psu_n_n_1*hermitian_transpose(c_kron)+ c_kron*Qv*hermitian_transpose(c_kron));
            hsu_n_n=hsu_n_n_1+k_gain_p*(y_tilde[i][p]-((c_kron*hsu_n_n_1).get_col(0))); // U did this bcoz c_kron*hsu_n_n_1 is a cvec of form cmat
            psu_n_n=( eye(Nt*l)-(k_gain_p*c_kron) )*psu_n_n_1;
			
            // populate arrays
            K_gain[p][i]=k_gain_p;
            h_hat[p]=hsu_n_n;
		}
	}

















return 0;
}
