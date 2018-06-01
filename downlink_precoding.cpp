#include <itpp/itcomm.h>
#include <iostream>
#include <complex>
#include <cstdlib>
#include<fstream>

using namespace itpp;
using namespace std;
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
	int iter=20;	  // No. of iterations
	int C=3;		  // # of cells in a network
	int U=3;	      // # of active users in each cell (worst case)
	int N=60;          // # of maximum transmitting antennas at each base station
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
	ivec Nt_vals = linspace_fixed_step(20, N, 10);                  //some change
	int Nt;	
	int no_iterations_in_antenna;
	no_iterations_in_antenna=((N-Nt_vals[0])/10)+1;
	vec Data_Rate[no_iterations_in_antenna];
	int count=0;
	//Nt=N; // look here loop terminating guy

	QPSK qpsk;                     //The QPSK modulator class   
	AWGN_Channel awgn_channel;     //The AWGN channel class

	/* Generation of H_bc/h_bcu  which is coming from cell to base station */
	RNG_randomize();
	double Ec=1.0;
	double Eb=Ec/2.0;
	vec EbN0dB,N0,EbN0;

	EbN0dB=linspace(0.0,20,10);
	EbN0=inv_dB(EbN0dB);
	N0=Eb/EbN0;   

	cout<<"No. of Cells Selected ="<<C<<endl;
	cout<<"No. of Users Selected ="<<U<<endl;
	
	//for(int Nt=20;Nt<=N;Nt=Nt+10){
	for (int Nt : Nt_vals) {
         final_sum_rate="0.0";
         cout<<"No. of antennas at each base station ="<<Nt<<endl;
         for (int it = 0; it < iter; ++it) {
			cout<<"The iteration number is ="<<it<<endl;     

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
			vec fDTs= 0.005 + (0.05-0.005)*randu(U);
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
	
	//cout<<"path loss"<<pl<<endl; // pl = 0.0466729
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
		//a_bb[i]=pow(alpha_bcu[0][0][i], v);
		a_bb[i]=zeros(l);
        a_bb[i](0)=alpha_bcu[0][0][i];
		//cout<<"a vec"<<a_bb[i]<<endl;
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

	cvec hsbu_n_n[U];   // hsbu_n_n means h stacked with n given n
	cvec hsbu_n_n_1[U]; // hsbu_n_n_1 means h stacked with n given n-1
	cmat psbu_n_n[U];
	cmat psbu_n_n_1[U];
	cmat p_bu[U];      // Storing error matrix in kalman equation
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
            hsu_n_n=hsu_n_n_1+k_gain_p*(y_tilde[i][p]-((c_kron*hsu_n_n_1).get_col(0))); // c_kron*hsu_n_n_1 is a cvec of form cmat
            psu_n_n=( eye(Nt*l)-(k_gain_p*c_kron) )*psu_n_n_1;
			
            // populate arrays
            K_gain[p][i]=k_gain_p;
            h_hat[p]=hsu_n_n;
            p_bu[p]=psu_n_n;
		}
	}
	
	cmat Fb;  // Precoding matrix for b-th cell
	Fb.set_size(Nt,U);
	Fb=zeros_c(Nt,U);
	cvec h_hat_vec[U];
	cmat P_bu[U]; // P_bu is extraction from stacked p_bu
	/*Extracting p_bu from stacked version of p_bu */
	for(int i=0;i<U;i++){									
		//cout<<"h_hat["<<i<<"] = "<<h_hat[i]<<endl;
		h_hat_vec[i]=h_hat[i];
		h_hat_vec[i]=c_kron*h_hat_vec[i];
		Fb.set_col(i,h_hat_vec[i]);
		//cout<<"h_hat_vec "<<h_hat_vec<<endl;
		//cout<<"Fb["<<i<<"] = "<<Fb<<endl;
		P_bu[i]=c_kron*p_bu[i];
		P_bu[i]=P_bu[i]*transpose(c_kron);
		//cout<<"P_bu "<<i<<"] = "<<P_bu[i]<<endl;
	}
	
	/* Finding lambda for precoding matrix */
	double lambda = 0*U;  // lagrangian lambda
	cvec f_bu[U];
	cmat H_pre[U];
	double lambda_b; // Normalizig factor to have avg transmit power constraint at b=0 base station
	double alpha_b=1.0;
	double pf=1.0;
	cvec y_downlink[U];
	cvec interference[C-1];
	cvec noise[U];
	vec sig_pow[U];
	vec int_pow;
	vec noise_pow[U];
	vec sinr[U];
	vec rate[U];
	vec sum_rate="0.0";
	cmat H_fbu[U]; // H_fbu= H_pre + transpose(H_pre)
	mat H_fbu_r[U];   // H_fbu_r = real(H_fbu)
	mat U_fbu[U],V_fbu[U];
	vec s_fbu[U];
	cmat UhhtransU[U];  // UhhtransU[U]= transpose(U)*h_hat_vec* hermitian(h_hat_vec)*U
	vec g[U];           // g[i] = real(diag(U'*h'h^H*U))
	for(int i=0;i<U;i++){
		H_pre[i]=outer_product(h_hat_vec[i],conj(h_hat_vec[i])) + P_bu[i];
		H_fbu[i]=H_pre[i] + transpose(H_pre[i]);
		H_fbu_r[i]=real(H_fbu[i]);	
		svd(H_fbu_r[i],U_fbu[i],s_fbu[i],V_fbu[i]);
		UhhtransU[i]=U_fbu[i]*outer_product(h_hat_vec[i],conj(h_hat_vec[i]))*U_fbu[i];
		g[i]=real(diag(UhhtransU[i]));
		//cout<<"side wala matrix= "<<g[i]<<endl;
		//cout<<"lamda wala="<<s_fbu[i]<<endl;
		// precoded vector 
		//f_bu[i]=2*inv(2*lambda*eye(Nt) + H_pre[i] + transpose(H_pre[i])) * h_hat_vec[i]; // alpha n b multiplication pending
	}
	
	/* Constructing the polynomial equation */
	// rel= (no. of terms to be convolved)(sequence length) - (no. of terms to be convolved-1)
	// root equation length(rel) or number of terms in the lagrange equation
	int maxcoeffeq = (Nt*U)*(3) - ((Nt*U)-1);  // length when all the terms are multiplied together
	vec eqcoeff[U][Nt]; 			// coeff of lagrangian polynomial equation 
	vec denpoly[U][Nt];
	vec rooteqcoeff=zeros(maxcoeffeq);
	rooteqcoeff(0)=1;

	for(int i=0;i<U;i++){
		for(int j=0;j<Nt;j++){	
			vec denpolyvec=zeros(3);
			if( abs(g[i](j)) < 1e-4){
				g[i](j)=0;
			}
			denpolyvec(0)=4;
			denpolyvec(1)=(4*g[i](j));	 //  second coefficient of the polynomial in denomiator
			denpolyvec(2)=(g[i](j)*g[i](j));	 // third coefficient of the polynomial in denomiator
			denpoly[i][j]=denpolyvec;
			cout<<"denpooly="<<denpoly[i][j]<<endl;
		}
	}
	
	for(int i=0;i<U;i++){
		for(int j=0;j<Nt;j++){
			eqcoeff[i][j]=zeros(maxcoeffeq);
			eqcoeff[i][j](0)=1;
			//cout<<"look i and j="<<i<<"  "<<j<<endl;
			for(int k=0;k<U;k++){
				for(int l=0;l<Nt;l++){
					if(!(i==k && j==l)){
						//cout<<"look at the sequence="<<k<<"  "<<l<<endl;
						eqcoeff[i][j]=filter(denpoly[k][l], "1", eqcoeff[i][j]);
					}
				}
			}
			eqcoeff[i][j]=s_fbu[i](j)*eqcoeff[i][j];
		}
	}
    
    for(int i=0;i<U;i++){
		for(int j=0;j<Nt;j++){	
			rooteqcoeff=filter(denpoly[i][j],"1",rooteqcoeff);
		}
	}
    
    for(int i=0;i<U;i++){
		for(int j=0;j<Nt;j++){	
			rooteqcoeff=rooteqcoeff - eqcoeff[i][j];
		}
	}
	
	double lambdaroots_real;
	cvec roots_of_lambda=roots(rooteqcoeff);
	//cout<<"roots "<<roots_of_lambda<<endl;	
	for(int i=0;i<(size(roots_of_lambda));i++){
		if(real(roots_of_lambda(i))>0 & imag(roots_of_lambda(i))==0){
			lambdaroots_real=real(roots_of_lambda(i));
		}
	}
	 
	cout<<"real roots="<<lambdaroots_real<<endl;
	

	lambda=lambdaroots_real;
	for(int i=0;i<U;i++){
		// precoded vector 
		f_bu[i]=2*inv(2*lambda*eye(Nt) +  H_pre[i] + transpose(H_pre[i])) * h_hat_vec[i] ; // alpha n b multiplication pending
	}
	
	
	
	
    /* Saving vectors in data file*/
   /*ofstream myfile ("pilot_contamination_data.txt");
	if(myfile.is_open()){
		myfile<<"N0. of antenas="<<Nt<<endl;
		for(int i=0;i<U;i++){
			myfile<<"H_pre["<<i<<"]=\n"<<H_pre[i]<<"\n";
			cout<<"H_pre["<<i<<"]=\n"<<H_pre[i]<<"\n";
		}
		for(int i=0;i<U;i++){
			myfile<<"H_pretranspose["<<i<<"]=\n"<<transpose(H_pre[i])<<"\n";
		}
		for(int i=0;i<U;i++){
			myfile<<"h_hat["<<i<<"]=\n"<<h_hat_vec[i]<<"\n";
		}
		myfile<<"End \n";
		myfile.close();
	}
	*/	 

	/* Saving vectors in csv file */
	/*ofstream myfile("H_pre_data.csv");
	if(myfile.is_open()){
		//myfile<<"N0. of antenas="<<Nt<<endl;
		for(int i=0;i<U;i++){
			cout<<"H_pre["<<i<<"]=\n"<<H_pre[i]<<"\n";
			for(int j=0;j<Nt;j++){
				for(int k=0;k<Nt;k++){
					myfile<<H_pre[i](j,k)<<",";
					
				}
				myfile<<"\n";
			}
		}
		for(int i=0;i<U;i++){
			myfile<<"H_pretranspose["<<i<<"]=\n"<<transpose(H_pre[i])<<"\n";
		}
		for(int i=0;i<U;i++){
			myfile<<"h_hat["<<i<<"]=\n"<<h_hat_vec[i]<<"\n";
		}
		myfile.close();
	}
	ofstream hfile("h_hat_data.csv");
	if(hfile.is_open()){
		for(int i=0;i<U;i++){
			cout<<"h_hat["<<i<<"]=\n"<<h_hat_vec[i]<<"\n";
			for(int j=0;j<Nt;j++){
				hfile<<h_hat_vec[i](j)<<",";
			}	
			hfile<<"\n";
		}
	hfile.close();
	}*/
    
    /* f_bu should be multiplied by normalizer*/

    for(int i=0;i<U;i++){
		//cout<<"f_bu["<<i<<"] = "<<f_bu[i]<<endl;
		Fb.set_col(i,f_bu[i]);
		//cout<<"Fb["<<i<<"] = "<<Fb<<endl;
	}
	
	cmat dumFb=Fb*hermitian_transpose(Fb);
	std::complex< double > trace_lambda;
	trace_lambda=trace(dumFb);
	lambda_b=1.0; // see here
	//lambda_b=1.0/real(trace_lambda);
	//cout<<"lambda_b "<<lambda_b<<endl;	
	
	// Multiplying with alpha_b*b=(sqrt(p_f * lambda_b)) to normalize f_bu
	Fb= sqrt(pf * lambda_b)*Fb; // didn't normalize in calculating f_bu 
	lambda_b=1.0/real(trace_lambda);//see here
	/* Calculating signal component in the downlink signal */
	for(int i=0;i<U;i++){
		y_downlink[i]="0+0i";
		//cout<<"h_pbcu[P-1][0][0][i] = "<<h_pbcu[P-1][0][0][i]<<endl;
		y_downlink[i]=sqrt(pf*lambda_b)*conj(h_pbcu[P-1][0][0][i])*Fb.get_col(i)*transmitted_symbols(i); 
	}
	
	/* Interference power */
	cmat Fc[C-1];
	vec lambda_c(C-1);
	cvec dummy_h_interfere;
	cmat dummy_Fc;
	dummy_Fc.set_size(Nt,U);

	for(int i=1;i<C;i++){
		for(int j=0;j<U;j++){
			dummy_Fc.set_col(j,h_pbcu[P-1][i][i][j]);		
		}
		Fc[i-1] =pl*dummy_Fc;
		//cout<<"Fc["<<i<<"-1] "<<Fc[i-1]<<endl;
		cmat dumFc=dummy_Fc*hermitian_transpose(dummy_Fc);
		std::complex< double > trace_lambda;
		trace_lambda=trace(dumFc);
		lambda_c(i-1)=1.0/real(trace_lambda);
		//cout<<"lambda_c "<<i<<"="<<lambda_c[i-1]<<endl;	
	}

	for(int i=0;i<C-1;i++){
		interference[i]	= "0+0i";	
		for(int j=0;j<U;j++){
			interference[i]	=interference[i]+sqrt(pf*lambda_c(i))*conj(Fc[i].get_col(j))*Fc[i].get_col(j)*transmitted_symbols(0); //picking first transmitted symbol
			//cout<<"inner interference[i] "<<interference[i]<<endl;	
		}	
	//cout<<"outer"<<endl;
	}
	
	int_pow="0.0";
	for(int i=0;i<C-1;i++){
		int_pow=int_pow+sqr(abs(interference[i]));
		//cout<<"interer "<<int_pow<<endl;
	}
	vec interference_pow;
	interference_pow=int_pow;
	//cout<<"interference_pow"<<interference_pow<<endl;
	/* Noise Power */
	for(int i=0;i<U;i++){
		noise[i]=randn_c(1);
	}
	int noise_power=1;
	int interference_power=1;
	
	/* Data Rate */
	for(int i=0;i<U;i++){
		sig_pow[i]=sqr(abs(y_downlink[i]));
		//cout<<"sig_power["<<i<<"] = "<<sig_pow[i]<<endl;
		noise_pow[i]=noise_power*sqr(abs(noise[i]));
		//cout<<"noise_power["<<i<<"] = "<<noise_pow[i]<<endl;
		//cout<<"int"<<interference_pow<<endl;
		sinr[i]=elem_div(sig_pow[i],(interference_power*interference_pow+noise_pow[i]));
		//cout<<"sinr["<<i<<"] "<<sinr[i]<<endl;
		rate[i]=log2(1+sinr[i]);
		//cout<<"rate["<<i<<"] "<<rate[i]<<endl;	
	}

	for(int i=0;i<U;i++){
		sum_rate=sum_rate+rate[i];
	}
	
	final_sum_rate+=sum_rate;
	//cout<<"Total sum rate = "<<sum_rate<<endl;
	}
	avg_sum_rate=final_sum_rate/iter;

                cout<<"average of sum rates for "<<Nt<<" antennas "<<avg_sum_rate<<endl;


                Data_Rate[count]=avg_sum_rate;
                ++count;
        }
        
        vec dr(no_iterations_in_antenna);
        for(int i=0;i<no_iterations_in_antenna;i++){
                dr(i)=Data_Rate[i][0];
                cout<<Data_Rate[i]<<"  "<<endl;
		}
		tt.toc();
		tt.toc_print();

return 0;
}
