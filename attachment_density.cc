// to compile : g++ -O2 -o AD attachment_density.cc /usr/lib64/liblapack.so.3
// to compile icpc attachment_density2.cc -L/opt/intel/composerxe/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
//g++ -O2 -o AD attachment_density.cc -L/usr/lib/lapack/liblapack.a -lblas -llapack

//compare the CI coefficients of two states
/* insert libraries defining standard functions */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include <iostream>
#include <fstream>
//#include "mkl.h"

using namespace std;

extern "C" void dsyev_(char&, char&, int&, double*, int&, double*, double*, int&, int&);   

int main(int argc, char *argv[]){
  
  if (argc != 5){
    cerr << "npi, nci, norb, exstate (0 to read from cin) " << endl;
    exit(1); 
  }
  
  int npi=atoi(argv[1]);
  int nci=atoi( argv[2] );
  int norb=atoi( argv[3] );
  int index=atoi( argv[4]);
  
  float CI[npi+1][npi+1];
  for(int i = 0 ; i <= npi ; i++){
    for(int j = 0 ; j <= npi ; j++){
      CI[i][j]=0.0;
    }
  }
  
  
  double Amath[norb*norb];
  double w[norb];
  double work[3*norb];
  int lwork = 3*norb;
  int info;
  int norb2=npi-norb;
  double w2[norb2];
  double Amate[norb2*norb2];
  double work2[3*norb2];
  int lwork2 = 3*norb2;
    
  ifstream inci("CI_coefs.dat");
  if(!inci){
    cerr << "ERROR: could not find file CI_coefs.dat" << endl;
    exit(1);
  }
  
  float ftemp;
  int itemp[5];
  int exstate;
  int iter;
  int nind=norb*norb;
  int ind[nind][2];
  int ind1,ind2;
  iter=0;
  for(int i = 1 ; i <= norb ; i++){
    for(int j = 1 ; j <= norb ; j++){
      ind[iter][0]=i;
      ind[iter][1]=j;
      iter++;
    }
  }
  
  int step=0;
  while(inci){
    exstate=index;
    if(index==0) cin >> exstate;
    
    for(int i = 0 ; i < npi ; i++){
      for(int j = 0 ; j < npi ; j++){
	CI[i][j]=0.0;
      }
    }
    
    for(int n = 1 ; n <= nci ; n++){
      for(int m = 0 ; m < nci ; m++){
	inci >> itemp[0] >> itemp[1] >> itemp[2] >> itemp[3] >> itemp[4]  >> ftemp;
	//cout << ftemp << endl;
	if(n==exstate){
	  CI[itemp[3]][itemp[4]]=ftemp;
	}
      }
    }
    if(inci){
      iter=0;
      
#pragma omp parallel for private(ftemp,ind1,ind2)
      for(int i = 0 ; i < nind ; i++){
	ind1=ind[i][0];
	ind2=ind[i][1];
	ftemp=0.0;
	for(int n = (norb+1) ; n <= npi ; n++){
	  ftemp+=CI[ind1][n]*CI[ind2][n];
	}
	Amath[i]=ftemp;
	//cout << ftemp << endl;
      }

      char v='V';
      char u='U'; 
      dsyev_(v,u,norb,Amath,norb,w,work,lwork,info); 

      for(int i = 1 ; i <= norb ; i++){
	cout << i << "  " << w[norb-i] << endl;
	for(int j = 0 ; j < norb ; j++){
	  cout << j+1 << "  " << Amath[(norb-i)*norb+j] << endl;
	}
      }
      
      iter=0;
      for(int i = (norb+1) ; i <= npi ; i++){
	for(int j = (norb+1) ; j <= npi ; j++){
	  ftemp=0.0;
	  for(int n = 1 ; n <= norb ; n++){
	    ftemp+=CI[n][i]*CI[n][j];
	  }
	  Amate[iter]=ftemp;
	  iter++;
	  
	  //cout << step << "  " << 1 << "  " << i-norb << "  " << j-norb << "  " << ftemp << endl;
	}
      }

      dsyev_(v,u,norb2,Amate,norb2,w,work,lwork,info);
      
      for(int i = 1 ; i <= norb2 ; i++ ){
	cout << i << "  " << w[norb2-i] << endl;
	for(int j = 0 ; j < norb2 ; j++){
	  cout << j+1 << "  " << Amate[(norb2-i)*norb2+j] << endl;
	}
      }
    }
    step++;
  }
  inci.close();
}
