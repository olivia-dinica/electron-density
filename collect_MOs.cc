/* insert libraries defining standard functions */
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<iostream>
#include<fstream>
using namespace std;

int main(int argc, char *argv[]){
  //this program needs three files 
  //MO_coefs.dat
  if (argc != 5){
    cerr << "Usage: npi,norb,step,index (0 data, 1 movie)" << endl;
    exit(1); 
  }
  
  int npi=atoi(argv[1]);
  int norb=atoi(argv[2]);
  int step=atoi(argv[3]);
  int index=atoi(argv[4]);
  
  float pos[npi][3];
  int itemp[6];
  float ftemp[3];
  float pmag[npi][2];
  
  ifstream inmo("MO_coefs.dat");
  if(!inmo){
    cerr << "***ERROR: could not find the file MO_coefs.dat ***" << endl;
    exit(1);
  }
  
  double holeval[norb];
  int nvac=npi-norb;
  double elecval[nvac];
  
  double holevec[norb][norb];
  double elecvec[nvac][nvac];
  
  float mocoef[npi][npi];
  
  float dr[2];
  float dr2[2];
  double holemag[npi][2];
  double elecmag[npi][2];

  for(int s = 0 ; s <= step ; s++){
    for(int i = 0 ; i < norb ; i++){
      cin >> itemp[0] >> holeval[i];    
      //cout << holeval[i] << "   s: " << s << "   iorb  " << i << endl;
      for(int j = 0 ; j < norb ; j++){
	cin >> itemp[0] >> holevec[i][j];
      }
    }

    for(int i = 0 ; i < nvac ; i++){
      cin >> itemp[0] >> elecval[i];
      for(int j = 0 ; j < nvac ; j++){
	cin >> itemp[0] >> elecvec[i][j];
      }
    }
  
    //read in MO coeficients
    for(int n = 0 ; n < npi ; n++){
      inmo >> itemp[0] >> itemp[1] >> pos[n][0] >> pos[n][1] >> pos[n][2];
	//cout << pos[n][0] << endl;
      for(int i = 0 ; i < npi ; i++){
	inmo >> mocoef[i][n];
      }
    }
    
    for(int n = 0 ; n < npi ; n++){
      holemag[n][0]=0.0;
      holemag[n][1]=0.0;
      elecmag[n][0]=0.0;
      elecmag[n][1]=0.0;
      for(int i = 0 ; i < norb ; i++){
	for(int j = 0 ; j < norb ; j++){
	  holemag[n][0]+=holeval[i]*holevec[i][j]*mocoef[j][n]; //gives an idea about the phase of the wavefunction. was spit out back when adam was trying to track diabatic states along a trajectory
	  holemag[n][1]+=holeval[i]*holevec[i][j]*holevec[i][j]*mocoef[j][n]*mocoef[j][n];//hole density 
	}
      }
      for(int i = 0 ; i < nvac ; i++){
	for(int j = 0 ; j < nvac ; j++){
	  elecmag[n][0]+=elecval[i]*elecvec[i][j]*mocoef[j+norb][n];
	  elecmag[n][1]+=elecval[i]*elecvec[i][j]*elecvec[i][j]*mocoef[j+norb][n]*mocoef[j+norb][n];
	}
      }
    }
    for(int i = 0 ; i < npi ; i++){
      cout << s << "  " << i << "  " << holemag[i][0] << "  " << elecmag[i][0] << "  " << holemag[i][1] << "  " << elecmag[i][1] << "  " << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << endl;
    }
  }
  inmo.close();
}
