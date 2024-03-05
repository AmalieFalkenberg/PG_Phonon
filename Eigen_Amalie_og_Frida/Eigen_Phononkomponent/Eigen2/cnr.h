/* Include file cnr.h - 
* Declarations for allocation routines of vectors and matrices 
* Last change: KL 01.08.02
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 2.4, August 2002
* 
============================================
*/
#pragma once
#include<complex.h>
  #define komplex double complex // Changed by Mads
  #define kvector cvector
  #define freekvector freecvector
  #define kmatrix cmatrix
  #define freekmatrix freecmatrix
  #define zero (komplex) 0.0 // Changed by Mads
  #define one (komplex) 1.0 // Changed by Mads
  #define Arg(a) (a==zero ? 0.0 : arg(a))
  #define real creal // Added by Mads
  #define imag cimag // Added by Mads
  #define arg(a) atan2(imag(a),real(a))
  #define abs(a) sqrt(real(a)*real(a) + imag(a)*imag(a))
  #define sqrabs(a) (real(a)*real(a) + imag(a)*imag(a))
  #define eksp(a) exp(real(a))*(cos(imag(a))+I*sin(imag(a))) 
  #define conj(a)  real(a)-I*imag(a)
  #define skrt(a) (sqrt(abs(a))*(cos(Arg(a)/2.0)+I*sin(Arg(a)/2.0) ) )
//#define skrt(a) sqrt(a)

long long *ivector(long long nl, long long nh);
void freeivector(long long*,long long,long long);
long long **imatrix(long long,long long,long long,long long);
void freeimatrix(long long**,long long,long long,long long,long long);

unsigned long long *lvector(long long nl, long long nh);
void freelvector(unsigned long long*,long long,long long);
unsigned long long **lmatrix(long long,long long,long long,long long);
void freelmatrix(unsigned long long**,long long,long long,long long,long long);

float *vector(long long nl, long long nh);
void freevector(float*,long long,long long);
float **matrix(long long, long long, long long, long long);
void freematrix(float**,long long,long long,long long,long long);

double *dvector(long long nl, long long nh);
void freedvector(double*,long long,long long);
double **dmatrix(long long, long long, long long, long long);
void freedmatrix(double**,long long,long long,long long,long long);

komplex* cvector(long long nl, long long nh);
void freecvector(komplex*,long long,long long);
komplex **cmatrix(long long, long long, long long, long long);
void freecmatrix(komplex**,long long,long long,long long,long long);

void nrerror(const char*);


