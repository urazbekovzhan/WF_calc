//
//  functions.cpp
//  WF_calc
//
//  Created by Bakytzhan on 28.01.2018.
//  Copyright © 2018 Bakytzhan. All rights reserved.
//

#include "FUNCTIONS.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_pow_int.h>
#include <math.h>
#include "CLASSES.hpp"
#include "GLOBALS.hpp"



double COERRELATION_FUNC_OF2DSET (double x, double y){
    double funct0=1;
//    double funct1=0;
//    double funct2=0;
    int j2;
//    double wij;
    for (int part=0; part<I.Partitions; part++) {
        I.SwitchPartitionTo(part);
        for (int i=0; i<I.Dimension(); i++) {
            for (int j=0; j<I.Dimension(); j++) {
                for (int j1=0; j1<=I.l()+I.lam(); j1++) {
                    j2=I.l()+I.lam()-j1;
     //               wij=I.mu(i)+I.mu(j)-0.25*pow(I.rho(i)+I.rho(j),2)/(I.nu(i)+I.nu(j));
                    funct0+= pow(x, 2*j1+2)  *pow(y, 2*j2+2)
                    * A(I, j1, j2, i, j)*A(I, j1, j2, i ,j)
                    *I.c(i)*I.c(j)
                    *exp(  -(I.w(i)+I.w(j)) *x*x  -(I.nu(i)+I.nu(j))*y*y   );
                    
                }
            
            }
        }
    }
    return funct0;
}


double COERRELATION_FUNC (double x, double y){
    double func0=0;
    for (int part=0; part<I.Partitions; part++) {
        I.SwitchPartitionTo(part);
        for (int i=0; i<I.Dimension(); i++) {
            for (int j=0; j<I.Dimension(); j++) {
                func0+=I.c(i)*I.c(j)  *pow(x, 2*I.lam()+2)  *pow(y, 2*I.l()+2)
                *exp(  -(I.a(i)+I.a(j))*x*x  -(I.b(i)+I.b(j))*y*y   );
                
            }
        }
    }
    return func0;
}


double INTEGRAL_FERMDIST(int l, double a){          // Integral[ x^(2+l)Exp(-a x^2 )  , x:0->Inf]
    return 0.5*gsl_sf_gamma((3+double(l))/2.0)*pow(a,-(3+double(l))/2.0);
}


double Pij (WaveFunction& I, WaveFunction& J, int i, int j){
    double pij;
    pij = pow(I.a(i),0.5*I.lam()+0.75);
    pij*= pow(J.a(j),0.5*J.lam()+0.75);
    pij*= pow(I.b(i),0.5*I.l()+0.75);
    pij*= pow(J.b(j),0.5*J.l()+0.75);
    return pij;
}

double Aij (WaveFunction& I, WaveFunction& J, double t, int i, int j){
    double aij;
    aij = pow(I.a(i)+J.a(j),0.5*I.lam()+0.5*J.lam()+0.5*t);
    return aij;
}

double Bij (WaveFunction& I, WaveFunction& J, double t, int i, int j){
    double bij;
    bij = pow(I.b(i)+J.b(j),0.5*I.l()+0.5*J.l()+0.5*t);
    return bij;
}

double overlap_ME(WaveFunction& I, WaveFunction& J, int i, int j){
    double overlap;
    overlap=pow(2, I.lam()+I.l()+3);
    overlap*=Pij(I,I, i, j);
    overlap=overlap/Aij(I, I, 3,  i, j);
    overlap=overlap/Bij(I, I, 3,  i, j);
    return overlap;
}

double norm (WaveFunction& I, int i){
    double norm0=1;
    norm0*=2*pow(I.a(i),I.lam()+1.5);
    norm0*=pow(I.b(i),I.l()+1.5);
    norm0*= 1. /3.14159265359 /gsl_sf_doublefact(2* I.l() +1 ) /gsl_sf_doublefact(2 *I.lam() +1);
    norm0=sqrt(norm0);
    norm0*=pow(2.,I.l()+I.lam()+3.);
    return norm0;
}


inline double doubler (double l){
    return 2*l+1;
}

double A (WaveFunction& I, double  j7, double j8, int i, int j){
    
    double a=-0.5*(I.rho(i)+I.rho(j)) /(I.nu(i)+I.nu(j));
    

    double t11=I.T[0][0]+a*I.T[0][1];
    double t12=I.T[0][1];
    double t21=I.T[1][0]+a*I.T[1][1];
    double t22=I.T[1][1];
    
    double j1,j2,j4,j5;
    double j3=I.lam();
    double j6=I.l();
    double j9=I.L();
    
    
    double A0;
    double A1;
    double sum=0;
    for( j1=0; j1<=j3; j1++){
        j2=j3-j1;
        for(  j4=0; j4<=j6; j4++){
            j5=j6-j4;
            A0=1;
            A0*=gsl_sf_coupling_9j(2*j1, 2*j2, 2*j3,
                                   2*j4, 2*j5, 2*j6,
                                   2*j7, 2*j8, 2*j9);
            if (A0==0) continue;
            
            A0*=gsl_sf_coupling_3j(2*j1, 2*j4, 2*j7, 0, 0, 0);
            if (A0==0) continue;
            
            A0*=gsl_sf_coupling_3j(2*j2, 2*j5, 2*j8, 0, 0, 0);
            if (A0==0) continue;
            
            A0*=gsl_sf_pow_int(-1, j3+j6);
            A0*=gsl_sf_pow_int(t11, j1) *gsl_sf_pow_int(t12, j2)
                *gsl_sf_pow_int(t21, j4) *gsl_sf_pow_int(t22, j5);
            
            A1=1.;
            A1*=gsl_sf_fact(doubler(j3)) *gsl_sf_fact(doubler(j6));
            A1*= doubler(j1) *doubler(j2) *doubler(j3) *doubler(j4)
            *doubler(j5) *doubler(j6) *doubler(j7) *doubler(j8);
            
            A1*=1./( gsl_sf_fact(doubler(j1)) *gsl_sf_fact(doubler(j2)) );
            A1*=1./( gsl_sf_fact(doubler(j4)) *gsl_sf_fact(doubler(j5)) );
            A1=sqrt(A1);
            sum+=A0*A1;
        }
        
    }
    return sum;
}

double MATTER_DIST_SUBALPHA(double R, double aij, double bij, int lx, int ly){
    double gamma0=0.7024;
    double rho0=0.4229;
    double matter_dist0= rho0
        *15.74960994572242
        *pow(aij,-1.5-lx)
        *pow(bij+gamma0,-1.5-ly)
        *gsl_sf_gamma(1.5+lx)
        *gsl_sf_doublefact(2*ly)
        *gsl_sf_laguerre_n(ly, 0.5, -R*R*gamma0*gamma0/(bij+gamma0))
        *exp(R*R*(-gamma0+gamma0*gamma0/(bij+gamma0)));
    return matter_dist0;
}

double MatterDistSubNucleon(double y, double aij, double bij, int lx, int ly){
    double dens0;
    dens0=INTEGRAL_FERMDIST(lx, aij)
    *pow(y, 2*ly+2)
    *exp(-bij*y*y);
    return dens0;
}



