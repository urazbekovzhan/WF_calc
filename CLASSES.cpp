//
//  classes.cpp
//  WF_calc
//
//  Created by Bakytzhan on 31.01.2018.
//  Copyright © 2018 Bakytzhan. All rights reserved.
//

#include "CLASSES.hpp"
#include "GLOBALS.hpp"
#include "FUNCTIONS.hpp"
#include "CLASSES.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_pow_int.h>
#include <math.h>

MatterDistribution::MatterDistribution(){
    
}

MatterDistribution::~MatterDistribution(){
    
}


double MatterDistribution::SubNucleonDensity(double x){
    double dens0=0;
    if (OwnShoulder==true) {
        for (int part=0; part<I.Partitions; part++) {
            I.SwitchPartitionTo(part);
            for (int i=0; i<I.Dimension(); i++) {
                for (int j=0; j<I.Dimension(); j++) {
                    dens0+=I.c(i)*I.c(j)
                    *MatterDistSubNucleon(x, I.a(i)+I.a(j), (I.b(i)+I.b(j)) /I.y0 /I.y0, I.lam(), I.l());
                }
            }
        }
        
    } else {
        int j2;
        for (int part=0; part<I.Partitions; part++) {
            I.SwitchPartitionTo(part);
            for (int i=0; i<I.Dimension(); i++) {
                for (int j=0; j<I.Dimension(); j++) {
                    for (int j1=0; j1<=I.l()+I.lam(); j1++) {
                        j2=I.l()+I.lam()-j1;
                        dens0+= pow(A(I, j1, j2, i, j),2)
                        *I.c(i)*I.c(j)
                        *MatterDistSubNucleon(x, I.w(i)+I.w(j), 1.5*(I.nu(i)+I.nu(j)) /I.y0 /I.y0, I.lam(),I.l());
                        
                    }
                }
            }
        }
    }
    
    
    
    return dens0;
}


void MatterDistribution::SwitchShoulder(){
    MatterDistribution::OwnShoulder=false;
}


double MatterDistribution::SubAlphaDensity(double x){
    double dens0=0;
    if(MatterDistribution::OwnShoulder){
        
        for(int part=0; part<I.Partitions; part++){
            I.SwitchPartitionTo(part);
            for (int i=0; i<I.Dimension(); i++){
                for (int j=0; j<I.Dimension(); j++){
                    dens0+= I.c(i)*I.c(j)*
                    MATTER_DIST_SUBALPHA(x, (I.a(i)+I.a(j)), (I.b(i)+I.b(j)) /I.y0 /I.y0, I.lam(),I.l());
                }
            }
        }
    
    
    }
    else{
        int j2;
        for (int part=0; part<I.Partitions; part++) {
            I.SwitchPartitionTo(part);
            for (int j1=0; j1<I.lam()+I.l(); j1++) {
                for (int i=0; i<I.Dimension(); i++) {
                    for (int j=0; j<I.Dimension(); j++) {
                        j2=I.lam()+I.l()-j1;
                        dens0+= pow(A(I, j1, j2, i, j),2)
                        *I.c(i)*I.c(j)
                        *MATTER_DIST_SUBALPHA(x, (I.w(i)+I.w(j)), (I.nu(i)+I.nu(j))  /I.y0 /I.y0, I.lam(),I.l());
                        
                    }
                }
            }
        }
        
        
    }
    
    return dens0;
    }






    WaveFunction::WaveFunction(){
        
    }

    WaveFunction::~WaveFunction() {
    }


void WaveFunction::Initialize(char *path){
    FILE *info=fopen(path, "rb");
    int k1,k2,k3,k4;
    double lf1,lf2,lf3,lf4;
    int i=0;
    int j=0;
    char line[1000];
    while (fgets(line, sizeof line, info)) {
        if (sscanf(line, "QUANTUM MOMENTS %2d%2d%2d%2d", &k1, &k2, &k3, &k4) == 4) {
            
            WaveFunction::QuantumMoments[j][0]=k1;
            WaveFunction::QuantumMoments[j][1]=k2;
            WaveFunction::QuantumMoments[j][2]=k3;
            WaveFunction::QuantumMoments[j][3]=k4;
            j++;
            Partitions=j;
            i=0;
        }
        else {
            sscanf(line, "%lf %lf %lf %lf", &lf1, &lf2, &lf3, &lf4);
            WaveFunction::a_data[i][j-1]=lf2;
            WaveFunction::b_data[i][j-1]=lf3;
            WaveFunction::c_data[i][j-1]=lf4;
            i++;
            WaveFunction::dimention_data[j-1]=i;
        }
        
        
    }
    
    
}

    double WaveFunction::lam(){
        return QuantumMoments[CurrentPartition][0];
    }
    
    double WaveFunction::l(){
        return QuantumMoments[CurrentPartition][1];
    }
    
    double WaveFunction::L(){
        return QuantumMoments[CurrentPartition][2];
    }
    
    double WaveFunction::Dimension(){
        return dimention_data[CurrentPartition];
    }
    
    void WaveFunction::SwitchPartitionTo(int new_partition){
        if (new_partition >= Partitions) {
            abort();
        }
        else {
        WaveFunction::CurrentPartition=new_partition;
        }
    }
    
    double WaveFunction::a(int i){
        return a_data[i][CurrentPartition];
    }
    
    double WaveFunction::b(int i){
       return b_data[i][CurrentPartition];
    }
    
    double WaveFunction::c(int i){
       return c_data[i][CurrentPartition];
    }

    double WaveFunction::nu(int i){
        return nu_data[i][CurrentPartition];
    }
    
    double WaveFunction::w(int i){
        return w_data[i][CurrentPartition];
    }

double WaveFunction::rho(int i){
    return rho_data[i][CurrentPartition];
}


double WaveFunction::mu(int i){
    return mu_data[i][CurrentPartition];
}
    
    void WaveFunction::InitNewSet(){
        
        T[0][0]= 1.0;
        T[0][1]= 0.0;
        T[1][0]= 0.0;
        T[1][1]= 1.0;
       
        

/*
 T[0][0]= -0.2;   // He, Li
 T[0][1]= -0.979796;
 T[1][0]= 0.979796;
 T[1][1]= 0.2;       
 
 T[0][0]= -0.8;  // Be
 T[0][1]= -0.6;
 T[1][0]= 0.6;
 T[1][1]= 0.8;
 
 */
        




/*

 T[0][0]= -m1 /(m3+m1);
 T[0][1]= -1.0;
 T[1][0]= m3 *(m1+m2+m3) /(m3+m1) /(m2+m3);
 T[1][1]= -m2 /(m2+m3);
 
 
 T[0][0]= -m1 /(m3+m1);
 T[0][1]= -1.0;
 T[1][0]= m3 *(m1+m2+m3) /(m3+m1) /(m2+m3);
 T[1][1]= -m2 /(m2+m3);
 
 T[0][0]= -m2 /(m2+m3);
 T[0][1]= 1.0;
 T[1][0]= -m3 *(m1+m2+m3) /(m3+m2) /(m1+m3);
 T[1][1]= -m1 /(m1+m3);
 
 */
        

        
        
        for(int part=0; part<I.Partitions; part++){
            SwitchPartitionTo(part);
            for (int i=0; i<Dimension(); i++) {
                mu_data[i][part]=T[0][0]*T[0][0]*a(i)+T[1][0]*T[1][0]*b(i);
                nu_data[i][part]=T[0][1]*T[0][1]*a(i)+T[1][1]*T[1][1]*b(i);
                rho_data[i][part]=2*T[0][0]*T[0][1]*a(i)+2*T[1][1]*T[1][0]*b(i);
                w_data[i][part]=mu_data[i][part]
                -0.25*rho_data[i][part]*rho_data[i][part] /nu_data[i][part];
            }
        }
    }



