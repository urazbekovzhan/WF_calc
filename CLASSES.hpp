//
//  classes.hpp
//  WF_calc
//
//  Created by Bakytzhan on 31.01.2018.
//  Copyright © 2018 Bakytzhan. All rights reserved.
//

#ifndef CLASSES_hpp
#define CLASSES_hpp
#include <gsl/gsl_matrix.h>
#include <stdio.h>


class MatterDistribution {
public:
    MatterDistribution();
    ~MatterDistribution();
    void SwitchShoulder();
    double SubAlphaDensity(double x);
    double SubNucleonDensity(double x);
private:
    bool OwnShoulder=true;
    
    
};




class WaveFunction {
public:
    WaveFunction();
    ~WaveFunction();
    void Initialize(char *path);
    double lam();
    double l();
    double L();
    double Dimension();
    void SwitchPartitionTo(int NewPartition);
    double a(int i);
    double b(int i);
    double c(int i);
    double nu(int i);
    double w(int i);
    double rho(int i);
    double mu(int i);
    void InitNewSet();
    int Partitions;
    double T[4][4]={};
    double y0;
    double m1;
    double m2;
    double m3;
private:
    static const int N = 49;
    static const int n = 7;
    int CurrentPartition=0;
    
    
    
    double a_data[100][10]={};
    double b_data[100][10]={};
    double c_data[100][10]={};
    double QuantumMoments[10][4]{};
    
    double dimention_data[10]={};
    
    double nu_data[100][10]={};
    double w_data[100][10]={};
    double mu_data[100][10]={};
    double rho_data[100][10]={};
};
#endif /* classes_hpp */
