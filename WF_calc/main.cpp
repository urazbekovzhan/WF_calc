#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "FUNCTIONS.hpp"
#include "GLOBALS.hpp"
#include "CLASSES.hpp"
using namespace std;


int main (void)
{
   
    I.Initialize("/Users/zhan/Desktop/MacDesktop/9be_exe/be2.wf");
    I.m1=1;
    I.m2=4;
    I.m3=4;
    I.y0= (I.m2+I.m3) /(I.m1+I.m2+I.m3);
    
    I.InitNewSet();
    MatterDistribution Nucleus;
    cout<<"SubAlpha Matter Distr calculation in progress..."<<endl;
    FILE *MDSNucl=fopen("/Users/zhan/Desktop/MacDesktop/9be_exe/be_subnucl.plot", "wb");
    for(double R=0.05; R<10.5; R=R+0.05){
        
        fprintf(MDSNucl,"%.3F\t%.5E\n",R,Nucleus.SubNucleonDensity(R));
        fprintf(stdout,"%.3F\t%.5E\n",R,Nucleus.SubNucleonDensity(R));
    }
    fclose(MDSNucl);
    
    
    
    
    I.y0= (I.m1+I.m3) /(I.m1+I.m2+I.m3);
    cout<<"SubNucleon Matter Distr calculation in progress..."<<endl;
    FILE *MDS=fopen("/Users/zhan/Desktop/MacDesktop/9be_exe/be_subalpha.plot", "wb");
    Nucleus.SwitchShoulder();
    for(double R=0.05; R<10.5; R=R+0.05){
        
        fprintf(MDS,"%.3F\t%.5E\n",R,Nucleus.SubAlphaDensity(R));
        fprintf(stdout,"%.3F\t%.5E\n",R,Nucleus.SubAlphaDensity(R));
    }
    fclose(MDS);

    
    
    double sum=0;
    for(int part1=0; part1<I.Partitions; part1++){
        I.SwitchPartitionTo(part1);
        for (int i=0; i<I.Dimension(); i++){
            for (int j=0; j<I.Dimension(); j++){
                sum+= INTEGRAL_FERMDIST(2*I.lam(),I.a(i)+I.a(j)) *INTEGRAL_FERMDIST(2*I.l(),I.b(i)+I.b(j))
                *I.c(j)  *I.c(i);
            }
        }
        
    }
    cout<<sum<<endl;
    sum=0;
    I.InitNewSet();
    int j2;
    double wij;
    for(int part1=0; part1<I.Partitions; part1++){
        I.SwitchPartitionTo(part1);
        for (int i=0; i<I.Dimension(); i++){
            for (int j=0; j<I.Dimension(); j++){
                for (int j1=0; j1<=I.l()+I.lam(); j1++) {
                    j2=I.l()+I.lam()-j1;
                    wij=I.mu(i)+I.mu(j)-0.25*pow(I.rho(i)+I.rho(j),2)/(I.nu(i)+I.nu(j));
                    sum+= INTEGRAL_FERMDIST(2*j1,wij) *INTEGRAL_FERMDIST(2*j2,(I.nu(i)+I.nu(j)))
                    *I.c(j)  *I.c(i)
                    * A(I, j1, j2, i, j)*A(I, j1, j2, i ,j);
                }
            }
        }
        
    }
    

    cout<<sum<<endl;
    
    
    
    FILE *out=fopen("/Users/zhan/Desktop/MacDesktop/9be_exe/3d.data", "wb");
    for (double x=0.5; x<10.5; x=x+0.5) {
        for (double y=0.5; y<10.5; y=y+0.5) {
            fprintf(out,"%.3F\t%.3F\t%.10E\n",x,y,COERRELATION_FUNC(x, y));
        }
    }
    
    /*
     
     int limit;
     double sum=0;
     for(int part1=0; part1<7; part1++){
     I.switch_partition_to(part1);
     limit=I.dimension();
     for (int i=0; i<limit; i++){
     for (int j=0; j<limit; j++){
     sum+= INTEGRAL_FERMDIST(I.lam(),I.a(i)+I.a(j)) *INTEGRAL_FERMDIST(I.l(),I.b(i)+I.b(j))
     *I.c(j)  *I.c(i);
     }
     }
     
     }
     
     
     
     int j2;
     double wij;
     for(int part1=0; part1<7; part1++){
     I.switch_partition_to(part1);
     limit=I.dimension();
     for (int i=0; i<limit; i++){
     for (int j=0; j<limit; j++){
     for (int j1=0; j1<=I.l()+I.lam(); j1++) {
     j2=I.l()+I.lam()-j1;
     wij=I.mu(i)+I.mu(j)-0.25*pow(I.rho(i)+I.rho(j),2)/(I.nu(i)+I.nu(j));
     sum+= INTEGRAL_FERMDIST(j1,wij) *INTEGRAL_FERMDIST(j2,(I.nu(i)+I.nu(j)))
     *I.c(j)  *I.c(i)
     * A(I, j1, j2, i, j)*A(I, j1, j2, i ,j);
     }
     }
     }
     
     }
     
     FILE *out=fopen("/Users/zhan/WF_9Be/plot.data", "wb");
     for (double x=0.2; x<7.2; x=x+0.2) {
     for (double y=0.2; y<8.2; y=y+0.2) {
     fprintf(out,"%.3F\t%.3F\t%.5E\n",x,y,COERRELATION_FUNC_OF2DSET(I, x, y));
     }
     }
     
     */
 
    
    return 0;
    
}
