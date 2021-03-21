#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <random>
#include <string>
#include <stdlib.h>
#include "BiCen_Fam.h"

int main(){
    int seed,nclust,nrep; long ngrid=20;
    int px = 4; int p = 9 + px; double rho = 0.15;
    
    string inputfile = "simu_data_int.dat";
    string out_beta_file="beta.dat";
    string out_lambdaT_file="lambdaT.dat";
    string out_lambdaD_file="lambdaD.dat";
    string out_file="res.dat";
    ofstream myfile_beta (out_beta_file);
    ofstream myfile_lambdaT (out_lambdaT_file);
    ofstream myfile_lambdaD (out_lambdaD_file);
    ofstream myfile (out_file);
    if (myfile.is_open())
    {   ifstream data_file(inputfile);
        if (! data_file.is_open()){
            cout << "Unable to open input file" << endl;
            exit(1);
        }
        int nsub=2613; MatrixXd dat(nsub,p);
        // read in the data
        for (int i=0; i<nsub;++i) { for (int j=0; j<p && data_file >> dat(i,j); ++j) {}}
        
        VectorXi cluster(nsub); VectorXi relation(nsub);
        VectorXi Rinf(nsub); VectorXi epsilon(nsub); VectorXi Gene(nsub);
        
        for (int i=0;i<nsub;++i){
            cluster(i) = dat(i,0); relation(i) = dat(i,2);
            Rinf(i) = dat(i,7); epsilon(i) = dat(i,4); Gene(i) = dat(i,8);
        }
        VectorXd L = dat.col(5); VectorXd R = dat.col(6); VectorXd U = dat.col(3);
        MatrixXd X = dat.block(0,9,nsub,px);
        
        Gene_Fam Fam_info(relation, cluster, Gene, rho);
        BiCen_Fam data(L, R, Rinf, U, epsilon, X, Fam_info, ngrid);
        data.solve();
        myfile_beta<<data.beta_.transpose()<<" "<<data.gamma_.transpose()<<" ";
        myfile_beta<<data.rho_.transpose()<<" "<<data.sigma2_.transpose()<<endl;
        
        for (int j=0;j<data.mT_;++j) myfile_lambdaT<<data.tT_(j)<<" "<<data.lambda_(j)<<endl;
        for (int j=0;j<data.mD_;++j) myfile_lambdaD<<data.tD_(j)<<" "<<data.alpha_(j)<<endl;
        
        myfile<<data.iter_<<" "<<data.logLik_MLE_i.sum()<<endl;
        
        data.est_sd_score(5.0); myfile_beta<<data.sd_score.transpose()<<endl;
        
        myfile_beta.close(); myfile_lambdaT.close(); myfile_lambdaD.close(); myfile.close();
        cout << "File written!";
    }
    else cout << "Unable to open file";
    return(0);
}
