#ifndef BiCens_Fam_mis_h
#define BiCens_Fam_mis_h
#include "BiCens_Fam.h"

BiCens_Fam simu_data_mis(long nclust, VectorXd beta, VectorXd gamma, VectorXd rho, VectorXd sigma2, double rho_pop, long ngrid, default_random_engine & generator){
    //Distributions
    bernoulli_distribution X1dist(0.5);
    normal_distribution<double> X2dist(0.0,1.0);
    bernoulli_distribution Bdist(0.5);
    normal_distribution<double> bdist(0.0,sqrt(sigma2(0)));
    normal_distribution<double> rdist(0.0,sqrt(0.96*sigma2(1)));
    uniform_real_distribution<double> udist(0.0,1.0);
    uniform_real_distribution<double> E1dist(0,3);
    uniform_real_distribution<double> E2dist(3,6);
    uniform_real_distribution<double> Cdist(7,10);
    
    //Patterns of Families
    //P, P+F/M, P+F+M, P+F/M+S1, P+F/M+S2, P+F+M+S1, P+F+M+S2
    //with probabilities 0.5, 0, 0.2, 0, 0, 0, 0.3
    VectorXd cump_fam (7); cump_fam << 0.5, 0.5, 0.7, 0.7, 0.7, 0.7, 1;
    VectorXi ni(nclust); VectorXi pattern(nclust); vector<VectorXi> relation_c(nclust);
    for (int i=0;i<nclust;++i){
        double u=udist(generator);
        for (int j=0;j<7;++j){
            if (u<=cump_fam(j)) {
                pattern(i) = j;
                ni(i) = j+1;
                if ((j==3)|(j==4)) ni(i) = j;
                if ((j==5)|(j==6)) ni(i) = j-1;
                break;
            }
        }
        relation_c[i].resize(ni(i));
        if (pattern(i)==0) relation_c[i] << 0;
        if (pattern(i)==1) relation_c[i] << 0,1;
        if (pattern(i)==2) relation_c[i] << 0,1,1;
        if (pattern(i)==3) relation_c[i] << 0,1,2;
        if (pattern(i)==4) relation_c[i] << 0,1,2,2;
        if (pattern(i)==5) relation_c[i] << 0,1,1,2;
        if (pattern(i)==6) relation_c[i] << 0,1,1,2,2;
    }
    long nsub = ni.sum();
    MatrixXd X(nsub,2);VectorXd L(nsub); VectorXd R(nsub); VectorXi Rinf(nsub);
    VectorXd U(nsub); VectorXi epsilon(nsub);
    VectorXi relation(nsub); VectorXi cluster(nsub); VectorXi Gene(nsub); Gene.setZero();
    
    //For proband: 1/3, 1/3, 1/3
    MatrixXi gene_proband (3,1); gene_proband << 0,1,2;
    VectorXd p_proband (3); p_proband << 1.0/3.0, 1.0/3.0, 1.0/3.0;
    gene_p proband(gene_proband, p_proband);
    
    int line = 0; double b, giF, giM, gij;
    double expb, T, D, E, C; VectorXd rij, A, rijstar;
    int nL=0; int nI = 0; int nR = 0; int nNo = 0; int nE0 = 0; int nE1 = 0; int nE2 = 0;
    for (int i=0;i<nclust;++i){
        b = bdist(generator);
        VectorXi giP = gen_gene (proband, generator);
        gene_p parent = parent_poss (giP(0),rho_pop);
        VectorXi giFM = gen_gene (parent, generator);
        giF = giFM(0); giM = giFM(1);
        bool par1 = false;
        cout<<ni(i)<<" subjects in cluster "<<i<<endl;
        int nint = 0;
        rij.resize(ni(i)); A.resize(ni(i)); rijstar.resize(ni(i));
        for (int j=0;j<ni(i);++j) A(j) = rdist(generator);
        rijstar(0) = A(0);
        if (ni(i)>=3){
            rijstar(1) = A(0)*23.0/48.0 + A(1) * 5 * sqrt(71.0)/48.0;
            rijstar(2) = A(0)*23.0/48.0 - A(1) * 125.0 * sqrt(71.0)/3408.0 + A(2) * 5.0 * sqrt(4899.0)/426.0;
        }
        if (ni(i)>=5){
            rijstar(3) = A(0)*7.0/32.0 + A(1) * 115.0 * sqrt(71.0)/2272.0 + A(2) * 5.0 * sqrt(4899.0)/568.0 + A(3) * 5.0/8.0;
            rijstar(4) = A(0)*7.0/32.0 + A(1) * 115.0 * sqrt(71.0)/2272.0 + A(2) * 5.0 * sqrt(4899.0)/568.0 - A(3) * 5.0/8.0;
        }
        for (int j=0;j<ni(i);++j) rij(j) = rijstar(j) + 0.2 * b * sqrt(sigma2(1) / sigma2(0));
        for (int j=0;j<ni(i);++j){
            X(line,0) = X1dist(generator);
            X(line,1) = X2dist(generator);
            relation(line) = relation_c[i](j);
            cluster(line) = i;
            if (relation(line)==0){
                gij = giP(0); Gene(line) = gij;
            } else if (relation(line)==1){ // parents
                if (par1) {
                    gij = giM;
                } else {
                    gij = giF; par1 = true;
                }
            } else {
                gene_p child = poss_parent (giF, giM);
                VectorXi gic = gen_gene (child, generator);
                gij = gic(0);
            }
            //simulation for T
            expb=exp(X.row(line) * beta.tail(2) + gij * beta(0) + b + rij(j));
            T = genlog1aT(0.8,expb,generator);
            //simulation for D
            expb=exp(X.row(line) * gamma.tail(2) + gij * gamma(0) + rho(0) * b + rho(1) * rij(j));
            D = genaT(0.2,expb,generator);
            
            if ((relation(line) ==0)|(relation(line) ==2)) E = E1dist(generator); else E = E2dist(generator);
            if (relation(line) ==0){
                C = Cdist(generator);
                if (T<=E) {U(line) = E; epsilon(line) = 2; nE2++;}
                else if (T<=C) {U(line) = T; epsilon(line) = 1;nE1++;}
                else { U(line) = C; epsilon(line) = 0;nE0++;}
                
                vector<double> CC=genC(C, 0, 2, generator, nint);
                long nC = CC.size();
                if (nC==0){ L(line)=0; R(line)=1; Rinf(line)=1; nNo++;}
                else {
                    if (D<CC[0]) {L(line)=0; R(line)=CC[0]; Rinf(line)=0; nL++;}
                    else {
                        int j=0;
                        while (j<nC){
                            if (D>CC[j]) {L(line)=CC[j]; R(line)=1; Rinf(line)=1;}
                            else {L(line)=CC[j-1]; R(line)=CC[j]; Rinf(line)=0; break;}
                            j++;
                        }
                        if (Rinf(line)==1) nR++; else nI++;
                    }
                }
            } else{
                U(line) = E;
                if (T<=E) { epsilon(line) = 2; nE2++; } else { epsilon(line) = 0; nE0++; }
                if (D<=E) { L(line) = 0; R(line) = E; Rinf(line) = 0; nL++;}
                else {L(line) = E; R(line) = 1; Rinf(line) = 1; nR++;}
            }
            line++;
        }
    }
    Gene_Fam Fam_info(relation, cluster, Gene, rho_pop);
    BiCens_Fam data(L, R, Rinf, U, epsilon, X, Fam_info, ngrid);
    return(data);
}

BiCens_Fam simu_data_mis_gamma(long nclust, VectorXd beta, VectorXd gamma, VectorXd rho, VectorXd sigma2, double rho_pop, long ngrid, default_random_engine & generator){
    //Distributions
    bernoulli_distribution X1dist(0.5);
    normal_distribution<double> X2dist(0.0,1.0);
    bernoulli_distribution Bdist(0.5);
    double gamma_alpha = 1.0/(exp(sigma2(0))-1.0);
    double gamma_beta = exp(sigma2(0)/2.0)*(exp(sigma2(0))-1.0);
    gamma_distribution<double> expbdist(gamma_alpha,gamma_beta);
    normal_distribution<double> rdist(0.0,sqrt(sigma2(1)));
    uniform_real_distribution<double> udist(0.0,1.0);
    uniform_real_distribution<double> E1dist(0,3);
    uniform_real_distribution<double> E2dist(3,6);
    uniform_real_distribution<double> Cdist(7,10);
    
    //Patterns of Families
    //P, P+F/M, P+F+M, P+F/M+S1, P+F/M+S2, P+F+M+S1, P+F+M+S2
    //with probabilities 0.5, 0, 0.2, 0, 0, 0, 0.3
    VectorXd cump_fam (7); cump_fam << 0.5, 0.5, 0.7, 0.7, 0.7, 0.7, 1;
    VectorXi ni(nclust); VectorXi pattern(nclust); vector<VectorXi> relation_c(nclust);
    for (int i=0;i<nclust;++i){
        double u=udist(generator);
        for (int j=0;j<7;++j){
            if (u<=cump_fam(j)) {
                pattern(i) = j;
                ni(i) = j+1;
                if ((j==3)|(j==4)) ni(i) = j;
                if ((j==5)|(j==6)) ni(i) = j-1;
                break;
            }
        }
        relation_c[i].resize(ni(i));
        if (pattern(i)==0) relation_c[i] << 0;
        if (pattern(i)==1) relation_c[i] << 0,1;
        if (pattern(i)==2) relation_c[i] << 0,1,1;
        if (pattern(i)==3) relation_c[i] << 0,1,2;
        if (pattern(i)==4) relation_c[i] << 0,1,2,2;
        if (pattern(i)==5) relation_c[i] << 0,1,1,2;
        if (pattern(i)==6) relation_c[i] << 0,1,1,2,2;
    }
    long nsub = ni.sum();
    MatrixXd X(nsub,2);VectorXd L(nsub); VectorXd R(nsub); VectorXi Rinf(nsub);
    VectorXd U(nsub); VectorXi epsilon(nsub);
    VectorXi relation(nsub); VectorXi cluster(nsub); VectorXi Gene(nsub); Gene.setZero();
    
    //For proband: 1/3, 1/3, 1/3
    MatrixXi gene_proband (3,1); gene_proband << 0,1,2;
    VectorXd p_proband (3); p_proband << 1.0/3.0, 1.0/3.0, 1.0/3.0;
    gene_p proband(gene_proband, p_proband);
    
    int line = 0; double b, riF, riM, rij, giF, giM, gij;
    double expb, T, D, E, C;
    int nL=0; int nI = 0; int nR = 0; int nNo = 0; int nE0 = 0; int nE1 = 0; int nE2 = 0;
    for (int i=0;i<nclust;++i){
        b = log(expbdist(generator)); riF = rdist(generator); riM = rdist(generator);
        VectorXi giP = gen_gene (proband, generator);
        gene_p parent = parent_poss (giP(0),rho_pop);
        VectorXi giFM = gen_gene (parent, generator);
        giF = giFM(0); giM = giFM(1);
        bool par1 = false;
        cout<<ni(i)<<" subjects in cluster "<<i<<endl;
        int nint = 0;
        for (int j=0;j<ni(i);++j){
            X(line,0) = X1dist(generator);
            X(line,1) = X2dist(generator);
            relation(line) = relation_c[i](j);
            cluster(line) = i;
            if (relation(line)==0){
                gij = giP(0); Gene(line) = gij;
                double B = Bdist(generator); rij = riF * B + riM * (1-B);
            } else if (relation(line)==1){ // parents
                if (par1) {
                    rij = riM; gij = giM;
                } else {
                    rij = riF; gij = giF; par1 = true;
                }
            } else {
                double B = Bdist(generator); rij = riF * B + riM * (1-B);
                gene_p child = poss_parent (giF, giM);
                VectorXi gic = gen_gene (child, generator);
                gij = gic(0);
            }
            //simulation for T
            expb=exp(X.row(line) * beta.tail(2) + gij * beta(0) + b + rij);
            T = genlog1aT(0.8,expb,generator);
            //simulation for D
            expb=exp(X.row(line) * gamma.tail(2) + gij * gamma(0) + rho(0) * b + rho(1) * rij);
            D = genaT(0.2,expb,generator);
            
            if ((relation(line) ==0)|(relation(line) ==2)) E = E1dist(generator); else E = E2dist(generator);
            if (relation(line) ==0){
                C = Cdist(generator);
                if (T<=E) {U(line) = E; epsilon(line) = 2; nE2++;}
                else if (T<=C) {U(line) = T; epsilon(line) = 1;nE1++;}
                else { U(line) = C; epsilon(line) = 0;nE0++;}
                
                vector<double> CC=genC(C, 0, 2, generator, nint);
                long nC = CC.size();
                if (nC==0){ L(line)=0; R(line)=1; Rinf(line)=1; nNo++;}
                else {
                    if (D<CC[0]) {L(line)=0; R(line)=CC[0]; Rinf(line)=0; nL++;}
                    else {
                        int j=0;
                        while (j<nC){
                            if (D>CC[j]) {L(line)=CC[j]; R(line)=1; Rinf(line)=1;}
                            else {L(line)=CC[j-1]; R(line)=CC[j]; Rinf(line)=0; break;}
                            j++;
                        }
                        if (Rinf(line)==1) nR++; else nI++;
                    }
                }
            } else{
                U(line) = E;
                if (T<=E) { epsilon(line) = 2; nE2++; } else { epsilon(line) = 0; nE0++; }
                if (D<=E) { L(line) = 0; R(line) = E; Rinf(line) = 0; nL++;}
                else {L(line) = E; R(line) = 1; Rinf(line) = 1; nR++;}
            }
            line++;
        }
    }
    Gene_Fam Fam_info(relation, cluster, Gene, rho_pop);
    BiCens_Fam data(L, R, Rinf, U, epsilon, X, Fam_info, ngrid);
    return(data);
}

BiCens_Fam_pred simu_data_pred_mis_gamma(long nclust, VectorXd beta, VectorXd gamma, VectorXd rho, VectorXd sigma2, double rho_pop, long ngrid, default_random_engine & generator, double t1, double t2){
    // All information up to time t1 to predict survival at time t2
    //Distributions
    bernoulli_distribution X1dist(0.5);
    normal_distribution<double> X2dist(0.0,1.0);
    bernoulli_distribution Bdist(0.5);
    double gamma_alpha = 1.0/(exp(sigma2(0))-1.0);
    double gamma_beta = exp(sigma2(0)/2.0)*(exp(sigma2(0))-1.0);
    gamma_distribution<double> expbdist(gamma_alpha,gamma_beta);
    normal_distribution<double> rdist(0.0,sqrt(sigma2(1)));
    uniform_real_distribution<double> udist(0.0,1.0);
    uniform_real_distribution<double> E1dist(0,3);
    uniform_real_distribution<double> E2dist(3,6);
    uniform_real_distribution<double> Cdist(7,10);
    
    //Patterns of Families
    //P, P+F/M, P+F+M, P+F/M+S1, P+F/M+S2, P+F+M+S1, P+F+M+S2
    //with probabilities 0.5, 0, 0.2, 0, 0, 0, 0.3
    VectorXd cump_fam (7); cump_fam << 0.5, 0.5, 0.7, 0.7, 0.7, 0.7, 1;
    VectorXi ni(nclust); VectorXi pattern(nclust); vector<VectorXi> relation_c(nclust);
    for (int i=0;i<nclust;++i){
        double u=udist(generator);
        for (int j=0;j<7;++j){
            if (u<=cump_fam(j)) {
                pattern(i) = j;
                ni(i) = j+1;
                if ((j==3)|(j==4)) ni(i) = j;
                if ((j==5)|(j==6)) ni(i) = j-1;
                break;
            }
        }
        relation_c[i].resize(ni(i));
        if (pattern(i)==0) relation_c[i] << 0;
        if (pattern(i)==1) relation_c[i] << 0,1;
        if (pattern(i)==2) relation_c[i] << 0,1,1;
        if (pattern(i)==3) relation_c[i] << 0,1,2;
        if (pattern(i)==4) relation_c[i] << 0,1,2,2;
        if (pattern(i)==5) relation_c[i] << 0,1,1,2;
        if (pattern(i)==6) relation_c[i] << 0,1,1,2,2;
    }
    long nsub = ni.sum();
    MatrixXd X(nsub,2);VectorXd L(nsub); VectorXd R(nsub); VectorXi Rinf(nsub);
    VectorXd U(nsub); VectorXi epsilon(nsub);
    VectorXi relation(nsub); VectorXi cluster(nsub); VectorXi Gene(nsub); Gene.setZero();
    
    //For proband: 1/3, 1/3, 1/3
    MatrixXi gene_proband (3,1); gene_proband << 0,1,2;
    VectorXd p_proband (3); p_proband << 1.0/3.0, 1.0/3.0, 1.0/3.0;
    gene_p proband(gene_proband, p_proband);
    
    int line = 0; double b, riF, riM, rij, giF, giM, gij;
    VectorXd T_prob(nclust); VectorXd PT_prob(nclust);
    VectorXi Ind_event_t1(nclust); VectorXi Ind_event_t2(nclust);
    double expb, T, D, E, C;
    for (int i=0;i<nclust;++i){
        b = log(expbdist(generator)); riF = rdist(generator); riM = rdist(generator);
        VectorXi giP = gen_gene (proband, generator);
        gene_p parent = parent_poss (giP(0),rho_pop);
        VectorXi giFM = gen_gene (parent, generator);
        giF = giFM(0); giM = giFM(1);
        bool par1 = false;
        cout<<ni(i)<<" subjects in cluster "<<i<<endl;
        int nint = 0;
        for (int j=0;j<ni(i);++j){
            X(line,0) = X1dist(generator);
            X(line,1) = X2dist(generator);
            relation(line) = relation_c[i](j);
            cluster(line) = i;
            if (relation(line)==0){
                gij = giP(0); Gene(line) = gij;
                double B = Bdist(generator); rij = riF * B + riM * (1-B);
            } else if (relation(line)==1){ // parents
                if (par1) {
                    rij = riM; gij = giM;
                } else {
                    rij = riF; gij = giF; par1 = true;
                }
            } else {
                double B = Bdist(generator); rij = riF * B + riM * (1-B);
                gene_p child = poss_parent (giF, giM);
                VectorXi gic = gen_gene (child, generator);
                gij = gic(0);
            }
            //simulation for T
            expb=exp(X.row(line) * beta.tail(2) + gij * beta(0) + b + rij);
            T = genlog1aT(0.8,expb,generator);
            if (relation(line) ==0) {
                T_prob(i) = T;
                if (T<=t1) Ind_event_t1(i) = 1; else Ind_event_t1(i) = 0;
                if (T<=t2) Ind_event_t2(i) = 1; else Ind_event_t2(i) = 0;
                PT_prob(i) = exp(-log(1.0+0.8*t2)*expb) / exp(-log(1.0+0.8*t1)*expb);
            }
            //simulation for D
            expb=exp(X.row(line) * gamma.tail(2) + gij * gamma(0) + rho(0) * b + rho(1) * rij);
            D = genaT(0.2,expb,generator);
            
            if ((relation(line) ==0)|(relation(line) ==2)) E = E1dist(generator); else E = E2dist(generator);
            if (relation(line) ==0){
                C = Cdist(generator);
                U(line) = t1; if (T<=t1)  epsilon(line) = 2; else epsilon(line) = 0; // know if there is event at time t1
                
                vector<double> CC=genC(C, 0, 2, generator, nint);
                long nC = CC.size();
                if (nC==0){ L(line)=0; R(line)=1; Rinf(line)=1;}
                else if (CC[0]>t1){
                    L(line)=0; R(line)=1; Rinf(line)=1;
                } else {
                    if (D<CC[0]) {L(line)=0; R(line)=CC[0]; Rinf(line)=0;}
                    else {
                        int j=0;
                        while ((j<nC)&(CC[j]<=t1)){
                            if (D>CC[j]) {L(line)=CC[j]; R(line)=1; Rinf(line)=1;}
                            else {L(line)=CC[j-1]; R(line)=CC[j]; Rinf(line)=0; break;}
                            j++;
                        }
                    }
                }
            } else{
                U(line) = E;
                    if (T<=E) { epsilon(line) = 2;} else { epsilon(line) = 0;}
                    if (D<=E) { L(line) = 0; R(line) = E; Rinf(line) = 0;}
                    else {L(line) = E; R(line) = 1; Rinf(line) = 1;}
            }
            line++;
        }
    }
    Gene_Fam Fam_info(relation, cluster, Gene, rho_pop);
    BiCens_Fam_pred data_pred(L, R, Rinf, U, epsilon, X, Fam_info, ngrid, T_prob, PT_prob,Ind_event_t1,Ind_event_t2);
    return(data_pred);
}

#endif /* BiCens_Fam_mis_h */
