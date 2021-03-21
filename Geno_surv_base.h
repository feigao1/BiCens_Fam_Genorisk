#ifndef Geno_surv_base_h
#define Geno_surv_base_h
#include <cmath>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

// Sort VectorXd increasingly
VectorXd VecSort(VectorXd x){
    long p=x.size(); double x1;
    for (int i=1;i<p;++i){
        for(int j=0;j<p-i;++j){
            if (x(j)>=x(j+1)){x1=x(j);x(j)=x(j+1);x(j+1)=x1;}
        }
    }
    return(x);
}

// Sort VectorXd increasingly and with only unique values
VectorXd VecSortUniq(VectorXd x){
    long p=x.size(); double x1;
    for (int i=1;i<p;++i){
        for(int j=0;j<p-i;++j) {
            if (x(j)>=x(j+1)) {
                x1=x(j); x(j)=x(j+1); x(j+1)=x1;
            }
        }
    }
    VectorXd z(1); z(0)=x(0); int linez=0;
    for(int i=1;i<p;++i){
        if (x(i)>z(linez)){
            VectorXd z1=z; linez=linez+1; z.resize(linez+1);
            for (int i=0;i<linez;++i) z(i)=z1(i);
            z(linez)=x(i);
        }
    }
    return(z);
}

//Get the distinct increasing jump points for interval-censored data
VectorXd get_distinct_time_i(VectorXd L, VectorXd R, VectorXi Rinf){
    long n1=L.size(); long n2=R.size();
    VectorXd x(n1+n2+1); x.head(n1)=L; int line=(int) n1;
    for (int i=0;i<n2;++i) {
        if (Rinf(i)==0) {
            x(line)=R(i); line++;
        }
    }
    x(line) = 0; line++;//add zero to the jump set
    VectorXd x1=x.head(line); return(VecSortUniq(x1));
}

//Get the distinct increasing jump points for right-censored data
VectorXd get_distinct_time_r(VectorXd Y, VectorXi Delta){
    long n=Y.size(); VectorXd x(n+1);int line=0;
    for (int i=0;i<n;++i) {
        if (Delta(i)==1) {
            x(line)=Y(i); line++;
        }
    }
    x(line) = 0; line++;//add zero to the jump set
    VectorXd x1=x.head(line); return(VecSortUniq(x1));
}

//Get the distinct increasing jump points for mixed-censored data
// Delta = 0, right; Delta = 1, exact; Delta = 2, left
VectorXd get_distinct_time_m(VectorXd Y, VectorXi Delta){
    long n=Y.size(); VectorXd x(n+1);int line=0;
    for (int i=0;i<n;++i) {
        if (Delta(i)>0) {
            x(line)=Y(i); line++;
        }
    }
    x(line) = 0; line++;//add zero to the jump set
    VectorXd x1=x.head(line); return(VecSortUniq(x1));
}

// Compute the grid and weight for Hermite integration
// Matrix: ngrid*2: weight, grid
MatrixXd Hermite_grid (long ngrid){
    MatrixXd J(ngrid,ngrid); J.setZero();
    for (int j=0;j<ngrid-1;++j){
        J(j,j+1)=pow(((double)j+1.0)/2.0,0.5);
        J(j+1,j)=J(j,j+1);
    }
    EigenSolver<MatrixXd> es(J);
    VectorXcd x=es.eigenvalues();
    MatrixXcd y=es.eigenvectors();
    
    MatrixXd grid(ngrid,2);
    for (int j=0;j<ngrid;++j){
        grid(j,0)=pow(y(0,j).real(),2)*pow(M_PI,0.5);
        grid(j,1)=x(j).real();
    }
    return(grid);
}

// generate a sequence of U with Uniform(c1,c2)
vector<double> genC(double tau, double c1, double c2, default_random_engine & generator,int& nint){
    uniform_real_distribution<double> gapdist(c1,c2);
    double C1=0; double C2=gapdist(generator);
    vector<double> C;
    if (C2<tau){
        while(C2<tau){
            C1=C2; C2 += 0.1+gapdist(generator);
            C.push_back(C1); nint++;
        }
    }
    return(C);
}

// generate failure time from Lambda = log(1+at)
double genlog1aT(double a, double expb, default_random_engine & generator){
    uniform_real_distribution<double> udist(0.0,1.0);
    double u=udist(generator); double Lambda = -log(u);
    return(exp(Lambda/expb)/a-1.0/a);
}

// generate failure time from Lambda = at
double genaT(double a, double expb, default_random_engine & generator){
    uniform_real_distribution<double> udist(0.0,1.0);
    double u=udist(generator); double Lambda = -log(u);
    return(Lambda/expb/a);
}

//Lambda -> lambda
VectorXd cumdiff(VectorXd x){ long n=x.size(); for (int i=(int)n-1;i>0;--i) x(i) -= x(i-1); return(x);}

//lambda -> Lambda
VectorXd cumsum(VectorXd x){ long n=x.size(); for (int i=1;i<n;++i) x(i) += x(i-1); return(x);}

//Numerical Integration, start from the first element of t
double NumInt (VectorXd t, VectorXd S, double tau){
    double value = 0; long step = t.size();
    for (int l=1;l<step;++l){
        if (t(l)<=tau) value += (t(l)-t(l-1))*S(l-1);
        else { value += (tau - t(l-1))*S(l-1); break;}
    }
    if (t(step-1)<=tau) value += (tau - t(step-1))*S(step-1);
    return (value);
}

//get line
VectorXi gettline(VectorXd t, VectorXd Y){
    long m=t.size(); long n = Y.size(); VectorXi Yline(n);
    for (int i=0;i<n;++i){
        int l=0;
        while (l<m){
            if (t(l)<=Y(i)) l++;
            else break;
        }
        Yline(i)=l-1;
    }
    return(Yline);
}

VectorXi gettline(VectorXd t, VectorXd Y, VectorXi start){
    long m=t.size(); long n = Y.size(); VectorXi Yline(n);
    for (int i=0;i<n;++i){
        int l=start(i);
        while (l<m){
            if (t(l)<=Y(i)) l++;
            else break;
        }
        Yline(i)=l-1;
    }
    return(Yline);
}

class gene_p{
private:
public:
    MatrixXi gene_; VectorXd p_; long ng_;
    gene_p() = default;
    gene_p(MatrixXi gene, VectorXd p){
        gene_ = gene; p_ = p; ng_ = p_.size();
    }
};

double gene_prob (int g, double rho){
    if (g==0) return(pow(1-rho,2));
    else if (g==1) return(2*rho*(1-rho));
    else return(pow(rho,2));
}

double cond_parent (int gF, int gM, int g){ // Pr(g | gF, gM)
    double p = 0;
    if (gF+gM==4){if (g==2) p=1;} //AA and AA --> AA
    else if (gF+gM==3){if (g>0) p=0.5;} //AA and Aa --> AA or Aa
    else if (gF+gM==1){if (g<2) p=0.5;} //Aa and aa --> Aa or aa
    else if (gF+gM==0){if (g==0) p=1;} //aa and aa --> aa
    else{ //in total 2 mutations
        if (gF*gM==0){if (g==1) p=1;} //AA and aa --> Aa
        else {if (g==1) p=0.5; else p=0.25;} // Aa and Aa --> AA, Aa, or AA
    }
    return(p);
}

double parent_cond (int g, int gF, int gM, double rho){ // Pr (gF, gM | g)
    double numer = gene_prob(gF,rho) * gene_prob(gM,rho) * cond_parent(gF, gM, g);
    double denom = gene_prob(g,rho);
    return(numer/denom);
}

gene_p parent_poss (int g,double rho){ // given proband, Gene matrix for parents
    long ng = 4; MatrixXi gene; gene.resize(ng,2);
    if (g==0) { gene.col(0)<<0,0,1,1; gene.col(1)<<0,1,0,1;}
    else if (g==1) { ng = 7; gene.resize(ng,2); gene.col(0)<<0,0,1,1,1,2,2; gene.col(1)<<1,2,0,1,2,0,1;}
    else { gene.col(0)<<1,1,2,2; gene.col(1)<<1,2,1,2;}
    VectorXd p(ng);
    for (int j=0;j<ng;++j){ p(j) = parent_cond(g,gene(j,0),gene(j,1),rho);}
    return(gene_p(gene,p));
}

gene_p poss_parent (int gF, int gM){ // given parents, Gene matrix for children
    long ng = 1; MatrixXi gene; gene.resize(ng,1);
    if (gF+gM==4) gene.col(0)<<2; //AA and AA --> AA
    else if (gF+gM==0) gene(0)=0; //aa and aa --> aa
    else if (gF+gM==3){ng = 2; gene.resize(ng,1); gene.col(0)<< 1,2;} //AA and Aa --> AA or Aa
    else if (gF+gM==1){ng = 2; gene.resize(ng,1); gene.col(0)<< 0,1;} //Aa and aa --> Aa or aa
    else{ //in total 2 mutations
        if (gF*gM==0) gene(0)=1; //AA and aa --> Aa
        else {ng = 3; gene.resize(ng,1); gene.col(0)<< 0,1,2;} // Aa and Aa --> AA, Aa, or AA
    }
    VectorXd p(ng);
    for (int j=0;j<ng;++j){ p(j) = cond_parent(gF, gM, gene(j,0));}
    return(gene_p(gene,p));
}

VectorXi gen_gene (gene_p genep, default_random_engine & generator){
    if (genep.ng_>1){
        uniform_real_distribution<double> udist(0.0,1.0);
        double u = udist(generator); double cut = 0;
        for (int j=0;j<genep.ng_;++j) {
            cut += genep.p_(j);
            if (u<=cut) return(genep.gene_.row(j));
        }
    }
    return(genep.gene_.row(genep.ng_-1));
}

// A class for the information of a family
// Values ordered y clust
class Gene_Fam {
private:
    void Cal_clust(VectorXi cluster){
        nclust_=1; cluster_.resize(nsub_); cluster_.setZero();
        VectorXi name(1); ni_.resize(1); name(0)=cluster(0); ni_(0) = 1;
        for(int ij=1;ij<nsub_;++ij){
            for (int i=0;i<nclust_;++i){
                if (cluster(ij)==name(i)) { ni_(i)++; cluster_(ij) = i; break;}
                else if ((i==nclust_-1)&(cluster(ij)!=name(i))) {
                    VectorXi name1 = name; name.resize(nclust_+1); name.head(nclust_) = name1;
                    VectorXi ni1 = ni_; ni_.resize(nclust_+1); ni_.head(nclust_) = ni1;
                    name(nclust_)=cluster(ij); cluster_(ij) = (int) nclust_;
                    ni_(nclust_)=1; nclust_++; break;
                }
            }
        }
        cluster_id_.resize(nclust_); proband.resize(nclust_);
        parent.resize(nclust_); sibling.resize(nclust_);
        npar_.resize(nclust_); nsib_.resize(nclust_); nchild_.resize(nclust_);
        npar_.setZero(); nsib_.setZero(); nchild_.setZero();
        for (int i=0;i<nclust_;++i) cluster_id_[i].resize(ni_(i));
        VectorXi line(nclust_); line.setZero();
        VectorXi pline(nclust_); pline.setZero();
        VectorXi sline(nclust_); sline.setZero();
        for (int ij=0;ij<nsub_;++ij){
            int i = cluster_(ij);
            cluster_id_[i](line(i)) = ij;
            if (relation_(ij)==0) proband(i) = line(i);
            else if (relation_(ij)==1) { npar_(i)++;
                if (pline(i)==0) {parent[i].resize(1); parent[i](0) = line(i);}
                else{
                    VectorXi parenti = parent[i]; parent[i].resize(pline(i)+1);
                    parent[i].head(pline(i)) = parenti; parent[i](pline(i)) = line(i);
                    relation_(ij)=3;
                }
                pline(i)++;}
            else if (relation_(ij)==2) { nsib_(i)++;
                if (sline(i)==0) {sibling[i].resize(1); sibling[i](0) = line(i);}
                else{
                    VectorXi siblingi = sibling[i]; sibling[i].resize(sline(i)+1);
                    sibling[i].head(sline(i)) = siblingi; sibling[i](sline(i)) = line(i);
                }
                sline(i)++;}
            else nchild_(i)++;
            line(i)++;
        }
    }
    void Cal_Gene(VectorXi Gene, double rho){ // Calculate the probability of gene following Mendelian law
        Gene_pro_.resize(nclust_);
        Gene_par_.resize(nclust_); P_par_.resize(nclust_); ng_par_.resize(nclust_);
        Gene_sib_.resize(nclust_); P_sib_.resize(nclust_); ng_sib_.resize(nclust_);
        for (int i=0;i<nclust_;++i){
            if (npar_(i)>2) cout<<"More than 2 parents found in cluster "<<i<<endl;
            else{
                // Assume no children, only nuclear family
                if (nchild_(i)>0) cout<<nchild_(i)<< "children found in cluster "<<i<<endl;
                else {
                    Gene_pro_(i) = Gene(cluster_id_[i](proband(i)));
                    gene_p gene_parent = parent_poss(Gene_pro_(i),rho); ng_par_(i) = (int) gene_parent.ng_;
                    Gene_par_[i] = gene_parent.gene_; P_par_[i] = gene_parent.p_;
                    if (nsib_(i)>0){
                        Gene_sib_[i].resize(ng_par_(i)); P_sib_[i].resize(ng_par_(i)); ng_sib_[i].resize(ng_par_(i));
                        for (int j=0;j<ng_par_(i);++j){
                            gene_p gene_sib = poss_parent(Gene_par_[i](j,0),Gene_par_[i](j,1));
                            Gene_sib_[i][j] = gene_sib.gene_.col(0);
                            P_sib_[i][j] = gene_sib.p_;
                            ng_sib_[i](j) = (int) gene_sib.ng_;
                        }
                    }
                }
            }
        }
    }
public:
    long nsub_,nclust_; // #subjects, #cluster
    VectorXi ni_, npar_,nsib_,nchild_; //#subjects, #parents, #sibling, #children in each cluster
    vector<VectorXi> parent; VectorXi proband; vector<VectorXi> sibling; // position of parents, proband, siblings, in each cluster
    VectorXi relation_;//0-proband, 1-parent1, 2-sibling, 3-parent2, 4-child
    VectorXi cluster_; // nsub, record the cluster num of each subject
    vector<VectorXi> cluster_id_; // nclust <ni>, record the subject ID in each cluster
    VectorXi Gene_pro_; // proband's genotype
    VectorXi ng_par_; //#possible combinations of parents
    vector<MatrixXi> Gene_par_; //parents' possible genotypes, nclust <ng_par * 2>
    vector<VectorXd> P_par_; // parents' genotype probabilites, nclust <ng_par>
    vector<VectorXi> ng_sib_; //#possible combinations for each combinations of parents
    vector<vector<VectorXi>> Gene_sib_; //child's possible genotypes, nclust < ng_par <ng_sib>>
    vector<vector<VectorXd>> P_sib_; // child's genotype probabilites, nclust <ng_par <ng_sib>>
    Gene_Fam() = default;
    Gene_Fam(VectorXi relation, VectorXi cluster, VectorXi Gene, double rho){
        relation_ = relation; nsub_ = relation_.size();
        //cout<<"nsub = "<<nsub_<<endl;
        Cal_clust(cluster);
        Cal_Gene(Gene,rho);
    }
    MatrixXd Gene_calr(int ij){// subject ij
        int rel = relation_(ij); MatrixXd rij_p (1,2);
        if ((rel==0)|(rel==2)){ rij_p.resize(2,2); rij_p << 0.5, 0, 0.5, 1; }
        else if (rel==1) rij_p << 1, 0;
        else if (rel==3) rij_p << 1, 1;
        return (rij_p);
    }
    MatrixXd Gene_calg(int ij, int kFM){ // subject ij, giF, giM at kFM line
        int i = cluster_(ij); int rel = relation_(ij); MatrixXd gij_p (1,2);
        if (rel==0) { gij_p << 1, (double) Gene_pro_(i); }
        else if (rel==1) gij_p << 1, (double) Gene_par_[i](kFM,0);
        else if (rel==2){
            gij_p.resize(ng_sib_[i](kFM),2);
            for (int j=0;j<ng_sib_[i](kFM);++j){
                gij_p(j,0) = P_sib_[i][kFM](j); gij_p(j,1) = (double) Gene_sib_[i][kFM](j,0);
            }
        }
        else if (rel==3) gij_p << 1, (double) Gene_par_[i](kFM,1);
        return (gij_p);
    }
};

class Gene_Fam_knownG {
private:
    void Cal_clust(VectorXi cluster){
        nclust_=1; cluster_.resize(nsub_); cluster_.setZero();
        VectorXi name(1); ni_.resize(1); name(0)=cluster(0); ni_(0) = 1;
        for(int ij=1;ij<nsub_;++ij){
            for (int i=0;i<nclust_;++i){
                if (cluster(ij)==name(i)) { ni_(i)++; cluster_(ij) = i; break;}
                else if ((i==nclust_-1)&(cluster(ij)!=name(i))) {
                    VectorXi name1 = name; name.resize(nclust_+1); name.head(nclust_) = name1;
                    VectorXi ni1 = ni_; ni_.resize(nclust_+1); ni_.head(nclust_) = ni1;
                    name(nclust_)=cluster(ij); cluster_(ij) = (int) nclust_;
                    ni_(nclust_)=1; nclust_++; break;
                }
            }
        }
        cluster_id_.resize(nclust_); proband.resize(nclust_);
        parent.resize(nclust_); sibling.resize(nclust_);
        npar_.resize(nclust_); nsib_.resize(nclust_); nchild_.resize(nclust_);
        npar_.setZero(); nsib_.setZero(); nchild_.setZero();
        for (int i=0;i<nclust_;++i) cluster_id_[i].resize(ni_(i));
        VectorXi line(nclust_); line.setZero();
        VectorXi pline(nclust_); pline.setZero();
        VectorXi sline(nclust_); sline.setZero();
        for (int ij=0;ij<nsub_;++ij){
            int i = cluster_(ij);
            cluster_id_[i](line(i)) = ij;
            if (relation_(ij)==0) proband(i) = line(i);
            else if (relation_(ij)==1) { npar_(i)++;
                if (pline(i)==0) {parent[i].resize(1); parent[i](0) = line(i);}
                else{
                    VectorXi parenti = parent[i]; parent[i].resize(pline(i)+1);
                    parent[i].head(pline(i)) = parenti; parent[i](pline(i)) = line(i);
                    relation_(ij)=3;
                }
                pline(i)++;}
            else if (relation_(ij)==2) { nsib_(i)++;
                if (sline(i)==0) {sibling[i].resize(1); sibling[i](0) = line(i);}
                else{
                    VectorXi siblingi = sibling[i]; sibling[i].resize(sline(i)+1);
                    sibling[i].head(sline(i)) = siblingi; sibling[i](sline(i)) = line(i);
                }
                sline(i)++;}
            else nchild_(i)++;
            line(i)++;
        }
    }
public:
    long nsub_,nclust_; // #subjects, #cluster
    VectorXi ni_, npar_,nsib_,nchild_; //#subjects, #parents, #sibling, #children in each cluster
    vector<VectorXi> parent; VectorXi proband; vector<VectorXi> sibling; // position of parents, proband, siblings, in each cluster
    VectorXi relation_;//0-proband, 1-parent1, 2-sibling, 3-parent2, 4-child
    VectorXi cluster_; // nsub, record the cluster num of each subject
    vector<VectorXi> cluster_id_; // nclust <ni>, record the subject ID in each cluster
    VectorXi Gene_; // genotype
    Gene_Fam_knownG() = default;
    Gene_Fam_knownG(VectorXi relation, VectorXi cluster, VectorXi Gene, double rho){
        relation_ = relation; nsub_ = relation_.size();
        //cout<<"nsub = "<<nsub_<<endl;
        Cal_clust(cluster); Gene_ = Gene;
    }
    MatrixXd Gene_calr(int ij){// subject ij
        int rel = relation_(ij); MatrixXd rij_p (1,2);
        if ((rel==0)|(rel==2)){ rij_p.resize(2,2); rij_p << 0.5, 0, 0.5, 1; }
        else if (rel==1) rij_p << 1, 0;
        else if (rel==3) rij_p << 1, 1;
        return (rij_p);
    }
};

#endif /* Geno_surv_base_h */
