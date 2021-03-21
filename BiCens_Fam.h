#ifndef BiCens_Fam_h
#define BiCens_Fam_h
#include "Geno_surv_base.h"
#include <time.h>

class BiCens_Fam {
private:
    //Parameters
    VectorXd beta0_, gamma0_, lambda0_, alpha0_;
    Vector2d rho0_, sigma20_;
    
    //Speicification
    long ngrid_; MatrixXd Ggrid_;
    const double eps=1E-3; const int maxiter_=10000;
    
    //Observation-related
    VectorXi Lline_, Rline_, Uline_, Rstar_;
    vector<MatrixXd> rmat_; vector<vector<MatrixXd>> gmat_; VectorXd riFM_;
    VectorXd pr_; vector<VectorXd> pg_;
    
    //Variables
    VectorXd b2_i, r2_iF, r2_iM; double nrij;
    VectorXd thetaT_ij, GthetaT_ij, G2thetaT_ij, thetaD_ij, GthetaD_ij, G2thetaD_ij;
    VectorXd bthetaD_ij, b2thetaD_ij, rthetaD_ij, r2thetaD_ij, brthetaD_ij, GbthetaD_ij, GrthetaD_ij;
    MatrixXd W_ijl, WG_ijl, V_ijl, VG_ijl;
    MatrixXd Vb_ijl, Vr_ijl;
    
    VectorXd Score_; MatrixXd Hessian_;
    double sumtheta, sumthetab, sumthetabb, sumthetar, sumthetarr, sumthetabr;
    VectorXd sumthetax, sumthetabx, sumthetarx; MatrixXd sumthetaxx;
    MatrixXd pI_score;
    
    //Parameter Initiation
    void Para_init (){
        tT_ = get_distinct_time_m(U_, epsilon_); tD_ = get_distinct_time_i(L_, R_, Rinf_);
        mT_ = tT_.size(); mD_ = tD_.size();
        lambda_.resize(mT_); lambda_.setOnes(); lambda_ *= 1/(double)mT_; Lambda_ = cumsum(lambda_);
        alpha_.resize(mD_); alpha_.setOnes(); alpha_ *= 1/(double)mD_; A_ = cumsum(alpha_);
        
        beta_.resize(p_+1); beta_.setZero(); gamma_.resize(p_+1); gamma_.setZero();
        rho_.setOnes(); sigma2_.setOnes();
        
    }
    void Other_init(){
        log_lik_.resize(nclust_);
        
        b2_i.resize(nclust_); r2_iF.resize(nclust_); r2_iM.resize(nclust_);
        thetaT_ij.resize(nsub_); GthetaT_ij.resize(nsub_); G2thetaT_ij.resize(nsub_);
        thetaD_ij.resize(nsub_); GthetaD_ij.resize(nsub_); G2thetaD_ij.resize(nsub_);
        bthetaD_ij.resize(nsub_); b2thetaD_ij.resize(nsub_);
        rthetaD_ij.resize(nsub_); r2thetaD_ij.resize(nsub_); brthetaD_ij.resize(nsub_);
        GbthetaD_ij.resize(nsub_); GrthetaD_ij.resize(nsub_);
        
        W_ijl.resize(nsub_,mT_); WG_ijl.resize(nsub_,mT_);
        V_ijl.resize(nsub_,mD_); VG_ijl.resize(nsub_,mD_);
        Vb_ijl.resize(nsub_,mD_); Vr_ijl.resize(nsub_,mD_);
        
        sumthetax.resize(p_+1); sumthetabx.resize(p_+1); sumthetarx.resize(p_+1); sumthetaxx.resize(p_+1,p_+1);
        
        rmat_.resize(nsub_); gmat_.resize(nsub_); pr_.resize(nsub_); pg_.resize(nsub_); riFM_.resize(2);
        for (int ij=0;ij<nsub_;++ij){
            int i = Fam_info_.cluster_(ij);
            rmat_[ij] = Fam_info_.Gene_calr(ij); pr_(ij) = rmat_[ij].rows();
            gmat_[ij].resize(Fam_info_.ng_par_(i)); pg_[ij].resize(Fam_info_.ng_par_(i));
            for (int kFM=0;kFM <Fam_info_.ng_par_(i);++kFM){
                gmat_[ij][kFM] = Fam_info_.Gene_calg(ij, kFM); pg_[ij](kFM) = gmat_[ij][kFM].rows();
            }
        }
        Uline_= gettline(tT_,U_);
        Lline_ = gettline(tD_,L_); Rline_ = gettline(tD_,R_,Lline_); Rstar_ = Rline_;
        for (int i=0;i<nsub_;++i){ if (Rinf_(i)==1) { Rline_(i) = (int) mD_; Rstar_(i) = Lline_(i);}}
        nrij = (double) nclust_; for (int i=0;i<nclust_;++i) {if (ni_(i)>1) nrij++;}
    }
    
    double cal_thetaT(int ij, double bi, double rij, int gij, VectorXd beta){
        double thetaT = beta(0) * (double) gij + X_.row(ij) * beta.tail(p_) + bi + rij;
        return(exp(thetaT));
    }
    
    double cal_thetaD(int ij, double bi, double rij, int gij, VectorXd gamma, Vector2d rho){
        double thetaD = gamma(0) * (double) gij + X_.row(ij) * gamma.tail(p_) + rho(0) * bi + rho(1) * rij;
        return(exp(thetaD));
    }
    
    double cal_lij(int ij, double bi, double rij, int gij, double thetaT, double thetaD, VectorXd lambda, VectorXd Lambda, VectorXd alpha, VectorXd A){
        double lij = exp(-A(Lline_(ij))*thetaD);
        if (Rinf_(ij)==0) lij -= exp(-A(Rline_(ij))*thetaD);
        if (epsilon_(ij)==2) lij *= (1-exp(-Lambda(Uline_(ij))*thetaT));
        else lij *= exp(-Lambda(Uline_(ij))*thetaT);
        if (epsilon_(ij)==1) lij *= lambda(Uline_(ij))*thetaT;
        return (lij);
    }
 
    void cal_likelihood(VectorXd beta, VectorXd gamma, Vector2d rho, Vector2d sigma2, VectorXd lambda, VectorXd Lambda, VectorXd alpha, VectorXd A){
        b2_i.setZero(); r2_iF.setZero(); r2_iM.setZero();
        thetaT_ij.setZero(); GthetaT_ij.setZero(); G2thetaT_ij.setZero();
        thetaD_ij.setZero(); GthetaD_ij.setZero(); G2thetaD_ij.setZero();
        bthetaD_ij.setZero(); b2thetaD_ij.setZero();
        rthetaD_ij.setZero(); r2thetaD_ij.setZero(); brthetaD_ij.setZero();
        GbthetaD_ij.setZero(); GrthetaD_ij.setZero();
        W_ijl.setZero(); WG_ijl.setZero();
        V_ijl.setZero(); VG_ijl.setZero(); Vb_ijl.setZero(); Vr_ijl.setZero();
        VectorXd EG_ij(nsub_); EG_ij.setZero();
        log_lik_.setZero(); double lik_FM, lij;
        double wb,wF,wM; double bi,riF,riM; double rij,gij; double wlik_rg, wlik, wlik_j;
        double ww,vv,denom;
        double thetaT,thetaD;
        for (int i=0;i<nclust_;++i){
            VectorXd thetat(ni_(i)); VectorXd gthetat(ni_(i)); VectorXd g2thetat(ni_(i));
            VectorXd thetad(ni_(i)); VectorXd gthetad(ni_(i)); VectorXd g2thetad(ni_(i));
            VectorXd bthetad(ni_(i)); VectorXd b2thetad(ni_(i));
            VectorXd rthetad(ni_(i)); VectorXd r2thetad(ni_(i)); VectorXd brthetad(ni_(i));
            VectorXd gbthetad(ni_(i)); VectorXd grthetad(ni_(i));
            MatrixXd w(ni_(i),mT_); MatrixXd wg(ni_(i),mT_);
            MatrixXd v(ni_(i),mD_); MatrixXd vg(ni_(i),mD_);
            MatrixXd vb(ni_(i),mD_); MatrixXd vr(ni_(i),mD_);
            VectorXd eg(ni_(i));
            for (int kb=0;kb<ngrid_;++kb){
                wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2(0),0.5);
                if (ni_(i)>1){
                    for (int kF=0;kF<ngrid_;++kF){
                        wF=Ggrid_(kF,0)*pow(M_PI,-0.5); riF=Ggrid_(kF,1)*pow(2.0*sigma2(1),0.5);
                        for (int kM=0;kM<ngrid_;++kM){
                            wM=Ggrid_(kM,0)*pow(M_PI,-0.5); riM=Ggrid_(kM,1)*pow(2.0*sigma2(1),0.5);
                            riFM_(0) = riF; riFM_(1) = riM;
                            for (int kFM=0;kFM <Fam_info_.ng_par_(i);++kFM){
                                lik_FM = 1;
                                thetat.setZero(); gthetat.setZero(); g2thetat.setZero();
                                thetad.setZero(); gthetad.setZero(); g2thetad.setZero();
                                bthetad.setZero(); b2thetad.setZero();
                                rthetad.setZero(); r2thetad.setZero(); brthetad.setZero();
                                gbthetad.setZero(); grthetad.setZero();
                                w.setZero(); wg.setZero();
                                v.setZero(); vg.setZero(); vb.setZero(); vr.setZero();
                                eg.setZero();
                                VectorXd lijm(ni_(i)); lijm.setOnes();
                                for (int j=0;j<ni_(i);++j){
                                    int ij = Fam_info_.cluster_id_[i](j);
                                    lij = 0;
                                    for (int r=0;r<pr_(ij);++r){
                                        for (int g=0;g<pg_[ij](kFM);++g){
                                            rij = riFM_((int)rmat_[ij](r,1)); gij = gmat_[ij][kFM](g,1);
                                            thetaT = cal_thetaT(ij, bi, rij, gij, beta);
                                            thetaD = cal_thetaD(ij, bi, rij, gij, gamma, rho);
                                            wlik_rg = rmat_[ij](r,0) * gmat_[ij][kFM](g,0) * cal_lij(ij, bi, rij, gij, thetaT, thetaD, lambda, Lambda, alpha, A);
                                            lij += wlik_rg;
                                            thetat(j) += wlik_rg * thetaT;
                                            gthetat(j) += wlik_rg * gij * thetaT; g2thetat(j) += wlik_rg * pow(gij,2) * thetaT;
                                            thetad(j) += wlik_rg * thetaD;
                                            gthetad(j) += wlik_rg * gij * thetaD; g2thetad(j) += wlik_rg * pow(gij,2) * thetaD;
                                            bthetad(j) += wlik_rg * bi * thetaD; b2thetad(j) += wlik_rg * pow(bi,2) * thetaD;
                                            rthetad(j) += wlik_rg * rij * thetaD; r2thetad(j) += wlik_rg * pow(rij,2) * thetaD;
                                            brthetad(j) += wlik_rg * bi * rij * thetaD;
                                            gbthetad(j) += wlik_rg * gij * bi * thetaD; grthetad(j) += wlik_rg * gij * rij * thetaD;
                                            eg(j) += wlik_rg * gij;
                                            
                                            if (epsilon_(ij)==2){
                                                denom = 1-exp(-Lambda(Uline_(ij))*thetaT);
                                                for (int l=0;l<=Uline_(ij);++l) {
                                                    ww = lambda(l) * thetaT / denom; if (denom==0) ww = 0;
                                                    w(j,l) += wlik_rg * ww; wg(j,l) += wlik_rg * ww * gij;
                                                }
                                            }
                                            if (Rinf_(ij)==0){
                                                denom = 1-exp(-(A(Rline_(ij))-A(Lline_(ij)))*thetaD);
                                                for (int l=(Lline_(ij)+1);l<=Rline_(ij);++l) {
                                                    vv = alpha(l) * thetaD / denom; if (denom==0) vv = 0;
                                                    v(j,l) += wlik_rg * vv; vg(j,l) += wlik_rg * vv * gij;
                                                    vb(j,l) += wlik_rg * vv * bi; vr(j,l) += wlik_rg * vv * rij;
                                                }
                                            }
                                        }
                                    }
                                    lik_FM *= lij;
                                    for (int jj=0;jj<ni_(i);++jj){if (jj!=j) lijm(jj) *= lij;}
                                }
                                wlik = wb * wF * wM * Fam_info_.P_par_[i](kFM) * lik_FM;
                                log_lik_(i) += wlik;
                                b2_i(i) += wlik * pow(bi,2); r2_iF(i) += wlik * pow(riF,2); r2_iM(i) += wlik * pow(riM,2);
                                
                                for (int j=0;j<ni_(i);++j){
                                    wlik_j = wb * wF * wM * Fam_info_.P_par_[i](kFM) * lijm(j);
                                    int ij = Fam_info_.cluster_id_[i](j);
                                    thetaT_ij(ij) += wlik_j * thetat(j);
                                    GthetaT_ij(ij) += wlik_j * gthetat(j); G2thetaT_ij(ij) += wlik_j * g2thetat(j);
                                    thetaD_ij(ij) += wlik_j * thetad(j);
                                    GthetaD_ij(ij) += wlik_j * gthetad(j); G2thetaD_ij(ij) += wlik_j * g2thetad(j);
                                    bthetaD_ij(ij) += wlik_j * bthetad(j); b2thetaD_ij(ij) += wlik_j * b2thetad(j);
                                    rthetaD_ij(ij) += wlik_j * rthetad(j); r2thetaD_ij(ij) += wlik_j * r2thetad(j);
                                    brthetaD_ij(ij) += wlik_j * brthetad(j);
                                    GbthetaD_ij(ij) += wlik_j * gbthetad(j); GrthetaD_ij(ij) += wlik_j * grthetad(j);
                                    W_ijl.row(ij) += wlik_j * w.row(j); WG_ijl.row(ij) += wlik_j * wg.row(j);
                                    V_ijl.row(ij) += wlik_j * v.row(j); VG_ijl.row(ij) += wlik_j * vg.row(j);
                                    Vb_ijl.row(ij) += wlik_j * vb.row(j); Vr_ijl.row(ij) += wlik_j * vr.row(j);
                                    EG_ij(ij) += wlik_j * eg(j);
                                }
                            }
                        }
                    }
                }
                else{
                    gij = (double) Fam_info_.Gene_pro_(i);
                    int ij = Fam_info_.cluster_id_[i](0);
                    for (int kF=0;kF<ngrid_;++kF){
                        wF = Ggrid_(kF,0)*pow(M_PI,-0.5); rij = Ggrid_(kF,1)*pow(2.0*sigma2(1),0.5);
                        thetaT = cal_thetaT(ij, bi, rij, gij, beta);
                        thetaD = cal_thetaD(ij, bi, rij, gij, gamma, rho);
                        wlik = wb * wF * cal_lij(ij, bi, rij, gij, thetaT, thetaD, lambda, Lambda, alpha, A);
                        
                        thetaT_ij(ij) += wlik * thetaT;
                        GthetaT_ij(ij) += wlik * gij * thetaT; G2thetaT_ij(ij) += wlik * pow(gij,2) * thetaT;
                        thetaD_ij(ij) += wlik * thetaD;
                        GthetaD_ij(ij) += wlik * gij * thetaD; G2thetaD_ij(ij) += wlik * pow(gij,2) * thetaD;
                        bthetaD_ij(ij) += wlik * bi * thetaD; b2thetaD_ij(ij) += wlik * pow(bi,2) * thetaD;
                        rthetaD_ij(ij) += wlik * rij * thetaD; r2thetaD_ij(ij) += wlik * pow(rij,2) * thetaD;
                        brthetaD_ij(ij) += wlik * bi * rij * thetaD;
                        GbthetaD_ij(ij) += wlik * gij * bi * thetaD; GrthetaD_ij(ij) += wlik * gij * rij * thetaD;
                        
                        if (epsilon_(ij)==2){
                            denom = 1-exp(-Lambda(Uline_(ij))*thetaT);
                            for (int l=0;l<=Uline_(ij);++l) {
                                ww = lambda(l) * thetaT / denom; if (denom==0) ww = 0;
                                W_ijl(ij,l) += wlik * ww; WG_ijl(ij,l) += wlik * ww * gij;
                            }
                        }
                        if (Rinf_(ij)==0){
                            denom = 1-exp(-(A(Rline_(ij))-A(Lline_(ij)))*thetaD);
                            for (int l=(Lline_(ij)+1);l<=Rline_(ij);++l) {
                                vv = alpha(l) * thetaD / denom; if (denom==0) vv = 0;
                                V_ijl(ij,l) += wlik * vv; VG_ijl(ij,l) += wlik * vv * gij;
                                Vb_ijl(ij,l) += wlik * vv * bi; Vr_ijl(ij,l) += wlik * vv * rij;
                            }
                        }
                        EG_ij(ij) += wlik * gij;
                        log_lik_(i) += wlik;
                        b2_i(i) += wlik * pow(bi,2);
                        r2_iF(i) += wlik * pow(rij,2);
                    }
                }
            }
            b2_i(i) /= log_lik_(i); r2_iF(i)/= log_lik_(i); r2_iM(i)/= log_lik_(i);
            for (int j=0;j<ni_(i);++j){
                int ij = Fam_info_.cluster_id_[i](j);
                thetaT_ij(ij) /= log_lik_(i); GthetaT_ij(ij) /= log_lik_(i); G2thetaT_ij(ij) /= log_lik_(i);
                thetaD_ij(ij) /= log_lik_(i); GthetaD_ij(ij) /= log_lik_(i); G2thetaD_ij(ij) /= log_lik_(i);
                bthetaD_ij(ij) /= log_lik_(i); b2thetaD_ij(ij) /= log_lik_(i);
                rthetaD_ij(ij) /= log_lik_(i); r2thetaD_ij(ij) /= log_lik_(i); brthetaD_ij(ij) /= log_lik_(i);
                GbthetaD_ij(ij) /= log_lik_(i); GrthetaD_ij(ij) /= log_lik_(i);
                
                W_ijl.row(ij) /= log_lik_(i); WG_ijl.row(ij) /= log_lik_(i);
                V_ijl.row(ij) /= log_lik_(i); VG_ijl.row(ij) /= log_lik_(i);
                Vb_ijl.row(ij) /= log_lik_(i); Vr_ijl.row(ij) /= log_lik_(i);
                EG_ij(ij) /= log_lik_(i);
                if (epsilon_(ij)==1) {
                    W_ijl.row(ij).setZero(); W_ijl(ij,Uline_(ij)) = 1;
                    WG_ijl.row(ij).setZero(); WG_ijl(ij,Uline_(ij)) = EG_ij(ij);
                }
            }
            log_lik_(i) = log(log_lik_(i));
        }
    }
    
    VectorXd update_T(VectorXd& beta){
        Score_.resize(p_+1); Hessian_.resize(p_+1,p_+1); Score_.setZero(); Hessian_.setZero();
        VectorXd lambda(mT_);
        VectorXd Xij; double sumw; VectorXd sumwx(p_+1);
        sumtheta=0; sumthetax.setZero(); sumthetaxx.setZero();
        for (int l=((int)mT_-1);l>=0;--l){
            for (int ij=0;ij<nsub_;++ij){
                if (l==Uline_(ij)){
                    Xij = X_.row(ij).transpose();
                    sumtheta += thetaT_ij(ij);
                    sumthetax(0) += GthetaT_ij(ij);
                    sumthetax.tail(p_) += thetaT_ij(ij) * Xij;
                    sumthetaxx(0,0) += G2thetaT_ij(ij);
                    sumthetaxx.block(0,1,1,p_) += GthetaT_ij(ij) * Xij.transpose();
                    sumthetaxx.block(1,0,p_,1) += GthetaT_ij(ij)* Xij;
                    sumthetaxx.block(1,1,p_,p_) += thetaT_ij(ij) * Xij * Xij.transpose();
                }
            }
            sumw = W_ijl.col(l).sum(); sumwx(0) = WG_ijl.col(l).sum();
            sumwx.tail(p_) = W_ijl.col(l).transpose() * X_;
            Score_ += sumwx - sumw * sumthetax/sumtheta;
            Hessian_ += sumw * (sumthetaxx/sumtheta - sumthetax*sumthetax.transpose()/pow(sumtheta,2));
            
            lambda(l) = sumw / sumtheta;
        }
        beta += Hessian_.inverse() * Score_;
        return(lambda);
    }
    
    VectorXd update_D(VectorXd& gamma, Vector2d& rho){
        Score_.resize(p_+3); Hessian_.resize(p_+3,p_+3); Score_.setZero(); Hessian_.setZero();
        VectorXd alpha(mD_);
        VectorXd Xij; double sumv, sumvb, sumvr; VectorXd sumvx(p_+1);
        sumtheta = 0; sumthetax.setZero(); sumthetaxx.setZero();
        sumthetab = 0; sumthetabx.setZero(); sumthetabb = 0;
        sumthetar = 0; sumthetarx.setZero(); sumthetarr = 0;
        sumthetabr = 0;
        for (int l=((int)mD_-1);l>=0;--l){
            for (int ij=0;ij<nsub_;++ij){
                if (l==Rstar_(ij)){
                    Xij = X_.row(ij).transpose();
                    sumtheta += thetaD_ij(ij);
                    sumthetax(0) += GthetaD_ij(ij);
                    sumthetax.tail(p_) += thetaD_ij(ij) * Xij;
                    sumthetaxx(0,0) += G2thetaD_ij(ij);
                    sumthetaxx.block(0,1,1,p_) += GthetaD_ij(ij) * Xij.transpose();
                    sumthetaxx.block(1,0,p_,1) += GthetaD_ij(ij)* Xij;
                    sumthetaxx.block(1,1,p_,p_) += thetaD_ij(ij) * Xij * Xij.transpose();
                    
                    sumthetab += bthetaD_ij(ij); sumthetabx(0) += GbthetaD_ij(ij);
                    sumthetabx.tail(p_) += bthetaD_ij(ij) * Xij;
                    sumthetabb += b2thetaD_ij(ij);
                    
                    sumthetar += rthetaD_ij(ij); sumthetarx(0) += GrthetaD_ij(ij);
                    sumthetarx.tail(p_) += rthetaD_ij(ij) * Xij;
                    sumthetarr += r2thetaD_ij(ij);
                    
                    sumthetabr += brthetaD_ij(ij);
                }
            }
            sumv = V_ijl.col(l).sum(); sumvx(0) = VG_ijl.col(l).sum();
            sumvx.tail(p_) = V_ijl.col(l).transpose() * X_;
            sumvb = Vb_ijl.col(l).sum(); sumvr = Vr_ijl.col(l).sum();
            
            Score_.head(p_+1) += sumvx - sumv * sumthetax / sumtheta;
            Score_(p_+1) += sumvb - sumv * sumthetab / sumtheta;
            Score_(p_+2) += sumvr - sumv * sumthetar / sumtheta;
            
            Hessian_.block(0,0,p_+1,p_+1) += sumv * (sumthetaxx / sumtheta - sumthetax*sumthetax.transpose()/pow(sumtheta,2));
            Hessian_.block(0,p_+1,p_+1,1) += sumv * (sumthetabx / sumtheta - sumthetax * sumthetab / pow(sumtheta,2));
            Hessian_.block(0,p_+2,p_+1,1) += sumv * (sumthetarx / sumtheta - sumthetax * sumthetar / pow(sumtheta,2));
            Hessian_(p_+1,p_+1) += sumv * (sumthetabb / sumtheta - pow(sumthetab / sumtheta,2));
            Hessian_(p_+2,p_+2) += sumv * (sumthetarr / sumtheta - pow(sumthetar / sumtheta,2));
            Hessian_(p_+1,p_+2) += sumv * (sumthetabr / sumtheta - sumthetab * sumthetar / pow(sumtheta,2));
            
            alpha(l) = sumv / sumtheta;
        }
        Hessian_.block(p_+1,0,2,p_+1) = Hessian_.block(0,p_+1,p_+1,2).transpose();
        Hessian_(p_+2,p_+1) = Hessian_(p_+1,p_+2);
        VectorXd step = Hessian_.inverse() * Score_;
        gamma += step.head(p_+1);
        rho += step.tail(2);
        return(alpha);
    }
    
    VectorXd update_lambda(){
        VectorXd lambda(mT_);
        for (int l=0;l<mT_;++l){
            sumtheta = 0;
            for (int ij=0;ij<nsub_;++ij){
                if (l<=Uline_(ij)) sumtheta += thetaT_ij(ij);
            }
            lambda(l) = W_ijl.col(l).sum() / sumtheta; if (sumtheta ==0) lambda(l)=0;
        }
        return(lambda);
    }
    
    VectorXd update_alpha(){
        VectorXd alpha(mD_);
        for (int l=0;l<mD_;++l){
            sumtheta = 0;
            for (int ij=0;ij<nsub_;++ij){
                if (l<=Rstar_(ij)) sumtheta += thetaD_ij(ij);
            }
            alpha(l) = V_ijl.col(l).sum() / sumtheta; if (sumtheta ==0) alpha(l)=0;
        }
        return(alpha);
    }
    
    Vector2d update_sigma(){
        Vector2d sigma2; sigma2 << b2_i.mean(), (r2_iF.sum() + r2_iM.sum())/nrij;
        return(sigma2);
    }
    
    double diff(VectorXd a, VectorXd b){
        VectorXd diff=a-b;
        if (b.norm()<=0.01) return(diff.norm());
        else return(diff.norm()/b.norm());
    }
    double diff(double a, double b){
        double diff=a-b;
        if (abs(b)<=0.01) return(abs(diff));
        else return(abs(diff)/abs(b));
    }
    
    VectorXd profilelogLik_i(VectorXd beta, VectorXd gamma, Vector2d rho, Vector2d sigma2){
        iter_=0;
        cout<<"beta="<<beta.transpose()<<endl;
        cout<<"gamma="<<gamma.transpose()<<endl;
        cout<<"rho = "<<rho.transpose()<<", sigma2 = "<<sigma2.transpose()<<endl;
        VectorXd lambda=lambda_; VectorXd lambda0, Lambda;
        VectorXd alpha=alpha_; VectorXd alpha0, A;
        while (iter_<=maxiter_){
            lambda0 = lambda; alpha0 = alpha;
            Lambda = cumsum(lambda); A = cumsum(alpha);
            cal_likelihood(beta, gamma, rho, sigma2,lambda, Lambda, alpha, A);
            lambda = update_lambda(); alpha = update_alpha();
            cout<<"iter="<<iter_<<", alphadiff="<<diff(lambda,lambda0)+diff(alpha,alpha0)<<endl;
            cout<<" logLik="<<log_lik_.mean()<<endl;
            iter_++;
            if (diff(lambda,lambda0)+diff(alpha,alpha0)<eps) break;
        }
        return(log_lik_);
    }
    
    void update_h(int j, double h, VectorXd& beta, VectorXd& gamma, Vector2d& rho, Vector2d& sigma2){
        if (j<p_+1) beta(j) += h;
        else{
            if (j<2*(p_+1)) gamma(j-(p_+1)) += h;
            else{
                if (j<2*(p_+2)) rho(j-2*(p_+1)) += h;
                else sigma2(j-2*(p_+2)) += h;
            }
        }
    }
    
    void update_h_sig(int j, double h, VectorXd& beta, VectorXd& gamma, Vector2d& rho, Vector2d& sigma2){
        if (j<p_+1) beta(j) += h;
        else{
            if (j<2*(p_+1)) gamma(j-(p_+1)) += h;
            else{
                if (j<2*(p_+2)) rho(j-2*(p_+1)) += h;
                else sigma2(j-2*(p_+2)) *= exp(h);
            }
        }
    }
    
    void Profile_I(double h){
        VectorXd beta, gamma; Vector2d rho, sigma2;
        long P = 2 * (p_+1) + 4;
        MatrixXd pL_i(nclust_,P);
        /* calculate the value of pli(beta+h e)*/
        for (int j=0;j<P;++j){
            beta=beta_; gamma = gamma_;
            rho = rho_; sigma2 = sigma2_;
            update_h(j,h,beta, gamma, rho, sigma2);
            pL_i.col(j)=profilelogLik_i(beta, gamma, rho, sigma2);
        }
        pI_score.resize(P,P); pI_score.setZero(); MatrixXd pL_i_c=pL_i;
        for (int j=0;j<P;++j) pL_i_c.col(j) = (pL_i_c.col(j)-logLik_MLE_i)/h;
        for (int i=0;i<nclust_;++i) pI_score += pL_i_c.row(i).transpose()*pL_i_c.row(i);
    }
    
    void Profile_I_sig(double h){
        VectorXd beta, gamma; Vector2d rho, sigma2;
        long P = 2 * (p_+1) + 4;
        MatrixXd pL_i(nclust_,P);
        /* calculate the value of pli(beta+h e)*/
        for (int j=0;j<P;++j){
            beta=beta_; gamma = gamma_;
            rho = rho_; sigma2 = sigma2_;
            update_h_sig(j,h,beta, gamma, rho, sigma2);
            pL_i.col(j)=profilelogLik_i(beta, gamma, rho, sigma2);
        }
        pI_score.resize(P,P); pI_score.setZero(); MatrixXd pL_i_c=pL_i;
        for (int j=0;j<P;++j) pL_i_c.col(j) = (pL_i_c.col(j)-logLik_MLE_i)/h;
        for (int i=0;i<nclust_;++i) pI_score += pL_i_c.row(i).transpose()*pL_i_c.row(i);
    }
    
public:
    //Parameters
    Vector2d rho_, sigma2_; // rho1, rho2, sigmab2, sigmar2
    VectorXd beta_, gamma_, tT_, tD_, lambda_, alpha_, Lambda_, A_;
    long mT_, mD_;
    
    //Observed Data
    VectorXd L_, R_, U_; VectorXi Rinf_, epsilon_; MatrixXd X_;
    Gene_Fam Fam_info_;
    long nsub_, nclust_, p_; VectorXi ni_;
    
    //Other output
    VectorXd log_lik_; int iter_;
    VectorXd logLik_MLE_i;
    MatrixXd Cov_score; VectorXd sd_score;
    
    BiCens_Fam() = default;
    BiCens_Fam(VectorXd L, VectorXd R, VectorXi Rinf, VectorXd U, VectorXi epsilon, MatrixXd X, Gene_Fam Fam_info, long ngrid){
        L_ = L; R_ = R; Rinf_ = Rinf; U_ = U; epsilon_ = epsilon;
        X_ = X; Fam_info_ = Fam_info;
        nsub_ = L_.size(); p_ = X_.row(0).size(); nclust_ = Fam_info_.nclust_; ni_ = Fam_info_.ni_;
        ngrid_ = ngrid; Ggrid_ = Hermite_grid (ngrid_);
        Para_init(); Other_init();
    }
    BiCens_Fam(VectorXd L, VectorXd R, VectorXi Rinf, VectorXd U, VectorXi epsilon, MatrixXd X, Gene_Fam Fam_info, long ngrid, VectorXd beta, VectorXd gamma, VectorXd tT, VectorXd tD, VectorXd lambda, VectorXd alpha, Vector2d rho, Vector2d sigma2){
        L_ = L; R_ = R; Rinf_ = Rinf; U_ = U; epsilon_ = epsilon; X_ = X; Fam_info_ = Fam_info;
        nsub_ = L_.size(); p_ = X_.row(0).size(); nclust_ = Fam_info_.nclust_; ni_ = Fam_info_.ni_;
        ngrid_ = ngrid; Ggrid_ = Hermite_grid (ngrid_);
        
        tT_ = tT; tD_ = tD; mT_ = tT_.size(); mD_ = tD_.size();
        lambda_ = lambda; Lambda_ = cumsum(lambda_);
        alpha_ = alpha; A_ = cumsum(alpha_);
        beta_ = beta; gamma_ = gamma; rho_ = rho; sigma2_ = sigma2;
        
        Other_init();
    }
    void solve(){
        iter_=0; double rhodiff, betadiff,alphadiff,sigmadiff;
        while (iter_<=maxiter_){
            beta0_=beta_;gamma0_=gamma_;
            lambda0_=lambda_; alpha0_=alpha_;
            rho0_ = rho_; sigma20_ = sigma2_;
            Lambda_ = cumsum(lambda_); A_ = cumsum(alpha_);
            cal_likelihood(beta_, gamma_, rho_, sigma2_, lambda_, Lambda_, alpha_, A_);
            lambda_ = update_T(beta_); alpha_ = update_D(gamma_,rho_);
            sigma2_ = update_sigma();
            rhodiff=diff(rho0_,rho_);
            betadiff=max(diff(beta_,beta0_),diff(gamma_,gamma0_));
            alphadiff=max(diff(lambda_,lambda0_),diff(alpha_,alpha0_));
            sigmadiff=diff(sigma20_,sigma2_);
            iter_++;
            if (rhodiff+betadiff+alphadiff+sigmadiff<eps) break;
            cout<<"iter= "<<iter_<<", beta = "<<beta_.transpose()<<endl;
            cout<<"gamma = "<<gamma_.transpose()<<endl;
            cout<<"rho = "<<rho_.transpose()<<", sigma2 = "<<sigma2_.transpose()<<endl;
            cout<<"lambda"<<lambda_.head(5).transpose()<<lambda_.tail(5).transpose();
            cout<<"alpha"<<alpha_.head(5).transpose()<<alpha_.tail(5).transpose();
            cout<<"rhodiff="<<rhodiff<<", betadiff="<<betadiff<<endl;
            cout<<"alphadiff="<<alphadiff<<", sigmadiff="<<sigmadiff<<endl;
            cout<<"logLik="<<log_lik_.sum()/(double)(nsub_)<<endl;
        }
        logLik_MLE_i=log_lik_;
        if (iter_<=maxiter_){
            cout<<"The Algorithm converges in "<<iter_<<" iterations."<<endl;
        }
        else {cout<<"Maximum Iteration Times ("<<maxiter_<<") Reached!"<<endl;}
    }
    MatrixXd cumu_T(VectorXd GX){
        VectorXd tTD_(mT_+mD_); tTD_.head(mT_) = tT_; tTD_.tail(mD_) = tD_;
        tTD_ = VecSortUniq(tTD_); long mTD_ = tTD_.size();
        VectorXi Tline_ = gettline(tT_,tTD_); VectorXi Dline_= gettline(tD_,tTD_);
        VectorXi TDline_= gettline(tT_,tD_);
        MatrixXd tcumuinc(mTD_,2); tcumuinc.setZero(); tcumuinc.col(0) = tTD_;
        
        double wb, wr, bi, ri, thetaT, thetaD;
        for (int kb=0;kb<ngrid_;++kb){
            wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
            for (int kr=0;kr<ngrid_;++kr){
                wr=Ggrid_(kr,0)*pow(M_PI,-0.5); ri=Ggrid_(kr,1)*pow(2.0*sigma2_(1),0.5);
                thetaT = exp(GX.transpose() * beta_ + bi + ri);
                thetaD = exp(GX.transpose() * gamma_ + rho_(0) * bi + rho_(1) * ri);
                for (int j=0;j<mTD_;++j){
                    tcumuinc(j,1) += wb * wr * (1-exp(-Lambda_(Tline_(j))*thetaT)) * exp(-A_(Dline_(j))*thetaD);
                    for (int jj=0;jj<=Dline_(j);++jj)
                        tcumuinc(j,1) += wb * wr * alpha_(jj) * thetaD * exp(-A_(jj)*thetaD) * (1-exp(-Lambda_(TDline_(jj))*thetaT));
                }
            }
        }
        return(tcumuinc);
    }
    MatrixXd cumu_D(VectorXd GX){
        VectorXd tTD_(mT_+mD_); tTD_.head(mT_) = tT_; tTD_.tail(mD_) = tD_;
        tTD_ = VecSortUniq(tTD_); long mTD_ = tTD_.size();
        VectorXi Tline_ = gettline(tT_,tTD_); VectorXi Dline_= gettline(tD_,tTD_);
        VectorXi DTline_= gettline(tD_,tT_);
        MatrixXd tcumuinc(mTD_,2); tcumuinc.setZero(); tcumuinc.col(0) = tTD_;
        
        double wb, wr, bi, ri, thetaT, thetaD;
        for (int kb=0;kb<ngrid_;++kb){
            wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
            for (int kr=0;kr<ngrid_;++kr){
                wr=Ggrid_(kr,0)*pow(M_PI,-0.5); ri=Ggrid_(kr,1)*pow(2.0*sigma2_(1),0.5);
                thetaT = exp(GX.transpose() * beta_ + bi + ri);
                thetaD = exp(GX.transpose() * gamma_ + rho_(0) * bi + rho_(1) * ri);
                for (int j=0;j<mTD_;++j){
                    tcumuinc(j,1) += wb * wr * (1-exp(-A_(Dline_(j))*thetaD)) * exp(-Lambda_(Tline_(j))*thetaT);
                    for (int jj=0;jj<=Tline_(j);++jj)
                        tcumuinc(j,1) += wb * wr * lambda_(jj) * thetaT * exp(-Lambda_(jj)*thetaT) * (1-exp(-A_(DTline_(jj))*thetaD));
                }
            }
        }
        return(tcumuinc);
    }
    MatrixXd surv_T(VectorXd GX){
        MatrixXd tcumuinc(mT_,2); tcumuinc.setZero(); tcumuinc.col(0) = tT_;
        
        double wb, wr, bi, ri, thetaT;
        for (int kb=0;kb<ngrid_;++kb){
            wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
            for (int kr=0;kr<ngrid_;++kr){
                wr=Ggrid_(kr,0)*pow(M_PI,-0.5); ri=Ggrid_(kr,1)*pow(2.0*sigma2_(1),0.5);
                thetaT = exp(GX.transpose() * beta_ + bi + ri);
                for (int j=0;j<mT_;++j){
                    tcumuinc(j,1) += wb * wr * exp(-Lambda_(j)*thetaT);
                }
            }
        }
        return(tcumuinc);
    }
    MatrixXd surv_D(VectorXd GX){
        MatrixXd tcumuinc(mD_,2); tcumuinc.setZero(); tcumuinc.col(0) = tD_;
        
        double wb, wr, bi, ri, thetaD;
        for (int kb=0;kb<ngrid_;++kb){
            wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
            for (int kr=0;kr<ngrid_;++kr){
                wr=Ggrid_(kr,0)*pow(M_PI,-0.5); ri=Ggrid_(kr,1)*pow(2.0*sigma2_(1),0.5);
                thetaD = exp(GX.transpose() * gamma_ + rho_(0) * bi + rho_(1) * ri);
                for (int j=0;j<mD_;++j){
                    tcumuinc(j,1) += wb * wr * exp(-A_(j)*thetaD);
                }
            }
        }
        return(tcumuinc);
    }
    
    MatrixXd surv_T_given_Fam(Gene_Fam Fam_info, VectorXd L, VectorXd R, VectorXi Rinf, VectorXd U, VectorXi epsilon, MatrixXd X){
        long nsub = Fam_info.nsub_; long nclust = Fam_info.nclust_;
        vector<MatrixXd> rmat (nsub); vector<vector<MatrixXd>> gmat(nsub);
        VectorXd pr(nsub); vector<VectorXd> pg(nsub); VectorXd riFM(2);
        for (int ij=0;ij<nsub;++ij){
            int i = Fam_info.cluster_(ij);
            rmat[ij] = Fam_info.Gene_calr(ij); pr(ij) = rmat[ij].rows();
            gmat[ij].resize(Fam_info.ng_par_(i)); pg[ij].resize(Fam_info.ng_par_(i));
            for (int kFM=0;kFM <Fam_info.ng_par_(i);++kFM){
                gmat[ij][kFM] = Fam_info.Gene_calg(ij, kFM); pg[ij](kFM) = gmat[ij][kFM].rows();
            }
        }
        
        MatrixXd tcumuinc(mT_,nclust); tcumuinc.setZero();
        VectorXi Uline = gettline(tT_,U);
        VectorXi Lline = gettline(tD_,L);
        VectorXi Rline = gettline(tD_,R);
        double denom;
        double wb, wF, wM, bi, riF, riM, rij, gij, thetaT, thetaD, lij, lij_rg, lik_FM; VectorXd survT(mT_);
        for (int i=0;i<nclust;++i){
            denom = 0;
            for (int kb=0;kb<ngrid_;++kb){
                wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
                for (int kF=0;kF<ngrid_;++kF){
                    wF=Ggrid_(kF,0)*pow(M_PI,-0.5); riF=Ggrid_(kF,1)*pow(2.0*sigma2_(1),0.5);
                    for (int kM=0;kM<ngrid_;++kM){
                        wM=Ggrid_(kM,0)*pow(M_PI,-0.5); riM=Ggrid_(kM,1)*pow(2.0*sigma2_(1),0.5);
                        riFM(0) = riF; riFM(1) = riM;
                        for (int kFM=0;kFM <Fam_info.ng_par_(i);++kFM){
                            lik_FM = 1; survT.setZero();
                            for (int j=0;j<Fam_info.ni_[i];++j){
                                int ij = Fam_info.cluster_id_[i](j);
                                lij = 0;
                                for (int r=0;r<pr(ij);++r){
                                    for (int g=0;g<pg[ij](kFM);++g){
                                        rij = riFM((int)rmat[ij](r,1)); gij = gmat[ij][kFM](g,1);
                                        thetaT = exp(beta_(0) * (double) gij + X.row(ij) * beta_.tail(p_) + bi + rij);
                                        thetaD = exp(gamma_(0) * (double) gij + X.row(ij) * gamma_.tail(p_) + rho_(0) * bi + rho_(1) * rij);
                                        lij_rg = exp(-A_(Lline(ij))*thetaD);
                                        if (Rinf(ij)==0) lij_rg -= exp(-A_(Rline(ij))*thetaD);
                                        
                                        if (j==0) {
                                            for (int l=0;l<mT_;++l){
                                                survT(l) += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg * exp(-Lambda_(l)*thetaT);
                                            }
                                        }
                                        
                                        if (epsilon(ij)==2) lij_rg *= (1-exp(-Lambda_(Uline(ij))*thetaT));
                                        else lij_rg *= exp(-Lambda_(Uline(ij))*thetaT);
                                        if (epsilon(ij)==1) lij_rg *= lambda_(Uline(ij))*thetaT;
                                        
                                        lij += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg;
                                    }
                                }
                                if (j==0) { if (lij >0 ) survT /= lij;}
                                lik_FM *= lij;
                            }
                            tcumuinc.col(i) += wb * wF * wM * Fam_info.P_par_[i](kFM) * survT * lik_FM;
                            denom += wb * wF * wM * Fam_info.P_par_[i](kFM) * lik_FM;
                        }
                    }
                }
            }
            tcumuinc.col(i) /= denom;
        }
        return(tcumuinc);
    }
    
    VectorXd surv_T_given_Fam_t(Gene_Fam Fam_info, VectorXd L, VectorXd R, VectorXi Rinf, VectorXd U, VectorXi epsilon, MatrixXd X, double t){
        long nsub = Fam_info.nsub_; long nclust = Fam_info.nclust_;
        vector<MatrixXd> rmat (nsub); vector<vector<MatrixXd>> gmat(nsub);
        VectorXd pr(nsub); vector<VectorXd> pg(nsub); VectorXd riFM(2);
        for (int ij=0;ij<nsub;++ij){
            int i = Fam_info.cluster_(ij);
            rmat[ij] = Fam_info.Gene_calr(ij); pr(ij) = rmat[ij].rows();
            gmat[ij].resize(Fam_info.ng_par_(i)); pg[ij].resize(Fam_info.ng_par_(i));
            for (int kFM=0;kFM <Fam_info.ng_par_(i);++kFM){
                gmat[ij][kFM] = Fam_info.Gene_calg(ij, kFM); pg[ij](kFM) = gmat[ij][kFM].rows();
            }
        }
        VectorXd tcumuinc(nclust); tcumuinc.setZero();
        cout<<"nsub="<<nsub<<endl;
        cout<<"nclust="<<nclust<<endl;
        VectorXi Uline = gettline(tT_,U);
        VectorXi Lline = gettline(tD_,L);
        VectorXi Rline = gettline(tD_,R);
        VectorXd tv(1); tv<<t;
        VectorXi tline = gettline(tT_,tv);
        double denom;
        double wb, wF, wM, bi, riF, riM, rij, gij, thetaT, thetaD, lij, lij_rg, lik_FM; double survT = 0;
        for (int i=0;i<nclust;++i){
            denom = 0;
            for (int kb=0;kb<ngrid_;++kb){
                wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
                for (int kF=0;kF<ngrid_;++kF){
                    wF=Ggrid_(kF,0)*pow(M_PI,-0.5); riF=Ggrid_(kF,1)*pow(2.0*sigma2_(1),0.5);
                    for (int kM=0;kM<ngrid_;++kM){
                        wM=Ggrid_(kM,0)*pow(M_PI,-0.5); riM=Ggrid_(kM,1)*pow(2.0*sigma2_(1),0.5);
                        riFM(0) = riF; riFM(1) = riM;
                        for (int kFM=0;kFM <Fam_info.ng_par_(i);++kFM){
                            lik_FM = 1; survT = 0;
                            for (int j=0;j<Fam_info.ni_[i];++j){
                                int ij = Fam_info.cluster_id_[i](j);
                                lij = 0;
                                for (int r=0;r<pr(ij);++r){
                                    for (int g=0;g<pg[ij](kFM);++g){
                                        rij = riFM((int)rmat[ij](r,1)); gij = gmat[ij][kFM](g,1);
                                        thetaT = exp(beta_(0) * (double) gij + X.row(ij) * beta_.tail(p_) + bi + rij);
                                        thetaD = exp(gamma_(0) * (double) gij + X.row(ij) * gamma_.tail(p_) + rho_(0) * bi + rho_(1) * rij);
                                        lij_rg = exp(-A_(Lline(ij))*thetaD);
                                        if (Rinf(ij)==0) lij_rg -= exp(-A_(Rline(ij))*thetaD);
                                        
                                        if (j==0) {
                                            survT += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg * exp(-Lambda_(tline(0))*thetaT);
                                        }
                                        
                                        if (epsilon(ij)==2) lij_rg *= (1-exp(-Lambda_(Uline(ij))*thetaT));
                                        else lij_rg *= exp(-Lambda_(Uline(ij))*thetaT);
                                        if (epsilon(ij)==1) lij_rg *= lambda_(Uline(ij))*thetaT;
                                        
                                        lij += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg;
                                    }
                                }
                                if (j==0) { if (lij >0 ) survT /= lij;}
                                lik_FM *= lij;
                            }
                            tcumuinc(i) += wb * wF * wM * Fam_info.P_par_[i](kFM) * survT * lik_FM;
                            denom += wb * wF * wM * Fam_info.P_par_[i](kFM) * lik_FM;
                        }
                    }
                }
            }
            tcumuinc(i) /= denom;
        }
        return(tcumuinc);
    }
    VectorXd surv_T_given_Fam_t(BiCens_Fam data, double t){
        long nsub = data.Fam_info_.nsub_; long nclust = data.Fam_info_.nclust_;
        vector<MatrixXd> rmat (nsub); vector<vector<MatrixXd>> gmat(nsub);
        VectorXd pr(nsub); vector<VectorXd> pg(nsub); VectorXd riFM(2);
        for (int ij=0;ij<nsub;++ij){
            int i = data.Fam_info_.cluster_(ij);
            rmat[ij] = data.Fam_info_.Gene_calr(ij); pr(ij) = rmat[ij].rows();
            gmat[ij].resize(data.Fam_info_.ng_par_(i)); pg[ij].resize(data.Fam_info_.ng_par_(i));
            for (int kFM=0;kFM <data.Fam_info_.ng_par_(i);++kFM){
                gmat[ij][kFM] = data.Fam_info_.Gene_calg(ij, kFM); pg[ij](kFM) = gmat[ij][kFM].rows();
            }
        }
        VectorXd tcumuinc(nclust); tcumuinc.setZero();
        cout<<"nsub="<<nsub<<endl;
        cout<<"nclust="<<nclust<<endl;
        VectorXi Uline = gettline(tT_,data.U_);
        VectorXi Lline = gettline(tD_,data.L_);
        VectorXi Rline = gettline(tD_,data.R_);
        VectorXd tv(1); tv<<t;
        VectorXi tline = gettline(tT_,tv);
        double denom;
        double wb, wF, wM, bi, riF, riM, rij, gij, thetaT, thetaD, lij, lij_rg, lik_FM; double survT = 0;
        for (int i=0;i<nclust;++i){
            denom = 0;
            for (int kb=0;kb<ngrid_;++kb){
                wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2_(0),0.5);
                for (int kF=0;kF<ngrid_;++kF){
                    wF=Ggrid_(kF,0)*pow(M_PI,-0.5); riF=Ggrid_(kF,1)*pow(2.0*sigma2_(1),0.5);
                    for (int kM=0;kM<ngrid_;++kM){
                        wM=Ggrid_(kM,0)*pow(M_PI,-0.5); riM=Ggrid_(kM,1)*pow(2.0*sigma2_(1),0.5);
                        riFM(0) = riF; riFM(1) = riM;
                        for (int kFM=0;kFM <data.Fam_info_.ng_par_(i);++kFM){
                            lik_FM = 1; survT = 0;
                            for (int j=0;j<data.Fam_info_.ni_[i];++j){
                                int ij = data.Fam_info_.cluster_id_[i](j);
                                lij = 0;
                                for (int r=0;r<pr(ij);++r){
                                    for (int g=0;g<pg[ij](kFM);++g){
                                        rij = riFM((int)rmat[ij](r,1)); gij = gmat[ij][kFM](g,1);
                                        thetaT = exp(beta_(0) * (double) gij + data.X_.row(ij) * beta_.tail(p_) + bi + rij);
                                        thetaD = exp(gamma_(0) * (double) gij + data.X_.row(ij) * gamma_.tail(p_) + rho_(0) * bi + rho_(1) * rij);
                                        lij_rg = exp(-A_(Lline(ij))*thetaD);
                                        if (data.Rinf_(ij)==0) lij_rg -= exp(-A_(Rline(ij))*thetaD);
                                        
                                        if (j==0) {
                                            survT += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg * exp(-Lambda_(tline(0))*thetaT);
                                        }
                                        
                                        if (data.epsilon_(ij)==2) lij_rg *= (1-exp(-Lambda_(Uline(ij))*thetaT));
                                        else lij_rg *= exp(-Lambda_(Uline(ij))*thetaT);
                                        if (data.epsilon_(ij)==1) lij_rg *= lambda_(Uline(ij))*thetaT;
                                        
                                        lij += rmat[ij](r,0) * gmat[ij][kFM](g,0) * lij_rg;
                                    }
                                }
                                if (j==0) { if (lij >0 ) survT /= lij;}
                                lik_FM *= lij;
                            }
                            tcumuinc(i) += wb * wF * wM * data.Fam_info_.P_par_[i](kFM) * survT * lik_FM;
                            denom += wb * wF * wM * data.Fam_info_.P_par_[i](kFM) * lik_FM;
                        }
                    }
                }
            }
            tcumuinc(i) /= denom;
        }
        return(tcumuinc);
    }
    void est_sd_score(double h){
        h *= pow((double)nclust_,-0.5); Profile_I(h);
        long P = 2*(p_+1)+4;
        Cov_score=pI_score.inverse(); sd_score.resize(P);
        for (int j=0;j<P;++j) sd_score(j)=sqrt(Cov_score(j,j));
    }
    
    void est_sd_score_sig(double h){
        h *= pow((double)nclust_,-0.5); Profile_I_sig(h);
        long P = 2*(p_+1)+4;
        Cov_score=pI_score.inverse(); sd_score.resize(P);
        for (int j=0;j<P;++j) sd_score(j)=sqrt(Cov_score(j,j));
        for (int j=0;j<2;++j) sd_score(2*(p_+1)+2+j) *= sigma2_(j);
    }
};

BiCens_Fam simu_data(long nclust, VectorXd beta, VectorXd gamma, VectorXd rho, VectorXd sigma2, double rho_pop, long ngrid, default_random_engine & generator){
    //Distributions
    bernoulli_distribution X1dist(0.5);
    normal_distribution<double> X2dist(0.0,1.0);
    bernoulli_distribution Bdist(0.5);
    normal_distribution<double> bdist(0.0,sqrt(sigma2(0)));
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
        b = bdist(generator); riF = rdist(generator); riM = rdist(generator);
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

class BiCens_Fam_pred {
private:
    long nsub_;
public:
    VectorXd L_, R_, U_; VectorXi Rinf_, epsilon_; MatrixXd X_; Gene_Fam Fam_info_; long ngrid_;
    VectorXd T_; VectorXd PT_;
    VectorXi Ind_event_t1_; VectorXi Ind_event_t2_;
    BiCens_Fam_pred() = default;
    BiCens_Fam_pred(VectorXd L, VectorXd R, VectorXi Rinf, VectorXd U, VectorXi epsilon, MatrixXd X, Gene_Fam Fam_info, long ngrid, VectorXd T, VectorXd PT, VectorXi Ind_event_t1, VectorXi Ind_event_t2){
        L_ = L; R_ = R; Rinf_ = Rinf; U_ = U; epsilon_ = epsilon; X_ = X; Fam_info_ = Fam_info; ngrid_ = ngrid;
        T_ = T; PT_ = PT; Ind_event_t1_ = Ind_event_t1; Ind_event_t2_ = Ind_event_t2;
        nsub_ = T_.size();
    }
    
    BiCens_Fam BiCens_Fam_data(){
        BiCens_Fam data(L_, R_, Rinf_, U_, epsilon_, X_, Fam_info_, ngrid_);
        return(data);
    }
    double C_index(VectorXd PT_est){
        double numer = 0; double denom = 0;
        for (int i=0;i<nsub_;++i){
            if ((Ind_event_t1_(i)==0)&(Ind_event_t2_(i)==1)){
                for (int j=0;j<nsub_;++j){
                    if ((Ind_event_t1_(j)==0)&(Ind_event_t2_(j)==0)){
                        denom ++;
                        if (PT_est(i)<PT_est(j)) numer += 1.0;
                        if (PT_est(i)==PT_est(j)) numer += 0.5;
                    }
                }
            }
        }
        numer /= denom;
        return(numer);
    }
    
    double Ratio(VectorXd PT_est){
        double numer = 0; double denom = 0;
        for (int i=0;i<nsub_;++i){
            if (Ind_event_t1_(i)==0){
                numer = numer + (double)Ind_event_t2_(i);
                denom = denom + (1.0 - PT_est(i));
            }
        }
        numer /= denom;
        return(numer);
    }
    double MSEP(VectorXd PT_est){
        double MSEP = 0; double denom = 0;
        for (int i=0;i<nsub_;++i){
            if (Ind_event_t1_(i)==0){
                MSEP += pow(PT_est(i)-PT_(i),2);
                denom ++;
            }
        }
        MSEP /= denom;
        return(MSEP);
    }
};
BiCens_Fam_pred simu_data_pred(long nclust, VectorXd beta, VectorXd gamma, VectorXd rho, VectorXd sigma2, double rho_pop, long ngrid, default_random_engine & generator, double t1, double t2){
    // All information up to time t1 to predict survival at time t2
    //Distributions
    bernoulli_distribution X1dist(0.5);
    normal_distribution<double> X2dist(0.0,1.0);
    bernoulli_distribution Bdist(0.5);
    normal_distribution<double> bdist(0.0,sqrt(sigma2(0)));
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
        b = bdist(generator); riF = rdist(generator); riM = rdist(generator);
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

#endif /* BiCens_Fam_h */
