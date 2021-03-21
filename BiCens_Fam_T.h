#ifndef BiCens_Fam_T_h
#define BiCens_Fam_T_h
#include "Geno_surv_base.h"

class BiCens_Fam_T {
private:
    //Parameters
    VectorXd beta0_, lambda0_;
    Vector2d sigma20_;
    
    //Speicification
    long ngrid_; MatrixXd Ggrid_;
    const double eps=1E-3; const int maxiter_=10000;
    
    //Observation-related
    VectorXi Yline_; VectorXi Uline_;
    vector<MatrixXd> rmat; vector<vector<MatrixXd>> gmat; VectorXd riFM;
    VectorXd pr; vector<VectorXd> pg;
    
    //Variables
    VectorXd b_i, b2_i, r2_iF, r2_iM;
    VectorXd G_ij, r_ij;
    VectorXd thetaT_ij, GthetaT_ij, G2thetaT_ij;
    MatrixXd W_ijl, WG_ijl;
    
    VectorXd Score; MatrixXd Hessian;
    double sumtheta, sumthetab, sumthetabb, sumthetar, sumthetarr, sumthetabr;
    VectorXd sumthetax, sumthetabx, sumthetarx; MatrixXd sumthetaxx;
    MatrixXd pI_score;
    
    //Parameter Initiation
    void Para_init(){
        tT_ = get_distinct_time_m(U_, epsilon_); mT_ = tT_.size();
        lambda_.resize(mT_); lambda_.setOnes(); lambda_ *= 1/(double)mT_; Lambda_ = cumsum(lambda_);
        
        beta_.resize(p_+1); beta_.setZero(); sigma2_.setOnes();
    }
    void Other_init(){
        log_lik_.resize(nclust_);
        b_i.resize(nclust_); b2_i.resize(nclust_); r2_iF.resize(nclust_); r2_iM.resize(nclust_);
        G_ij.resize(nsub_); r_ij.resize(nsub_);
        thetaT_ij.resize(nsub_); GthetaT_ij.resize(nsub_); G2thetaT_ij.resize(nsub_);
        
        W_ijl.resize(nsub_,mT_); WG_ijl.resize(nsub_,mT_);
        
        sumthetax.resize(p_+1); sumthetabx.resize(p_+1); sumthetarx.resize(p_+1); sumthetaxx.resize(p_+1,p_+1);
        
        rmat.resize(nsub_); gmat.resize(nsub_); pr.resize(nsub_); pg.resize(nsub_); riFM.resize(2);
        for (int ij=0;ij<nsub_;++ij){
            int i = Fam_info_.cluster_(ij);
            rmat[ij] = Fam_info_.Gene_calr(ij); pr(ij) = rmat[ij].rows();
            gmat[ij].resize(Fam_info_.ng_par_(i)); pg[ij].resize(Fam_info_.ng_par_(i));
            for (int kFM=0;kFM <Fam_info_.ng_par_(i);++kFM){
                gmat[ij][kFM] = Fam_info_.Gene_calg(ij, kFM); pg[ij](kFM) = gmat[ij][kFM].rows();
            }
        }
        Uline_= gettline(tT_,U_);
    }
    double cal_thetaT(int ij, double bi, double rij, int gij, VectorXd beta){
        double thetaT = beta(0) * (double) gij + X_.row(ij) * beta.tail(p_) + bi + rij;
        return(exp(thetaT));
    }
    
    double cal_lij(int ij, double bi, double rij, int gij, double thetaT, VectorXd lambda, VectorXd Lambda){
        double lij = 1;
        if (epsilon_(ij)==2) lij *= (1-exp(-Lambda(Uline_(ij))*thetaT));
        else lij *= exp(-Lambda(Uline_(ij))*thetaT);
        if (epsilon_(ij)==1) lij *= lambda(Uline_(ij))*thetaT;
        return (lij);
    }
 
    void cal_likelihood(VectorXd beta, Vector2d sigma2, VectorXd lambda, VectorXd Lambda){
        b_i.setZero(); b2_i.setZero(); r2_iF.setZero(); r2_iM.setZero();
        G_ij.setZero(); r_ij.setZero();
        thetaT_ij.setZero(); GthetaT_ij.setZero(); G2thetaT_ij.setZero();
        W_ijl.setZero(); WG_ijl.setZero();
        log_lik_.setZero(); double lik_FM, lij;
        double wb,wF,wM; double bi,riF,riM; double rij,gij; double wlik_rg, wlik, wlik_j;
        double thetaT;
        for (int i=0;i<nclust_;++i){
            VectorXd eg(ni_(i)); VectorXd er(ni_(i));
            VectorXd thetat(ni_(i)); VectorXd gthetat(ni_(i)); VectorXd g2thetat(ni_(i));
            MatrixXd w(ni_(i),mT_); MatrixXd wg(ni_(i),mT_);
            for (int kb=0;kb<ngrid_;++kb){
                wb=Ggrid_(kb,0)*pow(M_PI,-0.5); bi=Ggrid_(kb,1)*pow(2.0*sigma2(0),0.5);
                if (ni_(i)>1){
                    for (int kF=0;kF<ngrid_;++kF){
                        wF=Ggrid_(kF,0)*pow(M_PI,-0.5); riF=Ggrid_(kF,1)*pow(2.0*sigma2(1),0.5);
                        for (int kM=0;kM<ngrid_;++kM){
                            wM=Ggrid_(kM,0)*pow(M_PI,-0.5); riM=Ggrid_(kM,1)*pow(2.0*sigma2(1),0.5);
                            riFM(0) = riF; riFM(1) = riM;
                            for (int kFM=0;kFM <Fam_info_.ng_par_(i);++kFM){
                                lik_FM = 1;
                                eg.setZero(); er.setZero();
                                thetat.setZero(); gthetat.setZero(); g2thetat.setZero();
                                w.setZero(); wg.setZero();
                                VectorXd lijm(ni_(i)); lijm.setOnes();
                                for (int j=0;j<ni_(i);++j){
                                    int ij = Fam_info_.cluster_id_[i](j);
                                    lij = 0;
                                    for (int r=0;r<pr(ij);++r){
                                        for (int g=0;g<pg[ij](kFM);++g){
                                            rij = riFM((int)rmat[ij](r,1)); gij = gmat[ij][kFM](g,1);
                                            thetaT = cal_thetaT(ij, bi, rij, gij, beta);
                                            wlik_rg = rmat[ij](r,0) * gmat[ij][kFM](g,0) * cal_lij(ij, bi, rij, gij, thetaT, lambda, Lambda);
                                            lij += wlik_rg;
                                            eg(j) += wlik_rg * gij; er(j) += wlik_rg * rij;
                                            thetat(j) += wlik_rg * thetaT;
                                            gthetat(j) += wlik_rg * gij * thetaT; g2thetat(j) += wlik_rg * pow(gij,2) * thetaT;
                                            if (epsilon_(ij)==1) {
                                                w(j,Uline_(ij)) += wlik_rg; wg(j,Uline_(ij)) += wlik_rg * gij;
                                            }
                                            else if (epsilon_(ij)==2){
                                                VectorXd ww (mT_); ww.setZero();
                                                for (int l=0;l<=Uline_(ij);++l) {
                                                    ww(l) = lambda(l) * thetaT / (1-exp(-Lambda(Uline_(ij))*thetaT));
                                                    if (1-exp(-Lambda(Uline_(ij))*thetaT)==0) ww(l) = 0;
                                                }
                                                w.row(j) += wlik_rg * ww;
                                                wg.row(j) += wlik_rg * ww * gij;
                                            }
                                        }
                                    }
                                    lik_FM *= lij;
                                    for (int jj=0;jj<ni_(i);++jj){if (jj!=j) lijm(jj) *= lij;}
                                }
                                wlik = wb * wF * wM * Fam_info_.P_par_[i](kFM) * lik_FM;
                                log_lik_(i) += wlik;
                                b_i(i) += wlik * bi; b2_i(i) += wlik * pow(bi,2);
                                r2_iF(i) += wlik * pow(riF,2); r2_iM(i) += wlik * pow(riM,2);
                            
                                for (int j=0;j<ni_(i);++j){
                                    wlik_j = wb * wF * wM * Fam_info_.P_par_[i](kFM) * lijm(j);
                                    int ij = Fam_info_.cluster_id_[i](j);
                                    G_ij(ij) += wlik_j * eg(j); r_ij(ij) += wlik_j * er(j);
                                    thetaT_ij(ij) += wlik_j * thetat(j);
                                    GthetaT_ij(ij) += wlik_j * gthetat(j); G2thetaT_ij(ij) += wlik_j * g2thetat(j);
                                    W_ijl.row(ij) += wlik_j * w.row(j);
                                    WG_ijl.row(ij) += wlik_j * wg.row(j);
                                }
                            }
                        }
                    }
                } else {
                    gij = (double) Fam_info_.Gene_pro_(i);
                    int ij = Fam_info_.cluster_id_[i](0);
                    for (int kF=0;kF<ngrid_;++kF){
                        wF = Ggrid_(kF,0)*pow(M_PI,-0.5); rij = Ggrid_(kF,1)*pow(2.0*sigma2(1),0.5);
                        thetaT = cal_thetaT(ij, bi, rij, gij, beta);
                        wlik = wb * wF * cal_lij(ij, bi, rij, gij, thetaT, lambda, Lambda);
                        
                        G_ij(ij) += wlik * gij; r_ij(ij) += wlik * rij;
                        thetaT_ij(ij) += wlik * thetaT;
                        GthetaT_ij(ij) += wlik * gij * thetaT; G2thetaT_ij(ij) += wlik * pow(gij,2) * thetaT;
                    
                        if (epsilon_(ij)==1) {
                            W_ijl(ij,Uline_(ij)) += wlik; WG_ijl(ij,Uline_(ij)) += wlik * gij;
                        }
                        else if (epsilon_(ij)==2){
                            VectorXd ww (mT_); ww.setZero();
                            for (int l=0;l<=Uline_(ij);++l) {
                                ww(l) = lambda(l) * thetaT / (1-exp(-Lambda(Uline_(ij))*thetaT));
                                if (1-exp(-Lambda(Uline_(ij))*thetaT)==0) ww(l) = 0;
                            }
                            W_ijl.row(ij) += wlik * ww; WG_ijl.row(ij) += wlik * ww * gij;
                        }
                        log_lik_(i) += wlik;
                        b_i(i) += wlik * bi; b2_i(i) += wlik * pow(bi,2);
                        r2_iF(i) += wlik * pow(rij,2); r2_iM(i) += wlik * pow(rij,2);
                    }
                }
            }
            b_i(i) /= log_lik_(i); b2_i(i) /= log_lik_(i);
            r2_iF(i)/= log_lik_(i); r2_iM(i)/= log_lik_(i);
            for (int j=0;j<ni_(i);++j){
                int ij = Fam_info_.cluster_id_[i](j);
                G_ij(ij) /= log_lik_(i); r_ij(ij) /= log_lik_(i);
                thetaT_ij(ij) /= log_lik_(i); GthetaT_ij(ij) /= log_lik_(i); G2thetaT_ij(ij) /= log_lik_(i);
                
                W_ijl.row(ij) /= log_lik_(i); WG_ijl.row(ij) /= log_lik_(i);
            }
            log_lik_(i) = log(log_lik_(i));
        }
    }
    
    VectorXd update_T(VectorXd& beta){
        Score.resize(p_+1); Hessian.resize(p_+1,p_+1); Score.setZero(); Hessian.setZero();
        VectorXd lambda(mT_);
        VectorXd Xij; double sumw; VectorXd sumwx(p_+1);
        for (int l=0;l<mT_;++l){
            sumtheta=0; sumthetax.setZero(); sumthetaxx.setZero();
            for (int ij=0;ij<nsub_;++ij){
                if (l<=Uline_(ij)){
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
            sumthetax /= sumtheta; sumthetaxx /= sumtheta;
            Score += sumwx - sumw * sumthetax;
            Hessian += sumw * (sumthetaxx - sumthetax*sumthetax.transpose());
            
            lambda(l) = W_ijl.col(l).sum() / sumtheta; if (sumtheta ==0) lambda(l)=0;
        }
        beta = beta + Hessian.inverse() * Score;
        return(lambda);
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
    
    Vector2d update_sigma(){
        Vector2d sigma2; sigma2 << b2_i.mean(), (r2_iF.mean() + r2_iM.mean())/2;
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
    
    VectorXd profilelogLik_i(VectorXd beta, Vector2d sigma2){
        iter_=0;
        cout<<"beta="<<beta.transpose()<<endl;
        cout<<"sigma2 = "<<sigma2.transpose()<<endl;
        VectorXd lambda=lambda_; VectorXd lambda0, Lambda;
        while (iter_<=maxiter_){
            lambda0 = lambda; Lambda = cumsum(lambda);
            cal_likelihood(beta, sigma2,lambda, Lambda);
            lambda = update_lambda();
            cout<<"iter="<<iter_<<", lambdadiff="<<diff(lambda,lambda0)<<endl;
            cout<<" logLik="<<log_lik_.mean()<<endl;
            iter_++;
            if (diff(lambda,lambda0)<eps) break;
        }
        return(log_lik_);
    }
    
    void update_h(int j, double h, VectorXd& beta, Vector2d& sigma2){
        if (j<p_+1) beta(j) += h;
        else sigma2(j-(p_+1)) *= exp(h);
    }
    
    void Profile_I(double h){
        VectorXd beta; Vector2d sigma2;
        long P = (p_+1) + 2;
        MatrixXd pL_i(nclust_,P);
        // calculate the value of pli(beta+h e)
        for (int j=0;j<P;++j){
            beta=beta_; sigma2 = sigma2_;
            update_h(j,h,beta, sigma2);
            pL_i.col(j)=profilelogLik_i(beta, sigma2);
        }
        pI_score.resize(P,P); pI_score.setZero(); MatrixXd pL_i_c=pL_i;
        for (int j=0;j<P;++j) pL_i_c.col(j) = (pL_i_c.col(j)-logLik_MLE_i)/h;
        for (int i=0;i<nclust_;++i) pI_score += pL_i_c.row(i).transpose()*pL_i_c.row(i);
    }
    
public:
    //Parameters
    Vector2d sigma2_; // sigmab2, sigmar2
    VectorXd beta_, tT_, lambda_, Lambda_;
    long mT_;
    
    //Observed Data
    VectorXd U_; VectorXi epsilon_; MatrixXd X_;
    Gene_Fam Fam_info_;
    long nsub_, nclust_, p_; VectorXi ni_;
    
    //Other output
    VectorXd log_lik_; int iter_;
    VectorXd logLik_MLE_i;
    MatrixXd Cov_score; VectorXd sd_score;
    
    BiCens_Fam_T() = default;
    BiCens_Fam_T(VectorXd U, VectorXi epsilon, MatrixXd X, Gene_Fam Fam_info, long ngrid){
        U_ = U; epsilon_ = epsilon; X_ = X; Fam_info_ = Fam_info;
        nsub_ = U_.size(); p_ = X_.row(0).size(); nclust_ = Fam_info_.nclust_; ni_ = Fam_info_.ni_;
        ngrid_ = ngrid; Ggrid_ = Hermite_grid (ngrid_);
        Para_init(); Other_init();
    }
    BiCens_Fam_T(VectorXd U, VectorXi epsilon, MatrixXd X, Gene_Fam Fam_info, long ngrid, VectorXd beta, VectorXd tT, VectorXd lambda, Vector2d sigma2){
        U_ = U; epsilon_ = epsilon; X_ = X; Fam_info_ = Fam_info;
        nsub_ = U_.size(); p_ = X_.row(0).size(); nclust_ = Fam_info_.nclust_; ni_ = Fam_info_.ni_;
        ngrid_ = ngrid; Ggrid_ = Hermite_grid (ngrid_);
        
        tT_ = tT; mT_ = tT_.size();
        lambda_ = lambda; Lambda_ = cumsum(lambda_);
        beta_ = beta; sigma2_ = sigma2;
        
        Other_init();
    }
    void solve(){
        iter_=0; double betadiff,lambdadiff,sigmadiff;
        while (iter_<=maxiter_){
            beta0_=beta_; lambda0_=lambda_;
            sigma20_ = sigma2_;
            Lambda_ = cumsum(lambda_);
            cal_likelihood(beta_, sigma2_, lambda_, Lambda_);
            lambda_ = update_T(beta_);
            sigma2_ = update_sigma();
            betadiff=diff(beta_,beta0_);
            lambdadiff=diff(lambda_,lambda0_);
            sigmadiff=diff(sigma20_,sigma2_);
            iter_++;
            if (betadiff+lambdadiff+sigmadiff<eps) break;
            cout<<"iter= "<<iter_<<", beta = "<<beta_.transpose()<<endl;
            cout<<"sigma2 = "<<sigma2_.transpose()<<endl;
            cout<<"lambda"<<lambda_.head(5).transpose()<<lambda_.tail(5).transpose();
            cout<<"betadiff="<<betadiff<<endl;
            cout<<"lambdadiff="<<lambdadiff<<", sigmadiff="<<sigmadiff<<endl;
            cout<<"logLik="<<log_lik_.sum()/(double)(nsub_)<<endl;
        }
        logLik_MLE_i=log_lik_;
        if (iter_<=maxiter_){
            cout<<"The Algorithm converges in "<<iter_<<" iterations."<<endl;
        }
        else {cout<<"Maximum Iteration Times ("<<maxiter_<<") Reached!"<<endl;}
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
        VectorXi Uline = gettline(tT_,data.U_);
        VectorXd tv(1); tv<<t;
        VectorXi tline = gettline(tT_,tv);
        double denom;
        double wb, wF, wM, bi, riF, riM, rij, gij, thetaT, lij, lij_rg, lik_FM; double survT = 0;
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
                                        lij_rg = 1;
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
        long P = (p_+1)+2;
        Cov_score=pI_score.inverse(); sd_score.resize(P);
        for (int j=0;j<P;++j) sd_score(j)=sqrt(Cov_score(j,j));
        for (int j=0;j<2;++j) sd_score(j+(p_+1)) *= sigma2_(j);
    }
};

BiCens_Fam_T simu_data(long nclust, double betaT, VectorXd gammaT, VectorXd sigma2, double rho, long ngrid, default_random_engine & generator){
    //Distributions
    bernoulli_distribution Xdist(0.5);
    bernoulli_distribution Bdist(0.5);
    normal_distribution<double> bdist(0.0,sqrt(sigma2(0)));
    normal_distribution<double> rdist(0.0,sqrt(sigma2(1)));
    uniform_real_distribution<double> udist(0.0,1.0);
    uniform_real_distribution<double> Cdist(0,6);
    
    //Patterns of Families
    //P, P+F/M, P+F+M, P+F/M+S1, P+F/M+S2, P+F+M+S1, P+F+M+S2
    VectorXi nip(7); nip << 1,2,3,3,4,4,5;
    VectorXd cump_fam (7); cump_fam << 0.4, 0.6, 0.8, 0.9, 0.9, 1, 1;
    VectorXi ni(nclust); VectorXi pattern(nclust); vector<VectorXi> relation_c(nclust);
    for (int i=0;i<nclust;++i){
        double u=udist(generator);
        for (int j=0;j<7;++j){
            if (u<=cump_fam(j)) { pattern(i) = j; ni(i) = nip(j); break;}
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
    MatrixXd X(nsub,1); VectorXd U(nsub); VectorXi epsilon(nsub);
    VectorXi relation(nsub); VectorXi cluster(nsub); VectorXi Gene(nsub); Gene.setZero();
    
    //For proband: 1/3, 1/3, 1/3
    MatrixXi gene_proband (3,1); gene_proband << 0,1,2;
    VectorXd p_proband (3); p_proband << 1.0/3.0, 1.0/3.0, 1.0/3.0;
    gene_p proband(gene_proband, p_proband);
    
    int line = 0; double b, riF, riM, rij, giF, giM, gij;
    double expb, T, C;
    for (int i=0;i<nclust;++i){
        cout<<ni(i)<<" subjects in cluster "<<i<<endl;
        b = bdist(generator); riF = rdist(generator); riM = rdist(generator);
        VectorXi giP = gen_gene (proband, generator);
        gene_p parent = parent_poss (giP(0),rho);
        VectorXi giFM = gen_gene (parent, generator);
        giF = giFM(0); giM = giFM(1);
        bool par1 = false;
        for (int j=0;j<ni(i);++j){
            X(line,0) = Xdist(generator);
            relation(line) = relation_c[i](j);
            cluster(line) = i;
            if (relation(line)==0){ // proband
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
            expb=exp(X.row(line) * gammaT + gij * betaT + b + rij);
            T = floor(genaT(0.5,expb,generator) * 100 + 0.5)/100;
            C = floor(Cdist(generator) * 100 + 0.5)/100;
            
            if (T<=C) {U(line) = T; epsilon(line) = 1;}
            else{ U(line) = C; epsilon(line) = 0;}
            line++;
        }
    }
    Gene_Fam Fam_info(relation, cluster, Gene, rho);
    BiCens_Fam_T data(U, epsilon, X, Fam_info, ngrid);
    return(data);
}

#endif /* BiCens_Fam_T_h */
