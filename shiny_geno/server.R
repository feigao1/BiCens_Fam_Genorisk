library('statmod')
beta <- c(0.57554500,-0.20268200,-0.53772400,-0.02268650,-0.09567591/22)
gamma <- c(9.769920e-02,2.824580e-01,2.368890e-01,1.507100e-01,2.194409e-05/22)
rho <- c(-0.0464307,1.3846400)
sigma2 <- c(0.415472,0.458008)

t_AD <- c(0, 20, 46, 54, 55, 56, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105)

lambda_AD <- c(0.000173147, 0.00563039, 0.00463082, 0.008618, 0.00235976, 0.000134312, 0.000977178, 0.00135529, 0.000730329, 0.000204033, 0.000389595, 0.00929352, 1.65026e-08, 1.34271e-09, 0.00170903, 0.0055136, 0.00559892, 0.0018253, 0.0107833, 0.00772264, 0.00938111, 0.0195799, 0.0177279, 0.00879089, 0.0251789, 0.0333902, 0.0303299, 0.0396658, 0.0358045, 0.0392793, 0.0544901, 0.0558864, 0.0832322, 0.0514982, 0.0919443, 0.0770457, 0.0936562, 0.10851, 0.167261, 0.17401, 0.259495, 0.262975, 0.214607, 0.195381, 0.329084, 2.03363e-110, 0.25799, 0.339133, 9.2343e-121, 0.232609, 1.50676e-103)
t_CVD <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 109, 117, 119)

lambda_CVD <- c(0, 0.00050169, 0.000721833, 0.000860749, 0.000936183, 0.000954973, 0.00103157, 0.00110504, 0.00127093, 0.00135764, 0.00144769, 0.0004267, 0.00045714, 0.000501036, 0.000517444, 0.000550457, 0.000584753, 0.000703253, 0.000863106, 0.000963502, 0.00122042, 0.00156245, 0.00196113, 0.00139517, 0.00168276, 0.00115851, 0.000832545, 0.000966369, 0.000710788, 0.000797961, 0.000830301, 0.00114975, 0.00142503, 0.00161749, 0.00202366, 0.00167683, 0.0024152, 0.00214733, 0.00188905, 0.00260636, 0.00275323, 0.00318986, 0.00319438, 0.0029329, 0.00343072, 0.00381603, 0.00305634, 0.0032381, 0.00416669, 0.00299803, 0.00206538, 0.00100278, 0.00112332, 0.00096225, 0.000486085, 0.000446202, 0.000150933, 0.000150123, 0.000121572, 6.3855e-05, 4.11982e-05, 1.27488e-05, 1.24925e-05, 5.76312e-06, 3.90906e-06, 5.80142e-05, 0.000556405, 0.00418934, 0.0054613, 0.0108044, 0.0142426, 0.0131268, 0.0159906, 0.00439285, 0.027532, 0.0195162, 0.0220237, 0.0221983, 0.0129938, 0.0211283, 0.0273072, 0.0279725, 0.0326987, 0.0209407, 0.0429271, 0.0224766, 0.0450978, 0.045922, 0.0418027, 0.0563925, 6.60433e-07, 0.0981743, 0.017689, 0.031742, 0.046841, 0.00313794, 9.40855e-05, 0.0673834, 0.0421236, 0.0350673, 0.0602591, 2.34501e-32, 1.64182e-178, 2.24346e-142, 8.9407e-127, 8.42438e-162, 0, 0, 0, 0)

ngrid=20
grid = gauss.quad(ngrid,kind="hermite")
pred_CVD <- function(GX, age, prob_CVD_int, prob_CVD_delta = 0){
	t = t_AD; mT = length(t); survT = rep(0,mT); denom = 0
	Tline = sum(t_AD<=age); Lline = sum(t_CVD<=prob_CVD_int[1]); Rline = sum(t_CVD<=prob_CVD_int[2])
	AD_cum = cumsum(lambda_AD); CVD_cum = cumsum(lambda_CVD)
	for (kb in 1:ngrid){
		wb = grid$weights[kb]/sqrt(pi); bi = grid$nodes[kb]*sqrt(2*sigma2[1])
		for (kr in 1:ngrid){
			wr = grid$weights[kr]/sqrt(pi); ri = grid$nodes[kr]*sqrt(2*sigma2[2])
			thetaT = exp(sum(GX*beta)+ bi + ri); thetaD = exp(sum(GX*gamma) + rho[1] * bi + rho[2] * ri)
			H = exp(-CVD_cum[Lline]*thetaD)
			if (is.finite(prob_CVD_int[2])) H = H - exp(-CVD_cum[Rline]*thetaD)
			if (prob_CVD_delta!=0) H = H * lambda_CVD[Lline]*thetaD
			for (j in 1:mT) survT[j] = survT[j] + wb * wr * exp(-AD_cum[j]*thetaT) * H
			denom = denom + wb * wr * exp(-AD_cum[Tline]*thetaT) * H;
		}
	}
	survT = survT / denom
	result = data.frame(t=t,pred=1-survT); result = result[result$t>=age,]
	return(result)
}
# GX = c(0,0,0,0,10); prob_CVD_int = c(65,Inf); age = 65; tail(pred_CVD(GX,age,prob_CVD_int))
# GX = c(1,0,0,0,10); prob_CVD_int = c(65,Inf); age = 65; tail(pred_CVD(GX,age,prob_CVD_int))

gene_prob <- function(G,pop){ return((G==2)*(pop^2)+(G==1)*(2*pop*(1-pop))+(G==0)*((1-pop)^2))}
child_given_par <- function(GF, GM, GC){
	if (GF+GM==4) return(1*(GC==2)) # AA and AA --> AA
    else if (GF+GM==0) return(1*(GC==0)) # aa and aa --> aa
    else if (GF+GM==3) return(0.5*(GC>=1)) # AA and Aa --> AA or Aa
    else if (GF+GM==1) return(0.5*(GC<=1)) # Aa and aa --> Aa or aa
    else{
    		if (GF*GM==0) return(1*(GC==1)) # AA and aa --> Aa
        else return(0.25*(GC==1)+0.25) # Aa and Aa --> AA, Aa, or AA
    }
}
 
pred_relative <- function(GX,age,prob_CVD_int, prob_CVD_delta = 0,rel_info,pop){
	# rel_info should be a dataframe: relation, AD_L, AD_R, AD_delta, CVD_L, CVD_R, CVD_delta, G, X1-X4 (9-12)
	rel_info = rel_info[order(rel_info$relation),] # in the order of parent, sibling
	par_info = rel_info[rel_info$relation==1,]; npar = dim(par_info)[1]
	sib_info = rel_info[rel_info$relation==2,]; nsib = dim(sib_info)[1]
	nsub = npar + nsib + 1
	######## G combo and corresponding p ########
	nsub_g = 2 + nsib + 1 # two parents
	Gene_comb = NULL; Gene = seq(3)-1
	for (j in 1:nsub_g){
		Mat = kronecker(kronecker(rep(1,3^(j-1)),Gene),rep(1,3^(nsub_g-j)))
		Gene_comb = cbind(Gene_comb,Mat)
	}
	Gene_comb = Gene_comb[Gene_comb[,1]==GX[1],]
	if (npar>0){for (j in 1:npar){ if (!is.na(par_info$G[j])){ Gene_comb = Gene_comb[Gene_comb[,(j+1)]==par_info$G[j],]}}}
	if (nsib>0){for (j in 1:nsib){ if (!is.na(sib_info$G[j])){ Gene_comb = Gene_comb[Gene_comb[,(j+3)]==sib_info$G[j],]}}}
	Gene_comb = as.matrix(Gene_comb)
	ng = dim(Gene_comb)[1]; pg = rep(0,ng)
	for (k in 1:ng){
		GF = Gene_comb[k,2]; GM = Gene_comb[k,3]
		pg[k] = gene_prob(GF,pop) * gene_prob(GM,pop) * child_given_par(GF, GM, Gene_comb[k,1])
		if (nsib>0){for (j in 1:nsib) pg[k] = pg[k] * child_given_par(GF, GM, Gene_comb[k,(j+3)])}
	}
	pg = pg / sum(pg)
	######## r combo and corresponding p ########
	nsub_r = nsib + 1 # only vary children
	r_comb = NULL; r_seq = seq(2)
	for (j in 1:nsub_r){
		Mat = kronecker(kronecker(rep(1,2^(j-1)),r_seq),rep(1,2^(nsub_r-j)))
		r_comb = cbind(r_comb,Mat)
	}
	if (nsib>0) r_comb = cbind(r_comb[,1],rep(1,2^nsub_r),rep(2,2^nsub_r),r_comb[,2:nsub_r]) else r_comb = cbind(r_comb[,1],rep(1,2^nsub_r),rep(2,2^nsub_r))
	nr = dim(r_comb)[1]; pr = rep(1/nr,nr)
	AD_L = AD_R = AD_delta = CVD_L = CVD_R = CVD_delta = rep(0,nsub)
	AD_Lline = AD_Rline = CVD_Lline = CVD_Rline = rep(0,nsub)
	for (j in 1:nsub){
		if (j==1){
			AD_L[j] = age; AD_R[j] = Inf; AD_delta[j] = 0
			CVD_L[j] = prob_CVD_int[1]; CVD_R[j] = prob_CVD_int[2]
			CVD_delta[j] = prob_CVD_delta
		} else {
			AD_L[j] = rel_info$AD_L[j-1]; AD_R[j] = rel_info$AD_R[j-1]
			AD_delta[j] = rel_info$AD_delta[j-1]
			CVD_L[j] = rel_info$CVD_L[j-1]; CVD_R[j] = rel_info$CVD_R[j-1]
			CVD_delta[j] = rel_info$CVD_delta[j-1]
		}
		AD_Lline[j] = sum(t_AD<=AD_L[j]); AD_Rline[j] = sum(t_AD<=AD_R[j])
		CVD_Lline[j] = sum(t_CVD<=CVD_L[j]); CVD_Rline[j] = sum(t_CVD<=CVD_R[j])
	}
	AD_cum = cumsum(lambda_AD); CVD_cum = cumsum(lambda_CVD)
	mT = length(AD_cum); lik = 0; survT = rep(0,mT)
	for (kb in 1:ngrid){
		wb = grid$weights[kb]/sqrt(pi); bi = grid$nodes[kb]*sqrt(2*sigma2[1])
		for (kF in 1:ngrid){
			wF = grid$weights[kF]/sqrt(pi); riF = grid$nodes[kF]*sqrt(2*sigma2[2])
			for (kM in 1:ngrid){
				wM = grid$weights[kM]/sqrt(pi); riM = grid$nodes[kM]*sqrt(2*sigma2[2])
				riFM = c(riF,riM)
				for (r in 1:nr){
					for (g in 1:ng){
						survT_rg = rep(0,mT); lik_rg = 1;
						for (j in 1:nsub){
							if (j==1) {Xj = as.numeric(GX[2:5])} else {
								if (j<=npar+1) {Xj = as.numeric(par_info[(j-1),9:12])
								} else Xj = as.numeric(sib_info[(j-npar-1),9:12])}
							rij = riFM[r_comb[r,j]]; gij = Gene_comb[g,j]; GXj = c(gij,Xj)
							thetaT = exp(sum(beta*GXj) + bi + rij)
							thetaD = exp(sum(gamma*GXj)+ rho[1] * bi + rho[2] * rij)
							
							lij = exp(-CVD_cum[CVD_Lline[j]]*thetaD)
							if (is.finite(CVD_R[j])) lij = lij - exp(-CVD_cum[CVD_Rline[j]]*thetaD)
							if (CVD_delta[j]==1) lij = lij * lambda_CVD[CVD_Lline[j]]*thetaD
							lik_rg = lik_rg * lij
							lij = exp(-AD_cum[AD_Lline[j]]*thetaT)
							if (is.finite(AD_R[j]))  lij = lij - exp(-AD_cum[AD_Rline[j]]*thetaT)
							if (AD_delta[j]==1) lij = lij * lambda_AD[AD_Lline[j]]*thetaT
							lik_rg = lik_rg * lij
							if (j==1) { 
								for (l in AD_Lline[j]:mT) survT_rg[l] = exp(-AD_cum[l]*thetaT + AD_cum[AD_Lline[j]]*thetaT)
							}
						}
						lik = lik + wb * wF * wM * pr[r] * pg[g] * lik_rg;
						survT = survT + wb * wF * wM * pr[r] * pg[g] * lik_rg * survT_rg;
						if (is.na(survT[1])) {print(c(kb,kF,kM,r,g));break}
						# kb=15;kF=1;kM=20;r=2;g=1
					}
				}
			}
		}
	}
	survT = survT / lik; t = t_AD
	result = data.frame(t=t,pred=1-survT); result = result[result$t>=age,]
	return(result)
}

# GX = c(0,0,0,0,10); prob_CVD_int = c(65,Inf); prob_CVD_delta = 0; age = 65
# relation=1;AD_L = 0; AD_R=65; AD_delta=0; CVD_L=0; CVD_R=65; CVD_delta=0; G = NA;
 # X1=0; X2 = 0; X3=0; X4=10
# rel_info = data.frame(relation, AD_L, AD_R, AD_delta, CVD_L, CVD_R, CVD_delta, G, X1,X2,X3,X4)
# res = pred_relative(GX,age,prob_CVD_int, prob_CVD_delta = 0,rel_info,pop=0.15)


determine_time <- function(dis, dis_ind, dis_date, dis_int, cur_age){
	int = rep(0,2)
	if (dis =='yes'){
		if (dis_ind == 'yes'){
			delta = 1; int[1] = dis_date; int[2] = Inf
		} else { delta = 0; int=dis_int}
	} else {
		if (dis=='no'){ delta = 0; int = c(cur_age,Inf)}
		else {delta = 0; int = c(0,Inf)}
	}
	return(list(int = int, delta = delta))
}
inputtores <-function(input){
	GX = c(1*(input$APOEe4=='1')+2*(input$APOEe4=='2'), 1*(input$sex=='male'),1*(input$race=='white'),1*(input$race=='black'),input$educ)
	age = input$age
	prob_CVD_res = determine_time(input$CVD, input$CVD_ind, input$CVD_date, input$CVD_int, age)
	prob_CVD_delta = prob_CVD_res$delta
	prob_CVD_int = prob_CVD_res$int

	# if (input$CVD =='yes'){
		# if (input$CVD_ind == 'yes'){
			# prob_CVD_delta = 1; prob_CVD_int[1] = input$CVD_date; prob_CVD_int[2] = Inf
		# } else { prob_CVD_delta = 0; prob_CVD_int=input$CVD_int}
	# } else {
		# if (input$CVD=='no'){ prob_CVD_delta = 0; prob_CVD_int = c(age,Inf)}
		# else {prob_CVD_delta = 0; prob_CVD_int = c(0,Inf)}
	# }
	if (input$relinfo != 'yes'){ # no relative information
		return(pred_CVD(GX, age, prob_CVD_int, prob_CVD_delta))
	} else {
		# first relative
		if (input$Rel1 == 'parent') relation = 1 else relation = 2
		if (input$Rel1_APOEe4=='unknown') G=NA else {
			if (input$Rel1_APOEe4=='0') G=0 else {
				if (input$Rel1_APOEe4=='1') G=1 else G=2
			}
		}
		sex = 1*(input$Rel1_sex=='male')
		race_w = 1*(input$Rel1_race=='white'); race_b = 1*(input$Rel1_race=='black')
		educ = input$Rel1_educ
		AD_res = determine_time(input$Rel1_AD, input$Rel1_AD_ind, input$Rel1_AD_date, input$Rel1_AD_int, input$Rel1_AD_cens)
		AD_delta = AD_res$delta; AD_L = AD_res$int[1]; AD_R = AD_res$int[2]
		CVD_res = determine_time(input$Rel1_CVD, input$Rel1_CVD_ind, input$Rel1_CVD_date, input$Rel1_CVD_int, input$Rel1_CVD_cens)
		CVD_delta = CVD_res$delta; CVD_L = CVD_res$int[1]; CVD_R = CVD_res$int[2]
		if (input$Rel2_w=='yes'){  # two relatives
			if (input$Rel2 == 'parent') relation = c(relation,1) else relation = relation = c(relation,2)
			if (input$Rel2_APOEe4=='unknown') G_p=NA else {
				if (input$Rel2_APOEe4=='0') G_p=0 else {
					if (input$Rel2_APOEe4=='1') G_p=1 else G_p=2
				}
			}
			G = c(G,G_p); sex = c(sex, 1*(input$Rel2_sex=='male'))
			race_w = c(race_w, 1*(input$Rel2_race=='white')); race_b = c(race_b, 1*(input$Rel2_race=='black'))
			educ = c(educ, input$Rel2_educ)
			AD_res = determine_time(input$Rel2_AD, input$Rel2_AD_ind, input$Rel2_AD_date, input$Rel2_AD_int, input$Rel2_AD_cens)
			AD_delta = c(AD_delta, AD_res$delta); AD_L = c(AD_L, AD_res$int[1]); AD_R = c(AD_R, AD_res$int[2])
			CVD_res = determine_time(input$Rel2_CVD, input$Rel2_CVD_ind, input$Rel2_CVD_date, input$Rel2_CVD_int, input$Rel2_CVD_cens)
			CVD_delta = c(CVD_delta, CVD_res$delta); CVD_L = c(CVD_L, CVD_res$int[1]); CVD_R = c(CVD_R, CVD_res$int[2])
		}
		rel_info = data.frame(relation, AD_L, AD_R, AD_delta, CVD_L, CVD_R, CVD_delta, G, sex, race_w, race_b, educ)
		return(pred_relative(GX,age,prob_CVD_int, prob_CVD_delta,rel_info,pop=0.15))
	}
}

# Define server
server <- function(input, output) {
	dataInput <- eventReactive(input$go, {inputtores(input)})
	output$plotgraph <- renderPlot({
		t = dataInput()$t; pred = dataInput()$pred
		pred = pred[t<=90]; t = t[t<=90]
		par(family = 'sans')
		plot(t,pred,type = 'l', lwd=1.5, cex.main=1.5, cex.lab=1.5,ylim=c(0,1),xlab = 'Age',ylab = '‘Cumulative incidence of Alzheimer’s disease',main = 'Prediction of Future Disease Risk')
	})
	
	output$plottable <- renderTable({
		t_choice = 5*(11 + seq(7)) ### 60 - 90
		Age = round(dataInput()$t[dataInput()$t %in% t_choice])
		Cumulative.Incidence = dataInput()$pred[dataInput()$t %in% t_choice]
		table_data = data.frame(Age,Cumulative.Incidence)
		table_data
	})
}
