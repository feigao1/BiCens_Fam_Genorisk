# BiCens_Fam_Genorisk

This is the C++ code for "Semiparametric Regression Analysis of Bivariate Censored Events in a Family Study of Alzheimer’s Disease".

BiCens_Fam.h includes the main computation codes for the proposed bivariate regression analysis approach.
BiCens_Fam_T.h includes the computation codes for univariate regression analysis.
BiCens_Fam_mismodel.h includes the codes for generating data from mis-specified models.
Geno_surv_base.h includes supporting codes for computation.

data.cpp gives computation code for analyzing an data example given in simu_data_int.dat, with regression analysis results given in Data.res.pdf.

shiny_geno provides a shiny app that reproduce an interactive web application that predicts Alzheimer’s Disease risk.
