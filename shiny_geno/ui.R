ui <- fluidPage(
	# App title ----
	titlePanel("Risk Prediction for Alzheimer's Disease"),
	
	# Sidebar layout with input and output definitions ----
	sidebarLayout(
		# Sidebar panel for inputs ----
		sidebarPanel(
			# Input ----
  			selectInput('sex', 'Gender', c("male","female")),
  			selectInput('race', 'Race', c("white", "black", "other")),
  			sliderInput('age', 'Current age',  min = 0, max = 90, value = 50),
  			selectInput("APOEe4", "APOE epsilon4 allelle frequency",  c(0,1,2)),
  			sliderInput("educ", "Years of education",  min = 0, max = 25, value = 12),
  			selectInput('CVD', 'Have you developed any cardiovascular disease?', c("unknown","yes", "no")),
  			conditionalPanel(condition = "input.CVD == 'yes'", selectInput("CVD_ind", "Do you know the age of your first cardivascular disease onset?", c("","yes","no"))),
  			conditionalPanel(condition = "input.CVD_ind == 'yes'", sliderInput("CVD_date", "Your age of first cardivascular disease onset", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.CVD_ind == 'no'", sliderInput("CVD_int", "Your best estimate of an interval that brackets your age of first cardivascular disease onset", min = 0, max = 90, value = c(40,60))),
  			selectInput('relinfo', 'Do you have information on Alzheimer\'s disease or cardivascular disease history of your parents and/or siblings?', c("","yes","no")),
  			
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1", "Please input disease history information: what is his/her relation to you?", c("","parent","sibling"))),
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1_sex", "His/her gender", c("female", "male"))),
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1_race", "His/her race", c("black", "white", "other"))),
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1_APOEe4", "His/her APOE epsilon4 allelle frequency",  c('unknown','0','1','2'))),
  			conditionalPanel(condition = "input.relinfo == 'yes'", sliderInput("Rel1_educ", "His/her years of education",  min = 0, max = 25, value = 10)),
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1_AD", "Has he/she developed Alzheimer\'s disease?", c( "unknown","yes", "no"))),
  			conditionalPanel(condition = "input.Rel1_AD == 'yes'", selectInput("Rel1_AD_ind", "Do you know his/her age of Alzheimer's disease onset?", c("","yes","no"))),
  			conditionalPanel(condition = "input.Rel1_AD_ind == 'yes'", sliderInput("Rel1_AD_date", "His/her age of Alzheimer's disease onset", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.Rel1_AD_ind == 'no'", sliderInput("Rel1_AD_int", "Your best estimate of the interval that brackets his/her age of Alzheimer's disease onset", min = 0, max = 90, value = c(40,60))),
  			conditionalPanel(condition = "input.Rel1_AD == 'no'", sliderInput("Rel1_AD_cens", "What is his/her current age?", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.relinfo == 'yes'", selectInput("Rel1_CVD", "Has he/she developed any cardivascular disease?", c( "unknown","yes", "no"))),
  			conditionalPanel(condition = "input.Rel1_CVD == 'yes'", selectInput("Rel1_CVD_ind", "Do you know his/her age of first cardivascular disease onset?", c("","yes","no"))),
  			conditionalPanel(condition = "input.Rel1_CVD_ind == 'yes'", sliderInput("Rel1_CVD_date", "His/her age of first cardiovascular disease onset", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.Rel1_CVD_ind == 'no'", sliderInput("Rel1_CVD_int", "Your best estimate of the interval that brackets his/her age of first cardiovascular disease onset", min = 0, max = 90, value = c(40,60))),
  			conditionalPanel(condition = "input.Rel1_CVD == 'no'", sliderInput("Rel1_CVD_cens", "What is his/her current age?", min = 0, max = 90, value = 50)),
  			
  			conditionalPanel(condition = "input.Rel1 != ''", selectInput("Rel2_w", "Do you have Alzheimer\'s disease and/or cardivascular disease history of another parent and/or sibling?", c("","yes","no"))),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2", "Please input disease history information: what is his/her relation to you?", c("","parent","sibling"))),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2_sex", "His/her gender", c("female", "male"))),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2_race", "His/her race", c("black", "white", "other"))),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2_APOEe4", "His/her APOE epsilon4 allelle frequency",  c('unknown','0','1','2'))),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", sliderInput("Rel2_educ", "His/her years of education",  min = 0, max = 25, value = 10)),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2_AD", "Has he/she developed Alzheimer\'s disease?", c( "unknown","yes", "no"))),
  			conditionalPanel(condition = "input.Rel2_AD == 'yes'", selectInput("Rel2_AD_ind", "Do you know his/her age of Alzheimer's disease onset?", c("","yes","no"))),
  			conditionalPanel(condition = "input.Rel2_AD_ind == 'yes'", sliderInput("Rel2_AD_date", "His/her age of Alzheimer's disease onset", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.Rel2_AD_ind == 'no'", sliderInput("Rel2_AD_int", "Your best estimate of the interval that brackets his/her age of Alzheimer's disease onset", min = 0, max = 90, value = c(40,60))),
  			conditionalPanel(condition = "input.Rel2_AD == 'no'", sliderInput("Rel2_AD_cens", "What is his/her current age?", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.Rel2_w == 'yes'", selectInput("Rel2_CVD", "Has he/she developed any cardivascular disease?", c( "unknown","yes", "no"))),
  			conditionalPanel(condition = "input.Rel2_CVD == 'yes'", selectInput("Rel2_CVD_ind", "Do you know his/her age of first cardivascular disease onset?", c("","yes","no"))),
  			conditionalPanel(condition = "input.Rel2_CVD_ind == 'yes'", sliderInput("Rel2_CVD_date", "His/her age of first cardiovascular disease onset", min = 0, max = 90, value = 50)),
  			conditionalPanel(condition = "input.Rel2_CVD_ind == 'no'", sliderInput("Rel2_CVD_int", "Your best estimate of the interval that brackets his/her age of first cardiovascular disease onset", min = 0, max = 90, value = c(40,60))),
  			conditionalPanel(condition = "input.Rel2_CVD == 'no'", sliderInput("Rel2_CVD_cens", "What is his/her current age?", min = 0, max = 90, value = 50)),
  			
  			actionButton("go", h4("Predict"))
  		),
  		# Main panel for displaying outputs ----
  		mainPanel(
 			h4("Prediction of risk of Alzheimer's Disease given patient’s characteristics including genetic risk factors and medical history of CVD and his/her family members\' disease history."),
  			h5("Gao, F., Zeng, D., & Wang, Y. (2021). Semiparametric regression analysis of bivariate censored events in a family study of Alzheimer’s disease. Biostatistics.kxab014"),
  			fluidRow(
  				column(7, plotOutput('plotgraph')),#, height = 600, width = 600)),
  				column(5,offset = .5, tableOutput('plottable'))
  			),
  			h5("Reference is a 50-year-old white man with 12 years of education, no APOE epsilon4 allelle, no cardiovascular dissease history, and unknown family history of Alzheimer\'s disease or cardiovascular disease. Cum.Inc and Cum.Inc.ref and predicted cumulative incidence of the individual and the reference.")
  		)
  	)
)