setwd('/Users/rita-gaofei/Desktop/2020_Genotype_survival/shiny_geno')
library(shiny)
source('ui.R')
source('server.R')

# Create Shiny app ----
shinyApp(ui = ui, server = server)