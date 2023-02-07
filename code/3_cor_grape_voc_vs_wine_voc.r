#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:\\Users\\Bigghost\\Documents\\GitHub\\smoke_taint_prediction"
setwd(WORKING_DIR)

#Create the output directory for storing the output figures
dir.create("output")
dir.create("output/3_cor_grape_voc_vs_wine_voc")

#Clean the working space and import necessary libraries
rm(list=ls())


#Library for plotting
library(RColorBrewer)
library(circlize)           #(Color palette generator)
library(ComplexHeatmap)
library(ggplot2)

#Functions for loading data
source("code/functions_data_loading.r")



#Loading data
data = read.csv("data/dataset_for_analysis.csv",header=T,check.names = F,encoding="UTF-8",stringsAsFactors = F)
data$`Sensory_ash (smoke taint out of 100)` = as.numeric(data$`Sensory_ash (smoke taint out of 100)`)
grape_data = data[data$`Sample type` == "Grapes",]
wine_data = data[data$`Sample type` != "Grapes",]

#Constants of the data set (e.g., column indices)
colidx_sample_name = 1
colidx_voc = 20 : (ncol(data) - 1)
colidx_voc_with_binary_indicators = 9 : (ncol(data) - 1)
colidx_smoke_taint_idx = ncol(data)

#Get median and mean for grape and wine data respectively, 
#then merged them together for analysis

integrated_grape_wine_data = get_integrated_grape_wine_data(get_median_grape_data(grape_data),
                                                            get_mean_wine_data(wine_data))


#Cor (between grape feature and wine feature)
cor_matrix = cor(integrated_grape_wine_data[,1:20],integrated_grape_wine_data[,21:40])
cor_matrix = t(cor_matrix)
colnames(cor_matrix) = gsub("(grape)","",colnames(cor_matrix),fixed=T)
colnames(cor_matrix) = gsub("_"," ",colnames(cor_matrix),fixed=T)
rownames(cor_matrix) = gsub("_"," ",rownames(cor_matrix),fixed=T)
svg("output/3_cor_grape_voc_vs_wine_voc/cor_grape_wine_feature.svg",width=10,height=8)
h = Heatmap(cor_matrix,
            column_title = "(Grape VOCs)",
            row_title = "(Wine VOCs)",
            col = colorRamp2(c(-0.5, 0.5, 1.5), c("blue", "white", "red")),
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = T,
            show_row_names = T,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(formatC(cor_matrix[i,j], digits=2, format="f"), x, y)
            }, 
            column_title_side = "bottom",
            show_heatmap_legend = FALSE)
draw(h)
dev.off()

svg("output/3_cor_grape_voc_vs_wine_voc/grape_vs_wine_free_phenol.svg",width=5,height=5)
df = data.frame(x=integrated_grape_wine_data[,"(grape)free_phenol"],y=integrated_grape_wine_data[,"free_phenol"])
h = ggplot(df, aes(x=x,y=y)) +
    geom_point(size=3) + 
    theme_classic() + 
    ggtitle(paste0("Free Phenol Concentration\n(Grape VOCs VS Wine VOCs)(PCC = ",round(cor(df$x,df$y)*100)/100,")")) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Wine VOC concetration (ug/L)",x="Grape VOC concentration (ug/Kg)")
print(h)
dev.off()

svg("output/3_cor_grape_voc_vs_wine_voc/grape_vs_wine_free_guaiacol_o-cresol.svg",width=5,height=5)
df = data.frame(x=integrated_grape_wine_data[,"(grape)free_guaiacol"],y=integrated_grape_wine_data[,"free_o-cresol"])
h = ggplot(df, aes(x=x,y=y)) +
    geom_point(size=3) + 
    theme_classic() + 
    ggtitle(paste0("Concentration\n(Grape VOCs VS Wine VOCs)(PCC = ",round(cor(df$x,df$y)*100)/100,")")) + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12)) +
    labs(y="Free o-cresol (Wine) concetration (ug/L)",x="Free guaiacol (Grape) concentration (ug/Kg)")
print(h)
dev.off()