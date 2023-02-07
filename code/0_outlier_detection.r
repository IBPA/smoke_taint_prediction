#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:\\Users\\Bigghost\\Documents\\GitHub\\smoke_taint_prediction"
setwd(WORKING_DIR)

#Create the output directory for storing the output figures
dir.create("output")
dir.create("output/0_outlier_detection")

#Clean the working space and import necessary libraries
rm(list=ls())

#Library for plotting
library(RColorBrewer)
library(circlize)           #(Color palette generator)
library(ComplexHeatmap)
library(ggplot2)

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

#Function sets
getOutlierInd = function(x){
    q1 = quantile(x,0.25)
    q3 = quantile(x,0.75)
    iqr = q3 - q1
    
    ind = rep(FALSE, length(x))
    ind[which(x < (q1 - 1.5*iqr))] = T
    ind[which(x > (q3 + 1.5*iqr))] = T
    return(ind)
}

get_relative_std_rep_grape = function(grape_data){
    #Generate relative std (among triplicates) of grape VOC data
    #Input: Grape data
    #Output: Relative std (std/mean) among triplicates for each sample and for each VOCs (rows: samples, col: VOCs)
    relative_std_rep_grape = apply(grape_data[,colidx_voc],2,function(x){
        sample_name = unique(grape_data[,colidx_sample_name])
        results = numeric(length(sample_name))
        for (i in 1 : length(sample_name)){
            select_idx = grape_data[,colidx_sample_name] == sample_name[i]
            if (sum(abs(x[select_idx]))==0){
                results[i] = 0
            }else{
                results[i] = sd(x[select_idx])/mean(x[select_idx])
            }
        }
        names(results)= sample_name
        return(results)
    })
    return(relative_std_rep_grape)
}

get_relative_absdiff_rep_wine = function(wine_data){
    #Generate relative abs diff (among duplicates) of wine VOC data
    #Input: Wine data
    #Output: Relative abs diff (absdiff/mean) among duplicates for each sample and for each VOCs (rows: samples, col: VOCs)
    relative_absdiff_dup_wine = apply(wine_data[,c(colidx_voc)],2,function(x){
        sample_name = unique(paste0(wine_data[,colidx_sample_name],"_",wine_data$`Sample type`))
        results = numeric(length(sample_name))
        for (i in 1 : length(sample_name)){
            select_idx = paste0(wine_data[,colidx_sample_name],"_",wine_data$`Sample type`) == sample_name[i]
            if (is.na(sum(abs(x[select_idx])))){
                results[i] = NA
            }else{
                if (length(which(select_idx==T)) != 2){
                    stop()
                }else{
                    if (sum(abs(x[select_idx]))==0){
                        results[i] = 0
                    }else{
                        results[i] = abs(diff(x[select_idx]))/mean(x[select_idx])
                    }
                }
            }
        }
        names(results)= sample_name
        return(results)
    })
    return(relative_absdiff_dup_wine)
}



plot_grape_data_relative_std = function(filename, relative_std_rep_grape){
    #Show the heatmap of the relative std (std/mean) of triplicates of grape data
    #Input: Relative std (std/mean) of triplicates of grape data (Row: samples, Col: VOCs)
    #Output: The figures showing the relative std of triplicates of grape data
    
    mat_letters = apply(relative_std_rep_grape, 2, getOutlierInd)
    mat_letters = apply(mat_letters, 2, function(x){
        tmp = x
        tmp[x == T] = "*"
        tmp[x == F] = ""
        return(tmp)
    })
    
    dir.create("output")
    svg(filename,width=6,height=6)
    H = Heatmap(t(relative_std_rep_grape), name = "Std/Mean", column_title = "Relative standard deviation (Grapes)",
                show_column_dend = F, show_row_dend = F,
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    grid.text(t(mat_letters)[i, j], x, y, gp = gpar(fontsize=10))
                },
                row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8),
                heatmap_legend_param = list(direction="horizontal",title_position="leftcenter",x=unit(-5,"npc"), just = c("left", "bottom")))
    draw(H, heatmap_legend_side = "bottom")
    dev.off()
}


plot_wine_data_relative_absdiff = function(filename, relative_absdiff_dup_wine){
    #Show the heatmap of the relative std (std/mean) of triplicates of grape data
    #Input: Relative std (std/mean) of triplicates of grape data (Row: samples, Col: VOCs)
    #Output: The figures showing the relative std of triplicates of grape data
    
    mat_letters = apply(relative_absdiff_dup_wine, 2, getOutlierInd)
    mat_letters = apply(mat_letters, 2, function(x){
        tmp = x
        tmp[x == T] = "*"
        tmp[x == F] = ""
        return(tmp)
    })
    svg(filename,width=6,height=6)
    H = Heatmap(t(relative_absdiff_dup_wine[,1:(ncol(relative_absdiff_dup_wine))]), name = "Abs_Diff/Mean", column_title = "Relative Absolute Difference (Wines)",
                show_column_dend = F, show_row_dend = F,
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    grid.text(t(mat_letters)[i, j], x, y, gp = gpar(fontsize=10))
                },
                row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8),
                heatmap_legend_param = list(direction="horizontal",title_position="leftcenter",x=unit(-5,"npc"), just = c("left", "bottom")))
    draw(H, heatmap_legend_side = "bottom")
    dev.off()
}


#Scripts
relative_std_rep_grape = get_relative_std_rep_grape(grape_data)
relative_absdiff_dup_wine = get_relative_absdiff_rep_wine(wine_data)

plot_grape_data_relative_std("output/0_outlier_detection/grape_relative_std_triplicates.svg", relative_std_rep_grape)
plot_wine_data_relative_absdiff("output/0_outlier_detection/wine_relative_absdiff_duplicates.svg", relative_absdiff_dup_wine)

