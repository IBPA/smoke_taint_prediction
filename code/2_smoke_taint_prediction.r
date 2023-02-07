#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:\\Users\\Bigghost\\Documents\\GitHub\\smoke_taint_prediction"
setwd(WORKING_DIR)

#Please set the number of trials for cross-validation
N_Trial = 2

#Create the output directory for storing the output figures
dir.create("output")
dir.create("output/2_smoke_taint_prediction")

#Clean the working space and import necessary libraries
rm(list=ls())

#Library for prediction model
library(glmnet)             #(Lasso)
library(e1071)              #(SVM)


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


#Functions for model prediction
source("code/functions_smoke_taint_prediction.r")

#Functions for plotting

plot_cv_result = function(filename, width, height, cv_results, method_name){
    svg(filename, width, height)
    
    df = data.frame(x = cv_results$truth_Y,
                    y = cv_results$pred_Y,
                    fold = as.factor(cv_results$fold_id))
    names(df) = c("x","y","fold")
    point_size = 3
    
    pcc = cor(df$x, df$y)
    
    h = ggplot(df, aes(x = x, y = y, color = fold)) + 
        geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle(paste0("Smoke Taint Index Prediction\n(",method_name,", ",max(cv_results$fold_id),"-fold CV)")) + 
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        labs(y="Predicted Smoke Taint (out of 100)",x="Smoke Taint (out of 100)",color="Fold")+
        scale_color_brewer(palette="Set1", direction=-1) + 
        scale_y_continuous(breaks=seq(0,110,20),limits=c(0,110)) + 
        scale_x_continuous(breaks=seq(0,110,20),limits=c(0,110)) +  
        annotate(geom="text", x=0, y=100, label=paste0("N = ",length(cv_results$truth_Y)," (Wines)\nPCC = ",format(pcc,digits=2)), hjust=0,
                 color="black") +  
        theme(legend.position = c(0.9, 0.3))
    print(h)
    
    dev.off()
}


plot_lasso_loading = function(filename, lasso_coef, xmin,xmax,xstep){
    lasso_coef = lasso_coef[,1]
    lasso_coef = lasso_coef[lasso_coef != 0]
    lasso_coef = lasso_coef[order(abs(lasso_coef),decreasing = F)]
    
    y = names(lasso_coef)
    y = sapply(y, function(x){
        gsub("_"," ",x)
    })
    y = sapply(y, function(x){
        gsub(")",") ",x)
    })
    y = factor(y,levels = y)
    
    x = lasso_coef
    pos = as.factor(x >= 0)
    
    df = data.frame(x = x, y = y, pos = pos)
    
    svg(filename,width=6,height=5)
    p = ggplot(df, aes(x=x, y=y, fill=pos)) +
        ggtitle("Coefficients of Lasso Model") + 
        geom_col() + 
        theme_classic() + 
        scale_x_continuous(breaks=seq(xmin,xmax,xstep),limits=c(xmin,xmax)) +   
        geom_vline(xintercept = 0) + 
        xlab(paste0("Coefficients")) + 
        theme(legend.position = "none",
              plot.margin =  ggplot2::margin(7, 7, 7, 7, "pt"),
              plot.title = element_text(hjust = 0.5, size=14),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=12),
              axis.title.y = element_blank())
    
    print(p)
    dev.off()
}

#Linear regression
features = integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)]
scores = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
linear_model = lm(scores~features)

lm_predict = predict(linear_model)
cor(lm_predict,scores)

#Lasso regression
lambdas = 10^seq(3,-3,-0.0001)


lasso_result_list = list()
lasso_result_grape_list = list()
svr_result_list = list()
svr_result_grape_list = list()
rf_result_list = list()
rf_result_grape_list = list()

log_lasso_result_list = list()
log_lasso_result_grape_list = list()
log_svr_result_list = list()
log_svr_result_grape_list = list()
log_rf_result_list = list()
log_rf_result_grape_list = list()

for (i in 1 : N_Trial){
    print(i)
    lasso_result_list[[i]] = run_lasso(matrix(scores,ncol=1),features) 
    lasso_result_grape_list[[i]] = run_lasso(matrix(scores,ncol=1),features[,1:20]) 
    svr_result_list[[i]] = run_svr(matrix(scores,ncol=1),features)
    svr_result_grape_list[[i]] = run_svr(matrix(scores,ncol=1),features[,1:20])
    rf_result_list[[i]] = run_rf(matrix(scores,ncol=1),features)
    rf_result_grape_list[[i]] = run_rf(matrix(scores,ncol=1),features[,1:20])
    
    log_lasso_result_list[[i]] = run_lasso(matrix(scores,ncol=1),log(features+1)) 
    log_lasso_result_grape_list[[i]] = run_lasso(matrix(scores,ncol=1),log(1+features[,1:20])) 
    log_svr_result_list[[i]] = run_svr(matrix(scores,ncol=1),log(1+features))
    log_svr_result_grape_list[[i]] = run_svr(matrix(scores,ncol=1),log(1+features[,1:20]))
    log_rf_result_list[[i]] = run_rf(matrix(scores,ncol=1),log(1+features))
    log_rf_result_grape_list[[i]] = run_rf(matrix(scores,ncol=1),log(1+features[,1:20]))
}

pcc_lasso = sapply(lasso_result_list,function(x){x$lasso_performance[1]})
pcc_lasso_grape = sapply(lasso_result_grape_list,function(x){x$lasso_performance[1]})
pcc_svr = sapply(svr_result_list,function(x){x$svr_performance[1]})
pcc_svr_grape = sapply(svr_result_grape_list,function(x){x$svr_performance[1]})
pcc_rf = sapply(rf_result_list,function(x){x$rf_performance[1]})
pcc_rf_grape = sapply(rf_result_grape_list,function(x){x$rf_performance[1]})

pcc_log_lasso = sapply(log_lasso_result_list,function(x){x$lasso_performance[1]})
pcc_log_lasso_grape = sapply(log_lasso_result_grape_list,function(x){x$lasso_performance[1]})
pcc_log_svr = sapply(log_svr_result_list,function(x){x$svr_performance[1]})
pcc_log_svr_grape = sapply(log_svr_result_grape_list,function(x){x$svr_performance[1]})
pcc_log_rf = sapply(log_rf_result_list,function(x){x$rf_performance[1]})
pcc_log_rf_grape = sapply(log_rf_result_grape_list,function(x){x$rf_performance[1]})

get_spearman = function(input){
    raw_cor_val = sapply(input,function(x){cor(x$pred_Y,x$truth_Y,method="spearman")})
    return(raw_cor_val)
}


results_pcc = list(pcc_lasso = pcc_lasso,
                   pcc_lasso_grape = pcc_lasso_grape,
                   pcc_svr = pcc_svr,
                   pcc_svr_grape = pcc_svr_grape,
                   pcc_rf = pcc_rf,
                   pcc_rf_grape = pcc_rf_grape,
                   pcc_log_lasso = pcc_log_lasso,
                   pcc_log_lasso_grape = pcc_log_lasso_grape,
                   pcc_log_svr = pcc_log_svr,
                   pcc_log_svr_grape = pcc_log_svr_grape,
                   pcc_log_rf = pcc_log_rf,
                   pcc_log_rf_grape = pcc_log_rf_grape)

heatmap_matrix = matrix(0,nrow=2,ncol=6)
heatmap_matrix[1,1] = mean(results_pcc$pcc_lasso)
heatmap_matrix[1,2] = mean(results_pcc$pcc_svr)
heatmap_matrix[1,3] = mean(results_pcc$pcc_rf)

heatmap_matrix[1,4] = mean(results_pcc$pcc_log_lasso)
heatmap_matrix[1,5] = mean(results_pcc$pcc_log_svr)
heatmap_matrix[1,6] = mean(results_pcc$pcc_log_rf)

heatmap_matrix[2,1] = mean(results_pcc$pcc_lasso_grape)
heatmap_matrix[2,2] = mean(results_pcc$pcc_svr_grape)
heatmap_matrix[2,3] = mean(results_pcc$pcc_rf_grape)

heatmap_matrix[2,4] = mean(results_pcc$pcc_log_lasso_grape)
heatmap_matrix[2,5] = mean(results_pcc$pcc_log_svr_grape)
heatmap_matrix[2,6] = mean(results_pcc$pcc_log_rf_grape)

method_name = c("Lasso","SVR","RF","Lasso","SVR","RF")
feature_selection_method = c("(Use All Features)","(Grape Features Only)")
svg("output/2_smoke_taint_prediction/model_performance.svg",width=8,height=1.5)
h = Heatmap(heatmap_matrix,
        #column_title = "Model Performance (PCC) (5-Fold CV, 10 times)",
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
        column_split = c(rep("(Conc)",3),rep("(Log(Conc+1))",3)),
        row_split = c(1,2),
        show_row_dend = F,
        show_column_dend = F,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_title_gp = gpar(fontsize=12),
        column_names_gp = gpar(fontsize=10),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(formatC(heatmap_matrix[i,j], digits=2, format="f"), x, y)
        }, 
        top_annotation = HeatmapAnnotation(
            text = anno_text(method_name, rot = 0, just = "center"),
            annotation_height = max_text_height(method_name)
        ), 
        left_annotation = rowAnnotation(
            text = anno_text(feature_selection_method, rot = 0, just = "left"),
            annotation_width = max_text_width(feature_selection_method)
        ),
        row_title = NULL,
        column_title_side = "bottom",
        show_heatmap_legend = FALSE)
draw(h)
dev.off()


spearman_lasso = get_spearman(lasso_result_list)
spearman_lasso_grape = get_spearman(lasso_result_grape_list)
spearman_svr = get_spearman(svr_result_list)
spearman_svr_grape = get_spearman(svr_result_grape_list)
spearman_rf = get_spearman(rf_result_list)
spearman_rf_grape = get_spearman(rf_result_grape_list)

spearman_log_lasso = get_spearman(log_lasso_result_list)
spearman_log_lasso_grape = get_spearman(log_lasso_result_grape_list)
spearman_log_svr = get_spearman(log_svr_result_list)
spearman_log_svr_grape = get_spearman(log_svr_result_grape_list)
spearman_log_rf = get_spearman(log_rf_result_list)
spearman_log_rf_grape = get_spearman(log_rf_result_grape_list)

results_spearman = list(spearman_lasso = spearman_lasso,
                   spearman_lasso_grape = spearman_lasso_grape,
                   spearman_svr = spearman_svr,
                   spearman_svr_grape = spearman_svr_grape,
                   spearman_rf = spearman_rf,
                   spearman_rf_grape = spearman_rf_grape,
                   spearman_log_lasso = spearman_log_lasso,
                   spearman_log_lasso_grape = spearman_log_lasso_grape,
                   spearman_log_svr = spearman_log_svr,
                   spearman_log_svr_grape = spearman_log_svr_grape,
                   spearman_log_rf = spearman_log_rf,
                   spearman_log_rf_grape = spearman_log_rf_grape)

heatmap_matrix = matrix(0,nrow=2,ncol=6)
heatmap_matrix[1,1] = mean(results_spearman$spearman_lasso)
heatmap_matrix[1,2] = mean(results_spearman$spearman_svr)
heatmap_matrix[1,3] = mean(results_spearman$spearman_rf)

heatmap_matrix[1,4] = mean(results_spearman$spearman_log_lasso)
heatmap_matrix[1,5] = mean(results_spearman$spearman_log_svr)
heatmap_matrix[1,6] = mean(results_spearman$spearman_log_rf)

heatmap_matrix[2,1] = mean(results_spearman$spearman_lasso_grape)
heatmap_matrix[2,2] = mean(results_spearman$spearman_svr_grape)
heatmap_matrix[2,3] = mean(results_spearman$spearman_rf_grape)

heatmap_matrix[2,4] = mean(results_spearman$spearman_log_lasso_grape)
heatmap_matrix[2,5] = mean(results_spearman$spearman_log_svr_grape)
heatmap_matrix[2,6] = mean(results_spearman$spearman_log_rf_grape)

method_name = c("Lasso","SVR","RF","Lasso","SVR","RF")
feature_selection_method = c("(Use All Features)","(Grape Features Only)")
svg("output/2_smoke_taint_prediction/model_performance_spearman.svg",width=8,height=1.5)
h = Heatmap(heatmap_matrix,
            #column_title = "Model Performance (Spearman) (5-Fold CV, 10 times)",
            col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
            column_split = c(rep("(Conc)",3),rep("(Log(Conc+1))",3)),
            row_split = c(1,2),
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = F,
            show_row_names = F,
            column_title_gp = gpar(fontsize=12),
            column_names_gp = gpar(fontsize=10),
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(formatC(heatmap_matrix[i,j], digits=2, format="f"), x, y)
            }, 
            top_annotation = HeatmapAnnotation(
                text = anno_text(method_name, rot = 0, just = "center"),
                annotation_height = max_text_height(method_name)
            ), 
            left_annotation = rowAnnotation(
                text = anno_text(feature_selection_method, rot = 0, just = "left"),
                annotation_width = max_text_width(feature_selection_method)
            ),
            row_title = NULL,
            column_title_side = "bottom",
            show_heatmap_legend = FALSE)
draw(h)
dev.off()




plot_cv_result("output/2_smoke_taint_prediction/lasso_regression_results.svg",width=5,height=5,lasso_result_list[[1]],"Lasso")
plot_cv_result("output/2_smoke_taint_prediction/lasso_regression_results_grape_feature_only.svg",width=5,height=5,lasso_result_grape_list[[1]],"Lasso (Grape Features Only)")
plot_cv_result("output/2_smoke_taint_prediction/svr_regression_results.svg",width=5,height=5,svr_result_list[[1]],"SVR (Radial)")
plot_cv_result("output/2_smoke_taint_prediction/svr_regression_results_grape_feature_only.svg",width=5,height=5,svr_result_grape_list[[1]],"SVR (Radial) (Grape Features Only)")
plot_cv_result("output/2_smoke_taint_prediction/rf_regression_results.svg",width=5,height=5,rf_result_list[[1]],"RF")
plot_cv_result("output/2_smoke_taint_prediction/rf_regression_results_grape_feature_only.svg",width=5,height=5,rf_result_grape_list[[1]],"RF (Grape Features Only)")

plot_cv_result("output/2_smoke_taint_prediction/log_lasso_regression_results.svg",width=5,height=5,log_lasso_result_list[[1]],"Lasso (Log)")
plot_cv_result("output/2_smoke_taint_prediction/log_lasso_regression_results_grape_feature_only.svg",width=5,height=5,log_lasso_result_grape_list[[1]],"Lasso (Grape Features Only) (Log)")
plot_cv_result("output/2_smoke_taint_prediction/log_svr_regression_results.svg",width=5,height=5,log_svr_result_list[[1]],"SVR (Radial) (Log)")
plot_cv_result("output/2_smoke_taint_prediction/log_svr_regression_results_grape_feature_only.svg",width=5,height=5,log_svr_result_grape_list[[1]],"SVR (Radial) (Grape Features Only) (Log)")
plot_cv_result("output/2_smoke_taint_prediction/log_rf_regression_results.svg",width=5,height=5,log_rf_result_list[[1]],"RF (Log)")
plot_cv_result("output/2_smoke_taint_prediction/log_rf_regression_results_grape_feature_only.svg",width=5,height=5,log_rf_result_grape_list[[1]],"RF (Grape Features Only) (Log)")


plot_lasso_loading("output/2_smoke_taint_prediction/lasso_coef_log.svg",log_lasso_result_list[[1]]$coef_lasso,0,12,2)
plot_lasso_loading("output/2_smoke_taint_prediction/lasso_coef_log_grape.svg",log_lasso_result_grape_list[[1]]$coef_lasso,0,12,2)

save.image("output/2_smoke_taint_prediction/ModelResult.rdata")

