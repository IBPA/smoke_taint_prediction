#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:\\Users\\Bigghost\\Documents\\GitHub\\smoke_taint_prediction"
setwd(WORKING_DIR)

#Create the output directory for storing the output figures
dir.create("output")
dir.create("output/1_general_analysis")

#Clean the working space and import necessary libraries
rm(list=ls())

#Library for data analysis
library(preprocessCore)     #(Quantile Normalization for tsne)
library(tsne)
library(e1071)
library(glmnet)


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

#Function sets for visualization
plot_histogram_smoke_taint = function(filename, integrated_grape_wine_data){
    #Hist of wines
    svg(filename,width=6,height=6)
    dd<-data.frame(x=integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"])
    h = ggplot(dd, aes(x=x))+ geom_histogram(breaks=seq(0,70,10),fill="white",color="black",boundary=0) + theme_classic() +
        ggtitle(paste0("Distribution of Smoke Taint Index (N = ", nrow(integrated_grape_wine_data), ")")) + 
        #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        #scale_fill_manual(values=ColorSet)+
        labs(y="Count",x="Smoke Taint Index") + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 20))  + 
        scale_x_continuous(breaks=seq(0,70,10)) +
        stat_bin(breaks=seq(0,70,10), geom="text", color='black', aes(label=after_stat(count)),
                 position=position_stack(vjust = 0.5), boundary=0) 
    print(h)
    dev.off()
}


plot_pca_results = function(filename_pc1_pc2_smoke_taint,
                            filename_pc1_smoke_taint,
                            filename_abspc2_smoke_taint,
                            filename_pc1_loading,
                            filename_pc2_loading,
                            pca_result,
                            input_data,
                            xlim_loading1,
                            xlim_loading2){
    
    pc1 = pca_result$x[,1]
    pc2 = pca_result$x[,2]
    
    pc1_label = paste0("PC 1 (", format(100*pca_result$sdev[1]/sum(pca_result$sdev), digits=2), "% of variance explained)")
    pc2_label = paste0("PC 2 (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")
    abs_pc2_label = paste0("Abs(PC 2) (", format(100*pca_result$sdev[2]/sum(pca_result$sdev), digits=2), "% of variance explained)")

    score_tag = sapply(integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"], function(x){
        if (x <= 10) return("<= 10")
        if (x <= 20) return("10 < x <= 20")
        if (x <= 30) return("20 < x <= 30")
        if (x <= 40) return("30 < x <= 40")
        if (x <= 50) return("40 < x <= 50")
        if (x <= 60) return("50 < x <= 60")
        if (x <= 70) return("60 < x <= 70")
    })
    
    
    point_size = 3
    
    #SVM boundary (for PC1 VS PC2 plot)
    prepare_svm_grid_df = function(pc1, pc2){
        df = data.frame(pc1 = pc1, pc2 = pc2)
        high_low_index = rep("Low",nrow(df))
        high_low_index[integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"] > 30] = "High"
        df = cbind(df, high_low_index)
        svm_res = svm(high_low_index~pc1+pc2, data = df)
        
        grid_x = seq(-4, 18, length.out = 50)
        grid_y = seq(-7, 7, length.out = 50)
        grids = lapply(grid_x,function(x){
            cbind(rep(x,length(grid_y)),grid_y)
        })
        grids = do.call(rbind,grids)
        df_grids = data.frame(pc1 = grids[,1], pc2 = grids[,2])
        predict_grids = predict(svm_res, df_grids)
        predict_grids = as.character(predict_grids)
        predict_grids[predict_grids == "High"] = 1
        predict_grids[predict_grids == "Low"] = 0
        predict_grids = as.numeric(predict_grids)
        
        df_val = data.frame(x = grids[,1], y = grids[,2], zi = as.factor(predict_grids), z = predict_grids)
        
        return(list(df = df_val, support_vector_idx = svm_res$index))
    }
    
    svm_grid_df = prepare_svm_grid_df(pc1, pc2)
    
    df = data.frame(pc1 = pc1, pc2 = pc2, 
                    score_tag = score_tag, 
                    cont_score = integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"],
                    abs_pc2 = abs(pc2))
    
    
    #PC1 VS PC2
    svg(filename_pc1_pc2_smoke_taint,width=6,height=4)
    h = ggplot(df, aes(x = pc1, y = pc2, color = score_tag)) + geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle(paste0("PCA Plot (PC1 VS PC2) (N = ", nrow(integrated_grape_wine_data), ")")) + 
        #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        labs(y=pc2_label,x=pc1_label,color="Smoke Taint Index") +
        scale_color_brewer(palette="RdBu", direction=-1) + 
        stat_contour(data = svm_grid_df$df, aes(x=x,y=y,z=z), color="black", linetype=1,linewidth=1, breaks=c(0,1)) + 
        geom_point(data = df[svm_grid_df$support_vector_idx,], aes(x = pc1, y = pc2), color = "black", shape=4)
    print(h)
    dev.off()
    
    #PC1 VS Smoke Taint
    svg(filename_pc1_smoke_taint,width=4,height=4)
    h = ggplot(df, aes(x = pc1, y = cont_score)) + geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle(paste0("PC1 VS Smoke Taint Index (N = ", nrow(integrated_grape_wine_data), ")", "\n(PCC = ",format(cor(df$pc1,df$cont_score),digits=2) ,")")) + 
        #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        labs(y="Smoke Taint Index",x=pc1_label) +
        scale_color_brewer(palette="RdBu", direction=-1)
    print(h)
    dev.off()
    
    #ABS(PC2) VS Smoke Taint
    svg(filename_abspc2_smoke_taint,width=4,height=4)
    h = ggplot(df, aes(x = abs_pc2, y = cont_score)) + geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle(paste0("Abs(PC2) VS Smoke Taint Index (N = ", nrow(integrated_grape_wine_data), ")", "\n(PCC = ",format(cor(df$abs_pc2,df$cont_score),digits=2) ,")")) + 
        #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        labs(y="Smoke Taint Index",x=abs_pc2_label) +
        scale_color_brewer(palette="RdBu", direction=-1)
    print(h)
    dev.off()
    
    
    prepare_pca_loading_df = function(idx){
        order_entry = order(pca_result$rotation[,idx],decreasing=F)
        y = names(pca_result$rotation[,idx])[order_entry]
        y = sapply(y, function(x){
            gsub("_"," ",x)
        })
        y = sapply(y, function(x){
            gsub(")",") ",x)
        })
        y = factor(y,levels = y)
        
        x = pca_result$rotation[,idx][order_entry]
        pos = factor(pca_result$rotation[,idx][order_entry] >= 0, levels = c("TRUE","FALSE"))
        
        df = data.frame(x = x, y = y, pos = pos)
        return(df)
    }
    
    pc1_loading_df = prepare_pca_loading_df(1)
    pc2_loading_df = prepare_pca_loading_df(2)
    
    
    #PC1 Loading
    svg(filename_pc1_loading,width=6,height=8)
    p = ggplot(pc1_loading_df, aes(x=x, y=y, fill=pos)) +
        ggtitle("Loadings of PC 1") + 
        geom_col() + 
        theme_classic() + 
        xlim(xlim_loading1) + 
        geom_vline(xintercept = 0) + 
        xlab("Loading (PC 1)") + 
        ylab("VOC") + 
        theme(legend.position = "none") + 
        theme(plot.title = element_text(hjust = 0.5, size=14),
              axis.text.x = element_text(size=11),
              axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=11),
              axis.title.y = element_text(size=14))
    
    print(p)
    dev.off()
    
    
    #PC2 Loading
    svg(filename_pc2_loading,width=6,height=8)
    p = ggplot(pc2_loading_df, aes(x=x, y=y, fill=pos)) +
        ggtitle("Loadings of PC 2") + 
        geom_col() + 
        theme_classic() + 
        xlim(xlim_loading2) + 
        geom_vline(xintercept = 0) + 
        xlab("Loading (PC 2)") + 
        ylab("VOC") + 
        theme(legend.position = "none") + 
        theme(plot.title = element_text(hjust = 0.5, size=14),
              axis.text.x = element_text(size=11),
              axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=11),
              axis.title.y = element_text(size=14))
    
    print(p)
    dev.off()
}


plot_tsne = function(data, integrated_grape_wine_data,
                     approaches = c("none", "min_max","z","quantile"),
                     n_neighbors = c(2,5,10,20)){
    #TSNE
    find_region = function(x){
        x = strsplit(x,"|",fixed=T)[[1]][1]
        idx = which(x == data[,colidx_sample_name])
        if (length(idx) == 0) stop()
        return(data$`Region (AVA)`[idx[1]])
    }
    
    find_county = function(x){
        x = strsplit(x,"|",fixed=T)[[1]][1]
        idx = which(x == data[,colidx_sample_name])
        if (length(idx) == 0) stop()
        return(data$County[idx[1]])
    }
    
    score_tag = sapply(integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"], function(x){
        if (x <= 10) return("<= 10")
        if (x <= 20) return("10 < x <= 20")
        if (x <= 30) return("20 < x <= 30")
        if (x <= 40) return("30 < x <= 40")
        if (x <= 50) return("40 < x <= 50")
        if (x <= 60) return("50 < x <= 60")
        if (x <= 70) return("60 < x <= 70")
    })
    
    
    point_size = 3
    
    for (cur_approach in approaches){
        dir.create(paste0("output/1_general_analysis/tsne_",cur_approach))
        
        if (cur_approach == "none"){
            data_tsne = integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)]
            normalize_approach_descriptor = "(Not Normalizaed)"
        }else if (cur_approach == "min_max"){
            data_tsne = apply(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],2,function(x){
                return((x - min(x))/(max(x) - min(x)))  
            })
            normalize_approach_descriptor = "(Min-Max Normalizaed)"
        }else if (cur_approach == "z"){
            data_tsne = apply(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],2,function(x){
                return((x - mean(x))/sd(x))  
            })
            normalize_approach_descriptor = "(Z Normalizaed)"
        }else{
            new_data = t(normalize.quantiles(t(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)])))
            data_tsne = integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)]
            data_tsne[,] = new_data[,]
            normalize_approach_descriptor = "(Quantile Normalizaed)"
        }
        
        
        for (cur_n_neighbor in n_neighbors){
            cur_dir = paste0("output/1_general_analysis/tsne_",cur_approach,"/n_neighbor_",cur_n_neighbor)
            dir.create(cur_dir)
            
            n_neighbor_descriptor = paste0("(# of neighbor = ", cur_n_neighbor,")")
            
            tsne_title = paste0("TSNE (N = ", nrow(integrated_grape_wine_data), ")", n_neighbor_descriptor, "\n",normalize_approach_descriptor)
            
            tsne_res = tsne(data_tsne,perplexity = cur_n_neighbor,max_iter = 10000)
            
            df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], score_tag = score_tag)
            point_size = 3
            
            svg(paste0(cur_dir,"/tsne_index.svg"),width=6,height=4)
            h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = score_tag)) + geom_point(size=point_size) + 
                scale_size(guide="none")  +
                theme_classic() +
                ggtitle(tsne_title) + 
                #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      axis.title.y = element_text(size=12)) +
                labs(y="Dim 2",x="Dim 1",color="Smoke Taint Index")+
                scale_color_brewer(palette="RdBu", direction=-1)
            print(h)
            dev.off()
            
            
            
            
            df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                                 type = sapply(strsplit(rownames(integrated_grape_wine_data),"|",fixed=T),function(x){x[2]}))
            svg(paste0(cur_dir,"/tsne_Scale.svg"),width=6,height=4)
            h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
                scale_size(guide="none")  +
                theme_classic() +
                ggtitle(tsne_title) + 
                #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      axis.title.y = element_text(size=12)) +
                labs(y="Dim 2",x="Dim 1",color="Scale")+
                scale_color_brewer(palette="Set3")
            print(h)
            dev.off()
            
            df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                                 type = substr(rownames(integrated_grape_wine_data),3,4))
            svg(paste0(cur_dir,"/tsne_variety.svg"),width=6,height=4)
            h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
                scale_size(guide="none")  +
                theme_classic() +
                ggtitle(tsne_title) + 
                #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      axis.title.y = element_text(size=12)) +
                labs(y="Dim 2",x="Dim 1",color="Variety")
            print(h)
            dev.off()
            
            
            df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                                 type = sapply(rownames(integrated_grape_wine_data),find_region))
            svg(paste0(cur_dir,"/tsne_region.svg"),width=8,height=4)
            h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
                scale_size(guide="none")  +
                theme_classic() +
                ggtitle(tsne_title) + 
                #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      axis.title.y = element_text(size=12)) +
                labs(y="Dim 2",x="Dim 1",color="Region")
            print(h)
            dev.off()
            
            
            df_tsne = data.frame(dim1 = tsne_res[,1], dim2 = tsne_res[,2], 
                                 type = sapply(rownames(integrated_grape_wine_data),find_county))
            svg(paste0(cur_dir,"/tsne_county.svg"),width=6,height=4)
            h = ggplot(df_tsne, aes(x = dim1, y = dim2, color = type)) + geom_point(size=point_size) + 
                scale_size(guide="none")  +
                theme_classic() +
                ggtitle(tsne_title) + 
                #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
                theme(plot.title = element_text(hjust = 0.5),
                      axis.text.x = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.text.y = element_text(size=12),
                      axis.title.y = element_text(size=12)) +
                labs(y="Dim 2",x="Dim 1",color="County")+
                scale_color_brewer(palette="Set3")
            print(h)
            dev.off()
        }
    }
}

#Heatmap
plot_summary_heatmap = function(filename, integrated_grape_wine_data){
    find_region = function(x){
        x = strsplit(x,"|",fixed=T)[[1]][1]
        idx = which(x == data[,colidx_sample_name])
        if (length(idx) == 0) stop()
        return(data$`Region (AVA)`[idx[1]])
    }
    
    find_county = function(x){
        x = strsplit(x,"|",fixed=T)[[1]][1]
        idx = which(x == data[,colidx_sample_name])
        if (length(idx) == 0) stop()
        return(data$County[idx[1]])
    }
    
    prepareAnnotation = function(x, colorMapName, reverse=F, order_by_name = F){
        tab = table(x)
        tab = tab[order(tab,decreasing = T)]
        colors = brewer.pal(n = min(length(tab),brewer.pal.info[colorMapName,"maxcolors"]), colorMapName)
        if (reverse) colors = colors[length(colors):1]
        if (length(tab) > length(colors)){
            new_tab = tab[1:length(colors)]
            new_tab[length(new_tab)] = sum(tab) - sum(tab[1:(length(colors)-1)])
            names(new_tab)[length(new_tab)] = "others"
            
            new_x = x
            new_x[which(!(x %in% names(tab)[1:(length(colors)-1)]))] = "others"
        }else{
            new_tab = tab
            new_x = x
        }
        
        if (order_by_name){
            names(colors) = names(new_tab)[order(names(new_tab))]
        }else{
            names(colors) = names(new_tab)
        }
        
        return(list(x = new_x, colors = colors))
    }
    
    
    scale_annotation =  prepareAnnotation(sapply(strsplit(rownames(integrated_grape_wine_data),"|",fixed=T),function(x){x[2]}),"Pastel1")
    variety_annotation = prepareAnnotation(substr(rownames(integrated_grape_wine_data),3,4),"Set3")
    county_annotation = prepareAnnotation(sapply(rownames(integrated_grape_wine_data),find_county),"Pastel2")
    score_annotation = prepareAnnotation( sapply(integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"], function(x){
        if (x <= 10) return("<= 10")
        if (x <= 20) return("10 < x <= 20")
        if (x <= 30) return("20 < x <= 30")
        if (x <= 40) return("30 < x <= 40")
        if (x <= 50) return("40 < x <= 50")
        if (x <= 60) return("50 < x <= 60")
        if (x <= 70) return("60 < x <= 70")
    }), "RdBu", reverse = T, order_by_name = T)
    
    ha = HeatmapAnnotation(
        "Scale" = scale_annotation$x,
        "Variety" = variety_annotation$x,
        "County" = county_annotation$x,
        "Smoke Taint Index" = score_annotation$x,
        col = list("Scale" = scale_annotation$colors,
                   "Variety" = variety_annotation$colors,
                   "County" = county_annotation$colors,
                   "Smoke Taint Index" = score_annotation$colors),
        border = F
    )
    
    ###(TEST REGION)
    title=paste0("Heatmap (N = ", nrow(integrated_grape_wine_data), ")")
    testProfiles = apply(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],2,function(x){(x-min(x))/(max(x)-min(x))})
    svg(filename = filename, width = 10, height = 8)
    ht_list = Heatmap(t(testProfiles), 
                      col = colorRamp2(c(0, 0.5, 1), c("green", "white", "red")),
                      name = "Normalized Compound Concentration", clustering_distance_rows = "euclidean",
                      clustering_method_rows = "complete", column_gap = unit(1, "mm"), border=F,
                      column_title = title,show_column_names = F,show_row_names = T, top_annotation = ha, 
                      show_row_dend = F,
                      heatmap_legend_param = list(direction="horizontal",title_position="leftcenter",x=unit(0,"npc"), just = c("left", "bottom"))
    )
    draw(ht_list, merge_legend = F, heatmap_legend_side = "bottom", legend_grouping="original",
         annotation_legend_side = "right")
    dev.off()
}


analysis_hypothesis1 = function(integrated_grape_wine_data){
    #Hypothesis -- Does higher average VOCs (for each samples) yield higher smoke taint index?
    avg_val = apply(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],1,mean)
    pcc_val = cor(avg_val, integrated_grape_wine_data[,ncol(integrated_grape_wine_data)])
    t_score = (pcc_val * sqrt(length(avg_val)-2)) / sqrt(1 - pcc_val^2)
    p_val = pt(t_score, length(avg_val)-2, lower.tail = F)
    
    print("Hypothesis 1 -- Does higher average VOCs (for each samples) yield higher smoke taint index?")
    print(paste0("PCC (Avg VOC intensities VS smoke taint index) = ",format(pcc_val,digits=2), " (p = ",format(p_val,digits=2),")"))
}

analysis_hypothesis2 = function(integrated_grape_wine_data){
    #Hypothesis -- Are wines from Napa with higher smoke taint index?
    smoke_taint_index = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
    county_info = sapply(rownames(integrated_grape_wine_data),find_county)
    pval = t.test(smoke_taint_index[county_info=="Napa"],smoke_taint_index[county_info!="Napa"],alternative = "greater")$p.value
    print("Hypothesis 2 -- Are wines from Napa with higher smoke taint index?")
    print(paste0("Avg Smoke Taint Index (Napa) = ", format(mean(smoke_taint_index[county_info=="Napa"]), digits=4)))
    print(paste0("Avg Smoke Taint Index (Others) = ", format(mean(smoke_taint_index[county_info!="Napa"]), digits=4)))
    print(paste0("p = ",format(p_val,digits=2),""))
}


#Uni-variate analysis
plot_cor_between_features_and_smoke_taint = function(filename, integrated_grape_wine_data){
    pcc_smoke_idx = apply(integrated_grape_wine_data, 2, function(x){cor(x,integrated_grape_wine_data[,ncol(integrated_grape_wine_data)])})
    spearman_smoke_idx = apply(integrated_grape_wine_data, 2, function(x){cor(x,integrated_grape_wine_data[,ncol(integrated_grape_wine_data)],method = "spearman")})
    log_pcc_smoke_idx = apply(log(integrated_grape_wine_data+1), 2, function(x){cor(x,integrated_grape_wine_data[,ncol(integrated_grape_wine_data)])})
    
    correlation_matrix = cbind(pcc_smoke_idx, log_pcc_smoke_idx, spearman_smoke_idx)
    colnames(correlation_matrix) = c("PCC","PCC (Log value)","Spearman")
    feature_name = rownames(correlation_matrix)
    feature_name = sapply(feature_name, function(x){gsub("_"," ",x)})
    feature_name = sapply(feature_name, function(x){gsub(")",") ",x)})
    rownames(correlation_matrix) = feature_name
    correlation_matrix = correlation_matrix[-nrow(correlation_matrix),]
    
    
    svg(filename,width=7,height=8)
    h = Heatmap(correlation_matrix, 
            name = "correlation", 
            col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
            column_title = "Correlation (Features VS Smoke Taint Index)",
            column_title_gp = gpar(fontsize = 14),
            cluster_columns = F,
            cluster_rows = T,
            show_row_dend = F,
            row_names_side = "left",
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(formatC(correlation_matrix[i,j], digits=2, format="f"), x, y)
            }, 
            show_column_names = FALSE, 
            top_annotation = HeatmapAnnotation(
                text = anno_text(colnames(correlation_matrix), rot = 0, just = "center"),
                annotation_height = max_text_height(colnames(correlation_matrix))
            ), 
            show_heatmap_legend = FALSE)
    draw(h)
    dev.off()
    
}

plot_cor = function(filename, width, height, x,y, title, xlab){
    svg(filename, width, height)
    df = data.frame(x = x, 
                    y = y)
    point_size = 3
    pcc = cor(x,y)
    log_pcc = cor(log(x+1),y)
    spearman = cor(x,y,method="spearman")
    h = ggplot(df, aes(x = x, y = y)) + 
        geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle(title) + 
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        labs(y="Smoke Taint Index (out of 100)",x=xlab, color="value")+
        scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100)) + 
        annotate(geom="text", x=0, y=90, 
                 label=paste0("N = ",length(x)," (Wines)\nPCC = ",format(pcc,digits=2),
                              "\nPCC (Log(Conc+1)) = ",format(log_pcc,digits=2),
                              "\nSpearman = ",format(spearman,digits=2)), 
                 hjust=0,
                 color="black")
    print(h)
    dev.off()
}

plot_cor_log = function(filename, width, height, x,y){
    svg(filename, width, height)
    df = data.frame(x = log(x+1), 
                    y = y)
    point_size = 1
    pcc = cor(x,y)
    log_pcc = cor(log(x+1),y)
    spearman = cor(x,y,method="spearman")
    h = ggplot(df, aes(x = x, y = y)) + 
        geom_point(size=point_size) + 
        scale_size(guide="none")  +
        theme_classic() +
        ggtitle("Log(Conc+1) VS Index") + 
        theme(plot.title = element_text(hjust = 0.5,size=10),
              axis.text.x = element_text(size=10),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size=10),
              axis.title.y = element_blank()) +
        scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100))
    print(h)
    dev.off()
}


#The script part for general analysis

#Histogram of smoke taint index of wines
plot_histogram_smoke_taint("output/1_general_analysis/histogram_smoke_taint.svg",integrated_grape_wine_data)

#PCA Analysis
pca_result = prcomp(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],center = T,scale. = T)
pca_result_grape_voc = prcomp(integrated_grape_wine_data[,1:20],center = T,scale. = T)

plot_pca_results("output/1_general_analysis/pca_pc1_pc2.svg",
                 "output/1_general_analysis/pca_pc1_score.svg",
                 "output/1_general_analysis/pca_pc2_score.svg",
                 "output/1_general_analysis/pc1_loading.svg",
                 "output/1_general_analysis/pc2_loading.svg",
                 pca_result,
                 integrated_grape_wine_data,
                 c(-0.05,0.25),
                 c(-0.35,0.25)
                 )
plot_pca_results("output/1_general_analysis/pca_pc1_pc2_grape_voc_only.svg",
                 "output/1_general_analysis/pca_pc1_score_grape_voc_only.svg",
                 "output/1_general_analysis/pca_pc2_score_grape_voc_only.svg",
                 "output/1_general_analysis/pc1_loading_grape_voc_only.svg",
                 "output/1_general_analysis/pc2_loading_grape_voc_only.svg",
                 pca_result_grape_voc,
                 integrated_grape_wine_data,
                 c(0,0.4),
                 c(-0.40,0.40)
                 )


plot_tsne(data, integrated_grape_wine_data, n_neighbors = c(2,5))
plot_summary_heatmap("output/1_general_analysis/summary_heatmap.svg", integrated_grape_wine_data)

#Analysis
analysis_hypothesis1(integrated_grape_wine_data)
analysis_hypothesis2(integrated_grape_wine_data)

plot_cor_between_features_and_smoke_taint("output/1_general_analysis/heatmap_univariate_analysis.svg",integrated_grape_wine_data)
plot_cor("output/1_general_analysis/cor_free_m-cresol.svg",
         5,5,
         x = integrated_grape_wine_data[,"free_m-cresol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)],
         paste0("Free m-cresol VS Smoke Taint Index"),
         "Concentration (ug/L))"
         )
plot_cor_log("output/1_general_analysis/cor_free_m-cresol_log.svg",
         2,2,
         x = integrated_grape_wine_data[,"free_m-cresol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
)


plot_cor("output/1_general_analysis/cor_total_syringol_grape.svg",
         5,5,
         x = integrated_grape_wine_data[,"(grape)total_syringol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)],
         paste0("Total syringol (Grape) VS Smoke Taint Index"),
         "Concentration (ug/Kg))"
)
plot_cor_log("output/1_general_analysis/cor_total_syringol_grape_log.svg",
         2,2,
         x = integrated_grape_wine_data[,"(grape)total_syringol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
)


plot_cor("output/1_general_analysis/cor_total_4-methylsyringol_grape.svg",
         5,5,
         x = integrated_grape_wine_data[,"(grape)total_4-methylsyringol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)],
         paste0("Total 4-Methylsyringol (Grape) VS Smoke Taint Index"),
         "Concentration (ug/Kg))"
)
plot_cor_log("output/1_general_analysis/cor_total_4-methylsyringol_grape_log.svg",
         2,2,
         x = integrated_grape_wine_data[,"(grape)total_4-methylsyringol"], 
         y = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
)



#Linear regression

features = integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)]
scores = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]
linear_model = lm(scores~features)

lm_predict = predict(linear_model)
cor(lm_predict,scores)

#Lasso regression
lambdas = 10^seq(3,-3,-0.0001)

source("code/functions_smoke_taint_prediction.r")
N_Trial = 10
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
svg("output/1_general_analysis/model_performance.svg",width=8,height=1.5)
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
svg("output/1_general_analysis/model_performance_spearman.svg",width=8,height=1.5)
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

