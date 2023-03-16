#Smoke Taint Index Prediction Project

#Please set your working directory here.
WORKING_DIR = "C:/Users/Bigghost/Documents/GitHub/smoke_taint_prediction/"
setwd(WORKING_DIR)

#Create the output directory for storing the output figures
dir.create("output")
dir.create("output/1_general_analysis")
dir.create("output/1_general_analysis/1a_descriptions")

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
library(dendextend)

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
plot_histogram_smoke_taint = function(filename, integrated_grape_wine_data, width = 4, height = 4){
    #Hist of wines
    svg(filename,width,height)
    dd<-data.frame(x=integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"])
    h = ggplot(dd, aes(x=x))+ geom_histogram(breaks=seq(0,70,10),fill="white",color="black",boundary=0) + theme_classic() +
        ggtitle(paste0("Distribution of Smoke Taint Index")) + 
        #geom_hline(yintercept=median(Cor_T_Merged), color = "red") +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size=12),
              axis.title.x = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.title.y = element_text(size=12)) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
        #scale_fill_manual(values=ColorSet)+
        labs(y="Count",x="Smoke Taint Index") + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 20))  + 
        scale_x_continuous(expand = c(0, 0), breaks=seq(0,70,10)) +
        stat_bin(breaks=seq(0,70,10), geom="text", color='black', aes(label=after_stat(count)),
                 position=position_stack(vjust = 0.5), boundary=0)  +  
        annotate(geom="text", x=55, y=18, label=paste0("N = ",nrow(integrated_grape_wine_data)), hjust=0,
                 color="black")
    print(h)
    dev.off()
}



#Heatmap
find_region = function(x){
    x = strsplit(x,"|",fixed=T)[[1]][1]
    idx = which(x == data[,colidx_sample_name])
    if (length(idx) == 0) stop()
    return(data$`Region (AVA)`[idx[1]])
}

old_county_name = c("Yolo County","Lake County","Napa","Sonoma","Mendocino","Monterey ","Polk","Santa Barbara","Yamhill, OR")
new_county_name = c("Yolo","Lake","Napa","Sonoma","Mendocino","Monterey","Polk, OR","Santa Barbara","Yamhill, OR")



find_county = function(x){
    x = strsplit(x,"|",fixed=T)[[1]][1]
    idx = which(x == data[,colidx_sample_name])
    if (length(idx) == 0) stop()
    result = data$County[idx[1]]
    idx_replace = which(result == old_county_name)
    return(new_county_name[which(result == old_county_name)])
}

prepare_score_annotation = function(integrated_grape_wine_data){
    entries = sapply(integrated_grape_wine_data[,"Sensory_ash (smoke taint out of 100)"], function(x){
        if (x <= 10) return("<= 10")
        if (x <= 20) return("10 < x <= 20")
        if (x <= 30) return("20 < x <= 30")
        if (x <= 40) return("30 < x <= 40")
        if (x <= 50) return("40 < x <= 50")
        if (x <= 60) return("50 < x <= 60")
        if (x <= 70) return("60 < x <= 70")
    })
    colors = brewer.pal(7, "RdBu")
    colors = colors[length(colors):1]
    names(colors) = unique(entries)
    names(colors) = names(colors)[order(names(colors))]
    return (list(x = factor(entries,levels=names(colors)), colors = colors))
}


plot_summary_heatmap = function(filename, integrated_grape_wine_data){
    prepare_scale_fermentation_annotation = function(integrated_grape_wine_data){
        entries = sapply(strsplit(rownames(integrated_grape_wine_data),"|",fixed=T),function(x){x[2]})
        old_name = c("Wine_Bucket","Wine_TJ","Wine_Pilot")
        new_name = c("Bucket","T.J. Fermentor","Pilot")
        
        for (i in 1 : length(old_name)){
            entries = gsub(old_name[i],new_name[i],entries)
        }
        
        colors = c("#CC9900","#A6A6A6","#B4C7E7")
        names(colors) = new_name
        return (list(x = factor(entries,levels=names(colors)), colors = colors))
    }
    
    prepare_variety_annotation = function(integrated_grape_wine_data){
        entries = substr(rownames(integrated_grape_wine_data),3,4)
        old_name = c("CS","PN","TO","SB","CF","ME","PV")
        new_name = c("Cabernet Sauvignon","Pinot Noir","Torrontes","Sauvignon blanc","Cabernet Franc","Merlot","Petit Verdot")
        
        for (i in 1 : length(old_name)){
            entries[which(entries==old_name[i])] = new_name[i]
        }
        entries[which(!(entries %in% new_name))] = "Others"
        
        colors = c("#660033","#660066","#92D050","#C5E0B4","#666699","#990099","#CC0000","#D9D9D9")
        names(colors) = c(new_name,"Others")
        return (list(x = factor(entries,levels=names(colors)), colors = colors))
    }
    
    prepare_county_annotation = function(integrated_grape_wine_data){
        entries = sapply(rownames(integrated_grape_wine_data),find_county)
        colors = brewer.pal(9, "Set3")
        names(colors) = c("Yolo","Lake","Napa","Sonoma","Mendocino","Monterey","Santa Barbara","Polk, OR","Yamhill, OR")
        return (list(x = factor(entries,levels=names(colors)), colors = colors))
    }
    
    scale_annotation =  prepare_scale_fermentation_annotation(integrated_grape_wine_data)
    variety_annotation = prepare_variety_annotation(integrated_grape_wine_data)
    county_annotation = prepare_county_annotation(integrated_grape_wine_data)
    score_annotation = prepare_score_annotation(integrated_grape_wine_data)
    
    ha = HeatmapAnnotation(
        "County" = county_annotation$x,
        "Variety" = variety_annotation$x,
        "Scale" = scale_annotation$x,
        "Smoke Taint Index" = score_annotation$x,
        col = list("County" = county_annotation$colors,
                   "Variety" = variety_annotation$colors,
                   "Scale" = scale_annotation$colors,
                   "Smoke Taint Index" = score_annotation$colors),
        annotation_name_side = "left",
        border = F
    )
    
    testProfiles = apply(integrated_grape_wine_data[,1:(ncol(integrated_grape_wine_data)-1)],2,
                         function(x){
                             y = x
                             return((y-mean(y))/sd(y))
                             })
    #testProfiles = t(testProfiles)
    
    total_voc_type_tag = "(Total VOC)"
    free_voc_type_tag = "(Free VOC)"
    
    grape_testProfiles = testProfiles[,which(grepl("(grape)",colnames(testProfiles)))]
    colnames(grape_testProfiles) = gsub("(grape)","",colnames(grape_testProfiles),fixed=T)
    wine_testProfiles = testProfiles[,which(!grepl("(grape)",colnames(testProfiles)))]
    
    grape_voc_type_labels = rep(total_voc_type_tag, ncol(grape_testProfiles))
    grape_voc_type_labels[grep("free_",colnames(grape_testProfiles))] = free_voc_type_tag
    colnames(grape_testProfiles) = gsub("free_","",colnames(grape_testProfiles),fixed=T)
    colnames(grape_testProfiles) = gsub("total_","",colnames(grape_testProfiles),fixed=T)
    
    wine_voc_type_labels = rep(total_voc_type_tag, ncol(wine_testProfiles))
    wine_voc_type_labels[grep("free_",colnames(wine_testProfiles))] = free_voc_type_tag
    colnames(wine_testProfiles) = gsub("free_","",colnames(wine_testProfiles),fixed=T)
    colnames(wine_testProfiles) = gsub("total_","",colnames(wine_testProfiles),fixed=T)
    
    dend = as.dendrogram(hclust(dist(testProfiles, method="canberra")))
    dend = color_branches(dend, k = 2)
    
    hclust_col = hclust(dist(t(wine_testProfiles[,which(wine_voc_type_labels==total_voc_type_tag)]), method="canberra"))
    hclust_col_order = hclust_col$order
    
    
    wine_total_vocs = wine_testProfiles[,which(wine_voc_type_labels==total_voc_type_tag)]
    wine_total_vocs = wine_total_vocs[,hclust_col_order]
    
    grape_total_vocs = grape_testProfiles[,which(grape_voc_type_labels==total_voc_type_tag)]
    grape_total_vocs = grape_total_vocs[,colnames(wine_total_vocs)]
    
    wine_free_vocs = wine_testProfiles[,which(wine_voc_type_labels==free_voc_type_tag)]
    wine_free_vocs = wine_free_vocs[,colnames(wine_total_vocs)]
    
    grape_free_vocs = grape_testProfiles[,which(grape_voc_type_labels==free_voc_type_tag)]
    grape_free_vocs = grape_free_vocs[,colnames(wine_total_vocs)]
    
    
    wine_testProfiles = cbind(wine_total_vocs, wine_free_vocs)
    grape_testProfiles = cbind(grape_total_vocs, grape_free_vocs)
    wine_voc_type_labels = c(rep(total_voc_type_tag,ncol(wine_total_vocs)), rep(free_voc_type_tag, ncol(wine_free_vocs)))
    grape_voc_type_labels = c(rep(total_voc_type_tag,ncol(grape_total_vocs)), rep(free_voc_type_tag, ncol(grape_free_vocs)))
    
    wine_voc_type_labels = factor(wine_voc_type_labels, levels = unique(wine_voc_type_labels))
    grape_voc_type_labels = factor(grape_voc_type_labels, levels = unique(grape_voc_type_labels))
    
    ht_wine = Heatmap(t(wine_testProfiles), 
                      col = colorRamp2(seq(-4,4,length.out=12), hcl.colors(12,"viridis")),
                      name = "Z-Normalized Compound Concentration", 
                      cluster_columns = dend,
                      column_split = 2,
                      row_split = wine_voc_type_labels,
                      clustering_method_columns = "complete", column_gap = unit(1, "mm"), border=F,
                      cluster_rows = F,
                      show_column_names = F,
                      show_row_names = F, 
                      cluster_row_slices = F,
                      top_annotation = ha, 
                      show_row_dend = F,
                      show_column_dend = T,
                      column_title = " ",
                      row_names_side = "left",
                      heatmap_legend_param = list(plot=F,
                                                  direction="horizontal",
                                                  title_position="leftcenter",
                                                  x=unit(0,"npc"), 
                                                  just = c("left", "bottom")),
                      row_title = "Wine VOC",
                      left_annotation = rowAnnotation(foo = anno_block(gp = gpar(col=NA),
                                                                          labels=c("(Total VOC)","(Free VOC)"), 
                                                                          labels_gp = gpar(col="black")),
                                                      text = anno_text(rownames(t(wine_testProfiles)))))
    
    ht_grape = Heatmap(t(grape_testProfiles), 
                       col = colorRamp2(seq(-4,4,length.out=12), hcl.colors(12,"viridis")),
                       name = "Z-Normalized Compound Concentration", 
                       cluster_columns = dend,
                       column_split = 2,
                       row_split = wine_voc_type_labels,
                       clustering_method_columns = "complete", column_gap = unit(1, "mm"), border=F,
                       cluster_rows = F,
                       show_column_names = F,
                       show_row_names = F, 
                       cluster_row_slices = F,
                       show_row_dend = F,
                       show_column_dend = T,
                       column_title = " ",
                       row_names_side = "left",
                       heatmap_legend_param = list(direction="horizontal",title_position="leftcenter",x=unit(0,"npc"), just = c("left", "bottom")),
                       row_title = "Grape VOC",
                       left_annotation = rowAnnotation(foo = anno_block(gp = gpar(col=NA),
                                                                        labels=c("(Total VOC)","(Free VOC)"), 
                                                                        labels_gp = gpar(col="black")),
                                                       text = anno_text(rownames(t(wine_testProfiles)))))
    
    
    svg(filename = filename, width = 8, height = 10)
    h = draw(ht_wine %v% ht_grape, merge_legend = F, heatmap_legend_side = "bottom", legend_grouping="original",
         annotation_legend_side = "right")
    print(h)
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
    print(paste0("p = ",format(pval,digits=2),""))
}



#Histogram of smoke taint index of wines
plot_histogram_smoke_taint("output/1_general_analysis/1a_descriptions/histogram_smoke_taint.svg",integrated_grape_wine_data, 4, 4)

#Table for Pie chart
table_for_pie_chart = matrix(NA,nrow=nrow(integrated_grape_wine_data),ncol=3)
rownames(table_for_pie_chart) = rownames(integrated_grape_wine_data)
colnames(table_for_pie_chart) = c("scale","variety","county")

table_for_pie_chart[,"scale"] = sapply(strsplit(rownames(integrated_grape_wine_data),"|",fixed=T),function(x){x[2]})
table_for_pie_chart[,"variety"] = substr(rownames(integrated_grape_wine_data),3,4)
table_for_pie_chart[,"county"] = sapply(rownames(integrated_grape_wine_data),find_county)

write.csv(table_for_pie_chart, "output/1_general_analysis/1a_descriptions/table_for_pie_chart.csv")

#Heatmap
plot_summary_heatmap("output/1_general_analysis/1a_descriptions/summary_heatmap.svg", integrated_grape_wine_data)


#Boxplot of smoke taint index for binned by Scale/Variety/County
#County

smoke_taint_index = integrated_grape_wine_data[,ncol(integrated_grape_wine_data)]


county_color = brewer.pal(9, "Set3")
names(county_color) = c("Yolo","Lake","Napa","Sonoma","Mendocino","Monterey","Santa Barbara","Polk, OR","Yamhill, OR")
median_smoke_taint_index_county = sapply(names(county_color), function(x){
    county_list = sapply(rownames(integrated_grape_wine_data),find_county)
    return(median(smoke_taint_index[which(county_list==x)]))
})
county_color = county_color[order(median_smoke_taint_index_county,decreasing = F)]

df_county = data.frame(county = factor(sapply(rownames(integrated_grape_wine_data),find_county), levels=names(county_color)), 
                       smoke_taint_index = smoke_taint_index)
# Compute the analysis of variance
res.aov <- aov(smoke_taint_index ~ county, data = df_county)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov)



ggplot(data = df_county, aes(x=smoke_taint_index,y=county,fill=county)) + geom_boxplot() + 
    theme_classic() + 
    theme(
        legend.position = "none"
    ) +
    xlab("Smoke Taint Index") + 
    ylab("County") + 
    scale_fill_manual(values=county_color)



#Map Plot
find_gps = function(x){
    x = strsplit(x,"|",fixed=T)[[1]][1]
    idx = which(x == data[,colidx_sample_name])
    if (length(idx) == 0) stop()
    return(data$GPS[idx[1]])
}
gps_string = sapply(rownames(integrated_grape_wine_data),find_gps)
gps_string_N = sapply(gps_string,function(x){
    strsplit(x,"[N,]")[[1]][1]
})
gps_string_W = sapply(gps_string,function(x){
    tmp = strsplit(x,"[N,]")[[1]]
    return(tmp[length(tmp)])
})
parse_gps_string = function(x){
    x = gsub("[^0-9A-Z¢X'\" .]","",x)
    if (!is.na(as.numeric(x))) return(as.numeric(x))
    x = gsub("[NW]","",x)
    while(x != gsub("^ ","",x)) x = gsub("^ ","",x)
    strsplit_x = strsplit(x,"[¢X'\" .]")[[1]]
    strsplit_x = as.numeric(strsplit_x)
    return(strsplit_x[1] + strsplit_x[2]/60 + strsplit_x[3]/3600)
}
gps_string_N_numeric = sapply(gps_string_N,parse_gps_string)
gps_string_W_numeric = sapply(gps_string_W,parse_gps_string)

target_counties = map_data("county")
target_counties = subset(target_counties, region %in% c("california","oregon"))
target_county_list = sapply(names(county_color),function(x){
    return(tolower(strsplit(x,",",fixed=T)[[1]][1]))
})
target_counties$color = sapply(1:nrow(target_counties),function(i){
    x = target_counties$subregion[i]
    y = target_counties$region[i]
    
    if (x %in% target_county_list){
        if ((y == "oregon" && (x %in% c("yamhill","polk"))) || y == "california"){
            idx = which(x == target_county_list)
            return(names(target_county_list)[idx])
        }
        return(NA)
    }
    return(NA)
})
target_counties$color = factor(target_counties$color, levels=c("Yolo","Lake","Napa","Sonoma","Mendocino","Monterey","Santa Barbara","Polk, OR","Yamhill, OR"))

smoke_taint_index_annotation = prepare_score_annotation(integrated_grape_wine_data)

df_gps_smoke_taint = data.frame(long = -gps_string_W_numeric, lat = gps_string_N_numeric, smoke_taint_index = smoke_taint_index_annotation$x)
new_county_color = adjustcolor(county_color,alpha.f=0.4)
names(new_county_color) = names(target_county_list)

svg("output/1_general_analysis/1a_descriptions/map_large.svg",width=6,height=10)
ggplot(target_counties, aes(long, lat, group = group)) +
    geom_polygon(data = target_counties, aes(fill = color), color="#CCCCCC") + 
    geom_point(data = df_gps_smoke_taint, aes(x=long,y=lat,color=smoke_taint_index),size=2.8,
               inherit.aes = FALSE) + 
    geom_point(data = df_gps_smoke_taint, aes(x=long,y=lat,color="black"),size=3,shape=21,
               inherit.aes = FALSE) + 
    coord_quickmap() + 
    theme(legend.position = "right",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c(new_county_color)) +
    scale_color_manual(values=smoke_taint_index_annotation$colors) + 
    xlab("Longtitude") + 
    ylab("Latitude") + 
    labs(fill = "County") + 
    labs(color = "Smoke Taint Index")
dev.off()

svg("output/1_general_analysis/1a_descriptions/map_small.svg",width=3,height=3)
selected_target_counties = target_counties[which(target_counties$subregion %in% c("napa","sonoma","yolo","lake","mendocino")),]
ggplot(selected_target_counties, aes(long, lat, group = group)) +
    geom_polygon(data = selected_target_counties, aes(fill = color), color="#CCCCCC") + 
    geom_point(data = df_gps_smoke_taint, aes(x=long,y=lat,color=smoke_taint_index),size=3,
               inherit.aes = FALSE) + 
    coord_quickmap() + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c(new_county_color)) +
    scale_color_manual(values=smoke_taint_index_annotation$colors) + 
    xlab("Longtitude") + 
    ylab("Latitude") + 
    labs(fill = "County") + 
    labs(color = "Smoke Taint Index") + 
    xlim(c(-124,-121)) +
    ylim(c(38,40))

dev.off()



ca_counties_selected = ca_counties[which(ca_counties$subregion %in% c("napa","sonoma","yolo","lake","mendocino","monterey","santa barbara","polk")),]
or_counties_selected = or_counties[which(or_counties$subregion %in% c("yamhill")),]

data_selected = rbind(ca_counties,or_counties)





