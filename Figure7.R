#install.packages('reticulate')
#reticulate::install_miniconda()
#reticulate::conda_install('r-reticulate', 'python-kaleido')
#reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
#reticulate::use_miniconda('r-reticulate')

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(psych)
library(plotly)
library(eulerr)
library(ggrepel)
library(reticulate)

read_maxlfq = function(file_path){
  out_data = fread(file_path) %>%
    filter(V1 != "")
}

count_id_num = function(report_data, miss_num, col_name, method_name, level_name){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-147):col_num)
  
  col_name <- enquo(col_name)
  report_data_group = report_data_gather %>%
    filter(!is.na(quant)) %>%
    group_by(!!col_name) %>%
    summarise(n=n()) %>%
    ungroup()
  colnames(report_data_group)[1] = "id"
  
  report_data_group_out = report_data_group %>%
    filter(n >= miss_num)
  
  report_data_filter = report_data %>%
    filter(!!col_name %in% report_data_group$id)
  
  report_data_filter_gather = gather(report_data_filter, raw_files, quant, (col_num-147):col_num)
  report_data_filter_gather_group = report_data_filter_gather%>%
    filter(!is.na(quant)) %>%
    group_by(raw_files) %>%
    summarise(n=n()) %>%
    ungroup()
  colnames(report_data_filter_gather_group)[1] = "V1"
  report_data_filter_gather_group$V1 = gsub("\\\\", "", report_data_filter_gather_group$V1)
  report_data_filter_gather_group$V1 = str_replace_all(report_data_filter_gather_group$V1, "G:diaPASEF_PandeyRAW_DIA", "")
  report_data_filter_gather_group$method = method_name
  report_data_filter_gather_group$Cleavage = level_name
  return(report_data_filter_gather_group)
}

### Summarize average
average_summarize = function(report_data, one_class_data, miss_ratio){
  report_data_filter = missing_summarize_one_class(report_data, one_class_data, miss_ratio)
  report_data_used = report_data %>%
    filter(V1 %in% report_data_filter$V1)
  col_num = ncol(report_data_used)
  report_data_gather = gather(report_data_used, raw_files, quant, (col_num-147):col_num)
  report_data_group = report_data_gather %>%
    filter(!is.na(quant)) %>%
    group_by(raw_files) %>%
    summarise(n=n()) %>%
    ungroup()
  
  return(report_data_group)
}

### Missing values estimate
missing_summarize = function(report_data, method){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-147):col_num)
  
  report_data_group = report_data_gather %>%
    filter(!is.na(quant)) %>%
    group_by(V1) %>%
    summarise(n=n()) %>%
    ungroup()
  
  report_data_group_out_25 = report_data_group %>%
    filter(n >= 37)
  report_data_group_out_50 = report_data_group %>%
    filter(n >= 74)
  report_data_group_out_100 = report_data_group %>%
    filter(n == 148)
  
  report_data_group_out_25_pure = report_data_group_out_25 %>%
    filter(!V1 %in% report_data_group_out_50$V1)
  report_data_group_out_50_pure = report_data_group_out_50 %>%
    filter(!V1 %in% report_data_group_out_100$V1)
  report_data_group_out_other = report_data_group %>%
    filter(!V1 %in% report_data_group_out_25$V1)
  
  id_portion = c("100%", ">50%", ">25%", "<25%")
  id_num = c(nrow(report_data_group_out_100), nrow(report_data_group_out_50_pure), nrow(report_data_group_out_25_pure), nrow(report_data_group_out_other))
  
  out_data = data_frame(id_portion, id_num)
  out_data$method = method
  
  return(out_data)
}

### Missing values estimate
missing_summarize_one_class = function(report_data, class_file, thres_pass){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, sample_name, quant, (col_num-147):col_num)
  report_data_gather$sample_name = str_replace_all(report_data_gather$sample_name, " ", "_")
  
  report_data_gather_con = left_join(report_data_gather, class_file, by="sample_name")
  report_data_group = report_data_gather_con %>%
    filter(!is.na(quant)) %>%
    group_by(V1, condition) %>%
    summarise(n=n()) %>%
    ungroup()
  
  fig4_lib_exp_one_class_num = class_file%>%
    group_by(condition) %>%
    summarise(num=n())
  report_data_group_num = left_join(report_data_group, fig4_lib_exp_one_class_num, by="condition")
  report_data_group_num$ratio = report_data_group_num$n/report_data_group_num$num
  report_data_group_num_pass = report_data_group_num %>%
    filter(ratio >= thres_pass) %>%
    distinct(V1)
  

  return((report_data_group_num_pass))
}

# prepare analysist
generate_over_anno = function(exp_anno){
  exp_anno$condition = ifelse(str_detect(exp_anno$sample_name, "Epithelium"), "Epithelium", 
                              ifelse(str_detect(exp_anno$sample_name, "T_cell"), "T-cell zone", 
                                     ifelse(str_detect(exp_anno$sample_name, "mantle_zone") | str_detect(exp_anno$sample_name, "manle_zone"), "Mantel zone", "Germinal center")))
  
  return(exp_anno)
}

generate_detail_anno = function(exp_anno){
  exp_anno$condition = ifelse(str_detect(exp_anno$sample_name, "Epithelium"), "Epithelium", 
                              ifelse(str_detect(exp_anno$sample_name, "T_cell"), "T_cell", 
                                     ifelse(str_detect(exp_anno$sample_name, "dark_zone"), "dark_zone", 
                                            ifelse(str_detect(exp_anno$sample_name, "grey_zone"), "grey_zone", 
                                                   ifelse(str_detect(exp_anno$sample_name, "light_zone"), "light_zone", 
                                                          "mantle_zone")))))
  return(exp_anno)
}

fig3_lib_exp = fread("./lowInputData/Tonsil_lib_result/experiment_annotation.tsv") %>%
  filter(str_ends(file, ".d"))
fig3_lib_pg_matrix = fread("./lowInputData/Tonsil_lib_result/diann-output/report.pg_matrix.tsv")
fig3_lib_exp_over_anno = generate_over_anno(fig3_lib_exp)
fig3_lib_exp_detail_anno = generate_detail_anno(fig3_lib_exp)
write.table(fig3_lib_exp_over_anno, "./lowInputData/Tonsil_lib_result/experiment_annotation_over.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(fig3_lib_exp_detail_anno, "./lowInputData/Tonsil_lib_result/experiment_annotation_detail.tsv",
            quote = F, row.names = F, sep = "\t")
fig3_lib_exp_detail_anno_germ = fig3_lib_exp_detail_anno %>%
  filter(condition %in% c("dark_zone", "grey_zone", "light_zone"))
fig3_lib_pg_matrix_gather = gather(fig3_lib_pg_matrix, "file", "value", 6:153)
fig3_lib_pg_matrix_gather = fig3_lib_pg_matrix_gather %>%
  filter(file %in% fig3_lib_exp_detail_anno_germ$file)
fig3_lib_pg_matrix_gather_spread = spread(fig3_lib_pg_matrix_gather, file, value)
write.table(fig3_lib_exp_detail_anno_germ, "./lowInputData/Tonsil_lib_result/experiment_annotation_germ.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(fig3_lib_pg_matrix_gather_spread, "./lowInputData/Tonsil_lib_result/diann-output/report.pg_matrix_germ.tsv",
            quote = F, row.names = F, sep = "\t")

fig4_lib_exp = fread("./lowInputData/Direct_lib_result/experiment_annotation.tsv") %>%
  filter(str_ends(file, ".d"))
fig4_lib_exp_over_anno = generate_over_anno(fig4_lib_exp)
fig4_lib_exp_detail_anno = generate_detail_anno(fig4_lib_exp)
write.table(fig4_lib_exp_over_anno, "./lowInputData/Direct_lib_result/experiment_annotation_over.tsv",
            quote = F, row.names = F, sep = "\t")
write.table(fig4_lib_exp_detail_anno, "./lowInputData/Direct_lib_result/experiment_annotation_detail.tsv",
            quote = F, row.names = F, sep = "\t")

# missing value one class
single_cell_fig4_lib_pro_result = read_maxlfq("./lowInputData/Direct_lib_result/diann-output/protein_maxlfq.tsv")
single_cell_fig3_lib_pro_result = read_maxlfq("./lowInputData/Tonsil_lib_result/diann-output/protein_maxlfq.tsv")
single_cell_diann_pro_result = read_maxlfq("./lowInputData/DIANN_results_Fig4_original_study/protein_maxlfq.tsv")

#!!! The number changed since I update the annotation file. One file and sample is mapped to a wrong group in the begining.
fig4_lib_exp_one_class = fread("./lowInputData/Direct_lib_result/experiment_annotation_over.tsv") %>%
  select(sample_name, condition)
fig4_lib_exp_one_class$sample_name = paste("Bluto_230", fig4_lib_exp_one_class$sample_name, sep = "")

fig4_lib_miss_num = list()
fig3_lib_miss_num = list()
paper_miss_num = list()
non_miss_ratio = list()
for (num in 0:10) {
  non_miss_ratio = append(non_miss_ratio, paste(num*10, "%", sep=""))
  fig4_lib_miss_num = append(fig4_lib_miss_num, nrow(missing_summarize_one_class(single_cell_fig4_lib_pro_result, fig4_lib_exp_one_class, num/10)))
  fig3_lib_miss_num = append(fig3_lib_miss_num, nrow(missing_summarize_one_class(single_cell_fig3_lib_pro_result, fig4_lib_exp_one_class, num/10)))
  paper_miss_num = append(paper_miss_num, nrow(missing_summarize_one_class(single_cell_diann_pro_result, fig4_lib_exp_one_class, num/10)))
}
single_cell_missing_class_num = data.frame(`Anuar et al.` = array(unlist(paper_miss_num)), 
                                           `FragPipe high-input Lib` = array(unlist(fig3_lib_miss_num)),
                                           `FragPipe` = array(unlist(fig4_lib_miss_num)),
                                           ratio = array(unlist(non_miss_ratio)))
single_cell_missing_class_num_gather = gather(single_cell_missing_class_num, method, num, 1:3)
#single_cell_missing_class_num_gather$ratio = factor(single_cell_missing_class_num_gather$ratio, levels = c("50%", "60%", "70%", "80%", "90%", "100%"), ordered = TRUE)
single_cell_missing_class_num_gather$ratio = factor(single_cell_missing_class_num_gather$ratio, levels = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"), ordered = TRUE)
single_cell_missing_class_num_gather$small_num = single_cell_missing_class_num_gather$num/1000
single_cell_missing_class_num_full_plot = ggplot(single_cell_missing_class_num_gather, aes(x=ratio, y=small_num, group=method)) +
  geom_line(aes(color=method))+
  geom_point(aes(color=method, shape=method)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 3.6)) +
  scale_color_manual(values = c("#E54D37", "#5CBED3", "#06A088")) +
  ylab("# Proteins (x1000)") +
  xlab("Completeness") +
  theme_light() +
  theme(legend.position=c(0.33,0.15),
        legend.title = element_text(size=4.5, face="bold"),
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(0.1, "cm"),
        legend.background = element_blank()) +
  geom_vline(xintercept = 8, linetype="dashed", color = "black", size=0.2)+
  geom_vline(xintercept = 10, linetype="dashed", color = "red", size=0.2)+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))
single_cell_missing_class_num_full_plot
ggsave("./figures/Figure7a.pdf", single_cell_missing_class_num_full_plot, width=2.5, height = 2, units = c("in"), dpi=400)

# summarize numbers
single_cell_fig3_lib_pro_result_group = average_summarize(single_cell_fig3_lib_pro_result, fig4_lib_exp_one_class, 0)
single_cell_fig4_lib_pro_result_group = average_summarize(single_cell_fig4_lib_pro_result, fig4_lib_exp_one_class, 0)
single_cell_diann_pro_result_group = average_summarize(single_cell_diann_pro_result, fig4_lib_exp_one_class, 0)
mean(single_cell_fig3_lib_pro_result_group$n)
mean(single_cell_fig4_lib_pro_result_group$n)
mean(single_cell_diann_pro_result_group$n)

single_cell_fig3_lib_pro_result_group_70 = average_summarize(single_cell_fig3_lib_pro_result, fig4_lib_exp_one_class, 0.7)
single_cell_fig4_lib_pro_result_group_70 = average_summarize(single_cell_fig4_lib_pro_result, fig4_lib_exp_one_class, 0.7)
mean(single_cell_fig3_lib_pro_result_group_70$n)
mean(single_cell_fig4_lib_pro_result_group_70$n)

single_cell_fig3_lib_pro_result_group_90 = average_summarize(single_cell_fig3_lib_pro_result, fig4_lib_exp_one_class, 0.9)
single_cell_fig4_lib_pro_result_group_90 = average_summarize(single_cell_fig4_lib_pro_result, fig4_lib_exp_one_class, 0.9)
mean(single_cell_fig3_lib_pro_result_group_90$n)
mean(single_cell_fig4_lib_pro_result_group_90$n)

single_cell_70_non_miss_over = plot(venn(list(fig3_lib = missing_summarize_one_class(single_cell_fig3_lib_pro_result, fig4_lib_exp_one_class, 0.7)$V1, 
                fig4_lib = missing_summarize_one_class(single_cell_fig4_lib_pro_result, fig4_lib_exp_one_class, 0.7)$V1,
                paper = missing_summarize_one_class(single_cell_diann_pro_result, fig4_lib_exp_one_class, 0.7)$V1)), 
                labels = list(labels=c("FragPipe\nhigh-input\nLib", "FragPipe", "Anuar\n et al."), fontsize=4.5),
                quantities = list(fontsize = 5),
                fills = c("#2270B5", "#9ECAE1", "#6BAED6"))
single_cell_70_non_miss_over
ggsave("./figures/Figure7b.pdf", single_cell_70_non_miss_over, width=1.5, height = 1.5, units = c("in"), dpi=400)

# PCA.
plot_pca = function(imputed_data, annotation_data, padding){
  imputed_data_gather = gather(imputed_data, sample_name, value, 2:149)
  imputed_data_gather_spread = spread(imputed_data_gather, ProteinID, value)
  annotation_data_simple = annotation_data %>%
    select(sample_name, condition)
  annotation_data_simple$sample_name = paste(padding, annotation_data_simple$sample_name, "_imputed_intensity", sep="")
  imputed_data_gather_spread_cond = left_join(imputed_data_gather_spread, annotation_data_simple, by="sample_name")
  imputed_data_gather_spread_cond$sample_name = NULL
  imputed_data_gather_spread_cond = imputed_data_gather_spread_cond[order(imputed_data_gather_spread_cond$condition),]
  print(imputed_data_gather_spread_cond$condition)
  imputed_data_gather_spread_cond_pca = select(imputed_data_gather_spread_cond, -condition)
  pc_comp <- prcomp(imputed_data_gather_spread_cond_pca,
                    center = T,
                    scale. = F)
  components <- pc_comp[["x"]]
  components <- data.frame(components)
  components$PC2 <- -components$PC2
  components$PC3 <- -components$PC3
  components = cbind(components, imputed_data_gather_spread_cond$condition)
  
  fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~imputed_data_gather_spread_cond$condition, 
                 colors = c('#636EFA','#EF553B','#00CC96', "blue"),
                 marker = list(size = 9, line = list(color = 'black',
                                                     width = 1.5))) %>%
    add_markers(size=0)
  
  
  fig <- fig %>%
    layout(scene = list(bgcolor = "white", 
                        xaxis = list(title = list(text=paste('PC1=', summary(pc_comp)$importance[2,][1][[1]]*100, "%", sep=""), font=list(size=15)), 
                                     tickfont = list(size = 6)),
                        yaxis = list(title = list(text=paste('PC2=', summary(pc_comp)$importance[2,][2][[1]]*100, "%", sep=""), font=list(size=15)), 
                                     tickfont = list(size = 6)),
                        zaxis = list(title = list(text=paste('PC3=', summary(pc_comp)$importance[2,][3][[1]]*100, "%", sep=""), font=list(size=15)), 
                                     tickfont = list(size = 6)),
                        camera = list(eye = list(x = -1, y = -2.4, z = 1.5))))
  
  return(fig)
}
fig4_lib_exp_over_anno = fread("./lowInputData/Direct_lib_result/experiment_annotation_over.tsv")

fig3_lib_imputed_table = fread("./lowInputData/Tonsil_lib_result/fragpipe_analyst/Imputed_matrix.csv")
fig3_lib_pca = plot_pca(fig3_lib_imputed_table, fig4_lib_exp_over_anno, "30")
fig3_lib_pca

### Imputed_matrix.csv download from FragPipe-analysist.
fig4_lib_imputed_table = fread("./lowInputData/Direct_lib_result/fragpipe_analyst/Imputed_matrix.csv")
single_cell_fig4_lib_gg_result = read_maxlfq("./lowInputData/Direct_lib_result/diann-output/gene_maxlfq.tsv")
single_cell_fig4_lib_pro_result_70_non = missing_summarize_one_class(single_cell_fig4_lib_gg_result, fig4_lib_exp_one_class, 0.7)
fig4_lib_imputed_table_passed = fig4_lib_imputed_table %>%
  filter(ProteinID %in% single_cell_fig4_lib_pro_result_70_non$V1)

fig4_lib_pca = plot_pca(fig4_lib_imputed_table_passed, fig4_lib_exp_over_anno, "")
fig4_lib_pca # PCA plot save manually since need to modify 3D position.
reticulate::import("sys")
reticulate::import("plotly")
save_image(fig4_lib_pca, "./figures/pca.pdf", scale=2)


# gene expression non_imputed
get_expre_one_gene = function(report_data, class_file, gene_name_list){
  col_num = ncol(report_data)
  report_data_used = report_data %>%
    filter(V1 %in% gene_name_list)
  report_data_gather = gather(report_data_used, sample_name, quant, (col_num-147):col_num)
  
  report_data_gather_con = inner_join(report_data_gather, class_file, by="sample_name")
  
  return(report_data_gather_con)
}
get_expre_one_gene_imputed = function(report_data, class_file, gene_name_list){
  col_num = ncol(report_data)
  report_data_used = report_data %>%
    filter(ProteinID %in% gene_name_list)
  report_data_gather = gather(report_data_used, sample_name, quant, (col_num-147):col_num)
  
  report_data_gather_con = inner_join(report_data_gather, class_file, by="sample_name")
  
  return(report_data_gather_con)
}
genes_needed =c("CD19", "CD3D", "CDH1")
#genes_needed =c("IL16", "IL18", "STAT1", "CD19", "CD3D", "CDH1") All six genes have good correlation with original styudy. Only select three to show in main figure
fig4_lib_exp_one_class = fread("./lowInputData/Direct_lib_result/experiment_annotation_over.tsv") %>%
  select(sample_name, condition)
single_cell_fig4_lib_gg_result = read_maxlfq("./lowInputData/Direct_lib_result/diann-output/gene_maxlfq.tsv")
fig4_lib_exp_one_class$sample_name = paste("Bluto_230", fig4_lib_exp_one_class$sample_name, sep = "")
genes_needed_data = get_expre_one_gene(single_cell_fig4_lib_gg_result, fig4_lib_exp_one_class, genes_needed)
genes_needed_data$quant_log = log2(genes_needed_data$quant)

genes_needed_data_plot = ggplot(genes_needed_data, aes(x=condition, y=quant_log, color=condition)) +
  geom_boxplot(outlier.size = 0, size=0.1) +
  geom_point(position=position_jitterdodge(0.1), size=0.1, alpha=0.6)+
  scale_x_discrete(labels=c('Ep', 'GC', 'MZ', 'T cell')) + 
  scale_color_manual(values = c("#E54D37", "#5CBED3", "#06A088", "#3C5587")) +
  ylab("Protein Abundance (log2)") +
  xlab("Region") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "none") +
  theme(strip.text = element_text(size = 5, color = "black", margin = margin(b=1)),
        strip.background = element_blank())+
  facet_wrap(~ V1, ncol=1)
genes_needed_data_plot
ggsave("./figures/Figure7d.pdf", genes_needed_data_plot, width=1.2, height = 3, units = c("in"), dpi=400)


# Volcano plot
fig4_lib_de_table = fread("./lowInputData/Direct_lib_result/fragpipe_analyst/DE_results.csv")
fig4_lib_de_table_tc_mz = fig4_lib_de_table %>%
  select(`Gene Name`, `T_cell_vs_mantle_zone_log2 fold change`, `T_cell_vs_mantle_zone_p.adj`)
colnames(fig4_lib_de_table_tc_mz) = c("gene", "fold_change", "p_adj")
fig4_lib_de_table_tc_mz$label_name = ifelse((abs(fig4_lib_de_table_tc_mz$fold_change)>=1.5 & -log10(fig4_lib_de_table_tc_mz$p_adj) >20), fig4_lib_de_table_tc_mz$gene, NA)
fig4_lib_de_table_tc_mz$sig = ifelse((abs(fig4_lib_de_table_tc_mz$fold_change)>=0.6 & fig4_lib_de_table_tc_mz$p_adj <0.05), "sig", "nsg")
tc_mz_volcano = ggplot(fig4_lib_de_table_tc_mz, aes(x=fold_change, y=-log10(p_adj), color=sig)) +
  geom_point(size=0.2, alpha=0.8) +
  scale_color_manual(values=c("#999999", "black")) +
  ylab("Adjusted P-value (-log10)") +
  xlab("log2 Fold Change") +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_text_repel(label = fig4_lib_de_table_tc_mz$label_name, max.overlaps = Inf, color="black", size=1.5, segment.color = 'transparent', box.padding = 0.05) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "none")
tc_mz_volcano
ggsave("./figures/Figure7e.pdf", tc_mz_volcano, width=5.5, height = 3.5, units = c("in"), dpi=400)


