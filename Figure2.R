library(tidyverse)
library(data.table)

########################
### ID numbers plot
########################
count_id_num = function(report_data, miss_num, col_name, method_name, level_name){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-15):col_num)
  
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
  
  report_data_filter_gather = gather(report_data_filter, raw_files, quant, (col_num-15):col_num)
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
### Missing values estimate
missing_summarize = function(report_data, method){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-15):col_num)
  
  report_data_group = report_data_gather %>%
    filter(!is.na(quant)) %>%
    group_by(V1) %>%
    summarise(n=n()) %>%
    ungroup()
  
  report_data_group_out_25 = report_data_group %>%
    filter(n >= 4)
  report_data_group_out_50 = report_data_group %>%
    filter(n >= 8)
  report_data_group_out_100 = report_data_group %>%
    filter(n == 16)
  
  report_data_group_out_25_pure = report_data_group_out_25 %>%
    filter(!V1 %in% report_data_group_out_50$V1)
  report_data_group_out_50_pure = report_data_group_out_50 %>%
    filter(!V1 %in% report_data_group_out_100$V1)
  report_data_group_out_other = report_data_group %>%
    filter(!V1 %in% report_data_group_out_25$V1)
  
  id_portion = c("100%", ">50%", ">25%", ">0%")
  id_num = c(nrow(report_data_group_out_100), nrow(report_data_group_out_50_pure), nrow(report_data_group_out_25_pure), nrow(report_data_group_out_other))
  
  out_data = data_frame(id_portion, id_num)
  out_data$method = method
  
  return(out_data)
}
missing_summarize_sp_pro = function(spectronaut_report, method, col_name){
  col_name <- enquo(col_name)
  spectronaut_report_group = spectronaut_report %>%
    distinct(!!col_name, R.FileName) %>%
    group_by(!!col_name) %>%
    summarise(n=n()) %>%
    ungroup()
  spectronaut_report_group_out_25 = spectronaut_report_group %>%
    filter(n >= 4)
  spectronaut_report_group_out_50 = spectronaut_report_group %>%
    filter(n >= 8)
  spectronaut_report_group_out_100 = spectronaut_report_group %>%
    filter(n == 16)
  
  spectronaut_report_group_out_25_pure = spectronaut_report_group_out_25 %>%
    filter(!PG.GroupLabel %in% spectronaut_report_group_out_50$PG.GroupLabel)
  spectronaut_report_group_out_50_pure = spectronaut_report_group_out_50 %>%
    filter(!PG.GroupLabel %in% spectronaut_report_group_out_100$PG.GroupLabel)
  report_data_group_out_other = spectronaut_report_group %>%
    filter(!PG.GroupLabel %in% spectronaut_report_group_out_25$PG.GroupLabel)
  
  id_portion = c("100%", ">50%", ">25%", ">0%")
  id_num = c(nrow(spectronaut_report_group_out_100), nrow(spectronaut_report_group_out_50_pure), 
             nrow(spectronaut_report_group_out_25_pure), nrow(report_data_group_out_other))
  
  out_data = data_frame(id_portion, id_num)
  out_data$method = method
  
  return(out_data)
}
get_log_quant = function(diann_report, method){
  diann_report_gather = gather(diann_report, key="run", value="quanty", 2:17)
  diann_report_gather$log_quanty = log(diann_report_gather$quanty)
  diann_report_gather$method = method
  return(diann_report_gather)
}
read_maxlfq = function(file_path){
  out_data = fread(file_path) %>%
    filter(V1 != "")
}

spectronaut_report_num = fread("./revisionData/TNBCData/spectronaut/TNBC2020_4223_Spectronaut18-5_library_free_directDIA_IdentificationsOverview.tsv")
spectronaut_report = fread("./revisionData/TNBCData/spectronaut/20240921_233058_TNBC2020_4223_Spectronaut18-5_library_free_directDIA_ReportAll.tsv")
spectronaut_report_lib_based_num = fread("./revisionData/TNBCData/spectronaut/TNBC2020_4223_Spectronaut18-5_library_based_IdentificationsOverview.tsv")
spectronaut_report_lib_based = fread("./revisionData/TNBCData/spectronaut/20241002_140130_TNBC2020_4223_Spectronaut18-5_library_based_Report.tsv")


breast_cancer_diatracer_pg_processed = read_maxlfq("./revisionData/TNBCData/diaTracer_result_Frag22_methylthiolation/diann-output/protein_maxlfq.tsv")
breast_cancer_hybrid_pg_processed = read_maxlfq("./revisionData/TNBCData/diaTracer_hybrid_result_Frag22_methylthiolation/diann-output/protein_maxlfq.tsv")
breast_cancer_diann_pg_processed = read_maxlfq("./revisionData/TNBCData/14_MainSearch_4223_DIA-NN_1.8.1_library_free/14_MainSearch_4223_DIA-NN_1.8.1_library_free/protein_maxlfq.tsv")
breast_cancer_diann_pg_lib_based_processed = read_maxlfq("./revisionData/TNBCData/12_MainSearch_4223_DIA-NN_1.8.1_library_based/12_MainSearch_4223_DIA-NN_1.8.1_newSN16lib/protein_maxlfq.tsv")

breast_cancer_diatracer_pg_processed_num = count_id_num(breast_cancer_diatracer_pg_processed, 0, V1, "FragPipe ", "Tryptic")
breast_cancer_hybrid_pg_processed_num = count_id_num(breast_cancer_hybrid_pg_processed, 0, V1, "FragPipe", "Tryptic")
breast_cancer_diann_pg_processed_num = count_id_num(breast_cancer_diann_pg_processed, 0, V1, "DIA-NN\nlib-free", "Tryptic")

spectronaut_report_protein_num = spectronaut_report_num %>%
  select(FileName, Proteins)
colnames(spectronaut_report_protein_num) = c("V1", "n")
spectronaut_report_protein_num$method = "Spectronaut "
spectronaut_report_protein_num$Cleavage = "Tryptic"

spectronaut_report_lib_based_protein_num = spectronaut_report_lib_based_num %>%
  select(FileName, Proteins)
colnames(spectronaut_report_lib_based_protein_num) = c("V1", "n")
spectronaut_report_lib_based_protein_num$method = "Spectronaut"
spectronaut_report_lib_based_protein_num$Cleavage = "Tryptic"

breast_cancer_diatracer_pg_processed_num$type = "Direct DIA"
breast_cancer_diann_pg_processed_num$type = "Direct DIA"
spectronaut_report_protein_num$type = "Direct DIA"
breast_cancer_hybrid_pg_processed_num$type = "DDA/DIA hybrid library"
spectronaut_report_lib_based_protein_num$type = "DDA/DIA hybrid library"
breast_cancer_nums_for_plot_protein = bind_rows(breast_cancer_diatracer_pg_processed_num) %>%
  bind_rows(breast_cancer_hybrid_pg_processed_num) %>%
  bind_rows(breast_cancer_diann_pg_processed_num) %>%
  bind_rows(spectronaut_report_protein_num) %>%
  bind_rows(spectronaut_report_lib_based_protein_num)
breast_cancer_nums_for_plot_protein$small_n = breast_cancer_nums_for_plot_protein$n / 1000
breast_cancer_nums_for_plot_protein$method <- factor(breast_cancer_nums_for_plot_protein$method , levels=c("Spectronaut ", "DIA-NN\nlib-free", "Spectronaut",  "FragPipe ", "FragPipe"))
breast_cancer_nums_for_plot_protein$type = factor(breast_cancer_nums_for_plot_protein$type, levels=c('Direct DIA','DDA/DIA hybrid library'))
breast_cancer_pg_dis_plot = ggplot(data=breast_cancer_nums_for_plot_protein, aes(x=method, y=small_n, col=method)) +
  geom_boxplot(outlier.size = 0, size=0.1) +
  geom_point(position=position_jitterdodge(0.1), size=0.3, alpha=0.6)+
  scale_y_continuous( expand = c(0, 0), limits = c(5, 10.3)) +
  scale_color_brewer(palette="Dark2") +
  ylab("# Proteins (x1000)") +
  xlab("Method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "none",
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=6, colour = 'black')) +
  facet_grid(~type, scales="free", space = "free")
breast_cancer_pg_dis_plot

breast_cancer_diatracer_pg_processed_miss = missing_summarize(breast_cancer_diatracer_pg_processed, "FragPipe ")
breast_cancer_hybrid_pg_processed_miss = missing_summarize(breast_cancer_hybrid_pg_processed, "FragPipe")
breast_cancer_diann_pg_processed_miss = missing_summarize(breast_cancer_diann_pg_processed, "DIA-NN\nlib-free")
breast_cancer_spectronaut_pro_miss = missing_summarize_sp_pro(spectronaut_report, "Spectronaut ", PG.GroupLabel)
breast_cancer_spectronaut_lib_based_pro_miss = missing_summarize_sp_pro(spectronaut_report_lib_based, "Spectronaut", PG.GroupLabel)

breast_cancer_diatracer_pg_processed_miss$type = "Direct DIA"
breast_cancer_diann_pg_processed_miss$type = "Direct DIA"
breast_cancer_spectronaut_pro_miss$type = "Direct DIA"
breast_cancer_hybrid_pg_processed_miss$type = "DDA/DIA hybrid library"
breast_cancer_spectronaut_lib_based_pro_miss$type = "DDA/DIA hybrid library"
breast_cancer_num_pro_miss = bind_rows(breast_cancer_diatracer_pg_processed_miss) %>%
  bind_rows(breast_cancer_hybrid_pg_processed_miss) %>%
  bind_rows(breast_cancer_diann_pg_processed_miss) %>%
  bind_rows(breast_cancer_spectronaut_pro_miss) %>%
  bind_rows(breast_cancer_spectronaut_lib_based_pro_miss)
breast_cancer_num_pro_miss$small_id_num = breast_cancer_num_pro_miss$id_num/1000
breast_cancer_num_pro_miss$id_portion = factor(breast_cancer_num_pro_miss$id_portion, levels = c(">0%", ">25%", ">50%", "100%" ), ordered = TRUE)
breast_cancer_num_pro_miss$method = factor(breast_cancer_num_pro_miss$method, levels=c("Spectronaut ", "DIA-NN\nlib-free", "FragPipe ", "Spectronaut", "FragPipe"), ordered = TRUE)
breast_cancer_num_pro_miss$type = factor(breast_cancer_num_pro_miss$type, levels=c('Direct DIA','DDA/DIA hybrid library'))
breast_cancer_num_pro_miss_plot = ggplot(breast_cancer_num_pro_miss, aes(x=method, y=small_id_num, fill= id_portion)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) +
  scale_fill_brewer(palette="Blues", name="Completeness") +
  scale_y_continuous(expand = c(0.01, 0)) +
  ylab("# Proteins (x1000)") +
  xlab("Method") +
  theme_light() +
  labs(fill= "Non-missing\nvalue filter") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "right",
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=6, colour = 'black')) +
  facet_grid(~type, scales="free", space = "free")
breast_cancer_num_pro_miss_plot

breast_cancer_plot_all= ggarrange(breast_cancer_pg_dis_plot, breast_cancer_num_pro_miss_plot,
                       ncol = 1, nrow = 2, align="v", labels = c("a", "b"), font.label = list(size = 10))
ggsave("./figures/Figure2.pdf", breast_cancer_plot_all, width=4, height = 3, units = c("in"), dpi=400)

