library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(ggrepel)
library(openxlsx)

########################
### ID numbers plot
#######################

### Missing values estimate
missing_summarize = function(report_data, method){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-39):col_num)
  
  report_data_group = report_data_gather %>%
    filter(!is.na(quant)) %>%
    group_by(V1) %>%
    summarise(n=n()) %>%
    ungroup()
  
  report_data_group_out_25 = report_data_group %>%
    filter(n >= 10)
  report_data_group_out_50 = report_data_group %>%
    filter(n >= 20)
  report_data_group_out_100 = report_data_group %>%
    filter(n == 40)
  
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
read_maxlfq = function(file_path){
  out_data = fread(file_path) %>%
    filter(V1 != "")
}

# a. missing value comp with paper
plasma_diann_try_pro_processed = read_maxlfq("./plasmaData/DIA_NN_LF_result/protein_maxlfq.tsv")
plasma_dia_try_pro_processed = read_maxlfq("./plasmaData/DIA_lib_try_result/diann-output/protein_maxlfq.tsv")
plasma_dia_semi_pro_processed = read_maxlfq("./plasmaData/DIA_lib_semi_result/diann-output/protein_maxlfq.tsv")

plasma_diann_try_pro_processed_miss = missing_summarize(plasma_diann_try_pro_processed, "DIA-NN\nlib-free(Tr)")
plasma_dia_try_pro_processed_miss = missing_summarize(plasma_dia_try_pro_processed, "FP-diaTracer\n(Tr)")
plasma_dia_semi_pro_processed_miss = missing_summarize(plasma_dia_semi_pro_processed, "FP-diaTracer\n(Semi)")
plasma_num_pro_miss = bind_rows(plasma_diann_try_pro_processed_miss) %>%
  bind_rows(plasma_dia_try_pro_processed_miss) %>%
  bind_rows(plasma_dia_semi_pro_processed_miss)

plasma_num_pro_miss$small_id_num = plasma_num_pro_miss$id_num/1000
plasma_num_pro_miss$id_portion = factor(plasma_num_pro_miss$id_portion, levels = c("<25%", ">25%", ">50%", "100%" ), ordered = TRUE)
plasma_num_pro_miss$method = factor(plasma_num_pro_miss$method, levels = c("FP-diaTracer\n(Tr)", "FP-diaTracer\n(Semi)", "DIA-NN\nlib-free(Tr)"), ordered = TRUE)
plasma_num_pro_miss_plot = ggplot(plasma_num_pro_miss, aes(x=method, y=small_id_num, fill= id_portion)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) + 
  scale_fill_brewer(palette="Blues") +
  scale_y_continuous(expand = c(0.01, 0)) +
  ylab("# Proteins (x1000)") +
  xlab("Method") +
  theme_light() +
  labs(fill= "Non-missing\nvalue filter") +
  #theme(legend.position="none") +
  #annotate(geom="text", x=0.8, y=4, label="x 1000", size = 2) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "right",
        legend.title = element_text(size=4, face="bold"),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.15, "cm"),
        legend.margin = margin(-10, 2, 0, -10))
plasma_num_pro_miss_plot
plasma_num_pro_miss_sum = plasma_num_pro_miss %>%
  group_by(method) %>%
  summarise(n=sum(id_num))

plasma_dia_try_pre_processed = read_maxlfq("./plasmaData/DIA_lib_try_result/diann-output/precursor_maxlfq.tsv")
plasma_dia_semi_pre_processed = read_maxlfq("./plasmaData/DIA_lib_semi_result/diann-output/precursor_maxlfq.tsv")
plasma_diann_try_pre_processed = read_maxlfq("./plasmaData/DIA_NN_LF_result/precursor_maxlfq.tsv")

plasma_diann_try_pre_processed_miss = missing_summarize(plasma_diann_try_pre_processed, "DIA-NN\nlib-free(Tr)")
plasma_dia_try_pre_processed_miss = missing_summarize(plasma_dia_try_pre_processed, "FP-diaTracer\n(Tr)")
plasma_dia_semi_pre_processed_miss = missing_summarize(plasma_dia_semi_pre_processed, "FP-diaTracer\n(Semi)")
plasma_num_pre_miss = 
  bind_rows(plasma_diann_try_pre_processed_miss) %>%
  bind_rows(plasma_dia_try_pre_processed_miss) %>%
  bind_rows(plasma_dia_semi_pre_processed_miss)

plasma_num_pre_miss$small_id_num = plasma_num_pre_miss$id_num/1000
plasma_num_pre_miss$id_portion = factor(plasma_num_pre_miss$id_portion, levels = c("<25%", ">25%", ">50%", "100%" ), ordered = TRUE)
plasma_num_pre_miss$method = factor(plasma_num_pre_miss$method, levels = c("FP-diaTracer\n(Tr)", "FP-diaTracer\n(Semi)", "DIA-NN\nlib-free(Tr)"), ordered = TRUE)
plasma_num_pre_miss_plot = ggplot(plasma_num_pre_miss, aes(x=method, y=small_id_num, fill= id_portion)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) + 
  scale_fill_brewer(palette="Blues", name="Portion") +
  scale_y_continuous(expand = c(0.01, 0)) +
  ylab("# Precursors (x1000)") +
  xlab("Method") +
  theme_light() +
  labs(fill= "Non-missing\nvalue filter") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "right",
        legend.title = element_text(size=4, face="bold"),
        legend.text = element_text(size = 4),
        legend.key.size = unit(0.15, "cm"),
        legend.margin = margin(-10, 2, 0, -10))
plasma_num_pre_miss_plot
plasma_num_pre_miss_sum = plasma_num_pre_miss %>%
  group_by(method) %>%
  summarise(n=sum(id_num))

# gene expression
plasma_try_de_result = fread("./plasmaData/DIA_lib_try_result/frapipe_analyst/DE_results.csv")
plasma_semi_de_result = fread("./plasmaData/DIA_lib_semi_result/fragpipe_analyst/DE_results.csv")
plasma_de_result = list("Tryptic" = plasma_try_de_result, "Semi-Tryptic" = plasma_semi_de_result)
write.xlsx(plasma_de_result, file = './supplements/table_s3.xls')

plasma_try_de_result_sig = plasma_try_de_result %>%
  filter(significant == TRUE)
plasma_semi_de_result_sig = plasma_semi_de_result %>%
  filter(significant == TRUE)

plasma_semi_de_result_sig_unique = plasma_semi_de_result_sig %>%
  filter(!`Gene Name` %in% plasma_try_de_result_sig$`Gene Name`)

plasma_semi_de_result_cc = plasma_semi_de_result %>%
  select(`Gene Name`, `cancer_vs_control_log2 fold change`, cancer_vs_control_p.adj)
colnames(plasma_semi_de_result_cc) = c("gene", "fold_change", "p_adj")
plasma_semi_de_result_cc$label_name = ifelse((abs(plasma_semi_de_result_cc$fold_change)>=1 & -log10(plasma_semi_de_result_cc$p_adj) >1.6), plasma_semi_de_result_cc$gene, NA)
plasma_semi_de_result_cc$sig = ifelse((abs(plasma_semi_de_result_cc$fold_change)>=0.6 & plasma_semi_de_result_cc$p_adj <0.05), "sig", "nsg")
plasma_semi_de_result_cc_volcano = ggplot(plasma_semi_de_result_cc, aes(x=fold_change, y=-log10(p_adj), color=sig)) +
  geom_point(size=0.2) +
  ylab("Adjusted P-value (-log10)") +
  xlab("log2 Fold Change") +
  scale_color_manual(values=c("#999999", "black")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_text_repel(label = plasma_semi_de_result_cc$label_name, max.overlaps = Inf, color="black", size=1.5, segment.color = 'transparent', box.padding = 0.05) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "none")
plasma_semi_de_result_cc_volcano

plasma_semi_de_result_up = plasma_semi_de_result %>%
  filter(significant == TRUE) %>%
  filter(`cancer_vs_control_log2 fold change` > 0)
plasma_semi_de_result_down = plasma_semi_de_result %>%
  filter(significant == TRUE) %>%
  filter(`cancer_vs_control_log2 fold change` < 0)

# feature
plasma_annotation_data = fread("./plasmaData/experiment_annotation.txt") %>%
  select(sample_name, condition)
feature_data_try = plasma_dia_try_pro_processed %>%
  filter(V1 %in% c("P16104", "P31751"))
feature_data_try$V1 = str_replace_all(feature_data_try$V1, "P31751", "AKT2")
feature_data_try$V1 = str_replace_all(feature_data_try$V1, "P16104", "H2AX")
feature_data_try_gather = gather(feature_data_try, sample_name, quant, 2:41)
feature_data_try_gather_con = inner_join(feature_data_try_gather, plasma_annotation_data, by="sample_name")
feature_data_try_gather_con$quant_log = log2(feature_data_try_gather_con$quant)
feature_data_try_gather_con$type = "Tryptic"

feature_data_semi = plasma_dia_semi_pro_processed %>%
  filter(V1 %in% c("P16104", "P31751"))
feature_data_semi$V1 = str_replace_all(feature_data_semi$V1, "P31751", "AKT2")
feature_data_semi$V1 = str_replace_all(feature_data_semi$V1, "P16104", "H2AX")
feature_data_semi_gather = gather(feature_data_semi, sample_name, quant, 2:41)
feature_data_semi_gather_con = inner_join(feature_data_semi_gather, plasma_annotation_data, by="sample_name")
feature_data_semi_gather_con$quant_log = log2(feature_data_semi_gather_con$quant)
feature_data_semi_gather_con$type = "Semi"
feature_data_all = bind_rows(feature_data_try_gather_con, feature_data_semi_gather_con)

plasma_feature_plot = ggplot(feature_data_all, aes(x=condition, y=quant_log, color=condition)) +
  geom_boxplot(outlier.size = 0, size=0.1) +
  geom_point(position=position_jitterdodge(0.1), size=0.1, alpha=0.6)+
  scale_x_discrete(labels=c('NSCLC', 'Control')) + 
  scale_color_manual(values = c("#E54D37", "#5CBED3")) +
  stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5 , label.y = 13, size=2)+
  ylab("Protein Abun (log2)") +
  xlab("Group") +
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
  facet_grid(V1~type)
plasma_feature_plot


plasma_plot = ggarrange(ggarrange(plasma_num_pro_miss_plot, plasma_num_pre_miss_plot, widths = c(1.5, 1.5),
                                ncol = 2, nrow = 1, align="h", labels = c("a", "b"), font.label = list(size = 10), 
                                common.legend = T, legend = "right"),
                      ggarrange(plasma_semi_de_result_cc_volcano, plasma_feature_plot, widths = c(2, 1.5),
                                ncol = 2, nrow = 1, align="h", labels = c("c", "d"), font.label = list(size = 10)),
                      nrow = 2, ncol=1, heights = c(1.3,2.2))
plasma_plot
ggsave("./figures/Figure3.pdf", plasma_plot, width=5.3, height = 4, units = c("in"), dpi=400)



