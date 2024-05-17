library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(ggpubr)
library(ggplot2)
library(ggrepel)

########################
### ID numbers plot
########################
count_id_num = function(report_data, miss_num, col_name, method_name, level_name){
  col_num = ncol(report_data)
  report_data_gather = gather(report_data, raw_files, quant, (col_num-33):col_num)
  
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
  
  report_data_filter_gather = gather(report_data_filter, raw_files, quant, (col_num-33):col_num)
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
read_maxlfq = function(file_path){
  out_data = fread(file_path) %>%
    filter(V1 != "")
}

# a gg num
csf_diann_try_gg_processed = read_maxlfq("./CSFData/DIA_NN_LF_result/gene_maxlfq.tsv")
csf_dda_try_gg_processed = read_maxlfq("./CSFData/DDA_lib_try_result/diann-output/gene_maxlfq.tsv")
csf_dia_try_gg_processed = read_maxlfq("./CSFData/DIA_lib_try_result/diann-output/gene_maxlfq.tsv")
csf_dda_dia_semi_gg_processed = read_maxlfq("./CSFData/DDA_DIA_lib_semi_result/diann-output/gene_maxlfq.tsv")
csf_dda_semi_gg_processed = read_maxlfq("./CSFData/DDA_lib_semi_result/diann-output/gene_maxlfq.tsv")
csf_dia_semi_gg_processed = read_maxlfq("./CSFData/DIA_lib_semi_result/diann-output/gene_maxlfq.tsv")
csf_dda_dia_try_gg_processed = read_maxlfq("./CSFData/DDA_DIA_lib_try_result/diann-output/gene_maxlfq.tsv")

csf_dda_dia_try_gg_processed_num = count_id_num(csf_dda_dia_try_gg_processed, 0, V1, "FP-diaTracer\nhybrid", "Tryptic")
csf_diann_try_gg_processed_num = count_id_num(csf_diann_try_gg_processed, 0, V1, "DIA-NN\nlib-free", "Tryptic")
csf_dda_try_gg_processed_num = count_id_num(csf_dda_try_gg_processed, 0, V1, "FP-\nDDALib", "Tryptic")
csf_dia_try_gg_processed_num = count_id_num(csf_dia_try_gg_processed, 0, V1, "FP-\ndiaTracer", "Tryptic")
csf_dda_dia_semi_gg_processed_num = count_id_num(csf_dda_dia_semi_gg_processed, 0, V1, "FP-diaTracer\nhybrid", "Semi-tryptic")
csf_dda_semi_gg_processed_num = count_id_num(csf_dda_semi_gg_processed, 0, V1, "FP-\nDDALib", "Semi-tryptic")
csf_dia_semi_gg_processed_num = count_id_num(csf_dia_semi_gg_processed, 0, V1, "FP-\ndiaTracer", "Semi-tryptic")

csf_nums_for_plot_protein = bind_rows(csf_dda_dia_try_gg_processed_num) %>%
  bind_rows(csf_diann_try_gg_processed_num) %>%
  bind_rows(csf_dda_try_gg_processed_num) %>%
  bind_rows(csf_dia_try_gg_processed_num) %>%
  bind_rows(csf_dda_dia_semi_gg_processed_num) %>%
  bind_rows(csf_dda_semi_gg_processed_num) %>%
  bind_rows(csf_dia_semi_gg_processed_num)
csf_nums_for_plot_protein$small_n = csf_nums_for_plot_protein$n / 1000
csf_nums_for_plot_protein$method <- factor(csf_nums_for_plot_protein$method , levels=c("FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer", "DIA-NN\nlib-free"))
csf_gg_dis_plot = ggplot(csf_nums_for_plot_protein, aes(x=method, y=small_n, col=Cleavage)) +
  geom_boxplot(outlier.size = 0, size=0.1) +
  geom_point(position=position_jitterdodge(0.1), size=0.3, alpha=0.6)+
  scale_y_continuous(limits = c(0,1.7), expand = c(0, 0)) +
  scale_color_brewer(palette="Dark2") +
  ylab("# Protein (x1000)") +
  xlab("Method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.85, .85),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
csf_gg_dis_plot

csf_nums_for_plot_protein_avg = csf_nums_for_plot_protein %>%
  group_by(method, Cleavage) %>%
  summarise(avg = mean(n))

# b. precursor num
csf_diann_try_pr_processed = read_maxlfq("./CSFData/DIA_NN_LF_result/precursor_maxlfq.tsv")
csf_dda_try_pr_processed = read_maxlfq("./CSFData/DDA_lib_try_result/diann-output/precursor_maxlfq.tsv")
csf_dia_try_pr_processed = read_maxlfq("./CSFData/DIA_lib_try_result/diann-output/precursor_maxlfq.tsv")
csf_dda_dia_semi_pr_processed = read_maxlfq("./CSFData/DDA_DIA_lib_semi_result/diann-output/precursor_maxlfq.tsv")
csf_dda_semi_pr_processed = read_maxlfq("./CSFData/DDA_lib_semi_result/diann-output/precursor_maxlfq.tsv")
csf_dia_semi_pr_processed = read_maxlfq("./CSFData/DIA_lib_semi_result/diann-output/precursor_maxlfq.tsv")
csf_dda_dia_try_pr_processed = read_maxlfq("./CSFData/DDA_DIA_lib_try_result/diann-output/precursor_maxlfq.tsv")

csf_dda_dia_try_pr_processed_num = count_id_num(csf_dda_dia_try_pr_processed, 0, V1, "FP-diaTracer\nhybrid", "Tryptic")
csf_diann_try_pr_processed_num = count_id_num(csf_diann_try_pr_processed, 0, V1, "DIA-NN\nlib-free", "Tryptic")
csf_dda_try_pr_processed_num = count_id_num(csf_dda_try_pr_processed, 0, V1, "FP-\nDDALib", "Tryptic")
csf_dia_try_pr_processed_num = count_id_num(csf_dia_try_pr_processed, 0, V1, "FP-\ndiaTracer", "Tryptic")
csf_dda_dia_semi_pr_processed_num = count_id_num(csf_dda_dia_semi_pr_processed, 0, V1, "FP-diaTracer\nhybrid", "Semi-tryptic")
csf_dda_semi_pr_processed_num = count_id_num(csf_dda_semi_pr_processed, 0, V1, "FP-\nDDALib", "Semi-tryptic")
csf_dia_semi_pr_processed_num = count_id_num(csf_dia_semi_pr_processed, 0, V1, "FP-\ndiaTracer", "Semi-tryptic")

csf_nums_for_plot_precursor = bind_rows(csf_dda_dia_try_pr_processed_num) %>%
  bind_rows(csf_diann_try_pr_processed_num) %>%
  bind_rows(csf_dda_try_pr_processed_num) %>%
  bind_rows(csf_dia_try_pr_processed_num) %>%
  bind_rows(csf_dda_dia_semi_pr_processed_num) %>%
  bind_rows(csf_dda_semi_pr_processed_num) %>%
  bind_rows(csf_dia_semi_pr_processed_num)
csf_nums_for_plot_precursor$small_n = csf_nums_for_plot_precursor$n / 1000
csf_nums_for_plot_precursor$method <- factor(csf_nums_for_plot_precursor$method , levels=c("FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer", "DIA-NN\nlib-free"))
csf_pr_dis_plot = ggplot(csf_nums_for_plot_precursor, aes(x=method, y=small_n, col=Cleavage)) +
  geom_boxplot(outlier.size = 0, size=0.1) +
  geom_point(position=position_jitterdodge(0.1), size=0.3, alpha=0.6)+
  scale_y_continuous(limits = c(0,19), expand = c(0, 0)) +
  scale_color_brewer(palette="Dark2") +
  ylab("# Precursor (x1000)") +
  xlab("Method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.85, .85),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
csf_pr_dis_plot
csf_nums_for_plot_precursor_avg = csf_nums_for_plot_precursor %>%
  group_by(method, Cleavage) %>%
  summarise(avg = mean(n))

csf_nums_for_plot_precursor_avg_wt_diann <- csf_nums_for_plot_precursor_avg[-c(7), ]
csf_nums_for_plot_precursor_avg_wt_diann_spread = spread(csf_nums_for_plot_precursor_avg_wt_diann, Cleavage, avg)
csf_nums_for_plot_precursor_avg_wt_diann_spread$diff = (csf_nums_for_plot_precursor_avg_wt_diann_spread$`Semi-tryptic`- csf_nums_for_plot_precursor_avg_wt_diann_spread$Tryptic)/csf_nums_for_plot_precursor_avg_wt_diann_spread$Tryptic
mean(csf_nums_for_plot_precursor_avg_wt_diann_spread$diff)

# supplement Overall protein and precursor num
csf_method_list = c("DIA-NN\nlib-free", "FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer", "FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer")
csf_enzyme_list = c("Tryptic", "Tryptic", "Tryptic", "Tryptic", "Semi-tryptic", "Semi-tryptic", "Semi-tryptic")

csf_protein_overal_num = c(nrow(csf_diann_try_gg_processed), nrow(csf_dda_dia_try_gg_processed),
                           nrow(csf_dda_try_gg_processed), nrow(csf_dia_try_gg_processed),
                           nrow(csf_dda_dia_semi_gg_processed), nrow(csf_dda_semi_gg_processed), 
                           nrow(csf_dia_semi_gg_processed))
csf_protein_overal_num_plot_data = data.frame(csf_method_list, csf_enzyme_list, csf_protein_overal_num)
csf_protein_overal_num_plot_data$csf_method_list <- factor(csf_protein_overal_num_plot_data$csf_method_list , levels=c("FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer", "DIA-NN\nlib-free"))
csf_protein_overal_num_plot = ggplot(csf_protein_overal_num_plot_data, aes(x=csf_method_list, y= csf_protein_overal_num, fill=csf_enzyme_list)) +
  geom_bar(position=position_dodge(preserve = "single"), stat="identity", color="black", size=0.05, width = 0.8) +
  scale_fill_brewer(name = "Cleavage", palette="Dark2") +
  ylab("# Protein") +
  xlab("Method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.85, .85),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
csf_protein_overal_num_plot

csf_precursor_overal_num = c(nrow(csf_diann_try_pr_processed), nrow(csf_dda_dia_try_pr_processed),
                           nrow(csf_dda_try_pr_processed), nrow(csf_dia_try_pr_processed),
                           nrow(csf_dda_dia_semi_pr_processed), nrow(csf_dda_semi_pr_processed), 
                           nrow(csf_dia_semi_pr_processed))
csf_precursor_overal_num_plot_data = data.frame(csf_method_list, csf_enzyme_list, csf_precursor_overal_num)
csf_precursor_overal_num_plot_data$csf_method_list <- factor(csf_precursor_overal_num_plot_data$csf_method_list , levels=c("FP-diaTracer\nhybrid", "FP-\nDDALib", "FP-\ndiaTracer", "DIA-NN\nlib-free"))
csf_precursor_overal_num_plot = ggplot(csf_precursor_overal_num_plot_data, aes(x=csf_method_list, y= csf_precursor_overal_num, fill=csf_enzyme_list)) +
  geom_bar(position=position_dodge(preserve = "single"), stat="identity", color="black", size=0.05, width = 0.8) +
  scale_fill_brewer(name = "Cleavage", palette="Dark2") +
  ylab("# Precursor") +
  xlab("Method") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.85, .85),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
csf_precursor_overal_num_plot
ggsave("./supplements/CSF_overall_proteins.pdf", csf_protein_overal_num_plot, width=4.5, height = 4, units = c("in"), dpi=400)
ggsave("./supplements/CSF_overall_precursors.pdf", csf_precursor_overal_num_plot, width=4.5, height = 4, units = c("in"), dpi=400)


# d. running time
csf_time=c(661, 109.8, 661, 86.9, 1496.29)
csf_process=c("diaTracer", "FragPipe*", "diaTracer", "FragPipe*", "DIA-NN")
csf_method = c("FP-diaTracer\n(Semi)", "FP-diaTracer\n(Semi)", "FP-diaTracer\n(Tr)", 
               "FP-diaTracer\n(Tr)", "DIA-NN\nlib-free(Tr)")
csf_running_time = data_frame(csf_time, csf_process, csf_method)
csf_running_time$csf_process = factor(csf_running_time$csf_process, levels = c("DIA-NN", "FragPipe*", "diaTracer"))
csf_running_time$csf_hour = csf_running_time$csf_time/60
csf_running_time$csf_method = factor(csf_running_time$csf_method, levels = c("FP-diaTracer\n(Tr)", "FP-diaTracer\n(Semi)", "DIA-NN\nlib-free(Tr)"))
csf_running_time_plot = ggplot(csf_running_time, aes(x=csf_method, y=csf_hour, fill= csf_process)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) + 
  scale_fill_manual( name = "Process", values = c("#F39B7FB2", "#4DBBD5B2", "#00A087B2")) +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  ylab("Time (h)") +
  xlab("Method") +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 3.6),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.3, .85),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
  

csf_running_time_plot

# semi analysis USing DIA-lib data
hum_uniprot = fread("./database/uniprotkb_reviewed_true_AND_model_organ_2024_03_18.tsv")
hum_uniprot_sig = hum_uniprot %>%
  filter(`Signal peptide` != "")
hum_uniprot_sig_sep = separate(hum_uniprot_sig, `Signal peptide`, into = c("pos"), sep = ";", remove = T)
hum_uniprot_sig_sep_pos = separate(hum_uniprot_sig_sep, pos, into = c("name", "stop_pos"), sep="\\..")

csf_dia_lib_semi_psm = fread("./CSFData/DIA_lib_semi_result/psm.tsv")
csf_dia_lib_semi_psm_semi_pep_nterm = csf_dia_lib_semi_psm %>%
  filter(`Number of Enzymatic Termini` == 1) %>%
  filter(!`Prev AA` %in% c("K", "R", "-")) %>%
  distinct(Peptide, `Modified Peptide`, `Protein Start`, `Protein ID`)
csf_dia_lib_semi_psm_semi_pep_cterm = csf_dia_lib_semi_psm %>%
  filter(`Number of Enzymatic Termini` == 1) %>%
  filter(!`Next AA` %in% c("-")) %>%
  filter(!str_ends(Peptide, "K")) %>%
  filter(!str_ends(Peptide, "R")) %>%
  distinct(Peptide, `Modified Peptide`, `Protein Start`, `Protein ID`)

csf_dia_lib_semi_psm_semi_pep_allterm = rbind(csf_dia_lib_semi_psm_semi_pep_nterm, csf_dia_lib_semi_psm_semi_pep_cterm)
csf_dia_lib_semi_psm_semi_pep_allterm_num = csf_dia_lib_semi_psm_semi_pep_allterm %>%
  group_by(`Protein ID`) %>%
  summarise(semi_num=n())

### Semi/Try ratio
csf_dia_lib_tryp = csf_dia_lib_semi_psm %>%
  filter(`Prev AA` %in% c("K", "R", "-")) %>%
  filter(str_ends(Peptide, "K") | str_ends(Peptide, "R") | `Next AA` %in% c("-"))%>%
  distinct(Peptide, `Modified Peptide`, `Protein Start`, `Protein ID`)
csf_dia_lib_tryp_num = csf_dia_lib_tryp %>%
  group_by(`Protein ID`) %>%
  summarise(tryp_num=n())

csf_dia_lib_tryp_semi_num = left_join(csf_dia_lib_semi_psm_semi_pep_allterm_num, csf_dia_lib_tryp_num, by=c("Protein ID"))
csf_dia_lib_tryp_semi_num$ratio = csf_dia_lib_tryp_semi_num$semi_num/csf_dia_lib_tryp_semi_num$tryp_num
colnames(csf_dia_lib_tryp_semi_num)[1] = "Entry"
csf_dia_lib_tryp_semi_num_out = left_join(csf_dia_lib_tryp_semi_num, hum_uniprot, by="Entry")
write.table(csf_dia_lib_tryp_semi_num_out, "./supplements/table_s1.csv", row.names = F, quote = F, sep = "\t")

### Signal peptides
csf_dia_lib_semi_psm_semi_pep_allterm$preStart = csf_dia_lib_semi_psm_semi_pep_allterm$`Protein Start` - 1
csf_dia_lib_semi_psm_semi_pep_allterm_for_map = csf_dia_lib_semi_psm_semi_pep_allterm %>%
  distinct(`Protein ID`, preStart)
colnames(csf_dia_lib_semi_psm_semi_pep_allterm_for_map)  =c("Entry", "stop_pos")
hum_uniprot_sig_sep_pos$stop_pos = as.numeric(hum_uniprot_sig_sep_pos$stop_pos)
### csf_signal_result is the all proteins identidied after signal peptide
csf_signal_result = inner_join(hum_uniprot_sig_sep_pos, csf_dia_lib_semi_psm_semi_pep_allterm_for_map, by=c("Entry", "stop_pos"))
csf_dia_lib_semi_psm_semi_pep_sig = csf_dia_lib_semi_psm_semi_pep_allterm %>%
  filter(`Protein ID` %in% csf_signal_result$Entry) %>%
  filter(preStart < 50)
csf_dia_lib_semi_psm_semi_pep_sig_pro = csf_dia_lib_semi_psm_semi_pep_sig %>%
  group_by(`Protein ID`) %>%
  summarise(n=n())
write.table(csf_signal_result, "./supplements/table_s2.csv", row.names = F, quote = F, sep = "\t")


# Function to assign rows to peptides
assign_rows <- function(data) {
  rows <- rep(1, nrow(data))
  for (i in 2:nrow(data)) {
    for (j in 1:(i-1)) {
      if (is_overlap(data$Start[i], data$End[i], data$Start[j], data$End[j])) {
        rows[i] <- max(rows[i], rows[j] + 1)
      }
    }
  }
  rows
}
is_overlap <- function(start1, end1, start2, end2) {
  !(end1 <= start2 || start1 >= end2)
}
# Protein select
select_semi_protein = function(all_data, protein_id){
  all_data_filter = all_data %>%
    filter(str_equal(`Protein ID`,protein_id)) %>%
    distinct(Peptide, `Assigned Modifications`, `Protein Start`, `Protein End`, `Prev AA`, `Next AA`)
  all_data_filter$tryp = ifelse(all_data_filter$`Prev AA`%in% c("K", "R", "-") & 
                                              (str_ends(all_data_filter$Peptide, "K") | 
                                                 str_ends(all_data_filter$Peptide, "R") | 
                                                 all_data_filter$`Next AA` %in% c("-")), "y", "n")
  all_data_filter_plot = all_data_filter %>%
    select(Peptide, `Protein Start`, `Protein End`, tryp)
  colnames(all_data_filter_plot) = c("Peptide", "Start", "End", "tryp")
  all_data_filter_plot$Length = all_data_filter_plot$End - all_data_filter_plot$Start
  
  # Sort the data frame by the length of peptides
  all_data_filter_plot <- all_data_filter_plot[order(all_data_filter_plot$Length, decreasing = TRUE), ]
  print(head(all_data_filter_plot))
  
  all_data_filter_plot$Row <- assign_rows(all_data_filter_plot)
  return(all_data_filter_plot)
  
}


# c. Create the plot
csf_dia_lib_semi_psm_P02790_plot = select_semi_protein(csf_dia_lib_semi_psm, "P02790") %>%
  add_row(Peptide = "a", Start = 1, End=23, tryp = "s", Length = 23, Row = 0)%>%
  add_row(Peptide = "b", Start = 24, End=462, tryp = "o", Length = 439, Row = 0)

P02790_overlap_plot = ggplot(csf_dia_lib_semi_psm_P02790_plot, aes(x=Start, xend=End, y=Row, yend=Row, color=tryp)) +
  geom_segment(size=0.4) +
  theme_minimal() +
  labs(x=element_blank(), y=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank()) +
  scale_color_manual(values=c("y"="#D95F02", "n"="#42AE8D", "s"="red", "o"="#E6AB01","t"="orange", "i"="black")) +
  annotate(geom="text", x=11, y=-0.8, label="Signal", size = 2) +
  scale_y_continuous(limits = c(-1.5,12.5), expand = c(0, 0)) +
  theme(legend.position="none", plot.margin = unit(c(-10,0,-10,0), "mm")) +
  annotate(x=230, y=-1, label="HPX Protein Sequence",geom="text", color="Black", size=2)
  
P02790_overlap_plot

# E. mass offset result
csf_massoffset_summary = fread("./CSFData/mass_offset_result/ptm-shepherd-output/global.modsummary.tsv")
csf_massoffset_summary_used = csf_massoffset_summary %>%
  filter(!str_detect(Modification, "isotopic") & 
           !str_detect(Modification, "Unannotated mass-shift") & 
           !str_detect(Modification, "None") &
           `Mass Shift` >= -20 &
           `Mass Shift` <= 300)
colnames(csf_massoffset_summary_used)[2] = "mass_shift"
csf_massoffset_summary_used_label = csf_massoffset_summary_used %>%
  filter(Modification %in% c('Methyl', "Acetyl", "Phospho", "Hex", "Carbamyl" ,
                             "Oxidation", "Carbamidomethyl/Addition of G", 
                             "Didehydrobutyrine/Water loss", "Lysine not cleaved/Addition of K"))
csf_massoffset_summary_used_label$dataset01_percent_PSMs = csf_massoffset_summary_used_label$dataset01_percent_PSMs +1.4
csf_massoffset_summary_used_label$Modification = c("Carbamidomethyl\n(57.02)", "Carbamyl\n(43.00)", 
                                                   "Water loss\n(-18.01)", "Oxidation\n(15.99)", 
                                                   "Addition of K\n(128.095)", 'Methyl\n(14.015)', 
                                                   "Acetyl\n(42.01)", "Phospho\n(79.97)", "Hex\n(162.05)")

csf_massoffset_summary_used_plot = ggplot(csf_massoffset_summary_used, aes(x=mass_shift, ymax=dataset01_percent_PSMs, ymin=0)) +
  geom_linerange(linewidth = 0.1) +
  theme_light() +
  ylab("Percent of PSMs(%)") +
  xlab("Mass Shift") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 11)) +
  geom_text(data=csf_massoffset_summary_used_label,aes(x=mass_shift, y=dataset01_percent_PSMs, label=Modification), size=1.5) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 3),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "None") +
  theme(plot.margin = unit(c(0,0.2,0.1,0.5), "cm"))
csf_massoffset_summary_used_plot

figure2_bind_plot = ggarrange(ggarrange(csf_gg_dis_plot, csf_pr_dis_plot, widths = c(2.2, 2.2),
                                ncol = 2, nrow = 1, align="h", labels = c("a", "b"), font.label = list(size = 10)),
                      ggarrange(P02790_overlap_plot, csf_running_time_plot, widths = c(2, 1),
                                ncol = 2, nrow = 1, align="h", labels = c("c", "d"), font.label = list(size = 10)),
                      ggarrange(csf_massoffset_summary_used_plot, widths = c(1),
                                ncol = 1, nrow = 1, align="h", labels = c("e"), font.label = list(size = 10)),
                      nrow = 3, ncol=1, align = "v", heights = c(2, 1.5, 1))
figure2_bind_plot
ggsave("./figures/Figure2.pdf", figure2_bind_plot, width=5.3, height = 4.2, units = c("in"), dpi=400)

### Manually modification using AI later.


