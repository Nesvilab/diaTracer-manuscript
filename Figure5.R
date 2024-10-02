library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpmisc)
library(GGally)
library(ggplot2)

get_mod_pos = function(diann_pr, frag_pep, pep_max){
  diann_pr = diann_pr %>%
    filter(`Modified.Sequence` %in% pep_max$V1)
  phospho_pep=list()
  position= list()
  probability = list()
  for (one_pho in diann_pr$`STY:79.96633`) {
    phos_comp = strsplit(one_pho, "\\(|)")
    pos_index = 0
    one_part_index = 1
    for (one_part in phos_comp[[1]]) {
      if (one_part_index%%2 == 1) {
        pos_index = pos_index + str_length(one_part)
      } else{
        prob = as.numeric(one_part)
        phospho_pep = append(phospho_pep, one_pho)
        position = append(position, pos_index)
        probability = append(probability, prob)
      }
      one_part_index = one_part_index +1
    }
    
  }
  
  phos_pos_table = data.frame(array(unlist(phospho_pep)), array(unlist(position)), array(unlist(probability)))
  colnames(phos_pos_table) = c("phospho_pep", "position", "probability")

  diann_pr_pep = diann_pr %>%
    distinct(`Protein.Ids`, `Stripped.Sequence`, `STY:79.96633`)
  colnames(diann_pr_pep) = c("protein_id", "pep_seq", "phospho_pep")
  phos_pos_table_hig = phos_pos_table %>%
    filter(probability >= 0.75)
  diann_pos = inner_join(diann_pr_pep, phos_pos_table_hig, by=c("phospho_pep"))
  
  colnames(frag_pep) = c("pep_seq", "protein_start", "protein_id")
  out_data = inner_join(diann_pos, frag_pep, by=c("protein_id", "pep_seq"))
  out_data$mod_pos = out_data$protein_start + out_data$position -1
  
  return(out_data)
}

get_mod_full_pos = function(diann_pr, frag_pep, pep_max){
  diann_pr = diann_pr %>%
    filter(`Modified.Sequence` %in% pep_max$V1)
  phospho_pep=list()
  position= list()
  probability = list()
  for (one_pho in diann_pr$`STY:79.96633`) {
    phos_comp = strsplit(one_pho, "\\(|)")
    pos_index = 0
    one_part_index = 1
    for (one_part in phos_comp[[1]]) {
      if (one_part_index%%2 == 1) {
        pos_index = pos_index + str_length(one_part)
      } else{
        prob = as.numeric(one_part)
        phospho_pep = append(phospho_pep, one_pho)
        position = append(position, pos_index)
        probability = append(probability, prob)
      }
      one_part_index = one_part_index +1
    }
    
  }
  
  phos_pos_table = data.frame(array(unlist(phospho_pep)), array(unlist(position)), array(unlist(probability)))
  colnames(phos_pos_table) = c("phospho_pep", "position", "probability")
  
  diann_pr_pep = diann_pr %>%
    distinct(`Protein.Ids`, `Stripped.Sequence`, `STY:79.96633`)
  colnames(diann_pr_pep) = c("protein_id", "pep_seq", "phospho_pep")
  
  diann_pos = inner_join(diann_pr_pep, phos_pos_table, by=c("phospho_pep"))
  
  colnames(frag_pep) = c("pep_seq", "protein_start", "protein_id")
  out_data = inner_join(diann_pos, frag_pep, by=c("protein_id", "pep_seq"))
  out_data$mod_pos = out_data$protein_start + out_data$position -1
  
  out_data_unique = out_data %>%
    distinct(protein_id, mod_pos, probability) %>%
    group_by(protein_id, mod_pos) %>% 
    slice(which.max(probability)) 
  out_data_unique$class = ifelse(out_data_unique$probability >0.75, "Class I", ifelse(out_data_unique$probability >0.5, "Class II", "Class III"))
  
  return(out_data_unique)
}

get_class_num = function(out_data_unique, gradient){
  out_data_unique_class = out_data_unique %>%
    group_by(class) %>%
    summarise(n=n())
  out_data_unique_class$Gradient = gradient
  
  return(out_data_unique_class)
}

get_strip_pep = function(diann_pr, pep_max){
  diann_pr_pass = diann_pr %>%
    filter(`Modified.Sequence` %in% pep_max$V1)
  
  diann_pr_pass_strip_pep = diann_pr_pass %>%
    distinct(`Stripped.Sequence`)
  
  return(nrow(diann_pr_pass_strip_pep))

}

get_strip_pep_spec = function(spec_report){
  spec_report_strip_pep = spec_report %>%
    filter(str_detect(`EG.ModifiedPeptide`, "Phospho")) %>%
    distinct(`PEP.StrippedSequence`)
  
  return(nrow(spec_report_strip_pep))
  
}

phospho_7min_mod_pep = fread("./phosphoData/7min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_7min_pr = fread("./phosphoData/7min_result/diann-output/report.pr_matrix.tsv") 
phospho_7min_pr_pass = phospho_7min_pr %>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_7min_pep_table = fread("./phosphoData/7min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_7min_mod_sites = get_mod_pos(phospho_7min_pr_pass, phospho_7min_pep_table, phospho_7min_mod_pep)
phospho_7min_mod_sites_uni = phospho_7min_mod_sites %>%
  distinct(protein_id, mod_pos)

phospho_10min_mod_pep = fread("./phosphoData/10min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_10min_pr = fread("./phosphoData/10min_result/diann-output/report.pr_matrix.tsv")
phospho_10min_pr_pass = phospho_10min_pr %>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_10min_pep_table = fread("./phosphoData/10min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_10min_mod_sites = get_mod_pos(phospho_10min_pr_pass, phospho_10min_pep_table, phospho_10min_mod_pep)
phospho_10min_mod_sites_uni = phospho_10min_mod_sites %>%
  distinct(protein_id, mod_pos)

phospho_15min_mod_pep = fread("./phosphoData/15min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_15min_pr = fread("./phosphoData/15min_result/diann-output/report.pr_matrix.tsv")
phospho_15min_pr_pass = phospho_15min_pr %>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_15min_pep_table = fread("./phosphoData/15min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_15min_mod_sites = get_mod_pos(phospho_15min_pr_pass, phospho_15min_pep_table, phospho_15min_mod_pep)
phospho_15min_mod_sites_uni = phospho_15min_mod_sites %>%
  distinct(protein_id, mod_pos)

phospho_21min_mod_pep = fread("./phosphoData/21min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_21min_pr = fread("./phosphoData/21min_result/diann-output/report.pr_matrix.tsv")
phospho_21min_pr_pass = phospho_21min_pr %>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_21min_pep_table = fread("./phosphoData/21min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_21min_mod_sites = get_mod_pos(phospho_21min_pr_pass, phospho_21min_pep_table, phospho_21min_mod_pep)
phospho_21min_mod_sites_uni = phospho_21min_mod_sites %>%
  distinct(protein_id, mod_pos)

phospho_30min_mod_pep = fread("./phosphoData/30min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_30min_pr = fread("./phosphoData/30min_result/diann-output/report.pr_matrix.tsv") 
phospho_30min_pr_pass = phospho_30min_pr%>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_30min_pep_table = fread("./phosphoData/30min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_30min_mod_sites = get_mod_pos(phospho_30min_pr_pass, phospho_30min_pep_table, phospho_30min_mod_pep)
phospho_30min_mod_sites_uni = phospho_30min_mod_sites %>%
  distinct(protein_id, mod_pos)

phospho_60min_mod_pep = fread("./phosphoData/60min_result/diann-output/modified_sequence_maxlfq.tsv")
phospho_60min_pr = fread("./phosphoData/60min_result/diann-output/report.pr_matrix.tsv")
phospho_60min_pr_pass = phospho_60min_pr %>%
  filter(`STY:79.96633 Best Localization` > 0.75)
phospho_60min_pep_table = fread("./phosphoData/60min_result/peptide.tsv") %>%
  distinct(Peptide, `Protein Start`, `Protein ID`)
phospho_60min_mod_sites = get_mod_pos(phospho_60min_pr_pass, phospho_60min_pep_table, phospho_60min_mod_pep)
phospho_60min_mod_sites_uni = phospho_60min_mod_sites %>%
  distinct(protein_id, mod_pos)

# a. Stripped seq num
phos_strip_pep_num = c(get_strip_pep(phospho_7min_pr, phospho_7min_mod_pep),
                  get_strip_pep(phospho_10min_pr, phospho_10min_mod_pep),
                  get_strip_pep(phospho_15min_pr, phospho_15min_mod_pep),
                  get_strip_pep(phospho_21min_pr, phospho_21min_mod_pep),
                  get_strip_pep(phospho_30min_pr, phospho_30min_mod_pep),
                  get_strip_pep(phospho_60min_pr, phospho_60min_mod_pep))
gradient_list = c("7", "10", "15", "21", "30", "60")
phos_strip_pep_num_plot_data = data.frame(gradient_list, phos_strip_pep_num)
phos_strip_pep_num_plot_data$small_phos_strip_pep_num = phos_strip_pep_num_plot_data$phos_strip_pep_num/1000

spec_7_min_report = fread("./phosphoData/Spectronaut analysis files/20240305_112743_phospho_7min_Report.tsv")
spec_10_min_report = fread("./phosphoData/Spectronaut analysis files/20240305_153512_phospho_10min_Report.tsv")
spec_15_min_report = fread("./phosphoData/Spectronaut analysis files/20240305_153519_phospho_15min_Report.tsv")
spec_21_min_report = fread("./phosphoData/Spectronaut analysis files/20240305_153525_Phospho_21min_Report.tsv")
spec_30_min_report = fread("./phosphoData/Spectronaut analysis files/20240305_153530_Phospho_30min_Report.tsv")
spec_60_min_report = fread("./phosphoData/Spectronaut analysis files/20240226_113655_Phospho_60min_Report.tsv")

phos_strip_pep_num_spec = c(get_strip_pep_spec(spec_7_min_report),
                            get_strip_pep_spec(spec_10_min_report),
                            get_strip_pep_spec(spec_15_min_report),
                            get_strip_pep_spec(spec_21_min_report),
                            get_strip_pep_spec(spec_30_min_report),
                            get_strip_pep_spec(spec_60_min_report))
phos_strip_pep_num_plot_data_spec = data.frame(gradient_list, `Oliinyk et al.`= phos_strip_pep_num_spec)
colnames(phos_strip_pep_num_plot_data)[2] = "FP-diaTracer"
phos_strip_pep_num_plot_data_spec_Frag = inner_join(phos_strip_pep_num_plot_data, phos_strip_pep_num_plot_data_spec, by = "gradient_list")
phos_strip_pep_num_plot_data_spec_Frag$small_phos_strip_pep_num = NULL
phos_strip_pep_num_plot_data_spec_Frag_gather = gather(phos_strip_pep_num_plot_data_spec_Frag, 2:3, key="method", value="Num")
phos_strip_pep_num_plot_data_spec_Frag_gather$small_sum = phos_strip_pep_num_plot_data_spec_Frag_gather$Num/1000

phos_sites_num = c(nrow(phospho_7min_mod_sites_uni),
                   nrow(phospho_10min_mod_sites_uni),
                   nrow(phospho_15min_mod_sites_uni),
                   nrow(phospho_21min_mod_sites_uni),
                   nrow(phospho_30min_mod_sites_uni),
                   nrow(phospho_60min_mod_sites_uni))
phos_sites_num_plot_data = data.frame(gradient_list, phos_sites_num)
phos_sites_num_plot_data$small_phos_sites_num = phos_sites_num_plot_data$phos_sites_num/1000
phos_sites_num_plot_data$method = "FP-diaTracer"

phos_strip_pep_num_spec_Frag_plot = ggplot() +
  geom_bar(data=phos_strip_pep_num_plot_data_spec_Frag_gather, aes(x=reorder(gradient_list, Num), y=small_sum, fill=method), stat="identity", color="black", size=0.05, width = 0.8, position=position_dodge()) +
  ylab("#Phosphopeptide Sequence \n(X1000)") +
  xlab("Gradient (min)") +
  theme_light() +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        legend.position = c(.18, .9))
phos_strip_pep_num_spec_Frag_plot

# intensity comp
phospho_7min_pr_rep1_2 = phospho_7min_pr %>%
  filter(!is.na(`STY:79.96633`)) %>%
  filter(`STY:79.96633 Best Localization` >=0.75)
phospho_7min_pr_rep1_2 = phospho_7min_pr_rep1_2[,11:14]
colnames(phospho_7min_pr_rep1_2) = c("Rep 1", "Rep 4", "Rep 3", "Rep 2")
phospho_7min_pr_rep1_2 = na.omit(phospho_7min_pr_rep1_2)
phospho_7min_pr_rep1_2$`Rep 1` = log10(phospho_7min_pr_rep1_2$`Rep 1`)
phospho_7min_pr_rep1_2$`Rep 4` = log10(phospho_7min_pr_rep1_2$`Rep 4`)
phospho_7min_pr_rep1_2$`Rep 3` = log10(phospho_7min_pr_rep1_2$`Rep 3`)
phospho_7min_pr_rep1_2$`Rep 2` = log10(phospho_7min_pr_rep1_2$`Rep 2`)
phospho_7min_int_cor_plot = ggpairs(phospho_7min_pr_rep1_2, 
        upper = list(method="p", continuous = wrap('cor', size = 2, stars = F)),
        lower = list(continuous = wrap("points", size=0.03, alpha = 0.4),
                     combo = wrap("dot", alpha = 0.4, size=0.03)),
        diag = list(continuous = "densityDiag", size=0.2))+ 
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 0.05),
        strip.text = element_text(
          size = 5, margin = margin(0,0,0,0, "cm")))
pdf("./figures/Figure5b.pdf", width=2.7, height = 2)
phospho_7min_int_cor_plot
dev.off()

# b. Sites num
phos_sites_num = c(nrow(phospho_7min_mod_sites_uni),
                   nrow(phospho_10min_mod_sites_uni),
                   nrow(phospho_15min_mod_sites_uni),
                   nrow(phospho_21min_mod_sites_uni),
                   nrow(phospho_30min_mod_sites_uni),
                   nrow(phospho_60min_mod_sites_uni))
phos_sites_num_plot_data = data.frame(gradient_list, phos_sites_num)
phos_sites_num_plot_data$small_phos_sites_num = phos_sites_num_plot_data$phos_sites_num/1000
phos_sites_num_plot = ggplot(phos_sites_num_plot_data, aes(x=reorder(gradient_list, small_phos_sites_num), y=small_phos_sites_num)) +
  geom_bar(stat="identity", fill="#4292C6", color="black", size=0.05, width = 0.8) +
  ylab("# Phosphorylation sites\n(x1000)") +
  xlab("Gradient (min)") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))

bind_plot = ggarrange(phos_strip_pep_num_spec_Frag_plot, widths = c(2.2),
                                ncol = 1, nrow = 1, align="h", labels = c("a"), font.label = list(size = 10))
ggsave("./figures/Figure5a.pdf", bind_plot, width=2.5, height = 2, units = c("in"), dpi=400)

### Assemble all figures using AI

# Supplement
phospho_7min_mod_sites_full = get_mod_full_pos(phospho_7min_pr, phospho_7min_pep_table, phospho_7min_mod_pep)
phospho_10min_mod_sites_full = get_mod_full_pos(phospho_10min_pr_pass, phospho_10min_pep_table, phospho_10min_mod_pep)
phospho_15min_mod_sites_full = get_mod_full_pos(phospho_15min_pr_pass, phospho_15min_pep_table, phospho_15min_mod_pep)
phospho_21min_mod_sites_full = get_mod_full_pos(phospho_21min_pr_pass, phospho_21min_pep_table, phospho_21min_mod_pep)
phospho_30min_mod_sites_full = get_mod_full_pos(phospho_30min_pr_pass, phospho_30min_pep_table, phospho_30min_mod_pep)
phospho_60min_mod_sites_full = get_mod_full_pos(phospho_60min_pr_pass, phospho_60min_pep_table, phospho_60min_mod_pep)

phospho_7min_mod_sites_full_class = get_class_num(phospho_7min_mod_sites_full, "7")
phospho_10min_mod_sites_full_class = get_class_num(phospho_10min_mod_sites_full, "10")
phospho_15min_mod_sites_full_class = get_class_num(phospho_15min_mod_sites_full, "15")
phospho_21min_mod_sites_full_class = get_class_num(phospho_21min_mod_sites_full, "21")
phospho_30min_mod_sites_full_class = get_class_num(phospho_30min_mod_sites_full, "30")
phospho_60min_mod_sites_full_class = get_class_num(phospho_60min_mod_sites_full, "60")
phos_sites_num_full = phospho_7min_mod_sites_full_class %>%
  rbind(phospho_10min_mod_sites_full_class) %>%
  rbind(phospho_15min_mod_sites_full_class) %>%
  rbind(phospho_21min_mod_sites_full_class) %>%
  rbind(phospho_30min_mod_sites_full_class) %>%
  rbind(phospho_60min_mod_sites_full_class)

phos_sites_num_full$Gradient <- factor(phos_sites_num_full$Gradient , levels=c("7", "10", "15", "21", "30", "60"))
phos_sites_num_full_plot = ggplot(phos_sites_num_full, aes(x=Gradient, y=n, fill=class)) +
  geom_bar(position="stack", stat="identity", color="black", size=0.05, width = 0.8) +
  ylab("# Phosphorylation sites\n(x1000)") +
  xlab("Gradient (min)") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))
phos_sites_num_full_plot
ggsave("./supplements/FigureS7.pdf", phos_sites_num_full_plot, width=3, height = 1.5, units = c("in"), dpi=400)

# CV
get_spectro_pass = function(spectronaut_report){
  spectronaut_report_ptm = spectronaut_report%>%
    distinct(EG.PrecursorId, `EG.PTMProbabilities [Phospho (STY)]`)
  spectronaut_report_ptm_sep = separate_rows(spectronaut_report_ptm, `EG.PTMProbabilities [Phospho (STY)]`, sep=";")
  spectronaut_report_ptm_sep_pass = spectronaut_report_ptm_sep %>%
    filter(`EG.PTMProbabilities [Phospho (STY)]` >= 0.75) %>%
    distinct(EG.PrecursorId)
  spectronaut_report_out = spectronaut_report %>%
    filter(EG.PrecursorId %in% spectronaut_report_ptm_sep_pass$EG.PrecursorId) %>%
    select(R.FileName, EG.PrecursorId, `EG.TotalQuantity (Settings)`)
  colnames(spectronaut_report_out) = c("Run", "PrecursorID", "Quantity")
  return(spectronaut_report_out)
}

get_diann_pass = function(diann_report){
  diann_report_pass = diann_report %>%
    filter(!is.na(`STY:79.96633`)) %>%
    filter(`STY:79.96633 Best Localization` >=0.75)
  diann_report_pass = diann_report_pass[,10:14]
  colnames(diann_report_pass)[1] = "PrecursorID"
  diann_report_out = gather(diann_report_pass, key = "Run", "Quantity", 2:5)
  return(diann_report_out)
}

calculate_cv = function(ptm_report, method, gradient, num){
  ptm_report_num = ptm_report %>%
    filter(!is.na(Quantity)) %>%
    group_by(PrecursorID) %>%
    summarise(n=n())%>%
    ungroup()%>%
    filter(n>num)
  ptm_report_out = ptm_report %>%
    filter(!is.na(Quantity)) %>%
    filter(PrecursorID %in% ptm_report_num$PrecursorID) %>%
    group_by(PrecursorID) %>%
    summarise_at(vars(Quantity), list(sd=sd, mean=mean))
  ptm_report_out$cv = (ptm_report_out$sd/ptm_report_out$mean) * 100
  ptm_report_out$method = method
  ptm_report_out$gradient = gradient
  return(ptm_report_out)
}

spec_7_min_report_ptm = get_spectro_pass(spec_7_min_report)
phospho_7min_pr_ptm = get_diann_pass(phospho_7min_pr)
spec_10_min_report_ptm = get_spectro_pass(spec_10_min_report)
phospho_10min_pr_ptm = get_diann_pass(phospho_10min_pr)
spec_15_min_report_ptm = get_spectro_pass(spec_15_min_report)
phospho_15min_pr_ptm = get_diann_pass(phospho_15min_pr)
spec_21_min_report_ptm = get_spectro_pass(spec_21_min_report)
phospho_21min_pr_ptm = get_diann_pass(phospho_21min_pr)
spec_30_min_report_ptm = get_spectro_pass(spec_30_min_report)
phospho_30min_pr_ptm = get_diann_pass(phospho_30min_pr)
spec_60_min_report_ptm = get_spectro_pass(spec_60_min_report)
phospho_60min_pr_ptm = get_diann_pass(phospho_60min_pr)

cv_num = 1
phospho_7min_pr_cv = calculate_cv(phospho_7min_pr_ptm, "FP-diaTracer", "7min", cv_num)
spec_7_min_report_cv = calculate_cv(spec_7_min_report_ptm, "Spectronaut", "7min", cv_num)
phospho_10min_pr_cv = calculate_cv(phospho_10min_pr_ptm, "FP-diaTracer", "10min", cv_num)
spec_10_min_report_cv = calculate_cv(spec_10_min_report_ptm, "Spectronaut", "10min", cv_num)
phospho_15min_pr_cv = calculate_cv(phospho_15min_pr_ptm, "FP-diaTracer", "15min", cv_num)
spec_15_min_report_cv = calculate_cv(spec_15_min_report_ptm, "Spectronaut", "15min", cv_num)
phospho_21min_pr_cv = calculate_cv(phospho_21min_pr_ptm, "FP-diaTracer", "21min", cv_num)
spec_21_min_report_cv = calculate_cv(spec_21_min_report_ptm, "Spectronaut", "21min", cv_num)
phospho_30min_pr_cv = calculate_cv(phospho_30min_pr_ptm, "FP-diaTracer", "30min", cv_num)
spec_30_min_report_cv = calculate_cv(spec_30_min_report_ptm, "Spectronaut", "30min", cv_num)
phospho_60min_pr_cv = calculate_cv(phospho_60min_pr_ptm, "FP-diaTracer", "60min", cv_num)
spec_60_min_report_cv = calculate_cv(spec_60_min_report_ptm, "Spectronaut", "60min", cv_num)

cv_plot_data = phospho_7min_pr_cv %>%
  bind_rows(spec_7_min_report_cv) %>%
  bind_rows(phospho_10min_pr_cv) %>%
  bind_rows(spec_10_min_report_cv) %>%
  bind_rows(phospho_15min_pr_cv) %>%
  bind_rows(spec_15_min_report_cv) %>%
  bind_rows(phospho_21min_pr_cv) %>%
  bind_rows(spec_21_min_report_cv) %>%
  bind_rows(phospho_30min_pr_cv) %>%
  bind_rows(spec_30_min_report_cv) %>%
  bind_rows(phospho_60min_pr_cv) %>%
  bind_rows(spec_60_min_report_cv)

cv_plot_data$gradient <- factor(cv_plot_data$gradient , levels=c("7min", "10min", "15min", "21min", "30min", "60min"))
cv_plot = ggplot(cv_plot_data, aes(x= gradient, y=cv, fill=method)) + 
  geom_boxplot(size=0.1, width=0.5, alpha=0.9, position = position_dodge(0.9), outliers = F) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +
  ylab("Coefficient Variation (%)") +
  xlab("Gradient") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "top",
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
cv_plot
cv_plot_out = ggplot(cv_plot_data, aes(x= gradient, y=cv, fill=method)) + 
  geom_boxplot(size=0.1, width=0.5, alpha=0.9, position = position_dodge(0.9), outlier.size = 0.1) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +
  ylab("Coefficient Variation (%)") +
  xlab("Gradient") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = "top",
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
cv_plot_all= ggarrange(cv_plot_out, cv_plot,
                       ncol = 1, nrow = 2, align="v", labels = c("a", "b"), font.label = list(size = 10), 
                       common.legend = T, legend = "top")
ggsave("./supplements/FigureS8.pdf", cv_plot_all, width=3, height = 3, units = c("in"), dpi=400)

### We tried to compare the site level, but lots of things are different. Hard to compare. e.x: Protein inference; localization scoring






