library(tidyverse)
library(data.table)
library(ggpubr)
library(ggpattern)

get_charge_data = function(charge_data){
  spec_di_report_charge_pro = charge_data %>%
    group_by(Precursor.Charge) %>%
    summarise(n=n())
  colnames(spec_di_report_charge_pro)[1] = "charge"
  spec_di_report_charge_pro$portion = spec_di_report_charge_pro$n/nrow(charge_data)
  spec_di_report_charge_pro <- spec_di_report_charge_pro %>%
    arrange(desc(charge)) %>%
    mutate(lab.ypos = cumsum(portion) - 0.5*portion)
  spec_di_report_charge_pro$portion_label = as.character(round(spec_di_report_charge_pro$portion, 3)*100)
  spec_di_report_charge_pro$portion_label = paste(spec_di_report_charge_pro$portion_label, "%", "")
  spec_di_report_charge_pro$charge = as.character(spec_di_report_charge_pro$charge)
  return(spec_di_report_charge_pro)
}

hla_pr_processed = fread("./HLAData/DIA_lib_result/diann-output/precursor_maxlfq.tsv")
hla_report = fread("./HLAData/DIA_lib_result/diann-output/report.tsv") %>%
  distinct(`Stripped.Sequence`, `Precursor.Id` )
colnames(hla_report) = c("stripped_seq", "V1")
hla_pr_processed_with_seq = inner_join(hla_pr_processed, hla_report, by="V1")
hla_pr_processed_with_seq$length = str_length(hla_pr_processed_with_seq$stripped_seq)

# a. comparison with spectronaut
spec_di_report = fread("./HLAData/original_study/Spectronaut_direct/20230612_142629_directDIA_whi40_Report.tsv")
spec_di_report_pep = spec_di_report%>%
  distinct(PEP.StrippedSequence)
hla_pr_processed_with_seq_strip_pep = hla_pr_processed_with_seq %>%
  distinct(stripped_seq) %>%
  filter(str_length(stripped_seq) <=12)
hla_method = c("Wahle et al.", "FP-diaTracer")
pep_num = c(nrow(spec_di_report_pep)/1000, nrow(hla_pr_processed_with_seq_strip_pep)/1000)
hla_num_plot_data = data.frame(hla_method, pep_num)
hla_num_plot = ggplot(hla_num_plot_data, aes(x=hla_method, y=pep_num, fill=hla_method)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E"))+
  scale_y_continuous(expand = c(0.01, 0)) +
  ylab("# Immunopeptides\n(x1000)") +
  xlab("Method") +
  theme_light() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))
hla_num_plot

### NetMhcpan
hla_pr_processed_with_seq_for_net = hla_pr_processed_with_seq %>%
  filter(length <=14) %>%
  filter(length >=8) %>%
  distinct(stripped_seq)
write.table(hla_pr_processed_with_seq_for_net$stripped_seq, "./HLAData/DIA_lib_result/diann-output/pure_pep_8_14.txt", quote = F, row.names = F, col.names = F)

netmhcpan_result = fread("./HLAData/DIA_lib_result/diann-output/netMHCPan_format_result.txt")

netmhcpan_result$binder = if_else(netmhcpan_result$EL_Rank>2, "no_bin", if_else(netmhcpan_result$EL_Rank<=0.5, "strong", "weak"))
netmhcpan_result_bind = netmhcpan_result %>%
  filter(!str_equal(binder, "no_bin")) %>%
  filter(str_length(Peptide) >= 8) %>%
  filter(str_length(Peptide) <= 12) 

netmhcpan_result_bind_pep = netmhcpan_result_bind %>%
  distinct(Peptide)

netmhcpan_result_bind_num = netmhcpan_result_bind %>%
  group_by(MHC, binder) %>%
  summarise(n=n()) %>%
  ungroup()

netmhcpan_result_bind_overall = netmhcpan_result_bind %>%
  distinct(Peptide, binder)
netmhcpan_result_bind_overall_strong = netmhcpan_result_bind_overall %>%
  filter(str_equal(binder, "strong")) %>%
  distinct(Peptide)

netmhcpan_result_bind_num = netmhcpan_result_bind_num %>%
  add_row(MHC="All", binder="strong", n=nrow(netmhcpan_result_bind_overall_strong)) %>%
  add_row(MHC="All", binder="weak", n=nrow(netmhcpan_result_bind_pep)-nrow(netmhcpan_result_bind_overall_strong))

netmhcpan_result_bind_num$binder = factor(netmhcpan_result_bind_num$binder, levels = c("weak", "strong"), ordered = TRUE)
# b. binder numbers
netmhcpan_result_bind_num$smallN = netmhcpan_result_bind_num$n/1000
netmhcpan_result_bind_num_plot = ggplot(netmhcpan_result_bind_num, aes(x=MHC, y=smallN, fill=binder)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) +
  scale_fill_brewer(palette="Blues") + 
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=0.2)+
  geom_vline(xintercept = 3.5, linetype="dashed", color = "black", size=0.2)+
  geom_vline(xintercept = 5.5, linetype="dashed", color = "black", size=0.2)+
  ylab("# Predicted Binder\n(x1000)") +
  xlab("Allele") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.87, .8),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
netmhcpan_result_bind_num_plot

# c. length distribution
spec_di_netmhcpan = fread("./HLAData/original_study/Spectronaut_direct/spectronaut_netMHCpan_predictions.tsv")
spec_di_netmhcpan_pass = spec_di_netmhcpan %>%
  filter(Binder != "Non-binder") %>%
  distinct(Peptide)

pan_library_netmhcpan = fread("./HLAData/original_study/pan_library/whi40_pan_netMHCpan_predictions.tsv")
pan_library_netmhcpan_pass = pan_library_netmhcpan %>%
  filter(Binder != "Non-binder") %>%
  distinct(Peptide)

dda_library_netmhcpan = fread("./HLAData/original_study/dda_library/whi40_fragger_netMHCpan_predictions.tsv")
dda_library_netmhcpan_pass = dda_library_netmhcpan %>%
  filter(Binder != "Non-binder") %>%
  distinct(Peptide)

diaTracer_netmhcpan = fread("./HLAData/DIA_lib_result/diann-output/netMHCPan_format_result.txt")
diaTracer_netmhcpan_pass = diaTracer_netmhcpan %>%
  filter(EL_Rank<=2) %>%
  distinct(Peptide) %>%
  filter(str_length(Peptide) <= 12)

spec_di_netmhcpan_pass$length = str_length(spec_di_netmhcpan_pass$Peptide)
diaTracer_netmhcpan_pass$length = str_length(diaTracer_netmhcpan_pass$Peptide)
dda_library_netmhcpan_pass$length = str_length(dda_library_netmhcpan_pass$Peptide)
pan_library_netmhcpan_pass$length = str_length(pan_library_netmhcpan_pass$Peptide)
spec_di_netmhcpan_pass$class="directDIA"
diaTracer_netmhcpan_pass$class="FP-diaTracer"
dda_library_netmhcpan_pass$class="Experimental DDA Library"
pan_library_netmhcpan_pass$class = "panlibrary"

length_dist = spec_di_netmhcpan_pass%>%
  bind_rows(diaTracer_netmhcpan_pass)%>%
  bind_rows(dda_library_netmhcpan_pass) %>%
  bind_rows(pan_library_netmhcpan_pass)

length_dist_por = length_dist %>%
  group_by(length, class) %>%
  summarise(num=n()) %>%
  ungroup() %>%
  filter(length<13)
length_dist_por_total = length_dist_por %>%
  group_by(class) %>%
  summarise(total=sum(num))

length_dist_por = left_join(length_dist_por, length_dist_por_total, by="class")
length_dist_por$portion = length_dist_por$num/length_dist_por$total
colnames(length_dist_por)[2] = "Method"
length_dist_plot = ggplot(length_dist_por, aes(x=length, y=portion, fill=Method))+
  geom_bar(stat="identity", color="black", size=0.05, width = 0.8,position=position_dodge()) +
  ylab("Proportion") +
  xlab("Length") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = c("#1F77B4",  "#0FA18A", "#3D5587", "#FF7F0E"))+
  #scale_y_break(c(260, 2200 ), expand = c(0.01, 0)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.75, .8),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
length_dist_plot


# d.charge distribution
dda_library_report = fread("./HLAData/original_study/dda_library/whi40_fraggerlib_report.tsv")
pan_library_report = fread("./HLAData/original_study/pan_library/panlib_report.tsv")
diaTracer_report = fread("./HLAData/DIA_lib_result/diann-output/report.tsv") %>%
  filter(`Precursor.Id` %in% hla_pr_processed$V1)

spec_di_report_charge = spec_di_report %>%
  filter(PEP.StrippedSequence %in% spec_di_netmhcpan_pass$Peptide) %>%
  distinct(PEP.StrippedSequence, FG.Charge)
colnames(spec_di_report_charge) = c("Stripped.Sequence", "Precursor.Charge")
dda_library_report_charge = dda_library_report %>%
  filter(Stripped.Sequence %in% dda_library_netmhcpan_pass$Peptide) %>%
  distinct(Stripped.Sequence, Precursor.Charge)
pan_library_report_charge = pan_library_report %>%
  filter(Stripped.Sequence %in% pan_library_netmhcpan_pass$Peptide) %>%
  distinct(Stripped.Sequence, Precursor.Charge)
diaTracer_report_charge = diaTracer_report %>%
  filter(Stripped.Sequence %in% diaTracer_netmhcpan_pass$Peptide) %>%
  distinct(Stripped.Sequence, Precursor.Charge)

spec_di_report_charge_plot_data = get_charge_data(spec_di_report_charge)
dda_library_report_charge_plot_data = get_charge_data(dda_library_report_charge)
pan_library_report_charge_plot_data = get_charge_data(pan_library_report_charge)
diaTracer_report_charge_plot_data = get_charge_data(diaTracer_report_charge)

mycols <- c( "#0FA18A", "#FF7F0E",  "#3D5587")

spec_di_report_charge_plot = ggplot(spec_di_report_charge_plot_data, aes(x="", y=portion, fill=charge)) +
  geom_bar(width = 0.5, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(x= 1.2, y = lab.ypos, label = portion_label), color = "black", size=2)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  theme(legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(-10, 2, 0, -10))
spec_di_report_charge_plot

diaTracer_report_charge_plot = ggplot(diaTracer_report_charge_plot_data, aes(x="", y=portion, fill=charge)) +
  geom_bar(width = 0.5, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(x= 1.2, y = lab.ypos, label = portion_label), color = "black", size=2)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  theme(legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(-10, 2, 0, -10))
diaTracer_report_charge_plot


# Supplement Upset plot
spec_di_netmhcpan_pass$`Spectronaut directDIA` = 1
pan_library_netmhcpan_pass$PanLibrary = 1
dda_library_netmhcpan_pass$`experimental DDA library` = 1
diaTracer_netmhcpan_pass$`FP-diaTracer` = 1

hla_upset_data = spec_di_netmhcpan_pass %>%
  full_join(pan_library_netmhcpan_pass, by="Peptide") %>%
  full_join(dda_library_netmhcpan_pass, by="Peptide") %>%
  full_join(diaTracer_netmhcpan_pass, by="Peptide")

hla_upset_data[is.na(hla_upset_data)] <- 0
m = make_comb_mat(hla_upset_data)
cs = comb_size(m)
hla_immu_overlap = UpSet(m,
      top_annotation = HeatmapAnnotation(
        "Predicted binders\n Intersections" = anno_barplot(cs, 
                                                           ylim = c(0, max(cs)*1.1),
                                                           border = FALSE, 
                                                           gp = gpar(fill = "black"), 
                                                           height = unit(4, "cm")
        ), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),)

pdf("./supplements/hla_immu_overlap.pdf", width=7, height = 4)
hla_immu_overlap
dev.off()

# spectronaut result binder distribution
spec_di_netmhcpan$binder = if_else(spec_di_netmhcpan$EL_Rank>2, "no_bin", if_else(spec_di_netmhcpan$EL_Rank<=0.5, "strong", "weak"))
spec_di_netmhcpan_bind = spec_di_netmhcpan %>%
  filter(!str_equal(binder, "no_bin"))

spec_di_netmhcpan_bind_pep = spec_di_netmhcpan_bind %>%
  distinct(Peptide)

spec_di_netmhcpan_bind_num = spec_di_netmhcpan_bind %>%
  group_by(Allele, binder) %>%
  summarise(n=n()) %>%
  ungroup()

spec_di_netmhcpan_bind_overall = spec_di_netmhcpan_bind %>%
  distinct(Peptide, binder)
spec_di_netmhcpan_bind_overall_strong = spec_di_netmhcpan_bind_overall %>%
  filter(str_equal(binder, "strong")) %>%
  distinct(Peptide)

spec_di_netmhcpan_bind_num = spec_di_netmhcpan_bind_num %>%
  add_row(Allele="All", binder="strong", n=nrow(spec_di_netmhcpan_bind_overall_strong)) %>%
  add_row(Allele="All", binder="weak", n=nrow(spec_di_netmhcpan_bind_pep)-nrow(spec_di_netmhcpan_bind_overall_strong))

spec_di_netmhcpan_bind_num$binder = factor(spec_di_netmhcpan_bind_num$binder, levels = c("weak", "strong"), ordered = TRUE)
spec_di_netmhcpan_bind_num$smallN = spec_di_netmhcpan_bind_num$n/1000
spec_di_netmhcpan_bind_num_plot = ggplot(spec_di_netmhcpan_bind_num, aes(x=Allele, y=smallN, fill=binder)) +
  geom_bar(stat="identity",  color="black", size=0.05, width = 0.8) +
  scale_fill_brewer(palette="Blues") + 
  geom_vline(xintercept = 1.5, linetype="dashed", color = "black", size=0.2)+
  geom_vline(xintercept = 3.5, linetype="dashed", color = "black", size=0.2)+
  geom_vline(xintercept = 5.5, linetype="dashed", color = "black", size=0.2)+
  ylab("# Predicted Binder\n(x1000)") +
  xlab("Allele") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.87, .8),
        legend.title = element_text(size=5, face="bold"),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"))
spec_di_netmhcpan_bind_num_plot
ggsave("./supplements/spectronaut_binders.pdf", spec_di_netmhcpan_bind_num_plot, width=3, height = 3, units = c("in"), dpi=400)

bind_plot = ggarrange(ggarrange(hla_num_plot, netmhcpan_result_bind_num_plot, length_dist_plot, widths = c(1.2, 2.2, 2.2),
                                ncol = 3, nrow = 1, align="h", labels = c("a", "b", "c"), font.label = list(size = 10)),
                      ggarrange(
                        ggarrange(diaTracer_report_charge_plot, spec_di_report_charge_plot, widths = c(1,1),
                                  ncol = 1, nrow = 2, align = "v", labels = c("d", "e"), font.label = list(size = 10),
                                  common.legend = T, legend = "right"), 
                        widths = c(1), ncol = 3, nrow = 1, align="h"),
                      nrow = 2, ncol=1, heights = c(2,2))
bind_plot
ggsave("./figures/Figure4Part.pdf", bind_plot, width=5.3, height = 4, units = c("in"), dpi=400)

### Assemble all figures using AI
