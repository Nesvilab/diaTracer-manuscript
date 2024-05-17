library(tidyverse)
library(data.table)
library(ggpubr)

hla_pr_processed = fread("./HLAData/DIA_lib_result/diann-output/precursor_maxlfq.tsv")
hla_report = fread("./HLAData/DIA_lib_result/diann-output/report.tsv") %>%
  distinct(`Stripped.Sequence`, `Precursor.Id` )
colnames(hla_report) = c("stripped_seq", "V1")
hla_pr_processed_with_seq = inner_join(hla_pr_processed, hla_report, by="V1")
hla_pr_processed_with_seq$length = str_length(hla_pr_processed_with_seq$stripped_seq)

# a. comparison with spectronaut
hla_pr_processed_with_seq_strip_pep = hla_pr_processed_with_seq %>%
  distinct(stripped_seq)
hla_method = c("Maria et al.", "FP-diaTracer")
pep_num = c(2.13, nrow(hla_pr_processed_with_seq_strip_pep)/1000)
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
  filter(!str_equal(binder, "no_bin"))

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

# c.charge distribution
hla_pr_processed_with_seq$charge = str_sub(hla_pr_processed_with_seq$V1,-1,-1)
hla_pr_processed_with_seq_charge_num = hla_pr_processed_with_seq %>%
  group_by(charge) %>%
  summarise(n=n())
hla_pr_processed_with_seq_charge_num$smallN = hla_pr_processed_with_seq_charge_num$n/1000
hla_pr_processed_with_seq_charge_num_plot = ggplot(hla_pr_processed_with_seq_charge_num, aes(x=charge, y=smallN)) +
  geom_bar(stat="identity",  fill="#4DBBD5B2", color="black", size=0.05, width = 0.8) +
  ylab("# Immunopeptides\n(x1000)") +
  xlab("Charge") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05),
        legend.position = c(.85, .8))
hla_pr_processed_with_seq_charge_num_plot
# d. length distribution
hla_pr_processed_with_seq_strip_pep$length = str_length(hla_pr_processed_with_seq_strip_pep$stripped_seq)
hla_pr_processed_with_seq_strip_pep_len_num = hla_pr_processed_with_seq_strip_pep %>%
  group_by(length) %>%
  summarise(n=n())

hla_pr_processed_with_seq_strip_pep_len_num_short = hla_pr_processed_with_seq_strip_pep_len_num %>%
  filter(length <= 15)
hla_pr_processed_with_seq_strip_pep_len_num_short$smallN = hla_pr_processed_with_seq_strip_pep_len_num_short$n/1000
hla_pr_processed_with_seq_strip_pep_len_num_short_plot = ggplot(hla_pr_processed_with_seq_strip_pep_len_num_short, aes(x=length, y=smallN)) +
  geom_bar(stat="identity", fill="#4DBBD5B2", color="black", size=0.05, width = 0.8) +
  ylab("# Immunopeptides\n(x1000)") +
  xlab("Length") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  #scale_y_break(c(260, 2200 ), expand = c(0.01, 0)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))
hla_pr_processed_with_seq_strip_pep_len_num_short_plot

# Supplement
hla_pr_processed_with_seq_strip_pep_len_num$smallN = hla_pr_processed_with_seq_strip_pep_len_num$n/1000
hla_pr_processed_with_seq_strip_pep_len_num_plot = ggplot(hla_pr_processed_with_seq_strip_pep_len_num, aes(x=length, y=smallN)) +
  geom_bar(stat="identity", fill="#4DBBD5B2", color="black", size=0.05, width = 0.8) +
  ylab("# Immunopeptides\n(x1000)") +
  xlab("Length") +
  theme_light() +
  scale_y_continuous(expand = c(0.01, 0)) +
  #scale_y_break(c(260, 2200 ), expand = c(0.01, 0)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust=0.5, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.05))
hla_pr_processed_with_seq_strip_pep_len_num_plot
ggsave("./supplements/hla_full_length.pdf", hla_pr_processed_with_seq_strip_pep_len_num_plot, width=5.3, height = 4, units = c("in"), dpi=400)


bind_plot = ggarrange(ggarrange(hla_num_plot, netmhcpan_result_bind_num_plot, hla_pr_processed_with_seq_strip_pep_len_num_short_plot, widths = c(1.2, 2.2, 2.2),
                                ncol = 3, nrow = 1, align="h", labels = c("a", "b", "c"), font.label = list(size = 10)),
                      ggarrange(hla_pr_processed_with_seq_charge_num_plot, widths = c(1),
                                ncol = 3, nrow = 1, align="h", labels = c("d"), font.label = list(size = 10)),
                      nrow = 2, ncol=1, heights = c(2,2))
bind_plot
ggsave("./figures/Figure4Part.pdf", bind_plot, width=5.3, height = 4, units = c("in"), dpi=400)

### Assemble all figures using AI
