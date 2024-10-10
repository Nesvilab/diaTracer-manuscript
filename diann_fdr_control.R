library(iq)
library(stringr)
library(matrixStats)


generate_new_reports = function(path){
  
  globalPrecursorFDR <- 0.01
  globalProteinFDR <- 0.01
  runSpecificPrecursorFDR <- 0.01
  runSpecificProteinFDR <- 1.1

  out_path <- str_c(dirname(path), "/precursor_maxlfq.tsv")
  if (file.exists(out_path) == FALSE) {
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))
    
    df <- fast_read(path,
                    sample_id = "Run",
                    primary_id = "Precursor.Id",
                    secondary_id = "Precursor.Id",
                    intensity_col = "Precursor.Normalised",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                    filter_double_greater = NULL,
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "0")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
  }
  
  out_path <- str_c(dirname(path), "/modified_sequence_maxlfq.tsv")
  if (file.exists(out_path) == FALSE) {
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))
    
    df <- fast_read(path,
                    sample_id = "Run",
                    primary_id = "Modified.Sequence",
                    secondary_id = "Precursor.Id",
                    intensity_col = "Precursor.Normalised",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                    filter_double_greater = NULL,
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "0")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, log2_intensity_cutoff = 0, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
  }
  
  out_path <- str_c(dirname(path), "/protein_maxlfq.tsv")
  if (file.exists(out_path) == FALSE) {
    print(out_path)
    print(str_c("global precursor FDR = ", globalPrecursorFDR))
    print(str_c("global protein FDR = ", globalProteinFDR))
    print(str_c("run specific precursor FDR = ", runSpecificPrecursorFDR))
    print(str_c("run specific protein FDR = ", runSpecificProteinFDR))
    
    df <- fast_read(path,
                    sample_id = "Run",
                    primary_id = "Protein.Group",
                    secondary_id = "Precursor.Id",
                    intensity_col = "Precursor.Normalised",
                    annotation_col = NULL,
                    filter_string_equal = NULL,
                    filter_string_not_equal = NULL,
                    filter_double_less = c("Global.Q.Value" = globalPrecursorFDR, "Global.PG.Q.Value" = globalProteinFDR, "Q.Value" = runSpecificPrecursorFDR, "PG.Q.Value" = runSpecificProteinFDR),
                    filter_double_greater = NULL,
                    intensity_col_sep = NULL,
                    intensity_col_id = NULL,
                    na_string = "0")
    df_norm <- fast_preprocess(df$quant_table, median_normalization = FALSE, pdf_out = NULL) # Do not enable median normalization because the tools have their own normalization algorithm.
    df_maxlfq <- fast_MaxLFQ(df_norm, row_names = df$protein[, 1], col_names = df$sample) # MaxLFQ requires a precursor having non-zero intensities in at least two runs, which makes the "quantified" entries slightly fewer than those from the DIA-NN main report.
    maxlfq <- df_maxlfq$estimate
    maxlfq[maxlfq <= 0] <- NA
    write.table(2^maxlfq, out_path, quote = FALSE, sep = "\t", col.names = NA)
  }
  
}

### CSF
generate_new_reports("./CSFData/DDA_lib_try_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DDA_lib_semi_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DDA_DIA_lib_semi_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DDA_DIA_lib_try_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DIA_lib_semi_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DIA_lib_try_result/diann-output/report.tsv")
generate_new_reports("./CSFData/DIA_NN_LF_result/report.tsv")

### Plasma
generate_new_reports("./plasmaData/DIA_lib_semi_result/diann-output/report.tsv")
generate_new_reports("./plasmaData/DIA_lib_try_result/diann-output/report.tsv")
generate_new_reports("./plasmaData/DIA_NN_LF_result/report.tsv")

### HLA
generate_new_reports("./HLAData/DIA_lib_result/diann-output/report.tsv")

### PhosphoData
generate_new_reports("./phosphoData/7min_result/diann-output/report.tsv")
generate_new_reports("./phosphoData/10min_result/diann-output/report.tsv")
generate_new_reports("./phosphoData/15min_result/diann-output/report.tsv")
generate_new_reports("./phosphoData/21min_result/diann-output/report.tsv")
generate_new_reports("./phosphoData/30min_result/diann-output/report.tsv")
generate_new_reports("./phosphoData/60min_result/diann-output/report.tsv")

### Low input
generate_new_reports("./lowInputData/Tonsil_lib_result/diann-output/report.tsv")
generate_new_reports("./lowInputData/Direct_lib_result/diann-output/report.tsv")
generate_new_reports("./lowInputData/DIANN_results_Fig4_original_study/report.tsv")

### Breast cancer
generate_new_reports("./revisionData/TNBCData/14_MainSearch_4223_DIA-NN_1.8.1_library_free/14_MainSearch_4223_DIA-NN_1.8.1_library_free/report.tsv")
generate_new_reports("./revisionData/TNBCData/diaTracer_hybrid_result_Frag22_methylthiolation/diann-output/report.tsv")
generate_new_reports("./revisionData/TNBCData/diaTracer_result_Frag22_methylthiolation/diann-output/report.tsv")
generate_new_reports("./revisionData/TNBCData/12_MainSearch_4223_DIA-NN_1.8.1_library_based/12_MainSearch_4223_DIA-NN_1.8.1_newSN16lib/report.tsv")


