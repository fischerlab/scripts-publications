library(tidyverse)
library(GGally)
library(ggrepel)


# Input file has to have the following structure: Col1 - "Gene Name", Col2 - "sgRNA sequence", remaining columns should contain read number data. 
# Colnames should contain the following elements, separated by underscores: Experiment-name, Cell-line-name, Treatment, Selection-criteria (e.g. sorting-gate), replicate. "For example BCL6.reporter.screen_HEK293T_BI3802_High.gate_RepA"

# ==== Section 1: Defining functions ====

# Function 1: Reads in a table containing all gene names, sgRNA sequences and read counts for one experiment and selects the relevant conditions for further comparison
split_table_function <- function (experiment_path, input_file_name, screen_name) {
  input_data <- read_csv(paste0(experiment_path, "/", input_file_name))
  colnames(input_data)[c(1,2)] <- c("Gene", "Sequence")
  colnames_for_single_comparison <- colnames(input_data) %>% str_subset(screen_name)
  selected_input_data <- input_data %>% select(Gene, Sequence, one_of(colnames_for_single_comparison))
  return(selected_input_data)
}

split_output_table <- split_table_function(experiment_path)

# Function 2: Transforms the wide data into a long data format for normalization
long_data_transformation_function <- function(wide_selected_input_data, screen_name) {
  long_selected_input_data <- wide_selected_input_data %>% 
    gather(key = "Experiment", value = "Read_numbers", -Gene, -Sequence) %>% 
    separate(Experiment, into=c("Experiment",  "Cell_line", "Drug", "Gate", "Replicate"), sep="_")
  return(long_selected_input_data)}

# Function 3: Normalizes the read number to total read number across all selected conditions in Function 1
normalization_function <- function(long_selected_input_data) {
  # substitute 0 and NAs with 1
  long_selected_input_data [long_selected_input_data == 0 ] <- 1
  long_selected_input_data [is.na(long_selected_input_data)] <- 1
  
  # Read normalization (Every read number is divided by the mean of the reads within single sample followed by multiplication of the average read number for the whole experiment)
  normalized_reads <- long_selected_input_data %>% 
    group_by(Experiment, Cell_line, Drug, Gate, Replicate) %>% 
    mutate(Norm_read_numbers = Read_numbers/mean(Read_numbers)) %>% 
    ungroup() %>% 
    mutate(Norm_read_numbers = Norm_read_numbers * mean(Read_numbers))
  
  return(normalized_reads)
}

# Function 4: Transforms the long normalized reads table generated in Function 3 into a wide table for correlation plotting
wide_normalized_reads_function <- function(normalized_reads) {
  normalized_reads_wide <- normalized_reads %>%
    select(-Read_numbers) %>% 
    unite(col = "Experiment", Experiment, Cell_line, Drug, Gate, Replicate, sep="_") %>% 
    spread(key = Experiment, value = Norm_read_numbers)
  return(normalized_reads_wide)
}

# Function 5: Generates pair-wise correlations across all selected conditions in Function 1
read_correlation_function <- function(normalized_reads_wide) {
  read_count_log <- normalized_reads_wide %>% select_if(is.numeric)  %>% as.matrix() %>% log() %>% as.data.frame()
  correlation_key_word_plot <- 
    ggpairs(read_count_log, columns = 1:(ncol(read_count_log)),
            upper = list(continuous = wrap("density", alpha = 0.8, col = "steelblue"), combo = "box_no_facet"),
            lower = list(continuous = wrap("points", alpha = 0.2, col = "black", shape = 16, size = 0.8), combo = wrap("dot_no_facet", alpha = 0.4))) +
    theme_classic(base_size = 6) +
    theme(text = element_text(family = "Helvetica"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"))
  return(correlation_key_word_plot)
}


# Function 6: Calculates the ratio of gate 1 to gate 2 for each guide RNA and ranks them based on the ratio. Replicates are combined by summing up ranks 
# and by calculating the median ratio for each guide. Finally, guides are merged to genes by calculating the medium ranks and ratios.
rank_guides = function(df_norm, gate_list, replicate_list, screen_name) {
  ranked_df <- df_norm %>%
    filter(Replicate %in% replicate_list,
           Gate %in% gate_list) %>%
    # Combine gates
    group_by(Gene, Sequence, Replicate) %>%
    summarise(Ratio = Norm_read_numbers[Gate == gate_list[1]]/Norm_read_numbers[Gate == gate_list[2]]) %>%
    ungroup() %>%
    group_by(Replicate) %>%
    mutate(Rank = rank(desc(Ratio))) %>%
    ungroup() %>%
    # Combine replicates
    group_by(Gene, Sequence) %>%
    summarize(Rank = sum(Rank), Ratio = median(Ratio)) %>%
    ungroup() %>%
    # Combine sgRNAs
    group_by(Gene) %>%
    summarise(Rank = median(Rank), Ratio = median(Ratio))
  return(ranked_df)
}    

# Function 7: Simulates a p value distribution with guide RNAs that were assigned random ranks
simulate_p_value <- function(total_sgRNAs = 2852, n_iterations = 500, replicate_list = c("Rep1", "Rep2"), n_sgRNAs = 4) {
  null_rank_matrix <- matrix(nrow = total_sgRNAs, ncol = n_iterations)
  
  for (i in 1:n_iterations) {
    null_rank_matrix[,i] <- rowSums(cbind(replicate(length(replicate_list), sample(1:total_sgRNAs, total_sgRNAs, replace=FALSE))))
    
  }
  
  null_rank_df <- data.frame(null_rank_matrix) 
  
  null_rank_df <- null_rank_df %>%
    mutate(ID = rep(1:(nrow(null_rank_df)/n_sgRNAs), length.out = nrow(null_rank_df))) %>%
    group_by(ID) %>%
    summarise_all(median, na.rm = T)
  
  null_rank_vector <- data.matrix(null_rank_df[-1]) %>% as.vector()
  return(null_rank_vector)
}

# Function 8: Tests the calculated ranks (Function 6) against the simulated distribution (Function 7). Returns non-corrected p values and false-discovery rates (FDR)
calculate_p_value <- function(ranked_df, null_rank_vector, total_sgRNAs = 2852, n_iterations = 500, n_sgRNAs = 4, screen_name, FDR_thr = 0.1) {
  calculate_pVal <- function(i){
    sum(null_rank_vector <= as.numeric(ranked_df$Rank[i]))/(total_sgRNAs*(n_iterations/n_sgRNAs))
  }
  
  p_value <- map_dbl(1:(total_sgRNAs/n_sgRNAs), calculate_pVal)
  
  p_value_df <- data_frame(Gene = ranked_df$Gene, Ratio = ranked_df$Ratio, p_value) %>%
    mutate(p_value = ifelse(p_value < 0.5, p_value*2, (1 - p_value)*2),
           p_value = ifelse(p_value == 0.0, min(p_value[p_value > 0]), p_value)) %>%
    mutate(  FDR = p.adjust(p_value, method = "BH"),
             #q_value = qvalue(p_value)$qvalue,
             hits = ifelse(FDR < FDR_thr, T, F) ) %>%
    arrange(desc(p_value))
  return(p_value_df)
}

# Function 9: Generates two volcano plots for a selected display of hits or every hit that makes a specified FDR threshold
volcanoplot_screen <-  function(p_value_df, selected_hits, screen_name, gate_list = c("H", "L"), main_path = main_path) {
  experiment_path <- paste0(main_path, screen_name, "/")
  # Generates a volcano plot where x is the ratio of the sum of counts for each guide targeting a gene and y the p-value
  p = ggplot(p_value_df, aes(x = Ratio, y = -log10(p_value))) +
    geom_point(alpha = 0.4, shape = 16, size = 0.6, col = "darkgray") +
    scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) +
    scale_y_continuous(limits = c(0, -log10(min(p_value_df$p_value)) + 0.5)) +
    labs(title = screen_name, x = paste("Enrichment in", gate_list[1], "vs", gate_list[2]), y = expression('-log'[10]*'(p-value)')) +
    theme_classic(base_size = 6) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      title = element_text(size = 6),
      strip.background=element_blank())
  
  p1 = p + geom_point(data = filter(p_value_df, hits == T), alpha = 1, shape = 21, size = 1.5, fill = "#F45053") +
    geom_text_repel(
      data = filter(p_value_df, hits == T),
      aes(label = Gene),
      size = 2,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.3, "lines"))
  
  p2 = p + geom_point(data = filter(p_value_df, Gene %in% selected_hits), alpha = 1, shape = 21, size = 1.5, fill = "#F45053") +
    geom_text_repel(
      data = filter(p_value_df, Gene %in% selected_hits),
      aes(label = Gene),
      size = 2,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.3, "lines"))
  
  ggsave(paste0(experiment_path, "Volcano_", screen_name, ".pdf"), p1, width = 4.5, height = 4.5, units = "cm")
  ggsave(paste0(experiment_path, "Volcano_selected_", screen_name, ".pdf"), p2, width = 4.5, height = 4.5, units = "cm")
  
}

# Function 10: Generates a plot for read numbers of the individual guide RNAs targeting genes that make the FDR threshold
plot_individual_guides = function(df_norm, p_value_df, selected_hits, screen_name, gate_list = c("H", "L"), replicate_list = c("Rep1", "Rep2", "Rep3"), main_path = main_path) {
  experiment_path <- paste0(main_path, screen_name, "/")
  hit_list <- filter(p_value_df, hits == T)$Gene
  
  guides_df <- df_norm %>%
    filter(Replicate %in% replicate_list,
           Gate %in% gate_list) %>%
    # Combine replicates
    group_by(Gene, Sequence, Gate) %>%
    summarise(Count = mean(Norm_read_numbers)) %>%
    ungroup()
  
  # Arrange guides by mean read count in the high gate
  levels = filter(guides_df, Gate == gate_list[1]) %>%
    group_by(Gene) %>%
    summarize(Mean = mean(Count)) %>%
    arrange(desc(Mean)) %>%
    `$`(Gene) %>%
    as.character()
  
  
  posn_jd = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)
  posn_d = position_dodge(width = 0.4)
  
  p1 = ggplot(filter(guides_df, Gene %in% hit_list), aes(x = factor(Gene, levels = levels), y = Count, fill = Gate, group = Gate)) +
    stat_summary(geom = "crossbar", 
                 fun.ymin = mean,
                 fun.ymax = mean,
                 fun.y = mean,
                 width = 0.2,
                 size = 0.3,
                 position = posn_d,
                 fatten = 1,
                 color = "gray",
                 show.legend = F) +
    stat_summary(geom = "errorbar", 
                 fun.ymin = function(x) {mean(x)-sd(x)}, 
                 fun.ymax = function(x) {mean(x)+sd(x)},
                 width = 0.1,
                 size = 0.3,
                 position = posn_d,
                 color = "gray",
                 show.legend = F) +
    geom_point(position = posn_jd, shape = 21, size = 1) +
    labs(y = "Normalized read count") +
    scale_fill_manual("Gate", values = c("#3397A9", "darkred")) +
    theme_classic(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          text = element_text(family = "Helvetica")
    )
  
  p2 = ggplot(filter(guides_df, Gene %in% selected_hits), aes(x = factor(Gene, levels = levels), y = Count, fill = Gate, group = Gate)) +
    stat_summary(geom = "crossbar", 
                 fun.ymin = mean,
                 fun.ymax = mean,
                 fun.y = mean,
                 width = 0.2,
                 size = 0.3,
                 position = posn_d,
                 fatten = 1,
                 color = "gray",
                 show.legend = F) +
    stat_summary(geom = "errorbar", 
                 fun.ymin = function(x) {mean(x)-sd(x)}, 
                 fun.ymax = function(x) {mean(x)+sd(x)},
                 width = 0.1,
                 size = 0.3,
                 position = posn_d,
                 color = "gray",
                 show.legend = F) +
    geom_point(position = posn_jd, shape = 21, size = 1) +
    labs(y = "Normalized read count") +
    scale_fill_manual("Gate", values = c("#00B050", "red")) +
    theme_classic(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          text = element_text(family = "Helvetica")
          )
  
  
  
  ggsave(paste0(experiment_path, "Individual_guides_", screen_name, ".pdf"), p1, width = 8, height = 5, units = "cm")
  ggsave(paste0(experiment_path, "Individual_guides_selected_", screen_name, ".pdf"), p2, width = 8, height = 5, units = "cm")
}



# Wrapper function that runs all previous functions in sequence
run_pipeline <- function(main_path, input_file_name, screen_name, gate_list = c("H", "L"), replicate_list = c("Rep1", "Rep2", "Rep3"), 
                         n_iterations = 10, n_sgRNAs = 4, 
                         selected_hits, 
                         FDR_thr = 0.1) {
  
  experiment_path <- paste0(main_path, screen_name, "/")
  if (!dir.exists(experiment_path)){dir.create(experiment_path)}
  input_df <- read_csv(paste0(main_path, input_file_name))

  output_function_1 <- long_data_transformation_function(input_df, screen_name)
  total_sgRNAs <- unique(output_function_1$Gene) %>% length() * n_sgRNAs
  output_function_2 <- normalization_function (output_function_1)
  write_csv(output_function_2, paste0(experiment_path, "Normalized_reads_long", screen_name, ".csv"))
  output_function_3 <- wide_normalized_reads_function  (output_function_2)
  write_csv(output_function_3, paste0(experiment_path, "Normalized_reads_", screen_name, ".csv"))
  correlation_plots <- read_correlation_function(output_function_3)
  ggsave(paste0(experiment_path, "Correlation_", screen_name, ".png"), correlation_plots, width = 10, height = 10)
  
  output_function_4 <- rank_guides(df_norm = output_function_2, gate_list = gate_list, replicate_list = replicate_list, screen_name = screen_name)
  write_csv(output_function_4, paste0(experiment_path, "Ranked_", screen_name, ".csv"))
  output_function_5 <- simulate_p_value(total_sgRNAs = total_sgRNAs, n_iterations = n_iterations, replicate_list = replicate_list, n_sgRNAs = n_sgRNAs)
  output_function_6 <- calculate_p_value(ranked_df = output_function_4, null_rank_vector = output_function_5, total_sgRNAs = total_sgRNAs, n_iterations = n_iterations, n_sgRNAs = n_sgRNAs, screen_name = screen_name, FDR_thr = FDR_thr)
  write_csv(output_function_6, paste0(experiment_path, "Pval_", screen_name, ".csv"))
  
  volcanoplot_screen(p_value_df = output_function_6, selected_hits = selected_hits, screen_name = screen_name, gate_list = gate_list, main_path = main_path)
  plot_individual_guides(df_norm = output_function_2, p_value_df = output_function_6, selected_hits = selected_hits, screen_name = screen_name, gate_list = gate_list, replicate_list = replicate_list, main_path = main_path)
  
}

# ==== Section 2: Running functions on the data used in the manuscript ====

# Genome-scale eGFP-BCL6-FL reporter screen - For identification of gene-losses which prevent BI-3802-induced BCL6 degradation (Fig. 1e)
run_pipeline("/pathname", 
             screen_name = "GW.BCL6FL.deg_H293T_BI3802", 
             gate_list = c("H", "L"), 
             replicate_list = c("Rep1", "Rep2", "Rep3"), 
             n_iterations = 100, 
             n_sgRNAs = 4, 
             selected_hits = c("SIAH1", "SIAH2", "FBXO11"), 
             FDR_thr = 0.15)

# Genome-scale eGFP-BCL6-FL reporter screen - For identification of gene-losses which stabilize BCL6 in the absence of BI-3802 (Extended Data Fig. 2)
run_pipeline("/pathname", 
             screen_name = "201808.GW.BCL6FL.deg_H293T_DMSO",
             gate_list = c("H", "L"), 
             replicate_list = c("Rep1", "Rep2", "Rep3"), 
             n_iterations = 50, 
             n_sgRNAs = 4, 
             selected_hits = c("SIAH1", "SIAH2", "FBXO11"), 
             FDR_thr = 0.15)

# Genome-scale BI-3802 resistance screen in SuDHL4 cells - BI-3802-treated cells compared to DMSO treated cells at day 20 (Fig. 1g)
run_pipeline("/pathname", 
             screen_name = "201812.BCL6.GW_SuDHL4", 
             gate_list = c("BI3802.DAY20", "DMSO.DAY20"), 
             replicate_list = c("Rep1", "Rep2", "Rep3"), 
             n_iterations = 100, 
             n_sgRNAs = 4, 
             selected_hits = c("SIAH1", "SIAH2", "FBXO11"), 
             FDR_thr = 0.01)

# Genome-scale BI-3802 resistance screen in SuDHL4 cells - DMSO-treated cells at day 20 compared to day 5 (Extended Data Fig. 3)
run_pipeline("/pathname", 
             screen_name = "201812.BCL6.GW_SuDHL4", 
             gate_list = c("DMSO.DAY20", "NotTreated.Day5"), 
             replicate_list = c("Rep1", "Rep2", "Rep3"), 
             n_iterations = 100, 
             n_sgRNAs = 4, 
             selected_hits = c("SIAH1", "SIAH2", "FBXO11"), 
             FDR_thr = 0.01)

# Alanine scan of the BCL6-BTB domain in HEK293T cells - For identification of mutants that prvent BI-3802-induced BCL6 degradation (Extended Data Fig. 5)
run_pipeline("/pathname", 
             screen_name = "Ala.BCL6FL_H293T_1uM.BI3802", 
             gate_list = c("H", "L"), 
             replicate_list = c("Rep1", "Rep2", "Rep3"), 
             n_iterations = 1000,
             n_sgRNAs = 8,
             selected_hits = c("WT"), 
             FDR_thr = 0.001)