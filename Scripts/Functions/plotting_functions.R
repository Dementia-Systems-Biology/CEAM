make_stacked_bar <- function(cohort1_list,
                             cohort2_list,
                             cell_types,
                             filename = NULL,
                             cohort1_name = "BDR",
                             cohort2_name = "UKBBN2",
                             shared_cpg) 
  {
  library(tidyverse)
  library(ggplot2)
  
  # pre-define names
  unique1_name <- cohort1_name
  unique2_name <- cohort2_name
  
  
    
  # Helper to compute counts for a single cell type
  get_cpg_counts <- function(cell_type) {
    cpgs1 <- cohort1_list[[cell_type]]
    cpgs2 <- cohort2_list[[cell_type]]
    
    if (is.list(cpgs1) & is.list(cpgs2)) {
      shared  <- unlist(unique(c(
        intersect(unlist(cpgs1[1:3]), unlist(cpgs2[1:3])), 
                          intersect(intersect(unlist(cpgs1[4]), unlist(cpgs2[4])), shared_cpg)
        )))
                        
      unique1 <- unlist(unique(c(setdiff(unlist(cpgs1[1:3]), unlist(cpgs2[1:3])), setdiff(unlist(cpgs1[4]), shared_cpg))))
      unique2 <- unlist(unique(c(setdiff(unlist(cpgs1[1:3]), unlist(cpgs2[1:3])), setdiff(unlist(cpgs2[4]), shared_cpg))))
      
    } else {
      
      shared  <- intersect(cpgs1, cpgs2)
      unique1 <- setdiff(cpgs1, cpgs2)
      unique2 <- setdiff(cpgs2, cpgs1)
    }

      

    df <- tibble(
      Cohort = c("Both", cohort1_name, cohort2_name),
      CellType = cell_type,
      CpG_Type = c("Shared", unique1_name, unique2_name),
      Count = c(length(shared), length(unique1), length(unique2)),
      Coverage = paste("(",round((c(length(shared)/778549, length(unique1)/778549, length(unique2)/778549) *100), 1), "%", ")", sep = "") # same as count but divided by all probes in an .idat file
    )
    df$count_with_percentage <- paste(df$Count, df$Coverage, sep = "\n")
    return(df)
  }
  
  # combine results across cell types
  overlap_df <- map_dfr(cell_types, get_cpg_counts) %>%
    mutate(
      label_color = ifelse(CpG_Type == unique1_name, "black", "white"),
      prop = Count / sum(Count),
      vjust_val = ifelse((prop < 0.02 | Count < 3000) & CpG_Type == unique1_name, -0.6, 0.5)
    )
  
  # change factor levels for consistent fill order
  overlap_df$CpG_Type <- factor(
    overlap_df$CpG_Type,
    levels = c(unique1_name, unique2_name, "Shared")
  )
  

  
  # create stacked bar chart
  stacked_bar <- ggplot(overlap_df, aes(x = CellType, y = Count, fill = CpG_Type)) +
    geom_col(position = "stack") +
    geom_text(
      aes(label = count_with_percentage, color = label_color, vjust = vjust_val),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    scale_fill_manual(
      values = setNames(
        c("#DDAA33", "#BB5566", "#004488"),
        c(cohort1_name,
          cohort2_name,
          "Shared")
      )
    ) +
    labs(
      x = "Cell Type",
      y = "CpG Count",
      fill = "Cohort"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 12)
    )
  
  # (optional) save the plot to reuse and change later
  if (!is.null(filename)){
    saveRDS(stacked_bar, file = filename)
    message("Figure saved to: ", filename)
  }

  
  # return both the dataframe and plot
  list(data = overlap_df, plot = stacked_bar)
}

# function to create UpSet plots from cell type enrichment results
plot_CpG_UpSet <- function(ORA_result, min_set_size = 3, num_breaks = 3){
  # Prepare input list
  tmp <- strsplit(ORA_result$enriched_CpGs, split = ",")
  # change the names to the abbreviations
  names(tmp) <- sapply(rownames(ORA_result), function(x) {
    if (grepl("NeuN", x)) {
      "NEU"
    } else if (grepl("IRF8", x)) {
      "MG"
    } else if (grepl("SOX10", x)) {
      "OLIG"
    } else if (grepl("TN", x)) {
      "AST"
    } else {
      x  # fallback to original name if no match
    }
  })
  
  list_input <- UpSetR::fromList(tmp)
  colnames(list_input) <- names(tmp)
  rownames(list_input) <- unique(unlist(tmp))
  list_input <- as.data.frame(list_input == 1)
  
  # Plot
  ComplexUpset::upset(
    list_input,
    names(tmp),
    name = "CpG sets",
    base_annotations=list(
      'Intersection size'=intersection_size(counts=FALSE) # change to TRUE if values should be displayed
    ),
    min_size = min_set_size,
    set_sizes = (
      upset_set_size() +
        scale_y_reverse(
          breaks = scales::pretty_breaks(n = num_breaks) # makes nicer breaks for visualizing total overlap size
        ) +
        theme(axis.text.x = element_text(angle = 90))
    ),
    width_ratio = 0.3)
}

# function to manually annotate UpSet plots
annotate_sig <- function(upset, sig_label) {
  #' @param upset the UpSet plot object from ComplexUpset to be annotated
  #' @param sig_label a vector of characters used to indicate significance in the barplot, commonly "*"
  pbuilt <- ggplot_build(upset[[3]])
  bar_data <- pbuilt$data[[2]]  # Usually the first layer is the bars
  upset[[3]] <- upset[[3]] +  annotate("text", x = as.numeric(bar_data[,"x"]) -0.25, y = bar_data[,"count"] + 0.25*bar_data[,"count"], label = sig_label, size = 6)
  return(upset)
}

# function to quickly plot all the line plots with ribbons
plot_with_ribbon <- function(data,
                             x,
                             y,
                             ymin,
                             ymax,
                             group = NULL,
                             title = "",
                             x_label = "nr of DMPs",
                             y_label = "Odds Ratio/P-value",
                             colors = c("#DDAA33", "#004488"),
                             y_u_lim = 1) {
  
  if (!is.null(group)) {
    p <- ggplot(data, aes(x = !!sym(x),
                          y = !!sym(y),
                          color = !!sym(group),
                          fill = !!sym(group))) +
      geom_line(linewidth = 1.0) +
      geom_ribbon(aes(ymin = !!sym(ymin), ymax = !!sym(ymax)), colour = NA,
                  alpha = 0.2) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors)
  } else {
    p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
      geom_line(color = colors[1], linewidth = 1.0) +
      geom_ribbon(aes(ymin = !!sym(ymin), ymax = !!sym(ymax)),
                  fill = colors[1], alpha = 0.2)
  }
  
  # Shared styling
  p +
    theme_minimal(base_size = 8) +
    ylim(0, y_u_lim) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(
      axis.title = element_text(size = 7),
      title = element_text(size = 8),
      legend.title = element_text(size = 7),   # smaller legend title
      legend.text  = element_text(size = 6),   # smaller legend labels
      legend.key.size = unit(0.6, "lines")     # smaller legend boxes (optional)
    )
}

# function to apply colors to a list of plots generated from the simulations
apply_manual_color <- function(plot_list, color_values) {
  lapply(plot_list, function(row) {
    lapply(row, function(plot) {
      plot +
        scale_color_manual(values = color_values, name = "Correction",
                           labels = c("Uncorrected", "Corrected")) +
        scale_fill_manual(values = color_values, name = "Correction",
                          labels = c("Uncorrected", "Corrected")) +
        theme(legend.position = "none")
    })
  })
}

# plotting function for the OR permutation plots
plot_OR_dist_loop <- function(df, label_col = "label", value_col = "OR", bins = 20, ORA_res) {
  # change the ORA_res names so it is the same as the rest of the plotting function
  ORA_df <- data.frame(
    label = rownames(ORA_res),
    OR_value = ORA_res$Odds.ratio
  )
  
  # list to store individual plots
  plots <- list()
  
  # loop over each cell type
  for (lab in unique(df[[label_col]])) {
    # subset the data for each cell type
    df_sub <- df[df[[label_col]] == lab, ]
    OR <- ORA_df$OR_value[ORA_df$label == lab]
    
    # create the histogram
    p <- ggplot(df_sub, aes(x = .data[[value_col]])) +
      geom_histogram(bins = bins, fill = "#004488", color = "white", alpha = 0.8) +
      #geom_density(fill = "#004488", alpha = 0.5, color = "#004488") +
      geom_vline(xintercept = OR, color = "#BB5566", linetype = "dashed", size = 1) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste0(lab),
        x = "Odds Ratio (OR)",
        y = "Count"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold")
      )
    
    # Save plot in the list
    plots[[lab]] <- p
  }
  
  return(plots)
}