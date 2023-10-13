# peak_gene_eqtl
total_tissue = c("11_tissues","19_tissues") %>%
  data.frame() %>%
  dplyr::rename(total_tissue = 1) %>%
  arrange()

promoter_list = c("promoter_1","promoter_2","promoter_3") %>%
  data.frame() %>%
  dplyr::rename(promoter = 1) %>%
  arrange()

for (t in 1:(nrow(total_tissue)))
{
  tissue = total_tissue[t,1]
  
  for (p in 1:(nrow(promoter_list)))
  {
    promoter_name = promoter_list[p,1]
    
    if (tissue == "11_tissues") {
      total_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/", promoter_name,"/total_",tissue,"_",promoter_name,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/tissue_peak_gene_eqtl_function.csv")
    } else {
      total_path = paste0("D:/R_analysis/Epigenetics/e_", tissue,"/", tissue,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/tissue_peak_gene_eqtl_function.csv")
    }
    
    total = read.csv(total_path)
    
    ## classify the genetic_function_interpretation to 8 groups
    functional_gene = total %>% 
      select(tissue, PeakID, gene_id, enhancer_interpretation) %>%
      distinct() %>%
      group_by(tissue, PeakID) %>%
      dplyr::mutate(distinct_gene_per_peak = length(unique(gene_id)),
                    R_genetic_enhancer_per_peak = length(enhancer_interpretation[enhancer_interpretation == "enhancer"])/distinct_gene_per_peak,
                    R_genetic_silencer_per_peak = length(enhancer_interpretation[enhancer_interpretation == "silencer"])/distinct_gene_per_peak,
                    R_funtional_region_per_peak = R_genetic_enhancer_per_peak + R_genetic_silencer_per_peak) %>%
      dplyr::mutate(genetic_function_interpretation = case_when(distinct_gene_per_peak == 1 ~ paste0("one-to-one", " & ", enhancer_interpretation),
                                                                distinct_gene_per_peak >= 2 & R_funtional_region_per_peak == 0 ~ paste0("one-to-many", " & ", "gene_independent_undetermined"),
                                                                distinct_gene_per_peak >= 2 & R_funtional_region_per_peak == 1 & R_genetic_enhancer_per_peak == 1 & R_genetic_silencer_per_peak == 0 ~ paste0("one-to-many", " & ", "gene_independent_enhancer"),
                                                                distinct_gene_per_peak >= 2 & R_funtional_region_per_peak == 1 & R_genetic_enhancer_per_peak == 0 & R_genetic_silencer_per_peak == 1 ~ paste0("one-to-many", " & ", "gene_independent_silencer"),
                                                                distinct_gene_per_peak >= 2 & R_funtional_region_per_peak == 1 & R_genetic_enhancer_per_peak > 0 & R_genetic_enhancer_per_peak < 1 & R_genetic_silencer_per_peak > 0 & R_genetic_silencer_per_peak < 1 ~ paste0("one-to-many", " & ", "gene_dependent_functional_region"),
                                                                distinct_gene_per_peak >= 2 & R_funtional_region_per_peak > 0  & R_funtional_region_per_peak < 1  ~ paste0("one-to-many", " & ", "multifaceted_region"))) %>%
      ungroup %>%
      group_by(tissue, gene_id) %>%
      dplyr::mutate(distinct_peak_per_gene = length(unique(PeakID))) %>%
      ungroup 
    write.csv(functional_gene, function_gene_path, row.names = F)
    
    full_function = total %>%
      inner_join(functional_gene, ., by=c("tissue", "PeakID", "gene_id", "enhancer_interpretation"))
    write.csv(full_function, full_function_path, row.names = F)
    
  }
}


# plot

for (t in 1:(nrow(total_tissue)))
{
  tissue = total_tissue[t,1]
  
  for (p in 1:(nrow(promoter_list)))
  {
    promoter_name = promoter_list[p,1]
    
    if (tissue == "11_tissues") {
      total_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/", promoter_name,"/total_",tissue,"_",promoter_name,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/tissue_peak_gene_eqtl_function.csv")
      mylabel = c("Transverse Colon", "Thyroid Gland",                    
                 "Skeletal Muscle", "Sigmoid Colon",  "Prostate Gland",
                 "Pancreas", "Ovary", "Liver",             
                 "Heart",  "Aorta", "Adrenal Gland")
      
      figure_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/Relationship between peak and gene across human tissues.png")
      
    } else {
      total_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/",tissue,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/tissue_peak_gene_eqtl_function.csv")
      mylabel = c("Uterus", "Transverse Colon", "Tibial Nerve", 
                 "Tibial Artery", "Thyroid Gland",  "Testis",
                 "SAT", "Stomach", "Spleen",                     
                 "Skeletal Muscle", "Sigmoid Colon",  "Prostate Gland",
                 "Pancreas", "Ovary", "Liver",             
                 "Heart", "Coronary Artery",  "Aorta", "Adrenal Gland")
      figure_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/Relationship between peak and gene across human tissues.png")
    }
    
    average_gene_tissue = read.csv(function_gene_path) %>%
      dplyr::select(tissue, PeakID, distinct_gene_per_peak) %>%
      distinct() %>%
      group_by(tissue) %>%
      dplyr::mutate(average_gene = sum(distinct_gene_per_peak)/length(unique(PeakID)))%>%
      ungroup 
    
    average_peak_tissue = read.csv(function_gene_path) %>%
      dplyr::select(tissue, gene_id, distinct_peak_per_gene) %>%
      distinct() %>%
      group_by(tissue) %>%
      dplyr::mutate(average_peak = sum(distinct_peak_per_gene)/length(unique(gene_id)))%>%
      ungroup
    
    gene = average_gene_tissue %>%
      dplyr::select(tissue, average_gene) %>%
      distinct() %>% 
      dplyr::mutate(gene_peak = "gene") %>%
      dplyr::rename(values = "average_gene")
    
    peak = average_peak_tissue %>%
      dplyr::select(tissue, average_peak) %>%
      distinct() %>% 
      dplyr::mutate(gene_peak = "peak") %>%
      dplyr::rename(values = "average_peak")
    
    ## find the order
    order <- rev(levels(as.factor(peak$tissue)))
    breaks <- -3:3
    names(breaks) <- abs(breaks)
    
    # back to back bar chart
    p1 = rbind(gene, peak) %>% 
      dplyr::mutate(values = if_else(gene_peak == "peak", -values, values)) %>% 
      ggplot(aes(x = tissue, y = values, 
                 pattern= gene_peak)) +
      geom_bar_pattern(stat = "identity", 
                       width = 0.75,
                       color = "black", 
                       fill = "white",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6) +
      # geom_bar(stat = "identity", width = 0.75) +
      coord_flip() +
      scale_x_discrete(
                       limits = order,
                       labels = mylabel) +
      scale_y_continuous(position = "right",
                         limits = c(-3, 3),
                         breaks = breaks, labels = names(breaks)) +
      labs(x = "", y = "Average Number", pattern = NULL) +
      scale_pattern_manual(values = c(gene = "stripe", peak = "none"),
                           breaks=c("gene", "peak"),
                           labels=c("Corresponding Gene Numbers per Peak",
                                    "Corresponding Peak Numbers per Gene")) +
      theme_classic() +
      theme(legend.position="bottom",
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_text(size=14, face="bold"),
            axis.title.x= element_text(size=20, face="bold"),
            legend.text = element_text(size=10)) + 
      guides(pattern = guide_legend(reverse = TRUE,
                                    nrow=8,
                                    byrow=TRUE)) # reverse the order of items in legend
    
    genetic_function_of_peak_by_tissue = read.csv(function_gene_path) %>%
      dplyr::select(tissue, PeakID, genetic_function_interpretation) %>%
      distinct() %>%
      {table(.$genetic_function_interpretation, .$tissue)} %>%
      as.data.frame.matrix() %>%
      tibble::rownames_to_column("genetic_function_interpretation") %>%
      pivot_longer(-1,
                   names_to = "tissue", 
                   values_to = "values") %>%
      dplyr::mutate(genetic_function_interpretation = factor(genetic_function_interpretation,
                                                             levels= c("one-to-one & enhancer",
                                                                       "one-to-one & silencer", 
                                                                       "one-to-one & undetermined",
                                                                       "one-to-many & gene_independent_enhancer", 
                                                                       "one-to-many & gene_independent_silencer", 
                                                                       "one-to-many & gene_dependent_functional_region",
                                                                       "one-to-many & gene_independent_undetermined",
                                                                       "one-to-many & multifaceted_region"))) 
    
    levels(genetic_function_of_peak_by_tissue$genetic_function_interpretation)
    
    p2 = genetic_function_of_peak_by_tissue %>% 
      ggplot(aes(fill = genetic_function_interpretation,  
                 x = values, 
                 y = tissue)) + 
      geom_bar(position = "stack", stat = "identity", width = 0.75) + 
      # position = "fill" will present as the percentages
      scale_fill_manual(values = c("brown1",
                                   "dodgerblue1",
                                   "lightgoldenrod1",
                                   "brown3", 
                                   "dodgerblue3", 
                                   "lightgoldenrod3",
                                   "gray70",
                                   "gray90"),
                        # "hotpink3", "deepskyblue3", "lightgoldenrod3"
                        # "dodgerblue2", high = "brown2",
                        breaks = c("one-to-one & enhancer", 
                                   "one-to-one & silencer", 
                                   "one-to-one & undetermined",
                                   "one-to-many & gene_independent_enhancer",
                                   "one-to-many & gene_independent_silencer",
                                   "one-to-many & gene_dependent_functional_region",
                                   "one-to-many & gene_independent_undetermined",
                                   "one-to-many & multifaceted_region"),
                        labels = c("One-to-One & Enhancer", 
                                   "One-to-One & Silencer", 
                                   "One-to-One & Undetermined",
                                   "One-to-Many & Gene Independent Enhancer",
                                   "One-to-Many & Gene Independent Silencer",
                                   "One-to-Many & Gene Dependent Functional Region",
                                   "One-to-Many & Gene Independent Undetermined",
                                   "One-to-Many & Multifaceted Region")) +
      labs(x = "Peak Numbers", y = "", fill=NULL) + 
      scale_x_continuous(position = "top") +
      scale_y_discrete(limits = rev(levels(as.factor(genetic_function_of_peak_by_tissue$tissue)))) +
      theme_classic() +
      theme(legend.position="bottom",
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_blank(),
            axis.title.x= element_text(size=20, face="bold"),
            legend.text = element_text(size=10)) + 
      guides(fill=guide_legend(nrow=8,
                               byrow=TRUE))
    
    # "#d9c980", "#8f809f"/" #a889ab" 
    # lightpink2 steelblue2
    # bisque3,  darkseagreen4/ #617558
    # p = grid.arrange(p1,p2,p3,p4,p5, ncol=5, widths=c(0.37,0.25,0.25,0.25,0.25)) 
    # ggsave("overall epigenetic functionality across 19 human tissues.pdf", plot = p,width = 65, height = 30, units = "cm")
    
    p = grid.arrange(p1,p2, ncol=2, widths=c(0.5,0.5)) 
    ggsave(figure_path, plot = p,width = 45, height = 30, units = "cm")
    
  }
}

# select gene_dependent_functional_region

total_tissue = c("11_tissues","19_tissues") %>%
  data.frame() %>%
  dplyr::rename(total_tissue = 1) %>%
  arrange()

promoter_list = c("promoter_1","promoter_2","promoter_3") %>%
  data.frame() %>%
  dplyr::rename(promoter = 1) %>%
  arrange()

for (t in 1:(nrow(total_tissue)))
{
  tissue = total_tissue[t,1]
  
  for (p in 1:(nrow(promoter_list)))
  {
    promoter_name = promoter_list[p,1]
    
    if (tissue == "11_tissues") {
      total_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/", promoter_name,"/total_",tissue,"_",promoter_name,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/tissue_peak_gene_eqtl_function.csv")
      gene_dependent_functional_region_path = paste0("D:/R_analysis/Epigenetics/g_",tissue,"/",promoter_name,"/gene_dependent_functional_region.csv")
    } else {
      total_path = paste0("D:/R_analysis/Epigenetics/e_", tissue,"/", tissue,".csv")
      function_gene_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/functional_gene.csv")
      full_function_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/tissue_peak_gene_eqtl_function.csv")
      gene_dependent_functional_region_path = paste0("D:/R_analysis/Epigenetics/e_",tissue,"/gene_dependent_functional_region.csv")
    }
    
    gene_dependent_functional_region = read.csv(full_function_path) %>%
      filter(genetic_function_interpretation == "one-to-many & gene_dependent_functional_region") %>%
      dplyr::select(tissue,genetic_function_interpretation,PeakID,Chr,Start,End,Annotation,Gene_Description,Gene_Type,
                    distinct_gene_per_peak, Peak_Center_gene_tss, gene_id,gene_name, enhancer_interpretation, distinct_selected_rs_id) %>%
      distinct()
    
    write.csv(gene_dependent_functional_region, gene_dependent_functional_region_path, row.names = F)
    
  }
}
