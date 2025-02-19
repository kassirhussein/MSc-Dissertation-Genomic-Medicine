# COrrelations ------------------------------------------------------------


# Gene/DX + visualisation -------------------------------------------------




genopheno_genes_disease <- genopheno_genes %>%
  group_by( `Diagnosed Gene` , Disease) %>%
  tally() %>%
  arrange(-n) %>%
  inner_join(df_disease, by = c("Diagnosed Gene" = "Diagnosed Gene")) %>%
  mutate(percentage_disease = (n.x/n.y)*100) %>%
  select(1,2,6)



data <- data.frame(genopheno_genes_disease)

# for loop Function to plot diseases for each gene


plot_diseases_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Disease, -percentage_disease), x = percentage_disease, fill = Disease)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Disease Occurrences for Gene:", gene_data$Diagnosed.Gene[1]),
         x = "Disease",
         y = "Count (n)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Splitting the data by Gene and plotting each separately

unique_genes <- unique(data$Diagnosed.Gene)

for (gene in unique_genes) {
  gene_data <- filter(data, Diagnosed.Gene == gene)
  print(plot_diseases_by_gene(gene_data))
}




######For loop visulaisation of multiple Genes.




data <- data.frame(hpo_percentage_all)

data <- data %>%
  arrange(Diagnosed.Gene)

plot_hpo_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(HPO.terms.s. , percentage_hpo_all ) ,x = percentage_hpo_all, fill = HPO.terms.s.)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("HPO term Occurrences for: ", gene_data$Diagnosed.Gene[1]),
         x = "Frequency (%)",
         y = "HPO term") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Splitting the data by Gene and plotting each separately

unique_genes <- unique(data$Diagnosed.Gene)

for (gene in unique_genes) {
  gene_data <- filter(data, Diagnosed.Gene == gene)
  print(plot_hpo_by_gene(gene_data))
  
  
}
