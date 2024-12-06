
#Fisher test for each gene/hpo term correlations
######Fisher test with OR#



genopheno_all_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

genes <- unique(genopheno_all_separated$`Diagnosed Gene`)
hpo_terms <- unique(genopheno_all_separated$`HPO terms(s)`)



results_fisher_hpoterm <- data.frame(Gene = character(), 
                                     HPO_term = character(), 
                                     p_value = numeric(), 
                                     odds_ratio = numeric(),
                                     stringsAsFactors = FALSE)

# Loop through each gene and HPO term to perform Fisher's test and calculate Odds Ratio
for (gene in genes) {
  for (hpo in hpo_terms) {
    
    
    A <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` == gene & genopheno_all_separated$`HPO terms(s)` == hpo])
    B <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` == gene & genopheno_all_separated$`HPO terms(s)` != hpo])
    C <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` != gene & genopheno_all_separated$`HPO terms(s)` == hpo])
    D <- sum(genopheno_all_separated$n[genopheno_all_separated$`Diagnosed Gene` != gene & genopheno_all_separated$`HPO terms(s)` != hpo])
    
    contingency_table <- matrix(c(A, B, C, D), nrow = 2)
    
    
    fisher_test <- fisher.test(contingency_table)
    
    
    odds_ratio <- ifelse(!is.null(fisher_test$estimate), fisher_test$estimate, NA)
    
    
    results_fisher_hpoterm <- rbind(results_fisher_hpoterm, data.frame(`Diagnosed Gene` = gene, 
                                                                       `HPO terms(s)` = hpo, 
                                                                       p_value = fisher_test$p.value, 
                                                                       odds_ratio = odds_ratio))
  }
}

# Display the results
print(results_fisher_hpoterm)

results_fisher_hpoterm$p_adjusted <- p.adjust(results_fisher_hpoterm$p_value, method = "BH")

results_fisher_hpoterm <- results_fisher_hpoterm %>%
  filter(p_value<0.05)%>%
  filter(odds_ratio !="Inf")