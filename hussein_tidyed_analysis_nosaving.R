#################
.libPaths(c( .libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
library(tidyverse)



# import hpo file #---------------------------------------------------------



setwd("/nas/weka.gel.zone/re_gecip/shared_allGeCIPs/dsmedley/FUSIL/GenoPheno")

hpo1 = read_delim("hpo_ontology.txt")


#imort genopheno file

genopheno = read_delim("extract_allelic_series_30_08_24_plus_denovo.tsv") %>%
  unique()

genopheno <- genopheno %>%
  mutate(Disease = ifelse(Disease == "intellectual disability", "Intellectual disability", Disease))

genopheno <- genopheno %>%
  mutate(Disease = ifelse(Disease == "intellectual-disability", "Intellectual disability", Disease))

phenotype_300824 = read_delim("phenotype_300824.hpoa")

#Import annotations

setwd("/nas/weka.gel.zone/re_gecip/shared_allGeCIPs/dsmedley/FUSIL/GenoPheno/annotations")

panelapp_genes = read_delim("panelapp_by_gene_panel_aus_data_accessed_190824.txt")

gene_to_phenotype = read_delim("genes_to_phenotype_300824.txt")




dim(genopheno)
length(unique(genopheno$`GEL family ID`))

genopheno_dups = genopheno$`GEL family ID`[duplicated(genopheno$`GEL family ID`)]
genopheno_dups_check = genopheno %>%
  filter(`GEL family ID` %in% unique(genopheno_dups))

gene_dups <- genopheno_dups_check %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

### there are family ids with more than one diagnosed gene/variant


### 5735 cases in total

gene_consistency <- genopheno %>%
  group_by(`GEL family ID`) %>%
  summarise(Same_Gene = n_distinct(`Diagnosed Gene`) == 1)

##### There 142 GEL family ID that have different genes as Diagnosed Gene


gene_occurrence <- genopheno %>%  ## to export
  group_by(`GEL family ID`) %>%
  summarise(`Diagnosed Gene` = paste(unique(`Diagnosed Gene`), collapse = ", "))



##### gene_occurence summarise the genes for each GEL family id

genopheno %>%
  count(Disease == "Intellectual disability")

all_disease <-genopheno %>%
  group_by(Disease) %>%
  tally

### after correction 1313 cases of ID as disease, 

all_genes <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n)



id_genes <- genopheno %>%
  group_by(Disease,`Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n) %>%
  filter(Disease == "Intellectual disability") %>%
  filter(n >=10)




### 420 genes with ID as disease.

### when filtering >10 cases, we are left with 22 genes that have ID as disease




genopheno_genes <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$ `Diagnosed Gene`)

genopheno_familyid <- genopheno_genes %>%
  group_by(`GEL family ID`) %>%
  tally()%>%
  select(1)



#### we have 437 cases that have these 22 genes as diagnosed genes

genopheno_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease  == "Intellectual disability")

genopheno_no_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(!Disease  == "Intellectual disability")



### 345 of these cases have ID as disease


df_disease <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally()%>%
  arrange(-n)%>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n/sum(n)*100))

sum(df_disease$n) ###437 cases


ggplot(data=df_disease, aes(x= percentage, y = reorder(`Diagnosed Gene`,percentage) ,fill=`Diagnosed Gene`)) +
  geom_bar(stat = "identity", position = "dodge")+labs(title="Percenatge of ID genes",x="percentage", y= "Diagnosed_Gene")+
  theme(legend.position = "none")
##################################################


ggplot(data=df_disease, aes(x= n, y = reorder(`Diagnosed Gene`,percentage) ,fill=`Diagnosed Gene`)) +
  geom_bar(stat = "identity", position = "dodge")+
  labs(title="Number of diagnoses per gene (genes > 10 ID diagnoses)",x="number of diagnoses", y= "diagnosed sene") +
  theme(legend.position = "none")



##################################################

#### filtering for ID only###

df_disease_id <- genopheno %>%
  group_by( Disease, `Diagnosed Gene`) %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  tally()%>%
  arrange(-n)%>%
  mutate(n_id=n)%>%
  select(1,2,4)

ggplot(data=df_disease_id, aes(x= n_id , y = reorder(`Diagnosed Gene` , n_id), fill= `Diagnosed Gene`))+
  geom_bar(stat = "identity", position = "dodge")+labs(title="Number of patients with ID in each gene",x="Count", y= "Diagnosed Gene")

# Every other disease


df_disease_od <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease != "Intellectual disability") %>%
  group_by(`Diagnosed Gene`) %>%
  tally()%>%
  arrange(-n)%>%
  mutate(n_other_disease=n)%>%
  select(1,3)



#### combing the two data frames to compare
#### 2 visualzitiona methods


compare_table_genes <- df_disease_id %>%
  left_join(df_disease_od , by = c("Diagnosed Gene" = "Diagnosed Gene"))

data_long_genes <- compare_table_genes %>%
  pivot_longer(cols = c(n_id,n_other_disease),
               names_to = "Type",
               values_to = "Count")



ggplot(data = data_long_genes, aes(y = `Diagnosed Gene`, x = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~Type)


###################################################



# MOI overall###

genopheno_moi <- genopheno %>%
  group_by(MOI) %>%
  tally()

genopheno_moi <- genopheno_genes %>%
  group_by(MOI) %>%
  tally()


# NOw to see MOI for each gene:

genopheno_genes_moi <- genopheno_genes %>%
  group_by( `Diagnosed Gene` , MOI) %>%
  tally() %>%
  arrange(-n) %>%
  inner_join(df_disease, by = c("Diagnosed Gene" = "Diagnosed Gene")) %>%
  mutate(percentage_zygosity = (n.x/n.y)*100) %>%
  select(1,2,6)



ggplot(data=genopheno_genes_moi) +
  geom_bar(mapping = aes(y= `Diagnosed Gene` , x =  percentage_zygosity, fill=MOI ), stat = "identity", position = "dodge")



##### Variant overall

variant_all <- genopheno %>%
  group_by(`Diagnosed Variant(s)`) %>%
  tally()

variant_all <- genopheno_genes %>%
  group_by(`Diagnosed Variant(s)`) %>%
  tally()

variant_consequence <- genopheno %>%
  group_by(`Consequence(s)`) %>%
  tally()


denovo<- genopheno %>%
  group_by(`De novo`) %>%
  tally()

denovo<- genopheno_genes %>%
  group_by(`De novo`) %>%
  tally()





# COrrelations ------------------------------------------------------------


# Gene/DX + visualisation -------------------------------------------------




genopheno_genes_disease <- genopheno_genes %>%
  group_by( `Diagnosed Gene` , Disease) %>%
  tally() %>%
  arrange(-n) %>%
  inner_join(df_disease, by = c("Diagnosed Gene" = "Diagnosed Gene")) %>%
  mutate(percentage_disease = (n.x/n.y)*100) %>%
  select(1,2,6)

#export_1


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


##  those with ID only
## and those with ID + other disease categories

genes_only_id <- genopheno_genes_disease %>%
  group_by(`Diagnosed Gene`) %>%
  summarise(all_diseases = paste(unique(Disease), collapse = ", ")) %>%
  filter(all_diseases == "Intellectual disability")

### Seems only FOXP1/SHANK3/WAC have intellectual disability only as disease

#### these function is great!!
### need to highlight hte more pleiotropic genes

# HPO terms work ----------------------------------------------------------



#MOst occuring HPO terms between genes#

most_occuring_hpoterms <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`HPO terms(s)`, `Diagnosed Gene`) %>%
  tally()


most_occuring_hpoterms <- most_occuring_hpoterms %>%
  group_by(`HPO terms(s)`)%>%
  tally() %>%
  mutate(percentage_top_hpo = (n/22)*100)%>%
  filter(percentage_top_hpo>50)



# #####Gene/HPO term###### ------------------------------------------------



## HPO terms in each Gene occurence##


genopheno_all_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)


df_disease <- genopheno %>%
  group_by(`Diagnosed Gene`) %>%
  tally()%>%
  arrange(-n)%>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n/sum(n)*100))



hpo_percentage_all <- genopheno_all_separated %>%
  left_join(df_disease, by = "Diagnosed Gene") %>%
  mutate(percentage_hpo_all = (n.x/n.y)*100) %>%
  select(1,2,6)

###filterin for the hpo term occuring in >30% of patietns for each gene

hpo_percentage_all <- hpo_percentage_all %>%
  filter(percentage_hpo_all>30)




hpoterm_gene_all <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by (`HPO terms(s)`, `Diagnosed Gene`) %>%
  tally() %>%
  arrange(-n)


hpoterm_gene_all <- hpoterm_gene_all %>%
  left_join(df_disease, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x/n.y)*100) %>%
  select(1,2,6) %>%
  filter(Percentage>20)

hpoterm_gene_all_2 <- hpoterm_gene_all %>%
  group_by(`HPO terms(s)`)%>%
  tally()



##############NOw to see for each gene of the 24 genes, it occurs with what hpo term####
######### visualization


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

##### theser are great too 
##### but you need to specify these are percentages
##### are not counts, and how the percentages were computed
##### total number of diagnoses in that gene?????

################# Visualize each hpo term in what genes it occurs###


data <- data.frame(hpoterm_gene_all)

data <- data %>%
  arrange(Diagnosed.Gene)

plot_gene_by_hpo <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Diagnosed.Gene , -Percentage ) ,x = Percentage, fill = Diagnosed.Gene)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Gene Observed in HPO term :", gene_data$HPO.terms.s.[1]),
         x = "Percentage of Diagnoses in the Gene",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Splitting the data by Gene and plotting each separately

unique_hpoterms <- unique(data$HPO.terms.s.)

for (hpoterm in unique_hpoterms) {
  gene_data <- filter(data, HPO.terms.s. == hpoterm)
  print(plot_gene_by_hpo(gene_data))
}


#### here you need to decide which terms to plot
#### delayed motor development is a good one

# Matching for top level  -------------------------------------------------


#################

genopheno_all_separated2 <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo <- genopheno_all_separated2 %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1,2,17) %>%
  unique()


hpo_toplevel_df <- joined_hpo %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n)


df_disease2 <- genopheno %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`)%>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1,2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`)%>%
  tally()


hpo_toplevel_df <- hpo_toplevel_df %>%
  left_join(df_disease2, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x/n.y)*100)%>%
  filter(Percentage>20)

#Work below is to visualioze each gene and the top level hpo term occurence



data3 <- data.frame(hpo_toplevel_df)

plot_hpotoplevel_by_gene <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(hpo_ancestors_description ,  Percentage ) ,x = Percentage, fill = hpo_ancestors_description)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("HPO term Occurrences for:", gene_data$Diagnosed.Gene[1]),
         y = "HPO term top Level",
         x = "Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Splitting the data by Gene and plotting each separately

unique_genes3 <- unique(data3$Diagnosed.Gene)

for (gene in unique_genes3) {
  gene_data <- filter(data3, Diagnosed.Gene == gene)
  print(plot_hpotoplevel_by_gene(gene_data))
}

#####
example_ancestor <- joined_hpo %>%
  group_by(hpo_ancestors_description, `Diagnosed Gene`) %>%
  tally() %>%
  left_join(df_disease2, by = "Diagnosed Gene") %>%
  mutate(Percentage = (n.x/n.y)*100)%>%
  filter(Percentage>20)


data <- data.frame(example_ancestor)

data <- data %>%
  arrange(Diagnosed.Gene)

plot_gene_by_hpo <- function(gene_data) {
  ggplot(gene_data, aes(y = reorder(Diagnosed.Gene , -Percentage ) ,x = Percentage, fill = Diagnosed.Gene)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = paste("Genes Observed in HPO Ancestor term :", gene_data$hpo_ancestors_description[1]),
         x = "Percentage of Diagnoses in the Gene",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d(option = "plasma")
}

# Splitting the data by Gene and plotting each separately

unique_hpoterms <- unique(data$hpo_ancestors_description)

for (hpoterm in unique_hpoterms) {
  gene_data <- filter(data, hpo_ancestors_description == hpoterm)
  print(plot_gene_by_hpo(gene_data))
}




combined_df <- example_ancestor%>%
  group_by(hpo_ancestors_description) %>%
  summarize(`Diagnosed Gene` = paste(`Diagnosed Gene`, collapse = ", "))


nb_df <- example_ancestor %>%
  group_by(hpo_ancestors_description) %>%
  tally()

#######################################

# phenotype enrichmenetin genes with or without ID ------------------------



##### HPO term in ID as disease

genopheno_id_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease  == "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

df_disease_id <- genopheno %>%
  group_by( Disease, `Diagnosed Gene`) %>%
  tally()%>%
  arrange(-n)%>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease == "Intellectual disability") %>%
  mutate(percentage_id = (n/sum(n)*100))


hpo_percentage_id <- genopheno_id_separated %>%
  left_join(df_disease_id, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_Id = (n.x/n.y)*100) %>%
  select(1,2,4,7)%>%
  filter(Percentage_hpo_Id>10)


#### I think here is ok, because in df_disease_id, even when you
#### compute the percentages, you are not using them

#export_2

############ HPO terms filtering  ID disease ID OUT


genopheno_no_id_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(Disease  != "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`) %>%
  tally() %>%
  arrange(-n)

df_disease_no_id <- genopheno %>%
  filter(Disease != "Intellectual disability") %>%
  group_by(`Diagnosed Gene`) %>%
  tally()%>%
  arrange(-n)%>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  mutate(percentage = (n/sum(n)*100))



hpo_percentage_no_id <- genopheno_no_id_separated %>%
  left_join(df_disease_no_id, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_Other_Diseases = (n.x/n.y)*100) %>%
  select(1,2,6)%>%
  filter(Percentage_hpo_Other_Diseases>10)

#export_3


########################
#Compare hpo term when  for ID patients and non-intellectual disability
########################

compare_table_hpoterm <- full_join( hpo_percentage_id,hpo_percentage_no_id, by = c("Diagnosed Gene", "HPO terms(s)")) %>%
  select(1,2,4,5)%>%
  arrange(`Diagnosed Gene`)

df_compare_table_hpoterm <- data.frame(compare_table_hpoterm)
#to visualise the comparasion tbale 

data_long_genes_hpoterm <- df_compare_table_hpoterm %>%
  pivot_longer(cols = c(Percentage_hpo_Other_Diseases, Percentage_hpo_Id),
               names_to = "Type",
               values_to = "Percentage")

genes <- unique(df_compare_table_hpoterm$Diagnosed.Gene)

for (gene in genes) {
  
  gene_data <- subset(data_long_genes_hpoterm, Diagnosed.Gene == gene)
  
  data_long_genes_hpoterm <- df_compare_table_hpoterm %>%
    pivot_longer(cols = c(Percentage_hpo_Other_Diseases, Percentage_hpo_Id),
                 names_to = "Type",
                 values_to = "Percentage")
  
  
  p <- ggplot(gene_data, aes(y = reorder (HPO.terms.s. , -Percentage), x= Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
    scale_fill_manual(values = c("Percentage_hpo_Other_Diseases" = "blue", "Percentage_hpo_Id" = "red")) +
    labs(title = paste("HPO Terms for", gene),
         y = "HPO term",
         x = "Percentage",
         fill = "Percentage Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display the plot
  print(p)
}





# Gene Ancestor terms in ID and Other diseases ----------------------------



##
#Matching for top level ID
##

genopheno_all_separated_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(`Disease` == "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo_id <- genopheno_all_separated_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1,2,17) %>%
  unique()

df_disease_id2 <-genopheno_all_separated_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")%>%
  select(1,2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`) %>%
  tally()



hpo_toplevel_id <- joined_hpo_id %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n)%>%
  left_join(df_disease_id2, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_ancetor_ID = (n.x/n.y)*100) %>%
  filter(Percentage_hpo_ancetor_ID>10)


##
#Matching for top level other diseases
#

genopheno_all_separated_no_id <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  filter(`Disease` != "Intellectual disability") %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo_no_id <- genopheno_all_separated_no_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")%>%
  select(1,2,17) %>%
  unique()

df_disease_no_id2 <-genopheno_all_separated_no_id %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")%>%
  select(1,2) %>%
  unique() %>%
  group_by(`Diagnosed Gene`) %>%
  tally()

hpo_toplevel_no_id <- joined_hpo_no_id %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  left_join(df_disease_no_id2, by = "Diagnosed Gene") %>%
  mutate(Percentage_hpo_ancetor_Other_diseases = (n.x/n.y)*100)%>%
  filter(Percentage_hpo_ancetor_Other_diseases>10)



##### compare hpoterms in ID and Other diseases

compare_table_hpoterm_2 <- full_join( hpo_toplevel_id,hpo_toplevel_no_id ,by = c("Diagnosed Gene", "hpo_ancestors_description")) %>%
  arrange(`Diagnosed Gene`)


df_compare_table_hpoterm_2 <- data.frame(compare_table_hpoterm_2)
#to visualise the comparasion tbale 

data_long_genes_hpoterm_2 <- df_compare_table_hpoterm_2 %>%
  pivot_longer(cols = c(Percentage_hpo_ancetor_ID, Percentage_hpo_ancetor_Other_diseases),
               names_to = "Type",
               values_to = "Frequencies")

genes <- unique(df_compare_table_hpoterm_2$Diagnosed.Gene)

for (gene in genes) {
  
  gene_data <- subset(data_long_genes_hpoterm_2, Diagnosed.Gene == gene)
  
  data_long_genes_hpoterm_2 <- df_compare_table_hpoterm_2 %>%
    pivot_longer(cols = c(Percentage_hpo_ancetor_ID, Percentage_hpo_ancetor_Other_diseases),
                 names_to = "Type",
                 values_to = "Frequencies")
  
  
  p <- ggplot(gene_data, aes(y = reorder (hpo_ancestors_description , -Frequencies), x= Frequencies, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
    scale_fill_manual(values = c("Percentage_hpo_ancetor_ID" = "green4", "Percentage_hpo_ancetor_Other_diseases" = "orange")) +
    labs(title = paste("HPO Terms for", gene),
         x = "Frequencies %",
         y = "HPO Ancestor Term",
         fill = "Frequencies Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display the plot
  print(p)
}





####################################

# MOI and HPO term --------------------------------------------------------



##
genopheno_separated <- genopheno_genes %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")



id_genes_moi <- genopheno %>%
  group_by(`Diagnosed Gene`, MOI) %>%
  tally() %>%
  arrange(-n) %>%
  filter(`Diagnosed Gene` %in% c("SYNGAP1","SHANK3","SETD5","SCN2A","MED13L","MECP2","KDM5C","DDX3X","ANKRD11"))


grouped_data_disease_moi <- genopheno_separated%>%
  group_by(`Diagnosed Gene`, `HPO terms(s)`, MOI) %>%
  summarize(count = n(), .groups = 'drop') %>%
  filter(`Diagnosed Gene` %in% c("SYNGAP1","SHANK3","SETD5","SCN2A","MED13L","MECP2","KDM5C","DDX3X","ANKRD11"))


grouped_data_disease_moi <- grouped_data_disease_moi %>%
  full_join(id_genes_moi, by= c("Diagnosed Gene", "MOI"))%>%
  mutate(Percentage = (count/n)*100)%>%
  select(1,2,3,6)%>%
  filter(`Diagnosed Gene` %in% c("KDM5C", "MECP2", "SCN2A"))%>%
  filter(Percentage>15)




grouped_data_disease_moi <- data.frame(grouped_data_disease_moi)
## plotting each gene


unique_genes <- unique(grouped_data_disease_moi$`Diagnosed.Gene`)

for (gene in unique_genes) {
  
  gene_data <- subset(grouped_data_disease_moi, `Diagnosed.Gene` == gene)
  
  
  plot <- ggplot(gene_data, aes(y = `HPO.terms.s.`,x=Percentage, fill = MOI)) +
    geom_bar(stat="identity",position = "dodge") +
    labs(title = paste("Phenotypes Occurrence for", gene), x = "Percentage", y = "HPO terms") +
    theme_minimal()
  
  print(plot)
  
}



# Variant type and hpo term analysis --------------------------------------



variant_consequence <- genopheno %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`)%>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")%>%
  select(1,2,10,17) %>%
  unique() %>%
  group_by(`Diagnosed Gene`, `Consequence(s)`, `hpo_ancestors_description`) %>%
  tally()


df_variant <- genopheno %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`)%>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-") %>%
  select(1,2,10) %>%
  unique() %>%
  group_by(`Diagnosed Gene`, `Consequence(s)`)%>%
  tally()


variant_hpoterm <- variant_consequence %>%
  left_join(df_variant, by=c("Diagnosed Gene", "Consequence(s)"))%>%
  mutate(Percentage = (n.x/n.y)*100)

###PLotiing variant and hpo term for each gene

variant_hpoterm_df <- data.frame(variant_hpoterm)

variant_hpoterm_df <- variant_hpoterm_df %>%
  filter(n.y >3 & n.x >1)%>%
  select(1,2,3,6)



unique_genes <- unique(variant_hpoterm_df$`Diagnosed.Gene`)

for (gene in unique_genes) {
  
  gene_data <- subset(variant_hpoterm_df, `Diagnosed.Gene` == gene)
  
  
  plot <- ggplot(gene_data, aes(y = `hpo_ancestors_description`,x=Percentage, fill = Consequence.s.)) +
    geom_bar(stat="identity",position = "dodge") +
    labs(title = paste("Phenotypes Occurrence for", gene), x = "Percentage", y = "HPO ANcestor term") +
    theme_minimal()
  
  print(plot)
  
}





# PanelApp Analysis -------------------------------------------------------




panelapp_genes_idgenes <- panelapp_genes %>%
  filter(gene_symbol %in% id_genes$`Diagnosed Gene`) %>%
  group_by(gene_symbol) %>%
  tally()

ggplot(panelapp_genes_idgenes, aes(y = reorder( gene_symbol, n) , x = n, fill = gene_symbol)) +
  geom_bar(stat = "identity") +
  labs(title = "Gene Diagnosis Counts", x = "Gene", y = "Panel numbers") +
  theme_minimal()



panelapp_genes_idgenes_panels <- panelapp_genes %>%
  filter(gene_symbol %in% id_genes$`Diagnosed Gene`) %>%
  group_by(panel_name,gene_symbol)


# Statistical Analysis ----------------------------------------------------



#Statistical analysis of Diagnosed Gene with HPO terms####


####Chi square test for Diganosed Gene and HPO term contingnecy table as a whole


genopheno_all_separated <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  select(1,2,8)


contingency_table <- genopheno_all_separated%>%
  count(`Diagnosed Gene`, `HPO terms(s)`) %>%
  pivot_wider(names_from = `HPO terms(s)`, values_from = n, values_fill = list(n = 0))

contingency_matrix <- as.matrix(contingency_table[,-1])
rownames(contingency_matrix) <- contingency_table$`Diagnosed Gene`

chi_test <- chisq.test(contingency_matrix)

p_value <- chi_test$p.value

print(paste("P-Value:", p_value))

gene_hpo_combinations <- expand.grid(
  `Diagnosed Gene` = rownames(chi_test$observed),
  `HPO terms(s)` = colnames(chi_test$observed)
)

results_table <- data.frame(
  gene_hpo_combinations,
  Observed = as.vector(chi_test$observed),
  Expected = as.vector(chi_test$expected),
  Residuals = as.vector(chi_test$residuals),
  P_Value = p_value
)


results_table$StandardizedResiduals <- with(results_table, (Observed - Expected) / sqrt(Expected))

results_table$Significant <- abs(results_table$StandardizedResiduals) > 1.96

print(results_table)

results_table_chisq <- results_table %>%
  filter(Observed >= Expected)%>%
  filter(Significant=="TRUE")


# Load necessary libraries
library(ggplot2)
library(reshape2)

unique_genes <- unique(results_table_chisq$`Diagnosed.Gene`)

for (gene in unique_genes) {
  gene_data <- subset(results_table_chisq, `Diagnosed.Gene` == gene)
  
  heatmap_data <- dcast(gene_data, `Diagnosed.Gene` ~ `HPO.terms.s.`, value.var = "Residuals")
  
  heatmap_melted <- melt(heatmap_data, id.vars = "Diagnosed.Gene")
  
  p <- ggplot(heatmap_melted, aes(x = variable, y = `Diagnosed.Gene`, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Residuals") +
    labs(title = paste("Heatmap of Residuals for", gene),
         x = "HPO Terms(s)",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(p)
}

results_table_chisq <- results_table %>%
  filter(Observed >= Expected)%>%
  filter(Significant=="TRUE")%>%
  filter(Observed >1)





#Logistic regression for GENE HPO term
contingency_table_lr <- genopheno_all_separated%>%
  filter(`HPO terms(s)` %in% most_occuring_hpoterms$`HPO terms(s)`)%>%
     count(`Diagnosed Gene`, `HPO terms(s)`) %>%
     pivot_wider(names_from = `HPO terms(s)`, values_from = n, values_fill = list(n = 0))


long_data_lr <- contingency_table_lr %>%
     pivot_longer(
         cols = -`Diagnosed Gene`,                    
         names_to = "HPO_Term",           
         values_to = "Count"              
       )


long_data_lr <- long_data_lr %>%
  mutate(
    Presence = ifelse(Count > 0, 1, 0)  
  )

long_data_lr$`Diagnosed Gene` <- as.factor(long_data_lr$`Diagnosed Gene`)
long_data_lr$`HPO_Term` <- as.factor(long_data_lr$`HPO_Term`)

model <- glm(Presence ~ `Diagnosed Gene` + HPO_Term, data = long_data_lr, family = binomial)
summary(model)

summary_model <- summary(model)
print(summary_model$coefficients)
p_values <- summary_model$coefficients[, "Pr(>|z|)"]


results_table_lr <- data.frame(
  Term = rownames(summary_model$coefficients),
  Estimate = summary_model$coefficients[, "Estimate"],
  Std_Error = summary_model$coefficients[, "Std. Error"],
  z_value = summary_model$coefficients[, "z value"],
  P_Value = p_values
)


long_data_lr$Predicted_Probability <- predict(model, type = "response")


predictions <- predict(model, type = "response")
long_data_lr$Predicted_Probability <- predictions

print(long_data_lr)


library(pROC)

roc_curve <- roc(long_data_lr$Presence, long_data_lr$Predicted_Probability)

plot(roc_curve, main = "ROC Curve HPO term/Gene", col = "blue", lwd = 2)
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "blue", lwd = 2)


heatmap_data_lr <- long_data_lr %>%
  select(`Diagnosed Gene`, HPO_Term, Predicted_Probability) %>%
  spread(key = HPO_Term, value = Predicted_Probability)

# Plot heatmap
ggplot(melt(heatmap_data_lr), aes(y = variable, x = `Diagnosed Gene`, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(
    title = "Heatmap of Predicted Probabilities",
    x = "HPO Term",
    y = "Gene"
  )



#Fisher test for each gene/hpo term
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





##########Statistical analysis of Diagnosed Gene with HPO ancestor####


#Chi square test for Diganosed Gene and HPO ancestor contingnecy table as a whole

genopheno_all_separated2 <- genopheno %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`) %>%
  separate_rows(`HPO terms(s)`, sep = "\\|")

joined_hpo <- genopheno_all_separated2 %>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")

joined_hpo <- joined_hpo %>%
  select(1,2,17)%>%
  unique




contingency_table_ancestor <- joined_hpo%>%
  count(`Diagnosed Gene`, `hpo_ancestors_description`) %>%
  pivot_wider(names_from = `hpo_ancestors_description`, values_from = n, values_fill = list(n = 0))

contingency_matrix_ancestor <- as.matrix(contingency_table_ancestor[,-1])
rownames(contingency_matrix_ancestor) <- contingency_table_ancestor$`Diagnosed Gene`

chi_test <- chisq.test(contingency_matrix_ancestor)

p_value <- chi_test$p.value

print(paste("P-Value:", p_value))

gene_hpo_ancestor_combinations <- expand.grid(
  `Diagnosed Gene` = rownames(chi_test$observed),
  `hpo_ancestors_description` = colnames(chi_test$observed)
)

results_table_ancestor <- data.frame(
  gene_hpo_ancestor_combinations,
  Observed = as.vector(chi_test$observed),
  Expected = as.vector(chi_test$expected),
  Residuals = as.vector(chi_test$residuals),
  P_Value = p_value
)


results_table_ancestor$StandardizedResiduals <- with(results_table_ancestor, (Observed - Expected) / sqrt(Expected))

results_table_ancestor$Significant <- abs(results_table_ancestor$StandardizedResiduals) > 1.96

print(results_table_ancestor)





unique_genes <- unique(results_table_ancestor$`Diagnosed.Gene`)

for (gene in unique_genes) {
  gene_data <- subset(results_table_ancestor, `Diagnosed.Gene` == gene)
  
  heatmap_data <- dcast(gene_data, `Diagnosed.Gene` ~ `hpo_ancestors_description`, value.var = "Residuals")
  
  heatmap_melted <- melt(heatmap_data, id.vars = "Diagnosed.Gene")
  
  p <- ggplot(heatmap_melted, aes(x = variable, y = `Diagnosed.Gene`, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid  = "white", high = "red", midpoint = 0, name = "Residuals") +
    labs(title = paste("Heatmap of Residuals for", gene),
         x = "HPO Ancestor Terms",
         y = "Diagnosed Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

results_table_ancestor <- results_table_ancestor %>%
  filter(Significant=="TRUE")





##LOgistic regression HPO ancestor Gene

contingency_table_lr_ancestor <- joined_hpo%>%
  mutate(`HPO terms(s)`= hpo_ancestors_description)%>%
  select(1,2,4)%>%
  count(`Diagnosed Gene`, `HPO terms(s)`) %>%
  pivot_wider(names_from = `HPO terms(s)`, values_from = n, values_fill = list(n = 0))


long_data_lr_ancestor <- contingency_table_lr_ancestor %>%
  pivot_longer(
    cols = -`Diagnosed Gene`,                    
    names_to = "HPO_Term",           
    values_to = "Count"              
  )


long_data_lr_ancestor <- long_data_lr_ancestor %>%
  mutate(
    Presence = ifelse(Count > 0, 1, 0)  
  )

long_data_lr_ancestor$`Diagnosed Gene` <- as.factor(long_data_lr_ancestor$`Diagnosed Gene`)
long_data_lr_ancestor$`HPO_Term` <- as.factor(long_data_lr_ancestor$`HPO_Term`)

model2 <- glm(Presence ~ `Diagnosed Gene` + HPO_Term, data = long_data_lr_ancestor, family = binomial)
summary(model2)


coefficients <- summary(model2)$coefficients
p_values <- coefficients[, 4] 
coefficients_df <- as.data.frame(coefficients)

long_data_lr_ancestor$Predicted_Probability <- predict(model2, type = "response")


predictions <- predict(model2, type = "response")
long_data_lr_ancestor$Predicted_Probability <- predictions

print(long_data_lr_ancestor)


library(pROC)

roc_curve <- roc(long_data_lr_ancestor$Presence, long_data_lr_ancestor$Predicted_Probability)

plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), col = "blue", lwd = 2)


heatmap_data_lr_ancestor <- long_data_lr_ancestor %>%
  select(`Diagnosed Gene`, HPO_Term, Predicted_Probability) %>%
  spread(key = HPO_Term, value = Predicted_Probability)

# Plot heatmap
ggplot(melt(heatmap_data_lr_ancestor), aes(y = variable, x = `Diagnosed Gene`, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue",high = "red3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(
    title = "Heatmap of Predicted Probabilities",
    x = "HPO Ancestor Term",
    y = "Gene"
  )

#Fisher test for each gene/hpo ancestor
#########################################################################

hpo_toplevel_df <- joined_hpo %>%
  group_by(`Diagnosed Gene`, hpo_ancestors_description) %>%
  tally() %>%
  arrange(-n)

######Fisher test with OR##

genes <- unique(hpo_toplevel_df$`Diagnosed Gene`)
hpo_terms <- unique(hpo_toplevel_df$`hpo_ancestors_description`)



results_ancestor_fisher <- data.frame(Gene = character(), 
                                      HPO_term = character(), 
                                      p_value = numeric(), 
                                      odds_ratio = numeric(),
                                      stringsAsFactors = FALSE)

# Loop through each gene and HPO term to perform Fisher's test and calculate Odds Ratio
for (gene in genes) {
  for (hpo in hpo_terms) {
    
    ###gene = "ANKRD11"
    ###hpo = hpo_terms[1]
    
    A <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` == gene & hpo_toplevel_df$hpo_ancestors_description == hpo])
    B <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` == gene & hpo_toplevel_df$hpo_ancestors_description != hpo])
    C <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` != gene & hpo_toplevel_df$hpo_ancestors_description == hpo])
    D <- sum(hpo_toplevel_df$n[hpo_toplevel_df$`Diagnosed Gene` != gene & hpo_toplevel_df$hpo_ancestors_description != hpo])
    
    contingency_table <- matrix(c(A, B, C, D), nrow = 2)
    
    
    fisher_test <- fisher.test(contingency_table)
    
    
    odds_ratio <- ifelse(!is.null(fisher_test$estimate), fisher_test$estimate, NA)
    
    
    results_ancestor_fisher <- rbind(results_ancestor_fisher, data.frame(`Diagnosed Gene` = gene, 
                                                                         `hpo_ancestors_description` = hpo, 
                                                                         p_value = fisher_test$p.value, 
                                                                         odds_ratio = odds_ratio))
  }
}

# Display the results
print(results_ancestor_fisher)

results_ancestor_fisher$p_adjusted <- p.adjust(results_ancestor_fisher$p_value, method = "BH")



##Statistical test for ID/OD ####

compare_df_id_od <- genopheno_id_separated %>%
  left_join(genopheno_no_id_separated, by=c("Diagnosed Gene", "HPO terms(s)")) %>%
  mutate(n_id=n.x, n_od=n.y)%>%
  select(1,2,5,6) %>%
  replace(is.na(.), 0)
  
results_id_od <- data.frame(
  Diagnosed_Gene = character(),
  HPO_Term = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)


# Loop through each row in the data frame
for (i in 1:nrow(compare_df_id_od)) {
  # Extract the current row
  row <- compare_df_id_od[i, ]
  
  # Create a contingency table
  contingency_table <- matrix(
    c(row$n_id, row$n_od,
      sum(compare_df_id_od$n_id) - row$n_id,
      sum(compare_df_id_od$n_od) - row$n_od),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      Category = c("Intellectual_Disability", "Other_Disease"),
      HPO_Term_Presence = c("Present", "Absent")
    )
  )
  
  # Perform the Chi-Square test (use Fisher's exact test if needed)
  if (any(contingency_table < 5)) {
    test_result <- fisher.test(contingency_table)
  } else {
    test_result <- chisq.test(contingency_table)
  }
  
  # Append the results to the results data frame
  results_id_od <- rbind(
    results_id_od,
    data.frame(
      Diagnosed_Gene = row$`Diagnosed Gene`,
      HPO_Term = row$`HPO terms(s)`,
      P_Value = test_result$p.value,
      stringsAsFactors = FALSE
    )
  )
}


results_id_od$Adjusted_P_Value <- p.adjust(results_id_od$P_Value, method = "BH")


# Print the results
print(results_id_od)



##Statistical analysis for variant type and hpo term####


variant_consequence <- genopheno %>%
  separate_rows(`HPO terms(s)`, sep = "\\|") %>%
  filter(`Diagnosed Gene` %in% id_genes$`Diagnosed Gene`)%>%
  left_join(hpo1, by = c("HPO terms(s)" = "hpo_term_description")) %>%
  filter(top_level_phenotypic_abnormality != "-")%>%
  unique() %>%
  group_by(`Diagnosed Gene`, `Consequence(s)`, `hpo_ancestors_description`) %>%
  tally()



results_df_variant <- data.frame(
  Diagnosed_Gene = character(),
  Chi_Square_Statistic = numeric(),
  Degrees_of_Freedom = numeric(),
  P_Value = numeric(),
  Significant = logical(),
  stringsAsFactors = FALSE
)


genes <- unique(variant_consequence$`Diagnosed Gene`)

for (gene in genes) {
  
  gene_data <- variant_consequence %>% filter(`Diagnosed Gene` == gene)
  
  
  contingency_table <- xtabs(n ~ `Consequence(s)` + hpo_ancestors_description, data = gene_data)
  
  
  chi_test <- chisq.test(contingency_table)
  
  
  
  is_significant <- chi_test$p.value < 0.05
  
  
  results_df_variant <- rbind(
    results_df_variant,
    data.frame(
      Diagnosed_Gene = gene,
      Chi_Square_Statistic = chi_test$statistic,
      Degrees_of_Freedom = chi_test$parameter,
      P_Value = chi_test$p.value,
      Significant = is_significant
    )
  )
}


results_df_variant$Adjusted_P_Value <- p.adjust(results_df_variant$P_Value, method = "BH")





results_df_variant$Significant_Adjusted <- results_df_variant$Adjusted_P_Value < 0.05

print(results_df_variant)



##
#######Statistical analysis moi


moi_hpoterm <- genopheno_separated %>%
  group_by(`Diagnosed Gene`, MOI, `HPO terms(s)`) %>%
  tally() %>%
  filter(`Diagnosed Gene` %in% c("KDM5C", "MECP2" ,"SCN2A"))

data <-data.frame(moi_hpoterm)



# KDM5C WORK --------------------------------------------------------------



pivoted_data_moi_kdm5c<- data %>%
  group_by(Diagnosed.Gene, HPO.terms.s.) %>%
  summarise(
    XD = sum(n[MOI == "_XD"], na.rm = TRUE),
    XR = sum(n[MOI == "_XR"], na.rm = TRUE),
  )  %>%
  filter(!(XD == 0 & XR == 0))%>%
  filter(`Diagnosed.Gene` == "KDM5C") %>%
  select(2,3,4)



#Clustering Analysis for KDM5C

agg_data <- pivoted_data_moi_kdm5c %>%
  group_by(HPO.terms.s.) %>%
  summarise(XD = sum(XD), XR = sum(XR))


# Step 2: Normalize the data
scaled_data <- agg_data %>%
  mutate(XD_scaled = scale(XD), XR_scaled = scale(XR)) %>%
  select(XD_scaled, XR_scaled)


set.seed(123)  
kmeans_result <- kmeans(scaled_data, centers = 3)  


agg_data$Cluster <- as.factor(kmeans_result$cluster)


ggplot(agg_data, aes(x = XD, y = XR, color = Cluster)) +
  geom_point(size = 0.5) +
  geom_text(aes(label = HPO.terms.s., vjust = -0.5, hjust = 0.5)) +
  labs(title = "KDM5C Clustering of HPO Terms Based on MOI", 
       x = "XD Counts", y = "XR Counts") +
  theme_minimal()


print(agg_data)

##

# MECP2 WORK --------------------------------------------------------------


pivoted_data_moi_mecp2<- data %>%
  group_by(Diagnosed.Gene, HPO.terms.s.) %>%
  filter(`Diagnosed.Gene`== "MECP2") %>%
  summarise(
    XD = sum(n[MOI == "_XD"], na.rm = TRUE),
    XR = sum(n[MOI == "_XR"], na.rm = TRUE),
  )  %>%
  filter(!(XD == 0 & XR == 0))



#Clustering Analysis for MECP2

agg_data <- pivoted_data_moi_mecp2 %>%
  group_by(HPO.terms.s.) %>%
  summarise(XD = sum(XD), XR = sum(XR))


# Step 2: Normalize the data
scaled_data <- agg_data %>%
  mutate(XD_scaled = scale(XD), XR_scaled = scale(XR)) %>%
  select(XD_scaled, XR_scaled)


set.seed(123)  
kmeans_result <- kmeans(scaled_data, centers = 3)  


agg_data$Cluster <- as.factor(kmeans_result$cluster)


ggplot(agg_data, aes(x = XD, y = XR, color = Cluster)) +
  geom_point(size = 0.5) +
  geom_text(aes(label = HPO.terms.s., vjust = -0.5, hjust = 0.5)) +
  labs(title = "MECP2 Clustering of HPO Terms Based on MOI", 
       x = "XD Counts", y = "XR Counts") +
  theme_minimal()


print(agg_data)




# SN2A WORK ---------------------------------------------------------------



pivoted_data_moi_scn2a<- data %>%
  group_by(Diagnosed.Gene, HPO.terms.s.) %>%
  filter(`Diagnosed.Gene`== "SCN2A") %>%
  summarise(
    AD = sum(n[MOI == "_AD"], na.rm = TRUE),
    AR = sum(n[MOI == "_AR"], na.rm = TRUE),
  )  %>%
  filter(!(AD == 0 & AR == 0))



agg_data <- pivoted_data_moi_scn2a %>%
  group_by(HPO.terms.s.) %>%
  summarise(AD = sum(AD), AR = sum(AR))


# Step 2: Normalize the data
scaled_data <- agg_data %>%
  mutate(AD_scaled = scale(AD), AR_scaled = scale(AR)) %>%
  select(AD_scaled, AR_scaled)


set.seed(123)  
kmeans_result <- kmeans(scaled_data, centers = 3)  


agg_data$Cluster <- as.factor(kmeans_result$cluster)


ggplot(agg_data, aes(x = AD, y = AR, color = Cluster)) +
  geom_point(size = 4) +
  geom_text(aes(label = HPO.terms.s., vjust = -0.5, hjust = 0.5)) +
  facet_wrap(~ Cluster)
labs(title = "Clustering of HPO Terms Based on AD and AR", 
     x = "AD Counts", y = "AR Counts") +
  theme_minimal()


print(agg_data)
