
##Clustering analysis for Gene HPO term visualisation









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
