
##Clustering analysis code for Gene HPO term visualisation









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
