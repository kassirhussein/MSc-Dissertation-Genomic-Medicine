


#Logistic regression model for GENE HPO term correlations



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

