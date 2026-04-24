#### Lesion removed PCA 

## Remove lesion birds from dataset 
all_metrics_cnt_no <- all_metrics_inj %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes")
all_metrics_cnt_no <- all_metrics_cnt_no %>% filter(BirdID != "bl191or41" & BirdID != "bl10or50" & BirdID != "bl25or15")

###PCA on sleap metrics - whole exp 
  #Going with simple version - 4 metrics only 
  # Removing Lesion birds 
metrics_inj_simp <- all_metrics_cnt_no %>%
  select(
    mean_speed,
    head_var,
    concave_area,
    shannon_entropy,
  ) %>%
  drop_na()

pca_res_inj_simp <- prcomp(
  metrics_inj_simp,
  center = TRUE,
  scale. = TRUE
)

summary(pca_res_inj_simp)

pca_scores_inj_simp <- as.data.frame(pca_res_inj_simp$x) %>%
  bind_cols(
    all_metrics_cnt_no %>%
      select(BirdID, Condition, time_bin_3min, Inj_Lesion) %>%
      slice(rownames(metrics_inj_simp) |> as.integer())
  )

loadings_inj_simp <- as.data.frame(pca_res_inj_simp$rotation)
loadings_inj_simp
screeplot(pca_res_inj_simp, type = "lines", main = "Scree Plot (Base R)")

variance <- pca_res_inj_simp$sdev^2
pve <- variance / sum(variance)

scree_df <- data.frame(
  PC = 1:length(pve),
  variance = pve
)

ggplot(scree_df, aes(x = PC, y = variance)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    alpha = 0.8,
    fill = "steelblue"
  ) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 1:nrow(scree_df)) +
  labs(
    x = "Principal Component",
    y = "Proportion of Variance Explained",
    title = ""
  ) +
  theme_minimal()+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


# Install and load ggfortify if not already installed
install.packages("ggfortify")
library(ggfortify)

autoplot(pca_res_inj_simp, data = all_metrics_cnt_no,
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  geom_point(size = 1, alpha = 0.3)


#    ----------------------------------------------------------------------





pca_scores_inj_simp$time_bin_3min <- factor(pca_scores_inj_simp$time_bin_3min, ordered = TRUE)

Inj_PC1_model_simp <- glmmTMB(PC1 ~ Condition * time_bin_3min * Inj_Lesion + (1|BirdID),
                         family = gaussian(link = "identity"), data = pca_scores_inj_simp)
Anova(Inj_PC1_model_simp)

summary(Inj_PC1_model_simp)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_PC1_model_simp, plot = TRUE)

##Just IgG
pca_igg_simp <- pca_scores_inj_simp %>% filter(Inj_Lesion == "IgG-SAP_No")

pca_igg_simp$time_bin_3min <- factor(pca_igg_simp$time_bin_3min, ordered = TRUE)

IgG_PC1_model_simp <- glmmTMB(PC1 ~ Condition * time_bin_3min + (1 | BirdID),
                         family = gaussian(link = "identity"), data = pca_igg_simp)

Anova(IgG_PC1_model_simp)

summary(IgG_PC1_model_simp)

pairs(emmeans(IgG_PC1_model_simp, ~ Condition | time_bin_3min), adjust = "fdr")

#residuals
simulationOutput <- simulateResiduals(fittedModel = IgG_PC1_model, plot = TRUE)

##just anti-dbh_no 
pca_anti_no_simp <- pca_scores_inj_simp %>% filter(Inj_Lesion == "anti-DBH-SAP_No")
pca_anti_no_simp$time_bin_3min <- factor(pca_anti_no_simp$time_bin_3min, ordered = TRUE)

No_PC1_model_simp <- glmmTMB(PC1 ~ Condition * time_bin_3min + (1|BirdID),
                        family = gaussian(link = "identity"), data = pca_anti_no_simp)

Anova(No_PC1_model_simp)

summary(No_PC1_model_simp)

pairs(emmeans(No_PC1_model_simp, ~ Condition | time_bin_3min), adjust = "fdr")

#residuals
simulationOutput <- simulateResiduals(fittedModel = No_PC1_model_simp, plot = TRUE)

### Manip Raw points and model predictions plot ######### 
emm_PCA_inj_simp <- emmeans(
  Inj_PC1_model_simp,
  ~ Condition * time_bin_3min | Inj_Lesion
)


#to plot raw values in back (no log on data)
emm_df_PCA_inj_simp <- as.data.frame(emm_PCA_inj_simp) %>%
  mutate(
    area_hat = emmean,
    lower.CL = lower.CL,
    upper.CL = upper.CL
  )

ggplot(pca_scores_inj_simp,
       aes(x = time_bin_3min, y = PC1, color = Condition)) +
  
  # Raw data
  geom_jitter(
    aes(group = BirdID),
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  geom_line(
    aes(group = BirdID),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  
  # Model-predicted means (correct scale)
  geom_line(
    data = emm_df_PCA_inj_simp,
    aes(y = area_hat, group = Condition),
    linewidth = 1.2
  ) +
  
  # 95% CI ribbon
  geom_ribbon(
    data = emm_df_PCA_inj_simp,
    aes(
      y = area_hat,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = Condition,
      group = Condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  
  facet_wrap(~Inj_Lesion, nrow = 1, labeller = labeller(Inj_Lesion = c("IgG-SAP_No" = "IgG-SAP", 
                                                                       "anti-DBH-SAP_No" = "DBH-SAP")))+
  scale_color_manual(values = c("F"="purple", "N"="orange")) +
  scale_fill_manual(values = c("F"="purple", "N"="orange")) +
  labs(x = "Time bin (3 min)", y = "Movement Intensity (PC1)") +
  theme_minimal()+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )




##########.  FIRST 3 MINs ############### 
all_metrics_30sec_simp <- all_metrics_inj_30sec %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes")
all_metrics_30sec_simp <- all_metrics_30sec_simp %>% filter(BirdID != "bl191or41" & BirdID != "bl10or50" & BirdID != "bl25or15")

all_metrics_inj_simp_bin1_30 <- all_metrics_30sec_simp %>% filter(time_bin_30sec < 7)

###PCA on all sleap metrics 
metrics_30sec_simp <- all_metrics_inj_simp_bin1_30 %>%
  select(
    mean_speed,
    head_var,
    concave_area,
    shannon_entropy
  ) %>%
  drop_na()

pca_res_30sec_simp <- prcomp(
  metrics_30sec_simp,
  center = TRUE,
  scale. = TRUE
)

summary(pca_res_30sec_simp)

pca_scores_30sec_simp <- as.data.frame(pca_res_30sec_simp$x) %>%
  bind_cols(
    all_metrics_inj_simp_bin1_30 %>%
      select(BirdID, Condition, time_bin_30sec, Inj_Lesion) %>%
      slice(rownames(metrics_30sec_simp) |> as.integer())
  )

loadings_30sec_simp <- as.data.frame(pca_res_30sec_simp$rotation)
loadings_30sec_simp

screeplot(pca_res_30sec_simp, type = "lines", main = "Scree Plot (Base R)")



# first 3 min PCA tables --------------------------------------------------

# Extract PCA summary
pca_summary_3min <- summary(pca_res_30sec_simp)

pca_table <- data.frame(
  PC = paste0("PC", seq_along(pca_summary_3min$importance[1, ])),
  SD = round(pca_summary_3min$importance["Standard deviation", ], 3),
  Proportion_Variance = round(pca_summary_3min$importance["Proportion of Variance", ], 3),
  Cumulative_Variance = round(pca_summary_3min$importance["Cumulative Proportion", ], 3)
)

library(gt)

pca_table %>%
  gt() %>%
  tab_header(
    title = "Principal Component Analysis Summary (First 3 Min Data)",
    subtitle = "Standard deviation and variance explained by each component"
  ) %>%
  cols_label(
    PC = "Component",
    SD = "Std. Dev.",
    Proportion_Variance = "Proportion of Variance",
    Cumulative_Variance = "Cumulative Variance"
  ) %>%
  fmt_number(
    columns = c(SD, Proportion_Variance, Cumulative_Variance),
    decimals = 3
  )

# Extract loadings matrix
loadings_df <- as.data.frame(pca_res_30sec_simp$rotation)

# Add variable names as a column
loadings_df <- data.frame(
  Variable = rownames(loadings_df),
  loadings_df,
  row.names = NULL
)

# Round all PC columns
loadings_df[, -1] <- round(loadings_df[, -1], 3)

library(gt)
loadings_df %>%
  gt() %>%
  tab_header(
    title = "PCA Loadings (First 3 min)",
    subtitle = "Contribution of each variable to the principal components"
  ) %>%
  cols_label(
    Variable = "Variable"
  ) %>%
  fmt_number(
    columns = -1,  # all columns except Variable
    decimals = 3
  ) %>%
  # Optional: color-code by loading magnitude
  data_color(
    columns = -1,
    fn = scales::col_numeric(
      palette = c("tomato", "white", "lightgreen"),
      domain = c(-1, 1)
    )
  )






#    ----------------------------------------------------------------------



pca_scores_30sec_simp$time_bin_30sec <- factor(pca_scores_30sec_simp$time_bin_30sec, ordered = TRUE)
Inj_PC1_model_bin1_30_simp <- glmmTMB(PC1 ~ Condition * Inj_Lesion * time_bin_30sec + (1|BirdID),      
                                 family = gaussian(link = "identity"), data = pca_scores_30sec_simp)

Anova(Inj_PC1_model_bin1_30_simp)

summary(Inj_PC1_model_bin1_30_simp)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_PC1_model_bin1_30, plot = TRUE)


### UM Raw points and model predictions plot ######### 
emm_PCA_inj_30_simp <- emmeans(
  Inj_PC1_model_bin1_30_simp,
  ~ Condition * time_bin_30sec | Inj_Lesion
)


#to plot raw values in back (no log on data)
emm_df_PCA_inj_30_simp <- as.data.frame(emm_PCA_inj_30_simp) %>%
  mutate(
    area_hat = emmean,
    lower.CL = lower.CL,
    upper.CL = upper.CL
  )

pca_scores_30sec_simp <- pca_scores_30sec_simp %>%
  mutate(Inj_Lesion = factor(Inj_Lesion, levels = c("IgG-SAP_No", "anti-DBH-SAP_No")))
emm_df_PCA_inj_30_simp <- emm_df_PCA_inj_30_simp %>%
  mutate(Inj_Lesion = factor(Inj_Lesion, levels = c("IgG-SAP_No", "anti-DBH-SAP_No")))

ggplot(pca_scores_30sec_simp,
       aes(x = time_bin_30sec, y = PC1, color = Condition)) +
  
  # Raw data
  geom_jitter(
    aes(group = BirdID),
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  
  geom_line(
    aes(group = BirdID),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  
  # Model-predicted means (correct scale)
  geom_line(
    data = emm_df_PCA_inj_30_simp,
    aes(y = area_hat, group = Condition),
    linewidth = 1.2
  ) +
  
  # 95% CI ribbon
  geom_ribbon(
    data = emm_df_PCA_inj_30_simp,
    aes(
      y = area_hat,
      ymin = lower.CL,
      ymax = upper.CL,
      fill = Condition,
      group = Condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  
  facet_wrap(~Inj_Lesion, nrow = 1, labeller = labeller(Inj_Lesion = c("IgG-SAP_No" = "IgG", 
                                                                       "anti-DBH-SAP_No" = "DBH")))+
  scale_color_manual(values = c("F"="purple", "N"="orange")) +
  scale_fill_manual(values = c("F"="purple", "N"="orange")) +
  labs(x = "Time bin (30 sec)", y = "Movement Intensity (PC1)") +
  theme_minimal()+
  #scale_x_discrete(
  # labels = function(x) {
  #  x <- as.numeric(x)
  # paste0(((x - 1) * 30), "-", x * 30)
  #  }
  #)+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )




##Just IgG
pca_igg_bin1_30_simp <- pca_scores_30sec_simp %>% filter(Inj_Lesion == "IgG-SAP_No")

pca_igg_bin1_30_simp$time_bin_30sec <- factor(pca_igg_bin1_30_simp$time_bin_30sec, ordered = TRUE)

IgG_PC1_model_bin1_30_simp <- glmmTMB(PC1 ~ Condition * time_bin_30sec + (1 | BirdID),
                                 family = gaussian(link = "identity"), data = pca_igg_bin1_30_simp)

Anova(IgG_PC1_model_bin1_30_simp)

summary(IgG_PC1_model_bin1_30_simp)

pairs(emmeans(IgG_PC1_model_bin1_30_simp, ~ Condition | time_bin_30sec), adjust = "fdr")

#residuals
simulationOutput <- simulateResiduals(fittedModel = IgG_PC1_model_bin1_30, plot = TRUE)



##just anti-dbh_no 
pca_dbh_no_bin1_30_simp <- pca_scores_30sec_simp %>% filter(Inj_Lesion == "anti-DBH-SAP_No")

pca_dbh_no_bin1_30_simp$time_bin_30sec <- factor(pca_dbh_no_bin1_30_simp$time_bin_30sec, ordered = TRUE)

No_PC1_model_bin1_30_simp <- glmmTMB(PC1 ~ Condition * time_bin_30sec + (1 | BirdID),
                                family = gaussian(link = "identity"), data = pca_dbh_no_bin1_30_simp)

Anova(No_PC1_model_bin1_30_simp)

summary(No_PC1_model_bin1_30_simp)

pairs(emmeans(No_PC1_model_bin1_30_simp, ~ Condition | time_bin_30sec), adjust = "fdr")

#residuals
simulationOutput <- simulateResiduals(fittedModel = No_PC1_model_bin1_30_simp, plot = TRUE)



