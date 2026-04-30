# Principal Component Analysis of Pose-estimation metrics # 
library(ggfortify)

# PCA -----------------------------------------------------------
#load in all metrics df, removing lesion birds 
all_metrics_cnt_no <- all_metrics_inj %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes")

###PCA on sleap metrics - whole exp 
# Going with simple version - 4 metrics only 
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

variance <- pca_res_inj_simp$sdev^2
pve <- variance / sum(variance)

scree_df <- data.frame(
  PC = 1:length(pve),
  variance = pve
)

#Scree plot (Fig 6a)
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

#Vector bi-plot (fig 6b)
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


# Stats PCA ---------------------------------------------------------------
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

### Raw points and model predictions plot ######### 
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





# First 3 mins only PCA  --------------------------------------------------
  #Re-bin data and calculate metrics over 30-sec bins 
#Make 30 sec time bins 
Inj_Birds_all_30sec <- Inj_Birds_all %>%
  group_by(BirdID) %>%
  mutate(
    time_bin_30sec = cut(adj_frame, breaks = 30, labels = FALSE)  # Divide into 30 equal bins (30 sec)
  ) %>%
  ungroup()

#CV in 30 sec bins 
### Inj circ var 
circular_variance_inj_30sec <- Inj_Birds_all_30sec %>%
  group_by(BirdID, Condition, Inj_Lesion, time_bin_30sec) %>%
  summarise(head_var = var.circular(circular(head_angle), na.rm=TRUE), .groups = "drop")

#Area 30 sec bins 
Concave_hull_areas_3min_inj_30 <- Inj_Birds_all_30sec %>%
  group_by(BirdID, Condition, time_bin_30sec, Inj_Lesion) %>%
  summarise(area = concave_hull_area(across(c(x_head, y_head_neg))))  # Correct way to pass data

#Entropy 30 sec bins 
##Inj 
compute_spatial_entropy_time_inj <- function(data, grid_width, grid_height, grid_size = 100, jitter_amount = 0.0001) {
  results <- data %>%
    group_by(BirdID, Condition, time_bin_30sec, Inj_Lesion) %>%
    summarise(
      entropy = list({
        # Subset data for the current group
        data_group <- cur_data()
        
        # Jitter positions to avoid duplicate points
        data_group$x_jittered <- jitter(data_group$x_head, amount = jitter_amount)
        data_group$y_jittered <- jitter(data_group$y_head, amount = jitter_amount)
        
        #Ensure jitter is not above grid size
        data_group$x_jittered <- pmin(
          pmax(data_group$x_jittered, 0),
          grid_width
        )
        
        data_group$y_jittered <- pmin(
          pmax(data_group$y_jittered, 0),
          grid_height
        )
        
        
        # Define fixed grid (based on total data)
        x_breaks <- seq(0, grid_width, length.out = grid_size + 1)
        y_breaks <- seq(0, grid_height, length.out = grid_size + 1)
        
        # Assign each (x, y) to a grid cell
        x_bins <- cut(data_group$x_jittered, breaks = x_breaks,
                      labels = FALSE, include.lowest = TRUE)
        y_bins <- cut(data_group$y_jittered, breaks = y_breaks,
                      labels = FALSE, include.lowest = TRUE)
        grid_labels <- factor(paste(x_bins, y_bins, sep = "_"))  # Grid cell ID
        
        # Create full grid with all possible grid labels
        full_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size) %>%
          mutate(grid_id = factor(paste(x, y, sep = "_")))
        
        # Count visits to each grid cell
        grid_counts <- as.data.frame(table(grid_labels))
        colnames(grid_counts) <- c("grid_id", "count")
        
        # Merge with full grid to include unvisited cells
        grid_counts <- full_grid %>%
          left_join(grid_counts, by = "grid_id") %>%
          mutate(count = ifelse(is.na(count), 0, count))  # Set unvisited cells to zero
        
        # Compute probabilities
        grid_counts$prob <- grid_counts$count / sum(grid_counts$count)
        
        # Define spatial window
        win <- owin(xrange = c(0, grid_width), yrange = c(0, grid_height))
        
        # Create spatial point pattern with jittered points
        ppp_obj <- ppp(data_group$x_jittered, data_group$y_jittered, window = win, marks = grid_labels)
        
        # Compute Shannon entropy
        entropy_values <- shannon(ppp_obj)
        
        # Return entropy values
        tibble(
          shannon_entropy = entropy_values$shann,
          rel_shannon_entropy = entropy_values$rel.shann, 
          occupied_grid_cells = length(unique(grid_labels))  # This is I
        )
      }),
      .groups = "drop"
    ) %>%
    unnest(entropy)  # Unpack entropy columns
  
  return(results)
}

Inj_int_ent_30sec <- compute_spatial_entropy_time_inj(Inj_Birds_all_30sec, grid_width = 40, grid_height = 21, grid_size = grid_size)

#Average speed 30 sec bins 
Inj_binned_speed_30sec <- Inj_Birds_all_30sec %>%
  group_by(BirdID, Condition, time_bin_30sec, Inj_Lesion) %>%
  summarise(
    median_head_speed = median(speed_head_psec),
    mean_speed = mean(speed_head_psec),
    .groups = "drop"
  )

#Combine metric dataframes 
Inj_int_ent_30sec$ID <- paste(Inj_int_ent_30sec$BirdID, Inj_int_ent_30sec$Condition, Inj_int_ent_30sec$time_bin_30sec, sep = "_")
Concave_hull_areas_3min_inj_30$ID <- paste(Concave_hull_areas_3min_inj_30$BirdID, Concave_hull_areas_3min_inj_30$Condition, Concave_hull_areas_3min_inj_30$time_bin_30sec, sep = "_")
Concave_hull_areas_3min_inj_30 <- rename(Concave_hull_areas_3min_inj_30, concave_area = area)
circular_variance_inj_30sec$ID <- paste(circular_variance_inj_30sec$BirdID, circular_variance_inj_30sec$Condition, circular_variance_inj_30sec$time_bin_30sec, sep = "_")
Inj_binned_speed_30sec$ID <- paste(Inj_binned_speed_30sec$BirdID, Inj_binned_speed_30sec$Condition, Inj_binned_speed_30sec$time_bin_30sec, sep = "_")

all_metrics_inj_30sec <- merge(Inj_binned_speed_30sec, circular_variance_inj_30sec[, c("ID", "head_var")], by = "ID", all.x = TRUE)
all_metrics_inj_30sec <- merge(all_metrics_inj_30sec, Inj_int_ent_30sec[, c("ID", "shannon_entropy")], by = "ID", all.x = TRUE)
all_metrics_inj_30sec <- merge(all_metrics_inj_30sec, Concave_hull_areas_3min_inj_30[, c("ID", "concave_area")], by = "ID", all.x = TRUE)

#Subset to first 3 mins 
all_metrics_inj_bin1_30 <- all_metrics_inj_30sec %>% filter(time_bin_30sec < 7)
#Remove lesion birds 
all_metrics_30sec_simp <- all_metrics_inj_bin1_30 %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes")

#PCA 
###PCA on all sleap metrics 
metrics_30sec_simp <- all_metrics_30sec_simp %>%
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
    all_metrics_30sec_simp %>%
      select(BirdID, Condition, time_bin_30sec, Inj_Lesion) %>%
      slice(rownames(metrics_30sec_simp) |> as.integer())
  )

loadings_30sec_simp <- as.data.frame(pca_res_30sec_simp$rotation)
loadings_30sec_simp

screeplot(pca_res_30sec_simp, type = "lines", main = "Scree Plot (Base R)")

variance <- pca_res_30sec_simp$sdev^2
pve <- variance / sum(variance)

scree_df <- data.frame(
  PC = 1:length(pve),
  variance = pve
)
#Scree plot ()
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

#Vector bi-plot ()
autoplot(pca_res_30sec_simp, data = all_metrics_30sec_simp,
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


# First 3 min PCA stats ---------------------------------------------------
pca_scores_30sec_simp$time_bin_30sec <- factor(pca_scores_30sec_simp$time_bin_30sec, ordered = TRUE)
Inj_PC1_model_bin1_30_simp <- glmmTMB(PC1 ~ Condition * Inj_Lesion * time_bin_30sec + (1|BirdID),      
                                      family = gaussian(link = "identity"), data = pca_scores_30sec_simp)

Anova(Inj_PC1_model_bin1_30_simp)

summary(Inj_PC1_model_bin1_30_simp)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_PC1_model_bin1_30, plot = TRUE)


### Raw points and model predictions plot ######### 
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




