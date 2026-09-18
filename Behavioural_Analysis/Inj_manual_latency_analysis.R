# --- Injected birds manual latencies analysis --  


# Load df 
all_coh <- read_xlsx("/Users/michaelabierman/Downloads/LAB/Daria/Results/data/manual/Saporin_batches/familiar_vs_unfamiliar_saporin_allbatches.xlsx")

all_coh_6 <-all_coh %>% filter(Behaviour != "Eating")
all_coh_6 <-all_coh_6 %>% filter(Behaviour != "Grooming")

excluded_birds <- c("bl85gy195","bl21or11", "bl48pu128", "bl191or41", "bl10or50", "bl25or15")

all_coh_6 <- all_coh_6 %>% 
  filter(!(BirdID %in% excluded_birds))

#Nest injection and lesion 
all_coh_6$Inj_Lesion <- paste(all_coh_6$Injection, all_coh_6$Lesion, sep = "_")  

#Subset to movements and vocalizations 
Inj_lat_move <- all_coh_6 %>% filter(Behaviour == "Long Hop" | Behaviour == "Short Hop"| Behaviour == "Beak Swipe")
Inj_lat_vocal <- all_coh_6 %>% filter(Behaviour == "Long Call" | Behaviour == "Short Call"| Behaviour == "Singing")

#count stimuli 
print(sum(Inj_lat_move$Condition =="F" & Inj_lat_move$Inj_Lesion =="IgG-SAP_No" & Inj_lat_move$Behaviour =="Short Hop"))
print(sum(Inj_lat_move$Condition =="N" & Inj_lat_move$Inj_Lesion =="IgG-SAP_No" & Inj_lat_move$Behaviour =="Short Hop"))
print(sum(Inj_lat_move$Condition =="F" & Inj_lat_move$Inj_Lesion =="anti-DBH-SAP_No" & Inj_lat_move$Behaviour =="Short Hop"))
print(sum(Inj_lat_move$Condition =="N" & Inj_lat_move$Inj_Lesion =="anti-DBH-SAP_No" & Inj_lat_move$Behaviour =="Short Hop"))


#MOVEMENT STATS -----------------------
hist(Inj_lat_move$Latencies)
hist(log(Inj_lat_move$Latencies))

#Including Lesion birds 
Inj_lat_move_model <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move
)

#Excluding Lesion birds
Inj_lat_move_model <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes")
)

Anova(Inj_lat_move_model)
summary(Inj_lat_move_model)


#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_move_model, plot = TRUE)
plot(residuals(Inj_lat_move_model))
qqnorm(resid(Inj_lat_move_model))
qqline(resid(Inj_lat_move_model))

#MOVEMENTS POST-HOC 
# excluding lesions model 
# emm comparisons: two-way condition:behaviour 
emm_beh <- emmeans(
  Inj_lat_move_model,
  ~ Condition | Behaviour,
  type = "response"   # back-transforms from log scale
)
pairs(emm_beh)
emm_move_igg <- emmeans(Inj_lat_move_model, ~ Condition | Behaviour, type = "response")
confint(pairs(emm_move_igg, adjust = "fdr", type = "response"))


# Post hoc models (for model including lesion birds)
# IgG #
Inj_lat_move_igg <- Inj_lat_move %>% filter(Inj_Lesion == "IgG-SAP_No")

Inj_lat_move_model_igg <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move_igg
)

Anova(Inj_lat_move_model_igg)
summary(Inj_lat_move_model_igg)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_move_model_igg, plot = TRUE)
plot(residuals(Inj_lat_move_model_igg))
qqnorm(resid(Inj_lat_move_model_igg))
qqline(resid(Inj_lat_move_model_igg))
pairs(emmeans(Inj_lat_move_model_igg, ~ Condition | Behaviour, adjust = "fdr"))

# DBH NO #
Inj_lat_move_no <- Inj_lat_move %>% filter(Inj_Lesion == "anti-DBH-SAP_No")

Inj_lat_move_model_no <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move_no
)

Anova(Inj_lat_move_model_no)
summary(Inj_lat_move_model_no)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_move_model_no, plot = TRUE)
plot(residuals(Inj_lat_move_model_no))
qqnorm(resid(Inj_lat_move_model_no))
qqline(resid(Inj_lat_move_model_no))
pairs(emmeans(Inj_lat_move_model_no, ~ Condition | Behaviour, adjust = "fdr"))

# DBH YES #
Inj_lat_move_yes <- Inj_lat_move %>% filter(Inj_Lesion == "anti-DBH-SAP_Yes")

Inj_lat_move_model_yes <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move_yes
)

Anova(Inj_lat_move_model_yes)
summary(Inj_lat_move_model_yes)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_move_model_yes, plot = TRUE)
plot(residuals(Inj_lat_move_model_yes))
qqnorm(resid(Inj_lat_move_model_yes))
qqline(resid(Inj_lat_move_model_yes))
pairs(emmeans(Inj_lat_move_model_yes, ~ Condition | Behaviour, adjust = "fdr"))






# Vocalization STATS --------------------------
hist(Inj_lat_vocal$Latencies)
hist(log(Inj_lat_vocal$Latencies))

#Including Lesions 
Inj_lat_vocal_model <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal
)

#Excluding Lesions 
Inj_lat_vocal_model <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes")
)


Anova(Inj_lat_vocal_model)
summary(Inj_lat_vocal_model)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_vocal_model, plot = TRUE)
plot(residuals(Inj_lat_vocal_model))
qqnorm(resid(Inj_lat_vocal_model))
qqline(resid(Inj_lat_vocal_model))

#POST-HOC 
# Excluding lesions (two-way trend condition:Inj_Lesion) 
# and including lesions (three way)

# IgG #
Inj_lat_vocal_igg <- Inj_lat_vocal %>% filter(Inj_Lesion == "IgG-SAP_No")

Inj_lat_move_vocal_igg <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_igg
)

Anova(Inj_lat_move_vocal_igg)
summary(Inj_lat_move_vocal_igg)


#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_move_vocal_igg, plot = TRUE)
plot(residuals(Inj_lat_move_vocal_igg))
qqnorm(resid(Inj_lat_move_vocal_igg))
qqline(resid(Inj_lat_move_vocal_igg))
pairs(emmeans(Inj_lat_move_vocal_igg, ~ Condition | Behaviour, adjust = "fdr"))
emm_vocal_igg <- emmeans(Inj_lat_move_vocal_igg, ~ Condition | Behaviour, type = "response")
confint(pairs(emm_vocal_igg, adjust = "fdr", type = "response"))




# DBH NO #
Inj_lat_vocal_no <- Inj_lat_vocal %>% filter(Inj_Lesion == "anti-DBH-SAP_No")

Inj_lat_vocal_model_no <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_no
)

Anova(Inj_lat_vocal_model_no)
summary(Inj_lat_vocal_model_no)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_vocal_model_no, plot = TRUE)
plot(residuals(Inj_lat_vocal_model_no))
qqnorm(resid(Inj_lat_vocal_model_no))
qqline(resid(Inj_lat_vocal_model_no))
pairs(emmeans(Inj_lat_vocal_model_no, ~ Condition | Behaviour, adjust = "fdr"))
emm_vocal_dbhno <- emmeans(Inj_lat_vocal_model_no, ~ Condition | Behaviour, type = "response")
confint(pairs(emm_vocal_dbhno, adjust = "fdr", type = "response"))




# DBH YES #
Inj_lat_vocal_yes <- Inj_lat_vocal %>% filter(Inj_Lesion == "anti-DBH-SAP_Yes")

Inj_lat_vocal_model_yes <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_yes
)

Anova(Inj_lat_vocal_model_yes)
summary(Inj_lat_vocal_model_yes)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_vocal_model_yes, plot = TRUE)
plot(residuals(Inj_lat_vocal_model_yes))
qqnorm(resid(Inj_lat_vocal_model_yes))
qqline(resid(Inj_lat_vocal_model_yes))
pairs(emmeans(Inj_lat_vocal_model_yes, ~ Condition | Behaviour, adjust = "fdr"))






# lats capped at 3 mins ---------------------------------------------------
dodge <- position_dodge(width = 0.6)

Inj_lat_move_3min <- Inj_lat_move %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes") %>%
  mutate(Latencies = pmin(Latencies, 18))


Inj_lat_vocal_3min <- Inj_lat_vocal %>% filter(Inj_Lesion != "anti-DBH-SAP_Yes") %>%
  mutate(Latencies = pmin(Latencies, 18))


#Plot move 
ggplot(Inj_lat_move_3min, aes(x = Inj_Lesion, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "IgG-SAP_No" = "IgG",
      "anti-DBH-SAP_No" = "DBH",
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
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
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))


#Stats move 
Inj_lat_move_model_3min <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_move_3min
)


Anova(Inj_lat_move_model_3min)

summary(Inj_lat_move_model_3min)


emm_beh_3 <- emmeans(
  Inj_lat_move_model_3min,
  ~ Condition | Behaviour,
  type = "response"   # back-transforms from log scale
)

pairs(emm_beh_3)
confint(pairs(emm_beh_3, adjust = "fdr", type = "response"))




# Vocal 3 mins ------------------------------------------------------------
#Plot vocal
ggplot(Inj_lat_vocal_3min, aes(x = Inj_Lesion, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "IgG-SAP_No" = "IgG",
      "anti-DBH-SAP_No" = "DBH",
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
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
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))


#Stats vocal 
Inj_lat_vocal_model_3min <- glmmTMB(
  Latencies ~ Condition * Inj_Lesion * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_3min
)


Anova(Inj_lat_vocal_model_3min)

summary(Inj_lat_vocal_model_3min)

## Post hoc move 3 mins 
##### IgG ##########
Inj_lat_vocal_igg_3min <- Inj_lat_vocal_3min %>% filter(Inj_Lesion == "IgG-SAP_No")

Inj_lat_vocal_igg_3min_model <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_igg_3min
)

Anova(Inj_lat_vocal_igg_3min_model)

summary(Inj_lat_vocal_igg_3min_model)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_vocal_igg_3min_model, plot = TRUE)

plot(residuals(Inj_lat_vocal_igg_3min_model))
qqnorm(resid(Inj_lat_vocal_igg_3min_model))
qqline(resid(Inj_lat_vocal_igg_3min_model))
pairs(emmeans(Inj_lat_vocal_igg_3min_model, ~ Condition | Behaviour, adjust = "fdr"))
emm_vocal_igg_3min <- emmeans(Inj_lat_vocal_igg_3min_model, ~ Condition | Behaviour, type = "response")
confint(pairs(emm_vocal_igg_3min, adjust = "fdr", type = "response"))



######### DBH NO #############
Inj_lat_vocal_no_3min <- Inj_lat_vocal_3min %>% filter(Inj_Lesion == "anti-DBH-SAP_No")
hist(Inj_lat_vocal_no_3min$Latencies)

#False covergence
Inj_lat_vocal_model_no_3min <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = Inj_lat_vocal_no_3min
)

#Converges with identical values to truncated_nbinom2
Inj_lat_vocal_model_no_3min <- glmmTMB(
  Latencies ~ Condition * Behaviour + (1|BirdID),
  family = poisson(link = "log"),  
  data = Inj_lat_vocal_no_3min
)
Anova(Inj_lat_vocal_model_no_3min)

summary(Inj_lat_vocal_model_no_3min)

#residuals
simulationOutput <- simulateResiduals(fittedModel = Inj_lat_vocal_model_no_3min, plot = TRUE)









# Plots # -----------
# MOVEMENTS 
Inj_lat_move$Inj_Lesion <- factor(
  Inj_lat_move$Inj_Lesion,
  levels = c("IgG-SAP_No", "anti-DBH-SAP_No", "anti-DBH-SAP_Yes")
)


Inj_lat_move$Behaviour <- factor(
  Inj_lat_move$Behaviour,
  levels = c("Short Hop", "Long Hop", "Beak Swipe")
)

dodge <- position_dodge(width = 0.6)

#Fig 5a
ggplot(Inj_lat_move%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes"), aes(x = Inj_Lesion, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 95),
    breaks = seq(0, 90, by = 15)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "IgG-SAP_No" = "IgG",
      "anti-DBH-SAP_No" = "DBH",
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
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
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))




# VOCALIZATIONS 
Inj_lat_vocal$Inj_Lesion <- factor(
  Inj_lat_vocal$Inj_Lesion,
  levels = c("IgG-SAP_No", "anti-DBH-SAP_No", "anti-DBH-SAP_Yes")
)


Inj_lat_vocal$Behaviour <- factor(
  Inj_lat_vocal$Behaviour,
  levels = c("Short Call", "Long Call", "Singing")
)


dodge <- position_dodge(width = 0.6)

#Fig 5b 
ggplot(Inj_lat_vocal%>%filter(Inj_Lesion!="anti-DBH-SAP_Yes"), aes(x = Inj_Lesion, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 95),
    breaks = seq(0, 90, by = 15)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "IgG-SAP_No" = "IgG",
      "anti-DBH-SAP_No" = "DBH",
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
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
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))


# Day 1 Analysis Figure 5 -------------------------------------------------
day_1 <- read.csv("/Users/michaelabierman/Downloads/LAB/Daria/Results/data/manual/Saporin_batches/familiar_vs_unfamiliar_experimental_day1_v3.csv")

exclude_birds <- c("bl85gy195", "bl21or11", "bl48pu128", "bl191or41", "bl10or50", "bl25or15")
day_1 <- day_1 %>% filter(!(BirdID %in% exclude_birds))
day_1$Inj_Lesion <- paste(day_1$Injection, day_1$Lesion, sep = "_")

day_1 <- day_1 %>%
  left_join(all_coh_6 %>% distinct(BirdID, Condition), by = "BirdID")

day_1 <- rename(day_1, D1_latencies = "Latencies")

day_1 <- day_1 %>%
  left_join(all_coh_6 %>% distinct(BirdID, Latencies, Behaviour), by = c("BirdID", "Behaviour"))

day_1 <- day_1 %>% mutate(lat_diff = Latencies - D1_latencies)

day_1 <- day_1 %>% mutate(
  D1_latencies_3min = case_when(
    D1_latencies >17 ~ 18,
    TRUE             ~ D1_latencies),
  Latencies_3min = case_when(
    Latencies >17 ~ 18,
    TRUE             ~ Latencies)
)


day1_DBH_SAP_lesionNo<-day_1%>%filter(Inj_Lesion=="anti-DBH-SAP_No")
unique(day1_DBH_SAP_lesionNo$BirdID)
day1_control<-day_1%>%filter(Inj_Lesion=="control_No")
unique(day1_control$BirdID)
#All DBH-SAPNo, missing 1 control no day 1 video 
day_1_move <- day_1 %>% filter(Behaviour %in% c("Long Hop", "Short Hop", "Beak Swipe"))
day_1_vocal <- day_1 %>% filter(Behaviour %in% c("Long Call", "Short Call", "Singing"))

# Diff in day 1 lats bw manip ---------------------------------------------
day_1_move$Behaviour <- factor(
  day_1_move$Behaviour, levels = c("Short Hop", "Long Hop", "Beak Swipe")
)
day_1_move$Inj_Lesion <- factor(
  day_1_move$Inj_Lesion, levels = c("control_No", "anti-DBH-SAP_No")
)

day_1_vocal$Behaviour <- factor(
  day_1_vocal$Behaviour, levels = c("Short Call", "Long Call", "Singing")
)
day_1_vocal$Inj_Lesion <- factor(
  day_1_vocal$Inj_Lesion, levels = c("control_No", "anti-DBH-SAP_No")
)


# plots  ------------------------------------------------------------------
ggplot(day_1_long_vocal%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes"), aes(x=day, y=lat_3min, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    aes(group=Behaviour),
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group=Behaviour),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    aes(group=Behaviour),
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~interaction(Condition,Inj_Lesion), nrow=1) +
  #facet_grid(rows=vars(Condition), cols=vars(Behaviour))+
  scale_x_discrete(
    labels = c(
      "D1_latencies" = "Exp",
      "Latencies" = "Test")
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    #axis.text.x =element_text(angle=45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))

ggplot(day_1_long_move%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes"), aes(x=interaction(Inj_Lesion,day,Condition), y=lat_3min, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    aes(group=Behaviour),
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group=Behaviour),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    aes(group=Behaviour),
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  #facet_wrap(~interaction(Condition,Inj_Lesion), nrow=1) +
  #facet_grid(rows=vars(Condition), cols=vars(Behaviour))+
  scale_x_discrete(
    labels = c(
      "D1_latencies" = "Exp",
      "Latencies" = "Test")
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    #axis.text.x =element_text(angle=45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))

ggplot(day_1_long_move%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes"), aes(x=interaction(Condition,Inj_Lesion), y=lat_3min, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    aes(group=day),
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group=day),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    aes(group=day),
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  #facet_grid(rows=vars(Condition), cols=vars(Behaviour))+
  scale_x_discrete(
    labels = c(
      "control_No" = "CON",
      "anti-DBH-SAP_No" = "DBH-SAP")
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    #axis.text.x =element_text(angle=45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))


ggplot(day_1_long_vocal%>%filter(Inj_Lesion != "anti-DBH-SAP_Yes"), aes(x=interaction(Condition,Inj_Lesion), y=lat_3min, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    aes(group=day),
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group=day),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    aes(group=day),
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  #facet_grid(rows=vars(Condition), cols=vars(Behaviour))+
  scale_x_discrete(
    labels = c(
      "control_No" = "CON",
      "anti-DBH-SAP_No" = "DBH-SAP")
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    #axis.text.x =element_text(angle=45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))

# stats -------------------------------------------------------------------
#Full models 

#Movements 
move_full_model <- glmmTMB(
  lat ~ day * Inj_Lesion * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_move%>%filter(Inj_Lesion!="anti-DBH-SAP_Yes")
)
Anova(move_full_model)

#Vocalizations 
vocal_full_model <- glmmTMB(
  lat ~ day * Inj_Lesion * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion!="anti-DBH-SAP_Yes")
)
Anova(vocal_full_model)

#Control
vocal_full_model_igg <- glmmTMB(
  lat ~ day * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion=="control_No")
)
Anova(vocal_full_model_igg)
pairs(emmeans(vocal_full_model_igg, ~ day | Condition , adjust = "fdr"))

#DBH
vocal_full_model_dbh <- glmmTMB(
  lat ~ day * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion=="anti-DBH-SAP_No")
)
Anova(vocal_full_model_dbh)

#Movements 3mins 
move_full_model_3min <- glmmTMB(
  lat_3min ~ day * Inj_Lesion * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_move%>%filter(Inj_Lesion!="anti-DBH-SAP_Yes")
)
Anova(move_full_model_3min)
pairs(emmeans(move_full_model_3min, ~ day | Condition , adjust = "fdr"))



#Vocalizations 3mins 
vocal_full_model_3min <- glmmTMB(
  lat_3min ~ day * Inj_Lesion * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion!="anti-DBH-SAP_Yes")
)
Anova(vocal_full_model_3min)

#Control
vocal_full_model_3min_igg <- glmmTMB(
  lat_3min ~ day * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion=="control_No")
)
Anova(vocal_full_model_3min_igg)
pairs(emmeans(vocal_full_model_3min_igg, ~ day | Condition , adjust = "fdr"))

#DBH 
vocal_full_model_3min_DBH <- glmmTMB(
  lat_3min ~ day * Behaviour * Condition + (1|BirdID),
  family = truncated_nbinom2(link = "log"),  
  data = day_1_long_vocal%>%filter(Inj_Lesion=="anti-DBH-SAP_No")
)
Anova(vocal_full_model_3min_DBH)








# Supplemental lesion latencies sup fig 3 -------------------------------------

ggplot(Inj_lat_move%>%filter(Inj_Lesion=="anti-DBH-SAP_Yes"), aes(x = Condition, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 95),
    breaks = seq(0, 90, by = 15)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))

ggplot(Inj_lat_vocal%>%filter(Inj_Lesion=="anti-DBH-SAP_Yes"), aes(x = Condition, y = Latencies, fill = Condition, color = Condition)) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    position = dodge,
    alpha = 0.7,
    width = 0.5
  ) +
  
  # Error bars
  stat_summary(
    aes(group = Condition),
    fun.data = mean_se,
    geom = "errorbar",
    position = dodge,
    width = 0.2
  ) +
  
  # Raw data
  geom_jitter(
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.6
    ),
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 95),
    breaks = seq(0, 90, by = 15)
  )+
  theme_minimal() +
  theme(
    text = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines"),
    panel.grid = element_blank(),
  ) +
  
  labs(
    x = NULL,
    y = "Latencies (Trial #)",
    fill = "Condition"
  ) +
  facet_wrap(~Behaviour) +
  scale_x_discrete(
    labels = c(
      "anti-DBH-SAP_Yes" = "DBH-Yes"
    )
  )+
  theme(
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("N" = "orange", "F" = "purple")) +
  scale_color_manual(values = c("N" = "orange3", "F" = "purple"))









