#Script to summarize and analyse quantification of histology images in NCM and LC 

library(ggplot2)
library(readxl)
library(lme4)
library(lmerTest)  
library(emmeans)
library(readxl)
library(dplyr)
library(tidyr)
library(DHARMa)
library(coin)
library(car)
library(performance)
library(RVAideMemoire)
library(glmmTMB)

birds_to_exclude <- c("bl40pi140", "bl116pi160","bl85gy195",
                      "bl21or11", "bl104or24","bl23gr3", "bl186gr54", "bl48pu128",
                      "bl191or41", "bl25or15", "bl10or50")

#  NCM ---------------------------------------------------------------------####
#Import NCM batch1-6 file and remove birds with Lesion
df_NCM<- read_excel("DBH_NCM_batch123456.xlsx")%>%
  filter(!BirdID %in% birds_to_exclude)
df_NCM$Injection <- dplyr::recode(
  df_NCM$Injection,
  "Control" = "IgG-SAP"
)

#Remove birds with Lesion
df_NCM_no_lesion <- df_NCM %>%
  filter(Lesion == "No")

#correlation of Michaela and Jons counts
cor(df_NCM_no_lesion$`Jon counts`, df_NCM_no_lesion$`Michaela counts`, method = "spearman")
#0.82 STRONG

# Compute mean DBH+ innervation value per bird from all sections
bird_mean_NCM<- df_NCM_no_lesion %>%
  group_by(BirdID, Injection, Batch) %>%
  summarise(
    bird_mean = mean(`Averaged counts`, na.rm = TRUE),
    .groups = "drop"
  )

#to order control left and anti-DBH-SAP right
bird_mean_NCM$Injection <- factor(
  bird_mean_NCM$Injection,
  levels = c("IgG-SAP", "Anti-DBH-SAP")
)

#section counts per bird ID
df_NCM %>%
  distinct(BirdID, Section) %>%
  count(BirdID, name = "n_sections") %>%
  print(n=22)


#transform Michaela and Jon columns into rows and add to a new sheet
df_NCM_MJ <- df_NCM %>%
  pivot_longer(
    cols = c(`Jon counts`, `Michaela counts`),
    names_to = "Rater",
    values_to = "Count"
  ) %>%
  mutate(
    Rater = factor(Rater),
    Injection = factor(Injection, levels = c("IgG-SAP", "Anti-DBH-SAP")),
    BirdID = factor(BirdID),
    Section = factor(Section)
  )


#barplot IgG-SAP vs. Anti-DBH-SAP (No lesion) Fig 2C
ggplot(bird_mean_NCM, aes(x = Injection, y = bird_mean)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    fill = "grey80",
    color = "black",
    width = 0.65
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2
  ) +
  geom_jitter(
    width = 0.15,
    size = 2.2,
    alpha = 0.8,
    color = "black"
  ) +
  scale_x_discrete(labels = c("IgG-SAP" = "IgG", "Anti-DBH-SAP" = "DBH-SAP")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "DBH+ Innervation",
    title = "DBH+ Innervation in NCM, by Injection"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold") 
  )


####generalized mixed effects model on raw NCM data####
glmtmb_NCM <- glmmTMB(
  bird_mean ~ Injection + (1 | Batch),
  family = Gamma (link = "log"),
  data = bird_mean_NCM
)

Anova(glmtmb_NCM, type = 2)
emm <- emmeans(glmtmb_NCM, ~ Injection, type = "response")
emm

pairs(emm, adjust = "fdr")

qqnorm(residuals(glmtmb_NCM))
qqline(residuals(glmtmb_NCM))
simulationOutput <- simulateResiduals(fittedModel = glmtmb_NCM, plot = TRUE)
#Some flags from Dharma are ok due to small sample size






# LC------------------------------------------------------------------------####

#import LC batch1-6 file
df_LC <- read_excel("DBH_LC_batch123456.xlsx")%>%
  filter(!BirdID %in% birds_to_exclude)
#change control label to IgG-SAP
df_LC$Injection <- dplyr::recode(
  df_LC$Injection,
  "Control" = "IgG-SAP"
)

#remove birds with Lesion
df_LC_no_lesion <- df_LC %>%
  filter(Lesion == "No")

# Average DBH cell count per bird per injection
bird_mean_LC <- df_LC_no_lesion %>%
  group_by(BirdID, Injection, Batch, Lesion) %>%
  summarise(
    mean_DBH = mean(DBH_cell_count, na.rm = TRUE),
    .groups = "drop")

#to order control left and anti-DBH-SAP right
bird_mean_LC$Injection <- factor(
  bird_mean_LC$Injection,
  levels = c("IgG-SAP", "Anti-DBH-SAP")
)


#show number of birdID
df_LC %>%
  distinct(BirdID) %>%
  count(BirdID) %>%
  print(n = 22)

#shows section counts per bird ID
df_LC %>%
  distinct(BirdID, Section) %>%
  count(BirdID, name = "n_sections") %>%
  print(n = 22)

#add section number 
section_counts <- df_LC %>%
  distinct(BirdID, Section) %>%   # ensures unique sections
  count(BirdID, name = "n_sections")
bird_mean_LC <- bird_mean_LC %>%
  left_join(section_counts, by = "BirdID")


# barplot IgG-SAP vs. Anti-DBH-SAP Fig2E
ggplot(bird_mean_LC, aes(x = Injection, y = mean_DBH)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    fill = "grey80",
    color = "black",
    width = 0.65
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2
  ) +
  geom_jitter(
    width = 0.15,
    size = 2.2,
    alpha = 0.8,
    color = "black"
  ) +
  scale_x_discrete(labels = c("IgG-SAP" = "IgG", "Anti-DBH-SAP" = "DBH-SAP")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "# of DBH+ Neurons",
    title = "Number of DBH+ Neurons in LC, by Injection"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold") 
  )


#### generalized mixed effects model on LC count data ####
glmtmb_LC <- glmmTMB(
  mean_DBH ~ Injection + (1 | Batch),
  family = Gamma (link = "log"),
  data = bird_mean_LC
)

Anova(glmtmb_LC, type = 2)
emm <- emmeans(glmtmb_LC, ~ Injection, type = "response")
emm
pairs(emm, adjust = "fdr")
qqnorm(residuals(glmtmb_LC))
qqline(residuals(glmtmb_LC))
simulationOutput <- simulateResiduals(fittedModel = glmtmb_LC, plot = TRUE)
check_overdispersion(glmtmb_LC)



# LC AREA-------------------------------------------------------------------####

df_LC_area <- read_excel("LC_area.xlsx")

#add to excel with birds without Lesion
df_LC_no_lesion <- left_join(
  df_LC_no_lesion,
  df_LC_area,
  by = c("BirdID", "Section")
)


#average LC area per bird
bird_mean_area_LC <- df_LC_no_lesion %>%
  group_by(BirdID, Injection, Batch.x, Lesion) %>%
  summarise(
    mean_DBH_area = mean(area, na.rm = TRUE),
    .groups = "drop")

#to order control left and anti-DBH-SAP right
bird_mean_area_LC$Injection <- factor(
  bird_mean_area_LC$Injection,
  levels = c("IgG-SAP", "Anti-DBH-SAP")
)


#generalized mixed effects model for LC area on DBH-area averaged per bird
glmtmb_LC_area <- glmmTMB(
  mean_DBH_area ~ Injection + (1|Batch.x),
  family = Gamma (link = "log"),
  data = bird_mean_area_LC
)

Anova(glmtmb_LC_area, type = 2)
emm <- emmeans(glmtmb_LC_area, ~ Injection, type = "response")
emm
pairs(emm, adjust = "fdr")
qqnorm(residuals(glmtmb_LC_area))
qqline(residuals(glmtmb_LC_area))
simulationOutput <- simulateResiduals(fittedModel = glmtmb_LC_area, plot = TRUE)
check_overdispersion(glmtmb_LC_area)

# PLOT DBH AREA Fig2F
ggplot(bird_mean_area_LC, aes(x = Injection, y = mean_DBH_area)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    fill = "grey80",
    color = "black",
    width = 0.65
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2
  ) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.8) +
  scale_x_discrete(labels = c("IgG-SAP" = "IgG", "Anti-DBH-SAP" = "DBH-SAP")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "DBH+ Area (µm²)",
    title = "DBH+ Area in LC, by Injection"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold") 
  )
  




# NeuN ---------------------------------------------------------------------####
df_NeuN<- read_excel("all_batches_NeuN_NCM.xlsx")%>%
  filter(!BirdID %in% birds_to_exclude)

#compute per bird means for bars and error bars (since right now I have averaged DBH+ innervation per section but not PER BIRD)
bird_mean_NeuN<- df_NeuN %>%
  group_by(BirdID, Injection, Batch, Lesion) %>%
  summarise(
    bird_mean = mean(Area, na.rm = TRUE),
    .groups = "drop"
  )

#check if data is normally distributed
#check: rough symmetry, no extreme skew, no long tails
ggplot(bird_mean_NeuN, aes(x = bird_mean)) +
  geom_histogram(bins = 10, fill = "grey80", color = "black") +
  theme_classic(base_size = 14) +
  labs(x = "NeuN Area", y = "Count")
#noramllyish distributed?

#create PlotGrop
bird_mean_NeuN <- bird_mean_NeuN %>%
  mutate(
    PlotGroup = case_when(
      Injection == "IgG-SAP" ~ "IgG-SAP",
      Injection == "Anti-DBH-SAP" & Lesion == "No"  ~ "Anti-DBH-SAP (No lesion)",
      Injection == "Anti-DBH-SAP" & Lesion == "Yes" ~ "Anti-DBH-SAP (Lesion)"
    ),
    PlotGroup = factor(
      PlotGroup,
      levels = c(
        "IgG-SAP",
        "Anti-DBH-SAP (No lesion)",
        "Anti-DBH-SAP (Lesion)"
      )
    )
  )

glmtmb_NeuN <- glmmTMB(
  bird_mean ~ PlotGroup + (1 | Batch),
  family = tweedie,
  data = bird_mean_NeuN
)

Anova(glmtmb_NeuN, type = 2)
emm <- emmeans(glmtmb_NeuN, ~ PlotGroup, type = "response")
emm
pairs(emm, adjust = "fdr")
qqnorm(residuals(glmtmb_NeuN))
qqline(residuals(glmtmb_NeuN))
simulationOutput <- simulateResiduals(fittedModel = glmtmb_NeuN, plot = TRUE)
check_overdispersion(glmtmb_NeuN)


# barplot Fig S3B
ggplot(bird_mean_NeuN, aes(x = PlotGroup, y = bird_mean)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    fill = "grey80",
    color = "black",
    width = 0.65
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2
  ) +
  geom_jitter(
    width = 0.15,
    size = 2.2,
    alpha = 0.8,
    color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "NeuN Area (mm²)",
    title = "NeuN Area across Groups"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
    plot.title = element_text(hjust = 0.5)
  )