library(ggplot2)
library(scales)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(ggbeeswarm) # to make jitter plots 
library(patchwork) # for putting ggplot objects together 
library(readxl)
library(conflicted)
library(dplyr)
library(emmeans)
library(readxl)
library(car)
library(DHARMa)
library(multcomp)
library(multcompView)
library(patchwork)

library(RColorBrewer)


# 1. Prepare the data
# We set the "Down-regulated" counts as negative numbers so they plot below the x-axis.
df <- data.frame(
  Time = factor(c("5 min", "5 min", "10 min", "10 min", "30 min", "30 min"),
                levels = c("5 min", "10 min", "30 min")), 
  Regulation = factor(c("Up-regulated", "Down-regulated",
                        "Up-regulated", "Down-regulated",
                        "Up-regulated", "Down-regulated"),
                      levels = c("Up-regulated", "Down-regulated")), #
  Count = c(217, -1, 632, -41, 1239, -178)
)

# 2. Create the plot
ggplot(df, aes(x = Time, y = Count, fill = Regulation)) +
  # geom_col creates the bars; color="black" gives the crisp outlines
  geom_col(color = "black", width = 0.4) + 
  
  # Add the text labels above and below the bars
  # abs(Count) ensures we don't print negative signs for the down-regulated numbers
  geom_text(aes(label = abs(Count), 
                vjust = ifelse(Count > 0, -0.5, 1.5)), 
            size = 5.3) +
  
  # Customize colors (Red and Blue)
  scale_fill_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8")) +
  
  # Format the Y-axis: use abs() to hide negative signs on the axis ticks
  # Set the limits and breaks to match the specific ranges in your image
  scale_y_continuous(labels = abs, 
                     limits = c(-250, 1400), 
                     breaks = seq(-200, 1300, by = 300)) +
  
  # Add the bold zero-line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  
  # Clean up labels and apply a publication-ready theme
  labs(y = "Number of DEGs", x = NULL) +
  theme_classic() +
  theme(
    # Move the legend inside the plot area at the top
    legend.position = c(0.3, 0.9), 
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    legend.key.size = unit(0.7, "cm"),
    legend.direction = "vertical",
    
    # Adjust text sizes for manuscript readability
    axis.text = element_text(size = 15, color = "black"),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.line = element_line(color = "black"),
    
    # Remove the x-axis line since we drew our own at y=0
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank() 
  )
###W430, H310




######################### FigS2b_flowering
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=1)

p_values <- phenofill %>%
  group_by(geno) %>%
  summarise(
    p_value = t.test(height ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    # Format p-value
    p_label = case_when(
      p_value < 0.001 ~ "p < 0.001",
      p_value < 0.01 ~ sprintf("p = %.3f", p_value),
      p_value < 0.05 ~ sprintf("p = %.3f", p_value),
      TRUE ~ sprintf("p = %.2f", p_value)
    ),
    # Stars for significance
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    # Y position for label (adjust based on your data range)
    y_position = max(phenofill$height)+1
  )

print("P-values calculated:")
print(p_values)

phenofill$height <- as.numeric(phenofill$height)
means <- phenofill %>%
  group_by(geno, group) %>%
  summarise(mean_value = mean(height, na.rm = TRUE),
            .groups = "drop")

phenofill <- phenofill %>%
  mutate(group = factor(group, levels = c("UT", "T")))

phenofill <- phenofill %>%
  mutate(Treatment = factor(Treatment, levels = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")))

plot <- ggplot(phenofill, aes(x = geno, y = height, fill = Treatment, color=Treatment)) +
  
  # 1. BOXPLOT
  geom_boxplot(outlier.shape = NA,  # Don't show outliers (will use jitter)
               width = 0.6,
               color = "black",
               linewidth = 0.8,
               alpha = 0.9) +
  
  # 2. DATA POINTS (DOTS)
  geom_jitter(aes(fill=Treatment),
              position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.6),
              size = 2,          
              alpha = 0.6,
              shape = 21,
              color = "black") +
  
  # 4. SIGNIFICANCE BRACKETS
  geom_segment(data = p_values,
               aes(x = as.numeric(factor(geno)) - 0.2,
                   xend = as.numeric(factor(geno)) + 0.2,
                   y = y_position,
                   yend = y_position),
               inherit.aes = FALSE,
               size = 0.8,
               color = "black") +
  
  # 5. P-VALUE LABELS
  geom_text(data = p_values,
            aes(x = geno,
                y = y_position,
                label = p_label, vjust = -0.4),
            inherit.aes = FALSE,
            size = 4,
            fontface = "bold") +
  
  # 6. COLORS (ADJUST TO YOUR PREFERENCE)
  scale_fill_manual(values = c(
    "Col-0_UT"   = "#CAB2D6",
    "Col-0_T"    = "#6A3D9A",
    "rrtf1-1_UT" = "#B2DF8A",
    "rrtf1-1_T"  = "#33A02C",
    "rrtf1-2_UT" = "#B2DF8A",
    "rrtf1-2_T"  = "#33A02C"
  )) +
  scale_color_manual(values = c(
    "Col-0_UT"   = "#CAB2D6",
    "Col-0_T"    = "#6A3D9A",
    "rrtf1-1_UT" = "#B2DF8A",
    "rrtf1-1_T"  = "#33A02C",
    "rrtf1-2_UT" = "#B2DF8A",
    "rrtf1-2_T"  = "#33A02C"
  )) +
  
  # 7. Y-AXIS SCALE
  scale_y_continuous(
    limits = c(min(phenofill$height) - 2, max(phenofill$height) + 2),
    breaks = seq(0, ceiling(max(phenofill$height) +2 ), by = 3)
  ) +
  
  # 8. THEME
  theme_classic() +
  theme(
    # Axes
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 17, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 16, hjust = 0),
    # Panel
    panel.background = element_rect(fill = "white")
  ) +
  
  # 9. LABELS
  labs(x = "Genotype", y = "Flowering time (days)")
# Display the plot
print(plot)

ggsave("../final_boxplot.svg",
       plot,
       width = 5.0,
       height = 5.0,
       dpi = 300,
       bg = "white")



######FigS2C_Rosette leaf area
### stat
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(rosette ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = rosette, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-1" = "italic('rrtf1-1')",
      "rrtf1-2" = "italic('rrtf1-2')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$rosette)*1.1)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p = ", signif(p.value, 2))), size = 3.5, x = 1.5, y = 100, vjust = 1.75, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  labs(y = expression(bold("Rosette area ("*cm^2*")")),
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 17, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0),
    axis.title = element_text(size = 16, face = "bold")
  )

free_scale

ggsave("../free_scale.svg",plot = free_scale, height = 5, width = 5.1, bg = "white")

######FigS2D_ShoowFW
### stat
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(shootFW ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = shootFW, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-1" = "italic('rrtf1-1')",
      "rrtf1-2" = "italic('rrtf1-2')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$shootFW)*1.15)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p = ", signif(p.value, 2))), size = 3.5, x = 1.5, y = 100, vjust = 1.75, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  labs(y = "Shoot fresh weight (g)",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 17, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0),
    axis.title = element_text(size = 16, face = "bold")
  )

free_scale
ggsave("../free_scale.svg",plot = free_scale, height = 5, width = 5.1, bg = "white")

######FigS2E_rBr
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=2)
### stat
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(br ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))

phenofill <- phenofill %>%
  mutate(Treatment = factor(Treatment, levels = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = br, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-1" = "italic('rrtf1-1')",
      "rrtf1-2" = "italic('rrtf1-2')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$br)*1.15)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p = ", signif(p.value, 2))), size = 3.5, x = 1.5, y = 100, vjust = 1.75, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  labs(y = "Bud outgrowth (frequency)",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 17, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0),
    axis.title = element_text(size = 16, face = "bold")
  )

free_scale

######FigS3B_petiole
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=3)
### stat
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(petiole ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))

phenofill <- phenofill %>%
  mutate(Treatment = factor(Treatment, levels = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = petiole, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-1" = "italic('rrtf1-1')",
      "rrtf1-2" = "italic('rrtf1-2')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$petiole)*1.2)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T", "rrtf1-1_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  labs(y = "Petiole length (mm)",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 17, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0),
    axis.title = element_text(size = 16, face = "bold")
  )

free_scale
ggsave("../free_scale.svg",plot = free_scale, height = 5, width = 5.1, bg = "white")

ht.aov <- aov(petiole ~ group*geno, data = phenofill) #pick
summary(ht.aov)
shapiro.test(resid(ht.aov))
plot(ht.aov,1)
plot(ht.aov,2)
sr <- simulateResiduals(fittedModel = ht.aov, plot = T) #looks okay.
emm <- emmeans(ht.aov, ~group*geno, adjust = "mvt")
cld(emm, Letters = letters)



####################vFigS4_horomone####################v
# ============================================
horm <- read.csv("~/../thigmo_horm_tukey.csv",sep=",")
horm <- horm %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno))

summary_data <- horm %>%
  group_by(geno, time) %>%
  summarise(
    mean = mean(ABA, na.rm = TRUE),           # mean
    se = sd(ABA, na.rm = TRUE) / sqrt(n()),   # SE
    .groups = "drop")

print(summary_data)
summary_data$geno <- factor(summary_data$geno, 
                            levels = c("Col-0", "rrtf1"))
horm$time <- factor(horm$time, 
                    levels = c("0","2.5","5","10"))
# Tukey letters
horm.aov <- aov(ABA ~ geno*time, data = horm) #pick
summary(horm.aov)
shapiro.test(resid(horm.aov))
leveneTest(ABA ~ gbyt, data = horm)
plot(horm.aov,1)
plot(horm.aov,2)
sr <- simulateResiduals(fittedModel = horm.aov, plot = T) #looks okay.

emm <- emmeans(horm.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))

summary_data$time <- as.numeric(as.character(summary_data$time))
horm$time <- as.numeric(as.character(horm$time))

horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(data = horm,
             aes(x = time, y = ABA, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1)) + 
  geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean - se, mean + se), vjust = ifelse(geno == "Col-0", 1.5, -0.7)), size = 5, fontface = "bold", color = "black")  +
  
  geom_point(size = 5) +
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    
    "rrtf1" = "#33A02C"     
  ),
  labels = c(
    "Col-0" = expression(Col-0),
    "rrtf1" = expression(italic(rrtf1))
  )) +
  
  #axis
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10), labels = c("0", "2.5", "5", "7.5", "10")) +
  scale_y_continuous(limits = c(0, max(horm$ABA)*1.2)) +
  
  #Label
  labs(x = "Time (min)",
       y = "pmoles/gFW",
       color = "",
       title = "Abscisic acid (ABA)") + 
  theme(
    # Axes
    axis.text.x = element_text(size = 19, color = "black"),
    axis.text.y = element_text(size = 19, color = "black"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 23, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 25, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white"))

print(horm_plot)
ggsave("../final_boxplot.svg", horm_plot, width = 4.8, height = 6.08, dpi = 300, bg = "white")



# ============================================IAA============================================
horm <- read.csv("~/../thigmo_horm_tukey.csv",sep=",")
horm <- horm %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno))

summary_data <- horm %>%
  group_by(geno, time) %>%
  summarise(
    mean = mean(IAA, na.rm = TRUE),           # mean
    se = sd(IAA, na.rm = TRUE) / sqrt(n()),   # SE
    .groups = "drop")

print(summary_data)
summary_data$geno <- factor(summary_data$geno, 
                            levels = c("Col-0", "rrtf1"))
horm$time <- factor(horm$time, 
                    levels = c("0","2.5","5","10"))
# Tukey letters
horm.aov <- aov(IAA ~ geno*time, data = horm) #pick
summary(horm.aov)
shapiro.test(resid(horm.aov))
leveneTest(IAA ~ gbyt, data = horm)
plot(horm.aov,1)
plot(horm.aov,2)
sr <- simulateResiduals(fittedModel = horm.aov, plot = T) #looks okay.

emm <- emmeans(horm.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))

summary_data$time <- as.numeric(as.character(summary_data$time))
horm$time <- as.numeric(as.character(horm$time))

horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(data = horm,
             aes(x = time, y = IAA, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1)) + 
  geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean - se, mean + se), vjust = ifelse(geno == "Col-0", 1.5, -0.7)), size = 5, fontface = "bold", color = "black")  +
  
  geom_point(size = 5)+
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    
    "rrtf1" = "#33A02C"     
  ),
  labels = c(
    "Col-0" = expression(Col-0),
    "rrtf1" = expression(italic(rrtf1))
  )) +
  
  #axis
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10), labels = c("0", "2.5", "5", "7.5", "10")) +
  scale_y_continuous(limits = c(0, max(horm$IAA, na.rm = TRUE)*1.2)) +
  
  #Label
  labs(x = "Time (min)",
       y = "pmoles/gFW",
       color = "",
       title = "Indole-3-acetic acid (IAA)") + 
  theme(
    # Axes
    axis.text.x = element_text(size = 19, color = "black"),
    axis.text.y = element_text(size = 19, color = "black"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 23, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 25, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white"))

print(horm_plot)
ggsave("../final_boxplot.svg", horm_plot, width = 4.8, height = 6.08, dpi = 300, bg = "white")

# ============================================GAs 4, 1, 20, 8============================================
horm <- read.csv("~/../thigmo_horm_tukey.csv",sep=",")
horm <- horm %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno))

summary_data <- horm %>%
  group_by(geno, time) %>%
  summarise(
    mean = mean(GA20, na.rm = TRUE),           # mean
    se = sd(GA20, na.rm = TRUE) / sqrt(n()),   # SE
    .groups = "drop")

print(summary_data)
summary_data$geno <- factor(summary_data$geno, 
                            levels = c("Col-0", "rrtf1"))
horm$time <- factor(horm$time, 
                    levels = c("0","2.5","5","10"))
# Tukey letters
horm.aov <- aov(GA20 ~ geno*time, data = horm) #pick
summary(horm.aov)
shapiro.test(resid(horm.aov))
leveneTest(GA20 ~ gbyt, data = horm)
plot(horm.aov,1)
plot(horm.aov,2)
sr <- simulateResiduals(fittedModel = horm.aov, plot = T) #looks okay.

emm <- emmeans(horm.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))
tukey_all

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))

summary_data$time <- as.numeric(as.character(summary_data$time))
horm$time <- as.numeric(as.character(horm$time))

horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(data = horm,
             aes(x = time, y = GA20, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1)) + 
  geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean - se, mean + se), vjust = ifelse(geno == "Col-0", 1.5, -0.7)), size = 5, fontface = "bold", color = "black")  +
  
  geom_point(size = 5)+
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    
    "rrtf1" = "#33A02C"     
  ),
  labels = c(
    "Col-0" = expression(Col-0),
    "rrtf1" = expression(italic(rrtf1))
  )) +
  
  #axis
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10), labels = c("0", "2.5", "5", "7.5", "10")) +
  scale_y_continuous(limits = c(0, max(horm$GA20, na.rm = TRUE)*1.2)) +
  
  #Label
  labs(x = "Time (min)",
       y = "pmoles/gFW",
       color = "",
       title = expression("Gibberellin A20 (GA"[20]*")")) + 
  theme(
    # Axes
    axis.text.x = element_text(size = 19, color = "black"),
    axis.text.y = element_text(size = 19, color = "black"),
    axis.title = element_text(size = 19, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 23, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 25, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white"))

print(horm_plot)
ggsave("../final_boxplot.svg", horm_plot, width = 4.8, height = 6.08, dpi = 300, bg = "white")



##########################################FigS5a###############################################################
library("glmmTMB")
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=4)

pheno <- pheno %>%
  mutate(treat = factor(treat, levels = c("Mock", "MeJA")))
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

set.glmm <- glmmTMB(bw ~ geno * group * treat, 
                    family = tweedie(link = "log"), data = pheno)
shapiro.test(resid(set.glmm))

sr <- simulateResiduals(fittedModel = set.glmm, plot = T) #looks okay.

emm <- emmeans(set.glmm, ~geno*group*treat, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all

label_pos <- pheno %>%
  group_by(geno, group, treat) %>%
  summarise(y_pos = mean(bw, na.rm = TRUE) + sd(bw, na.rm = TRUE)/sqrt(n()) + max(pheno$bw) * 0.17,
            .groups = "drop") #max(pheno$bw) * 0.05 

tukey_all <- left_join(tukey_all, label_pos, by = c("geno", "group", "treat"))

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "Col-0", ]

phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

phenofill <- phenofill %>%
  mutate(Treatment = factor(Treatment, levels = c("Col-0_UT_Mock", "Col-0_T_Mock","Col-0_UT_MeJA","Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_tT_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")))

##group+treat
free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = bw, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("Col-0" = "Col-0")
    ,label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$bw)*1.1)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "Col-0", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black") +
  
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = "Leaf length-to-width ratio",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )
free_scale

##############rrtf1######################################################################
phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "rrtf1", ]
phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

free_scale1 <- phenofill %>% 
  ggplot(aes(x = group, y = bw, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("rrtf1" = expression(italic("RRTF1")))
    ,label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$bw)*1.2)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "rrtf1", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black")+
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = "",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_blank(),  # 
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )

combined_plot <- free_scale + free_scale1 + 
  plot_layout(ncol = 2, widths = c(1, 1))
combined_plot
ggsave("../free_scale.svg",plot = combined_plot, width = 5.6, height = 4.0, bg = "white")


##########################################FigS5b##############################################################
library("glmmTMB")
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=4)

pheno <- pheno %>%
  mutate(treat = factor(treat, levels = c("Mock", "MeJA")))
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

set.glmm <- glmmTMB(avgleaf ~ geno * group * treat, 
                    family = tweedie(link = "log"), dispformula = ~ geno * group * treat, data = pheno)
shapiro.test(resid(set.glmm))
sr <- simulateResiduals(fittedModel = set.glmm, plot = T) #looks okay.

emm <- emmeans(set.glmm, ~geno*group*treat, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all

label_pos <- pheno %>%
  group_by(geno, group, treat) %>%
  summarise(y_pos = mean(avgleaf, na.rm = TRUE) + sd(avgleaf, na.rm = TRUE)/sqrt(n()) + max(pheno$avgleaf) * 0.21,
            .groups = "drop") #max(pheno$bw) * 0.05 

tukey_all <- left_join(tukey_all, label_pos, by = c("geno", "group", "treat"))

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "Col-0", ]

phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

##group+treat
free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = avgleaf, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("Col-0" = "Col-0")
    ,label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$avgleaf)*1.1)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "Col-0", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black") +
  
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = expression(bold("Average rosette area ("*cm^2*")")),
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )
free_scale

##############rrtf1######################################################################
phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "rrtf1", ]
phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

free_scale1 <- phenofill %>% 
  ggplot(aes(x = group, y = avgleaf, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("rrtf1" = expression(italic("RRTF1")))
    ,label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$avgleaf)*1.2)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "rrtf1", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black")+
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = "",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_blank(),  # 
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )

combined_plot <- free_scale + free_scale1 + 
  plot_layout(ncol = 2, widths = c(1, 1))
combined_plot
ggsave("../free_scale.svg",plot = combined_plot, width = 5.6, height = 4.0, bg = "white")


######################### FigS5c_flowering
library(glmmTMB)
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=4)

pheno <- pheno %>%
  mutate(treat = factor(treat, levels = c("Mock", "MeJA")))
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

set.glmm <- glmmTMB(flower ~ geno * group * treat, 
                    family = (link = "log"), data = pheno)
shapiro.test(resid(set.glmm))
sr <- simulateResiduals(fittedModel = set.glmm, plot = T) #looks okay.

emm <- emmeans(set.glmm, ~geno*group*treat, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all

label_pos <- pheno %>%
  group_by(geno, group, treat) %>%
  summarise(y_pos = mean(flower, na.rm = TRUE) + sd(flower, na.rm = TRUE)/sqrt(n()) + max(pheno$flower) * 0.15,
            .groups = "drop") #max(pheno$flower) * 0.05 

tukey_all <- left_join(tukey_all, label_pos, by = c("geno", "group", "treat"))

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "Col-0", ]
phenofill <- droplevels(phenofill)
phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = flower, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("Col-0" = "Col-0")
    ,label_parsed)) + coord_cartesian(ylim = c(min(phenofill$flower)*0.9, max(phenofill$flower)*1.1)) +
  geom_boxplot(aes(fill = Treatment, color = Treatment), outlier.shape = NA,  # Don't show outliers (will use jitter)
               width = 0.6,
               color = "black",
               linewidth = 0.8,
               alpha = 0.9)  +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "Col-0", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black") +
  # 7. Y-AXIS SCALE
  scale_y_continuous(
    limits = c(min(phenofill$flower)*0.9, max(phenofill$flower)*1.1),
    breaks = seq(20,33, by = 2)
  ) +
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = "Flowering time (days)",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )
free_scale

##############rrtf1######################################################################
phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "rrtf1", ]
phenofill$group <- factor(phenofill$group, levels = c("UT", "T"))
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))
phenofill <- phenofill[phenofill$geno == "rrtf1", ]
phenofill <- droplevels(phenofill)

free_scale1 <- phenofill %>% 
  ggplot(aes(x = group, y = flower, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("rrtf1" = expression(italic("RRTF1")))
    ,label_parsed)) + coord_cartesian(ylim = c(min(phenofill$flower)*0.9, max(phenofill$flower)*1.1)) +
  geom_boxplot(aes(fill = Treatment, color = Treatment), outlier.shape = NA,  # Don't show outliers (will use jitter)
               width = 0.6,
               color = "black",
               linewidth = 0.8,
               alpha = 0.9)  +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = tukey_all[tukey_all$geno == "rrtf1", ],
            aes(x = group, y = y_pos, label = .group),
            inherit.aes = FALSE,
            size = 5, fontface = "bold", color = "black") +
  # 7. Y-AXIS SCALE
  scale_y_continuous(
    limits = c(min(phenofill$flower)*0.9, max(phenofill$flower)*1.1),
    breaks = seq(20, 33, by = 2)
  ) +
  scale_fill_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT_Mock"   = "#CAB2D6",
      "Col-0_T_Mock"    = "#6A3D9A",
      "Col-0_UT_MeJA" = "#FB9A84",
      "Col-0_T_MeJA"  = "#E31A1C",
      "rrtf1_UT_Mock"   = "#B2DF8A",
      "rrtf1_T_Mock"    = "#33A02C",
      "rrtf1_UT_MeJA" = "#FDBF6F",
      "rrtf1_T_MeJA"  = "#FF7F00"
    ),
    breaks = c("Col-0_UT_Mock", "Col-0_T_Mock", "Col-0_UT_MeJA", "Col-0_T_MeJA", "rrtf1_UT_Mock", "rrtf1_T_Mock", "rrtf1_UT_MeJA", "rrtf1_T_MeJA")
  ) +
  labs(y = "",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_blank(),  # 
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )

combined_plot <- free_scale + free_scale1 + 
  plot_layout(ncol = 2, widths = c(1, 1))
combined_plot

ggsave("../free_scale.svg",plot = combined_plot, width = 5.6, height = 4.0, bg = "white")

################################################FigS6_CBF#################################################
pheno <- read_excel("~/../Touch/Remake/supp/FigS_data.xlsx", sheet=5)
pheno %>%
  filter(geno == "Col-0" & group == "UT") 
pheno %>%
  filter(geno == "rrtf1" & group == "UT")
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))
phenofill <- phenofill %>%
  mutate(group = factor(group, levels = c("UT", "T")))


free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = mean, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1" = "italic('rrtf1')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, max(phenofill$mean)*1.35)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.8, color = c(
    "Col-0_UT"   = "#CAB2D6",
    "Col-0_T"    = "#B2DF8A",
    "rrtf1_UT" = "#6A3D9A", #It only works if #Col-0_T and rrtf1_UT are swapped.
    "rrtf1_T"  = "#33A02C"
  ),alpha=0.7)  +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1_UT" = "#B2DF8A",
      "rrtf1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1_UT", "rrtf1_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1_UT" = "#B2DF8A",
      "rrtf1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1_UT", "rrtf1_T")
  ) +
  labs(y = "Relative expression",
       x = "Touch",
       title = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    title = element_text(size = 14, face = "bold"),
    plot.caption= element_text(size = 12, hjust = 0)
  )

free_scale

ggsave("../free_scale.svg",plot = free_scale,width = 3, height = 4, bg = "white")


####FigS7b
##makeup missing
setwd("~/../Thigmo_exp/2025_RNA_seq/rrtf1_RNAseq_033125/Cluster")
merge_top <- read.csv("merge_topwrky.csv", sep=",")
rownames(merge_top) <- merge_top[,2]
merge_top <- merge_top[,-1:-2]
merge_top <- merge_top[,c(2,1)]
colnames(merge_top) <- c( "Col-0_10m","rrtf1_10m")

####WRKY expression comparison
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(0, 5, 10), c("darkblue", "white", "red"))
row_labels <- rownames(merge_top)
col_fun <- colorRamp2(seq(0, 10, length.out = 9), brewer.pal(9, "YlGnBu"))

Heatmap(merge_top,
        name = "LOG2",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontface = "bold.italic"),
        column_names_rot = 45,
        heatmap_legend_param = list(title = "LOG2"),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "black", fill = NA, lwd = 1.5))
        }
)
ggsave("../free_scale.svg",plot = free_scale,width = 2.7, height = 5, bg = "white")
