library(ggplot2)
library(scales)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(ggbeeswarm) 
library(patchwork) 
library(readxl)
library(dplyr)

library(emmeans) 

####Fig.1C
library(scales) #for pretty_breaks
GO_TF <- read.csv(file = "~/../TF_enrichment/GO_TF.csv",header=T, stringsAsFactors = F,sep =",")
GO_gp <- GO_TF
# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(GO_TF)
# Transform the column 'Gene_number' into a numeric variable
GO_TF$Gene_number <- as.numeric(GO_TF$Gene_number)
GO_TF$Pathway <- as.factor(GO_TF$Pathway)
# Transform FDR values by -log10('FDR values')
GO_TF$'|log10(FDR)|' <- -(log10(GO_TF$FDR))

# Change factor order
GO_TF$Group<- factor(GO_TF$Group,levels = c("5min","10min","30min"))
GO_TF$Pathway <- factor(GO_TF$Pathway,levels=c("ERF", "WRKY"))
group.labs <- c(`ERF` = "ERF TF",
                `WRKY` = "WRKY TF")
GO_TF <- GO_TF[1:4,]
# Draw the plot in facets by group with ggplot2
# to represent -log10(FDR), Number of genes and 
# Fold enrichment of each GO biological process per group (Figure 3)
#-------------------------------------------------------------------
ggplot(GO_TF, aes(x = Group, y = Pathway)) +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "azure4", linewidth=.5)+
  geom_point(data=GO_TF,aes(x=Group, y=Fold_enrichment,size = Gene_number, colour = `|log10(FDR)|`), alpha=.9)+
  scale_color_distiller(palette = "YlGnBu", direction = 1,
                        name = '|log10(FDR)|',
                        limits = c(0, NA)) +
  scale_y_continuous(breaks= pretty_breaks(), limits=c(1, 8))+
  scale_size_continuous(range = c(3, 7))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"), size = 12),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"), size = 12),
        axis.text = element_text(color = "black", face="bold"),
        axis.title.x = element_text(color = "black", face = "bold", size = 13),
        axis.title.y = element_text(color = "black", face = "bold", size = 13),
        panel.grid.minor = element_blank(),
        legend.title.align=0.5)+
  xlab("Time after touch")+
  ylab("Fold Enrichment")+
  labs(color="-log10(FDR)", size="Number\nof genes", face="bold")+
  facet_wrap(~Pathway,ncol=2,labeller=as_labeller(group.labs))+ #after "ncol=", specify the number of groups you have
  facet_grid(~Pathway, scales = "free", labeller=as_labeller(group.labs))+ 
  theme(strip.background = element_rect(color="black", fill="black", size=1, linetype="solid"), 
        strip.text.x = element_text(size = 10, color = "white", face = "bold"))+
  #scales each level of x but with different scales of Y-axis chosen by R freely
  guides(size = guide_legend(order=2),
         colour = guide_colourbar(order=1))

##############################Fig 1D
##makeup missing
merge_top <- read.csv("~/../Merge_bam/merge_top20each.csv", sep=",")

merge_top <- merge_top[,2:5]
rownames(merge_top) <- merge_top[,1]
merge_top <- merge_top[,-1]
colnames(merge_top) <- c("a_5min", "b_10min", "c_30min")
at_metadata <- data.frame(
  Time = c(5, 10, 30),
  row.names = c("a_5min", "b_10min", "c_30min")
)
gene_info <- read.csv("~/../Touch/Remake/gene_info.csv", sep=",")

merge_top$ensembl_gene_id <- ""
merge_top$ensembl_gene_id <- rownames(merge_top)
merge_top <- merge(merge_top, gene_info, by.x = "ensembl_gene_id", all.x = T)
rownames(merge_top) <- merge_top$external_gene_name
merge_top <- merge_top[,2:4]


col_fun <- colorRamp2(c(0, 5, 10), 
                      c(ylgnbu[1], ylgnbu[5], ylgnbu[9]))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
ylgnbu <- brewer.pal(9, "YlGnBu")
col_fun <- colorRamp2(c(0, 5, 10), c(ylgnbu[1], ylgnbu[5], ylgnbu[9]))
row_labels <- rownames(merge_top)
gene_colors <- ifelse(row_labels %in% c("RRTF1", "AT4G27654", "DDF1", "CBF4"), "#6A3D9A",
                      ifelse(row_labels == "bHLH19", "blue", "black"))
colnames(merge_top) <- c("5 min", "10 min", "30 min")
merge_top <- as.matrix(merge_top)
Heatmap(merge_top,
        name = "LOG2",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(col = gene_colors, fontface = "bold.italic"),
        column_names_rot = 45,
        heatmap_legend_param = list(title = "LOG2"),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "black", fill = NA, lwd = 1.5))
        }
)

##############################Fig2

pheno <- read_excel("~/../Touch/Remake/Fig2_data.xlsx", sheet=1)
### stat
library(dplyr)
pheno %>%
  filter(geno == "Col-0" & group == "UT") %>%
  summarise(mean = mean(height))
pheno %>%
  filter(geno == "rrtf1-1" & group == "UT") %>%
  summarise(mean = mean(height))
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(height ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = height, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-1" = "italic('rrtf1-1')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, 120)) +
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
      "rrtf1-1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T")
  ) +
  labs(y = "Shoot height (mm)",
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
#####################Fig2C#####################
pheno %>%
  filter(geno == "Col-0" & group == "UT") %>%
  summarise(mean = mean(rosette))
pheno %>%
  filter(geno == "rrtf1-1" & group == "UT") %>%
  summarise(mean = mean(rosette))
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
      "rrtf1-1" = "italic('rrtf1-1')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, 50)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p = ", signif(p.value, 2))), size = 3.5, x = 1.5, y = 50, vjust = 1.75, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-1_UT" = "#B2DF8A",
      "rrtf1-1_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-1_UT", "rrtf1-1_T")
  ) +
  labs(y = expression(bold("Rosette area ("*cm^2*")")),
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


##################.           TUKEY - Fig2C/2D        ##################
ht.aov <- aov(rosette ~ group*geno, data = pheno) #pick
summary(ht.aov)
shapiro.test(resid(ht.aov))
plot(ht.aov,1)
plot(ht.aov,2)
sr <- simulateResiduals(fittedModel = ht.aov, plot = T) #looks okay.

emm <- emmeans(ht.aov, ~group*geno, adjust = "mvt")
cld_result <- cld(emm, Letters = letters)
cld_result

#####################Fig2E#####################
pheno <- read_excel("~/../Touch/Remake/Fig2_data.xlsx", sheet=2)
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))

model_conc <- lm(height ~ geno * group, data = pheno)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ group | geno)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, sep = "_"))

free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = height, fill=Treatment)) +
  facet_wrap(. ~ geno, labeller = as_labeller(
    c(
      "Col-0"   = "Col-0",
      "rrtf1-2" = "italic('rrtf1-2')"
    ), label_parsed)) + coord_cartesian(ylim = c(0, 100)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p= ", signif(p.value, 2))), size = 3.5, x = 1.5, y = 100, vjust = 1.75, inherit.aes = FALSE) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  scale_colour_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "rrtf1-2_UT" = "#B2DF8A",
      "rrtf1-2_T"  = "#33A02C"
    ),
    breaks = c("Col-0_UT", "Col-0_T", "rrtf1-2_UT", "rrtf1-2_T")
  ) +
  labs(y = "Shoot height (mm)",
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
#W400/H350
ggsave("../free_scale.svg",plot = free_scale, height = 5, width = 5.7, bg = "white")

#######Fig2F_% relative height reduction#########################
Relative_height

library(tidyverse)
library(RColorBrewer)
library(ggbeeswarm)
library(patchwork)

#####Fig2F
per_height
sample_sizes <- per_height %>%
  group_by(set, geno) %>%
  summarise(n = n(),  .groups = "drop") %>%
  group_by(set) %>%
  summarise(n = nth(n,1), .groups = "drop")
per_height <- per_height %>%
  left_join(sample_sizes, by = "set") %>%
  mutate(set_label = paste0(set, "\n(n=", n, ")"))
# Set order
per_height$set_label <- factor(per_height$set_label,
                               levels = unique(per_height$set_label))
p_values <- per_height %>%
  group_by(set_label) %>%
  summarise(
    p_value = t.test(raw ~ geno)$p.value,
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
    y_position = max(per_height$percent) * 1.1
  )

print("P-values calculated:")
print(p_values)

means <- per_height %>%
  group_by(set_label, geno) %>%
  summarise(mean_value = mean(percent, na.rm = TRUE),
            .groups = "drop")
# ============================================
# CREATE THE PLOT
# ============================================

plot <- ggplot(per_height, aes(x = set_label, y = percent, fill = geno, color=geno)) +
  
  # 1. BOXPLOT
  geom_boxplot(outlier.shape = NA,  # Don't show outliers (will use jitter)
               width = 0.6,
               color = "black",
               size = 0.8,
               alpha = 0.9) +
  
  # 2. DATA POINTS (DOTS)
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.6),
              size = 2,           # Adjust dot size
              alpha = 0.6, color="black") +
  
  # 4. SIGNIFICANCE BRACKETS
  geom_segment(data = p_values,
               aes(x = as.numeric(set_label) - 0.2,
                   xend = as.numeric(set_label) + 0.2,
                   y = y_position,
                   yend = y_position),
               inherit.aes = FALSE,
               size = 0.8,
               color = "black") +
  
  # 5. P-VALUE LABELS
  geom_text(data = p_values,
            aes(x = set_label,
                y = y_position + max(per_height$percent) * 0.03,
                label = p_label, vjust = -0.05),
            inherit.aes = FALSE,
            size = 4,
            fontface = "bold") +
  
  # 6. COLORS (ADJUST TO YOUR PREFERENCE)
  scale_fill_manual(values = c(
    "Col-0" = "#6A3D9A",      
    "rrtf1-1" = "#33A02C"     
  )) +
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",      
    "rrtf1-1" = "#33A02C"     
  )) +
  
  # 7. Y-AXIS SCALE
  scale_y_continuous(
    limits = c(0, max(per_height$percent) * 1.2),
    breaks = seq(0, ceiling(max(per_height$percent) * 1.2), by = 25)
  ) +
  
  # 8. THEME
  theme_classic() +
  theme(
    # Axes
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    # Panel
    panel.background = element_rect(fill = "white")
  ) +
  
  # 9. LABELS
  labs(x = "Experiment", y = "Touched plant height (% of control)")
# Display the plot
print(plot)

# 3. mean
geom_text(data = means,
          aes(x = set_label, y = mean_value, 
              label = sprintf("%.1f", mean_value)),
          position = position_dodge(width = 0.6))       
# ============================================
# SAVE THE PLOT
# ============================================

ggsave("../final_boxplot.svg",
       plot,
       width = 3.9,
       height = 6.45,
       dpi = 300,
       bg = "white")




########### ============================================ ###########
########### gFig3B JA hormone levels              ###########
########### ============================================ ###########                            
# ============================================
# Add Tukey Letters to Line Plot
# ============================================

library(ggplot2)
library(dplyr)
library(multcomp)
library(multcompView)
library(emmeans)   # For pairwise comparisons
library(car)
library(DHARMa)

library(ggplot2)
horm <- read.csv("~/../thigmo_horm_tukey.csv",sep=",")
horm <- horm %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno))

summary_data <- horm %>%
  group_by(geno, time) %>%
  summarise(
    mean = mean(JA, na.rm = TRUE),           
    se = sd(JA, na.rm = TRUE) / sqrt(n()),   
    .groups = "drop")

print(summary_data)
summary_data$geno <- factor(summary_data$geno, 
                            levels = c("Col-0", "rrtf1"))
horm$time <- factor(horm$time, 
                    levels = c("0","2.5","5","10"))
# Tukey letters
horm.aov <- aov(JA ~ geno*time, data = horm) #pick
summary(horm.aov)
shapiro.test(resid(horm.aov))
leveneTest(JA ~ gbyt, data = horm)
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
             aes(x = time, y = JA, color = geno),
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
  scale_y_continuous(limits = c(0, max(horm$JA, na.rm = TRUE)*1.2)) +
  
  #Label
  labs(x = "Time (min)",
       y = "pmoles/gFW",
       color = "",
       title = "Jasmonic acid (JA)") + 
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


#########################################################################################
########################Fig4ABC################################################
# ============================================
# Add Tukey Letters to Line Plot
# ============================================
library(ggplot2)
library(dplyr)
library(multcomp)  # For Tukey HSD test
library(emmeans)   # For pairwise comparisons
library(car)
library(DHARMa)

expr <- read.csv(file="~/../thigmo_qpcr_tukey.csv", sep=",")

expr$time <- as.factor(expr$time)
gene.aov <- aov(CT ~ geno*time, data = expr[expr$gene == "RRTF1", ]) #pick
summary(gene.aov)
shapiro.test(resid(gene.aov))
leveneTest(expr$CT[expr$gene == "RRTF1"] ~ gbyt, data = expr[expr$gene == "RRTF1", ])

plot(gene.aov,1)
plot(gene.aov,2)
sr <- simulateResiduals(fittedModel = gene.aov, plot = T) #looks okay.

emm <- emmeans(gene.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))

#read excel_relatvie expression file for plot
expr1 <- read_excel("~/../Touch/Remake/Fig4_data.xlsx", sheet=1)
expr1$relative <- as.numeric(as.character(expr1$relative))
expr1$se <- as.numeric(as.character(expr1$se))

summary_data <- expr1[expr1$gene == "RRTF1", 2:5]

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno)) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))
print(summary_data)

summary_data$time <- as.numeric(as.character(summary_data$time))
colnames(summary_data)[3] <- c("mean")
# ==========================
#           Figure
# ==========================
summary_data$time <- as.numeric(as.character(summary_data$time))
horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(size = 5) + geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean + se, mean - se), vjust = ifelse(geno == "Col-0", -.8, 1)), size = 5, fontface = "bold", color = "black")  +
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",   
    "rrtf1" = "#33A02C"     
  )) +
  
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 2000)) +
  
  #Label
  labs(x = "Time (mins)",
       y = "Relative expression",
       color = "",
       title = expression(italic("RRTF1"))) + 
  theme(
    # Axes
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 20, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1))

print(horm_plot)

library(patchwork)
horm_plot1 <- ggplot(summary_data, aes(x = time, y = mean, 
                                       color = geno,
                                       group = geno)) +
  # Line (both line)
  geom_line(size = 1.2) + geom_point(size = 5) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                                                               width = 0.3, size = 0.8) +
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    
    "rrtf1" = "#33A02C"     
  )) +
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(limits = c(0, 2000)) +
  labs(x = "",
       y = "",
       color = "")+
  theme(
    # Axes
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 1, alpha = 0.3,
             position = position_jitter(width = 0.1))
print(horm_plot1)

p_inset <- horm_plot1 +
  coord_cartesian(ylim = c(0, 15), xlim = c(0, 3)) +
  theme_classic(base_size = 4) +                 
  theme(
    legend.position = "none",                    
    axis.text  = element_text(size = 12),         
    axis.title = element_text(size = 9),
    axis.line  = element_line(linewidth = 0.8),  
  )

p_final <- horm_plot +
  inset_element(
    p_inset,
    left = 0.05, bottom = 0.3, right = 0.55, top = 1,  
    align_to = "panel"
  )

ggsave("../final_boxplot.svg", p_final, width = 6, height = 4.5, dpi = 300, bg = "white")

########################################################################
#############                         JAZ7                              #######
########################################################################
expr$time <- as.factor(expr$time)
gene.aov <- aov(CT ~ geno*time, data = expr[expr$gene == "JAZ7", ]) #pick

summary(gene.aov)
shapiro.test(resid(gene.aov))
leveneTest(expr$CT[expr$gene == "JAZ7"] ~ gbyt, data = expr[expr$gene == "JAZ7", ])

plot(gene.aov,1)
plot(gene.aov,2)
sr <- simulateResiduals(fittedModel = gene.aov, plot = T) #looks okay.

emm <- emmeans(gene.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))

#read excel_relatvie expression file for plot
expr1 <- read_excel("~/../Touch/Remake/Fig4_data.xlsx", sheet=1)
expr1$relative <- as.numeric(as.character(expr1$relative))
expr1$se <- as.numeric(as.character(expr1$se))

summary_data <- expr1[expr1$gene == "JAZ7", 2:5]

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno)) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))
print(summary_data)

summary_data$time <- as.numeric(as.character(summary_data$time))
colnames(summary_data)[3] <- c("mean")
# ==========================
#           Figure
# ==========================
summary_data$time <- as.numeric(as.character(summary_data$time))
horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(size = 5) + geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean - se, mean + se), vjust = ifelse(geno == "Col-0", 1.5, -0.5)), size = 5, fontface = "bold", color = "black")  +
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    
    "rrtf1" = "#33A02C"     
  )) +
  
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 30)) +
  
  #Label
  labs(x = "Time (mins)",
       y = "Relative expression",
       color = "",
       title = expression(italic("JAZ7"))) + 
  theme(
    # Axes
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 20, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1))

print(horm_plot)

library(patchwork)
horm_plot1 <- ggplot(summary_data, aes(x = time, y = mean, 
                                       color = geno,
                                       group = geno)) +
  # Line (both line)
  geom_line(size = 1.2) + geom_point(size = 5) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                                                               width = 0.3, size = 0.8) +
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    # 
    "rrtf1" = "#33A02C"     # 
  )) +
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(x = "",
       y = "",
       color = "")+
  theme(
    # Axes
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 1, alpha = 0.3,
             position = position_jitter(width = 0.1))
print(horm_plot1)

p_inset <- horm_plot1 +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 5)) +
  theme_classic(base_size = 4) +                
  theme(
    legend.position = "none",                    
    axis.text  = element_text(size = 12),       
    axis.title = element_text(size = 9),
    axis.line  = element_line(linewidth = 0.8),  
  )

p_final <- horm_plot +
  inset_element(
    p_inset,
    left = 0.05, bottom = 0.3, right = 0.55, top = 1,  # location (0~1, based on panel)
    align_to = "panel"
  )
print(p_final)
ggsave("../final_boxplot.svg", p_final, width = 6, height = 4.5, dpi = 300, bg = "white")



########################################################################
#############                         ERF19                     #######
########################################################################
expr$time <- as.factor(expr$time)
gene.aov <- aov(CT ~ geno*time, data = expr[expr$gene == "ERF19", ]) 

summary(gene.aov)
shapiro.test(resid(gene.aov))
leveneTest(expr$CT[expr$gene == "ERF19"] ~ gbyt, data = expr[expr$gene == "ERF19", ])

plot(gene.aov,1)
plot(gene.aov,2)
sr <- simulateResiduals(fittedModel = gene.aov, plot = T) 
emm <- emmeans(gene.aov, ~geno*time, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)
tukey_all$time <- as.numeric(as.character(tukey_all$time))

#read excel_relatvie expression file for plot
expr1 <- read_excel("~/../Touch/Remake/Fig4_data.xlsx", sheet=1)
expr1$relative <- as.numeric(as.character(expr1$relative))
expr1$se <- as.numeric(as.character(expr1$se))

summary_data <- expr1[expr1$gene == "ERF19", 2:5]

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno)) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))
print(summary_data)

summary_data$time <- as.numeric(as.character(summary_data$time))
colnames(summary_data)[3] <- c("mean")
# ==========================
#           Figure
# ==========================
summary_data$time <- as.numeric(as.character(summary_data$time))
horm_plot <- ggplot(summary_data, aes(x = time, y = mean, 
                                      color = geno,
                                      group = geno)) +
  
  # Line (both line)
  geom_line(size = 1.2) +
  
  # DOT
  geom_point(size = 5) + geom_text(aes(label = letter, y = ifelse(geno == "Col-0", mean - se, mean + se), vjust = ifelse(geno == "Col-0", 1.5, -0.5)), size = 5, fontface = "bold", color = "black")  +
  
  # error bars
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.3, size = 0.8) +
  
  # color
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",    #
    "rrtf1" = "#33A02C"     #
  )) +
  
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(limits = c(0, 1000)) +
  
  #Label
  labs(x = "Time (mins)",
       y = "Relative expression",
       color = "",
       title = expression(italic("ERF19"))) + 
  theme(
    # Axes
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 20, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    plot.title = element_text(size = 20, hjust = 0.5),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 2, alpha = 0.3,
             position = position_jitter(width = 0.1))

print(horm_plot)

library(patchwork)
horm_plot1 <- ggplot(summary_data, aes(x = time, y = mean, 
                                       color = geno,
                                       group = geno)) +
  # Line (both line)
  geom_line(size = 1.2) + geom_point(size = 5) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                                                               width = 0.3, size = 0.8) +
  scale_color_manual(values = c(
    "Col-0" = "#6A3D9A",   
    "rrtf1" = "#33A02C"    
  )) +
  #axis
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(x = "",
       y = "",
       color = "")+
  theme(
    # Axes
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Panel
    panel.background = element_rect(fill = "white")) +
  geom_point(data = summary_data,
             aes(x = time, y = mean, color = geno),
             size = 1, alpha = 0.3,
             position = position_jitter(width = 0.1))
print(horm_plot1)

p_inset <- horm_plot1 +
  coord_cartesian(ylim = c(0, 30), xlim = c(0, 5)) +
  theme_classic(base_size = 4) +                
  theme(
    legend.position = "none",                    
    axis.text  = element_text(size = 12),      
    axis.title = element_text(size = 9),
    axis.line  = element_line(linewidth = 0.8),  
  )

p_final <- horm_plot +
  inset_element(
    p_inset,
    left = 0.05, bottom = 0.3, right = 0.55, top = 1,
    align_to = "panel"
  )
print(p_final)
ggsave("../final_boxplot.svg", p_final, width = 6, height = 4.5, dpi = 300, bg = "white")


#############======================================================##################################
#############=====================    Fig3A   =====================##################################
#############======================================================##################################
gene <- read_excel("~/../Touch/Remake/Fig3_data.xlsx", sheet=1)

####p-value
model_conc <- lm(ct ~ Treatment, data = gene)
anova(model_conc)
est_conc <- emmeans(model_conc, pairwise ~ Treatment)

conc_contrast <- est_conc$contrasts %>% 
  as.data.frame()
conc_contrast


gene <- gene[c(1,4),c(1,3,4)]
gene$Treatment <- factor(gene$Treatment, levels = c("Mock (H2O)", "MeJA 10 mM"))
free_scale <- gene %>% 
  ggplot(aes(x = Treatment, y = mean, fill=Treatment)) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.8, color = c(
    "Mock (H2O)"   = "#CAB2D6",
    "MeJA 10 mM"    = "#FB9A84"
  ),alpha=0.7)+
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  geom_text(data = conc_contrast, aes(label = paste0(
    "p = ", signif(p.value, 2))), x = 1.5, y = 50, vjust = -7, inherit.aes = FALSE, size = 5, fontface = "bold", color = "black")+
  scale_fill_manual(
    values = c(
      "Mock (H2O)"   = "#CAB2D6",
      "MeJA 10 mM"    = "#FB9A84"
    ),
    breaks = c("Mock (H2O)", "MeJA 10 mM")
  ) +
  scale_color_manual(
    values = c(
      "Mock (H2O)"   = "#CAB2D6",
      "MeJA 10 mM"    = "#FB9A84"
    ),
    breaks = c("Mock (H2O)", "MeJA 10 mM")
  )+
  labs(y = "Relative expression",
       x = "Col-0",
       title = expression(paste(italic("RRTF1"), " expression"))) +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.caption= element_text(size = 12, hjust = 0),
    legend.text = element_text(size = 15),
  )


########################################################################
#############               Fig3E                     #######
########################################################################
pheno <- read_excel("~/../Touch/Remake/Fig3_data.xlsx", sheet=2)
pheno <- pheno %>%
  mutate(group = factor(group, levels = c("UT", "T")))
library("glmmTMB")

set.glmm <- glmmTMB(height ~ geno * group * treat, 
                    family = tweedie(link = "log"), data = pheno)
shapiro.test(resid(set.glmm))

plot(set.glmm,1)
plot(set.glmm,2)
sr <- simulateResiduals(fittedModel = set.glmm, plot = T) #looks okay.

emm <- emmeans(set.glmm, ~geno*group*treat, adjust = "mvt")
cld_result <- cld(emm, Letters = letters, adjust = "mvt")
tukey_all <- as.data.frame(cld_result)
tukey_all$.group <- trimws(tukey_all$.group)

#read excel_relatvie expression file for plot
expr1 <- read_excel("~/../Touch/Remake/Fig4_data.xlsx", sheet=1)

summary_data <- expr1[expr1$gene == "RRTF1", 2:5]

tukey_all$time <- as.factor(tukey_all$time)
summary_data$time <- as.factor(summary_data$time)
summary_data <- summary_data %>%
  left_join(
    tukey_all %>% dplyr::select(geno, time, .group), by = c("geno", "time")) %>%
  rename(letter = `.group`) %>%
  mutate(geno = if_else(geno == "col", "Col-0", geno)) %>%
  mutate(geno = factor(geno, levels = c("Col-0", "rrtf1")))
print(summary_data)

summary_data$time <- as.numeric(as.character(summary_data$time))
colnames(summary_data)[3] <- c("mean")

library(patchwork) 

phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "Col-0", ]

phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

##group+treat
free_scale <- phenofill %>% 
  ggplot(aes(x = group, y = height, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("Col-0" = "Col-0")
    ,label_parsed)) + coord_cartesian(ylim = c(0, 100)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
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
  labs(y = "Shoot height (mm)",
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



##############rrtf1######################################################################
phenofill <- pheno %>%
  mutate(Treatment = paste(geno, group, treat, sep = "_"))
phenofill<- phenofill[phenofill$geno == "rrtf1", ]
phenofill$treat <- factor(phenofill$treat, levels = c("Mock", "MeJA"))

free_scale1 <- phenofill %>% 
  ggplot(aes(x = group, y = height, fill=Treatment)) +
  facet_grid(geno ~ treat, labeller = labeller(
    Treatment = c("Mock" = "Mock", "MeJA" = "MeJA"),
    geno = c("rrtf1" = expression(italic("RRTF1")))
    ,label_parsed)) + coord_cartesian(ylim = c(0, 100)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.5, linewidth = 0.8) +
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
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
ggsave("../free_scale.svg",plot = combined_plot, width = 5.6, height = 3.5, bg = "white")

####################################################################################
##########################################Fig8F / Fig3F####################################
per_height <- read_excel("~/../Touch/Remake/Fig3_data.xlsx", sheet=3)

per_height <- per_height %>%
  mutate(group = paste(geno, treat, sep = "_"))

# Set order
per_height$group <- factor(per_height$group, levels = c("Col-0_Mock", "Col-0_MeJA", "rrtf1_Mock","rrtf1_MeJA"))

p_values <- per_height %>%
  group_by(geno) %>%
  summarise(
    p_value = t.test(raw ~ treat)$p.value,
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
    y_position = max(per_height$percent) * 1.1
  )

print("P-values calculated:")
print(p_values)

# ============================================
# CREATE THE PLOT
# ============================================
plot <- ggplot(per_height, aes(x = geno, y = percent, fill = group, color = group)) +
  
  # 1. BOXPLOT
  geom_boxplot(outlier.shape = NA,  # Don't show outliers (will use jitter)
               width = 0.6,
               color = "black",
               size = 0.8,
               alpha = 0.9) +
  
  # 2. DATA POINTS (DOTS)
  geom_jitter(aes(color=group),
              position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.6),
              size = 2,          
              alpha = 0.5) +
  
  # 4. SIGNIFICANCE BRACKETS
  geom_segment(data = p_values,
               aes(x = as.numeric(as.factor(geno)) - 0.2,
                   xend = as.numeric(as.factor(geno)) + 0.2,
                   y = y_position,
                   yend = y_position),
               inherit.aes = FALSE,
               size = 0.8,
               color = "black") +
  
  # 5. P-VALUE LABELS
  geom_text(data = p_values,
            aes(x = geno,
                y = y_position + max(per_height$percent) * 0.03,
                label = p_label, vjust = -0.05),
            inherit.aes = FALSE,
            size = 5,
            fontface = "bold") +
  
  # 6. COLORS
  scale_fill_manual(values = c(
    "Col-0_Mock" = "#6A3D9A",      
    "Col-0_MeJA" = "#E31A1C",
    "rrtf1_Mock" = "#33A02C",
    "rrtf1_MeJA" = "#FF7F00"
  )) +
  scale_color_manual(values = c(
    "Col-0_Mock" = "#3D1F6B",      
    "Col-0_MeJA" = "#8B0000",
    "rrtf1_Mock" = "#1A5C16",
    "rrtf1_MeJA" = "#CC5500"
  )) +
  
  # 7. Y-AXIS SCALE
  scale_y_continuous(
    limits = c(0, max(per_height$percent) * 1.2),
    breaks = seq(0, ceiling(max(per_height$percent) * 1.2), by = 25)
  ) +
  
  # 8. THEME
  theme_classic() +
  theme(
    # Axes
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    # Legend
    legend.position = "top",
    text = element_text(size = 16, color = "black"),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    plot.caption= element_text(size = 12, hjust = 0),
    # Panel
    panel.background = element_rect(fill = "white")
  ) +
  
  # 9. LABELS
  labs(x = "Genotype", y = "Touched plant height (% of control)")
# Display the plot
print(plot)

########################################################################
#############                       Fig3G                     #######
########################################################################

gene <- read_excel("~/../Touch/Remake/Fig3_data.xlsx", sheet=4)

####p-value
est_conc <- pairwise.t.test(gene$ct, gene$Treatment, 
                            p.adjust.method = "bonferroni")

p_val <- est_conc$p.value[1, 1]
conc_contrast <- data.frame(
  p.value = p_val
)

gene <- gene[c(1,2,5,6),c(1,2,4,5)]


genefill <- gene %>%
  mutate(Treatment = paste(Treatment, Time, sep = "_"))
genefill$Treatment <- factor(genefill$Treatment, levels = c("Col-0_UT","aos_UT","Col-0_T","aos_T"))
genefill$Time <- factor(genefill$Time, levels = c("UT","T"))

genefill <- genefill %>%
  mutate(
    geno = ifelse(grepl("Col-0", Treatment), "Col-0", "aos")
  )
genefill$geno <- factor(genefill$geno, levels = c("Col-0", "aos"))


free_scale <- genefill %>% 
  ggplot(aes(x = geno, y = mean, fill=Treatment)) +
  coord_cartesian(ylim = c(0, 3200)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.7, linewidth = 0.8) + 
  facet_wrap(~Time, labeller = as_labeller(c("UT" = "Untouched","T" = "Touched")), scales = "free_x") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.8, color = c(
    "Col-0_UT"   = "#CAB2D6",
    "Col-0_T"    = "#6A3D9A",
    "aos_UT"    = "#A6CEE3",
    "aos_T"    = "#1F78B4"
  ),alpha=0.7)+
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "aos_UT"    = "#A6CEE3",
      "aos_T"    = "#1F78B4"
    )
  ) +
  scale_color_manual(
    values = c(
      "Col-0_UT"   = "#CAB2D6",
      "Col-0_T"    = "#6A3D9A",
      "aos_UT"    = "#A6CEE3",
      "aos_T"    = "#1F78B4"
    )
  )+
  labs(y = "Relative expression",
       x = "",
       title = expression(paste(italic("RRTF1"), " expression"))) +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    panel.spacing = unit(0.5, "lines"),
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.caption= element_text(size = 12, hjust = 0),
    legend.text = element_text(size = 15),
  )

free_scale
ggsave("../free_scale.svg",plot = free_scale, width = 4.8, height = 6, bg = "white") 
ggsave("../free_scale.svg",plot = free_scale, width = 3, height = 6, bg = "white") 

########################################################################
#############                       Fig3H                     #######
########################################################################

gene <- read_excel("~/../Touch/Remake/Fig3_data.xlsx", sheet=5)

####p-value
est_conc <- pairwise.t.test(gene$ct, gene$Treatment, 
                            p.adjust.method = "bonferroni")

p_val <- est_conc$p.value[1, 1]
conc_contrast <- data.frame(
  p.value = p_val
)

gene <- gene[c(1,2,5,6),c(1,2,4,5)]


genefill <- gene %>%
  mutate(Treatment = paste(Treatment, Time, sep = "_"))
genefill$Treatment <- factor(genefill$Treatment, levels = c("Mock (H2O)_UT","Jarin-1 50 µM_UT","Mock (H2O)_T","Jarin-1 50 µM_T"))
genefill$Time <- factor(genefill$Time, levels = c("UT","T"))

genefill <- genefill %>%
  mutate(
    geno = gene$Treatment)
genefill$geno <- factor(genefill$geno, levels = c("Mock (H2O)", "Jarin-1 50 µM"))


free_scale <- genefill %>% 
  ggplot(aes(x = geno, y = mean, fill=Treatment)) +
  coord_cartesian(ylim = c(0, 2200)) +
  geom_bar(stat = "summary", fill = NA, aes(color = Treatment), 
           alpha = 0.8, fun = mean, width = 0.7, linewidth = 0.8) + 
  facet_wrap(~Time, labeller = as_labeller(c("UT" = "Untouched","T" = "Touched")), scales = "free_x") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, linewidth = 0.8, color = c(
    "Mock (H2O)_UT"   = "#CAB2D6",
    "Mock (H2O)_T"    = "#6A3D9A",
    "Jarin-1 50 µM_UT"    = "#EEC229FF",
    "Jarin-1 50 µM_T"    = "#B15928"
  ),alpha=0.7)+
  ggbeeswarm::geom_quasirandom(
    shape = 21, size = 2.5, alpha = 0.8, color = "white",
    aes(fill = Treatment), width = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Mock (H2O)_UT"   = "#CAB2D6",
      "Mock (H2O)_T"    = "#6A3D9A",
      "Jarin-1 50 µM_UT"    = "#EEC229FF",
      "Jarin-1 50 µM_T"    = "#B15928"
    )
  ) +
  scale_color_manual(
    values = c(
      "Mock (H2O)_UT"   = "#CAB2D6",
      "Mock (H2O)_T"    = "#6A3D9A",
      "Jarin-1 50 µM_UT"    = "#EEC229FF",
      "Jarin-1 50 µM_T"    = "#B15928"
    )
  )+
  labs(y = "Relative expression",
       x = "",
       title = expression(paste(italic("RRTF1"), " expression"))) +
  theme_classic() +
  theme(
    legend.position = "top",
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    panel.spacing = unit(0.5, "lines"),
    title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.caption= element_text(size = 12, hjust = 0),
    legend.text = element_text(size = 15),
  )

free_scale

#####################################################Fig 4F###############################################################
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
core_wt10_rrtf10_rrtf0 <- read.csv(file = "~/../core_touch_compare/core_wt10_rrtf10_rrtf0.csv", sep=",") #core_wt10_rrtf10_rrtf0 = 679 genes
core_wt10_rrtf10_676 <- read.csv(file = "~/../Touch/Remake/supp/core_wt10_rrtf10_676.csv", sep=",")

rownames(core_wt10_rrtf10_676) <- core_wt10_rrtf10_676[,1]
core_wt10_rrtf10_676 <- core_wt10_rrtf10_676[,2:3]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

myCol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(9)

#Group order
g1 <- rownames(heat_capped)[heat_capped[, "rrtf1_T"] == 0 & heat_capped[, "Col-0_T"] != 0]  
g2 <- rownames(heat_capped)[heat_capped[, "Col-0_T"] == 0 & heat_capped[, "rrtf1_T"] != 0]  
g3 <- setdiff(rownames(heat_capped), c(g1, g2))                                               

# Reorder
g1 <- g1[order(heat_capped[g1, "Col-0_T"], decreasing = TRUE)]
g2 <- g2[order(heat_capped[g2, "rrtf1_T"], decreasing = TRUE)]
g3 <- g3[order(heat_capped[g3, "Col-0_T"], decreasing = TRUE)]

gene_order <- c(g1, g3, g2)

# Group label
group_label <- setNames(
  c(rep("Col-0 T-specific", length(g1)),
    rep("Both", length(g3)),
    rep("rrtf1 T-specific", length(g2))),
  gene_order
)

free_scale <- reshape2::melt(heat_capped) %>%
  setNames(c("gene", "sample", "value")) %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    sample = factor(sample, levels = c("Col-0_T", "rrtf1_T")),
    group = group_label[as.character(gene)],
    value2 = case_when(value >= 10 ~ 10, TRUE ~ value)
  ) %>%
  ggplot(aes(x = sample, y = gene)) +
  geom_tile(aes(fill = value2)) +
  scale_fill_gradientn(
    colors = myCol,
    breaks = c(0, 4, 8),
    labels = c("<0", "4", ">8")
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Touch", y = "Genes", fill = "Z-score",
       title = "676 DEGs: Col-0 T vs rrtf1 T") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0)
  )