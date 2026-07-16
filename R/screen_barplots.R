# Figure 1d,e and Supplemental Figure 1
# By Nolan Origer

# Load data
counts <- read.csv("~/supp_table_1.csv")

counts <- counts[order(counts$Count_inhib, decreasing = T),]
counts$Gene_Symbol <- factor(counts$Gene_Symbol, levels = unique(counts$Gene_Symbol))

adh.counts <- filter(counts, Gene_Symbol %in% adhesome$gene)
ggplot(adh.counts, aes(x = Viral_titer)) +
  geom_histogram()
mean(adh.counts$Viral_titer) + 3*sd(adh.counts$Viral_titer)
mean(adh.counts$Viral_titer) - 3*sd(adh.counts$Viral_titer)
adh.counts[which(adh.counts$Viral_titer < 1e+07),]
adh.counts <- adh.counts[-which(adh.counts$Viral_titer < 1e+07),]


control.counts <- filter(counts, Condition == "Control")
ggplot(control.counts, aes(x = Viral_titer)) +
  geom_histogram()
mean(control.counts$Viral_titer) + 3*sd(control.counts$Viral_titer)
mean(control.counts$Viral_titer) - 3*sd(control.counts$Viral_titer)
control.counts[which(control.counts$Viral_titer > 0 & control.counts$Viral_titer < 1e+07),]
control.counts <- control.counts[-which(control.counts$Viral_titer > 0 & control.counts$Viral_titer < 1e+07),]

neg.ctrl.counts <- filter(control.counts, Gene_Symbol %in% c("NEG.CTRL","EMPTY"))

control.counts <- filter(control.counts, Gene_Symbol %in% c("RFP", "LUCIFERASE", "GFP", "lacZ"))
control.counts$Gene_Symbol <- factor(control.counts$Gene_Symbol, levels = c("RFP", "LUCIFERASE", "GFP", "lacZ"))


# Calculate mean counts
{
  mean.counts <- rbind(adh.counts, control.counts) %>% 
    group_by(Gene_Symbol) %>% 
    summarise("Count", Mean = mean(Count))
  colnames(mean.counts)[2] <- "Stat"
  mean.counts[,2] <- rep("Count", length(mean.counts[,2]))
  mean.counts$log2fc <- log2(mean.counts$Mean/mean(control.counts$Count))
  mean.counts <- mean.counts[order(mean.counts$log2fc, decreasing = T),]
  
  mean.counts.inhib <- rbind(adh.counts, control.counts) %>% 
    group_by(Gene_Symbol) %>% 
    summarise("Count_inhib", Mean = mean(Count_inhib))
  colnames(mean.counts.inhib)[2] <- "Stat"
  mean.counts.inhib[,2] <- rep("Count", length(mean.counts.inhib[,2]))
  mean.counts.inhib$log2fc_inhib <- log2(mean.counts.inhib$Mean/mean(control.counts$Count_inhib))
  mean.counts.inhib <- mean.counts.inhib[order(mean.counts.inhib$log2fc_inhib, decreasing = T),]
}

# Calculate mean percent area
{
  mean.area <- rbind(adh.counts, control.counts) %>% 
    group_by(Gene_Symbol) %>% 
    summarise("pArea", Mean = mean(pArea))
  colnames(mean.area)[2] <- "Stat"
  mean.area[,2] <- rep("Area", length(mean.area[,2]))
  mean.area$log2fc <- log2(mean.area$Mean/mean(control.counts$pArea))
  mean.area <- mean.area[order(mean.area$log2fc, decreasing = T),]
  
  mean.area.inhib <- rbind(adh.counts, control.counts) %>% 
    group_by(Gene_Symbol) %>% 
    summarise("pArea_inhib", Mean = mean(pArea_inhib))
  colnames(mean.area.inhib)[2] <- "Stat"
  mean.area.inhib[,2] <- rep("Area", length(mean.area.inhib[,2]))
  mean.area.inhib$log2fc_inhib <- log2(mean.area.inhib$Mean/mean(control.counts$pArea_inhib))
  mean.area.inhib <- mean.area.inhib[order(mean.area.inhib$log2fc_inhib, decreasing = T),]
}

# Merge and set factor by count ranking
mean.combined <- rbind(mean.counts, mean.area)
mean.combined.inhib <- rbind(mean.counts.inhib, mean.area.inhib)

mean.combined.both <- merge(mean.combined, mean.combined.inhib, by = c("Gene_Symbol", "Stat")) %>%
  dplyr::filter(Stat == "Count") %>%
  dplyr::select("Gene_Symbol", "log2fc", "log2fc_inhib")
mean.combined.both$Label <- ifelse((mean.combined.both$log2fc > 2.77 | 
                                      mean.combined.both$log2fc_inhib > 1.7), 
                                   yes = "Top", 
                                   no = ifelse(mean.combined.both$Gene_Symbol %in% control.counts$Gene_Symbol, 
                                               yes = "Control", 
                                               no = "None"))

mean.combined <- filter(mean.combined, Gene_Symbol %in% adh.counts$Gene_Symbol)
mean.combined$Gene_Symbol <- factor(mean.combined$Gene_Symbol, levels = filter(mean.counts, Gene_Symbol %in% adh.counts$Gene_Symbol)$Gene_Symbol)
mean.combined.inhib <- filter(mean.combined.inhib, Gene_Symbol %in% adh.counts$Gene_Symbol)
mean.combined.inhib$Gene_Symbol <- factor(mean.combined.inhib$Gene_Symbol, levels = filter(mean.counts.inhib, Gene_Symbol %in% adh.counts$Gene_Symbol)$Gene_Symbol)

# Plot for Figure 1D (counts only, with inhibitors)
ggplot(filter(mean.combined.inhib, Stat == "Count"), aes(x = Gene_Symbol, y = log2fc_inhib, fill = Stat)) + 
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_hline(yintercept = 0, color = "black") + 
  theme_classic() +
  theme(axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(color="black", size = 10, 
                                   angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") + 
  labs(x = NULL, y = "Log Fold Change Colony Count", title = "Figure 1D") +
  scale_fill_manual(values = c("Count" = "#4e69b0"))

library(ggrepel)
ggplot(mean.combined.both, aes(x = log2fc, y = log2fc_inhib, color = Label)) + 
  geom_point() +
  geom_text_repel(data = filter(mean.combined.both, Label %in% c("Top", "Control")),
                  aes(x = log2fc, y = log2fc_inhib, label = Gene_Symbol),
                  max.overlaps = Inf) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") + 
  theme_classic() +
  theme(axis.text = element_text(color="black", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none") + 
  labs(x = "Log Fold Change (no inhibitors)", 
       y = "Log Fold Change (LSDi/ROCKi)", 
       title = "Figure 1E") +
  scale_color_manual(values = c("Top" = "red",
                                "None" = "black",
                                "Control" = "#979797")) + 
  geom_smooth(mapping = aes(x = log2fc, y = log2fc_inhib, color = "black"),
              method = 'lm', formula = y~x)

summary(rr <- rlm(log2fc_inhib ~ log2fc, data = mean.combined.both))
summary(ols <- lm(log2fc_inhib ~ log2fc, data = mean.combined.both))

# Plot for Supp Figure 1A/B
ggplot(mean.combined, aes(x = Gene_Symbol, y = log2fc, fill = Stat)) + 
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_hline(yintercept = 0, color = "black") + 
  theme_classic() +
  theme(axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(color="black", size = 10, 
                                   angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = c(0.95,0.75)) + 
  labs(x = NULL, y = "Log (Fold Change)", 
       title = "Dynamic Adhesome Reprogramming (no inhibitors)") +
  scale_fill_manual(values = c("Count" = "#4e69b0",
                               "Area" = "#bebebe"),
                    name = NULL)

ggplot(mean.combined.inhib, aes(x = Gene_Symbol, y = log2fc_inhib, fill = Stat)) + 
  geom_bar(stat = "identity", position = "dodge", color="black") +
  geom_hline(yintercept = 0, color = "black") + 
  theme_classic() +
  theme(axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(color="black", size = 10, 
                                   angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = c(0.95,0.75)) + 
  labs(x = NULL, y = "Log (Fold Change)", 
       title = "Dynamic Adhesome Reprogramming (LSDi / ROCKi)") +
  scale_fill_manual(values = c("Count" = "#4e69b0",
                               "Area" = "#bebebe"),
                    name = NULL)


library(ggbreak)

ggplot(counts, aes(x = Gene_Symbol, y = Count_inhib, fill = Gene_Symbol)) + 
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(color="black"), 
        axis.text.x = element_text(angle = 90),
        legend.position = "none") + 
  labs(x = NULL, y = "Count", title = "Control Colony Count Comparison")

mean(neg.ctrl.counts$Count_inhib)

s3.counts <- filter(counts, Gene_Symbol == "SHROOM3")
s3.counts$Gene_Symbol <- c("SHROOM3 #1", "SHROOM3 #2", "SHROOM3 #3")
s3.counts$Gene_Symbol <- factor(s3.counts$Gene_Symbol, levels = c("SHROOM3 #3", "SHROOM3 #2", "SHROOM3 #1"))


combined.counts <- rbind(control.counts, s3.counts)

mean.counts <- combined.counts %>% 
  group_by(Gene_Symbol) %>% 
  summarise("Count_inhib",mean=mean(Count_inhib), sd=sd(Count_inhib))

# Plot for Supp Figure 1C/D
ggplot(combined.counts, aes(x = Gene_Symbol, y = mean, fill = Gene_Symbol)) + 
  geom_bar(data = mean.counts, stat = "identity") +
  geom_dotplot(aes(x = combined.counts$Gene_Symbol, y = combined.counts$Count_inhib),
               stackdir = "center", 
               binaxis = "y", 
               binpositions = "all",
               stackratio = 2,
               dotsize = 0.25,
               fill = "black") + 
  geom_errorbar(data= mean.counts,
                aes(ymin=mean-sd, ymax=mean+sd), width=.1,
                position=position_dodge(.9)) +
  geom_hline(yintercept = mean(neg.ctrl.counts$Count_inhib),
             linetype = "dashed") + 
  geom_hline(yintercept = mean(control.counts$Count_inhib),
             linetype = "dashed", color = "red") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1200, by = 100)) +
  scale_y_break(c(300, 600), scales = 0.5) + 
  theme(axis.text = element_text(color="black", size = 14),
        axis.title=element_text(size=14, face="bold"),
        legend.position = "none") + 
  labs(x = NULL, y = "Colony Counts (w/ inhibs)", title = NULL) + 
  scale_fill_viridis_d()

mean.counts <- combined.counts %>% 
  group_by(Gene_Symbol) %>% 
  summarise("Count",mean=mean(Count), sd=sd(Count))

ggplot(combined.counts, aes(x = Gene_Symbol, y = mean, fill = Gene_Symbol)) + 
  geom_bar(data = mean.counts, stat = "identity") +
  geom_dotplot(aes(x = combined.counts$Gene_Symbol, y = combined.counts$Count),
               stackdir = "center", 
               binaxis = "y", 
               binpositions = "all",
               stackratio = 2,
               dotsize = 0.4,
               binwidth = 5,
               fill = "black") + 
  geom_errorbar(data= mean.counts,
                aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_hline(yintercept = mean(neg.ctrl.counts$Count),
             linetype = "dashed") + 
  geom_hline(yintercept = mean(control.counts$Count),
             linetype = "dashed", color = "red") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1200, by = 25)) +
  scale_y_break(c(80, 250), scales = 0.5) + 
  theme(axis.text = element_text(color="black", size = 14),
        axis.title=element_text(size=14, face="bold"),
        legend.position = "none") + 
  labs(x = NULL, y = "Colony Counts (no inhibitors)", title = NULL) + 
  scale_fill_viridis_d()
