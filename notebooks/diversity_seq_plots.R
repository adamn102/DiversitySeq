### Diversity Sequencing Plots for Dissertation 

###
# Read in packages 
###

library(tidyverse) 
library(readxl)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(survival)
library(survminer)


###
# Set Directory 
###
setwd("~/Dropbox (Partners HealthCare)/Harvard/DAC/Dissertation/analysis/notebooks")

###
# Set Parameters
###
group.colors <- c(JRCSF = "#7570b3", NL43 = "#1b9e77", REJOc ="#d95f02")
escape_group.colors_annotation <- c(Suppressed = "#FB9A99", Escaped = "#E31A1C")

###
# Functions
###

# Shannon Entropy - Scaled by Window Size
hsn <- function(freq,window=100){
  return(-1*sum(freq*log(freq))/window)
}

# bootstrap wrapper for Shannon Entropy 
hsn_bootstratp <- function(freq,bootstraps=10,sample_size=100){
  out <- c()
  for (i in 1:bootstraps){
    filter_df <- sample(freq,sample_size)
    diversity <- hsn(filter_df,window=2562) / sample_size
    out <- append(out,diversity)
  }
  return(out)
  
}

# Hill diversity function 
hill <- function(freq,q){
  return(sum(freq^q))
}


###
# Load data
###

#sample meta data
group_SNV_df <- read_csv("../data/mouse_exp/grouped_SNV_data.csv")
group_SNV_df <- group_SNV_df %>% select(c(sample,snv_count,mouse_number,ref,filtered_reads,exp,week,mouse,average_entropy))
group_SNV_df$sample = paste(group_SNV_df$exp,group_SNV_df$mouse_number,sep="-")


#variant data 
snv_df <- read_csv("../data/mouse_exp/var_data.csv")
snv_df <- snv_df %>% dplyr::select(c(pos,ref,alt,pos_alt,freq,aa_pos,aa_ref,aa_alt,SynNonSyn,MutType,exp,sample,virus))
snv_df$mouse_number <- snv_df$sample
snv_df$sample = paste(snv_df$exp,snv_df$mouse_number,sep="-")

#viral load data
VL_df <- read_csv("../data/mouse_exp/merged_viral_load_data.csv")
VL_df$sample <- VL_df$id
# compute avearge VL at week 3
vl_df <- VL_df %>% 
  filter(week == 3) %>%
  group_by(sample) %>%
  summarise(vl = mean(viral_load))

# compute average entropy from variant data
computed_diversity <- snv_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(Hsn = hsn(freq,window=2563),
                   hill_1 = hill(freq,q=1),
                   hill_2 = hill(freq,q=2),
                   hill_3 = hill(freq,q=3))

df <- full_join(group_SNV_df,computed_diversity)

# average entropy from Syn variant data
computed_syn_diversity <- snv_df %>%
  dplyr::filter(SynNonSyn == "Syn") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(Syn_Hsn = hsn(freq,window=2563))

df <- full_join(df,computed_syn_diversity)
df <- left_join(df,vl_df, by="sample")


#survival data
survival_df <- read_excel("../data/mouse_exp/survival_data.xlsx",sheet="Sheet1")
colnames(survival_df) <-c("sample","Time","Antibody","Status")

#merge data frames 
group_df <- dplyr::inner_join(survival_df, group_SNV_df, by = "sample") %>% dplyr::distinct()
group_df <-  left_join(group_df,vl_df, by="sample")


div_df <- dplyr::inner_join(df,survival_df, by="sample") %>% dplyr::distinct()

var_df <- dplyr::inner_join(survival_df, snv_df, by = "sample") %>% dplyr::distinct()

### Add average mutation frequency
boot_df <- var_df %>% 
  group_by(sample) %>%
  summarise(avg_freq = mean(freq))

### merge with survival df
div_df <- dplyr::inner_join(boot_df,div_df, by="sample") %>% dplyr::distinct()


###
# Set plotting parameters
####


###
# Figure 2.3.A Plot Average Entropy by Escape Status for REJO.c 
###

status_table <- c("Suppressed", "Escaped")
names(status_table) <- c(0,1)

dff <- group_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=average_entropy)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) +
  theme_prism() + 
  xlab("") + 
  ylab("Average Entropy (Hs)") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_diversity_vs_escape.png",units ="cm",width=15,height=10)


###
# Figure 2.3.B Plot Number of Called SNVs by Escape Status for REJO.c 
###

dff <- group_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=snv_count)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() +
  xlab("") + 
  ylab("Number of Called SNVs") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_SNVs_vs_escape.png",units ="cm",width=15,height=10)


###
# Figure 2.3.C Plot Average Entropy of Syn Mutation by Escape Status for REJO.c 
###

status_table <- c("Suppressed", "Escaped")
names(status_table) <- c(0,1)

dff <- div_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status)) %>% 
  dplyr::filter(Syn_Hsn > 0)
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=Syn_Hsn)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() + 
  scale_y_log10(
    limits = c(5e-6, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(6:3, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Synonymous  Mutation Entropy (Hs)") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_syn_diversity_vs_escape.png",units ="cm",width=15,height=10)



###
# Figure 2.3.D Plot Sequencing Coverage by Escape Status for REJO.c 
###

status_table <- c("Suppressed", "Escaped")
names(status_table) <- c(0,1)

dff <- div_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=filtered_reads)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() + 
  scale_y_log10(
    limits = c(1e3, 1e6), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Coverage (Filtered Reads)") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_seq_coverage_vs_escape.png",units ="cm",width=15,height=10)

###
# Figure 2.3.E Average Frequency of Variants within each sample
###

status_table <- c("Suppressed", "Escaped")
names(status_table) <- c(0,1)

dff <- var_df %>%
  dplyr::filter(virus ==  "REJOc") %>% 
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))

boot_df <- dff %>% 
  group_by(sample) %>%
  summarise(avg_freq = mean(freq))

### merge with survival df
dff <- dplyr::inner_join(boot_df,survival_df, by="sample") %>% dplyr::distinct()

  
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=avg_freq)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() + 
  scale_y_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Mutation Frequency") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_seq_mut_var_freq_vs_escape.png",units ="cm",width=15,height=10)

###
# Figure 2.3.F dN/dS 
###

#### Compute Synonymous and NonSynonymous sites

REJOc_possible_mut_df <- read_csv("../data/mouse_exp/REJOc_possible_mut_df.csv")
REJOc_possible_mut_df$syn <- REJOc_possible_mut_df$ref_AA == REJOc_possible_mut_df$mut_AA

syn_rejoc <- REJOc_possible_mut_df %>% 
  count(syn)
syn_count <- syn_rejoc %>% filter(syn == TRUE) %>% select(n)
syn_count <- syn_count$n
nonsyn_count <- syn_rejoc %>% filter(syn == FALSE) %>% select(n)
nonsyn_count <- nonsyn_count$n

status_table <- c("Suppressed", "Escaped")
names(status_table) <- c(0,1)

dff <- var_df %>%
  dplyr::filter(virus ==  "REJOc") %>% 
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))

syn_df <- dff %>% 
  group_by(sample,SynNonSyn) %>%
  count(SynNonSyn) %>%
  spread(SynNonSyn,n)

syn_df[is.na(syn_df)] <- 0

syn_df$ratio <- (syn_df$NonSyn/nonsyn_count) / (syn_df$Syn/syn_count) 



### merge with survival df
dff <- dplyr::inner_join(syn_df,survival_df, by="sample") %>% dplyr::distinct()


# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=ratio)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() + 
  xlab("") + 
  scale_y_log10(
    limits = c(0.1, 10), 
    minor_breaks = rep(1:9, 4)*(10^rep(-1:1, each = 9)),
    guide = "prism_minor"
  ) +
  ylab("dN/dS") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_seq_dNdS_vs_escape.png",units ="cm",width=15,height=10)


###
# Figure 2.3.G viral load
###

#### Compute Synonymous and NonSynonymous sites

dff <- group_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))
# 
VL_df$sample <- VL_df$id
vl_df <- VL_df %>% filter(week == 3) %>%
  group_by(sample) %>%
  summarise(vl = mean(viral_load))


dff <- left_join(dff,vl_df)

my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = Status ,y=vl)) + 
  geom_boxplot(aes(color=Status)) +
  geom_point(aes(color=Status)) +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  theme_prism() + 
  scale_y_log10(
    limits = c(1e4, 1e7), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Viral Load") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_VL_vs_escape.png",units ="cm",width=15,height=10)







###
# Test if input data correlates with observed outcomes 
###

dff <- div_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status)) 
# 
my_comparisons <- list(c("Suppressed", "Escaped"))
dff$Status <- as.factor(dff$Status)
dff$Status <- status_table[as.character(dff$Status)]

f <- ggplot(dff,aes(x = filtered_reads,y=average_entropy)) + 
  geom_point(aes(color=Status),size=8) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Number of Input Reads") + 
  ylab("Average Entropy (Hs)") +
  #ggtitle("Diversity of REJO.c Escape") + 
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:2, each = 9)),
    guide = "prism_minor"
  ) +
  scale_x_log10(
    limits = c(1e3, 1e6), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_input_reads_vs_diversity.png",units ="cm",width=15,height=10)

cor.test(dff$vl,dff$average_entropy)
cor.test(dff$filtered_reads,dff$average_entropy)
cor.test(dff$snv_count,dff$average_entropy)
cor.test(dff$avg_freq,dff$average_entropy)

cor.test(dff$snv_count,dff$avg_freq)

###
# SNV vs Hsn
###
f <- ggplot(dff,aes(x = snv_count,y=average_entropy)) + 
  geom_point(aes(color=Status),size=8) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Number of Called SNVs") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1, 100), 
    minor_breaks = rep(1:9, 4)*(10^rep(0:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_snpcount_vs_diversity.png",units ="cm",width=15,height=10)

###
# SNV vs Average Freq
###
f <- ggplot(dff,aes(x = snv_count,y=avg_freq)) + 
  geom_point(aes(color=Status),size=6) +
  #geom_smooth(method = "lm",aes(color=Status)) +
  theme_prism() +
  xlab("Number of Called SNVs") + 
  ylab("Mutation Frequency") +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1, 100), 
    minor_breaks = rep(1:9, 4)*(10^rep(0:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_snpcount_vs_freq.png",units ="cm",width=15,height=10)

###
# regression to explain this patterns
###
summary(lm(avg_freq ~ snv_count + log10(filtered_reads) , data = dff))
summary(lm(avg_freq ~ snv_count + log10(filtered_reads) + , data = dff))

###
# VL vs Hsn
###
f <- ggplot(dff,aes(x = vl,y=average_entropy)) + 
  geom_point(aes(color=Status),size=8) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Viral Load") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1e4, 1e7), 
    minor_breaks = rep(1:9, 4)*(10^rep(4:7, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_VL_vs_diversity.png",units ="cm",width=15,height=10)


###
# VL vs Hsn
###
f <- ggplot(dff,aes(x = avg_freq,y=average_entropy)) + 
  geom_point(aes(color=Status),size=8) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Average Mutation Frequency") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_mut_freq_vs_diversity.png",units ="cm",width=15,height=10)


###
# Logistic Regression of Input Variables 
###

dff <- div_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Antibody == 'VRC07') %>% 
  dplyr::filter(!is.na(Status))

logistic_regression <- glm(Status ~ log10(average_entropy) + log10(filtered_reads)  + log10(vl) , data = dff, family = binomial(link = "logit"))
logistic_regression <- glm(Status ~ log10(average_entropy) + log10(filtered_reads)  + log10(vl) , data = dff, family = binomial(link = "logit"))
summary(logistic_regression)
confint(logistic_regression)

### Scale all of the input parameters 
dff$scale_avg_freq <- scale(dff$avg_freq)
dff$scale_average_entropy <- scale(dff$average_entropy)
dff$scale_filtered_reads <- scale(log10(dff$filtered_reads))
dff$scale_snv_count<- scale(dff$snv_count)
dff$scale_vl <- scale(log10(dff$vl))

logistic_regression <- glm(Status ~ scale_average_entropy  + scale_filtered_reads  + scale_vl, data = dff, family = binomial(link = "logit"))
logistic_regression <- glm(Status ~ scale_average_entropy  + scale_filtered_reads  + scale_vl + scale_snv_count +  scale_avg_freq, data = dff, family = binomial(link = "logit"))

summary(logistic_regression)
confint(logistic_regression)

#### Relationship between Hsn, SNV count and average freq
model <- glm(Hsn ~ snv_count + avg_freq +  log10(vl) + log10(filtered_reads), data=dff)
summary(model)

#summary(glm(Status ~ log10(filtered_reads) + snv_count + log10(average_entropy), data=dff,family = "binomial"))
model <- glm(Status ~ log10(filtered_reads) + log10(snv_count) + log10(average_entropy), data=dff,family = "binomial")
summary(model)
confint(model)

###
# Figure 2.7.A - REJO.c KM Survival Plot 
###
dff <- group_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>% 
  dplyr::filter(Antibody %in% c("VRC07","Luciferase")) %>% 
  dplyr::filter(!is.na(Status))

fit <- survfit(Surv(Time, Status) ~ Antibody, data = dff)

survp <- ggsurvplot(
  fit, 
  data=dff, 
  conf.int = FALSE,
  pval = FALSE,
  palette = c("black",'#d95f02'),
  ggtheme = theme_bw(),
  ncensor.plot = TRUE,
  ylab = "Fraction Suppressed",
  xlab = "Time Weeks Post Treatment",
  risk.table = FALSE
  #legend.labs=c("Luciferase", "VRC07")
)
print(survp$plot)
ggsave("../plots/mouse_exp/REJOc_survival_plot.png", units ="cm",width=25,height=10)

###
# Figure 2.7.B - REJO.c Cox Regression 
###

dff <- div_df %>%
  dplyr::filter(ref ==  "REJOc") %>% 
  dplyr::filter(week == 3) %>% 
  dplyr::filter(Antibody %in% c("VRC07","Luciferase")) %>% 
  dplyr::filter(!is.na(Status))

res.cox  <- coxph(Surv(Time,Status) ~ as.factor(Antibody) + log10(Hsn) + log10(filtered_reads) + log10(vl), data = dff)
summary(res.cox)
int <- summary(res.cox)$conf.int[,c(3,4)]
coef <- summary(res.cox)$coefficients[,c(2,5)]
cox_df <- as.data.frame(cbind(coef,int))
rownames(cox_df) <- c("VRC07","diversity", "reads","viral load")
colnames(cox_df) <- c("HR","p-value","CI-lower.95","CI-higher.95")
cox_df

summary(lm(average_entropy ~ filtered_reads, data = dff))
summary(lm(snv_count ~ filtered_reads, data = dff))
summary(lm(average_entropy ~ snv_count, data = dff))
#summary(lm(average_entropy ~ filtered_reads + ref , data = group_df))


###
# Figure 2.4 Compare differences in mutation composition between suppressed and escaped REJOc samples
###

plot_SNP_type_all <- function(){
  dff <- var_df %>%
    dplyr::filter(virus ==  "REJOc") %>% 
    dplyr::filter(Antibody %in% c("VRC07")) %>% 
    dplyr::filter(!is.na(Status))
  
  # Categorized Mutations
  dff <- dff %>% group_by(sample,ref) %>%
    dplyr::mutate(snp_type =  paste(ref, "->" ,alt,sep="")) %>%
    as_tibble()
  
  dff$Status <- as.factor(dff$Status)
  
  # count type specifc mutations 
  cov_summary_df <- dff %>%
    dplyr::filter(ref %in% c("A","T","G","C")) %>%
    dplyr::group_by(Status,snp_type) %>% 
    dplyr::count(name='muts') %>%
    as_tibble()
  
  #add missing vals
  #Status <- factor(c("1"))
  #snp_type <- c("C->G")
  #muts <- as.integer(c(0))
  #val <- data.frame(Status,snp_type,muts)
  
  #cov_summary_df <- rbind(cov_summary_df,val)
  #cov_summary_df$muts <- as.integer(cov_summary_df$muts)
  # count total mutations 
  total_muts_summary_df <- dff %>%
    dplyr::filter(ref %in% c("A","T","G","C")) %>%
    dplyr::group_by(Status) %>% 
    dplyr::count(name='total') %>%
    as_tibble()
  
  total_muts_summary_df$Status <- as.factor(total_muts_summary_df$Status)
  cov_summary_df$Status <- as.factor(cov_summary_df$Status)
  
  dff <- full_join(cov_summary_df,total_muts_summary_df,on=Status)
  dff$norm <- dff$muts / dff$total
  
  dff$snp_type <- factor(dff$snp_type,levels= c("A->C","A->G","A->T",
                                                "C->A","C->G","C->T",
                                                "G->A","G->C","G->T",
                                                "T->A","T->C","T->G",
                                                "DEL","IN"))
  my_comparisons <- list(c("0", "1"))
  dff$Status <- as.factor(dff$Status)
  dff$Status <- status_table[as.character(dff$Status)]
  #breaks=c("0","1"), labels=c("Suppressed", "Escaped")
  
  # plot 
  f <- ggplot(dff,aes(x = snp_type ,y=norm))  +
    geom_bar(stat="identity",aes(fill = Status), position=position_dodge()) + 
    theme_prism(border = TRUE) +  
    xlab("SNP Type") + 
    scale_color_manual(values=escape_group.colors_annotation) +
    scale_fill_manual(values=escape_group.colors_annotation) +
    ylab("Mutation Fraction") 
  f
  ggsave("../plots/mouse_exp/REJOc_Mutation_Type_Fraction.png",units ="cm",width=25,height=10)
  
  ### Chi Sq
  sup <- dff %>% filter(Status == "Suppressed") %>% select(norm,muts)
  esc <- dff %>% filter(Status == "Escaped") %>% select(norm,muts)
  
  dff <- data.frame(sup$muts,esc$muts)
  test <- chisq.test(dff)
  test
  print(summary(test))
  print(test)
}

plot_SNP_type_all()

plot_SNP_type_syn <- function(){
  dff <- var_df %>%
    dplyr::filter(virus ==  "REJOc") %>% 
    dplyr::filter(Antibody %in% c("VRC07")) %>% 
    dplyr::filter(SynNonSyn == "Syn") %>%
    dplyr::filter(!is.na(Status))
  
  # Categorized Mutations
  dff <- dff %>% group_by(sample,ref) %>%
    dplyr::mutate(snp_type =  paste(ref, "->" ,alt,sep="")) %>%
    as_tibble()
  
  dff$Status <- as.factor(dff$Status)
  
  # count type specific mutations 
  cov_summary_df <- dff %>%
    dplyr::filter(ref %in% c("A","T","G","C")) %>%
    dplyr::group_by(Status,snp_type) %>% 
    dplyr::count(name='muts') %>%
    as_tibble()
  
  #add missing vals
  Status <- factor(c("0","0","0","1","1"))
  snp_type <- c("C->G","G->C","T->G","C->G","G->C")
  muts <- as.integer(c(0))
  val <- data.frame(Status,snp_type,muts)
  
  cov_summary_df <- rbind(cov_summary_df,val)
  cov_summary_df$muts <- as.integer(cov_summary_df$muts)
  
  # count total mutations 
  total_muts_summary_df <- dff %>%
    dplyr::filter(ref %in% c("A","T","G","C")) %>%
    dplyr::group_by(Status) %>% 
    dplyr::count(name='total') %>%
    as_tibble()
  
  #total_muts_summary_df$Status <- as.factor(total_muts_summary_df$Status)
  cov_summary_df$Status <- as.factor(cov_summary_df$Status)
  
  dff <- full_join(cov_summary_df,total_muts_summary_df,on=Status)
  dff$norm <- dff$muts / dff$total
  
  dff$snp_type <- factor(dff$snp_type,levels= c("A->C","A->G","A->T",
                                                "C->A","C->G","C->T",
                                                "G->A","G->C","G->T",
                                                "T->A","T->C","T->G",
                                                "DEL","IN"))
  my_comparisons <- list(c("0", "1"))
  dff$Status <- as.factor(dff$Status)
  dff$Status <- status_table[as.character(dff$Status)]
  
  # plot 
  f <- ggplot(dff,aes(x = snp_type ,y=norm))  +
    geom_bar(stat="identity",aes(fill = Status), position=position_dodge()) + 
    theme_prism(border = TRUE) +  
    xlab("SNP Type") + 
    scale_color_manual(values=escape_group.colors_annotation) +
    scale_fill_manual(values=escape_group.colors_annotation) +
    ylab("Syn Mutation Fraction") 
  f
  ggsave("../plots/mouse_exp/REJOc_Mutation_Type_Sun_Fraction.png",units ="cm",width=25,height=10)
  
  ### Chi Sq
  sup <- dff %>% filter(Status == "Suppressed") %>% select(norm,muts)
  esc <- dff %>% filter(Status == "Escaped") %>% select(norm,muts)
  df <- data.frame(sup$muts,esc$muts)
  
  test <- chisq.test(x = df$sup.muts, y = df$esc.muts,simulate.p.value=TRUE,B = 10000)
  print(test)
  #print(summary(test))
}

plot_SNP_type_syn()




###
# Figure 2.6 Site Specific Diversity Plots - sliding window
###

#siding window function 

hsn <- function(freq,window=100){
  return(-1*sum(freq*log(freq))/window)
}

# site specific diversity function 
sliding_div <- function(dff,window=100){
  

  max_site <- 2563 #max(dff$pos)
  min_site <- 0#min(dff$pos)
  
  #scale by number of samples
  scale <- length(unique(dff$sample))
  
  start_site <- round((min_site + window/2))
  end_site <- round((max_site - window/2))
  
  output_df <- data.frame(site=integer(),
                          diversity=integer()) 
  
  for (site in seq(start_site,end_site)){
    #print(site)
    dff_window <- dff %>% filter(pos < site + round(window/2)) %>% filter(pos > site - round(window/2))
    diversity <-hsn(dff_window$freq,window) / scale #-1*sum(dff_window$freq*log(dff_window$freq)) / window
    tmp_df <- data.frame(site,diversity)
    output_df <- rbind(output_df,tmp_df)
    
  }
  return(output_df)
}

# site specific diversity function with bootstraps 
sliding_div_boot <- function(dff,window=100,bootstraps=100){
  
  max_site <- 2563 #max(dff$pos)
  min_site <- 0#min(dff$pos)
  
  #scale by number of samples
  scale <- length(unique(dff$sample))
  
  start_site <- round((min_site + window/2))
  end_site <- round((max_site - window/2))
  
  output_df <- data.frame(site=integer(),
                          diversity=integer()) 
  
  for (site in seq(start_site,end_site)){
    #print(site)
    dff_window <- dff %>% filter(pos < site + round(window/2)) %>% filter(pos > site - round(window/2))
    for (b in seq(0,bootstraps)){
      sample <- dff_window %>% sample_frac(0.1,replace=TRUE) 
      diversity <-hsn(sample$freq,window*0.1) / scale 
      tmp_df <- data.frame(site,diversity,b)
      output_df <- rbind(output_df,tmp_df)
    }
    
  }
  return(output_df)
}


## Suppressed Diversity 
Suppressed <-var_df %>%
  dplyr::filter(Antibody == "VRC07") %>%
  dplyr::filter(virus == "REJOc") %>%
  dplyr::filter(Status == "0")

slide_sup <- sliding_div_boot(Suppressed,window=100,bootstraps=10)
slide_sup$Status <- "Suppressed"

slide_sup_average <- slide_sup %>% group_by(site) %>%
  dplyr::summarise(average=mean(diversity),
                   error = qnorm(0.975)*mean(diversity)/sqrt(n()))
slide_sup_average$Status <- "Suppressed"
slide_sup_average$upper <- slide_sup_average$average + slide_sup_average$error
slide_sup_average$lower <- slide_sup_average$average - slide_sup_average$error


## Escaped Diversity 
Escaped <-var_df %>%
  dplyr::filter(Antibody == "VRC07") %>%
  dplyr::filter(virus == "REJOc") %>%
  dplyr::filter(Status == "1")

slide_exp <- sliding_div_boot(Escaped,window=100,bootstraps=10)

slide_exp$Status <- "Escaped"
slide_exp_average <- slide_exp %>% group_by(site) %>%
  dplyr::summarise(average=mean(diversity),
                   error = qnorm(0.975)*mean(diversity)/sqrt(n()))
slide_exp_average$Status <- "Escaped"
slide_exp_average$upper <- slide_exp_average$average + slide_exp_average$error
slide_exp_average$lower <- slide_exp_average$average - slide_exp_average$error


slide_df <- rbind(slide_exp_average,slide_sup_average)
#lide_df <- rbind(slide_exp,slide_sup)
#slide_df %>% group_by(site)


ggplot(slide_df, aes(x=average,fill=Status)) + geom_histogram(bins = 50)
#ggqqplot(slide_df$diversity)
#shapiro.test(slide_df$diversity[0:1000])

#plot 1
f <- ggplot(slide_df,aes(x = site ,y=average))  +
  geom_ribbon(aes(fill=Status,
    ymin = lower, ymax = upper,
    alpha = 0.1
  )) + 
  geom_line(aes(color = Status))  + 
  theme_prism(border = TRUE) +  
  scale_color_manual(values=escape_group.colors_annotation) +
  scale_fill_manual(values=escape_group.colors_annotation) +
  xlab("Position") + 
  ylab("Average Diversity") +
  
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 3)*(10^-rep(4:2, each = 9)),
    guide = "prism_minor"
  ) 
f

ggsave("../plots/mouse_exp/REJOc_Mutation_Site_Diversity.png",units ="cm",width=25,height=10)

## Corr Plot for each site
spread_df <- slide_df %>% tidyr::spread(Status, diversity)
spread_df$FC <- (spread_df$Escaped - spread_df$Suppressed)
cor(spread_df$Escaped,spread_df$Suppressed)
f <- ggplot(spread_df,aes(x=site,y=FC,color=site)) + 
  theme_prism(border = TRUE) + 
  geom_point() 
f
ggsave("../plots/mouse_exp/REJOc_Mutation_Site_Diversity_Corr_plot.png",units ="cm",width=25,height=10)


### 
# JR-CSF vs REJO.c Diversity 
###

#sample meta data
#group_SNV_df <- read_csv("../data/mouse_exp/grouped_SNV_data.csv")
#group_SNV_df <- group_SNV_df %>% select(c(sample,snv_count,mouse_number,ref,filtered_reads,exp,week,mouse,average_entropy))
#group_SNV_df$sample = paste(group_SNV_df$exp,group_SNV_df$mouse_number,sep="-")

#variant data 
#snv_df <- read_csv("../data/mouse_exp/var_data.csv")
#snv_df <- snv_df %>% dplyr::select(c(pos,ref,alt,pos_alt,freq,aa_pos,aa_ref,aa_alt,SynNonSyn,MutType,exp,sample,virus))
#snv_df$mouse_number <- snv_df$sample
#snv_df$sample = paste(snv_df$exp,snv_df$mouse_number,sep="-")

#survival data
#survival_df <- read_excel("../data/mouse_exp/survival_data.xlsx",sheet="Sheet1")
#colnames(survival_df) <-c("sample","Time","Antibody","Status")

#merge data frames 
#group_df <- dplyr::inner_join(survival_df, group_SNV_df, by = "sample") %>% dplyr::distinct()
#group_df <- dplyr::left_join(group_df,vl_df,by="sample")
#var_df <- dplyr::inner_join(survival_df, snv_df, by = "sample") %>% dplyr::distinct()

# include both VRC07 and Luciferase Ab groups 
dff <- group_df %>%
  dplyr::filter(ref %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) 



####
#### Figure 2.8
####


###
# Figure 2.8.A Plot Average Entropy by Virus Strains 
###

my_comparisons <- list(c("REJOc","JRCSF"))
f <- ggplot(dff,aes(x = ref ,y=average_entropy)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) +
  theme_prism() + 
  xlab("") + 
  ylab("Average Entropy (Hs)") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_JRCSF_diversity.png",units ="cm",width=15,height=10)


###
# Figure 2.8.C Plot SNV calls by Virus Strains 
###

my_comparisons <- list(c("REJOc","JRCSF"))
f <- ggplot(dff,aes(x = ref ,y=snv_count)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  ylab("Number of Called SNVs") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_JRCSF_SNV.png",units ="cm",width=15,height=10)


###
# Figure 2.8.B Plot Average Entropy of Syn Mutation - JRCSF vs REJOc
###

dff <- div_df %>%
  dplyr::filter(ref %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(Syn_Hsn > 0)
# 
my_comparisons <- list(c("REJOc","JRCSF"))

f <- ggplot(dff,aes(x = ref ,y=Syn_Hsn)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  scale_y_log10(
    limits = c(5e-6, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(6:3, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Synonymous  Mutation Entropy (Hs)") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_vs_JRCSF_syn_diversity.png",units ="cm",width=15,height=10)



###
# Figure 2.8.D Plot Average Coverage  - JRCSF vs REJOc
###

dff <- div_df %>%
  dplyr::filter(ref %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(filtered_reads > 3e3)
# 
my_comparisons <- list(c("REJOc","JRCSF"))

f <- ggplot(dff,aes(x = ref ,y=filtered_reads)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  theme_prism() + 
  scale_y_log10(
    limits = c(1e3, 1e6), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) +
  ylab("Coverage (Filtered Reads)") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_vs_JRCSF_coverage.png",units ="cm",width=15,height=10)


###
# Figure 2.8.E Plot Average Mutation Frequency  - JRCSF vs REJOc
###

dff <- div_df %>%
  dplyr::filter(ref %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) 
# 
my_comparisons <- list(c("REJOc","JRCSF"))

f <- ggplot(dff,aes(x = ref ,y=avg_freq)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  theme_prism() + 
  scale_y_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) +
  xlab("") + 
  ylab("Mutation Frequency") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_vs_JRCSF_mut_freq.png",units ="cm",width=15,height=10)


###
# Figure 2.8.E Plot dN/dS  - JRCSF vs REJOc
###

#### Compute Synonymous and NonSynonymous sites

REJOc_possible_mut_df <- read_csv("../data/mouse_exp/REJOc_possible_mut_df.csv")
REJOc_possible_mut_df$syn <- REJOc_possible_mut_df$ref_AA == REJOc_possible_mut_df$mut_AA

JRCSF_possible_mut_df <- read_csv("../data/mouse_exp/JRCSF_possible_mut_df.csv")
JRCSF_possible_mut_df$syn <- JRCSF_possible_mut_df$ref_AA == JRCSF_possible_mut_df$mut_AA

sites_df <- rbind(REJOc_possible_mut_df,JRCSF_possible_mut_df) 

syn_site <- sites_df %>% 
  group_by(virus) %>%
  count(syn) %>%
  spread(syn,n)
colnames(syn_site) <- c('virus','NonSynSite','SynSite')


dff <- var_df %>%
  dplyr::filter(virus %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(!is.na(Status))

syn_df <- dff %>% 
  group_by(sample,virus,SynNonSyn) %>%
  count(SynNonSyn) %>%
  spread(SynNonSyn,n)

syn_df[is.na(syn_df)] <- 0

#add site data
syn_df <- full_join(syn_df,syn_site,by="virus")

syn_df$ratio <- (syn_df$NonSyn/syn_df$NonSynSite) / (syn_df$Syn/syn_df$SynSite) 



### merge with survival df
dff <- dplyr::inner_join(syn_df,survival_df, by="sample") %>% 
  dplyr::distinct() %>%
  filter(Syn >0  ) 



f <- ggplot(dff,aes(x = virus ,y=ratio)) + 
  geom_boxplot(aes(color=virus)) +
  geom_point(aes(color=virus)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  scale_y_log10(
    limits = c(0.1, 10), 
    minor_breaks = rep(1:9, 4)*(10^rep(-1:1, each = 9)),
    guide = "prism_minor"
  ) +
  ylab("dN/dS") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_vs_JRCSF_seq_dNdS_vs_escape.png",units ="cm",width=15,height=10)

dff %>% 
  filter(Syn >0  ) %>% 
  group_by(virus) %>% 
  dplyr::summarise(mean(log10(ratio)))

###
# Viral Load
###

dff <- div_df %>%
  dplyr::filter(ref %in% c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3)
# 
my_comparisons <- list(c("REJOc","JRCSF"))

f <- ggplot(dff,aes(x = ref ,y=vl)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  theme_prism() + 
  scale_y_log10(
    limits = c(1e3, 1e7), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) +
  ylab("Viral Load") +
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_vs_JRCSF_viral_load.png",units ="cm",width=15,height=10)


###
# Mut Freq vs SNV calls
###
c(JRCSF = "#7570b3", REJOc ="#d95f02")

f <- ggplot(dff,aes(x = snv_count,y=avg_freq)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm",aes(color=ref)) +
  theme_prism() +
  xlab("Number of Called SNVs") + 
  ylab("Mutation Frequency") +
  scale_color_manual(values=c(JRCSF = "#7570b3", REJOc ="#d95f02")) +
  scale_fill_manual(values=c(JRCSF = "#7570b3", REJOc ="#d95f02")) +
  scale_y_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1, 100), 
    minor_breaks = rep(1:9, 4)*(10^rep(0:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_mut_freq_vs_diversity.png",units ="cm",width=15,height=10)

### corr data
df_cor <- dff %>% filter(ref == "REJOc")
cor(df_cor$avg_freq,df_cor$snv_count)

df_cor <- dff %>% filter(ref == "JRCSF")
cor(df_cor$avg_freq,df_cor$snv_count)

summary(lm(avg_freq ~ snv_count +ref, data=dff ))














###
# Test if input data correlates with observed outcomes 
###
group.colors <- c(JRCSF = "#7570b3", REJOc ="#d95f02")
dff <- div_df %>%
  dplyr::filter(ref %in%  c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) %>%
  dplyr::filter(!is.na(Status)) 
# 

f <- ggplot(dff,aes(x = filtered_reads,y=average_entropy)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Number of Input Reads") + 
  ylab("Average Entropy (Hs)") +
  #ggtitle("Diversity of REJO.c Escape") + 
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:2, each = 9)),
    guide = "prism_minor"
  ) +
  scale_x_log10(
    limits = c(1e3, 1e6), 
    minor_breaks = rep(1:9, 4)*(10^rep(3:6, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_input_reads_vs_diversity.png",units ="cm",width=15,height=10)

cor.test(dff$vl,dff$average_entropy)
cor.test(dff$filtered_reads,dff$average_entropy)
cor.test(dff$snv_count,dff$average_entropy)
cor.test(dff$avg_freq,dff$average_entropy)

cor.test(dff$snv_count,dff$avg_freq)


###
# SNV vs Hsn
###
f <- ggplot(dff,aes(x = snv_count,y=average_entropy)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Number of Called SNVs") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1, 100), 
    minor_breaks = rep(1:9, 4)*(10^rep(0:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_snpcount_vs_diversity.png",units ="cm",width=15,height=10)

###
# SNV vs Average Freq
###
f <- ggplot(dff,aes(x = snv_count,y=avg_freq)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm",aes(color=Status)) +
  theme_prism() +
  xlab("Number of Called SNVs") + 
  ylab("Mutation Frequency") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1, 100), 
    minor_breaks = rep(1:9, 4)*(10^rep(0:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_snpcount_vs_freq.png",units ="cm",width=15,height=10)

###
# regression to explain this patterns
###
summary(lm(avg_freq ~ snv_count + log10(filtered_reads) , data = dff))
summary(lm(avg_freq ~ snv_count + log10(filtered_reads) , data = dff))

###
# VL vs Hsn
###
f <- ggplot(dff,aes(x = vl,y=average_entropy)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm",aes(color=ref)) +
  theme_prism() +
  xlab("Viral Load") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1e4, 1e7), 
    minor_breaks = rep(1:9, 4)*(10^rep(4:7, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_VL_vs_diversity.png",units ="cm",width=15,height=10)


###
# Average freq vs Hsn
###
f <- ggplot(dff,aes(x = avg_freq,y=average_entropy)) + 
  geom_point(aes(color=ref),size=6) +
  #geom_smooth(method = "lm") +
  theme_prism() +
  xlab("Average Mutation Frequency") + 
  ylab("Average Entropy (Hs)") +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:3, each = 9)),
    guide = "prism_minor"
  ) + 
  scale_x_log10(
    limits = c(1e-3, 1e0), 
    minor_breaks = rep(1:9, 4)*(10^-rep(4:1, each = 9)),
    guide = "prism_minor"
  ) 

f 
ggsave("../plots/mouse_exp/REJOc_JRCSF_mut_freq_vs_diversity.png",units ="cm",width=15,height=10)


###
# Logistic Regression of Input Variables 
###

dff <- div_df %>%
  dplyr::filter(ref %in% c("JRCSF","REJOc")) %>% 
  dplyr::filter(week == 3)
dff$ref <- factor(dff$ref,levels=c("JRCSF","REJOc"))
#logistic_regression <- glm(ref ~ log10(average_entropy) + log10(filtered_reads)  + log10(vl) , data = dff, family = binomial(link = "logit"))
logistic_regression <- glm(ref ~ log10(average_entropy) + log10(filtered_reads)  + log10(vl) + snv_count + avg_freq, data = dff, family = binomial(link = "logit"))
summary(logistic_regression)
confint(logistic_regression)




##### 
##### OLD
#####
# mean diversity 
dff %>% 
  group_by(ref) %>%
  dplyr::summarise(mean(average_entropy))

###
# Figure 2.8.C Plot SNV calls by Virus Strains 
###
dff <- var_df %>%
  dplyr::filter(virus %in%c("REJOc",'JRCSF'))

boot_df <- dff %>% 
  group_by(sample) %>%
  summarise(avg_freq = mean(freq))

group_SNV_df
### merge with survival df
dff <- dplyr::inner_join(boot_df,group_SNV_df, by="sample") %>% dplyr::distinct()


my_comparisons <- list(c("REJOc","JRCSF"))
f <- ggplot(dff,aes(x = ref ,y=avg_freq)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  ylab("Average Freq of SNVs") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
#ggsave("../plots/mouse_exp/REJOc_JRCSF_SNV.png",units ="cm",width=15,height=10)


###
# Figure 2.8.B Plot SNV freq calls by Virus Strains 
###

my_comparisons <- list(c("REJOc","JRCSF"))
f <- ggplot(dff,aes(x = ref ,y=snv_count)) + 
  geom_boxplot(aes(color=ref)) +
  geom_point(aes(color=ref)) +
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  theme_prism() + 
  xlab("") + 
  ylab("Number of Called SNVs") +
  # ggtitle("Diversity of REJO.c Escape") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test",label = "p.signif",comparisons = my_comparisons)
ggsave("../plots/mouse_exp/REJOc_JRCSF_SNV.png",units ="cm",width=15,height=10)




###
# 2.8.A JR-CSF vs REJO.c Plot Hs Diversity
###

f <- ggplot(dff,aes(x = ref ,y=average_entropy,color = ref)) + 
  geom_boxplot() +
  geom_point() +
  scale_color_manual(values=group.colors) + 
  theme_prism() +
  xlab("") + 
  ylab("Average Entropy (Hs)") +
  scale_x_discrete(breaks=c("JRCSF","REJOc"),
                   labels=c("JR-CSF", "REJO.c")) + 
  ggtitle("Diversity of JR-CSF vs REJO.c") + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test")
ggsave("../plots/mouse_exp/REJOc_JRCSF_diversity.png",units ="cm",width=15,height=10)

# mean diversity 
dff %>% 
  group_by(ref) %>%
  dplyr::summarise(mean(average_entropy))

###
# 2.8.B JR-CSF vs REJO.c Plot SNV counts
###
f <- ggplot(dff,aes(x = ref ,y=snv_count,color = ref)) + 
  geom_boxplot() +
  geom_point() +
  scale_color_manual(values=group.colors) + 
  theme_prism() +
  xlab("") + 
  ylab("Number of Called SNVs") +
  ggtitle("SNV Calls of JR-CSF vs REJO.c") + 
  scale_x_discrete(breaks=c("JRCSF","REJOc"),
                   labels=c("JR-CSF", "REJO.c")) + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test")
ggsave("../plots/mouse_exp/REJOc_JRCSF_SNV.png",units ="cm",width=15,height=10)

# mean SNVs 
dff %>% 
  group_by(ref) %>%
  dplyr::summarise(mean(snv_count))

summary(glm(log10(average_entropy) ~ log10(filtered_reads) + ref + log10(snv_count),data = dff))


###
# 2.8.C JR-CSF vs REJO.c Plot Average SNV frequency
### 
snv_freq_df <- snv_df %>% 
  group_by(sample) %>%
  summarise(avg_freq = mean(freq))

### merge with survival df
dff <- dplyr::inner_join(snv_freq_df,group_df, by="sample") %>% dplyr::distinct()


f <- ggplot(dff,aes(x = ref ,y=avg_freq,color = ref)) + 
  geom_boxplot() +
  geom_point() +
  scale_color_manual(values=group.colors) + 
  theme_prism() +
  xlab("") + 
  ylab("Frequency of SNVs") +
  ggtitle("SNV Calls of JR-CSF vs REJO.c") + 
  scale_x_discrete(breaks=c("JRCSF","REJOc"),
                   labels=c("JR-CSF", "REJO.c")) + 
  theme(legend.position = "none") 

f + stat_compare_means(method="wilcox.test")







###
# Plot site specific entropy
###

## Suppressed Diversity 
R_site <-var_df %>%
  dplyr::filter(virus == "REJOc") 

slide_R <- sliding_div(R_site,window=100)
slide_R$virus <- "REJOc"

## Escaped Diversity 
J_site <-var_df %>%
  dplyr::filter(virus == "JRCSF")

slide_J <- sliding_div(J_site,window=100)
slide_J$virus <- "JRCSF"

slide_df <- rbind(slide_J,slide_R)

# set colors
group.colors <- c(JRCSF = "#7570b3", REJOc ="#d95f02")

# plot 
f <- ggplot(slide_df,aes(x = site ,y=diversity))  +
  geom_line(aes(color = virus))  + 
  theme_prism(border = TRUE) +  
  scale_color_manual(values=group.colors) +
  scale_fill_manual(values=group.colors) +
  xlab("Position") + 
  ylab("Average Diversity") +
  
  scale_y_log10(
    limits = c(1e-4, 1e-2), 
    minor_breaks = rep(1:9, 3)*(10^-rep(4:2, each = 9)),
    guide = "prism_minor"
  ) 
f

ggsave("../plots/mouse_exp/REJOc_JRCSF_site_entropy.png",units ="cm",width=25,height=10)


### corr plot

## Corr Plot for each site
spread_df <- slide_df %>% tidyr::spread(virus, diversity)
cor(spread_df$JRCSF,spread_df$REJOc)
f <- ggplot(spread_df,aes(x=JRCSF,y=REJOc,color=site)) + 
  theme_prism(border = TRUE) + 
  geom_point() + 
  geom_smooth(method='lm')
f
ggsave("../plots/mouse_exp/REJOc_JRCSF_site_entropy_Corr_plot.png",units ="cm",width=25,height=10)



summary(lm(diversity ~ virus,data=slide_df))
###
# REJOc vs JRCSF Mutational Profile
###

dff <- var_df %>%
  dplyr::filter(virus %in% c("JRCSF","REJOc")) 

# Categorized Mutations
dff <- dff %>% group_by(sample,ref) %>%
  dplyr::mutate(snp_type =  paste(ref, "->" ,alt,sep="")) %>%
  as_tibble()


# count type specifc mutations 
cov_summary_df <- dff %>%
  dplyr::filter(ref %in% c("A","T","G","C")) %>%
  dplyr::group_by(virus,snp_type) %>% 
  dplyr::count(name='muts') %>%
  as_tibble()

total_muts_summary_df <- dff %>%
  dplyr::filter(ref %in% c("A","T","G","C")) %>%
  dplyr::group_by(virus) %>% 
  dplyr::count(name='total') %>%
  as_tibble()


dff <- full_join(cov_summary_df,total_muts_summary_df,on=virus)
dff$norm <- dff$muts / dff$total

dff$snp_type <- factor(dff$snp_type,levels= c("A->C","A->G","A->T",
                                              "C->A","C->G","C->T",
                                              "G->A","G->C","G->T",
                                              "T->A","T->C","T->G",
                                              "DEL","IN"))
dff$virus <- as.factor(dff$virus)
my_comparisons <- list(c("REJOc", "JRCSF"))
group.colors <- c(JRCSF = "#7570b3", REJOc ="#d95f02")

# plot 
f <- ggplot(dff,aes(x = snp_type ,y=norm,fill=virus))  +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_prism(border = TRUE) +  
  scale_fill_manual(values=group.colors) + 
  xlab("SNP Type") + 
  ylab("Mutation Fraction") 
f
ggsave("../plots/mouse_exp/REJOc_JRCSF_Mutation_Type_Fraction.png",units ="cm",width=25,height=10)

R <- dff %>% filter(virus == "REJOc") %>% select(norm,muts)
J <- dff %>% filter(virus == "JRCSF") %>% select(norm,muts)
dff <- data.frame(R$muts,J$muts)

test <- chisq.test(dff)
#test <- chisq.test(R)

test 
summary(test)


###
# Synonymouse Mutation SNP Types
###

dff <- var_df %>%
  dplyr::filter(virus %in% c("JRCSF","REJOc"))  %>% 
  dplyr::filter(SynNonSyn == "Syn")

# Categorized Mutations
dff <- dff %>% group_by(sample,ref) %>%
  dplyr::mutate(snp_type =  paste(ref, "->" ,alt,sep="")) %>%
  as_tibble()


# count type specifc mutations 
cov_summary_df <- dff %>%
  dplyr::filter(ref %in% c("A","T","G","C")) %>%
  dplyr::group_by(virus,snp_type) %>% 
  dplyr::count(name='muts') %>%
  as_tibble()

total_muts_summary_df <- dff %>%
  dplyr::filter(ref %in% c("A","T","G","C")) %>%
  dplyr::group_by(virus) %>% 
  dplyr::count(name='total') %>%
  as_tibble()

#add missing vals
virus <- factor(c("REJOc","REJOc"))
snp_type <- c("C->G","G->C")
muts <- as.integer(c(0))
val <- data.frame(virus,snp_type,muts)

cov_summary_df <- rbind(cov_summary_df,val)
cov_summary_df$muts <- as.integer(cov_summary_df$muts)


dff <- full_join(cov_summary_df,total_muts_summary_df,on=virus)
dff$norm <- dff$muts / dff$total

dff$snp_type <- factor(dff$snp_type,levels= c("A->C","A->G","A->T",
                                              "C->A","C->G","C->T",
                                              "G->A","G->C","G->T",
                                              "T->A","T->C","T->G",
                                              "DEL","IN"))
dff$virus <- as.factor(dff$virus)
my_comparisons <- list(c("REJOc", "JRCSF"))
group.colors <- c(JRCSF = "#7570b3", REJOc ="#d95f02")

# plot 
f <- ggplot(dff,aes(x = snp_type ,y=norm,fill=virus))  +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_prism(border = TRUE) +  
  scale_fill_manual(values=group.colors) + 
  xlab("SNP Type") + 
  ylab("Syn Mutation Fraction") 
f
ggsave("../plots/mouse_exp/REJOc_JRCSF_Syn_Mutation_Type_Fraction.png",units ="cm",width=25,height=10)

R <- dff %>% filter(virus == "REJOc") %>% select(norm,muts)
J <- dff %>% filter(virus == "JRCSF") %>% select(norm,muts)
dff <- data.frame(R$muts,J$muts)

test <- chisq.test(dff)
test <- chisq.test(R)

test 
summary(test)


###
# REJOc vs JRCSF survival analysis
###

###
# KM Survival Plot - REJO.c vs JRCSF
###
dff <- group_df %>%
  dplyr::filter(ref %in%  c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) %>% 
  dplyr::filter(Antibody %in% c("VRC07")) %>% 
  dplyr::filter(!is.na(Status))

fit <- survfit(Surv(Time, Status) ~ ref, data = dff)

survp<- ggsurvplot(
  fit, 
  data=dff, 
  conf.int = FALSE,
  pval = TRUE,
  palette = c("#7570b3",'#d95f02'),
  ggtheme = theme_bw(),
  ncensor.plot = TRUE,
  ylab = "Fraction Suppressed",
  xlab = "Time Weeks",
  risk.table = FALSE
  
  
  #legend.labs=c("Luciferase", "VRC07")
)

print(survp$plot)
ggsave("../plots/mouse_exp/REJOc_JRCSF_survival_plot.png", units ="cm",width=25,height=10)


###
# Cox Regression 
###

dff <- group_df %>%
  dplyr::filter(ref %in%  c("REJOc","JRCSF")) %>% 
  dplyr::filter(week == 3) %>% 
  dplyr::filter(!is.na(Status))

#re-order reference 
dff$ref <- factor(dff$ref,levels = c("REJOc","JRCSF"))

res.cox  <- coxph(Surv(Time,Status) ~  log10(average_entropy) + log10(vl) + log10(filtered_reads) , data = dff)
summary(res.cox)
int <- summary(res.cox)$conf.int[,c(3,4)]
coef <- summary(res.cox)$coefficients[,c(2,5)]
cox_df <- as.data.frame(cbind(coef,int))
rownames(cox_df) <- c("JRCSF","hsn diversity", "reads")
colnames(cox_df) <- c("HR","p-value","CI-lower.95","CI-higher.95")
cox_df

model <- glm(ref ~ log10(filtered_reads) + log10(snv_count) + log10(average_entropy), data=dff,family = "binomial")
summary(model)
confint(model)

summary(glm(average_entropy ~ filtered_reads + ref, data = dff))
summary(lm(snv_count ~ filtered_reads, data = dff))
summary(lm(average_entropy ~ snv_count, data = dff))
