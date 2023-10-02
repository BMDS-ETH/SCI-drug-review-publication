################################################################################
# SCI - Drug review - Numbers reported in the paper
# L. Bourguignon & L.P. Lukas
# First version : 26.04.2023
# Last update : 25.07.2023
################################################################################

# Load data

raw_data_review <- read.csv('./data/data-extracted-20230818.csv')

################################################################################
# Libraries

library(treemapify)
library(ggplot2)
library(waffle)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(devtools)
library(dplyr)
library(tidyr)
library(ggridges)
library(stringr)
library(gridExtra)
library(grid)
library(lattice)
library(forcats)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(ggstats)
library(svglite)

################################################################################
# Functions

unique_doi <- function(data){
  new_data <- data %>% 
    group_by(DOI.or.PMID) %>% 
    filter(rank(DOI.or.PMID, ties.method="first")==1) %>% # keep the first instance when duplicate DOI
    distinct %>%
    dplyr::ungroup() # need to ungroup for future steps
  return(new_data)
}

waffle_plot <- function(data, column){
  number_colors <- length(levels(factor(data[[column]])))
  plot <- waffle(table(data[[column]]), 
                 colors = brewer.pal(n = number_colors, name = "Paired"))
  return(plot)
}

barplot_plot <- function(data, column){
  plot <- ggplot(data = data, aes_string(x = column, y = 'n')) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9),
              vjust = -0.25, size = 6) +
    labs(x = column) +
    theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold")) +
    theme_half_open() +
    background_grid()
  return(plot)
}

levels_prevalence <- function(data, column){
  data[[column]] <- reorder(data[[column]], data[[column]], FUN = length)
  data[[column]] <- fct_rev(data[[column]])
  return(data[[column]])
}


counts_table <- function(data, column){
  counts <- data %>% 
    dplyr::count(across(column)) %>%
    mutate(freq = n / sum(n))
  return(counts)
}

plots_column <- function(data, column){
  data[[column]] <- levels_prevalence(data, column)
  counts <- counts_table(data, column)
  waffle <- waffle_plot(data, column)
  bar_plot <- barplot_plot(counts, column)
  return(list(counts, waffle, bar_plot))
}

harmonised_mixed_effects <- function(data){
  data$Drug.effect.on.functional.assessment <- as.character(data$Drug.effect.on.functional.assessment)
  data$Drug.effect.on.functional.assessment[grepl('mixed', data$Drug.effect.on.functional.assessment)] <- 'mixed'
  data$Drug.effect.on.functional.assessment[data$Drug.effect.on.functional.assessment == ''] <- NA
  data$Drug.effect.on.functional.assessment[data$Drug.effect.on.functional.assessment == '-'] <- NA
  return (data)
}

################################################################################
## Publications over time

test_doi <- unique_doi(raw_data_review)
print(paste('Number of unique DOI full-text publications studied:', dim(test_doi)[1]))


included_review <- raw_data_review %>% filter(Included.exclude == 'included')
included_unique_doi_all <- unique_doi(included_review)
print(paste('Number of unique publications included:', dim(included_unique_doi_all)[1]))
print(paste('Number of experiments included:', dim(included_review)[1]))
print(paste('Number of publications that contribute to more than 1 experiment:', 
            sum(table(included_review$DOI.or.PMID) > 1))) 
print(paste('% of publications that contribute to more than 1 experiment:', 
            100*sum(table(included_review$DOI.or.PMID) > 1)/dim(included_unique_doi_all)[1]))

excluded_review <- raw_data_review %>% filter(Included.exclude == 'excluded')
excluded_unique_doi_all <- unique_doi(excluded_review)
# table(excluded_unique_doi_all$Reason.for.exclusion)
# sum(table(excluded_unique_doi_all$Reason.for.exclusion.1))
# 
# intersect(included_unique_doi_all$DOI.or.PMID, excluded_unique_doi_all$DOI.or.PMID)

################################################################################
### Supplementary Figure 1A ###
# Number of unique publications over time

count_year_included <- included_unique_doi_all %>% dplyr::count(Year)
count_year_included$Year <- as.numeric(count_year_included$Year)
barplot_years <- ggplot(data = count_year_included, 
                        aes(x = Year, y = n)) +
  geom_bar(stat = "identity") +
  xlab('Year of publication') +
  ylab('Number of publications') +
  theme_half_open() +
  background_grid() + 
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'), 
        legend.key = element_rect(fill = "black", color = NA))
barplot_years

count_year_included_after2010 <- count_year_included %>% filter(Year >= 2010)

print(paste('Number of publications published in or after 2010:', sum(count_year_included_after2010$n)))
print(paste('% of publications that contribute to more than 1 experiment:', 
            100*sum(count_year_included_after2010$n)/dim(included_unique_doi_all)[1]))

#ggsave("./figures/paper-preparation/SupplementaryFigure1.png", barplot_years, width = 6, height=7, dpi=300, units = "in")
################################################################################

print('Languages')
table(included_unique_doi_all$Language)

included_review_animals <- included_review[included_review$Species != 'human', ]
included_review_unique_animals <- included_unique_doi_all[included_unique_doi_all$Species != 'human', ]

print(paste('Number of unique publications included animals:', dim(included_review_unique_animals)[1]))
print(paste('Number of unique experiments included animals:', dim(included_review_animals)[1]))

included_review_humans <- included_review[included_review$Species == 'human', ]
included_review_unique_humans <- included_unique_doi_all[included_unique_doi_all$Species == 'human', ]

print(paste('Number of unique publications included human:', dim(included_review_unique_humans)[1]))
print(paste('Number of unique experiments included human:', dim(included_review_humans)[1]))

################################################################################
## Drugs studied

print(paste('Number of drugs or combinations tested:', dim(table(included_review$Drug.name.harmonized))))

drugs_list <- data.frame(unclass(table(included_review$Drug.name.harmonized)))
#write.csv(drugs_list, file='./tables/TableS1.csv', row.names=TRUE)

count_drugs_comb <- dplyr::filter(included_review, grepl('+', Drug.name.harmonized, fixed = TRUE)) %>% dplyr::count(Drug.name.harmonized)
print(paste('Number of drug combination tested:', dim(count_drugs_comb)[1]))

drugs_interest <- read.csv('./data/drugs-of-interest-0205.csv')
drugs_interest_tested <- drugs_interest %>%
  filter(Was.this.drug.tested.in.any.of.the.papers.extracted...binary..1.yes. == 1)
name_drugs_interest_tested <- drugs_interest_tested$n

################################################################################
## Animal studies

table(included_review_animals$Species)
# other: 'Yucatan miniature pigs' (n=2) 
# and 'yellow eel Anguilla anguilla L.' (n=1)

################################################################################
### Supplementary Figure 1 ###
# Number of unique publications over time
library(scales)
colours <- viridis_pal()(8)

count_year_included <- included_review_animals %>% dplyr::count(Year, Species)
count_year_included$Year <- as.numeric(count_year_included$Year)
count_year_included$Species <- factor(count_year_included$Species, 
                                      levels = c("multiple", "other", 
                                                 "guinea pig", "dogs", 'cats',
                                                 'rabbit', 'mice', 'rats'))

# Number of unique publications over time 
barplot_years <- ggplot(data = count_year_included, 
                                aes(x = Year, y = n)) +
  geom_bar(stat = "identity") +
  xlab('Year of publication') +
  ylab('Number of experiments') +
  theme_half_open() +
  background_grid()+
  guides(fill = guide_legend(reverse = TRUE))+ 
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'),
        legend.text = element_text(color="black", size=10, family = 'Avenir'),
        legend.title = element_text(color="black", size=12, family = 'Avenir'))
barplot_years
ggsave("./figures/SuppFigure1/SupplementaryFigure1A.svg", barplot_years, width = 7, height=5, dpi=300, units = "in")


# Number of unique publications over time coloured by species
barplot_years_species <- ggplot(data = count_year_included, 
                        aes(x = Year, y = n, fill = Species)) +
  scale_fill_manual(values = c("#949798", "#690608","#B6503C", "#EB9E4B","#98A285", 
                               "#44809A", "#00568C")) +
  geom_bar(stat = "identity") +
  xlab('Year of publication') +
  ylab('Number of experiments') +
  theme_half_open() +
  background_grid()+
  guides(fill = guide_legend(reverse = TRUE))+ 
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'),
        legend.text = element_text(color="black", size=10, family = 'Avenir'),
        legend.title = element_text(color="black", size=12, family = 'Avenir'))
barplot_years_species
ggsave("./figures/SuppFigure1/SupplementaryFigure1B.svg", barplot_years_species, width = 7, height=5, dpi=300, units = "in")

# Number of unique publications over time coloured by sex
count_year_included_sex <- included_review_animals %>% dplyr::count(Year, Sex)
count_year_included_sex$Year <- as.numeric(count_year_included_sex$Year)
count_year_included_sex$Sex <- factor(count_year_included_sex$Sex, 
                                      levels = c("mixed","not reported",
                                                 "female", "male"))
barplot_years_sex <- ggplot(data = count_year_included_sex, 
                        aes(x = Year, y = n, fill = Sex)) +
  scale_fill_manual(values = c("#c86243","#d8b161", "#949798","#00547B")) +
  geom_bar(stat = "identity") +
  xlab('Year of publication') +
  ylab('Number of experiments') +
  theme_half_open() +
  background_grid()+
  guides(fill = guide_legend(reverse = TRUE))+ 
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'),
        legend.text = element_text(color="black", size=10, family = 'Avenir'),
        legend.title = element_text(color="black", size=12, family = 'Avenir'))
barplot_years_sex
ggsave("./figures/SuppFigure1/SupplementaryFigure1C.svg", barplot_years_sex, width = 7, height=5, dpi=300, units = "in")

barplot_years_sex_prop <- ggplot(data = count_year_included_sex, 
                                 aes(x = Year, fill = factor(Sex))) +
  scale_fill_manual(values = c("#902C25","#EB9E4B", "#98A285","#00547B")) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(.5)) + 
  xlab('Year of publication') +
  ylab('Number of experiments') +
  theme_half_open() +
  background_grid()+
  guides(fill = guide_legend(reverse = TRUE))
barplot_years_sex_prop

subset_sex_year <- included_review_animals %>% select(Year, Sex)
ggplot(data = subset_sex_year, aes(x=Year)) +
  geom_histogram(aes(color=Sex, fill=Sex),
                 position = "identity", bins = 30, alpha = 0.2) +
  scale_color_manual(values = c("#902C25","#EB9E4B", "#98A285","#00547B")) +
  scale_fill_manual(values = c("#902C25","#EB9E4B", "#98A285","#00547B"))

count_year_included_sex_complete <- count_year_included_sex %>% complete(Year, Sex)
ggplot(count_year_included_sex_complete,
       aes(x = factor(Year),
           y = n,
           fill = Sex)) +
  geom_bar(stat = "identity",
           position = "dodge")+
  scale_fill_manual(values = c("#902C25","#EB9E4B", "#98A285","#00547B"))

################################################################################

included_unique_doi_bounded.sample.size <- dplyr::filter(included_review_animals, grepl('>', Count..n))
included_unique_doi_no.sample.size <- included_review_animals[included_review_animals$Count..n == 'Not reported',]
included_unique_doi_ranges <- dplyr::filter(included_review_animals, grepl('-', Count..n))
included_unique_doi_num.sample.size <- included_review_animals[!is.na(as.numeric(included_review_animals$Count..n)), ]
included_unique_doi_num.sample.size$Count..n <- as.numeric(included_unique_doi_num.sample.size$Count..n)

print(paste('Number of experiments with sample size (partially) missing:', 
            dim(included_unique_doi_bounded.sample.size)[1] + 
              dim(included_unique_doi_ranges)[1] +
              dim(included_unique_doi_no.sample.size)[1]))

included_unique_doi_no.age <- included_review_animals[included_review_animals$Age..comment. == 'Not reported',]
print(paste('Number of publications with no age reported:', dim(included_unique_doi_no.age)[1]))
included_unique_doi_adult <- included_review_animals[included_review_animals$Age..comment. == 'Adult',]
print(paste('Number of publications with adult reported, note that some have both adult and specific age:', dim(included_unique_doi_adult)[1]))
included_unique_doi_age.range <- included_review_animals %>% filter(!is.na(Age..min.))
print(paste('Number of publications with age range reported:', dim(included_unique_doi_age.range)[1]))
included_unique_doi_age.mean <- included_review_animals %>% filter(!is.na(Age..mean.))
print(paste('Number of publications with age mean reported:', dim(included_unique_doi_age.mean)[1]))

included_unique_doi_adult_only <- included_review_animals %>% 
  filter(Age..comment. %in% c('Adult', 'Young')) %>% 
  filter(Age..units. == '')
  
print(paste('Number of publications with adult or young reported only:',
            dim(included_unique_doi_adult_only)[1]))


results_sex <- plots_column(included_review_animals, 'Sex')
results_sex[[1]]
results_sex[[2]]
results_sex[[3]]

included_mixed_sex_animals <- included_review_animals %>% filter(Sex == 'mixed')
table(included_mixed_sex_animals$Sex.....male.)

included_noinfo_exclusion_animals <- included_review_animals %>% 
  filter(Count..n.died.in.control.group == '',
         Count..n.excluded.in.control.group == '',
         Count..n.died.in.treatment.group == '',
         Count..n.excluded.in.treatment.group == '')

print(paste('Number of publications reporting exclusion/death animals:',
            dim(included_review_animals)[1] - dim(included_noinfo_exclusion_animals)[1]))

mean(included_unique_doi_num.sample.size$Count..n)
sd(included_unique_doi_num.sample.size$Count..n)
median(included_unique_doi_num.sample.size$Count..n)
quantile(included_unique_doi_num.sample.size$Count..n)

################################################################################
### Supplementary Table 3 ###


included_num.sample.size <- included_review_animals[!is.na(as.numeric(included_review_animals$Count..n)), ]
included_num.sample.size$Count..n <- as.numeric(included_num.sample.size$Count..n)
sub_included_num.sample.size <- included_num.sample.size[,c('Count..n', 'Species')]
stats_sample.size.per.species <- aggregate(. ~ Species, data.frame(sub_included_num.sample.size), 
                                           function(x) c(mean = round(mean(x), digits=2), 
                                                         sd = round(sd(x), digits=2),
                                                         median = round(quantile(x, probs = c(0.5)), digits=2),
                                                         Q1 = round(quantile(x, probs = c(0.25)), digits=2),
                                                         Q3 = round(quantile(x, probs = c(0.75)), digits=2)))
stats_sample.size.per.species

#write.table(stats_sample.size.per.species, file = "./figures/paper-preparation/SupplementaryTable3.txt", sep = ",", quote = FALSE, row.names = F)

################################################################################


included_age.range.week <- included_review_animals[included_review_animals$Age..units. == 'Weeks',]
included_age.range.days <- included_review_animals[included_review_animals$Age..units. == 'Days',]
included_age.range.months <- included_review_animals[included_review_animals$Age..units. == 'Months',]
sub_included_num.age.range.week <- included_age.range.week[,c('Age..min.', 'Age..max.', 'Species')]
sub_included_num.age.mean.week <- included_age.range.week[,c('Age..mean.', 'Age..SD.', 'Species')]

# Rats
sub_included_num.age.mean.week.rats <- sub_included_num.age.mean.week[sub_included_num.age.mean.week$Species == 'rats',]
dim(sub_included_num.age.mean.week.rats)[1]
mean(sub_included_num.age.mean.week.rats$Age..mean., na.rm = T)
dim(sub_included_num.age.mean.week.rats)[1] - sum(is.na(sub_included_num.age.mean.week.rats$Age..mean.))

sub_included_num.age.range.week.rats <- sub_included_num.age.range.week[sub_included_num.age.range.week$Species == 'rats',]
mean(sub_included_num.age.range.week.rats$Age..min., na.rm = T)
mean(sub_included_num.age.range.week.rats$Age..max., na.rm = T)
dim(sub_included_num.age.range.week.rats)[1] - sum(is.na(sub_included_num.age.range.week.rats$Age..min.))

# Mice
sub_included_num.age.mean.week.mice <- sub_included_num.age.mean.week[sub_included_num.age.mean.week$Species == 'mice',]
dim(sub_included_num.age.mean.week.mice)[1]
mean(sub_included_num.age.mean.week.mice$Age..mean., na.rm = T)
dim(sub_included_num.age.mean.week.mice)[1] - sum(is.na(sub_included_num.age.mean.week.mice$Age..mean.))

sub_included_num.age.range.week.mice <- sub_included_num.age.range.week[sub_included_num.age.range.week$Species == 'mice',]
mean(sub_included_num.age.range.week.mice$Age..min., na.rm = T)
mean(sub_included_num.age.range.week.mice$Age..max., na.rm = T)
dim(sub_included_num.age.range.week.mice)[1] - sum(is.na(sub_included_num.age.range.week.mice$Age..min.))


# mean(sub_included_num.age.mean.week.rats$Age..mean., na.rm = T)
# stats_age_weeks <- aggregate(. ~ Species, data.frame(sub_included_num.age.mean.week), 
#                              function(x) c(mean = mean(x), 
#                                            sd = sd(x),
#                                            median = quantile(x, probs = c(0.5)),
#                                            Q1 = quantile(x, probs = c(0.25)),
#                                            Q3 = quantile(x, probs = c(0.75)),
#                                            min = min(x),
#                                            max = max(x)))
# stats_age_weeks
# 
# included_review.week <- included_review_animals %>% 
#   filter(Age..units. == 'Weeks', Species == 'rats')
# included_review.days <- included_review_animals %>% 
#   filter(Age..units. == 'Days', Species == 'rats')
# included_review.months <- included_review_animals %>% 
#   filter(Age..units. == 'Months', Species == 'rats')
# 
# included_review.week %>% select(Age..min., Age..max., Age..mean.)
# mean(included_review.week$Age..min., na.rm = T)

################################################################################
# Injury characteristics

results_inj.mechanism <- plots_column(included_review_animals, 'Injury.mechanism')
results_inj.mechanism[[1]]
results_inj.mechanism[[2]]
results_inj.mechanism[[3]]

results_inj.severity <- plots_column(included_review_animals, 'Injury.severity')
results_inj.severity[[1]]
results_inj.severity[[2]]
results_inj.severity[[3]] + scale_x_discrete(labels=c("Not reported" = "Not\nspecifically\nreported",
                                                      "moderate - severe" = "moderate\nto\nsevere"))

#Level of injury

all_levels <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
                "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13",
                "L1", "L2", "L3", "L4", "L5")

# for (level in c("C1", "C3", "C8", "T2", "T5", "L5")){
#   levels(included_review$Injury.level)[length(levels(included_review$Injury.level))+1] <- level
# }

count_injury.level <- included_review_animals %>% dplyr::count(Injury.level)
count_injury.level$y <- 1
count_injury.level_unique <- count_injury.level %>% mutate(nb_ch = nchar(Injury.level)) %>% dplyr::filter(nb_ch <= 3)
count_injury.level_unique$Injury.level <- factor(count_injury.level_unique$Injury.level, levels = all_levels)

count_injury.level_unique <- complete(count_injury.level_unique, Injury.level, fill=list(n = 0, y = 1, nb_ch = 0))
count_injury.level_unique <- drop_na(count_injury.level_unique)

# ggplot(count_injury.level_unique, aes(x = Injury.level, y = y, color = n))+
#   geom_point() +
#   scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1) +
#   theme_dark()
gg.unique <- ggplot(data=count_injury.level_unique, aes(x = Injury.level, y = n)) +
  geom_bar(stat="identity") +
  #scale_fill_gradient2(low = "red", high = "blue", midpoint = 30) +
  theme_minimal() + 
  coord_flip(ylim = c(0,125)) +
  scale_x_discrete(limits = rev(all_levels))+
  xlab("") +
  theme(plot.margin=unit(c(0.5,-0.8,0.5,0.5), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  ylab("Number of experiments")+
  ggtitle("Unique level reported")

included_unique_doi_range.levels <- included_review_animals %>% 
  dplyr::filter(grepl(" - ", Injury.level))
included_unique_doi_range.levels <- included_unique_doi_range.levels %>% 
  dplyr::filter(!grepl("lumbar",Injury.level))
included_unique_doi_range.levels <- included_unique_doi_range.levels %>% 
  dplyr::filter(!grepl("group",Injury.level))

vec_range <- c()
vec_severity <- c()
vec_mechanism <- c()
vec_species <- c()
for (i in c(1:dim(included_unique_doi_range.levels)[1])){
  start <- str_replace_all(substr(included_unique_doi_range.levels$Injury.level[i], start = 1, stop = 3), " ", "")
  end <- str_replace_all(str_sub(included_unique_doi_range.levels$Injury.level[i], start = -3), " ", "")
  idx_start <- match(start, all_levels)
  idx_end <- match(end, all_levels)
  values <- all_levels[idx_start:idx_end]
  vec_range <- append(vec_range, values)
  vec_severity <- append(vec_severity, rep(included_unique_doi_range.levels$Injury.severity[i], length(values)))
  vec_mechanism <- append(vec_mechanism, rep(included_unique_doi_range.levels$Injury.mechanism[i], length(values)))
  vec_species <- append(vec_species, rep(included_unique_doi_range.levels$Species[i], length(values)))
}

df_range <- data.frame(Injury.level  = vec_range,
                       Injury.severity = vec_severity,
                       Injury.mechanism = vec_mechanism,
                       Species = vec_species)

count_injury.level_range <- vec_range %>%
  data.frame(Injury.level = .) %>%
  dplyr::count(Injury.level)
count_injury.level_range$Injury.level <- factor(count_injury.level_range$Injury.level, levels = all_levels)
count_injury.level_range <- complete(count_injury.level_range, Injury.level, fill=list(n = 0))

gg.range <- ggplot(data=count_injury.level_range, aes(x = Injury.level, y = n)) +
  geom_bar(stat="identity") +
  #scale_fill_gradient2(low = "red", high = "blue", midpoint = 30) +
  theme_minimal() +
  coord_flip(ylim = c(max(count_injury.level_range$n)+5,0)) +
  scale_x_discrete(limits = rev(all_levels)) + scale_y_reverse(limits=c(max(count_injury.level_range$n)+5,0)) + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin=unit(c(0.5,-0.8,0.5,0.5), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  xlab("Injury level(s)")+
  ylab("Number of experiments")+
  ggtitle("Ranges reported")

levels.overall <- grid.arrange(gg.range,
                               gg.unique,
                               widths=c(0.4,0.4),
                               ncol=2,
                               top = textGrob("Number of experiments per level",gp=gpar(fontsize=14)))

################################################################################

for (doi in as.vector(included_review_animals[['DOI.or.PMID']])){
  if (!(doi %in% as.vector(included_unique_doi_no.age[['DOI.or.PMID']]))){
    if (!(doi %in% included_unique_doi_adult$DOI.or.PMID)){
      if (!(doi %in% included_unique_doi_age.range$DOI.or.PMID)){
        if (!(doi %in% included_unique_doi_age.mean$DOI.or.PMID)){
          print(doi)
        }
      }
    }
  }
}

################################################################################
# Figure assessments animals

df_ass_animals <- read.csv('./data/animals_assessments_long.csv')

df_ass_animals <- df_ass_animals %>% 
  mutate(groups=
           case_when(
             Assessment %in% c("BBB",
                               "BMS",
                               "beam walk test",
                               "gait analysis",
                               "grid walking test",
                               "ladder walk test",
                               "inclined plane test",
                               "Tarlov scale",
                               "footprint analysis",
                               "swimming",
                               "locomotor (other)") ~ "locomotion", 
             Assessment %in% c("grip strength",
                              "rearing",
                              "reaching or retrieval",
                              "grooming") ~ "forelimb\nfunction",
             Assessment %in% c("mechanical reactivity",
                               "thermal reactivity",
                               "pain",
                               "toe spread test",
                               "other reflexes") ~ "sensory\nfunction", 
             Assessment %in% c("motor evoked potentials",
                               "somatosensory evoked potentials",
                               "spinal cord evoked potentials",
                               "electrophysiology (other)") ~ 'electro-\nphysiology',
             Assessment %in% c("Gale scale",
                               "composite scores",
                               "hindfoot bar grab test",
                               "spinal cord blood flow",
                               "trunk muscle activity",
                               "micturition") ~ 'other\nfunctional'
             ))

table(df_ass_animals$groups, useNA = 'always')

levels(factor(df_ass_animals$Assessment))

table(table(df_ass_animals$DOI.or.PMID))/sum(table(table(df_ass_animals$DOI.or.PMID)))

out <- df_ass_animals %>% 
  group_by(groups, Assessment) %>% 
  summarise(freq = n(), .groups = 'drop') %>% 
  mutate(Assessment = fct_reorder(Assessment, freq, .desc = F))

temp <- df_ass_animals %>% 
  group_by(DOI.or.PMID, groups) %>% 
  summarise(freq = n(), .groups = 'drop')
n_occur <- data.frame(table(temp$DOI.or.PMID))
dim(n_occur[n_occur$Freq > 1,])

levels(factor(out$Assessment))

out$groups <- factor(out$groups, levels = c("locomotion",
                                            "forelimb\nfunction",
                                            "sensory\nfunction",
                                            "electro-\nphysiology",
                                            'other\nfunctional'))

# p_assessment_v1 <-ggplot(data=out, aes(x=Assessment, y=freq, fill = groups)) +
#   geom_bar(stat="identity", colour="black") +
#   ylab('Number of experiments') +
#   scale_fill_brewer(palette = "Blues") + 
#   theme_half_open() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   background_grid()+
#   guides(fill = guide_legend(reverse = TRUE)) +
#   geom_text(aes(label=freq), position=position_dodge(width=0.9), vjust=-0.25) +
#   guides(fill=guide_legend(title="Assessment categories"))
# p_assessment_v1

p_assessment_v2 <- ggplot(data=out, aes(x=Assessment, y=freq, fill = groups)) +
  scale_fill_manual(values = c('#00547B', "#98A285", "#D9B061", "#D97847", "#690608", "grey")) +
  geom_bar(stat="identity", colour="black") +
  coord_flip() +
  ylab('Number of experiments') +
  theme_half_open() +
  background_grid()+
  facet_grid(groups ~ ., scales = "free_y", space = "free") +
  geom_text(aes(label=freq), position=position_dodge(width=1), hjust = -0.25,
            size=4, family = 'Avenir') +
  theme(axis.title.x = element_text(size = 14, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 16, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=12, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=12, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=14, family = 'Avenir',face='bold'),
        legend.text = element_text(color="black", size=12, family = 'Avenir'),
        legend.title = element_text(color="black", size=14, family = 'Avenir'),
        text = element_text(color="black", size=14, family = 'Avenir'),
        legend.position = "none")
p_assessment_v2

ggsave("./figures/Fig2/plot-assessments.png", p_assessment_v2,
       width = 12, height=10, dpi=300, units = "in",
       bg = "white")
ggsave(file="./figures/Fig2/plot-assessments.svg", plot=p_assessment_v2, width=10, height=8, dpi=300, units = "in")

################################################################################
## Drugs studied in animals

print(paste('Number of drugs or combinations tested:', dim(table(included_review_animals$Drug.name.harmonized))))

count_drugs_comb <- dplyr::filter(included_review_animals, grepl('+', Drug.name.harmonized, fixed = TRUE)) %>% dplyr::count(Drug.name.harmonized)
print(paste('Number of drug combination tested:', dim(count_drugs_comb)[1]))

table_sorted <- table(included_review_animals$Drug.name.harmonized)[order(table(included_review_animals$Drug.name.harmonized), decreasing = TRUE)]
sum(table(included_review_animals$Drug.name.harmonized))

dim(table_sorted[table_sorted > 4])
table_sorted_5 <- table_sorted[table_sorted > 4]
t(table_sorted_5)

plot_drugs_animals <- ggplot(data = as.data.frame(table_sorted_5), aes_string(x = 'Var1', y = 'Freq')) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.9),
            vjust = 0.5, hjust = -0.5, size = 6) +
  labs(x = 'Drug tested',
       y = 'Number of experiments') +
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) +
  theme_half_open() +
  background_grid() + coord_flip()
plot_drugs_animals

included_review_animals_harmomixed <- harmonised_mixed_effects(included_review_animals)

df_plot_results <- included_review_animals_harmomixed %>% 
  group_by(Drug.name.harmonized,
           Drug.effect.on.functional.assessment) %>%
  count() %>%
  filter(Drug.name.harmonized %in% as.data.frame(table_sorted_5)[['Var1']])

df_plot_results$Drug.name.harmonized <- factor(df_plot_results$Drug.name.harmonized, 
                                               levels = levels(factor(as.data.frame(table_sorted_5)[['Var1']])))

# Stacked
ggplot(df_plot_results, aes(fill=Drug.effect.on.functional.assessment, y=n, x=Drug.name.harmonized)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = 'Drug tested',
       y = 'Number of experiments') +
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) +
  theme_half_open() +
  background_grid() + coord_flip() +
  scale_fill_manual(values = c("#E31A1C","#FB9A99", 
                               "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3")) +
  guides(fill=guide_legend(title="Drug effect reported"))

test <- included_review_animals_harmomixed %>% 
  group_by(Drug.name.harmonized,
           Drug.effect.on.functional.assessment) %>%
  count() %>% 
  filter(Drug.effect.on.functional.assessment == 'mixed')
sum(test['freq'])

## Drugs studied in animals - results reported balloon plot

df_balloon <- df_plot_results %>%
  #select(Drug.effect.on.functional.assessment, Drug.name.harmonized) %>%
  #group_by(Drug.name.harmonized,
  #         Drug.effect.on.functional.assessment) %>%
  #count() %>%
  spread(Drug.effect.on.functional.assessment, n)

df_balloon <- cbind(Drug.name.harmonized = df_balloon[, 1],
                    df_balloon[, -1]/rowSums(df_balloon[, -1], na.rm=T))
df_balloon <- df_balloon %>% replace(is.na(.), 0)
df_balloon <- df_balloon %>%
  remove_rownames %>%
  column_to_rownames(var="Drug.name.harmonized")
df_balloon <- df_balloon[,order(colSums(-df_balloon,na.rm=TRUE))]

my_cols <- c("#00547B", "#ffffff")
balloon_plot_all_effects <- ggballoonplot(df_balloon, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
balloon_plot_all_effects

included_review_copy_only.mixed <- included_review_animals[grepl('mixed', included_review_animals$Drug.effect.on.functional.assessment),]
included_review_copy_only.mixed$Drug.name.harmonized <- factor(included_review_copy_only.mixed$Drug.name.harmonized,
                                                levels = levels(factor(as.data.frame(table_sorted_5)[['Var1']])))
df_balloon_mixed <- included_review_copy_only.mixed %>%
  select(Drug.effect.on.functional.assessment, Drug.name.harmonized) %>%
  filter(Drug.name.harmonized %in% as.data.frame(table_sorted_5)[['Var1']]) %>%
  group_by(Drug.name.harmonized,
           Drug.effect.on.functional.assessment) %>%
  count() %>%
  spread(Drug.effect.on.functional.assessment, n) %>%
  replace(is.na(.), 0) %>%
  remove_rownames %>%
  column_to_rownames(var="Drug.name.harmonized")
df_balloon_mixed <- df_balloon_mixed[,order(colSums(-df_balloon_mixed,na.rm=TRUE))]

balloon_plot_mixed_effects <- ggballoonplot(df_balloon_mixed, fill = "value")+
  scale_fill_gradientn(colors = my_cols)


## Drugs studied in animals - results reported balloon plot
## New balloon plot: now based on geom_point
## to colour points based on the proportion of a certain effect by drug

df_plot_results_points <- df_plot_results %>%
  group_by(Drug.name.harmonized) %>%
  mutate(freq = n / sum(n))

df_plot_results_points$Drug.name.harmonized <- fct_rev(df_plot_results_points$Drug.name.harmonized)
df_plot_results_points$Drug.effect.on.functional.assessment <- factor(
  df_plot_results_points$Drug.effect.on.functional.assessment, 
  levels = c("positive effect", "mixed", "no effect", "no stats", "not reported",
             "negative effect"))

balloon <- ggplot(df_plot_results_points, aes(x = Drug.effect.on.functional.assessment, 
                                   y = Drug.name.harmonized)) +
  geom_point(shape = 21, color = 'black',
             aes(size = n, fill = freq))+
                 #color = freq)) +
  scale_size(range = c(2, 9)) +
  scale_fill_gradientn(colors = c("#F0F921", '#D9B061', "#C86343", '#1c98c4', '#156f9c', "#00547B"),
                       name = 'Frequency of effect\n reported by drug',
                       limits=c(0, 0.8)) +
  # scale_color_gradientn(colors = c("#F0F921", '#D9B061', "#C86343", '#1c98c4', '#156f9c', "#00547B"),
  #                      name = 'Frequency of effect\n reported by drug',
  #                      limits=c(0, 0.8)) +
  labs(y = 'Drug name', x = 'Effect reported',
       size = "Number of experiments") + 
  theme_minimal() +
  annotate("text",
         x = 1:length(table(df_plot_results_points$Drug.effect.on.functional.assessment)),
         y = length(table(df_plot_results_points$Drug.name.harmonized)),
         label = paste0('n = ', aggregate(df_plot_results_points$n, 
                           by=list(Category=df_plot_results_points$Drug.effect.on.functional.assessment), 
                           FUN=sum)$x),
         col = "black",
         vjust = - 1.5,
         family = 'Avenir') +
  annotate("text",
           y = 1:length(table(df_plot_results_points$Drug.name.harmonized)),
           x = length(table(df_plot_results_points$Drug.effect.on.functional.assessment)) +0.35,
           label = paste0('n = ', aggregate(df_plot_results_points$n, 
                             by=list(Category=df_plot_results_points$Drug.name.harmonized), 
                             FUN=sum)$x),
           col = "black",
           vjust = 0.4,
           family = 'Avenir') + 
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir'), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'), 
        legend.text = element_text(size=10, family = 'Avenir'),
        legend.title = element_text(size=12, family = 'Avenir'),
        line = element_line(linewidth = 0.3, linetype = 1),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))

filter(df_plot_results_points, Drug.effect.on.functional.assessment == 'positive effect' & freq >0.50)

ggsave(file="./figures/Figure2B.svg", plot=balloon, width=10, height=9, dpi=300, units = "in")

##################

df_plot_results_points_mixed <- included_review_copy_only.mixed %>%
  select(Drug.effect.on.functional.assessment, Drug.name.harmonized) %>%
  filter(Drug.name.harmonized %in% as.data.frame(table_sorted_5)[['Var1']]) %>%
  droplevels() %>%
  group_by(Drug.name.harmonized,
           Drug.effect.on.functional.assessment) %>%
  count() %>%
  ungroup() %>%
  group_by(Drug.name.harmonized) %>%
  mutate(freq = n / sum(n))

df_plot_results_points_mixed$Drug.name.harmonized <- fct_rev(df_plot_results_points_mixed$Drug.name.harmonized)
df_plot_results_points_mixed$Drug.effect.on.functional.assessment <- factor(
  df_plot_results_points_mixed$Drug.effect.on.functional.assessment, 
  levels = names(df_balloon_mixed))

balloon_mixed <- ggplot(df_plot_results_points_mixed, aes(x = Drug.effect.on.functional.assessment, 
                                              y = Drug.name.harmonized)) +
  geom_point(shape = 21, color = 'black',
             aes(size = n, fill = freq))+
  #color = freq)) +
  scale_size(range = c(2, 9)) +
  scale_fill_gradientn(colors = c("#F0F921", '#D9B061', "#C86343", '#1c98c4', '#156f9c', "#00547B"),
                       name = 'Frequency of effect\n reported by drug',
                       limits=c(0, 0.8)) +
  # scale_color_gradientn(colors = c("#F0F921", '#D9B061', "#C86343", '#1c98c4', '#156f9c', "#00547B"),
  #                      name = 'Frequency of effect\n reported by drug',
  #                      limits=c(0, 0.8)) +
  labs(y = 'Drug name', x = 'Effect reported',
       size = "Number of experiments") + 
  theme_minimal() +
  annotate("text",
           x = 1:length(table(df_plot_results_points_mixed$Drug.effect.on.functional.assessment)),
           y = length(table(df_plot_results_points_mixed$Drug.name.harmonized)),
           label = paste0('n = ', aggregate(df_plot_results_points_mixed$n,
                                            by=list(Category=df_plot_results_points_mixed$Drug.effect.on.functional.assessment),
                                            FUN=sum)$x),
           col = "black",
           vjust = - 1.5,
           family = 'Avenir') +
  annotate("text",
           y = 1:length(table(df_plot_results_points_mixed$Drug.name.harmonized)),
           x = length(table(df_plot_results_points_mixed$Drug.effect.on.functional.assessment)) +0.35,
           label = paste0('n = ', aggregate(df_plot_results_points_mixed$n,
                                            by=list(Category=df_plot_results_points_mixed$Drug.name.harmonized),
                                            FUN=sum)$x),
           col = "black",
           vjust = 0.4,
           family = 'Avenir') +
  theme(axis.title.x = element_text(size = 12, family = 'Avenir', face='bold'),
        plot.title =  element_text(size = 14, family = 'Avenir', face='bold'),
        axis.text.x = element_text(color="black", size=10, family = 'Avenir', angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(color="black", size=10, family = 'Avenir'), 
        axis.title.y  = element_text(color="black", size=12, family = 'Avenir'), 
        legend.text = element_text(size=10, family = 'Avenir'),
        legend.title = element_text(size=12, family = 'Avenir'),
        line = element_line(linewidth = 0.3, linetype = 1),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))
balloon_mixed
ggsave(file="./figures/SuppFigure2/SuppFigure2.svg", plot=balloon_mixed, width=10, height=9, dpi=300, units = "in")

################################################################################
## Funnel chert

# Number of experiments
dim(included_review_animals)[1]

# Number of experiments with MP
included_review_animals_MP <- included_review[included_review$Drug.name.harmonized == 'methylprednisolone',]
dim(included_review_animals_MP)[1]

# Number of experiments with MP in rats
included_review_animals_MP_rats <- included_review_animals_MP[included_review_animals_MP$Species == 'rats',]
dim(included_review_animals_MP_rats)[1]

# Number of experiments with MP in rats and thoracic injuries
remove_if = c("C", "L")
included_review_animals_MP_rats_thoracic <- included_review_animals_MP_rats %>%
  filter(str_detect(Injury.level, paste(remove_if, collapse = "|"), negate = TRUE))
dim(included_review_animals_MP_rats_thoracic)[1]

# Number of experiments with MP in rats and thoracic injuries and BBB
included_review_animals_MP_rats_thoracic_BBB <- included_review_animals_MP_rats_thoracic %>%
  filter(str_detect(Name.type.of.asessement, 'BBB'))
dim(included_review_animals_MP_rats_thoracic_BBB)[1]

# Number of experiments with MP in rats and thoracic injuries and BBB at day 28
included_review_animals_MP_rats_thoracic_BBB_28d <- included_review_animals_MP_rats_thoracic_BBB[
  included_review_animals_MP_rats_thoracic_BBB$Assessment.on.day.28..yes.no. == 'Yes',]
dim(included_review_animals_MP_rats_thoracic_BBB_28d)[1]

# Number of experiments with MP in rats and thoracic injuries and BBB at day 28 with MP given in first 24h
included_review_animals_MP_rats_thoracic_BBB_28d_24h <- included_review_animals_MP_rats_thoracic_BBB_28d[
  !grepl("reported", included_review_animals_MP_rats_thoracic_BBB_28d$Time..minutes.pre.injury..minutes.post.injury.),]
dim(included_review_animals_MP_rats_thoracic_BBB_28d_24h)[1]

# Number of experiments with MP in rats and thoracic injuries and BBB at day 28 with MP given in first 24h, with 30 mg/kg
included_review_animals_MP_rats_thoracic_BBB_28d_24h_30mg <- included_review_animals_MP_rats_thoracic_BBB_28d_24h[
  included_review_animals_MP_rats_thoracic_BBB_28d_24h$Dose..absolute.dose.or.mg.kg. == '30 mg/kg',]
dim(included_review_animals_MP_rats_thoracic_BBB_28d_24h_30mg)[1]

#funnel chart :MP, level of injury 1, injury severity 2, timing when drug given [>24 versus <24] 4, dose 5, timing of assessment 3) (Lucie)


included_review_copy <- included_review
included_review_copy <- harmonised_mixed_effects(included_review_copy)

harmonised_mixed_effects <- function(data){
  data$Drug.effect.on.functional.assessment <- as.character(data$Drug.effect.on.functional.assessment)
  data$Drug.effect.on.functional.assessment[grepl('mixed', data$Drug.effect.on.functional.assessment)] <- 'mixed'
  data$Drug.effect.on.functional.assessment[data$Drug.effect.on.functional.assessment == ''] <- NA
  data$Drug.effect.on.functional.assessment[data$Drug.effect.on.functional.assessment == '-'] <- NA
  return (data)
}

included_review_positive <- included_review_copy %>%
  filter(Drug.effect.on.functional.assessment %in% c('positive effect', 'mixed', 'negative effect'))
length(unique(included_review_positive$Drug.name.harmonized))
