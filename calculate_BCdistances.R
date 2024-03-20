library(tidyverse)
setwd("~/Desktop/GF_16s_1")
getwd()
path <- "/Users/charlotte.kenneally/Desktop/GF_16s_1"


### reading in a tibble that contains my sample names, the day, and the OTUs
### i am going to filter out the in between days to only keep days at the beginning and at the end since we are measuring diversity distances 
### also filtering out OTUs that are 0
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  select(Group, day, everything()) ### this line is just reorganizing table to the order of columns i specify (group, day, everything else)


## copy and paste from up there ^^ so i dont delete it, adding days I want in this
## we want the first week, and the last 2 weeks
days_wanted <- c(1:7, 49:56) ## specify what days, run this

read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  select(Group, day, everything()) %>%
  filter(day %in% days_wanted) ## filter by column name "day" for what is within "days_wanted" you just created 

## next pivot data frame longer instead of wider so we can count the total number of sequences in each sample 
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(total = sum(value, na.rm = TRUE)) %>%
  arrange(total) %>%
  print(n = 20) 
## this now displays a tibble that is 126 rows by two columns, titled group and total in ascending order 
## only samples from days 1-7 and days 49 and 56 are in this
## future me: if the first row shows 1D-3-D56 with total OTU count as 987 then you did it correctly 
## purpose: see the lowest number of OTU for a sample so we know where to cut off for calculating distances 

## after reading thru the top 20, pick a threshold where there is a good sequence mass beginning 
## I think it could be just 2179 because that's already a high number of sequences (not compared to the rest of the group but it is high)
## another thought to picking this threshold could be where the samples start to "stabilize" and show less difference in numbers
## in this case that would be about 10,000 because there are more than one sample with around 10k and after this increase in sequence numbers is a slower progression 
## mutate the table by cutting out sequences lower than my designated threshold (2000)
## in this case I only lose 1 sample with this cut off -- can always come back and cut more if wanted 

## take the same lines and change the "summarize" into a mutate 
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 2000) 

## now we want to remove any OTUs that are 0
## this isn't causing any problem, it just makes the data set larger and will be easier if they are removed 
## add a new "group_by" and "name" so we group it by the OTU column 
## then summarize together by the value 
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 2000) %>%
  group_by(name) %>%
  summarize(total = sum(value)) %>% filter(total == 0)
## this now showed me there are 662 OTU values that do not have any counts 

## now lets take them out 
## change summarize to mutate again 
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 2000) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>% 
  filter(total != 0)
## now we removed any OTUs that had absolutely no sequences 

## clean up the data set now 
read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 2000) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Group)
## returns a wider data frame that has 125 samples total with only OTUs that are actually used in at least one of the samples
## removed the "total" column because we do not need to see that anymore 
## 125 samples is GOOD - double check with the "days wanted" to check the total samples there
## mine had 126, and remember we removed one sample when cleaning for OTU total counts, so 125 is where it should be


### now need to convert the tibble to a data frame
## also save the pipeline as "dt_start_end_BC"
dt_start_end_BC <- read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_wanted) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 2000) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Group) %>%
  as.data.frame()

## make the rownames the same 
rownames(dt_start_end_BC) <- dt_start_end_BC$Group
## remove the group Column
dt_start_end_BC <- dt_start_end_BC[, -1]
## turn the data frame into a matrix because vegan needs a matrix as input for its functions 
dt_start_end_BC <- as.matrix(dt_start_end_BC)


### NOW WE CALCULATE DISTANCES USING VEGAN 
library(vegan)

# make object called dist (for distances) and using vegdist specify you want bray curtis distances 
set.seed(2427)
dist <- vegdist(dt_start_end_BC, method = "bray")
nmds <- metaMDS(dist)

scores(nmds)
## ^^ this is being stored as a matrix
## in order to plot this with ggplot, we need to make it a tibble or a data frame 
scores(nmds) %>%
  as.tibble(rownames = "Group") %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point()
## points below 0 are from the earlier time points, points above 0 are from the later time points 

## alt option 
set.seed(2427)
dist <- avgdist(dt_start_end_BC, dmethod = "bray", sample = 2000)

set.seed(2427)
nmds <- metaMDS(dist)


### pcoa plots
library(glue)
pcoa <- cmdscale(dist, eig = TRUE, add = TRUE)
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
pretty_pe <- round(percent_explained[1:2], digits = 1)
labs <- c(glue("PCO 1 ({pretty_pe[1]}%)" ),
          glue("PCO 2 ({pretty_pe[2]}%)"))

positions %>% 
  as_tibble(rownames = "Samples") %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color = "Samples"))+
  geom_point()+
  labs(x=labs[1], y=labs[2])+
  theme_classic()

### using the entire data set now 
dt_ps3 %>%
  select(Sample, starts_with("OTU")) %>%
  mutate(day = str_replace(Sample, ".*D", "")) %>%
  filter(day %in% all_days) %>%
  select(Sample, day, everything()) %>%
  pivot_longer(-Sample)

## decide where to cut it now based on sequence size 
all_days <- c(1:56)
round1 <- c("1A", "1B", "1C", "1D")
dt_ps3 %>%
  select(Sample, starts_with("1")) %>%
  mutate(round1 = str_replace(round, "1", "")) %>%
  filter(round1 %in% round1)

days_1_7 <- c(1:7)

read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  filter(day %in% days_1_7) %>%
  select(Group, day, everything()) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(total = sum(value)) %>% 
  arrange(total) %>%
  print(n = 20)
## mainly the D0s with little to no sequences present 
## cut it off at 1600 
## this loses 15 samples 

days_1_7 <- c(1:7)

## create dataframe of the first week of samples and don't include pooled 2 because we only have few samples for pooled 2 in this run 
week1 <- read_tsv("OTU_day_name_BC.txt") %>%
  #select(Group, starts_with("OTU")) %>%
  mutate(Day = str_replace(Day, ".*D", "")) %>%
  filter(Day %in% days_1_7) %>%
  select(Group, Day, everything()) %>%
  filter(donor != "Pooled2")
  
 # select(-Day) %>% # remove the column day from the tibble 
 # pivot_longer(-Group) %>%
#  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 10000) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Group) %>%
  as.data.frame()


##rownames(week1) <- week1$Group
## week1 <- week1[,-1]
##week1 <- as.matrix(week1)

## run these since i did not convert to data frame in previous step 
## line up group names with the OTU matrix 
OTUw1 <- week1 %>% select(starts_with("OTU")) %>% as.matrix()
rownames(OTUw1) <- week1$Group

## keep the set seed 
set.seed(2427)
##dist_week1 <- vegdist(week1, method = "bray")
## calculate distances 
nmds_week1 <- metaMDS(OTUw1, distance = 'bray')
class(nmds_week1)
class(week1)
class(BC_doc_6000)

## plot by colors day and shape donors 
cbind.data.frame(scores(nmds_week1, 'sites'), "Day" = week1$Day, "Donor"=week1$donor) %>%
  as_tibble(rownames = "Group") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Day, shape=Donor)) +
  geom_point(alpha = 6) +
  stat_ellipse() +
  labs(title = "All stool days 1-7")

wk1_plots <- cbind.data.frame(scores(nmds_week1, 'sites'), "Day" = week1$Day, "Donor"=week1$donor)



## plot by donor only 
cbind.data.frame(scores(nmds_week1, 'sites'), "Day" = week1$Day, "Donor"=week1$donor) %>% 
  as_tibble(rownames = "Group") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Donor)) +
  geom_point(alpha = 6) +
  stat_ellipse() +
  theme_classic()+
  labs(title = "Week 1 Bray Curtis distribution")

## dont filter out by days here yet, just cut samples down based on otu level (6000) 
BC_doc_6000 <- read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  mutate(day = str_replace(Group, ".*D", "")) %>%
  select(Group, day, everything()) %>%
  select(-day) %>% # remove the column day from the tibble 
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  summarize(total = sum(value)) %>% 
  arrange(total) %>%
  print(n = 20)

### START HERE WHEN OPENING BACK UP 
## cut 6000 and lower 
BC_doc_6000 <- read_tsv("OTU_day_name_BC.txt") %>%
  select(Group, starts_with("OTU")) %>%
  ##select(Group, everything()) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>% 
  filter(total > 6000) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() #%>%
  #select(-total) %>%
  #pivot_wider(Group) %>%
 # as.data.frame()

class(BC_doc_6000)
## ^^ this is now my BC doc that has every sample which had 6000 reads or HIGHER
### use this object now for plots and filter which samples i want through that for each plot 
try1 <- BC_doc_6000 %>% select(starts_with("OTU")) %>% as.matrix()
rownames(try1) <- BC_doc_6000$Group
set.seed(2427)
## make the rownames the same 
#rownames(BC_doc_6000) <- BC_doc_6000$Group
## remove the group Column
#BC_doc_6000 <- BC_doc_6000[, -1]
## turn the data frame into a matrix because vegan needs a matrix as input for its functions 
#BC_matrix_6000 <- as.matrix(BC_doc_6000)
## plot the last week for all stools
##day56_allstools <- BC_doc_6000 %>% select(starts_with("OTU")) %>% as.matrix()
## calculate nmds scores 
nmds <- metaMDS(try1, distance = 'bray')

meta_nmds <- cbind.data.frame(scores(nmds, 'sites'), "Day" = BC_doc_6000$Day, "Donor"=BC_doc_6000$donor)


class(week1)
class(BC_doc_6000)
### try next time is says error memory exhausted 
Sys.setenv('R_MAX_VSIZE'=32000000000)

#set.seed(2427)
##dist_week1 <- vegdist(week1, method = "bray")
#nmds_day56 <- metaMDS(day56_allstools, distance = 'bray')



