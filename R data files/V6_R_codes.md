### Analysis of the microbiome data in R
##### Files to upload in R
###### A. ASV count table
###### B. Taxonomy metadata 
###### Metatable 

#### Load the following libraries
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2","readxl","openxlsx","vegan","ggplot2","dplyr","tidyverse","phyloseq","gridExtra","grid","ggpubr","agricolae","microbiome","MicrobeR","microbiomeSeq","apeglm"))
library(DESeq2)
library(readxl)
library(openxlsx)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(gridExtra)
library(grid)
library(ggpubr)
library(agricolae)
library(microbiome)
#library(MicrobeR)
library(pheatmap)
library("apeglm")
library(ranacapa)
library(devtools)
library(phylosmith)
library(VennDiagram)
library(UpSetR)
library(coin)
library(MicrobiotaProcess)

```
#### Load all the ASV counttables and metatable, taxonomy table and combine all into one big table in R
```
setwd("~/Desktop/microbiome/VIMS16S/mvpics/tsv/new/")
```
###### ASV Counttable
```
abund_table <-read.xlsx("~/Documents/sequence/All_data.xlsx", sheet="table", rowNames = TRUE)
abund_table<- as.data.frame(t(abund_table))
abund_table1 <- cbind(SampleID = rownames(abund_table), abund_table)
rownames(abund_table1) <- NULL
abun_count <-gather(abund_table1, FeatureID, Value, ASV_1:ASV_4155)
#add the rownames as a proper column
```
###### Metatable
```
meta_table <-read.xlsx("~/Documents/sequence/All_data.xlsx", sheet="meta_table", rowNames = TRUE)
meta_table1 <- cbind(SampleID = rownames(meta_table ), meta_table)
rownames(meta_table1 ) <- NULL
meta_table1$Treatment <- as.character(meta_table1$Treatment)
meta_table1$Treatment <- factor(meta_table1$Treatment, levels=unique(meta_table1$Treatment))
meta_table1$ProjectName <- as.character(meta_table1$ProjectName)
meta_table1$Tanksystem <- factor(meta_table1$ProjectName, levels=unique(meta_table1$ProjectName))
```
###### Read the taxonomy table
```
taxa1<-read_excel("~/Documents/sequence/All_data.xlsx", sheet = "metadata")
tax_table1 <-do.call(rbind, strsplit(as.character(taxa1$Taxon),";"))
colnames(tax_table1) <-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax_table1)<-taxa1$FeatureID
taxa_A <- as.data.frame(tax_table1)
taxa_A1  <- cbind(FeatureID = rownames(taxa_A ), taxa_A)
rownames(taxa_A1) <- NULL
```
#### Convert to phyloseq
```
Abund_table<- as.data.frame(t(abund_table))
A2 <-as.matrix(Abund_table)
# First we need to create a phyloseq object using our un-transformed count table
count_tab_pyA <- otu_table(A2, taxa_are_rows=T)
# Read in the taxonomy 
taxa_A2 <-taxa_A1 %>% remove_rownames %>% column_to_rownames(var="FeatureID")
taxa_A2 <-as.matrix(taxa_A2)
tax_tab_phy <- tax_table(taxa_A2)
# Read in the metatable
meta_table2 <-Meta_T %>% remove_rownames %>% column_to_rownames(var="SampleID")
sample_info_tab_phy <- sample_data(meta_table2)
All_S4 <-phyloseq(count_tab_pyA, tax_tab_phy, sample_info_tab_phy)
All_S4
```
#### Remove any taxa associated with chloroplast and mitochondria
```
justbacteria <- All_S4_prune %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   
      Family  != "mitochondria" &           
      Class   != "Chloroplast" &
      Order   != "mitochondria" &
      Order   != "Chloroplast"
  )
```

#### Join all the data
```
sample_sum_df <- data.frame(sum = sample_sums(justbacteria ))
sum(sample_sum_df[, 'sum'])
Sample_reads <- cbind(SampleID = rownames(sample_sum_df), sample_sum_df)
rownames(Sample_reads) <- NULL
#Join all the data
Meta_T <-full_join(meta_table1, Sample_reads)
M1count<-full_join(taxa_A1, abun_count)
combine_b <-full_join(M1count,Meta_T)
All_samples <- combine_b%>% group_by(SampleID,Order,Treatment)%>%summarise(avg_abundance = mean(sum))
```
#### Basic Statistics
```
newdata <-combine_b %>% group_by(SampleID,FeatureID) %>% summarise(Value=mean(Value))
newdata1 <-spread(newdata, SampleID, Value) ##spread the data into a wide format
x<-dim(newdata1)[2] # number of columns
newdata1[2:x] <- lapply(newdata1[2:x], function(x) as.numeric(as.character(x)))
```

#### Get basic stats and generate a table called "dataset_info"
```
seq_total<-apply(newdata1[2:x],2,sum) #number of sequences per sample
ASV_count<-colSums(newdata1[2:x]>0) # OTUs per sample
ASV_single<-colSums(newdata1[2:x]==1) # OTUs with only 1 seq
ASV_double<-colSums(newdata1[2:x]==2) # OTUs that have only 2 seqs
ASV_true<-colSums(newdata1[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,ASV_count)
#add the rownames as a proper column
dataset_info <- cbind(SampleID = rownames(dataset_info), dataset_info)
dataset_info
rownames(dataset_info) <- NULL
#combine the metadata with dataset_info
data_meta <- meta_table1 [c(1,2,4)]
basic_stat <- full_join(data_meta, dataset_info)
basic_stat
#write.csv(basic_stat,"basic_stat.csv")
``` 
#### Subset by larvae
```
phy <- subset_samples(justbacteria, Sym=="Y")
phy <- subset_samples(justbacteria, All=="IN")
Table1 <- otu_table(phy)
SD <- sample_data(phy) %>%
    data.frame() %>%
    select("Treatment","Hatchery", "Year","Month","Season","season","ProjectName", "Num","Type", "Name","All","TankLocation","Location","Env") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
sample_sum_df <- data.frame(sum = sample_sums(ps2 ))
sum(sample_sum_df[, 'sum'])
```
#### Data Filtering
##### Step1. Evaluate Amplicon Sequence Variants (ASV) summary statistics
####### Remove taxa with zero counts 
```
summary(taxa_sums(justbacteria))
phy <- prune_taxa(taxa_sums(justbacteria) > 0, justbacteria)
```
#### Step2. Factor reordering and renaming
##### Reorder Trials
```
levels(sample_data(phy)$Trial)
sample_data(phy)$Trial <- factor(sample_data(phy)$Trial, levels = c("T6","T8"))
levels(sample_data(phy)$Trial)

# Relabel Water Type
sample_data(phy)$Water <- factor(sample_data(phy)$Water, labels = c("T6NUV","T6UV","T8NUV","T8UV"))
levels(sample_data(phy)$Water)

# Reorder Hatchery labels
levels(sample_data(phy)$Hatchery)
sample_data(phy)$Hatchery <- factor(sample_data(phy)$Hatchery, levels = c("MOOK", "RWU"))
levels(sample_data(phy)$Hatchery)

# Create a new data frame of the sorted row sums, a column of sorted values from 1 to the total number of individuals/counts for each ASV and a categorical variable stating these are all ASVs.
readsumsdf <- data.frame(nreads = sort(taxa_sums(phy), decreasing = TRUE), 
                        sorted = 1:ntaxa(phy), type = "ASVs")
# Make a data frame with a column for the read counts of each sample for histogram production
sample_sum_df <- data.frame(sum = sample_sums(phy))
```

#### Step3. Make plots for number of reads
###### Generates a bar plot with # of reads (y-axis) for each taxa. Sorted from most to least abundant
```
p.reads = ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle("ASV Assessment") +
  scale_y_log10() +
  facet_wrap(~type, scales = "free") +
  ylab("# of Sequences")
 ``` 
###### Generates a second bar plot with # of reads (y-axis) per sample at various read depths. Sorted from most to least
```
p.reads.hist <- ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "firebrick3", binwidth = 150) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("# of Samples")
# Final plot, side-by-side
grid.arrange(p.reads, p.reads.hist, ncol = 2)
```

#### Step 4: Detect and remove outlier samples
##### Format a data table to combine sample summary data with sample variable data
```
ss <- sample_sums(phy)
sd <- as.data.frame(sample_data(phy))
ss.df <- merge(sd, data.frame("ASV" = ss), by ="row.names")
# Plot the data by the treatment variable
y = 10000 # Set a threshold for the minimum number of acceptable reads. 
x = "SampleName" # Set the x-axis variable you want to examine
label = "Row.names" # This is the label you want to overlay on the points
p.ss.boxplot <- ggplot(ss.df, aes_string(x, y = "ASV", color = "Treatment")) + 
  geom_boxplot(outlier.colour="NA", position = position_dodge(width = 0.8)) +
  geom_jitter(size = 2, alpha = 0.6) +
  scale_y_log10() +
  facet_wrap(~Water) +
  geom_hline(yintercept = y, lty = 2) +
  geom_text(aes_string(label = label), size = 3, nudge_y = 0.05, nudge_x = 0.05)
p.ss.boxplot
```
###### Remove samples that are outliers
```
nsamples(phy)
ps1 <- phy2 %>% subset_samples(SampleName != "VIMS4_S4A_3" & SampleName != "UVL_S4A_3"& SampleName != "VIMS1_CONA_1")
nsamples(ps1)
```

#### Step 5: Taxon cleaning
###### The following R chunk removes taxa not-typically part of a bacterial microbiome analysis.
```
get_taxa_unique(phy, "Kingdom")
get_taxa_unique(phy, "Class")
phy # Check the number of taxa prior to removal
ntaxa(phy)
ps2 <- phy %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    Family  != "mitochondria" &
    Class   != "Chloroplast" &
    Phylum != "Cyanobacteria/Chloroplast"
  )
ps2 # Confirm that the taxa were removed
get_taxa_unique(ps2, "Kingdom")
get_taxa_unique(ps2, "Class")
get_taxa_unique(ps2, "Phylum")
```
#### Step 6: Data Subsetting
###### All samples
```
ntaxa(ps2)
ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)
ntaxa(ps2)
# Trial1
ps2.T1 <- subset_samples(ps2, ProjectName == "T1")
any(taxa_sums(ps2.T1) == 0) # In this case it is TRUE, so remove the zero's
ps2.T1 <- prune_taxa(taxa_sums(ps2.T1) > 0, ps2.T1)
any(taxa_sums(ps2.T1) == 0) # It should now be false
# Trial2
ps2.T2 <- subset_samples(ps2, ProjectName == "T2")
any(taxa_sums(ps2.T2) == 0) # In this case it is TRUE, so remove the zero's
ps2.T2 <- prune_taxa(taxa_sums(ps2.T2) > 0, ps2.T2)
any(taxa_sums(ps2.T2) == 0) # It should now be false
# Trial3
ps2.T3 <- subset_samples(ps2, ProjectName == "T3")
any(taxa_sums(ps2.T3) == 0) # In this case it is TRUE, so remove the zero's
ps2.T3 <- prune_taxa(taxa_sums(ps2.T3) > 0, ps2.T3)
any(taxa_sums(ps2.T3) == 0) # It should now be false
# Trial4
ps2.T4 <- subset_samples(ps2, ProjectName == "T4")
any(taxa_sums(ps2.T4) == 0) # In this case it is TRUE, so remove the zero's
ps2.T4 <- prune_taxa(taxa_sums(ps2.T4) > 0, ps2.T4)
any(taxa_sums(ps2.T4) == 0) # It should now be false
# Trial5
ps2.T5 <- subset_samples(ps2, ProjectName == "T5")
any(taxa_sums(ps2.T5) == 0) # In this case it is TRUE, so remove the zero's
ps2.T5 <- prune_taxa(taxa_sums(ps2.T5) > 0, ps2.T5)
any(taxa_sums(ps2.T5) == 0) # It should now be false
# Trial6
ps2.T6 <- subset_samples(ps2, Name == "T6")
any(taxa_sums(ps2.T6) == 0) # In this case it is TRUE, so remove the zero's
ps2.T6 <- prune_taxa(taxa_sums(ps2.T6) > 0, ps2.T6)
any(taxa_sums(ps2.T6) == 0) # It should now be false
# Trial7
ps2.T7 <- subset_samples(ps2, Name == "T7")
any(taxa_sums(ps2.T7) == 0) # In this case it is TRUE, so remove the zero's
ps2.T7 <- prune_taxa(taxa_sums(ps2.T7) > 0, ps2.T7)
any(taxa_sums(ps2.T7) == 0) # It should now be false
# Trial8
ps2.T8 <- subset_samples(ps2, Name == "T8NUV")
any(taxa_sums(ps2.T8) == 0) # In this case it is TRUE, so remove the zero's
ps2.T8 <- prune_taxa(taxa_sums(ps2.T8) > 0, ps2.T8)
any(taxa_sums(ps2.T8) == 0) # It should now be false
```
###### Taxa Prevalence estimation
```
# Calculate ASVs  prevalence across the data set
prevdf <- apply(X = otu_table(ps2),MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to prevdf
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps2), tax_table(ps2))
```
######Prevalence plot
```
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))
p.prevdf1 <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps2),color=Family)) +
geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
geom_point(size = 3, alpha = 0.7) +
scale_x_log10() +
xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
facet_wrap(~Phylum) +
theme(legend.position="none") +
 ggtitle("Phylum Prevalence in All Samples\nColored by Family")
p.prevdf1
```	
	
#### The R chunk below implements several commonly used transformations in microbiome research and plots their results. Similar to outlier removal and prevalence filtering, your choice should be justified, tested and documented.
###### Transform to Relative abundances
```
ps2.ra <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
```
###### Log transformation moves to a more normal distribution
```
ps2.log <- transform_sample_counts(ps2, function(x) log(1 + x))
```
###### Histograms of the non-transformed data vs. the transformed data can address the shift to normality
```
p.no.log <- qplot(rowSums(otu_table(ps2))) + ggtitle("Raw Counts") +
  theme_bw() +
  xlab("Row Sum") +
  ylab("# of Samples")
p.log <- qplot(log10(rowSums(otu_table(ps2)))) +
  ggtitle("log10 transformed counts") +
  theme_bw() +
  xlab("Row Sum") +
  ylab("# of Samples")
ggarrange(p.no.log, p.log, ncol = 2, labels = c("A)", "B)"))
```
#### Rarefaction plot
####### We can plot the rarefaction curve for the observed ASVs in the entire data set. This is a way to check how has the richness captured in the sequencing effort.
```
Raf <-ggrare(ps2, step = 50, color="Trial",  se = TRUE)
```
### Alpha diversity
###### Alpha diversity is a standard tool researchers can use to calculate the number of bacterial taxa present in a study or study group and the relationships between relative abundance and how evenly taxa are distributed. These are classic representations of species number and diversity in a study which provide useful summary information about the numbers and relative abundances of bacterial taxa within your study.
```
SD <- sample_data(ps2) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","Hatchery","Location","season","Year","Month","Season","TankLocation") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
All_S4_prune_rarefy <- rarefy_even_depth(ps2, rngseed= 81, sample.size = 10000)
#Table1 <- otu_table(ps2)
adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(All_S4_prune_rarefy, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(All_S4_prune_rarefy, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(All_S4_prune_rarefy, measures = "Simpson"),
  "Chao1"= phyloseq::estimate_richness(All_S4_prune_rarefy, measures = "Chao1"),
  "Status" = phyloseq::sample_data(All_S4_prune_rarefy)$Treatment)
head(adiv)
adiv$ProjectName <- SD$ProjectName
adiv$Treatment <- SD$Treatment
adiv$Name <- SD$Name
adiv$Num <- SD$Num
adiv$Hatchery <- SD$Hatchery
adiv$Season <- SD$Season
adiv$Month <- SD$Month
adiv$Year <- SD$Year
adiv$TankLocation <- SD$TankLocation
adiv$Status <-as.character(adiv$Status)
newSTorder = c("Control",  "S4")
sample_types <- c("Control", "S4")
sample_labels <- c("C", "S4")
adiv$Hatchery <- factor(adiv$Hatchery, levels=c("VIMS","Mat","MOOK","RWU"))
adiv$ProjectName <- factor(adiv$ProjectName, levels=c("T1","T2","T3","T4","T5","T6","T7","T8"))
```
##### Simpson diversity
```
Div_SIM <- ggplot(adiv, aes(as.factor(Name), Simpson)) +
  facet_grid(.~Hatchery, scales="free",space="free_x",switch="y")+
  geom_boxplot(aes(fill=factor(Treatment)),width=0.5)+
  #geom_signif(comparisons = list(c("A", "B"), c("B","C")),
              #map_signif_level=TRUE)+
  scale_fill_manual(values=c("blue", "orange","blue", "orange")) +
  labs(y="Simpson's Diversity", x="Trial", fill="Treatment")+theme_bw()+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(strip.background = element_rect(fill="gray85"))+
  theme(axis.text.y = element_text(size="14", color="black"), axis.title.y = element_text(face="bold",size="14", color="black"))+theme(axis.text.x = element_text(size="12", color="black", angle = 45), axis.title.x.bottom = element_text(face="bold",size="14", color="black"))+theme(axis.text.x.top= element_text(face="bold",size="14", color="black"),
        legend.position = "bottom",panel.spacing = unit(0.2, "mm"),               
        panel.border = element_rect(colour="grey"))+
  scale_y_continuous(limits=c(0,1))
Div_SIM
## Statistics
Simp.aov <- aov(Simpson ~ ProjectName:Treatment+ProjectName+Treatment, data = adiv)
summary(Simp.aov)
TukeyHSD(Simp.aov)
tukey_Trial2 <-HSD.test(Simp.aov, "ProjectName", group = TRUE)
tukey_Trial2
kruskal.test(Simpson ~ ProjectName, data=adiv)
```
##### Richness
```
Div_rich <- ggplot(adiv, aes(as.factor(ProjectName), Chao1.Chao1)) +
  facet_grid(.~Hatchery, scales="free",space="free_x",switch="y")+
  geom_boxplot(aes(fill=factor(Treatment)),width=0.5)+
  #geom_signif(comparisons = list(c("A", "B"), c("B","C")),
              #map_signif_level=TRUE)+
  scale_fill_manual(values=c("blue", "orange","blue", "orange")) +
  labs(y="Chao1 Richness", x="Trial", fill="Treatment")+theme_bw()+
  theme(strip.text.x = element_text(size=12,face="bold"))+
  theme(strip.background = element_rect(fill="gray85"))+
  theme(axis.text.y = element_text(size="14", color="black"), axis.title.y = element_text(face="bold",size="14", color="black"))+theme(axis.text.x = element_text(size="14", color="black", angle = 45), axis.title.x.bottom = element_text(face="bold",size="14", color="black"))+theme(axis.text.x.top= element_text(face="bold",size="14", color="black"),
        legend.position = "bottom",panel.spacing = unit(0.2, "mm"),               
        panel.border = element_rect(colour="grey"))+
  scale_y_continuous(limits=c(0,300))
Div_rich
##Statistics
Rich.aov <- aov(Chao1.Chao1 ~ ProjectName:Treatment+ProjectName, data = adiv)
summary(Rich.aov)
tukey_Trial2 <-HSD.test(Rich.aov, "ProjectName", group = TRUE)
TukeyHSD(Rich.aov)
tukey_Trial2
kruskal.test(Chao1.Chao1 ~ Treatment, data=adiv)
kruskal.test(Chao1.Chao1 ~ ProjectName, data=adiv)
```
#### Shared and unique ASVs across samples
```
#By Trial
ps2A <- prune_taxa(taxa_sums(ps2) > 0, ps2)
upsetda <- get_upset(obj=ps2A, factorNames="ProjectName")
U <- upset(upsetda, sets = c("T1", "T2", "T3", "T4", "T5", 
    "T6", "T7", "T8"), sets.bar.color = "#56B4E9",
      order.by = "freq", keep.order = TRUE, queries = list(list(query = intersects, params = list("T1", "T2", "T3", "T4", "T5", 
    "T6", "T7", "T8"), color = "orange", active = T),
      empty.intersections = "on",group.by = "sets", cutoff = 5))
upset(upsetda, sets.bar.color = "darkblue",nintersects = 50, sets = c("T1", "T2", "T3", "T4", "T5","T6", "T7", "T8"), queries = list(list(query = intersects, params = list("T1", "T2", "T3", "T4", "T5", 
    "T6", "T7", "T8"), color = "orange", active = T)),order.by = "freq", empty.intersections = "on",keep.order = TRUE,number.angles = 0, point.size = 2.5,line.size = 1, mb.ratio = c(0.55, 0.45),
    mainbar.y.label = "ASVs Intersections", sets.x.label = "Total number of ASV per trial",text.scale = c(2.2, 2, 1.3, 1.5, 1.7, 1.3))
```    
```
#By Treatment    
upsetda <- get_upset(obj=ps2A, factorNames="Treatment")
U <- upset(upsetda, sets = c("S4", "C"), sets.bar.color = "#56B4E9",
      order.by = "freq", keep.order = TRUE,
      empty.intersections = "on",group.by = "sets", cutoff = 5)
upset(upsetda, sets.bar.color = "darkblue",nintersects = 40, sets = c("S4", "C"),order.by = "freq", empty.intersections = "on",keep.order = TRUE,number.angles = 0, point.size = 3.5,line.size = 1, 
    mainbar.y.label = "ASVs Intersections", sets.x.label = "ASV Per treatment in all trials",text.scale = c(2.2, 2, 1.5, 1.4, 2.7, 2.5))
```

#### NMDS Plot

###### Log transformation moves to a more normal distribution
```
#GPA1 <- ps2 %>% transform_sample_counts(function(x) log(1 + x))%>% psmelt()
all_order <- GPA1 %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order)) 
SD 
```
##### By Hatchery
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
          Treatment=as.factor(SD[,1]),Hatchery=as.factor(SD[,6]))

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Hatchery),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Hatchery)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Hatchery==g,], veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Hatchery=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Hatchery),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Hatchery, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Hatchery))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Hatchery), size=1)+
  scale_linetype_manual(values=c("VIMS"="solid", "MOOK"="longdash","Mat"="dotdash","RWU"="dashed"))+
  scale_colour_manual(values=c("VIMS"="gold2", "MOOK"="hotpink4","Mat"="cyan","RWU"="darksalmon"))+
  ggtitle("By Hatchery")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=Treatment), size=3) + scale_shape_manual(values=c(15,16)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Hatchery <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Hatchery
```
##### By Location
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
          Treatment=as.factor(SD[,1]),Location=as.factor(SD[,7]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Location),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Location)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Location==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Location=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Location),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Location, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Location))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Location), size=1)+
  scale_linetype_manual(values=c("Virginia"="solid","Maine"="dotted", "RhodeIsland"="longdash"))+
  scale_colour_manual(values=c("Virginia"="darkgreen","Maine"="purple", "RhodeIsland"="red"))+
  ggtitle("By Location")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=Treatment), size=3) + scale_shape_manual(values=c(15,16)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Location <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Location
```
##### By Year
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
         ProjectName=as.factor(SD[,2]),Year=as.factor(SD[,9]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Year),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Year)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Year==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Year=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Year),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Year, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Year))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Year), size=1)+
  scale_linetype_manual(values=c("2019"="solid","2020"="dotted", "2021"="longdash"))+
  scale_colour_manual(values=c("2019"="gray67","2020"="skyblue", "2021"="royalblue4"))+
  ggtitle("By Year")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=ProjectName), size=3) + scale_shape_manual(values=c(15,16,17,18,19,20,7,8,6)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Year <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Year
```
##### By Month
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
         Name=as.factor(SD[,5]),Month=as.factor(SD[,10]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Month),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Month)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Month==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Month=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Month),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Month, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Month))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Month), size=1)+
  scale_linetype_manual(values=c("January"="solid","May"="dotted", "June"="longdash","July"="dashed"))+
  scale_colour_manual(values=c("January"="coral","May"="gold", "June"="deeppink","July"="orange"))+
  ggtitle("By Month")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=Name), size=3) + scale_shape_manual(values=c(15,16,17,18,19,20,7,8,6)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Month <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Month
```
##### By Season
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
    ProjectName=as.factor(SD[,2]),Season=as.factor(SD[,11]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Season),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Season)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Season==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Season=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Season),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Season, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Season))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Season), size=1)+
  scale_linetype_manual(values=c("Spring"="solid","Summer"="dotted"))+
  scale_colour_manual(values=c("Spring"="darkgreen","Summer"="purple"))+
  ggtitle("By Season")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=ProjectName), size=3) + scale_shape_manual(values=c(15,16,17,18,19,20,7,8,6)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Season <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Season
```
##### By trial
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
         ProjectName=as.factor(SD[,2]), Treatment=as.factor(SD[,1]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$ProjectName),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$ProjectName)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$ProjectName==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),ProjectName=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$ProjectName),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~ProjectName, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDS_t1 <- cbind(SampleID = rownames(NMDS_t1 ), NMDS_t1)
rownames(SD1 ) <- NULL
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=ProjectName))+
  geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=3,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 #annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
           label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=ProjectName), size=1)+
  scale_linetype_manual(values=c("T1"="solid","T2"="dotted", "T3"="longdash","T4"="dashed","T5"="dotdash","T6"="solid","T7"="twodash","T8"="dotted"))+
  scale_colour_manual(values=c("T1"="darkgreen","T2"="purple", "T3"="red","T4"="Orange","T5"="chocolate4","T6"="blue","T7"="cyan","T8"="gray0"))+
  ggtitle("By Trial")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=Treatment), size=3) + scale_shape_manual(values=c(15,16,17,18,19,20,7,8,6)) 
  # reorder the legend
  guides(color = guide_legend(order=1), shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Trial <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Trial
```
##### By treatment
```
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(SD[,1]),ProjectName=as.factor(SD[,2]))
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                    veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1-1),y=min(NMDS_t1$NMDS2-1),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="#009E73","S4"="#652926"))+
  ggtitle("By Treatment")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=Treatment), size=3) + scale_shape_manual(values=c(15,16)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_Treat <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 16, face = "bold", hjust = 0.5))
NMDSplot_Treatment
```
#### By each trial

##### VIMS_Trial1
```
T1 <- subset_samples(ps2, ProjectName=="T1")
Table1 <- otu_table(T1)
T1_SD <- sample_data(T1) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T1_SD1 <- cbind(SampleID = rownames(T1_SD ), T1_SD)
GPA1 <- T1 %>% transform_sample_counts(function(x) log(1 + x))%>% psmelt()
all_order <- GPA1 %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T1_SD[,1]),TankLocation=as.factor(T1_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T1_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T1_SD, by=NULL,method="bray", k=2)
#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
#geom_segment(data=df_biofit, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),arrow = arrow(length = unit(0.4, "cm")),color="black",alpha=0.05)+
 # geom_text(data=as.data.frame(df_biofit*1.1),aes(NMDS1, NMDS2, label = rownames(df_biofit)),color="black",alpha=0.5)+
  #annotate("text",x=min(NMDS_t1$NMDS1-1),y=min(NMDS_t1$NMDS2-1), 
          # label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial1,VIMS June 2019")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,18,7,11,13,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T1 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T1 
```
##### VIMS_Trial2
```
T2 <- subset_samples(ps2, ProjectName=="T2")
Table1 <- otu_table(T2)
T2_SD <- sample_data(T2) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T2_SD1 <- cbind(SampleID = rownames(T2_SD ), T2_SD)
GPA <- T2 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)

NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T2_SD[,1]),TankLocation=as.factor(T2_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T2_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T2_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  #geom_segment(data=sig.spp.scrs, aes(xend = NMDS1, yend = NMDS2),arrow = arrow(length = unit(0.4, "cm")),color="black",alpha=0.05)+
  #geom_text(data=as.data.frame(sig.spp.scrs*1.1),aes(NMDS1, NMDS2, label = rownames(sig.spp.scrs)),color="black",alpha=0.5)+
  #annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           #label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial2,VIMS July 2019")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,15,16,16,16,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T2 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T2
```
##### VIMS_Trial3
```
T3 <- subset_samples(ps2, ProjectName=="T3")
Table1 <- otu_table(T3)
T3_SD <- sample_data(T3) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T3_SD1 <- cbind(SampleID = rownames(T3_SD ), T3_SD)
GPA <- T3 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)

NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T3_SD[,1]),TankLocation=as.factor(T3_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())

plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T3_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T3_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  #annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5), 
         #  label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial3,VIMS May 2020")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,18,16,16,16,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T3 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T3
```
##### VIMS_Trial4
```
T4 <- subset_samples(ps2, ProjectName=="T4")
Table1 <- otu_table(T4)
T4_SD <- sample_data(T4) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T4_SD1 <- cbind(SampleID = rownames(T4_SD ), T4_SD)
GPA <- T4 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T4_SD[,1]),TankLocation=as.factor(T4_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T4_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T4_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
           label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1-1),y=min(NMDS_t1$NMDS2-1), 
          # label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial4,VIMS June 2020")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,18,16,16,16,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T4 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T4
```

##### MOOK_Trial7
```
T7 <- subset_samples(ps2, Name=="T7")
Table1 <- otu_table(T7)
T7_SD <- sample_data(T7) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T7_SD1 <- cbind(SampleID = rownames(T7_SD ), T7_SD)
GPA <- T7 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T7_SD[,1]),TankLocation=as.factor(T7_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())

plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T7_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T7_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
  annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
           label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  #annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5), 
          # label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial7,MOOK June 2021")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,18,7,11,13,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T7 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T7
```
##### RWU_Trial8
```
T8 <- subset_samples(ps2, ProjectName=="T8")
Table1 <- otu_table(T8)
T8_SD <- sample_data(T8) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T8_SD1 <- cbind(SampleID = rownames(T8_SD ), T8_SD)
GPA <- T8 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
            Name=as.factor(T8_SD[,5]),TankLocation=as.factor(T8_SD[,6]),Treatment=as.factor(T8_SD[,1]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())

plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T8_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T8_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
  #annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
         #  label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  #annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5), 
           #label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="blue","S4"="orange"))+
  ggtitle("Trial8,RWU July 2021")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=8) + scale_shape_manual(values=c(15,16,17,18,16,16,16,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T8 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="18", color="black"), axis.title.y = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x = element_text(size="18", color="black"), axis.title.x.bottom = element_text(face="bold",size="18", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="20", color="black"))+theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
NMDSplot_T8
```
#### Abundance plot
##### phylum abundance
```
SD <- sample_data(ps2) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation","sum","Hatchery") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
Orderabundance <- ps2 %>%
tax_glom(taxrank = "Phylum") %>% # agglomerate at phylum level
transform_sample_counts(function(x) {x/sum(x)*100}) %>% # Transform to rel. abundance
  psmelt() %>% 
  #filter(Abundance > 0.01)%>% # Melt to long format
  arrange(Phylum)
head(Orderabundance)
all_order <- Orderabundance %>%group_by(Sample,Phylum)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="Phylum")
Order1c<- as.data.frame(t(Order))
Order1b <-Order1c[,order(colSums(Order1c),decreasing=TRUE)]
#Extract list of top 12 Taxa and assign the remaining as "Others"
N<-4 # can change number of taxa to plot here
taxa_list<-colnames(Order1b)[1:N]
#Generate a new table with everything added to Others
Order1d <-data.frame(Order1b[,colnames(Order1b) %in% taxa_list], Others=rowSums(Order1b[,!colnames(Order1b) %in% taxa_list]))
Order1d_long <-gather(Order1d, SampleID, Value)
Order1d_long$Taxa <-factor(Order1d_long$SampleID, levels = unique(Order1d_long$SampleID))
#Make a new data fram called meta_taxa and merge the metadata and long taxa data
Order1d_count<-as.data.frame(c(SD1, Order1d_long))
###Get the colours to use for plot
my_colors <- c("#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A","#F6AE2D","#86BBD8")
Order1d_count$Name <- as.character(Order1d_count$Name)
Order1d_count$Name <- factor(Order1d_count$Name, levels=c("T1","T2","T3","T4","T5","T6","T7","T8NUV","T8UV"))
Order1d_count$Treatment <- as.character(Order1d_count$Treatment)
Order1d_count$Treatment <- factor(Order1d_count$Treatment, levels=unique(Order1d_count$Treatment))
Order1d_count$Type <- as.character(Order1d_count$Type)
# plot the percent ASV abundances
Order_plot <- ggplot(Order1d_count)+
  geom_col(mapping = aes(x = SampleID, y = Value, fill = Taxa), position = "fill", show.legend = TRUE)+
  ylab("Percent Phylum Abundance") +
  # facet_wrap(. ~ProjectName+Type, scales = "free_x") +
  theme(strip.background = element_rect(colour="black", size=1))+
  theme(strip.text.x = element_text(size=30,face="bold"))+
  facet_grid(.~Name+Treatment+Hatchery, drop=TRUE,scale="free", space="free",switch="x")+
  scale_fill_manual(values = my_colors) +
  xlab("ProjectName")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme_minimal()+
   theme(axis.text.x = element_blank()) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.2,"lines"),axis.text.y.left = element_text(colour="black", size=14, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        title = element_text(size = 16, face="bold"))+ theme(legend.position="none")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
  labs(x="ProjectName", y="Percentage Phylum Abundance")

# plot the number of reads per sample
Read_plot <-ggplot(SD1,aes(SampleID,sum, fill=Treatment))+
  # plot using columns
  geom_col()+theme_grey()+
  # facet by type, then trial, then day, then treatment group.
  facet_grid(.~Name+Treatment+Hatchery, drop=TRUE, scale="free_x")+
  # label the axes
  labs(y="Read Abundance", x=NULL)+
  # change the fill color to Control= light blue, Treatment=dark red.
  scale_fill_manual(values=c("gray67","gray5","gray67","gray5","gray67","gray5"))+
  scale_y_continuous(expand = c(0,0),labels=scales::scientific_format(digits=1))+
  # edit the facets a little
  theme(strip.text = element_blank(),legend.position="none",
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.4,"lines"),axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size="14", color="black",face = "bold"), axis.title.y = element_text(face="bold",size="16", color="black"))

# get the legend for each plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
ASV_legend<-g_legend(ggplot(Order1d_count,aes(SampleID,Value,fill=Taxa))+                    geom_col()+scale_fill_manual(values= my_colors)+labs(fill="Taxa"))
read_legend<-g_legend(ggplot(SD1,aes(SampleID,sum, fill=Treatment))+
                     geom_col()+scale_fill_manual(values=c("#00545b","#ff856d","#640025","#3ddda5","#cdffaa","#150e00","#bae278"))+labs(fill="Treatment"))

# plot the figures and legends together

Order_plot1 <-grid.arrange(Read_plot, read_legend, Order_plot, ASV_legend, nrow=2, ncol=2, heights=c(1,2.5), widths=c(12.5,12.5))
```

##### Order abundance
```
SD <- sample_data(ps2) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation", "Hatchery","sum") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
Orderabundance <- ps2 %>%
tax_glom(taxrank = "Order") %>% # agglomerate at phylum level
transform_sample_counts(function(x) {x/sum(x)*100}) %>% # Transform to rel. abundance
  psmelt() %>% 
  #filter(Abundance > 0.01)%>% # Melt to long format
  arrange(Order)
head(Orderabundance)
all_order <- Orderabundance %>%group_by(Sample,Order)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="Order")
Order1c<- as.data.frame(t(Order))
Order1b <-Order1c[,order(colSums(Order1c),decreasing=TRUE)]
#Extract list of top 12 Taxa and assign the remaining as "Others"
N<-12 # can change number of taxa to plot here
taxa_list<-colnames(Order1b)[1:N]
#Generate a new table with everything added to Others
Order1d <-data.frame(Order1b[,colnames(Order1b) %in% taxa_list], Others=rowSums(Order1b[,!colnames(Order1b) %in% taxa_list]))
Order1d_long <-gather(Order1d, SampleID, Value)
Order1d_long$Taxa <-factor(Order1d_long$SampleID, levels = unique(Order1d_long$SampleID))
#Make a new data fram called meta_taxa and merge the metadata and long taxa data
Order1d_count<-as.data.frame(c(SD1, Order1d_long))
###Get the colours to use for plot
my_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#808080","#94FFB5","#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
"#F6AE2D","#86BBD8")
Order1d_count$Name <- as.character(Order1d_count$Name)
Order1d_count$Name <- factor(Order1d_count$Name, levels=c("T1","T2","T3","T4","T5","T6","T7","T8NUV","T8UV"))
Order1d_count$Treatment <- as.character(Order1d_count$Treatment)
Order1d_count$Treatment <- factor(Order1d_count$Treatment, levels=unique(Order1d_count$Treatment))
Order1d_count$Type <- as.character(Order1d_count$Type)
# plot the percent ASV abundances
Order_plot <- ggplot(Order1d_count)+
  geom_col(mapping = aes(x = SampleID, y = Value, fill = Taxa), position = "fill", show.legend = TRUE)+
  ylab("Percent Order Abundance") +
  # facet_wrap(. ~ProjectName+Type, scales = "free_x") +
  theme(strip.background = element_rect(colour="black", size=1))+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  facet_grid(.~Name+Treatment+Hatchery, drop=TRUE,scale="free", space="free",switch="x")+
  scale_fill_manual(values = my_colors) +
  xlab("ProjectName")+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme_minimal()+
   theme(axis.text.x = element_blank()) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.2,"lines"),axis.text.y.left = element_text(colour="black", size=14, face="bold"),
        axis.title.y = element_text(size = 14, face="bold"),
        title = element_text(size = 14, face="bold"))+ theme(legend.position="bottom")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
  labs(x="ProjectName", y="Percentage Phylum Abundance")

# plot the number of reads per sample
Read_plot <-ggplot(SD1,aes(SampleID,sum, fill=Treatment))+
  # plot using columns
  geom_col()+theme_grey()+
  # facet by type, then trial, then day, then treatment group.
  facet_grid(.~Name+Treatment+Hatchery, drop=TRUE, scale="free_x")+
  # label the axes
  labs(y="Read Abundance", x=NULL)+
  # change the fill color to Control= light blue, Treatment=dark red.
  scale_fill_manual(values=c("gray67","gray5","#640025","#3ddda5","#cdffaa","#150e00","#bae278"))+
  scale_y_continuous(expand = c(0,0),labels=scales::scientific_format(digits=1))+
  # edit the facets a little
  theme(strip.text = element_blank(),legend.position="none",
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.4,"lines"),axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size="14", color="black",face = "bold"), axis.title.y = element_text(face="bold",size="14", color="black"))

# get the legend for each plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
ASV_legend<-g_legend(ggplot(Order1d_count,aes(SampleID,Value,fill=Taxa))+                    geom_col()+scale_fill_manual(values= my_colors)+labs(fill="Taxa"))
read_legend<-g_legend(ggplot(SD1,aes(SampleID,sum, fill=Treatment))+
                     geom_col()+scale_fill_manual(values=c("gray67","gray5","#640025","#3ddda5","#cdffaa","#150e00","#bae278"))+labs(fill="Treatment"))

# plot the figures and legends together
Order_plot1 <-grid.arrange(Read_plot, read_legend, Order_plot, ASV_legend, nrow=2, ncol=2, heights=c(1,2.5), widths=c(9,3.5))

Order_plot1 <-grid.arrange(Read_plot, Order_plot, nrow=2, heights=c(1,2.5), widths=c(9,3.5))
```

#### Differential Abundance
##### By Trial
```
xt <- transform(ps2, 'hellinger')
Orderabundance <- xt %>% tax_glom(taxrank = "Order")
deres <- diff_analysis(obj =Orderabundance, classgroup = "Name",
                       mlfun = "lda",alltax = TRUE,
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 4,
                       subclwilc = TRUE,
                       secondalpha = 0.05,
                       lda=2)
es_pb <- ggeffectsize(obj=deres, 
                     lineheight=0.3,
                     linewidth=0.2,   pointsize = 4.5,) + 
        scale_color_manual(values=c("grey","blue","orange","red", "chocolate4","gray0", "cyan", "green", "pink1","#599861"))+theme(axis.text.x = element_text(colour="black", size=20, face="bold")) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.2,"lines"),axis.text.y.left = element_text(colour="black", size=20, face="bold"),axis.text.y.right = element_text(colour="black", size=20, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),axis.title.x = element_text(size = 20, face="bold", colour = "black"),
        title = element_text(size = 20, face="bold"))+labs(title = "Trial")+ylab("Taxa")+theme(legend.position="right")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
es_pb
```
##### By Treatment
```
xt <- transform(ps2, 'hellinger')
Orderabundance <- xt %>% tax_glom(taxrank = "Genus")
deres <- diff_analysis(obj =xt, classgroup = "Treatment",
                       mlfun = "lda",alltax = TRUE,
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.01,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 4,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=2)
es_pb <- ggeffectsize(obj=deres, 
                     lineheight=0.3,
                     linewidth=0.2,   pointsize = 4.5,) + 
        scale_color_manual(values=c("#009E73","#652926",
                                    "orange","red", "chocolate4","gray0", "cyan", "green", "pink1","#599861"))+theme(axis.text.x = element_text(colour="black", size=20, face="bold")) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.2,"lines"),axis.text.y.left = element_text(colour="black", size=20, face="bold"),axis.text.y.right = element_text(colour="black", size=20, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),axis.title.x = element_text(size = 20, face="bold", colour = "black"),
        title = element_text(size = 20, face="bold"))+labs(title = "Trial")+ylab("Taxa")+theme(legend.position="bottom")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
es_pb

diffboxT1 <- ggdiffbox(obj=deres, box_notch=FALSE, lineheight=0.8,
                     linewidth=0.6,   pointsize = 3.5,
             colorlist=c("#009E73","#652926"), l_xlabtext="relative abundance",l_ylabtext="Taxa")+ theme(text=element_text())
diffboxT1

+ theme_set(theme_gray())+theme(axis.text.x = element_text(colour="black", size=20, face="bold")) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.2,"lines"),axis.text.y.left = element_text(colour="black", size=20, face="bold"),axis.text.y.right = element_text(colour="black", size=20, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),axis.title.x = element_text(size = 20, face="bold", colour = "black"),
        title = element_text(size = 20, face="bold"))+labs(title = "Trial5")+ylab("Taxa")+theme(legend.position="bottom")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
```
##### Wilconxin rank sum tes
```
ps2.log <- transform_sample_counts(ps2, function(x) log(1 + x))
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps2.log))))
ps_wilcox$Status <- phyloseq::sample_data(ps2.log)$Treatment
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Status, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Status, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = ASV, value = abund, -Status) %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))    
#Show results
head(wilcox_results)
head(wilcox_results$data[[1]])
wilcox_results$wilcox_test[[1]]
#Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(ASV, p_value) %>%
  unnest()
alpha <- 0.01
wilcox_results1 <- wilcox_results[which(wilcox_results$p_value < alpha), ]
#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps2.log))
taxa_info <- taxa_info %>% rownames_to_column(var = "ASV")
#Computing FDR corrected p-values
wilcox.T1 <- wilcox_results1 %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.01) %>%
  dplyr::select(ASV, p_value, BH_FDR, everything())
write.csv(wilcox.T1, "Treatment_wilcox.csv", sep=",", row.names=FALSE, col.names=FALSE)
set.seed(1024)
```

#### Effect of Probiotics on Vibrio Community
##### Vibrio Diversity
```
(vibrio1 <- subset_taxa(ps2, Order== "Vibrionales"))
vibrio1_rarefy <- rarefy_even_depth(vibrio1, rngseed= 81, sample.size = min(sample_sums(vibrio1)))
sample_data(vibrio1)$ProjectName <- factor(sample_data(vibrio1)$ProjectName, levels = c("T1", "T2", "T3", "T4","T5","T6","T7","T8"))
levels(sample_data(vibrio1)$ProjectName)
tax_table(vibrio1)
sample_data(vibrio1)
summary(sample_sums(vibrio1))
reads_sample <- readcount(vibrio1)
vibrio2 <-psmelt(vibrio1) %>%as_tibble
vibrio3 <- vibrio2 %>% group_by(Sample,OTU) %>% summarise(Value=mean(Abundance))
Md <-spread(vibrio3, Sample, Value) ##spread the data into a wide format
#Me <- column_to_rownames(Md)
Me <-Md %>% remove_rownames %>% column_to_rownames(var="OTU")
Me <-as.data.frame(t(Me))
unique_taxa(vibrio1, treatment = 'Treatment')
summary(sample_sums(vibrio1))
plot_taxa_prevalence(vibrio1, "OTU")
#SD$SimpsonsVib <-diversity(Me, index="simpson")
Vibrio_rarefy <- rarefy_even_depth(vibrio1, rngseed= 81, sample.size = min(sample_sums(vibrio1)))
##Order the x-axis
newSTorder = c("ConS", "S4S","ConF","S4F")
sample_types <- c("ConS", "S4S","ConF","S4F")
sample_labels <- c("Control", "S4","Control", "S4")
adiv <- data.frame(
  "Chao1" = phyloseq::estimate_richness(Vibrio_rarefy, measures = "Chao1"),
  "Shannon" = phyloseq::estimate_richness(Vibrio_rarefy, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(Vibrio_rarefy, measures = "Simpson"),
  "Status" = phyloseq::sample_data(Vibrio_rarefy)$Treatment)
adiv$ProjectName <-as.character(SD$ProjectName)
adiv$Treatment <-as.character(SD$Treatment)
adiv$Num <-as.character(SD$Num)
#adiv$Type <- factor(adiv$Status, levels = newSTorder)

DIV <- adiv %>% 
  ggplot(aes(x=ProjectName,y=Simpson, label = Num, color=Status))+
  geom_boxplot(width=.5)+
  #facet_wrap(~ System, scales = "free")+
  #facet_grid(~Status, scales = "free") +
  # jittered text with geom_text
  #geom_text(check_overlap = FALSE,
           # position=position_jitter(width=0.05))+
  labs(y="Vibrio Simpson's Diversity Index ", x="Trial")+
scale_colour_manual(values=c("red", "purple", "green","purple"))+   #set fill colors
  #scale_x_discrete(                                                  #change x-axis labels
   # breaks = sample_types, 
   # labels = sample_labels)+ 
  scale_y_continuous(limits=c(0,1))
DIV_S <- DIV +theme(axis.text.y = element_text(size="14", color="black",face = "bold"), axis.title.y = element_text(face="bold",size="14", color="black"), axis.text.x = element_text(angle = 45,hjust=1, size="14", color="black",face = "bold"), axis.title.x = element_text(face="bold",size="14", color="black"),legend.position = "bottom")

aov.simpsons.vibrio = aov(adiv$Simpson ~Treatment+ProjectName+Treatment:ProjectName, data=adiv)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.simpsons.vibrio)
#using tukeyHDS test to between groups and correct for multiple  comparison
tukey_vib <-HSD.test(aov.simpsons.vibrio, "Treatment:ProjectName", group = TRUE)
tukey_vib
kruskal.test(Simpson ~ Treatment, data=adiv)
kruskal.test(Simpson ~ ProjectName, data=adiv)
```
##### Vibrio Abundance
```
SD <- sample_data(ps2.T1) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation","sum") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
Orderabundance <- ps2.T1 %>%
tax_glom(taxrank = "Order") %>% # agglomerate at phylum level
transform_sample_counts(function(x) {x/sum(x)*100}) %>% # Transform to rel. abundance
  psmelt() %>%     # Melt to long format
  arrange(Order)
head(Orderabundance)
all_order <- Orderabundance %>%group_by(Sample,Order,)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="Order")
#Transpose the data to have sampleID as rows
Order1c<- as.data.frame(t(Order))
Order1c<- cbind(SampleID = rownames(Order1c ), Order1c)
rownames(Order1c) <- NULL
trialorder_A <-as.data.frame(c(SD1, Order1c))
trialorder1 <- group_by(trialorder_A,Treatment)
trial_Vibrionales_means <- dplyr::summarise(trialorder1,
                         meanvib = mean(Vibrionales, na.rm=TRUE),
                         stdevvib = sd(Vibrionales, na.rm=TRUE))

Vibrionales.aov <- aov(trialorder1$Vibrionales ~ Treatment,  data = trialorder1 )
summary(Vibrionales.aov)
tukey_result2 <- HSD.test(Vibrionales.aov, "Treatment", group = TRUE)
print(tukey_result2)
group_data2 <- tukey_result2$groups[order(rownames(tukey_result2$groups)),]

plotVibrionales <-ggplot(trial_Vibrionales_means, aes(x = Treatment, y = meanvib, by=Treatment))+
  geom_bar(stat="identity", aes(fill=Treatment), position="dodge")+
  geom_errorbar(aes(ymin=pmax(meanvib-stdevvib,0), ymax=meanvib+stdevvib),
                width=0.1, position=position_dodge(.5))+
  geom_text(data = data.frame(),
            aes(x = rownames(group_data2), y = max(trialorder1$Vibrionales), label = group_data2$groups),
            col = 'black',
            size = 6) +
  labs(y="Relative Percent Vibrionales  Abundance", x="Treatment")+
  scale_fill_manual(values=c("darkgreen", "red","darkgreen", "red","purple"))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),strip.text = element_text(size = 12),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 12, face = "bold"),legend.position = "none")+
  scale_y_continuous(limits=c(0,100))
plotVibrionales
```
#### NMDS Plot
###### Trial1
```
Vib1 <- subset_samples(vibrio1, Name=="T1")
VibT1 <- otu_table(Vib1)
T1_SD <- sample_data(Vib1) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T1_SD1 <- cbind(SampleID = rownames(T1_SD ), T1_SD)

GPA <- VibT1 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T1_SD[,1]),TankLocation=as.factor(T1_SD[,6]))

library(vegan)
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())

plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T1_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T1_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5), 
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="darkblue","S4"="darkred"))+
  ggtitle("Trial1,VIMS June 2019")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=5) + scale_shape_manual(values=c(15,16,17,18,7,11,13,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T1 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 12, face = "bold", hjust = 0.5))
NMDSplot_T1 
```
###### Trial2
```
Vib2 <- subset_samples(vibrio1, Name=="T2")
VibT2 <- otu_table(Vib2)
T2_SD <- sample_data(Vib2) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name","TankLocation") %>%
    mutate_if(is.factor,as.character)
T2_SD1 <- cbind(SampleID = rownames(T2_SD ), T2_SD)
GPA <- VibT2 %>% transform_sample_counts(function(x) {x/sum(x)})%>% psmelt()
all_order <- GPA %>%group_by(Sample,OTU)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="OTU")
Abun <- as.data.frame(t(Order))
sol_t1 <-metaMDS(Abun ,distance = "bray", k = 2, trymax = 50)
NMDS_t1=data.frame(NMDS1=sol_t1$point[,1],NMDS2=sol_t1$point[,2],
                Treatment=as.factor(T2_SD[,1]),TankLocation=as.factor(T2_SD[,6]))
# define shapes for plots
shape_values<-seq(1,11)
# set theme for following plots
theme_set(theme_gray())
plot.new()
ord_t1<-ordiellipse(sol_t1, as.factor(T2_SD$Treatment),
                    display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Generate ellipse points based on 95% confidence (SD) intervals
df_ell_t1 <- data.frame()
for(g in levels(NMDS_t1$Treatment)){
  if(g!="" && (g %in% names(ord_t1))){
    df_ell_t1 <- rbind(df_ell_t1,
    cbind(as.data.frame(with(NMDS_t1[NMDS_t1$Treatment==g,],                            veganCovEllipse(ord_t1[[g]]$cov,ord_t1[[g]]$center,ord_t1[[g]]$scale))),Treatment=g))}}
#Store the center of the ellipses (mean NMDS)
NMDS.mean_t1=aggregate(NMDS_t1[,1:2],list(group=NMDS_t1$Treatment),mean)
#Calculate p-value:
adon_t1<-adonis2(Abun ~Treatment, data=T2_SD, by=NULL,method="bray", k=2)

#Plot NMDS
NMDSplot_t3<-ggplot(data=NMDS_t1,aes(NMDS1,NMDS2,col=Treatment))+
  #geom_text(data=NMDS_t1,aes(NMDS1,NMDS2,label=SampleID),size=5,vjust=0)+
  # label the middle of the ellipses with the name of the grouping factor
 # annotate("text",x=NMDS.mean_t1$x,y=NMDS.mean_t1$y,
          # label=NMDS.mean_t1$group,size=5)+
  # add the p-value in the bottom right corner
  annotate("text",x=min(NMDS_t1$NMDS1),y=min(NMDS_t1$NMDS2-0.5),
           label=paste("p= ", adon_t1$`Pr(>F)`[1]),size=3)+
  # draw the ellipses. define color based on the grouping factor
  geom_path(data=df_ell_t1, aes(x=NMDS1, y=NMDS2, linetype=Treatment), size=1)+
  scale_linetype_manual(values=c("C"="solid","S4"="dotted"))+
  scale_colour_manual(values=c("C"="darkblue","S4"="darkred"))+
  ggtitle("Other trials")+
  # add the points per sample. define shape based on Day
  geom_point(aes(shape=TankLocation), size=5) + scale_shape_manual(values=c(15,16,17,15,16,16,16,25,22,1,5,7,9,6,14)) +
  # reorder the legend
  guides(color = guide_legend(order=1),lty= guide_legend(order=1),
         shape = guide_legend(order=2), legend.position = "NONE")
NMDSplot_T2 <- NMDSplot_t3 + theme(axis.text.y = element_text(size="12", color="black"), axis.title.y = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x = element_text(size="12", color="black"), axis.title.x.bottom = element_text(face="bold",size="12", color="black"))+
theme(axis.text.x.top= element_text(face="bold",size="12", color="black"))+theme(plot.title=element_text(size = 12, face = "bold", hjust = 0.5))
NMDSplot_T2
```

##### Relative abundance of Vibrio
```
SD <- sample_data(vibrio1) %>%
    data.frame() %>%
    select("Treatment","ProjectName", "Num","Type", "Name") %>%
    mutate_if(is.factor,as.character)
SD1 <- cbind(SampleID = rownames(SD ), SD)
rownames(SD1 ) <- NULL
Orderabundance <- vibrio1 %>%
tax_glom(taxrank = "OTU") %>% # agglomerate at phylum level
transform_sample_counts(function(x) {x/sum(x)*100}) %>% # Transform to rel. abundance
  psmelt() %>% 
  #filter(Abundance > 0.01)%>% # Melt to long format
  arrange(Order)
head(Orderabundance)
all_order <- Orderabundance %>%group_by(Sample,Order,)%>%summarise(avg_abundance = mean(Abundance))
Oder <-spread(all_order, Sample,avg_abundance) ##spread the data into a wide format
Order <-Oder %>% remove_rownames %>% column_to_rownames(var="Order")
#Order <-column_to_rownames(Oder)
#rownames(Oder) <- NULL
#Transpose the data to have sampleID as rows
Order1c<- as.data.frame(t(Order))
Order1b <-Order1c[,order(colSums(Order1c),decreasing=TRUE)]
#Extract list of top 12 Taxa and assign the remaining as "Others"
N<-12 # can change number of taxa to plot here
taxa_list<-colnames(Order1b)[1:N]
#Generate a new table with everything added to Others
Order1d <-data.frame(Order1b[,colnames(Order1b) %in% taxa_list], Others=rowSums(Order1b[,!colnames(Order1b) %in% taxa_list]))
Order1d_long <-gather(Order1d, SampleID, Value)
Order1d_long$Taxa <-factor(Order1d_long$SampleID, levels = unique(Order1d_long$SampleID))
#Make a new data fram called meta_taxa and merge the metadata and long taxa data
Order1d_count<-as.data.frame(c(SD1, Order1d_long))
###Get the colours to use for plot
my_colors <- c("brown4","brown2","blue4","blue","chocolate2","aquamarine4","burlywood1","darkgoldenrod1","gold1","hotpink1","#808080","#94FFB5","#C84248","#8A7C64", "#599861","gray")
Order_plot <- ggplot(Order1d_count)+
  geom_col(mapping = aes(x = SampleID, y = Value, fill = Taxa), position = "fill", show.legend = TRUE)+
  ylab("Percent Order Abundance") +
   facet_wrap(.~ ProjectName+Type, scales = "free_x") +
  theme(strip.text.x = element_text(size=12,face="bold"))+
  labs(title = "Abundant Order (> 1%)")+
  #facet_grid(.~ProjectName+Treatment+Type, drop=TRUE,scale="free", space="free",switch="x")+
  scale_fill_manual(values = my_colors) +
  xlab(NULL)+
  scale_y_continuous(expand = c(0,0),labels=scales::percent)+
  theme_minimal()+
   theme(axis.text.x = element_blank()) +
  theme(strip.background = element_rect(fill="gray85"),panel.spacing = unit(0.1,"lines"),axis.text.y.left = element_text(colour="black", size=12),
        axis.title.y = element_text(size = 14, face="bold"),
        title = element_text(size = 14, face="bold"))
Order_plot+ theme(legend.position="bottom")+ theme(legend.title = element_text(colour="black", size=10))+ theme(legend.text = element_text(colour="black", size=14, face="bold"))
```
