# Set up environment and get dataset into R
library(tidyverse);library(vegan);library(patchwork)
setwd("~/Dropbox/Documents/Postdoc/MethodsOptimization/Data/")
options(scipen = 999)

#Need a Key for samples
SampleKey = read.csv("MetaDataKey.csv") %>% #Because this is all of the samples we have
  mutate(Sample_ID = case_when(str_detect(Sample_ID, "SISPApolyTC\\-") ~"SISPApolyTCneg",
                               str_detect(Sample_ID, "SISPApolyTC\\+") ~"SISPApolyTCpos",
                               str_detect(Sample_ID, "SISPARPC\\-") ~"SISPARPCneg",
                               str_detect(Sample_ID, "SISPARPC\\+") ~"SISPARPCpos",
                               TRUE ~Sample_ID),
         SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+"))) %>% 
  mutate(SampDigit = factor(as.character(str_match_all(Sample_ID, "\\d+"))),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"), #Spiked samples are in the series 10, 11, 12
         ExpGroup = case_when(SampDigit %in% c("7", "10") ~"Amicon",
                              SampDigit %in% c("8", "11") ~"Amicon+CsCl",
                              SampDigit %in% c("9", "12") ~"CsCl+Pellet",
                              str_detect(Sample_ID, "(pos)") ~"Controlpos",
                              str_detect(Sample_ID, "(neg)") ~"Controlneg",
                              str_detect(Sample_ID, "TARA") ~"TARA",
                              str_detect(Sample_ID, "99C") ~"Denature")) %>% 
  mutate(Experiment = ifelse(SampDigit %in% c("7", "8", "9", "10", "11", "12"), 2, 1),
         Prep = case_when(str_detect(Prep, "SISPA") & Spikein == "UnSpiked" ~"SISPA-RP",
                          str_detect(Prep, "SISPA") & Spikein == "Spiked" ~"SISPA-polyT",
                          TRUE ~Prep))


read.csv("Norm_by_mapped_renameheaders.csv") %>% head()

# DNA sequencing stats ----------------------------------------------------
DNAseq_df = read.csv("DNA_Sequencing_Stats.csv") %>% 
  mutate(Sample_ID = factor(Sample_ID, levels = read.csv("DNA_Sequencing_Stats.csv") %>% arrange(Total_Reads) %>% pull(Sample_ID))) %>% 
  mutate(ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit))
#Number of reads
DNAseq_df %>% 
  ggplot(., aes(x = Sample_ID, y = Total_Reads, fill = ExpGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "",
                    values = c("grey","steelblue","black")) +
  theme_minimal() 


DNAseq_df %>% 
  ggplot(., aes(x = ExpGroup, y = Total_Reads, fill = ExpGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "",
                    values = c("grey","steelblue","black")) +
  labs(y = "Reads",
       x = "") +
  theme_minimal() +
  facet_wrap(~ExtMethod) 

DNAseq_df %>% 
  ggplot(., aes(x = ExpGroup, y = Totalbp_Q30, fill = ExpGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "",
                    values = c("grey","steelblue","black")) +
  labs(y = "Bases > Q30",
       x = "") +
  theme_minimal() +
  facet_wrap(~ExtMethod)

# Assembly stuff ----------------------------------------------------------
HeadersComplete = read.table("ContigHeaders_All.txt", header = F) %>% 
  separate(V1, sep = "_", into = c("A", "B", "C", "D", "E")) %>% 
  unite(Sample_ID,B,C) %>% 
  unite(SeqID,D,E, remove = T) %>% 
  select(-A) %>% 
  mutate(Len = as.numeric(str_replace_all(V4, "\\w+\\=", "")),
         Multi = as.numeric(str_replace_all(V3, "\\w+\\=","")),
         Flag = str_replace_all(V2, "\\w+\\=",""),
         Sample_ID = case_when(str_detect(Sample_ID, regex("(7A)_S.")) ~"SISPA7A",
                               str_detect(Sample_ID, regex("(7D)_S.")) ~"SISPA7D",
                               str_detect(Sample_ID, regex("(8A)_S.")) ~"SISPA8A",
                               str_detect(Sample_ID, regex("(8D)_S.")) ~"SISPA8D",
                               str_detect(Sample_ID, regex("(9A)_S.")) ~"SISPA9A",
                               str_detect(Sample_ID, regex("(9D)_S.")) ~"SISPA9D",
                               TRUE ~Sample_ID)) %>% 
  select(Sample_ID, SeqID, Flag, Multi, Len) %>% 
  mutate(LibPrep = case_when(str_detect(Sample_ID, regex("SISPA[789][ABDE]")) ~"SISPApolyT",
                             TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"),
         ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               TRUE ~Sample_ID))

#Conbine these replicates with the old samples and write as a new csv file.
TableFromHeadersMay2025 = read.csv("extable.csv") %>% 
  rename(LibPrep = Group, Spikein = SampSpikein) %>% 
  mutate(Spikein = str_replace_all(Spikein, regex("nospike", ignore_case = T), "UnSpiked")) %>% 
  rbind(HeadersComplete)
#write.csv(TableFromHeadersMay2025, file = "TableFromHeadersMay2025.csv")

# ReAnalyze the unspiked SISPA samples with the additional replicates.
#Summary stats for lengths

#First plot alll of the data for each experiment type
TableFromHeadersMay2025 %>% 
  filter(str_detect(Sample_ID, "SISPA[789][ABDE]"),
         Spikein == "UnSpiked")  %>% 
  mutate(LibPrep = "SISPApolyT") %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit),
         SequenceGroup = ifelse(str_detect(Sample_ID, "\\d[AD]"),"new","old")) %>%  
  ggplot(., aes(x=ExpGroup , y = Len)) +
  geom_boxplot(outlier.size = .3, lwd = .3) +
  stat_summary(fun = "mean", shape = 23, size = .1, color = "red") +
  labs(x = "",
       y = "Contig Length") +
  facet_wrap(~ExtMethod) +
  theme_minimal()

#Creating a table with summary statistics for each Sample type.
Len_summary_table = TableFromHeadersMay2025 %>% 
  filter(str_detect(Sample_ID, "SISPA[789][ABDE]"),
         Spikein == "UnSpiked")  %>% 
  mutate(LibPrep = "SISPApolyT") %>% 
  group_by(Sample_ID, LibPrep,SampDigit,Spikein,ExtMethod) %>% 
  summarise(N_contigs = n(),
            MinLen = min(Len),
            MaxLen = max(Len),
            MeanLen = mean(Len),
            MedLen = median(Len),
            q1Len = quantile(Len, probs = 0.25),
            q3Len = quantile(Len, probs = 0.75),
            SdLen = sd(Len)) %>% 
  ungroup() %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit),
         SequenceGroup = ifelse(str_detect(Sample_ID, "\\d[AD]"),"new","old")) %>% 
  pivot_longer(cols = c(N_contigs, contains("Len")), names_to = "variable", values_to = "value")

# Compare The lengths and number of contigs from the different SampDigit groups
interestgroup = c("N_contigs", "MaxLen", "MedLen")
Len_summary_table %>% 
  filter(variable %in% interestgroup) %>% 
  #ggplot(., aes(x = ExpGroup, y = value)) +
  ggplot(., aes(x = ExpGroup, y = value, color = SequenceGroup)) +
  geom_boxplot() + 
  #geom_jitter(width = .3) +
  labs(y = "Value",
       x = "") +
  theme_minimal() +
  facet_wrap(~variable, scale = "free")

#Plot again to examine differences from sequencing runs
Len_summary_table %>% 
  filter(variable %in% interestgroup) %>% 
  #ggplot(., aes(x = ExpGroup, y = value)) +
  ggplot(., aes(x = ExpGroup, y = value, color = SequenceGroup)) +
  geom_boxplot() + 
  #geom_jitter(width = .3) +
  labs(y = "Value",
       x = "") +
  theme_minimal() +
  facet_wrap(~variable, scale = "free")

# Compare extraction methodologies
Len_summary_table %>% 
  filter(variable %in% interestgroup) %>% 
  ggplot(., aes(x=ExtMethod, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scale = "free") +
  labs(x = "") +
  theme_minimal()



#Boxplot lengths
HeadersComplete %>% 
  ggplot(., aes(y = Sample_ID, x = Len)) +
  geom_boxplot(outlier.size = .3, lwd = .3) +
  stat_summary(fun = "mean", shape = 23, size = .1) +
  labs(x = "Length",
       y = "Sample") +
  theme_minimal()

#Performed an anova and tukey but needs to be vetted
HeadersComplete %>% 
  aov(.$Len ~ .$Sample, data = .) %>%
  TukeyHSD(ordered = TRUE)


# RdRP contigs assemblies -------------------------------------------------
HeadersRdRP_table = read.table("curated_RdRP_contig_headersraw.txt", header = F) %>% 
  separate(V1, sep = "_", into = c("A", "B", "C", "D", "E")) %>% 
  unite(Sample_ID,B,C) %>% 
  unite(SeqID,D,E, remove = T) %>% 
  select(-A) %>% 
  mutate(Len = as.numeric(str_replace_all(V4, "\\w+\\=", "")),
         Multi = as.numeric(str_replace_all(V3, "\\w+\\=","")),
         Flag = str_replace_all(V2, "\\w+\\=",""),
         Sample_ID = case_when(str_detect(Sample_ID, regex("(7A)_S.")) ~"SISPA7A",
                               str_detect(Sample_ID, regex("(7D)_S.")) ~"SISPA7D",
                               str_detect(Sample_ID, regex("(8A)_S.")) ~"SISPA8A",
                               str_detect(Sample_ID, regex("(8D)_S.")) ~"SISPA8D",
                               str_detect(Sample_ID, regex("(9A)_S.")) ~"SISPA9A",
                               str_detect(Sample_ID, regex("(9D)_S.")) ~"SISPA9D",
                               TRUE ~Sample_ID)) %>% 
  select(Sample_ID, SeqID, Flag, Multi, Len) %>% 
  mutate(LibPrep = case_when(str_detect(Sample_ID, regex("SISPA[789][ABDE]")) ~"SISPApolyT",
                             TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"),
         ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               TRUE ~Sample_ID))

#Contrast the lenghs produced from each extraction method
HeadersRdRP_table %>% 
  group_by(Sample_ID, LibPrep,SampDigit,Spikein,ExtMethod) %>% 
  summarise(N_contigs = n(),
            MinLen = min(Len),
            MaxLen = max(Len),
            MeanLen = mean(Len),
            MedLen = median(Len),
            q1Len = quantile(Len, probs = 0.25),
            q3Len = quantile(Len, probs = 0.75),
            SdLen = sd(Len)) %>% 
  ungroup() %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit)) %>% 
  pivot_longer(cols = c(N_contigs, contains("Len")), names_to = "variable", values_to = "value") %>% 
  filter(variable == "N_contigs") %>% 
  #ggplot(., aes(x = ExpGroup, y = value)) +
  ggplot(., aes(x = ExpGroup, y = value)) +
  geom_boxplot() + stat_summary(fun = "mean", shape = 23, size = .1, color = "red") +
  geom_point(size = .8, aes(color = ExtMethod)) +
  #geom_jitter(width = .3) +
  labs(y = "N Contigs",
       x = "") +
  theme_minimal() +
  facet_wrap(~ExtMethod, scale = "free")

#Contrast the lenghs produced from each method
HeadersRdRP_table %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit)) %>%  
  ggplot(., aes(x=ExpGroup , y = Len)) +
  geom_boxplot(outlier.size = .3, lwd = .3) + 
  #geom_jitter(alpha = .1) +
  stat_summary(fun = "mean", shape = 23, size = .1, color = "red") +
  labs(x = "",
       y = "Contig Length") +
  facet_wrap(~ExtMethod) +
  theme_minimal()

# Contrast the number of contigs produced from each method
HeadersRdRP_table %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit)) %>% 
  count(ExtMethod,ExpGroup) %>% 
  ggplot(., aes(x = ExpGroup, y = n, fill = ExtMethod)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(name = "",
                    values = c("grey","steelblue")) +
  labs(y = "N contigs",
       x = "") +
  theme_minimal()


###########################################################################
# Diversity_Calculations --------------------------------------------------

RNAv_OTU_DF = read.csv(file = "Norm_by_mapped_renameheaders.csv", header = T) %>% 
  mutate(vOTU.ID = str_c("vOTU",row_number(), sep = "_"),
         Contig = str_c(vOTU.ID,Contig,sep=":")) %>% 
  select(vOTU.ID, everything()) 

Contig_pairs = RNAv_OTU_DF$Contig
RNAv_OTU_DF = RNAv_OTU_DF %>% select(-Contig) %>% column_to_rownames("vOTU.ID")

#Calculate diversity stats.
library(vegan)
div_InvSimp = RNAv_OTU_DF %>% 
  diversity(., index = "invsimpson", MARGIN = 2) %>% 
  data.frame("InvSimp" = .) %>% 
  rownames_to_column("Sample_ID")

#Calculate Sannon
div_Shannon = RNAv_OTU_DF %>% 
  diversity(., index = "shannon", MARGIN = 2) %>% 
  data.frame("Shannon" = .) %>% 
  rownames_to_column("Sample_ID")

#Calculate Species Richness
div_richness = RNAv_OTU_DF %>% rownames_to_column("ID") %>%  
  pivot_longer(cols = -ID, values_to = "values", names_to = "variables") %>% 
  arrange(variables) %>% mutate(bin = ifelse(values > 0,1,0)) %>% 
  group_by(variables) %>% 
  summarise("Richness" = sum(bin)) %>% 
  rename("Sample_ID" = variables)

TotalDiversity = div_richness %>% pull(Richness) %>% sum()

#Calculate Sequence Abundance
NormvOTU_abundances = RNAv_OTU_DF %>% 
  colSums(.) %>% 
  data.frame("Abundance" = .) %>% 
  rownames_to_column("Sample_ID")

## Create a diveristy tables
Diversity_Table = div_richness %>% 
  left_join(., div_InvSimp) %>% 
  left_join(., div_Shannon) %>% 
  left_join(., NormvOTU_abundances) %>% 
  data.frame() %>% left_join(SampleKey)

Treatments = c("Amicon", "Amicon+CsCl", "CsCl+Pellet", "TARA")

#Long version of the table for data exploration and plotting
Diversity_Table.long = Diversity_Table %>% left_join(AlignmentRate_DF %>% select(., Sample_ID, AlignmentRate)) %>% 
  #drop_na() %>% 
  pivot_longer(cols = c(Richness, InvSimp, Shannon, Abundance, AlignmentRate), names_to = "variable", values_to = "value") %>% 
  mutate(variable = factor(variable, levels = c("Shannon", "InvSimp", "Richness", "Abundance", "AlignmentRate")))
  
write.csv(Diversity_Table.long, file = "FullDipurity-Table1.csv")
MethodsToCompare = c("7","8","9","10","11","12", "006", "178", "208", "123") # Filter by these groups
getwd()
pd <- position_dodge(width = 0.7)
#Plot diversity estimates
Diversity_Table.long %>% 
  filter(ExpGroup %in% Treatments, 
         SampDigit %in% MethodsToCompare,
         Experiment == 2,
         Spikein == "UnSpiked",
         variable != "InvSimp") %>% #group_by(variable) %>% arrange(desc(value)) %>% select(-Experiment, SampDigit) %>% write_csv(.,path = "DiversityValues.csv")
  group_by(Extraction, ExpGroup, variable) %>% 
  summarise(MedValue = median(value),
            SdValue = sd(value)) %>% 
  ggplot(., aes(x = ExpGroup, y = MedValue, fill = Extraction)) +
  geom_bar(stat = "identity", position = pd, show.legend = T, color = "black") +
  geom_errorbar(position = pd, aes(ymin = MedValue - SdValue, ymax = MedValue + SdValue), 
                width = 0.3, 
                color = "black") +
  scale_fill_manual(name = "",
                    values = c("grey60","steelblue")) +
  labs(x = "",
       y = "Value") +
  facet_wrap(~variable, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Question 1 : SISPA bias and low yields
Diversity_Table.long %>% 
  filter(ExpGroup %in% c("Amicon", "Amicon+CsCl", "CsCl+Pellet"),
         #SampDigit %in% MethodsToCompare,
         #Experiment == 1,
         Spikein == "Spiked",
         #variable != "InvSimp",
         variable %in% c("Shannon", "InvSimp", "Richness", "AlignmentRate")) %>% 
  unite(Combo, Prep, Extraction, sep = "_", remove = F) %>% 
  mutate(Combo = factor(Combo, levels = c("Ovation_PS","SISPA-RP_PS","Ovation_QVR","SISPA-RP_QVR"))) %>% 
  #print(n = 100) 
  ggplot(., aes(x = ExpGroup, y = value, color = Combo)) +
  geom_boxplot(size = .3) +
  scale_color_manual(name = "",
                     values = c("grey30", "steelblue4","grey","steelblue1")) +
  facet_wrap(~variable, scale = "free") +
  ggthemes::theme_fivethirtyeight() + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(color = "grey60", size = 12))

## Here I contrasted the available prep_extraction permutations (Ovation vs SispaRP) available for Spiked samples to see if Sispa is actually more biased than other methods.

## The flipside in Unspiked samples...
Diversity_Table.long %>% 
  filter(ExpGroup %in% c("Amicon", "Amicon+CsCl", "CsCl+Pellet"),
         #SampDigit %in% MethodsToCompare,
         #Experiment == 1,
         Spikein == "UnSpiked",
         #variable != "InvSimp",
         variable %in% c("Shannon", "InvSimp", "Richness", "AlignmentRate")) %>% 
  unite(Combo, Prep, Extraction, sep = "_", remove = F) %>% 
  #mutate(Combo = factor(Combo, levels = c("Ovation_PS","SISPA-RP_PS","Ovation_QVR","SISPA-RP_QVR"))) %>% 
  mutate(Combo = factor(Combo, levels = c("Ovation_PS","SISPA-polyT_PS","Ovation_QVR","SISPA-polyT_QVR"))) %>% 
  #print(n = 100) 
  ggplot(., aes(x = ExpGroup, y = value, color = Combo)) +
  geom_boxplot(size = .3) +
  scale_color_manual(name = "",
                     values = c("grey30", "steelblue4","grey","steelblue1")) +
  facet_wrap(~variable, scale = "free") +
  ggthemes::theme_fivethirtyeight() + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(color = "grey60", size = 12))
theme_minimal()

## Question 2 & 3: Higher denature temperature increases the recovery of DS viruses with RP.
Diversity_Table.long %>% 
  filter(#ExpGroup %in% Treatments, 
    #SampDigit %in% MethodsToCompare,
    Extraction == "QVR",
    Experiment == 2,
    Spikein == "UnSpiked",
    variable %in% c("Shannon", "InvSimp", "Richness", "AlignmentRate")) %>%  filter(str_detect(ExpGroup, "Denat"))
  mutate(Prep = factor(Prep, levels = c("SISPA-polyT", "Ovation", "SISPA-RP")),
         ExpGroup = factor(ExpGroup, levels = c("Amicon+CsCl", "CsCl+Pellet", "Amicon", "TARA", "Denature", "Controlpos", "Controlneg"))) %>% 
  ggplot(., aes(y = value, x =  ExpGroup, color = Prep)) +
  geom_boxplot(show.legend = T) +
  scale_color_manual(name = "",
                     values = c("black", "#FFCE1B","steelblue")) +
  labs(x = "",
       y = "Value") +
  facet_wrap(~variable, scale = "free") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#Question 3 Highest richness in polyT-SISPA + CsCl gradient + QVR.
Diversity_Table.long %>% 
  filter(#ExpGroup %in% Treatments[1:3], 
    #SampDigit %in% MethodsToCompare,
    #Extraction == "QVR",
    Experiment == 2,
    Spikein == "UnSpiked",
    variable %in% c("Shannon", "InvSimp", "Richness", "AlignmentRate")) %>% # filter(str_detect(ExpGroup, "Denat"))
  mutate(Prep = factor(Prep, levels = c("SISPA-polyT", "Ovation", "SISPA-RP")),
         ExpGroup = factor(ExpGroup, levels = c("Amicon+CsCl", "CsCl+Pellet", "Amicon", "TARA", "Denature", "Controlpos", "Controlneg"))) %>% 
  ggplot(., aes(y = value, x =  ExpGroup, color = Prep, fill = Extraction)) +
  geom_boxplot(size = .3, show.legend = T) +
  scale_color_manual(name = "",
                     values = c("black", "#FFCE1B","steelblue")) +
  labs(x = "",
       y = "Value") +
  facet_wrap(~variable, scale = "free") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Question 4 PS vs QVR: “QVR consistently outperformed PS in all measurements and treatments.
Diversity_Table.long %>% #drop_na() %>% 
  filter(Spikein == "UnSpiked") %>% 
  ggplot(., aes(x = ExpGroup, y = value, fill = Extraction)) +
  geom_boxplot(outlier.size = .3, outlier.shape = 3, outlier.color = "red", size = .3) +
  scale_fill_manual(name = "",
                    values = c("grey60","steelblue")) +
  facet_wrap(~variable, scale = "free") +
  #theme_minimal() +
  ggthemes::theme_fivethirtyeight() + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(color = "grey60", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
## Question 5 - CsCl gradient on spikeins ###
Diversity_Table.long %>% filter() 
  filter(Spikein == "Spiked")

## Alignment rate stuff
AlignmentRate_DF = read.csv("~/Downloads/AlignmentRate.csv") %>% left_join(., SampleKey)
write.csv(AlignmentRate_DF, file = "AlignmentRates-Table1.csv")

AlignmentRate_DF %>% 
  drop_na() %>% 
  filter(Spikein == "UnSpiked",
         ExpGroup %in% Treatments) %>% view()
  group_by(Extraction, Prep, ExpGroup) %>% 
  summarise(MedRate = median(AlignmentRate),
            SdRate = sd(AlignmentRate)) %>% print(n=38) %>% 
  ggplot(., aes(x = ExpGroup, y = MedRate, fill = Extraction)) +
  geom_bar(stat = "identity", position = pd, show.legend = T, color = "black") +
  geom_errorbar(position = pd, aes(ymin = MedRate - SdRate, ymax = MedRate + SdRate), 
                width = 0.3, 
                color = "black") +
  scale_fill_manual(name = "",
                    values = c("grey60","steelblue")) +
  labs(x = "",
       y = "Percent Viral Reads") +
  facet_wrap(~Prep, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  
  
library(rstatix)
## Wilcox tables- Shannon
Shannon_pairPs = Diversity_Table.long %>% 
  filter(Extraction == "QVR"| Extraction == "PS",
         variable == "Shannon") %>% 
  split(.$Extraction) %>% 
  map(., ~pairwise_wilcox_test(value ~ ExpGroup, data = .)) #missing


##One plot of the old sample vOTU dataset (Find this otu table)
vOTU_table %>% 
  select(matches(regex("SISPA[789]"))) %>% 
  diversity(., index = "invsimpson", MARGIN = 2) %>% 
  data.frame("InvSimp" = .) %>% 
  rownames_to_column("Sample_ID") %>% left_join(.,vOTU_table %>% 
                                                  diversity(., index = "shannon", MARGIN = 2) %>% 
                                                  data.frame("Shannon" = .) %>% 
                                                  rownames_to_column("Sample_ID")) %>% 
  left_join(.,vOTU_table %>% rownames_to_column("ID") %>%  
              pivot_longer(cols = -ID, values_to = "values", names_to = "variables") %>% 
              arrange(variables) %>% mutate(bin = ifelse(values > 0,1,0)) %>% 
              group_by(variables) %>% 
              summarise("Richness" = sum(bin)) %>% 
              rename("Sample_ID" = variables)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"),
         ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = factor(SampDigit, levels = c("7","8","9")),
         ExpGroup = case_when(SampDigit == "7" ~"Amicon",
                              SampDigit == "8" ~"Amicon+CsCl",
                              SampDigit == "9" ~"CsCl+Pellet",
                              TRUE ~ SampDigit)) %>% 
  pivot_longer(cols = c(Richness, InvSimp,Shannon), names_to = "variable", values_to = "values") %>% 
  mutate(variable = factor(variable, levels = c("Shannon", "InvSimp", "Richness"))) %>% 
  ggplot(., aes(x = ExpGroup, y = values, fill = ExpGroup)) +
  geom_bar(stat = "identity", show.legend = F) +
  scale_fill_manual(name = "",
                    values = c("black","steelblue","grey")) +
  labs(x = "",
       y = "Value") +
  facet_wrap(~variable, scale = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
SampleBases = read.csv("SampleNumberBases.csv") %>% #
  filter(str_detect(Sample_ID, "(SISPA(7|8|9|10)")) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         SampDigit = factor(SampDigit, levels = c("7","8","9", "10")),
         SeqGroup = ifelse(str_detect(Sample_ID, "\\w+\\d[ADC]"), "New", "Old"))
SampleBases %>% 
  ggplot(., aes(x = SeqGroup, y = N.bases)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal()
NewOvations = c("")
### Round 3. Counting the number of bases from three previous rounds
# I would like to make apples to apples comparisons so the comparisons will be made among sample types, ie New Sispa to Old Sispa; New Ovation to old Ovation
# Define the samples for this, selecting only the samples from the old that have a New mate. For example, the Old set has more samples and types of sample prep than the New. And so filtering out the ones that arent present in both old and new—an inner_join() type of situation.

BasesAnalysisTable = read.csv("SampleNumberBases.csv") %>% 
  filter(str_detect(Sample_ID, regex("(sispa\\d)|(RNA[\\dC])", ignore_case = T))) %>% # get sispa and rna samples, nothing else.
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         SampDigit = factor(SampDigit, levels = c("7","8","9", "10", "11", "12")), #I want these as a symbol, not a digit, in the correct order.
         SeqGroup = ifelse(Target_bases > 4, "New", "Old"),
         PrepGroup = ifelse(str_detect(Sample_ID, "SISPA"), "SISPA", "Ovation")) %>% 
  unite(Group,SeqGroup, PrepGroup, remove = F) %>% 
  mutate(Group = factor(Group, levels = c("New_SISPA", "Old_SISPA", "New_Ovation", "Old_Ovation")))
#Manually curated the table to change RNA D and A to be Old.

BasesAnalysisTable %>% 
  ggplot(., aes(x=Group, y=Actual_bases, fill = factor(Target_bases))) +
  geom_boxplot(show.legend = T) +
  geom_jitter(width = .1, show.legend = F, alpha = .8, aes(color = factor(Target_bases))) +
  scale_fill_manual(name = "Target bases (Gbp)",
                    values = c("black","steelblue","grey")) +
  scale_color_manual(name = "Target bases (Gbp)",
                    values = c("black","steelblue","grey")) +
  labs(y = "Number of Raw Fastq Bases",
       x = "") +
  theme_minimal()

# As the rough model predicted, 5 Gbp was still much higher than the original (Old) samples
# This information is good enough to make a model


#save.image(file = "MethodsOptimizationReplicates.Rdata")
#load("MethodsOptimizationReplicates.Rdata")
##Species abundance curves
NormRNAv_OTU_DF %>% 
  pivot_longer(cols = everything(), names_to = "variables", values_to = "values") %>% 
  mutate(scaled_values = values * 10^19) %>% #pull(scaled_values) %>% summary()
  filter(!scaled_values < 1) %>% 
  

write.csv(NormRNAv_OTU_DF, file = "forgpt.csv")
##### FIGURES 

BasesAnalysisTable %>% 
  mutate()
