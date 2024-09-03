###################### SCRIPT INFO ######################
# This R_script corresponds to Fig. 1i, regarding WGS profiling of off-target SNVs.
library(tidyverse)

# Load csv file.
df_raw <- read_csv("../data/Input_files_Script5/OUT_CONCAT.csv", col_types = "ffdccdd") 

# Extract SNPs.
df1 <- filter(df_raw, grepl("^.$", df_raw$REF) & grepl("^.$", df_raw$ALT))

# Remove POS already present in the strain, absent plasmids.
list_p0 <- df1$POS[df1$Plasmid==0]
df2 <- df1 %>% filter(!POS %in% list_p0)

# Remove POS regions manipulated to create strain HBs30.
list_manipulation <- data.frame(start = c(1986236, 364749, 65780,890460), stop = c(1986449, 364750, 71351,890500),reason = c("araF","lacZ","araABCD","KanR_insert") )
df3 <- df2 %>% filter(!sapply(POS, function(region) { any( sapply( seq_len(nrow(list_manipulation)), function(i) {region >= list_manipulation$start[i] && region <= list_manipulation$stop[i]} ) ) } ) ) 

# Remove records, where AF < 0.25
df4 <- df3 %>% filter(AF>=0.25)

# Determine lowest Q1-DP (34), then remove records where DP is less than 34, unless other records at that POS are greater than 34.
RAW_DP <- read_csv("~/Desktop/HIRIdata/20240628align/2/List_DP.csv", col_types = c("f","i"))
DP_Q1 <- RAW_DP %>% group_by(FILE) %>% summarize(Q1= quantile(DP, 0.25)) %>% filter(FILE != "p0c0")
min_Q1 <- min(DP_Q1$Q1); print(min_Q1)
rm(RAW_DP, DP_Q1)
df5 <- df4 %>% filter(DP >= min_Q1 | POS %in% df4$POS[df4$DP>=min_Q1][which(df4$POS[df4$DP>=min_Q1] %in% df4$POS[df4$DP<min_Q1])])

# Establish SNV type.
df6 <- df5 %>% mutate(SNV = paste(REF,ALT, sep=""), SNV = case_match(SNV, "AC" ~ "TG","AG" ~ "TC","AT" ~ "TA","GT" ~ "CA","GC" ~ "CG","GA" ~ "CT", .default=SNV), .keep="unused") %>% 
  filter(Plasmid %in% c(11,167,190)) %>% mutate(Plasmid=factor(Plasmid, levels=c(11,167,190)), Colony=factor(Colony, levels=c(1,2,3)), SNV=factor(SNV, levels=c("TA","TC","TG","CA","CG","CT")))

# Tally SNVs
df7 <- df6 %>% count(Plasmid, Colony, SNV) %>% complete(Plasmid, Colony, SNV, fill=list(n=0))

# Generate summary data
df8 <- df7 %>% group_by(Plasmid, SNV) %>%
  summarize(n_mean = mean(n),
         n_sd = sd(n))

## Heatmap of individual samples
plotA <- ggplot(data=df7, aes(x=Colony, y=SNV, fill=n)) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#449EDF", limits=c(0,10)) +
  facet_grid(~ Plasmid, labeller = labeller( Plasmid = c("11"="nCas9", "167"="DarT(G49D)-nCas9", "190"="APOBEC-nCas9-UGI") ) )+
  geom_tile( color = "black", lwd = 0.5, linetype = 1) +
  geom_text(aes(label = n), color = "black", size = 6) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15))+
  theme( axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))+
  theme_minimal()
  


## Add motif information. Each POS was manually searched against NCBI Reference Sequence NC_000913.3 to extract the motif, with the SNV centered on character-3.
library(gt)
paste_motif <- strsplit("ATTG,ATCG,AGTG,GCTG,GTCG,ATTG,TCCA,GTCG,CTCA,TTTC,TTTC,ATCG,TCCA,CTCA,TGTG,GTCG,TGTG,CTCT,CTCT,CTCC,ATTG,GCTG,CTCC,GTCC,TTCC,GTCA,AATG,GTCG,TTCA,CGCC,ATCG,AATG,TTCC,TTCC,TTCC", split = ",")
df10 <- cbind(df6, paste_motif); colnames(df10)[7] <- "Motif"
df11 <- df10 %>% mutate(Editor = factor(Plasmid, levels=c(11,167,190),labels=c("nCas9","DarT(G49D)-nCas9","APOBEC-nCas9")),
                        Colony = Colony,
                        "Position" = POS,
                        "Read Depth" = DP,
                        "Frequency" = AF,
                        SNV = str_c(str_sub(SNV, 1, 1), ">", str_sub(SNV, 2, 2)), .keep="unused") %>% 
  select(c("Editor", "SNV", "Motif", "Colony", "Position", "Read Depth", "Frequency")) %>% arrange(Editor, SNV)

tableA <- gt(df11) ; tableA
###################### END OF SCRIPT ######################