###################### SCRIPT INFO ######################
# This R_script corresponds to Fig. 1h and Fig. S3, data regarding editing without Kanamycin selection.
library(GenomicAlignments)
library(tidyverse)

function_bam <- function(bamFile){
  bn <- substr(basename(bamFile),1,nchar(basename(bamFile))-4)
  POS <- "KanR_site:745-752"
  stack <- stackStringsFromBam(bamFile, param=GRanges(POS))
  stack_table <- table(stack)
  st1 <- as.data.frame(stack_table)
  st2 <- cbind(st1, bn)
  colnames(st2) <- c("allele", "obs_count", "sample_name")
  assign( paste("stack_",bn, sep=""), st2, envir=.GlobalEnv ) }

bam_file_list <- list.files(path="~/Desktop/HIRIdata/20240705align/1/BAM", pattern="\\.bam$", full.names = TRUE)
lapply(bam_file_list, function_bam)
list_of_stacks1 <- grep( "^stack_", ls(), value=TRUE )
list_of_stacks2 <- mget(list_of_stacks1)
df_combined <- do.call(rbind, list_of_stacks2); row.names(df_combined) <- NULL
rm( list = list_of_stacks1 )

## Determine fraction of cell culture edited, unedited, ambiguous, at KanR site
df1 <- df_combined %>% filter(sample_name != "p167p31c4"); df1$allele <- factor( df1$allele, levels = c("ACAAGTCT", "GTAGGAAA", "AMB"), labels = c("ALT", "REF", "AMB") )

df2 <- df1 %>%
  separate(col = sample_name, into = c("blank","plasmid_editor","plasmid_guide","culture"), sep = "p|c" ) %>%
  select(-c(blank)) %>%
  mutate(ID = paste(plasmid_editor, plasmid_guide, sep="x"), .keep="unused") %>%
  group_by(ID, culture, allele) %>%  summarize(obs_count = sum(obs_count)) %>% ungroup() %>%
  mutate( allele = as.character(allele),
          allele = ifelse(is.na(allele), "AMB", allele),
          allele = factor(allele, levels=c("REF","ALT","AMB")),
          ID = factor(ID, levels=c("167x32","167x31","168x32","168x31")),
          culture = factor(culture)) %>%
  complete(ID, culture, allele, fill = list(obs_count = 0)) %>%
  group_by(ID, culture) %>% mutate(frac = obs_count / sum(obs_count)) %>% ungroup()

df3 <- df2 %>%
  group_by(ID,allele) %>%
  summarize(frac_mean = mean(frac),
            frac_sd = sd(frac))

plotA <- ggplot() +
  geom_col( data=df3, aes(x=ID, y=frac_mean, fill=allele), position=position_dodge(width=1) ) +
  geom_errorbar( data=df3, aes(x=ID, ymin=(frac_mean-frac_sd), ymax=(frac_mean+frac_sd), fill=allele ), position=position_dodge(width=1), width=0.25 ) +
  geom_point(data=df2, aes(x=ID, y=frac, fill=allele), position=position_dodge(width=1) )

plotB <- plotA +
  scale_y_continuous( limits=c(0,1), breaks=seq(0,1, by= 0.1) ) +
  labs(y="Fraction of all sequencing reads",
       title="Non-selective editing at the KanR* locus") +
  scale_x_discrete(labels=c(
    "167x31"="DarT(G49D)-\nnCas9 \nT sgRNA",
    "167x32"="DarT(G49D)-\nnCas9 \nNT sgRNA",
    "168x31"="dDarT(D)-\nnCas9 \nT sgRNA",
    "168x32"="dDarT(D)-\nnCas9 \nNT sgRNA") ) +
  scale_fill_discrete( labels=c(
    "ALT"="Edited",
    "AMB"="Ambiguous",
    "REF"="Unedited") ) +
  theme_minimal() ; plotB


### Determine base mutation percentage at the ADPr-site.
df4 <- df_combined %>% 
  filter(sample_name != "p167p31c4") %>% 
  mutate(allele = as.character(allele)) %>%
  filter(str_detect(allele, "GTAGG[CGT]AA")) %>%
  mutate(allele = substr(allele,6,6)) %>%
  separate(sample_name,into = c("blank","plasmid_editor","plasmid_guide","culture"), sep = "p|c"  ) %>%
  select(-c("blank")) %>%
  mutate(ID = paste(plasmid_editor, plasmid_guide, sep = "x"), .keep = "unused")

df5 <- df2 %>% group_by(ID, culture) %>% summarize(obs_total = sum(obs_count))

df6 <- merge.data.frame(df4, df5, by = c("ID", "culture")) %>% 
  mutate(percent = 100 * obs_count / obs_total) %>% 
  mutate(ID = factor(ID, levels = c("167x32", "167x31", "168x32", "168x31") ),
         allele = factor(allele, levels = c("T","G","C"), labels = c("A","C","G") ))

df7 <- df6 %>% complete(ID, culture, allele, fill = list(percent=0))

plotC <- ggplot(data=df7, aes(x=culture, y=allele, fill=percent))+
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "#449EDF", limits=c(0,0.3), breaks=seq(0,0.3, by=0.1))+
  geom_tile( color = "black", lwd = 0.5, linetype = 1)+
  facet_wrap(~ID, labeller = labeller(ID = c(
    "167x31"="DarT(G49D)-\nnCas9 \nT sgRNA",
    "167x32"="DarT(G49D)-\nnCas9 \nNT sgRNA",
    "168x31"="dDarT(D)-\nnCas9 \nT sgRNA",
    "168x32"="dDarT(D)-\nnCas9 \nNT sgRNA")) )+
  guides(fill = guide_colourbar(barwidth = 1, barheight = 15))+
  theme( axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12));plotC
###################### END OF SCRIPT ######################