###################### SCRIPT INFO ######################
# This R_script corresponds to  Fig. S5, regarding growth curves of additional point mutations to DarT(D) in the delta-recA strain.

library(tidyverse)

## Load data ====
df_raw <- read_csv("../data/Input_files_Script6/20230826_toxicity_assay.csv")

df1 <- df_raw
df_names1 <- colnames(df1)
df_names2 <- sub("r.*", "", df_names1)
df_names3 <- sub("p","",df_names2); df_names3[1] <- c("Time")

## Plots of Growth curves. ====
df4 <- df1; colnames(df4) <- df_names3
BLK_mean <- mean(unlist(df4[1,colnames(df4) == "BLK"]))
df5 <- df4 %>% pivot_longer(cols = -c(Time), names_to = "Plasmid", values_to = "OD") %>% mutate(ODc = OD - BLK_mean) %>% filter(str_detect(Plasmid, "[0-9]"))
df6 <- df5 %>% group_by(Time, Plasmid) %>% summarize(OD_mean = mean(ODc), OD_sd = sd(ODc)) %>% ungroup() %>% mutate(Hours = as.numeric(Time) / 3600, Plasmid=factor(Plasmid, levels = c(167,195,182,180,181,196,21,28,168)))
custom_labels <- c("167" = "G49D", "195"="G49D, R57A", "182"="G49D, M84L", "180"="G49D, M86L", "181"="G49D, R92A", "196"="G49D, R166A", "21"="G49D, R193A", "28"="G49D, M86L, R92A, R193A", "168"="G49D, E170A")
plotA <- ggplot(data = df6, aes(x=Hours, y=OD_mean)) +
  geom_point(color="black") +
  geom_ribbon(aes(ymin = OD_mean - OD_sd, ymax = OD_mean + OD_sd), alpha = 0.25) +
  facet_wrap(vars(Plasmid), labeller = labeller(Plasmid = custom_labels)) +
  scale_x_continuous(breaks = seq(0,12, by= 2)) +
  scale_y_continuous(lim = c(0,0.6), breaks = seq(0, 0.6, by = 0.1)) +
  labs(y= "OD600") +
  theme(axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        strip.background = element_blank()) ;plotA

## Plot final OD for each sample. ====
Final_Time <- as.numeric(df1[nrow(df1),"Time"])
df8 <- df5 %>% filter(Time == Final_Time)
df9 <- df8 %>% group_by(Plasmid) %>% summarize(OD_mean = mean(ODc), OD_sd = sd(ODc)) %>% mutate(Plasmid = factor(Plasmid, levels = c(167,195,182,180,181,196,21,28,168)))
plotB <- ggplot() +
  geom_col(data=df9, aes(x=Plasmid, y=OD_mean), width=.5) +
  geom_errorbar(data=df9, aes(x=Plasmid, ymin=(OD_mean - OD_sd), ymax=(OD_mean + OD_sd)), width=0.25 ) +
  geom_point(data=df8, aes(x=Plasmid, y=ODc), position=position_jitter(), size=5 ) +
  labs(y="Final growth (OD600)") +
  scale_y_continuous(limits=c(0,0.6))+
  scale_x_discrete(labels = custom_labels)+
  theme_minimal();plotB

# Save plots ====
output_folder <- "../outputs/Output_files_Script6"
dir.create(output_folder, recursive = TRUE)

filename <- paste0(output_folder, "/PlotA_growth_curves.pdf")
ggplot2::ggsave(
  filename = filename,
  plot     = plotA,
  width    = 10,
  height   = 7
)

filename <- paste0(output_folder, "/PlotB_final_OD.pdf")
ggplot2::ggsave(
  filename = filename,
  plot     = plotB,
  width    = 12,
  height   = 7
)
###################### END OF SCRIPT ######################