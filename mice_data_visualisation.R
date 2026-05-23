library(tidyverse)
library(karyoploteR)

# ---- 1. ADATBETÖLTÉS ----
setwd("//wsl$/Ubuntu/home/melinda/projects/mice_test_dataset/results/coverage") # mappakönyvtár
files <- list.files(pattern="coverage.txt$") # lista a fájlok neveivel
coverage_data <- map_df(files, function(f) {
  df <- read_tsv(f, col_names = FALSE)
  colnames(df) <- c(
    "chr","start","end",
    "read_count","bases_covered",
    "bin_size","coverage_fraction"
  )
  df$sample <- f # új oszlop ($) beszúrása a fájl nevével
  df
})

# ---- 2. TISZTÍTÁS ----
coverage_clean <- coverage_data %>%
  filter(chr %in% c(as.character(1:22), "X", "Y", "MT")) %>%
  mutate(sample_clean = str_remove(sample, "_100kbp_coverage.txt")) %>%
  separate(sample_clean,
           into = c("genotype","treatment","replicate"),
           sep = "_",
           fill = "right")

# chr formátum plothoz
coverage_clean <- coverage_clean %>%
  mutate(chr_plot = case_when(
    chr == "MT" ~ "chrM",
    TRUE ~ paste0("chr", chr)
  ))

# ---- MT vs NUCLEAR ARÁNY ----
dna_plot <- coverage_clean %>%
  mutate(type = ifelse(chr == "MT", "mt", "nuclear")) %>%
  group_by(sample, type) %>%
  summarise(reads = sum(read_count), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(total = sum(reads),
         fraction = reads / total) %>%
  ungroup()

# százalék számítás felirathoz
dna_plot <- dna_plot %>%
  mutate(percent = round(fraction * 100, 1))

# minták sorszámozása
dna_plot <- dna_plot %>%
  mutate(sample_id = as.numeric(factor(sample)))
# legend beállítása
dna_plot <- dna_plot %>%
  mutate(type = recode(type,
                       "mt" = "Mitokondriális",
                       "nuclear" = "Sejtmagi"))

# plot
ggplot(dna_plot, aes(x = factor(sample_id), y = fraction, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  # geom_text(
  #   data = dna_plot %>% filter(type == "Sejtmagi"),
  #   aes(label = paste0(percent, "%")),
  #   position = position_stack(vjust = 0.5),
  #   size = 3
  # )+
  scale_fill_manual(values = c("Mitokondriális" = "red", "Sejtmagi" = "dodgerblue")) +
  scale_x_discrete(expand = c(0,01, 3))+
  labs(
  # title = "Mitokondriális és sejtmagi DNS arány mintánként",
    x = "Minta sorszám",
    y = "Arány (össz = 1)",
    fill = "DNS típusa"
  ) +
  theme_classic()
  # theme(
  #   plot.title = element_text(hjust = 0.5, size = 16, face = "bold")  # középre + nagyobb
  # )
#ggsave("mt_nuclear_ratio.png", width = 4, height = 4)
