library(tidyverse)
library(karyoploteR)
library(ggpubr)

# genotípus
chooseGenotype = "P53KO"

# ---- 1. ADATBETÖLTÉS ----
setwd("//wsl$/Ubuntu/home/melinda/projects/cytodna/results/coverage") # mappakönyvtár
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
  scale_fill_manual(values = c("Mitokondriális" = "red", "Sejtmagi" = "dodgerblue")) +
  labs(
    title = "Mitokondriális és sejtmagi DNS arány mintánként",
    x = "Minta sorszám",
    y = "Arány (össz = 1)",
    fill = "DNS típusa"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")  # középre + nagyobb
  )

# ---- 3. PCA (QC – összes minta) ----
pca_input <- coverage_clean %>%
  mutate(sample_id = paste(genotype, treatment, replicate, sep="_")) %>% # mintainfo.-k összevonása új oszlopba
  select(chr, start, end, sample_id, coverage_fraction) %>% # oszlopok kiválasztása
  pivot_wider(names_from = sample_id, values_from = coverage_fraction) # tábla átalakítása

pca_input[is.na(pca_input)] <- 0 # NA coverage bin-ek kezelése
pca_matrix <- pca_input %>%
  select(-chr, -start, -end) %>%
  as.matrix() # mátrixszá konvertálás
pca_matrix <- t(pca_matrix) # mátrix transzponálása
pca_matrix <- pca_matrix[, apply(pca_matrix, 2, var) != 0] # PCA cleanup - null varianciájú oszlopok kezelése

pca_res <- prcomp(pca_matrix, scale. = TRUE) # PCA futtatása (mean=0, sd=1 standardizál)
# PCA eredmények kinyerése
pca_df <- as.data.frame(pca_res$x) # koordiná
pca_df$sample <- rownames(pca_df) # mintanevek visszarakása
pca_df <- pca_df %>%
  separate(sample, into=c("genotype","treatment","replicate"), sep="_") # visszabontás
# variancia kiszámítása, hogy PC1 és PC2 mekkora varianciát magyaráz
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
# PCA plot
ggplot(pca_df, aes(PC1, PC2, color=treatment, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (", round(var_explained[1]*100,1), "% magyarázott variancia)")) +
  ylab(paste0("PC2 (", round(var_explained[2]*100,1), "% magyarázott variancia)")) +
  labs(color="Kezelés", shape="Genotípus") +
  scale_color_discrete(labels=c(
    "Aphidicolin"="Afidikolin",
    "Aphidicolin-ATRi"= "Afidikolin-ATRi",
    "Aphidicolin-Ciszplatin" = "Afidikolin-Ciszplatin",
    "untreated" = "Kezeletlen"
  )) +
  theme_classic() +
  theme(
    legend.background = element_rect(color = "black", fill = "white"),
    #legend.box.background = element_rect(color = "black", fill = "white")
  )

# ---- 4. KARYOPLOT – QC (replikátumok) ----
plot_replicates <- function(data, chr=NULL, col1, col2, title=NULL){
  
  if(!is.null(chr)){
    chr_plot <- ifelse(chr=="MT","chrM",paste0("chr", chr))
    kp <- plotKaryotype(genome="hg38", chromosomes=chr_plot)
  } else {
    kp <- plotKaryotype(genome="hg38")
  }
  
  reps <- unique(data$replicate)
  
  df1 <- data %>% filter(replicate==reps[1])
  df2 <- data %>% filter(replicate==reps[2])
  
  # felső görbe
  kpLines(kp,
          chr = df1$chr_plot,
          x = (df1$start + df1$end)/2,
          y = df1$coverage_fraction,
          col = col1,
          r0 = 0.5, r1 = 1)
  
  # alsó görbe
  kpLines(kp,
          chr = df2$chr_plot,
          x = (df2$start + df2$end)/2,
          y = df2$coverage_fraction,
          col = col2,
          r0 = 0, r1 = 0.5)
  
  if(!is.null(title)){
    kpAddMainTitle(kp, title, cex=1)
  }
  legend("topright",
         inset = c(-0.11, -0.085),
         legend = c("Rep1", "Rep2"),
         col = c("black", "grey"),
         lty = 1,
         lwd = 3,
         seg.len = 2.75,
         bty = "n")
}

# untreated QC
untreated <- coverage_clean %>%
  filter(genotype==chooseGenotype, treatment=="untreated")

plot_replicates(untreated, col1="navy", col2="dodgerblue", title="P53KO-BRCA1KO kezeletlen replikátumok")
plot_replicates(untreated, chr="21", col1="navy", col2="dodgerblue")
plot_replicates(untreated, chr="9", col1="navy", col2="dodgerblue") ###???választani chr-t

# treated QC
treated <- coverage_clean %>%
  filter(genotype==chooseGenotype, treatment=="Aphidicolin-ATRi")

plot_replicates(treated, col1="black", col2="grey", title="P53KO-BRCA1KO Afidikolin-ATRi replikátumok")
plot_replicates(treated, chr="9", col1="black", col2="grey")
plot_replicates(treated, chr="21", col1="black", col2="grey")

# ---- 5. REPLIKÁTUMOK ÖSSZEVONÁSA ----
coverage_merged <- coverage_clean %>%
  group_by(chr, start, end, genotype, treatment) %>%
  summarise(read_count = sum(read_count), .groups="drop")

# ---- 6. RPKM NORMALIZÁLÁS ----
total_reads <- coverage_merged %>%
  group_by(genotype, treatment) %>%
  summarise(total_reads = sum(read_count), .groups="drop")

coverage_norm <- coverage_merged %>%
  left_join(y=total_reads, by=c("genotype","treatment")) %>%
  mutate(RPKM = (read_count * 1e9) / (total_reads * (end-start)))

coverage_norm <- coverage_norm %>%
  mutate(chr_plot = case_when(
    chr=="MT" ~ "chrM",
    TRUE ~ paste0("chr", chr)
  ))

# ---- 7. KARYOPLOT – ÖSSZEHASONLÍTÁS ----
plot_data_norm <- coverage_norm %>%
  filter(genotype==chooseGenotype,
         treatment %in% c("untreated","Aphidicolin-ATRi"))

kp <- plotKaryotype(genome="hg38")

# untreated
kpLines(kp,
        chr = plot_data_norm$chr_plot[plot_data_norm$treatment=="untreated"],
        x = (plot_data_norm$start + plot_data_norm$end)[plot_data_norm$treatment=="untreated"]/2,
        y = plot_data_norm$RPKM[plot_data_norm$treatment=="untreated"],
        col="dodgerblue",
        r0=0.5, r1=1)

# treated
kpLines(kp,
        chr = plot_data_norm$chr_plot[plot_data_norm$treatment=="Aphidicolin-ATRi"],
        x = (plot_data_norm$start + plot_data_norm$end)[plot_data_norm$treatment=="Aphidicolin-ATRi"]/2,
        y = plot_data_norm$RPKM[plot_data_norm$treatment=="Aphidicolin-ATRi"],
        col="black",
        r0=0, r1=0.5)

legend("topright",
       inset = c(-0.13, -0.085),
       legend = c("Kezeletlen", "Afidikolin-ATRi"),
       col = c("dodgerblue", "black"),
       lty = 1,
       lwd = 3,
       seg.len = 2.75,
       bty = "n")
kpAddMainTitle(kp, "P53KO-BRACA1 Kezeletlen és Afidikolin-ATRi", cex=1)

# ---- 8. KROMOSZÓMA-SPECIFIKUS PLOT ----
plot_chr <- function(chr_name){
  
  chr_plot <- ifelse(chr_name=="MT","chrM",paste0("chr",chr_name))
  
  kp <- plotKaryotype(genome="hg38", chromosomes=chr_plot)
  
  kpLines(kp,
          chr = plot_data_norm$chr_plot[plot_data_norm$treatment=="untreated"],
          x = (plot_data_norm$start + plot_data_norm$end)[plot_data_norm$treatment=="untreated"]/2,
          y = plot_data_norm$RPKM[plot_data_norm$treatment=="untreated"],
          col="dodgerblue",
          r0=0.5, r1=1)
  
  kpLines(kp,
          chr = plot_data_norm$chr_plot[plot_data_norm$treatment=="Aphidicolin-ATRi"],
          x = (plot_data_norm$start + plot_data_norm$end)[plot_data_norm$treatment=="Aphidicolin-ATRi"]/2,
          y = plot_data_norm$RPKM[plot_data_norm$treatment=="Aphidicolin-ATRi"],
          col="black",
          r0=0, r1=0.5)
}

plot_chr("3")

# ---- 9. ENRICHMENT (MAD) ----
comp <- plot_data_norm %>%
  select(chr, start, end, treatment, RPKM) %>%
  pivot_wider(names_from=treatment, values_from=RPKM)

comp[is.na(comp)] <- 0 # NA értékek nullázása

comp <- comp %>%
  mutate(diff = `Aphidicolin-ATRi` - untreated)

# MAD és threshold
mad_value <- mad(comp$diff)
threshold <- 3 * mad(comp$diff)

# enriched binek
enriched <- comp %>%
  filter(abs(diff) > threshold)
n_enriched <- nrow(enriched)

# label szöveg
label_text <- paste0(
  chooseGenotype,
  " kezeletlen és Afidikolin-ATRi"
)
stats_text <- paste0(
  "MAD: ", round(mad_value, 4), "\n",
  "Küszöbérték: ", round(threshold, 4), "\n",
  "Dúsult régiók: ", n_enriched
)

# ---- 10. ENRICHED RÉGIÓK PLOT ----
kp <- plotKaryotype(genome="hg38")

kpPoints(kp,
         chr = paste0("chr", enriched$chr),
         x = (enriched$start + enriched$end)/2,
         y = enriched$diff,
         col="red",
         cex=0.5)

kpAddMainTitle(kp, label_text, cex=1.3)
par(xpd=NA)
text(
  x = grconvertX(0.85, from="npc", to="user"),
  y = grconvertY(0.9,  from="npc", to="user"),
  labels = stats_text,
  adj = c(0,1),
  cex = 1.2
)

# Legnagyobb eltérést mutató régiók kiválasztása
top5 <- enriched %>%
  arrange(desc(abs(diff))) %>% # különbségek abszolút értékének csökkenő sorrendje
  head(50)
ggtexttable(top5)
