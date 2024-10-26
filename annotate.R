## Annotating mutserve output with MITOMAP
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

## read in sample file
input_file=args[1]
output_file=args[2]
anno_dir=args[3]

mutserve_input = read.delim(input_file)

### can consider saving the combined df so dont have to run this part every time??
## read in MITOMAP annos
mitomap_disease_codingcontrol = read.delim(file=paste0(anno_dir, "/MutationsCodingControl_MITOMAP.csv"), sep=",")
mitomap_disease_RNA = read.delim(file=paste0(anno_dir, "/MutationsRNA_MITOMAP.csv"), sep=",")
mitomap_variants_codingRNA = read.delim(file=paste0(anno_dir, "/VariantsCoding_MITOMAP.csv"), sep=",")
mitomap_variants_control = read.delim(file=paste0(anno_dir, "/VariantsControl_MITOMAP.csv"), sep=",")

protein_coding_genes = c("MT-ATP6", "MT-ATP8", "MT-ATP8/6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6")
rRNA_genes = c("MT-RNR1", "MT-RNR2", "MT-RNR3")
tRNA_genes = c("MT-TA", "MT-TR", "MT-TN", "MT-TD", "MT-TC", "MT-TE", "MT-TQ", "MT-TG", "MT-TH", "MT-TI", "MT-TL1", "MT-TL2", "MT-TK", "MT-TM", "MT-TF", "MT-TP", "MT-TS1", "MT-TS2", "MT-TT", "MT-TW", "MT-TY", "MT-TV", "MT-TS1 precursor")
control_region = c("MT-CR")
noncoding_regions = c("MT-OLR", "MT-TER", "MT-NC1", "MT-NC2", "MT-NC3", "MT-NC4", "MT-NC5", "MT-NC6", "MT-NC7", "MT-NC8", "MT-NC9", "MT-NC10")
peptides = c("MT-Hum", "MT-MOTSc", "MT-SHLP1", "MT-SHLP2", "MT-SHLP3", "MT-SHLP4", "MT-SHLP5", "MT-SHLP6")
valid_gene_labels = c(protein_coding_genes, rRNA_genes, tRNA_genes, control_region)

## standardise the individual dfs before merging
mitomap_disease_codingcontrol = mitomap_disease_codingcontrol %>%
  select(-c(Allele)) %>%
  separate(NucleotideChange, into = c("REF", "ALT"), sep = "-") %>%
  separate(`Plasmy.Reports.Homo.Hetero.`, into = c("Homoplasmy", "Heteroplasmy"), sep = "/") %>%
  mutate(Amino.Acid.Change =`Amino.AcidChange`) %>% select(-c(`Amino.AcidChange`)) %>%
  separate(`GB.Freq..FL..CR...`, into = c("GB.Freq.FL", "GB.Freq.CR"), sep = "\\(") %>%
  mutate(GB.Freq.CR = str_replace(GB.Freq.CR, "\\)", "")) %>% 
  separate(`GB.Seqs.FL..CR..`, into = c("GB.Seqs.FL", "GB.Seqs.CR"), sep = "\\(", convert=TRUE) %>%
  mutate(GB.Seqs.CR = as.integer(str_replace(GB.Seqs.CR, "\\)", ""))) %>%
  mutate(Gene.Name = Locus) %>% select(-Locus) %>%
  # If Gene.Name is in valid_gene_labels, keep it, otherwise move it to Additional.Annotations
  mutate(Additional.Annotations = if_else(!(Gene.Name %in% valid_gene_labels), Gene.Name, NA),
         Gene.Name = if_else(Gene.Name %in% valid_gene_labels, Gene.Name, NA))

mitomap_disease_RNA = mitomap_disease_RNA %>%
  separate(Allele, into = c("REF", "ALT"), sep = "\\d+") %>%
  mutate(MitoTIP = `MitoTIP.`, Amino.Acid.Change = RNA) %>% select(-c(RNA, `MitoTIP.`)) %>%
  separate(`GB.Freq..FL..CR...`, into = c("GB.Freq.FL", "GB.Freq.CR"), sep = "\\(") %>%
  mutate(GB.Freq.CR = str_replace(GB.Freq.CR, "\\)", "")) %>% 
  separate(`GB.Seqs.FL..CR..`, into = c("GB.Seqs.FL", "GB.Seqs.CR"), sep = "\\(", convert=TRUE) %>%
  mutate(GB.Seqs.CR = as.integer(str_replace(GB.Seqs.CR, "\\)", ""))) %>%
  mutate(Gene.Name = Locus) %>% select(-Locus) %>%
  mutate(Additional.Annotations = if_else(!(Gene.Name %in% valid_gene_labels), Gene.Name, NA),
         Gene.Name = if_else(Gene.Name %in% valid_gene_labels, Gene.Name, NA))

mitomap_variants_codingRNA = mitomap_variants_codingRNA %>%
  select(-c(Codon.Number, Codon.Position)) %>%
  separate(Nucleotide.Change, into = c("REF", "ALT"), sep = "-") %>%
  mutate(GB.Freq.FL = `GB.Freq.`, GB.Seqs.FL = as.integer(GB.Seqs), GB.Freq.CR = "0.000%", GB.Seqs.CR = 0, References = `Curated.References`) %>%
  select(-c(`GB.Freq.`, `GB.Seqs`, `Curated.References`))  %>%
  mutate(Gene.Name = Locus) %>% select(-Locus) %>%
  mutate(Additional.Annotations = if_else(!(Gene.Name %in% valid_gene_labels), Gene.Name, NA),
         Gene.Name = if_else(Gene.Name %in% valid_gene_labels, Gene.Name, NA))
  
mitomap_variants_control = mitomap_variants_control  %>%
  separate(Nucleotide.Change, into = c("REF", "ALT"), sep = "-") %>%
  mutate(References = `Curated.References`) %>% select(-c(`Curated.References`)) %>%
  separate(`GB.FreqFL..CR...`, into = c("GB.Freq.FL", "GB.Freq.CR"), sep = "\\(") %>%
  mutate(GB.Freq.CR = str_replace(GB.Freq.CR, "\\)", "")) %>% 
  mutate(`GB.Seqstotal..FL.CR..` = if_else(`GB.Seqstotal..FL.CR..` == "0", "0 (0/0)" , `GB.Seqstotal..FL.CR..`)) %>%
  separate(`GB.Seqstotal..FL.CR..`, into = c(NA, "GB.Seqs.Ind"), sep = "\\(") %>%
  separate(GB.Seqs.Ind, into = c("GB.Seqs.FL", "GB.Seqs.CR"), sep = "\\/", convert=TRUE) %>%
  mutate(GB.Seqs.CR = as.integer(str_replace(GB.Seqs.CR, "\\)", "")))  %>%
  mutate(Amino.Acid.Change = "noncoding") %>%
  mutate(Locus = if_else(Locus == "Control Region", "MT-CR", Locus)) %>%
  mutate(Gene.Name = Locus, Gene.Type = "control region") %>% select(-Locus) %>%
  mutate(Additional.Annotations = if_else(!(Gene.Name %in% valid_gene_labels), Gene.Name, NA),
         Gene.Name = "MT-CR")

## join dfs
disease_df = bind_rows(mitomap_disease_codingcontrol, mitomap_disease_RNA) %>%
  mutate(Gene.Type = case_when(
           Gene.Name %in% control_region ~ "control region",
           Gene.Name %in% protein_coding_genes ~ "protein coding",
           Gene.Name %in% rRNA_genes ~ "rRNA",
           Gene.Name %in% tRNA_genes ~ "tRNA",
           Gene.Name %in% noncoding_regions ~ "noncoding",
           Gene.Name %in% peptides ~ "small peptide",
           TRUE ~ "Unknown"),
         GB.Freq.FL = str_extract(GB.Freq.FL, "^[^%]*%"),
         GB.Freq.CR = str_extract(GB.Freq.CR, "^[^%]*%"))

variants_df = bind_rows(mitomap_variants_codingRNA, mitomap_variants_control) %>%
  group_by(Position, REF, ALT) %>% summarise(across(everything(), ~paste(unique(na.omit(.)), collapse = ",")), .groups = 'drop') %>%
  mutate(Gene.Name = if_else(Gene.Name == "", Additional.Annotations, Gene.Name)) %>%
  mutate(Gene.Type = case_when(
          str_detect(Gene.Name, str_c(control_region, collapse = "|")) ~ "control region",
          str_detect(Gene.Name, str_c(protein_coding_genes, collapse = "|")) ~ "protein coding",
          str_detect(Gene.Name, str_c(rRNA_genes, collapse = "|")) ~ "rRNA",
          str_detect(Gene.Name, str_c(tRNA_genes, collapse = "|")) ~ "tRNA",
          str_detect(Gene.Name, str_c(noncoding_regions, collapse = "|")) ~ "noncoding",
          str_detect(Gene.Name, str_c(peptides, collapse = "|")) ~ "small peptide",
          TRUE ~ "Unknown")) %>%
  mutate(GB.Seqs.FL = as.integer(GB.Seqs.FL),
         GB.Seqs.CR = as.integer(GB.Seqs.CR),
         GB.Freq.FL = str_extract(GB.Freq.FL, "^[^%]*%"),
         GB.Freq.CR = str_extract(GB.Freq.CR, "^[^%]*%"), 
         References = as.integer(References))

combined_df <- bind_rows(
  disease_df %>%
    mutate(Source = "disease_df") %>% 
    select(Position, REF, ALT, Gene.Name, Gene.Type, Amino.Acid.Change, GB.Freq.FL, GB.Freq.CR, GB.Seqs.FL, GB.Seqs.CR, everything()),

  variants_df %>%
    mutate(Source = "variant_df") %>%  # Add a source column for variants
    select(Position, REF, ALT, Gene.Name, Gene.Type, Amino.Acid.Change, GB.Freq.FL, GB.Freq.CR, GB.Seqs.FL, GB.Seqs.CR, everything())
)

collapsed_df = combined_df %>%
  group_by(Position, REF, ALT) %>%
  summarise(across(everything(), ~paste(unique(na.omit(.)), collapse = ",")), .groups = 'drop')


## merge anno df with input df 
input_anno = mutserve_input %>% left_join(collapsed_df, by = c("Pos" = "Position", "Ref" = "REF", "Variant" = "ALT"))
write_delim(input_anno, output_file, delim="\t")


## save heteroplasmy plot to a file

# ggplot(input_anno, mapping = aes(Pos, VariantLevel)) +
#   geom_point(aes(color = Source), size = 3, alpha = 0.6) +
#   theme_bw() + theme(panel.grid = element_blank()) +
#   geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(-0.1,1.1)) +
#   scale_x_continuous(limits = c(0, 16569), breaks = (c(0, 5000, 10000, 15000))) +
#   ylab("Heteroplasmy") + xlab("Position in mtDNA genome")





