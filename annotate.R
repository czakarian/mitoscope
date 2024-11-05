## Annotating mutserve output with MITOMAP

## generalize to accept a vcf since annos done by postion, REF, ALT (to accept Mutect2 input)

library(tidyverse)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)
input_file=args[1]
anno_file=args[2]

input_prefix = file.path(dirname(input_file), gsub("\\.txt$", "", basename(input_file)))

## read in sample file and anno file
mutserve_input = read.delim(input_file)
anno_df = read.delim(anno_file, sep=",")

## merge anno df with input df 
input_anno = mutserve_input %>% left_join(anno_df, by = c("Pos" = "Position", "Ref" = "REF", "Variant" = "ALT"))

input_anno = input_anno %>%
  mutate(DiseaseVariantStatus = case_when(
    Source %in% c("variant_df", "disease_df,variant_df") ~ "Common Variant", 
    Source == "disease_df" ~ "Disease Variant",
    is.na(Source) ~ "Unknown Variant")
  )

write_delim(input_anno, paste0(input_prefix, ".annotated.txt"), delim="\t")

## generate annotated heteroplasmy plot
heteroplasmy_plot = ggplot(input_anno, mapping = aes(Pos, VariantLevel)) +
  geom_point(aes(color = DiseaseVariantStatus), size = 6, alpha = 0.7) +
  theme_bw() + theme(panel.grid = element_blank(),
    axis.title = element_text(size=24),
    axis.text = element_text(size=20),
    legend.text = element_text(size=18),
    legend.title = element_blank()) +
#  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0), limits = c(-0.1,1.1)) +
  scale_x_continuous(limits = c(0, 16569), breaks = (c(0, 4000, 8000, 12000, 16000))) +
  scale_color_manual(values = c("Common Variant" = "#00BFC4", "Disease Variant" = "#F8766D", "Unknown Variant" = "grey")) +
  ylab("Variant Level") + xlab("Position in mtDNA genome") + 
  geom_text_repel(data = subset(input_anno, Source %in% c("disease_df", NA)),
     aes(label = paste(Ref, Pos, Variant, "\n", VariantLevel*100, "%", sep="")), size = 7, color = "black", nudge_y = 0.08)

ggsave(paste0(input_prefix, ".heteroplasmy.png"), heteroplasmy_plot, width = 12, height= 8)



