library(tidyverse)

traits <- readxl::read_xlsx("fungal traits.xlsx")

tax$genus[1:10]
tax$primary_lifestyle <- traits[match(tax$genus, traits$GENUS), 8]
tax$secondary_lifestyle <- traits[match(tax$genus, traits$GENUS), 9]
tax <- as.data.frame(tax)

write.csv(tax, "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2020.csv")

write.csv(tax, "data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_FUNGI_2021.csv")


write.csv(tax, "other data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_2020.csv", quote = F, row.names = F)

write.csv(tax, "other data tables/mm_Uv9_ssu-its-lsu_SR_taxonomy_2021.csv", quote = F, row.names = F)

unique(tax$primary_lifestyle)


# to reassign secondary to primary

reassign <- which(tax$primary_lifestyle == "animal_parasite")
newnames <- tax$secondary_lifestyle[reassign]
tax$primary_lifestyle[reassign] <- newnames

