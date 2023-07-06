library(readxl)
library(dplyr)


### GO sets
KPC <- read_excel("C:/Users/edmondsonef/Desktop/KPC DSP GENE LIST.xlsx", sheet = "KPC GENE LIST 07.06.22_comps_MH")
KPC <- read_excel("C:/Users/edmondsonef/Desktop/KPC DSP GENE LIST.xlsx", sheet = "Significant")
#CosMX <- read_excel("C:/Users/edmondsonef/Desktop/KPC DSP GENE LIST.xlsx", sheet = "CosMX mouse neuro")
#axon <- read_excel("C:/Users/edmondsonef/Desktop/GO_term_summary_20230601_102507.xlsx")
FDA <- read_excel("C:/Users/edmondsonef/Desktop/FDA Drug Targets.xlsx")

#axon$SYMBOL <- axon$Symbol
#axon <- distinct(axon, SYMBOL, .keep_all= TRUE)
new <- dplyr::left_join(KPC, CosMX, by = "SYMBOL")

write.csv(new, "C:/Users/edmondsonef/Desktop/data.csv")


#### HWANG Paper
KPC <- read_excel("C:/Users/edmondsonef/Desktop/KPC DSP GENE LIST.xlsx", sheet = "KPC GENE LIST 07.06.22_comps_MH")
Hwang <- read_excel("C:/Users/edmondsonef/Desktop/KPC DSP GENE LIST.xlsx", sheet = "Hwang lists")

KPC$SYMBOL_L <- tolower(KPC$SYMBOL)
FDA$SYMBOL_L <- tolower(FDA$SYMBOL)

#SYMBOL <- Hwang$`Neural-like progenitor`
#SYMBOL <- Hwang$`Acinar-like`
#SYMBOL <- Hwang$`Classical-like`
#SYMBOL <- Hwang$Basaloid
#SYMBOL <- Hwang$Squamoid
#SYMBOL <- Hwang$Mesenchymal
SYMBOL <- Hwang$`Neuroendocrine-like`

LIST <- data.frame(SYMBOL)

LIST$SYMBOL_L <- tolower(LIST$SYMBOL)
rm(new)
new <- dplyr::left_join(KPC, FDA, by = "SYMBOL_L")

write.csv(new, "C:/Users/edmondsonef/Desktop/data.csv")



#Remove Duplicates
new <- distinct(KPC, SYMBOL, .keep_all= TRUE)
