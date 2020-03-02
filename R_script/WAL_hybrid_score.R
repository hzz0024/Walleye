################################################################################################################################carp
######## Scoring genotypes ratios for NEWHYBRIDS evaluation ########

# The purpose of this R script is to score three genotype ratios (homozygous AA and BB, heterozygous AB) for each examined fish.
# From your genotype file, copy your genotypes and transpose (lengthwise) them into a csv file containing the sample IDs and the 68 SNPs
# Follow the instructions below to run the script with your data and the output will be a file with the scored markers and percentages

## Part 1: Set working directory ##
# e.g., setwd("Your working directory")

## Part 2: Read into the reference and your sample files ##
# Make sure to have the "Reference Genotype" file on your working directory
# Change the second variable to reflect the correct file and name of your data
# Read into the reference and your sample files

refsnp <- read.csv(file="example_ref.csv")
genotype <- read.csv(file="example_genotypes.csv", na.string="")

## Part 3: Change the length of your samples ##
# This number should reflect the number of samples you will be running

length_samples <- 60
out <- vector(mode="list", length=length_samples)

## Part 4: For loop to calculate AA/BB/AB ##

for (i in 1:length_samples) {
  
  AA_cnt <- sum(as.character(refsnp[,3])==as.character(genotype[,i]), na.rm=TRUE)
  BB_cnt <- sum(as.character(refsnp[,2])==as.character(genotype[,i]), na.rm=TRUE)
  total <- sum(!is.na(genotype[,i]))
  he1 <- paste(refsnp[,2], refsnp[,3], sep="")
  he2 <- paste(refsnp[,3], refsnp[,2], sep="")
  he1calc <- sum(he1==as.character(genotype[,i]), na.rm=TRUE)
  he2calc <- sum(he2==as.character(genotype[,i]), na.rm=TRUE)
  hetot <- sum(he1calc+he2calc, na.rm=TRUE)
  
  A_freq <- ((AA_cnt+(0.5*hetot))/total)
  B_freq <- ((BB_cnt+(0.5*hetot))/total)
  HE <- (hetot/total)
  HOMO <- (1-HE)
  AA_freq <- (AA_cnt/total)
  BB_freq <- (BB_cnt/total)
  H_AB <- (hetot/total)
  
  out[[i]] <- data.frame(AA_cnt = AA_cnt, BB_cnt = BB_cnt, AB_cnt = hetot, Total_cnt = total, A_freq = A_freq, B_freq = B_freq, Heter_ratio = HE, Homo_ratio = HOMO, AA_freq = AA_freq, BB_freq = BB_freq, AB_freq = H_AB )
  
}

## Part 5: Store variable as your data ##

results <- do.call(rbind, out)

## Part 6: Write the data into a csv file ##
# This file will appear into your working directory

write.csv(results, file="results.csv", row.names = TRUE)
