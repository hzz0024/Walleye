# Walleye R script for hybridization analysis

The purpose of this R script is to score three genotype ratios (homozygous AA and BB, heterozygous AB) for each examined fish.. If you find this code useful in your research, please consider citing (manuscript under revision).

To download the R script and example files, please do
```
git clone git@github.com:hzz0024/Walleye.git
```

## Dependencies
- R 3.6
- RStudio 1.2

## How to use

- copy your genotypes into a csv file containing the sample IDs and the target SNPs (see example [here](https://github.com/hzz0024/Walleye/blob/master/Example_data/example_genotypes.csv))

- create a reference genotype file containing the target SNPs (have the same order as genotype file, see example [here](https://github.com/hzz0024/Walleye/blob/master/Example_data/example_ref.csv)) 

- Follow the instructions below to run the script in R or RStudio and the output will be a file with the scored markers and percentages

1. Set working directory
```
setwd("Your working directory")
```
2. Read into the reference and your genotype files
```
# Make sure to have the "Reference Genotype" file on your working directory
refsnp <- read.csv(file="example_ref.csv")
genotype <- read.csv(file="example_genotypes.csv", na.string="")
```
3. Change the length of your samples
```
# This number should reflect the number of samples you will be running
length_samples <- 60
out <- vector(mode="list", length=length_samples)
```
4. For loop to calculate AA/BB/AB
```
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
```
5. Store variable as your data
```
results <- do.call(rbind, out)
```
6. Write the data into a csv file
```
# This file will appear into your working directory
write.csv(results, file="results.csv", row.names = TRUE)
```
Questions or comments please send to: 
hzz0024@auburn.edu or hz269@cornell.edu
