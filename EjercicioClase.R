library(Biostrings)

##### version chida
### dividiendo el genoma 
ecoli <- readDNAStringSet("GCF_000005845.2_ASM584v2_genomic(1).fna.gz")
ecoli
ecoli2 <- subseq(ecoli, start = width(ecoli) - 1000)
ecoli2

####generando reads
pedazos2 <- c()
for(i in 1:10){
  inicio <- sample(1:width(ecoli2), 1)
  final <- inicio + sample(4:5, 1)
  pedazos <- subseq(ecoli2, start= inicio, end = final)
  pedazos2 <- c(pedazos2, pedazos)
}
pedazos2

pedazos2 <- as.data.frame(pedazos2)
dim(pedazos2)

#### pasandolos a string
pedazos3 <- c()
cont2 <- 1
for(i in 1:dim(pedazos2)[2]){
  pedazos3 <- c(pedazos3, toString(pedazos2[cont2]))
  cont2 <- cont2 + 1
}
pedazos3
pedazos3[1]

##### conteo de patrones
cont3 <- 1
conteos <- c()
for(i in 1:dim(pedazos2)[2]){
conteos <- c(conteos,vcountPattern(pedazos3[cont3], ecoli2))
  cont3 <- cont3 + 1
  
}

conteos2 <- matrix(conteos, ncol =1)
conteos2 <- as.numeric(conteos2)
pedazitos <- matrix(pedazos3, ncol= 1)
pedazitos
datardos <- cbind(pedazitos, conteos2)
datardos
names <- c("reads", "apariciones")
colnames(datardos) <- names
datardos <- as.data.frame(datardos)
plot(datardos$apariciones)
print(paste(conteos, pedazos3))

barplot(as.numeric(datardos$apariciones), names = datardos$reads, col = "blue"
        , horiz = T, las=1)


