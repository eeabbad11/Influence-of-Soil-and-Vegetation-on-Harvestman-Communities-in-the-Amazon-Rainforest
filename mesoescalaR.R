#mesoescala!

envHar<-read.csv("mesoEscala.csv")

envHar$Parcela


envPPbio<-read.csv("originalData/PPBio2014_Environment.RFAD.csv")

envPPbio$Parcela<-paste(envPPbio$Trail,envPPbio$Plot)

envHar$Parcela

envHar$Trail<-gsub("DU_(.)O([0-9])_(.*)","\\1\\2",envHar$Parcela)
envHar$Plot<-as.numeric(gsub("DU_(.)O([0-9])_(.*)","\\3",envHar$Parcela))
envHar$Parcela2<-paste(envHar$Trail,envHar$Plot)

envHar$Parcela2

envPPbioReordered<-envPPbio[match(envHar$Parcela2,envPPbio$Parcela),]


dim(envHar)


spHar<-envHar[,grepl("_",colnames(envHar))]

dim(spHar)
dim(envPPbioReordered)


## Filter only RFAD

head(envPPbioReordered)


spHarDucke<-spHar[!is.na(envPPbioReordered$Location),]
envPPbioReorderedDucke<-envPPbioReordered[!is.na(envPPbioReordered$Location),]

####

# number of species in the whole study area
sum(colSums(spHarDucke>0)>0)


N<-rowSums(spHarDucke)
S<-rowSums(spHarDucke>0)

library(ggplot2)

#areia 
plot(N~envPPbioReorderedDucke$Sand)
plot(S~envPPbioReorderedDucke$Sand) 
ggplot(envPPbioReorderedDucke, aes(x = Sand, y + N, color = factor(Sand))) + 
  geom_point() +
  scale_color_brewer(palette = "Set1") + #5B8E7D#8CB369#F4A259

#fósforo
plot(N~envPPbioReorderedDucke$P)
plot(S~envPPbioReorderedDucke$P)

#numero de arvores 
plot(N~envPPbioReorderedDucke$Vegetation.structure)
plot(S~envPPbioReorderedDucke$Vegetation.structure)

#número de palmeiras
plot(N~envPPbioReorderedDucke$Vegetaion.structure2)
plot(S~envPPbioReorderedDucke$Vegetaion.structure2)


### Multiple regression using abundance and richness vs. predictor variables

#summary(lm(N~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2))

#
summary(glm(N~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2,family="poisson"))

plot(N~envPPbioReorderedDucke$Vegetaion.structure2)
m0<-lm(N~envPPbioReorderedDucke$Vegetaion.structure2)
abline(m0)


#summary(lm(S~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2))

summary(glm(S~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2,family="poisson"))













