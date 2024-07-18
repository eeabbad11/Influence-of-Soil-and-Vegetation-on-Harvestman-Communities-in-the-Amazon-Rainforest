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

ggplot(envPPbioReorderedDucke, aes(x = Sand, y = N, color = factor(Sand))) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +  #8CB369
  labs(title = "Relationship between N and sand quantity",
       x = "Sand Quantity",
       y = "Abundance (N)")
#fósforo
# Scatter plot of S against phosphorus quantity with attractive colors
ggplot(envPPbioReorderedDucke, aes(x = P, y = S, color = factor(P))) +
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Relationship between S and phosphorus quantity",
       x = "Phosphorus Quantity",
       y = "Richness (S)")

#numero de arvores 
# Scatter plot of N and S against the number of trees
ggplot(envPPbioReorderedDucke, aes(x = Vegetation.structure, y = N)) +
  geom_point(color = "blue") +
  geom_point(aes(y = S), color = "red") +
  labs(title = "Abundance and Richness vs. Number of Trees",
       x = "Number of Trees",
       y = "Abundance/Richness") +
  theme_minimal()

#número de palemiras 
# Scatter plot of N and S against the number of palm trees
ggplot(envPPbioReorderedDucke, aes(x = Vegetaion.structure2, y = N)) +
  geom_point(color = "green") +
  geom_point(aes(y = S), color = "orange") +
  labs(title = "Abundance and Richness vs. Number of Palm Trees",
       x = "Number of Palm Trees",
       y = "Abundance/Richness") +
  theme_minimal()


### Multiple regression using abundance and richness vs. predictor variables

#summary(lm(N~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2))

#
summary(glm(N~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2,family="poisson"))

plot(N~envPPbioReorderedDucke$Vegetaion.structure2)
m0<-lm(N~envPPbioReorderedDucke$Vegetaion.structure2)
abline(m0)


#summary(lm(S~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2))

summary(glm(S~envPPbioReorderedDucke$Sand+envPPbioReorderedDucke$P+envPPbioReorderedDucke$Vegetation.structure+envPPbioReorderedDucke$Vegetaion.structure2,family="poisson"))

