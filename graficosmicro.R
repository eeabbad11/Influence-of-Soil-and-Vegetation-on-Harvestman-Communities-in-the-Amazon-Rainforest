####Análises Artigo1#### organizadoooooo
microData<-read.csv("microData.csv")
#### organização da tabela ####
head(microData)
microDataHab<-microData[,c(5,6,8,9,10,11)]
microDataHabPA<-ifelse(microDataHab=="A"|microDataHab==""|microDataHab=="a",0,1)
microDataHabPA[(microDataHabPA[,4]>0)&(microDataHabPA[,1]>0),1]<-0
microData$Habitat<-unlist(apply(microDataHabPA,1,function(x){
  resu<-c(1:6)[x>0]
  if(length(resu)==0){
    resu<-NA
  }
  resu
}))
head(microData)
microData$HabitatFac<-factor(colnames(microDataHab)[microData$Habitat])
microData$N<-1
#criou uma coluna N na microData e atribuiu o valor 1 pra ela
microData<-microData[!is.na(microData$HabitatFac),]
#a variável "microSpAb" contém uma tabela resumindo os dados da variável "microData", indicando a soma dos valores da coluna "N" para cada combinação de "Parcela_HabitatFac" e "Taxon".
microSpAb<-tapply(microData$N,list(paste0(microData$Parcela,"_",microData$HabitatFac),microData$Taxon),sum,default = 0)
#colocar no chatGPT
microEnv<-tapply(microData$HabitatFac,paste0(microData$Parcela,"_",microData$HabitatFac),function(x)as.character(x[1]))

####teste 
#ANÁLISES#####
library(ggplot2)
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")
N<-rowSums(microSpAb)#abundancia
S<-rowSums(microSpAb>0)#riqueza
#ANOVA abundância e microhabitat
boxplot(N~microEnv)
# como a abundância total está distribuída entre os diferentes tipos de habitat.

# Defina um tom de verde mais sólido
verde_solido <- "#1f7a1f"

# Defina os novos nomes dos níveis da variável microEnv
novos_nomes <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Seu código ggplot
ggplot(data.frame(N = N, microEnv = factor(microEnv)), aes(x = microEnv, y = N)) +
  geom_boxplot(fill = verde_solido, color = "#36454F", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = verde_solido) +
  scale_x_discrete(labels = novos_nomes) +
  labs(x = "Habitat Type", y = "Abundance", title = "Abundance Distribution by Habitat Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))

#Riqueza S
# Carregar a biblioteca ggplot2
library(ggplot2)

# Carregar os dados
microData <- read.csv("microData.csv")

# Organização da tabela
microDataHab <- microData[, c(5, 6, 8, 9, 10, 11)]
microDataHabPA <- ifelse(microDataHab == "A" | microDataHab == "" | microDataHab == "a", 0, 1)
microDataHabPA[(microDataHabPA[, 4] > 0) & (microDataHabPA[, 1] > 0), 1] <- 0
microData$Habitat <- unlist(apply(microDataHabPA, 1, function(x) {
  resu <- c(1:6)[x > 0]
  if (length(resu) == 0) {
    resu <- NA
  }
  resu
}))
microData$HabitatFac <- factor(colnames(microDataHab)[microData$Habitat])
microData$N <- 1
microData <- microData[!is.na(microData$HabitatFac),]
microSpAb <- tapply(microData$N, list(paste0(microData$Parcela, "_", microData$HabitatFac), microData$Taxon), sum, default = 0)
microEnv <- tapply(microData$HabitatFac, paste0(microData$Parcela, "_", microData$HabitatFac), function(x) as.character(x[1]))

# Cálculo de riqueza
S <- rowSums(microSpAb > 0)

# Definição das cores
verde_solido <- "#1f7a1f"
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")

# Definição dos nomes dos habitats
novos_nomes <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Criação do gráfico de boxplot
ggplot(data.frame(S = S, microEnv = factor(microEnv)), aes(x = microEnv, y = S)) +
  geom_boxplot(fill = verde_solido, color = "#36454F", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = verde_solido) +
  scale_x_discrete(labels = novos_nomes) +
  labs(x = "Habitat Type", y = "Species Richness", title = "Species Richness by Habitat Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))
#PCOA 
# Define the new levels and a subtle palette of green shades
new_levels <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")
subtle_green_palette <- c("#BFD8B8", "#A6C7A3", "#8EB793", "#76A983", "#5E9A73", "#467963")

# Create a boxplot for PCoA1 with new levels and the subtle green palette
ggplot(data_pcoa, aes(x = microEnv, y = PCOA1, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = subtle_green_palette) +
  scale_x_discrete(labels = new_levels) +  # Use new levels for x-axis labels
  xlab("MicroHabitat") +
  ylab("PCoA1") +
  ggtitle("PCoA Analysis") +
  theme_minimal() +
  theme(legend.position = "none")

#mapa
install.packages("ggmap")
library(ggplot2)
library(ggmap)
install.packages("readxl")

# Carregar a biblioteca readxl
library(readxl)
install.packages("leaflet")
library(leaflet)
install.packages("magrittr")
library(magrittr)

coords<-read.csv("bioData.csv")
library(rgdal)
lat<-coords$Latitude
long<-coords$Longitude
latlong<-cbind(long,lat)

library(leaflet)

coords <- read.csv("bioData.csv")

m <- leaflet(data = coords) %>%
  addTiles() %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude, radius = 5, color = "blue")

m

install.packages(c("ggplot2", "leaflet", "magrittr", "rgdal", "readxl"))
library(ggplot2)
library(leaflet)
library(magrittr)
library(rgdal)
library(readxl)

# Carregar os dados
coords <- read.csv("bioData.csv")

brasil_map <- leaflet() %>%
  addTiles() %>%
  setView(lng = -51.9253, lat = -14.2350, zoom = 4)  # Coordenadas centrais do Brasil e zoom
brasil_map <- brasil_map %>%
  addCircleMarkers(data = coords, lng = ~Longitude, lat = ~Latitude, radius = 5, color = "blue")
# Destacar o estado do Amazonas (usando a sigla "AM") com a cor #57cc99
amazonas_coords <- coords[coords$Estado == "AM", ]
brasil_map <- brasil_map %>%
  addCircleMarkers(data = amazonas_coords, lng = ~Longitude, lat = ~Latitude, radius = 5, color = "#57cc99")

brasil_map
  
# Instalar e carregar as bibliotecas necessárias
install.packages(c("ggplot2", "rgdal", "readxl"))
library(ggplot2)
library(rgdal)
library(readxl)

# Carregar os dados
coords <- read.csv("bioData.csv")
readOGR("shapes")
mapa<-readOGR("shapes")
mapa
