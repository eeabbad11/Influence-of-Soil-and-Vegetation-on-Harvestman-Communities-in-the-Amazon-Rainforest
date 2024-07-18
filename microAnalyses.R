####Análises Artigo1####
microData<-read.csv("microData.csv")
####organização da tabela####
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

#ANÁLISES#####
N<-rowSums(microSpAb)#abundancia
S<-rowSums(microSpAb>0)#riqueza
#ANOVA abundância e microhabitat
boxplot(N~microEnv)
summary(aovN<-aov(N~microEnv)) #vai dar os valores da análise
TukeyHSD(aovN) #compara os boxplots entre si: valor da diferença
# Define as cores para cada nível de microEnv
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")

#site coolors pra copiar e colar
# Cria o gráfico boxplot
boxplot(N ~ microEnv, col = cores, 
        xlab = "MicroHabitat", ylab = "Abundância", 
        main = "Distribuição de Abundância por MicroHabitat",cex.main = 1.0,border="black")
box(bty = "n")
library(ggplot2)
# Calcula a abundância e a riqueza
N <- rowSums(microSpAb)  # Abundância
S <- rowSums(microSpAb > 0)  # Riqueza
head(microSpAb)
# Realiza a análise de variância
aovN <- aov(N ~ microEnv)
tukey <- TukeyHSD(aovN)
tukey
S==N #?

# Cria o gráfico utilizando ggplot2
data <- data.frame(N = N, microEnv = factor(microEnv))

ggplot(data, aes(x = microEnv, y = N, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("Abundância") +
  ggtitle("Distribuição de Abundância por MicroHabitat") +
  theme_minimal() +
  theme(legend.position = "none")

#ANOVA riqueza e microhabitat
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9","#c8d5b9")
lm(S ~ microEnv)
summary(lm(S ~ microEnv))
# Cria o gráfico boxplot
boxplot(S ~ microEnv, col = cores, 
        xlab = "MicroHabitat", ylab = "Riqueza", 
        main = "Distribuição de Riqueza por MicroHabitat",cex.main = 1.0)
box(bty = "n")
boxplot(S~microEnv)
summary(aovS<-aov(S~microEnv))
TukeyHSD(aovS)

tapply(S, microEnv,mean)

#####ggplot####
# Cria o gráfico boxplot
data <- data.frame(S = S, microEnv = factor(microEnv))

ggplot(data, aes(x = microEnv, y = S, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("Riqueza") +
  ggtitle("Distribuição de Riqueza por MicroHabitat") +
  theme_minimal() +
  theme(legend.position = "none")


#pacote pra rodar a PCOA
install.packages("vegan")
library(vegan)
jac<-vegdist(microSpAb>0,"jac")
pcoa<-prcomp(jac)

#PCOA 
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9")

# Cria o gráfico boxplot


boxplot(scores(pcoa)[,1]~microEnv,las=2,col=cores, xlab = "MicroHabitat", ylab = "PCOA", 
        main = "PCOA",cex.main = 1.0)
box(bty = "n")
summary(aovPcoa<-aov(scores(pcoa)[,1]~microEnv))
lm(aovPcoa)
summary(lm(aovPcoa))
TukeyHSD(aovPcoa)
#COMPARAÇÃO DOS DIFERENTES HABITATS SEM RIQUEZA(NÚMERO ESPÉCIE) E SEM ABUNDÂNCIA (NÚMERO INDIVÍDUO)
####ggplot####
install.packages("vegan")
library(vegan)

# Cálculo da dissimilaridade de Jaccard
jac <- vegdist(microSpAb > 0, "jaccard")

# Realiza a PCoA
pcoa <- cmdscale(jac, k = 2)

# Cria o gráfico de boxplot
data <- data.frame(PCOA = pcoa[, 1], microEnv = factor(microEnv))

cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")
library(ggplot2)
ggplot(data, aes(x = microEnv, y = PCOA, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("PCOA") +
  ggtitle("PCOA") +
  theme_minimal() +
  theme(legend.position = "none")

summary(pcoa)
summary(jac)
data
pcoa





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


#Nesse exemplo, personalizamos as cores das caixas e dos pontos usando os argumentos fill e color, ajustamos o ângulo do texto do eixo x para facilitar a leitura com axis.text.x, e centralizamos o título do gráfico com plot.title. Sinta-se à vontade para ajustar esses parâmetros para atender às suas preferências visuais.


# Supondo que você tenha um vetor chamado "especies" com os nomes das espécies coletadas
# e um vetor chamado "microhabitats" com os micro-hábitats que cada espécie ocupou




library(tidyr)
tabela_frequencia <- pivot_wider(bioData, names_from = microDataHab, values_from = especie, values_fn = length)
print(tabela_frequencia)
cj

####Análises Artigo1####
# Carregar os dados
microData <- read.csv("microData.csv")

#### Organização da tabela ####
head(microData)

# Selecionar colunas dos habitats para organização
microDataHab <- microData[, c(5, 6, 8, 9, 10, 11)]

# Criação da variável microDataHabPA
microDataHabPA <- ifelse(microDataHab == "A" | microDataHab == "" | microDataHab == "a", 0, 1)
microDataHabPA[(microDataHabPA[, 4] > 0) & (microDataHabPA[, 1] > 0), 1] <- 0

# Função para determinar o tipo de habitat
microData$Habitat <- unlist(apply(microDataHabPA, 1, function(x) {
  resu <- c(1:6)[x > 0]
  if (length(resu) == 0) {
    resu <- NA
  }
  resu
}))

# Criação da coluna HabitatFac como fator
microData$HabitatFac <- factor(colnames(microDataHab)[microData$Habitat])

microData$N <- 1

# Filtrar linhas com HabitatFac não nulo
microData <- microData[!is.na(microData$HabitatFac),]

# Criação da tabela de resumo microSpAb
microSpAb <- tapply(microData$N, list(paste0(microData$Parcela, "_", microData$HabitatFac), microData$Taxon), sum, default = 0)

# Criação da variável microEnv
microEnv <- tapply(microData$HabitatFac, paste0(microData$Parcela, "_", microData$HabitatFac), function(x) as.character(x[1]))

# ANÁLISES
N <- rowSums(microSpAb) # Abundância
S <- rowSums(microSpAb > 0) # Riqueza

# ANOVA abundância e microhabitat
boxplot(N ~ microEnv)
summary(aovN <- aov(N ~ microEnv))
TukeyHSD(aovN)

# Definir cores para cada nível de microEnv
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")

# Criação do gráfico utilizando ggplot2
data <- data.frame(N = N, microEnv = factor(microEnv))

library(ggplot2)
ggplot(data, aes(x = microEnv, y = N, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("Abundância") +
  ggtitle("Distribuição de Abundância por MicroHabitat") +
  theme_minimal() +
  theme(legend.position = "none")

# ANOVA riqueza e microhabitat
boxplot(S ~ microEnv, col = cores,
        xlab = "MicroHabitat", ylab = "Riqueza",
        main = "Distribuição de Riqueza por MicroHabitat", cex.main = 1.0)
summary(aovS <- aov(S ~ microEnv))
TukeyHSD(aovS)
tapply(S, microEnv, mean)

# PCOA
library(vegan)
jac <- vegdist(microSpAb > 0, "jaccard")
pcoa <- prcomp(jac)

# Criação do gráfico boxplot PCOA
boxplot(scores(pcoa)[, 1] ~ microEnv, las = 2, col = cores, xlab = "MicroHabitat", ylab = "PCOA",
        main = "PCOA", cex.main = 1.0)
summary(aovPcoa <- aov(scores(pcoa)[, 1] ~ microEnv))
TukeyHSD(aovPcoa)

# Cálculo da dissimilaridade de Jaccard e criação do gráfico de boxplot PCOA com ggplot
jac <- vegdist(microSpAb > 0, "jaccard")
pcoa <- cmdscale(jac, k = 2)

data <- data.frame(PCOA = pcoa[, 1], microEnv = factor(microEnv))

ggplot(data, aes(x = microEnv, y = PCOA, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("PCOA") +
  ggtitle("PCOA") +
  ####Análises Artigo1####
# ... (código anterior até o ponto onde o erro ocorreu)

# Criação do gráfico utilizando ggplot2
data <- data.frame(N = N, microEnv = factor(microEnv))

library(ggplot2)
ggplot(data, aes(x = microEnv, y = N, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("Abundância") +
  ggtitle("Distribuição de Abundância por MicroHabitat") +
  theme_minimal() +
  theme(legend.position = "none")


# ANOVA riqueza e microhabitat
boxplot(S ~ microEnv, col = cores,
        xlab = "MicroHabitat", ylab = "Riqueza",
        main = "Distribuição de Riqueza por MicroHabitat", cex.main = 1.0)
summary(aovS <- aov(S ~ microEnv))
TukeyHSD(aovS)
tapply(S, microEnv, mean)

# PCOA
library(vegan)
jac <- vegdist(microSpAb > 0, "jaccard")
pcoa <- prcomp(jac)

# Criação do gráfico boxplot PCOA
boxplot(scores(pcoa)[, 1] ~ microEnv, las = 2, col = cores, xlab = "MicroHabitat", ylab = "PCOA",
        main = "PCOA", cex.main = 1.0)
summary(aovPcoa <- aov(scores(pcoa)[, 1] ~ microEnv))
TukeyHSD(aovPcoa)

# Cálculo da dissimilaridade de Jaccard e criação do gráfico de boxplot PCOA com ggplot
jac <- vegdist(microSpAb > 0, "jaccard")
pcoa <- cmdscale(jac, k = 2)

data <- data.frame(PCOA = pcoa[, 1], microEnv = factor(microEnv))
new_levels <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")

# Gráfico de boxplot PCOA com ggplot
ggplot(data, aes(x = microEnv, y = PCOA, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("PCOA") +
  ggtitle("PCOA") +
  theme_minimal() +
  theme(legend.position = "none")
####Análises Artigo1####
# ... (previous code until this point)

# Abundance and microhabitat ANOVA
boxplot(N ~ microEnv)
# How total abundance is distributed among different habitat types.

# Define a more solid green color
solid_green <- "#1f7a1f"

# Define the new levels for the microEnv variable
new_levels <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Your ggplot code
library(ggplot2)
ggplot(data.frame(N = N, microEnv = factor(microEnv)), aes(x = microEnv, y = N)) +
  geom_boxplot(fill = solid_green, color = "#36454F", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = solid_green) +
  scale_x_discrete(labels = new_levels) +
  labs(x = "Habitat Type", y = "Abundance", title = "Abundance Distribution by Habitat Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))

#### Análises Artigo1 ####
# ... (previous code until this point)

# PCoA analysis
install.packages("vegan")
library(vegan)

# Calculate Jaccard dissimilarity
jac <- vegdist(microSpAb > 0, "jaccard")

# Perform PCoA
pcoa <- cmdscale(jac, k = 2)

# Create a data frame for PCoA plot
data_pcoa <- data.frame(PCOA1 = pcoa[, 1], microEnv = factor(microEnv))

#### Análises Artigo1 ####
# ... (previous code until this point)

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
####Análises Artigo1####
microData<-read.csv("microData.csv")
# lê o arquivo csv e armazena os dados no Dataframe "microData"

####organização da tabela####
head(microData)
#pegou somente as colunas dos habitats (que queremos arrumar)
#ver as primeiras linhase colunas- esboço tabela
#exibe as primeiras linhas do microData pra ver se foi carregado corretamente

microDataHab<-microData[,c(5,6,8,9,10,11)]
#só as colunas pra arrumar
#se o valor em "microDataHab" for igual a "A" (maiúsculo), ou se estiver em branco (representado por ""), ou se for igual a "a" (minúsculo), a nova variável "microDataHabPA" receberá o valor 0. Caso contrário, ou seja, se o valor em "microDataHab" for diferente dessas condições, a nova variável "microDataHabPA" receberá o valor 1.
microDataHabPA<-ifelse(microDataHab=="A"|microDataHab==""|microDataHab=="a",0,1)
#Quando as duas partes da condição são verdadeiras (ou seja, quando os valores da coluna 4 e da coluna 1 são maiores que zero), o código atribui o valor 0 à primeira coluna da variável "microDataHabPA" para aquelas observações que satisfazem essa condição.
microDataHabPA[(microDataHabPA[,4]>0)&(microDataHabPA[,1]>0),1]<-0
# função "apply" está sendo utilizada para aplicar uma função anônima a cada linha da variável "microDataHabPA". A função anônima recebe cada linha da variável "microDataHabPA" como um vetor "x", verifica quais elementos de "x" são maiores que zero e retorna um vetor "resu" com os índices desses elementos.

#A linha "resu <- c(1:6)[x > 0]" cria o vetor "resu" com os índices dos elementos de "x" que são maiores que zero. Esses índices representam o tipo de habitat em que cada observação foi encontrada.

#A linha "if(length(resu)==0){ resu<-NA }" verifica se o vetor "resu" está vazio (ou seja, se não foram encontrados habitats para essa observação). Se o vetor estiver vazio, a linha atribui o valor "NA" a "resu">0.

microData$Habitat<-unlist(apply(microDataHabPA,1,function(x){
  resu<-c(1:6)[x>0]
  if(length(resu)==0){
    resu<-NA
  }
  resu
}))
#esse trecho de código cria uma nova coluna chamada "Habitat" no DataFrame "microData", que contém informações sobre o tipo de habitat em que cada observação foi encontrada, com base nas condições definidas no DataFrame "microDataHabPA".

head(microData)

#Dessa forma, a nova coluna "HabitatFac" é uma versão categórica da coluna "Habitat", que indica o tipo de habitat em que cada observação foi encontrada, utilizando os nomes das colunas da variável "microDataHab" como níveis do fator. Verificar alterações.

microData$HabitatFac<-factor(colnames(microDataHab)[microData$Habitat])
#a nova coluna "HabitatFac" contém os nomes das colunas da variável "microDataHab" que correspondem aos habitats indicados pela coluna "Habitat" do DataFrame "microData". Isso cria uma representação categórica dos tipos de habitat em que cada observação foi encontrada.


microData$N<-1
#criou uma coluna N na microData e atribuiu o valor 1 pra ela

microData<-microData[!is.na(microData$HabitatFac),]
# "microData" terá apenas as linhas em que a coluna "HabitatFac" não é NA, ou seja, as linhas onde os valores de habitat foram definidos e não são ausentes. Isso pode ser útil para remover observações sem informações relevantes.
#a variável "microSpAb" contém uma tabela resumindo os dados da variável "microData", indicando a soma dos valores da coluna "N" para cada combinação de "Parcela_HabitatFac" e "Taxon".

microSpAb<-tapply(microData$N,list(paste0(microData$Parcela,"_",microData$HabitatFac),microData$Taxon),sum,default = 0)
# função tapply() para calcular a soma da coluna "N" (que contém o valor 1) do "microData" com base em combinações de "Parcela_HabitatFac" e "Taxon". microData$N: Isso acessa a coluna "N" do "microData", que contém o valor 1 para cada observação.list(...): Aqui está criando uma lista de fatores para usar como índices.paste0(microData$Parcela, "_", microData$HabitatFac): Esta parte cria uma combinação de "Parcela" e "HabitatFac" separados por "_". Isso cria um fator composto que representa uma combinação única de parcela e habitat.microData$Taxon: Isso é usado como o segundo nível de índice. Ele representa a classificação taxonômica.sum: A função sum é aplicada para somar os valores da coluna "N" que correspondem a cada combinação única de "Parcela_HabitatFac" e "Taxon".default = 0: Isso define o valor padrão como 0 para as combinações onde não há correspondência, ou seja, onde não há observações correspondentes.Portanto, o resultado "microSpAb" é uma tabela que resume os dados do DataFrame "microData", mostrando a soma dos valores da coluna "N" para cada combinação única de "Parcela_HabitatFac" e "Taxon". Isso pode ser útil para análises de abundância de espécies em diferentes habitats e parcelas.

microEnv<-tapply(microData$HabitatFac,paste0(microData$Parcela,"_",microData$HabitatFac),function(x)as.character(x[1]))
#taply() para criar um vetor que representa o tipo de habitat (HabitatFac) associado a cada combinação única de "Parcela" e "HabitatFac". O resultado do código é um vetor "microEnv" que contém os tipos de habitat associados a cada combinação única de "Parcela_HabitatFac". Esse vetor pode ser usado em análises subsequentes para entender como os diferentes habitats estão distribuídos entre as parcelas.


#ANÁLISES#####
library(ggplot2)
N<-rowSums(microSpAb)#abundancia
S<-rowSums(microSpAb>0)#riqueza

#ANOVA abundância e microhabitat
boxplot(N~microEnv)
# como a abundância total está distribuída entre os diferentes tipos de habitat.

# Defina os novos nomes dos níveis da variável microEnv
# Define a solid shade of green
solid_green <- "#1f7a1f"

# Defina os novos nomes dos níveis da variável microEnv
new_names <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Define a solid shade of green
solid_green <- "#1f7a1f"

# Define a solid shade of black
solid_black <- "#000000"

# Seu ggplot code
ggplot(data.frame(N = N, microEnv = factor(microEnv)), aes(x = microEnv, y = N)) +
  geom_boxplot(fill = solid_green, color = solid_black, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = solid_black) +
  scale_x_discrete(labels = new_names) +
  labs(x = "Habitat Type", y = "Abundance", title = "Abundance Distribution by Habitat Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))


#Nesse exemplo, personalizamos as cores das caixas e dos pontos usando os argumentos fill e color, ajustamos o ângulo do texto do eixo x para facilitar a leitura com axis.text.x, e centralizamos o título do gráfico com plot.title. Sinta-se à vontade para ajustar esses parâmetros para atender às suas preferências visuais.



#agora ANOVA com a riqueza: 
boxplot(S~microEnv)
#Para analisar a variação da riqueza de espécies em diferentes micro-habitats usando o valor F de Fisher e o valor p, você pode realizar uma análise de variância (ANOVA) sobre a riqueza (número de espécies) em cada micro-habitat. A ANOVA testará se há diferenças significativas na média da riqueza entre os diferentes grupos de micro-habitats. Os valores F e p fornecerão informações sobre a significância dessas diferenças.
# Defina os novos nomes dos níveis da variável microEnv
new_names <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Define a solid shade of green
solid_green <- "#1f7a1f"

# Define a solid shade of black
solid_black <- "#000000"

# Seu ggplot code para a riqueza (S)
ggplot(data.frame(S = S, microEnv = factor(microEnv)), aes(x = microEnv, y = S)) +
  geom_boxplot(fill = solid_green, color = solid_black, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = solid_black) +
  scale_x_discrete(labels = new_names) +
  labs(x = "Habitat Type", y = "Richness", title = "Richness Distribution by Habitat Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))

# Realizar a ANOVA para Riqueza por MicroHabitat
aovRiqueza <- aov(S ~ microEnv)

# Obter o valor F e o valor p
summary(aovRiqueza)
#summary(aovRiqueza) exibirá um resumo da análise de variância, incluindo o valor F e o valor p. O valor F mede a variação entre as médias dos grupos em relação à variação dentro dos grupos. O valor p indica a probabilidade de obter uma diferença tão extrema quanto a observada, supondo que as médias reais sejam iguais.

# Executar o teste de Tukey para comparações múltiplas
tukeyRiqueza <- TukeyHSD(aovRiqueza)

# Resumo do teste de Tukey
tukeyRiqueza





######no ggplot####
# Calcula a abundância e a riqueza
N <- rowSums(microSpAb)  # Abundância
S <- rowSums(microSpAb > 0)  # Riqueza
head(microSpAb)
# Realiza a análise de variância
aovN <- aov(N ~ microEnv)
tukey <- TukeyHSD(aovN)
tukey
S==N #?


#ANOVA riqueza e microhabitat
#precisa ainda?????????
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9")
lm(S ~ microEnv)
summary(lm(S ~ microEnv))
# Cria o gráfico boxplot
boxplot(S ~ microEnv, col = cores, 
        xlab = "MicroHabitat", ylab = "Riqueza", 
        main = "Distribuição de Riqueza por MicroHabitat",cex.main = 1.0)
box(bty = "n")
boxplot(S~microEnv)
summary(aovS<-aov(S~microEnv))
TukeyHSD(aovS)

tapply(S, microEnv,mean)

#####ggplot####
# Cria o gráfico boxplot
data <- data.frame(S = S, microEnv = factor(microEnv))

ggplot(data, aes(x = microEnv, y = S, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("Riqueza") +
  ggtitle("Distribuição de Riqueza por MicroHabitat") +
  theme_minimal() +
  theme(legend.position = "none")

# ANOVA e Tukey para Abundância por MicroHabitat
aovN <- aov(N ~ microEnv)
tukeyN <- TukeyHSD(aovN)
tukeyN

# ANOVA e Tukey para Riqueza por MicroHabitat
aovS <- aov(S ~ microEnv)
tukeyS <- TukeyHSD(aovS)
tukeyS


#pacote pra rodar a PCOA
install.packages("vegan")
library(vegan)
jac<-vegdist(microSpAb>0,"jac")
pcoa<-prcomp(jac)
pcoa
jac
summary(pcoa)
summary(jac)
#PCOA 
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9")

# Cria o gráfico boxplot


boxplot(scores(pcoa)[,1]~microEnv,las=2,col=cores, xlab = "MicroHabitat", ylab = "PCOA", 
        main = "PCOA",cex.main = 1.0)
box(bty = "n")
summary(aovPcoa<-aov(scores(pcoa)[,1]~microEnv))
lm(aovPcoa)
summary(lm(aovPcoa))
TukeyHSD(aovPcoa)



#COMPARAÇÃO DOS DIFERENTES HABITATS SEM RIQUEZA(NÚMERO ESPÉCIE) E SEM ABUNDÂNCIA (NÚMERO INDIVÍDUO)
####ggplot####
install.packages("vegan")
library(vegan)

# Cálculo da dissimilaridade de Jaccard
jac <- vegdist(microSpAb > 0, "jaccard")

# Realiza a PCoA
pcoa <- cmdscale(jac, k = 2)

# Cria o gráfico de boxplot
data <- data.frame(PCOA = pcoa[, 1], microEnv = factor(microEnv))

cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")
library(ggplot2)
ggplot(data, aes(x = microEnv, y = PCOA, fill = microEnv)) +
  geom_boxplot() +
  scale_fill_manual(values = cores) +
  xlab("MicroHabitat") +
  ylab("PCOA") +
  ggtitle("PCOA") +
  theme_minimal() +
  theme(legend.position = "none")

summary(pcoa)
summary(jac)
data
pcoa

library(ggplot2)

# Defina as cores para cada nível de microEnv
verde_solido <- "#1f7a1f"
cores <- c("#344e41", "#3a5a40", "#588157", "#a3b18a", "#c8d5b9", "#c8d5b9")

# Defina os novos nomes dos níveis da variável microEnv
novos_nomes <- c("Shrub", "Tree", "Herbaceous", "Leaf litter", "Palm", "Decaying Log")

# Função para criar um gráfico de boxplot com jitter
create_boxplot_with_jitter <- function(data, x_var, y_var, x_label, y_label, title) {
  ggplot(data, aes(x = get(x_var), y = get(y_var))) +
    geom_boxplot(fill = verde_solido, color = "#36454F", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, color = verde_solido) +
    scale_x_discrete(labels = novos_nomes) +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5))
}

# Gráfico de Abundância por MicroHabitat
abundancia_plot <- create_boxplot_with_jitter(
  data = microData,
  x_var = "microEnv",
  y_var = "N",
  x_label = "Habitat Type",
  y_label = "Abundance",
  title = "Abundance Distribution by Habitat Type"
)

# Gráfico de Riqueza por MicroHabitat
riqueza_plot <- create_boxplot_with_jitter(
  data = microData,
  x_var = "microEnv",
  y_var = "S",
  x_label = "Habitat Type",
  y_label = "Riqueza",
  title = "Riqueza Distribution by Habitat Type"
)

# Gráfico de PCOA por MicroHabitat
pcoa_plot <- create_boxplot_with_jitter(
  data = data.frame(PCOA = pcoa[, 1], microEnv = factor(microEnv)),
  x_var = "microEnv",
  y_var = "PCOA",
  x_label = "Habitat Type",
  y_label = "PCOA",
  title = "PCOA Distribution by Habitat Type"
)

# Visualizar os gráficos
gridExtra::grid.arrange(abundancia_plot, riqueza_plot, pcoa_plot, ncol = 2)
