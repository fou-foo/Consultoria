setwd('C:\\Users\\fou-f\\Desktop\\decidata')#es una buena practica ahislar el enviroment de R desde el incicio asi se tiene control sobre donde esta el working directory
#por cuestiones de compatibilidad entre sistemas opte por copiar los datos de la hoja 'health' y guardarlos como .csv
#si todo jala bien las variables los numericas 
datos <- read.csv(file='datos_decidata.csv', header=TRUE, stringsAsFactors = TRUE)
summary(datos)#un breve analisis descritico de las variables 
colnames(datos)[1] <- 'Id'
datos.originales <- datos 
#si bien varias varibles son del tipo ordinales y categoricas las codifico como dummies, se pierde informacion
#pero estoy contra tiempo, esto lo puedo mejorar con analisis de correspondencia canonica 
library(reshape)
dummies <- model.matrix(Id~.-1,data=datos[, c(1,2,4,6:26 )]) #genero las variables dummies necesarias
#agrego la informacion de las variables continuas
dummies <- cbind(datos[, c(1, 3,5, 27:30)], dummies)
#notamos con una inspeccion visual que la cariable 'turno' en su parte FULL se codifica de tres maneras distintas haciendo referencia al mismo hecho 
#consolidamos la informacion 
library(dplyr)
dummies <- dummies%>% group_by(Id) %>%mutate(turnoFULL= sum(`turnoFULL `, `turnoFULL   ` ))
dummies$`turnoFULL ` <- dummies$`turnoFULL   ` <- NULL
#ahora todas nuestras variables son numericas excepto el Id
index.id <- colnames(dummies)=='Id'
  #en vista de que se crearon var¿iables dummies no es prudente usar PCA sobre la matriz de correlaciones 
  #al usar la matriz de covarianzas se refleja nuestro conocimiento apriori de la 'calidad' de las caracteriaticas a explicar 
PCA <- princomp(dummies[, -index.id], cor = FALSE)
(varianza.muestral <- PCA$sdev**2) #notemos que la primer componente explica 54%
cumsum(varianza.muestral)/sum(varianza.muestral) #la primer componente principal explica poco más del 54% de la varianza total, por lo que es un buen candidado a indice
#esto tambien se puede mejorar sí conocemos más acerca de la importancia de cada variable indivualmente
#entonces mi propuesta de indice es la rotación que obtiene de los datos y proyectada sobre la primer componente principal 
#un indice más completo consideraria más componentes o bien las componentes resultado 
#de un metodo basado en kernels como PCA, aunque se desempeñan mejor su interpretabilidad es 
#menos sencilla que la de PCA clasico. NOTEMOS QUE LA VARIABLE MAS IMPORTANTE REGLEJADA
# EN EL INDICADOR ES la presion sistolica
datos.rotados <- as.matrix(dummies[, -index.id])%*%PCA$loadings
datos$Indicador_salud <- datos.rotados[,'Comp.1']
#con solo mirar los pesos de la primer componente principal, nuestro indicador de salud
#vemos que la variable 'mas' importante en nuestro enfoque es 'pa_sis' la presión sistolica
#lo cual es ad hoc con la teoria medica, en contraste las demas variables son poco influyentes en comparacion
hist(datos$Indicador_salud) # la distribucion de nuestro indicador parece uniforme
#sin embargo el test de bondad de ajuste más noble rechaza la hipotesis anterior
ks.test(datos$Indicador_salud, "punif", min(datos$Indicador_salud), max(datos$Indicador_salud) )
#CON FINES DE INTERPRETABILIDAD PROPONGO QUE EL INDICE 'Indice_salud' se mueva en un rango de 0 a 100
summary(datos$Indicador_salud)
datos$Indicador_salud <- 100*(datos$Indicador_salud-min(datos$Indicador_salud))/(range(datos$Indicador_salud)[2]-range(datos$Indicador_salud)[1])
hist(datos$Indicador_salud)#lo cual no altera drasticamente su distribución y lo hace mas facil de transmitir 
library(ggplot2)
(m <- median(datos$Indicador_salud))
ggplot(datos, aes(x=Indicador_salud, fill=I('lightblue')))+geom_density()+
  ggtitle('Distribución del Indicador de salud')+
  xlab('valordel indicadoren %')+ylab('frecuencia dentro de la empresa')+theme_minimal()+
geom_rug(aes(x =m, color=I('red'), size=3))  
library(cluster)
index.id <-  which(colnames(dummies)=='Id')
set.seed(0) #fijamos semilla con fines de reproducibilidad
#aprovechando la naturaleza mixta de los tipos de datos que tenemos podemos usar metodos 
#de clustering, que aunque no usan kernels, utilizan la información de las variables nominales
#utilizamos el clasico metodo de cluster basado en la metrica de 'gower'
disimilaridad <- daisy(x = datos.originales[, -index.id],metric = "gower")#notese que no utilizo el indicador que cree en la sección pasada pues ello implicaría problemas de multicolinealidad
#optamos por un metodo de cluster approximado como k.means, para tener información acerca
#del numero prudente de segmentos a diseñar
#dejamos para despues los metodos basados en clusterizacion 
#jerarquica que son mas demandantes en recursos y en complejidad computacional
#utilizare una versión de kmeans por practicidad-todo esto se puede mejorar
#con al metrica inducida por la norma L_2
#exploremos un numero prudente de segmentos maximizando la varianza entre grupos diferentes
var.fuera <- lapply(1:20, function(x){
  kmeans(dummies[, -index.id], centers = x, nstart = 100)$betweenss
})
plot(unlist(var.fuera)) #podemos conjeturar que 8 y 9 son numeros adecuados de segmentos
#fijamos el numero de segmentos a 8, de inicio y realicemos un clustering jerarquico 
#basado en la disimilaridad que ya calculamos
segmentacion <- agnes(x = disimilaridad, method = 'ward')
datos$segmento <- cutree(as.hclust(segmentacion), k = 9)
library(dendextend) #pintamos el arbol con los clusters
colores <- c('red', 'green', 'green4', 'pink4', 
             'yellow', 'purple', 'navy',
             'orange', 'magenta')
arbol.feliz <- segmentacion %>% as.dendrogram %>%
  set("branches_k_color", k=9, colores) %>% set("branches_lwd", 1.2) %>%
  set("labels_colors", rep(colores, each=50)) %>%
  set("labels_cex", c(.01)) 
plot(arbol.feliz)
#Asi que 9 parece ser un numero adecuado de segmentos 
datos <- datos[order(datos$Indicador_salud, decreasing = TRUE),]#ordenamos los datos en funcion del indice que construimos
plot(density(datos$Indicador_salud)) #en vista de la distribucion de nuestro indicador
#es simetrica y ademas no es unimodal propongo como evaluacion global del estado de la empresa
# a la mediana del indicador (el valor que acomula el 50% de probabilidad)
# pude usar la media pero este estimador es en general menos robusto 
# si se presentan casos futuros con indices muy altos o muy bajos la mediana esta dada por:
(m <- median(datos$Indicador_salud))
rug(m, col = 'red', ticksize = 0.1) #en color rojo 
#con fines ilustratimos muestro la media tambien
m1 <- mean(datos$Indicador_salud)
rug(m1, col = 'green', ticksize = 0.1) #que en este caso son muy parecidas por la simetria de la distribucion
#del indicador 
#Nota, si bien la segmentación relizada con más información de las variables originales,
#es más informativa que el indicador que se construyo a partir de variables dummie
#la segmentación no refleja a primera vista lo que dice el indicador que se correlaciona negativamente
#con la variable 'pa_sis' es decir que mientras mejor sea el puntaje del indicador se tiene 'mejor salud'
#asi la mediana refleja la buena salud de la empresa, mietras que los grupos 
#que se formaron son grupos de empleados con caracteristicas parecidas pero diferentes a 
#la buena salud. Es posible construir el indicador de manera que se adapte a los grupos 
#por ejemplo usando Kernel PCa en lugar de PCA, sin embargo tendriamos información repetida
#o que refleja lo mismo.
#realizo un test para validar la correlacion negativa del indicaor yla variable 'pa_sis'
cor.test(datos$Indicador_salud, datos$pa_sis)
#finalmente guardo los resultados en formato.csv para mayor portabilidad
write.csv(datos, file='resultado_de_analisis_salud.csv', row.names=FALSE )
