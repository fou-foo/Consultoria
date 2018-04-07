setwd(dir = 'C:\\Users\\fou-f\\Desktop\\Consultoria\\ROLOCTS\\100\\')
#primero dejamos fijo el nombre del archivo a examinar
#archivo <- 'DET.TXT'
archivo <- 'E543.TXT'
#archivo <-'E87113.TXT'
#archivo <-'EXAMPLE1.TXT'
#archivo <- 'GEN.TXT'
archivo1 <- paste0(archivo, '_lambda.txt')
datos <- read.table(archivo1, header=TRUE, sep ='\t') #lectura de todas las soluciones encontradas
library(ggplot2)
p1<- ggplot(datos, aes(x=expected_value, y =risk, alpha=.0001 )) +
  geom_point()+theme_minimal()+
  theme(legend.position="none")+
  ggtitle(archivo)
p1 #graficamos los puntos de todas las soluciones encontradas
######################## lectura y graficacion de todos los puntos que son candidatos a formar la frontera de pareto
archivo2 <- paste0(archivo, 'dominantes.txt')
#archivo2 <- 'C:\\Users\\fou-f\\Downloads\\GEN.TXTdominantes.txt'
datos.verificacion <- read.table(archivo2, header=TRUE, sep ='\t')
p1+geom_point(data = datos.verificacion, aes(x=expected_value, y =risk,
                  colour=I('green') ))  
##########################
datos.verificacion <- datos.verificacion[order(datos.verificacion$expected_value,
                                               datos.verificacion$risk),]#ordenammos las posibles soluciones que forman la frontera por expected_value y para romper empates usamos el risk
# fijamos el minimo de pareto y el punto maximo en la expected_value que es facil de reconocer
mini <- datos.verificacion[1,]
maxi <- subset(datos.verificacion, risk==0)
maxi <- maxi[1,]
filtro <- data.frame() #estructura para guardar los puntos de la frontera
filtro <- rbind(filtro, mini)
actual <- mini
for(i in 1:dim(datos.verificacion)[1])
{
  #nos movemos en la dirección de la derivada y checamos puntualmente
  candidato <- datos.verificacion[i, ]
  if(candidato[,'expected_value'] >= actual[,'expected_value'] && (candidato[, 'expected_value']<=maxi[, 'expected_value']))
  {
    if( (candidato[, 'risk'] <= actual[,'risk']) )
    {
      filtro <- rbind(filtro, actual) # si encontramos un punto en direccion de la derivada con mayor funcion objetivo está en la frontera de pareto
      actual <- candidato
    }
  }
}
filtro <- unique(filtro)
#graficamos la frontera 
p1 + geom_point(data = filtro,
                 aes(x=expected_value, y =risk,
                  colour=I('purple') ) )
write.csv(filtro, paste0(archivo, 'fronteraPareto.csv'))
