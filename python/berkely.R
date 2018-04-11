library(jsonlite)
salidaResultados <- 'C:\\Users\\fou-f\\Desktop\\Consultoria\\python'
setwd(salidaResultados)# nos cambiamos de directorio
out <- 'C:\\Users\\fou-f\\Desktop\\Consultoria\\python\\datosdemanda'
entrada <- 'C:\\Users\\fou-f\\Desktop\\Consultoria\\python\\datosdemanda.zip'
archivosLista <- dir(out) #los archivos con los cuales trabajaremos
archivosLista <- sort(archivosLista)
datos <- lapply(1:length(archivosLista), FUN = function(i){
  pathArchivo <- paste0(out,'\\', archivosLista[i])
  archivoNojason  = fromJSON(pathArchivo)
  archivoNojason
})
data <- do.call("rbind", datos)
data$valorEnlace <- NULL
data <- apply(data, 2, function(x){
  as.numeric(x)
})
data <- as.data.frame(data)
summary(data)
class(data)
dias <- sapply(strsplit(archivosLista,"[.]"), `[`, 1)
dias <- as.numeric(dias)
dias <- rep(dias, each = 24)
meses <- sapply(strsplit(archivosLista,"[.]"), `[`, 2)
meses <- as.numeric(meses)
meses <- rep(meses, each = 24)
data$dia <- dias
data$mes <- meses
data$anio <- 2016
semana_mayo <-  c('Su', 'Mo', 'Tu', 'We', 'Th', 'Fr', 'Sa')
semana_abril <- c('Fr', 'Sa', 'Su', 'Mo', 'Tu', 'We', 'Th')
data$diaSemana <- 'dia de la semana'
for ( i in 1:dim(data)[1])
{
  if(data[i, 'mes'] == 5)
  {
    data$diaSemana[i] <- semana_mayo[data$dia[i]%%7+1]
  } else
  {
    data$diaSemana[i] <- semana_abril[data$dia[i]%%7+1]
  } 
}
table(data$diaSemana)
#las cuentas estan en python 