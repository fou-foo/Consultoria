library(readxl)     #paquete para leer archivos de excel
library(xlsx) # paquete de excel que puede crear problemas porque depende de java

base.datos <- read.xlsx(file = 'C:/Users/fou-f/Desktop/Ambiental/INFORMACION-COLEGA-SAENZ.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                        sheetName = 'GENERO PINO', startRow = 2, header = TRUE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
datos.importantes <- data.frame(medidas=base.datos$X.DN)
hist(datos.importantes$medidas) #este es el dafault de R
cortes <- base.datos$NA..1
cortes <- cortes[!is.na(cortes)]
cortes <- cortes[-1]
cortes <-unique(as.numeric(as.character(cortes))) 
datos.importantes$corte <- cut(datos.importantes$medidas, breaks = cortes )
frecuencias <- table(datos.importantes$corte)
frecuencias <- as.data.frame(frecuencias)
frecuencias$relativas <- frecuencias$Freq/sum(frecuencias$Freq)
library(ggplot2)
ggplot(data = frecuencias, aes(x = factor(Var1), y = Freq, fill = factor(Var1))) + #indicas la variable de agrupamiento -debe ser tipo factor- y la de interes, tambien requieres llenar el parametro 'fill' con la variable con la que agrupa 
  geom_bar(stat="identity") +   # se construyen las barras
  xlab('medida') + theme_minimal() +  ggtitle('GENERO PINO')+ ylab('Frecuencias totales') +  guides(fill=FALSE) #quitar leyenda
write.xlsx(frecuencias, file ='C:/Users/fou-f/Desktop/Ambiental/frecuencias.xlsx',
           sheetName = 'frecuencias')# tendras que cambiar la ruta de escritura
