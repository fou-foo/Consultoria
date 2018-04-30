library(readxl)     #paquete para leer archivos de excel
library(xlsx) # paquete de excel que puede crear problemas porque depende de java
library(dplyr)
raiz.de.arbol <- 1+.02 # dato sujeto a cambio: parametro 
CONVERSION.CARBONO.A.CO2 <- 3.67
cortes1 <- seq(5, 90, by =5)
cortes2 <- seq(7.4, 92.4, by = 5)
cortes.frecuencias <- c()
for(i in 1:length(cortes1))
{
  cortes.frecuencias <- c(cortes.frecuencias,cortes1[i], cortes2[i]  )
}
setwd('C:/Users/fou-f/Desktop/Consultoria/Ambiental2/29_de_4_2019/')
######################################### Analisis de cortes en funcion diametro del arbol
base.datos <- read_xlsx(path = 'Base_de_datos_de_trabajo_ANALISANDO_final_copya.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                        sheet = '2015_BD',  col_names = TRUE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
base.datos <- base.datos %>% select(`39.ESP`, `41.DN(cm)`, Vol., DENSIDAD, FEB, `%CC`)

base.datos2 <- read_xlsx(path = 'C:/Users/fou-f/Desktop/Consultoria/Ambiental2/29_de_4_2019/Base_de_datos_de_trabajo_ANALISANDO_final_copya.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                         sheet = '2015 CONSERVADA',  col_names = TRUE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
base.datos2 <- base.datos2 %>% select(`39.ESP`, `41.DN ( cm)`, Vol, DENSIDAD, FEB, `%CC`)

base.datos4 <- read_xlsx(path = 'C:/Users/fou-f/Desktop/Consultoria/Ambiental2/29_de_4_2019/Base_de_datos_de_trabajo_ANALISANDO_final_copya.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                         sheet = '2016BD',  col_names = TRUE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
base.datos4 <- base.datos4 %>% select(ESPECIE, `41.DN(cm)`, VOL., DENSIDAD, FEB, `%CC`)

base.datos5 <- read_xlsx(path = 'C:/Users/fou-f/Desktop/Consultoria/Ambiental2/29_de_4_2019/Base_de_datos_de_trabajo_ANALISANDO_final_copya.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                         sheet = '2017BD',  col_names = TRUE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
base.datos5 <- base.datos %>% select(`39.ESP`, `41.DN(cm)`, Vol., DENSIDAD, FEB, `%CC`)
names(base.datos) <- names(base.datos2) <- names(base.datos5) <- names(base.datos4)
base.datos$ESPECIE <- gsub("juniperus deppeana", "Juniperus deppeana", base.datos$ESPECIE )#por si las dudas


base.datos$`%CC` <- 0.531
base.datos$`%CC`[ base.datos$ESPECIE %in% c('Pinus arizonica', 'Pinus durangensis',
                   'Pinus leiophylla', 'Juniperus deppeana',
                   'Pinus ayacahuite', 'Pinus engelmannii', 'Pseudotsuga menziesii')  ] <- 0.521  # porque son pinos y para ellos el porcentaje de corte fue diferente
base.datos <- base.datos %>% mutate( BIOMASA = (VOL. * DENSIDAD * FEB   ), 
                                        CC = (BIOMASA* `%CC`), 
                                        CO2 = VOL.*DENSIDAD*FEB*raiz.de.arbol*`%CC`*CONVERSION.CARBONO.A.CO2 ) #RECALCULAMOS LOS CALCULOS DE jaime, para verificarlos
cuadro.especie <- base.datos %>% group_by(ESPECIE) %>% summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                  BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                  CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                  CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                  VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                  BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                  CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                  CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
library(ggplot2)
cuadro.especie <- head(cuadro.especie, 17)
#generamos graficos
ggplot(cuadro.especie, aes(x = ESPECIE, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
      coord_flip() +theme_minimal() +ggtitle("Volumen área de corta 2015")+ xlab('')
ggplot(cuadro.especie, aes(x = ESPECIE, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Biomasa área de corta 2015")+ xlab('')
ggplot(cuadro.especie, aes(x = ESPECIE, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de carbonó área de corta 2015")+ xlab('')
ggplot(cuadro.especie, aes(x = ESPECIE, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de CO2 área de corta 2015")+ xlab('')
#################La otra hoja 
base.datos2$`%CC`[ base.datos2$ESPECIE %in% c('Pinus arizonica', 'Pinus durangensis',
                                                   'Pinus leiophylla', 'Juniperus deppeana',
                                                   'Pinus ayacahuite', 'Pinus engelmannii') ] <- 0.521  # porque son pinos y para ellos el porcentaje de corte fue diferente
base.datos2 <- base.datos2 %>% mutate( BIOMASA = (VOL. * DENSIDAD * FEB   ), 
                                           CC = (BIOMASA* `%CC`), 
                                           CO2 = VOL.*DENSIDAD*FEB*raiz.de.arbol*`%CC`*CONVERSION.CARBONO.A.CO2 ) #RECALCULAMOS LOS CALCULOS DE jaime, para verificarlos
cuadro.especie2 <- base.datos2 %>% group_by(ESPECIE) %>% summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                     BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                     CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                     CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                     VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                     BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                     CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                     CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.especie2 <- head(cuadro.especie2, 12)
#generamos graficos
ggplot(cuadro.especie2, aes(x = ESPECIE, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Volumen área de corta 2015 (control)")+ xlab('')
ggplot(cuadro.especie2, aes(x = ESPECIE, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Biomasa área de corta 2015 (control)")+ xlab('')
ggplot(cuadro.especie2, aes(x = ESPECIE, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de carbonó área de corta 2015 (control)")+ xlab('')
ggplot(cuadro.especie2, aes(x = ESPECIE, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de CO2 área de corta 2015 (control)")+ xlab('')
############################### EL MISMO ANALISIS PERO PARA EL ANO 2016
base.datos4$`%CC` <- 0.531
base.datos4$`%CC`[ base.datos4$ESPECIE %in% c('Pinus arizonica', 'Pinus durangensis',
                                                   'Pinus leiophylla', 'Juniperus deppeana',
                                                   'Pinus ayacahuite', 'Pinus engelmannii', 'Pseudotsuga menziesii')  ] <- 0.521  # porque son pinos y para ellos el porcentaje de corte fue diferente
base.datos4 <- base.datos4 %>% mutate( BIOMASA = (VOL. * DENSIDAD * FEB   ), 
                                           CC = (BIOMASA* `%CC`), 
                                           CO2 = VOL.*DENSIDAD*FEB*raiz.de.arbol*`%CC`*CONVERSION.CARBONO.A.CO2 ) #RECALCULAMOS LOS CALCULOS DE jaime, para verificarlos
cuadro.especie4 <- base.datos4 %>% group_by(ESPECIE) %>% summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                     BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                     CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                     CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                     VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                     BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                     CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                     CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.especie4 <- head(cuadro.especie4, 16)
#generamos graficos
ggplot(cuadro.especie4, aes(x = factor(ESPECIE), y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Volumen área de corta 2016")+ xlab('')
ggplot(cuadro.especie4, aes(x = factor(ESPECIE), y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Biomasa área de corta 2016")+ xlab('')
ggplot(cuadro.especie4, aes(x = factor(ESPECIE), y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de carbonó área de corta 2016")+ xlab('')
ggplot(cuadro.especie4, aes(x = factor(ESPECIE), y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de CO2 área de corta 2016")+ xlab('')
#################La otra hoja 
base.datos5$ESPECIE <- gsub("juniperus deppeana", "Juniperus deppeana", base.datos5$ESPECIE )#por si las dudas
base.datos5$`%CC` <- 0.531
base.datos5$`%CC`[ base.datos5$ESPECIE %in% c('Pinus arizonica', 'Pinus durangensis',
                                                   'Pinus leiophylla', 'Juniperus deppeana',
                                                   'Pinus ayacahuite', 'Pinus engelmannii') ] <- 0.521  # porque son pinos y para ellos el porcentaje de corte fue diferente
base.datos5 <- base.datos5 %>% mutate( BIOMASA = (VOL. * DENSIDAD * FEB   ), 
                                           CC = (BIOMASA* `%CC`), 
                                           CO2 = VOL.*DENSIDAD*FEB*raiz.de.arbol*`%CC`*CONVERSION.CARBONO.A.CO2 ) #RECALCULAMOS LOS CALCULOS DE jaime, para verificarlos
cuadro.especie5 <- base.datos5 %>% group_by(ESPECIE) %>% summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                      BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                      CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                      CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                      VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                      BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                      CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                      CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.especie5 <- head(cuadro.especie5, 14)
cuadro.especie5$Especie <- factor(cuadro.especie5$ESPECIE)
#generamos graficos
ggplot(cuadro.especie5, aes(x = Especie, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Volumen área de corta 2017")+ xlab('')
ggplot(cuadro.especie5, aes(x = Especie, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Biomasa área de corta 2017")+ xlab('')
ggplot(cuadro.especie5, aes(x = Especie, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de carbonó área de corta 2017")+ xlab('')
ggplot(cuadro.especie5, aes(x = Especie, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  coord_flip() +theme_minimal() +ggtitle("Contenido de CO2 área de corta 2017")+ xlab('')
######################################################################
cuadro.especie5$Especie <- NULL
cuadro.especie$anio <- '2015'
cuadro.especie2$anio <- '2015 control'
cuadro.especie4$anio <- '2016'
cuadro.especie5$anio <- '2017'
base.datos$anio <- '2015'
base.datos2$anio <- '2015 control'
base.datos4$anio <- '2016'
base.datos5$anio <- '2017'
names(cuadro.especie) <- names(cuadro.especie2) <- names(cuadro.especie5) <-   names(cuadro.especie4)
totales.ano <- rbind(cuadro.especie, cuadro.especie2, cuadro.especie4, cuadro.especie5)
names(base.datos5) <- names(base.datos) <- names(base.datos2) <- names(base.datos4)
diametros <- rbind(base.datos, base.datos2, base.datos4, base.datos5)
diametros$topes <- cut(diametros$`41.DN(cm)`, breaks = cortes.frecuencias, include.lowest = TRUE )
cuadro.frecuencias <- diametros %>% filter(anio == '2015') %>% group_by(topes)  %>%summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                                             BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                                             CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                                             CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                                             VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                                             BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                                             CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                                             CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.frecuencias <- na.omit(cuadro.frecuencias)
ggplot(cuadro.frecuencias, aes(x = topes, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Volumen área de corta 2015")+ xlab('')
ggplot(cuadro.frecuencias, aes(x = topes, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Biomasa área de corta 2015")+ xlab('')
ggplot(cuadro.frecuencias, aes(x = topes, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
 theme_minimal() +ggtitle("Contenido de carbonó área de corta 2015")+ xlab('')
ggplot(cuadro.frecuencias, aes(x = topes, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
 theme_minimal() +ggtitle("Contenido de CO2 área de corta 2015")+ xlab('')
## 2015 area conservada 
cuadro.frecuencias2 <- diametros %>% filter(anio == '2015 control') %>% group_by(topes)  %>%summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                                             BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                                             CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                                             CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                                             VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                                             BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                                             CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                                             CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.frecuencias2 <- na.omit(cuadro.frecuencias2)
ggplot(cuadro.frecuencias2, aes(x = topes, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Volumen área de corta 2015 control")+ xlab('')
ggplot(cuadro.frecuencias2, aes(x = topes, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Biomasa área de corta 2015")+ xlab('')
ggplot(cuadro.frecuencias2, aes(x = topes, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de carbonó área de corta 2015")+ xlab('')
ggplot(cuadro.frecuencias2, aes(x = topes, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de CO2 área de corta 2015")+ xlab('')
##
cuadro.frecuencias3 <- diametros %>% filter(anio == '2016') %>% group_by(topes)  %>%summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                                                      BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                                                      CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                                                      CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                                                      VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                                                      BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                                                      CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                                                      CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.frecuencias3 <- na.omit(cuadro.frecuencias3)
ggplot(cuadro.frecuencias3, aes(x = topes, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Volumen área de corta 2016")+ xlab('')
ggplot(cuadro.frecuencias3, aes(x = topes, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Biomasa área de corta 2016")+ xlab('')
ggplot(cuadro.frecuencias3, aes(x = topes, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de carbonó área de corta 2016")+ xlab('')
ggplot(cuadro.frecuencias3, aes(x = topes, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de CO2 área de corta 2016")+ xlab('')

###
cuadro.frecuencias4 <- diametros %>% filter(anio == '2017') %>% group_by(topes)  %>%summarise(VOLUMEN.ESPECIE.suma = sum(VOL., na.rm = TRUE),
                                                                                              BIOMASA.ESPECIE.suma = sum(BIOMASA, na.rm = TRUE),
                                                                                              CC.ESPECIE.suma = sum(CC, na.rm =TRUE),
                                                                                              CO2.ESPECIE.suma = sum(CO2, na.rm = TRUE),
                                                                                              VOLUMEN.ESPECIE.promedio = mean(VOL., na.rm = TRUE),
                                                                                              BIOMASA.ESPECIE.promedio = mean(BIOMASA, na.rm = TRUE),
                                                                                              CC.ESPECIE.promedio = mean(CC, na.rm =TRUE),
                                                                                              CO2.ESPECIE.promedio = mean(CO2, na.rm = TRUE))
cuadro.frecuencias4 <- na.omit(cuadro.frecuencias4)
ggplot(cuadro.frecuencias4, aes(x = topes, y = VOLUMEN.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Volumen área de corta 2017")+ xlab('')
ggplot(cuadro.frecuencias4, aes(x = topes, y = BIOMASA.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Biomasa área de corta 2017")+ xlab('')
ggplot(cuadro.frecuencias4, aes(x = topes, y = CC.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de carbonó área de corta 2017")+ xlab('')
ggplot(cuadro.frecuencias4, aes(x = topes, y = CO2.ESPECIE.suma  )) + geom_bar(stat="identity") + 
  theme_minimal() +ggtitle("Contenido de CO2 área de corta 2017")+ xlab('')
write.xlsx(cuadro.especie, file ='cuadros.xlsx', 
           sheetName = '2015')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.especie2, file ='cuadros.xlsx', append=TRUE,
           sheetName = '2015 control')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.especie4, file ='cuadros.xlsx', append=TRUE,
           sheetName = '2016')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.especie5, file ='cuadros.xlsx', append=TRUE,
           sheetName = '2017')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.frecuencias, file ='cuadros.xlsx', append=TRUE,
           sheetName = '2015_f')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.frecuencias2, file ='cuadros.xlsx', append = TRUE,
           sheetName = '2015_control_frec')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.frecuencias3, file ='cuadros.xlsx',append=TRUE,
           sheetName = '2016_frecuencias')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro.frecuencias4, file ='cuadros.xlsx', append=TRUE, 
           sheetName = '2017_frecuencias')# tendras que cambiar la ruta de escritura

