library(readxl)     #paquete para leer archivos de excel
library(dplyr)      #paquete para manipular datos estilo SQL
base.datos <- read_xlsx(path = 'C:/Users/fou-f/Desktop/ta Garcia/MEMORIA-DE-CALCULO.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                  sheet='BASE DE DATOS', skip = 1, col_names = FALSE) #se lee la 'Hoja base de datos' que planeamos en un futuro sea un archivo .csv
colnames(base.datos) <- c("RODAL", "SITIO", "No", "ESPECIE",
                      "NOMBRE.CIENTIFICO", "DIAMETRO", "ALTURA",
                      "EDAD", "T.P", "AB", "VOL.UNITARIO", 
                      "GPO.SILVICOLA")          #le cambio el nombre a las columnas para que sea mas facil la programación
base.superficie <- read_xlsx(path = 'C:/Users/fou-f/Desktop/ta Garcia/MEMORIA-DE-CALCULO.xlsx', #CAMBIAR RUTA DE ARCHIVO A LA TUYA
                        sheet='SUP RODAL', skip = 2, col_names = TRUE) #se lee la que contiene la superficie por RODAL
names(base.superficie) <-  c("RODAL", "Superficie")
base.datos <- merge(base.datos, base.superficie, by.x = 'RODAL', by.y = 'RODAL', all.x =TRUE) #cruzamos las tablas para pegarle a la base da datos original la superficie el RODAL
base.datos <- base.datos[-dim(base.datos)[1],] # quito en ultimo reglon que es un total
############ se construye el cuadro 3
cuadro3.1 <- base.datos %>% select(RODAL, SITIO,NOMBRE.CIENTIFICO, AB, 
                                   VOL.UNITARIO,GPO.SILVICOLA, Superficie) 
cuadro3.1$c5a15 <- cuadro3.1$c20a25 <- cuadro3.1$c30a40 <- NA
parser <-list("5a15" = 10, "20a25" = 9, "30a40"=8 ) #tener cuidado con el numeor de columnas asociado
    #esta parte convierte la columna de GPO.SILVIOLA en varias con su respectivo VOL.UNITARIO, si no usamos excel se podria hacer sin el 'for'
columnas <- unlist(parser[cuadro3.1$GPO.SILVICOLA]) #convertimos usando una lista 
for(i in 1:length(columnas))
{
  cuadro3.1[i, columnas[i]] <- cuadro3.1[i, 'VOL.UNITARIO'] #asi cada registro tiene el VOl.UNITARIO en la columna de GRP.SILVICOLA que le corresponde
}
      #primer suma acumulado por RODAL y especie del volumen
cuadro3 <- cuadro3.1 %>% group_by(RODAL, NOMBRE.CIENTIFICO ) %>%
            summarise(c20a25 = sum(c20a25, na.rm=TRUE),
                      c30a40 = sum(c30a40, na.rm=TRUE),
                      c5a15 = sum(c5a15, na.rm=TRUE), Total.general= sum(c20a25, c30a40, c5a15),
                      Suma.de.AB = sum(AB, na.rm=TRUE))
        #calculo del numero de sitios
cuadro3.2 <- cuadro3.1 %>% select(RODAL,SITIO) %>%group_by(RODAL) %>% summarise(No.Sitios=length(unique(SITIO)) /10 ) #en porcentaje
    #se actualizan suma acumulado por RODAL y especie del volumen entre No.sitios
cuadro3 <- merge(cuadro3, cuadro3.2, all.x =TRUE)
cuadro3 <- cuadro3 %>% mutate(c20a25 = c20a25/ No.Sitios,
                              c30a40 = c30a40/ No.Sitios,
                              c5a15 = c5a15/ No.Sitios)
cuadro3$Total.general <- cuadro3$c20a25 + cuadro3$c30a40 + cuadro3$c5a15
cuadro3 <- cuadro3 %>% rename(Especie = NOMBRE.CIENTIFICO, er.entre.ha =Total.general  ) # se cambio nombre de columna
cuadrox <- cuadro3 #hago una copia para futuros calculos
cuadro3$c20a25 <- cuadro3$c30a40 <- cuadro3$c5a15 <- NULL                  #se elimnan columnas no cesarias 
cuadro3 <- merge(cuadro3, base.superficie, all.x=TRUE)
ic <- .3#se planea que el IC se pueda cambiar en el futuro
cuadro3 <- cuadro3 %>% mutate(er.totales = er.entre.ha*Superficie,
                   ab.entre.ha=Suma.de.AB/No.Sitios,
                   IC = ic) #se planea que el IC se pueda cambiar en el futuro 
cuadro3 <- cuadro3%>% mutate( M3.VTA.por.ha.verde = er.entre.ha*(1-ic),
                   Area.basal.M2.por.ha = ab.entre.ha*(1-ic),
                   M3.VTA.por.ha = er.entre.ha*IC,
                   M3.VTA.por.UMM = M3.VTA.por.ha*Superficie
                   )
cuadro3$SITIO <- NULL
cuadro3 <- cuadro3 %>% select(RODAL, Superficie, Especie, er.entre.ha, 
                             er.totales, ab.entre.ha, IC,
                             M3.VTA.por.ha.verde, Area.basal.M2.por.ha,
                             M3.VTA.por.ha, M3.VTA.por.UMM) #se reordenan columnas
##############Se termina cuadro 3
############comienzan cuadros de la hoja 4 nom
cuadro4.1 <- cuadro3%>%select(Especie, er.entre.ha, er.totales, ab.entre.ha, M3.VTA.por.ha.verde, 
                              Area.basal.M2.por.ha, M3.VTA.por.ha, M3.VTA.por.UMM ) %>%
                      group_by(Especie) %>%
                        summarise(Suma.de.erentre.ha = sum(er.entre.ha) ,
                                  Suma.de.er.totales = sum(er.totales),
                                  Suma.de.ab.entre.ha = sum(ab.entre.ha),
                                  Suma.de.RES.M3.VTA.por.ha = sum(M3.VTA.por.ha.verde), 
                                  Suma.de.RES.Area.basal.M2.por.ha = sum(Area.basal.M2.por.ha),
                                  Suma.de.POS.M3.VTA.por.ha = sum(M3.VTA.por.ha),
                                  Suma.de.POS.M3.VTA.por.UMM = sum(M3.VTA.por.UMM))
cuadro4.2 <- cuadro3%>%select(Especie, er.totales, M3.VTA.por.UMM ) %>%
  group_by(Especie) %>%
  summarise(Suma.de.er.totales = sum(er.totales),
            Suma.de.POS.M3.VTA.por.UMM = sum(M3.VTA.por.UMM))
#se terminan cuadros hoja 4 nom
#######comienzan cuadro de la hoja 4 y calculo de ICA
calculo.ica.1 <- base.datos%>% select(RODAL, T.P,EDAD, DIAMETRO) %>%
                                group_by(RODAL)%>%
                                summarise(  T.P =sum(T.P, na.rm= TRUE),
                                            EDAD = sum(EDAD, na.rm =TRUE),
                                            DIAMETRO =max(DIAMETRO, na.rm=TRUE))  #se construyen estadisticos 

calculo.ica.2 <- cuadro3 %>% select(RODAL, Especie,er.totales)%>%filter(Especie=='Pinus engelmannii') # se filtran los datos calculado en el cuadro 3 para una especie
calculo.ica.2 <- calculo.ica.2 %>% rename(ERT=er.totales)
calculo.ica <- merge(calculo.ica.1, calculo.ica.2)
calculo.ica <- calculo.ica %>% mutate(ICA = (ERT*10)/(T.P*DIAMETRO),
                       ICA_porcentaje = (ICA/ERT)*100,
                       IMA=ERT/EDAD)
cuadro4 <- cuadro3%>% select(RODAL,Especie,er.entre.ha, er.totales, ab.entre.ha,
                             IC,M3.VTA.por.ha.verde, Area.basal.M2.por.ha,
                             M3.VTA.por.ha, M3.VTA.por.UMM)
cuadro4 <- cuadro4%>% mutate(Sin.nombre= (er.totales- M3.VTA.por.UMM ))
cuadro4.1.1 <- cuadrox %>%select(c20a25, c30a40, c5a15, er.entre.ha)
cuadro4.1.1 <- cuadro4.1.1%>%rename(Total = er.entre.ha )
cuadro4.1.1 <- cbind(cuadro4, cuadro4.1.1)
cuadro4.1.1 <- merge(cuadro4.1.1, base.superficie) 
cuadro4.1.1 <- cuadro4.1.1 %>% mutate(ER.entre.TOTALES.M = Total*Superficie)
cuadro4.1.1<- merge(cuadro4.1.1, calculo.ica, all.x=TRUE)
cuadro4.1.1$ICA[is.na(cuadro4.1.1$ICA)] <- 0
cuadro4 <- cuadro4.1.1%>%select(Especie, er.entre.ha, er.totales, ab.entre.ha, IC, 
                  M3.VTA.por.ha.verde, Area.basal.M2.por.ha,
                  M3.VTA.por.ha, M3.VTA.por.UMM, Sin.nombre, c20a25,         
                  c30a40, c5a15, Total, ER.entre.TOTALES.M,
                  ICA)
cuadr4x <- cuadro4.1.1 #copio este cuadro para tomar valores de el y construir el cuadro 5
### se termina cuadro 4
####se comienza a construir cuadro5
cuadro5 <- base.superficie[order(base.superficie$RODAL),] #tomo la superficie de la hoja 'SUP RODAL' y la ordeno
    #esto es lo que esta en la hoja 'ABxHA'
cuadro5.1 <- base.datos%>% select(RODAL, No, NOMBRE.CIENTIFICO,AB)%>%
              group_by(RODAL, NOMBRE.CIENTIFICO)%>%
              summarise(Suma.de.AB=sum(AB), Suma.de.No=sum(No))
cuadro5.2 <- base.datos %>% select(RODAL,SITIO) %>%group_by(RODAL) %>% summarise(No.Sitios=length(unique(SITIO)) /10 ) #en porcentaje
cuadro5.1 <- merge(cuadro5.1, cuadro5.2, all.x = TRUE)
cuadro5.1 <- cuadro5.1 %>% mutate(AB.entre.HA=Suma.de.AB/No.Sitios,
                     NARBOL.entre.HA=Suma.de.No/No.Sitios )
    #esto es lo de la hoja 11
cuadro5.1 <- cuadro5.1 %>% group_by(RODAL)%>% 
                          summarise(Suma.de.AB.entre.HA= sum(AB.entre.HA),
                                    Suma.de.NARBOL.entre.HA= sum(NARBOL.entre.HA))
cuadro5 <- merge(cuadro5, cuadro5.1, all.y =TRUE)
cuadro5 <- cuadro5%>%rename(No.de.Arboles.entre.ha=Suma.de.NARBOL.entre.HA,
                            Area.Basal=Suma.de.AB.entre.HA)
cuadro5.3 <- calculo.ica%>% select(RODAL, T.P, IMA) 
cuadro5 <- merge(cuadro5, cuadro5.3, all.x = TRUE)
cuadro5 <- cuadro5 %>% rename(Tiempo.de.paso=T.P)
cuadro5 <- cuadro5 %>% mutate(ICA.m3.ha.ano.=(Tiempo.de.paso*10)/(Area.Basal*Superficie))
cuadro5 <- cuadro5%>%select(RODAL, Superficie, No.de.Arboles.entre.ha,
                 Area.Basal, Tiempo.de.paso, ICA.m3.ha.ano.,IMA)
cuadro5 <- cuadro5[order(cuadro5$RODAL),] #tenemos diferencias entre los tiempos de paso abria que checarlo 
####se termina cuadro 5
#####se comienza cuadro 6
cuadro6 <- cuadro3 %>% select(RODAL, Superficie, Especie, M3.VTA.por.UMM)
cuadro6 <- cuadro6 %>%rename(No= RODAL, Genero = Especie, m3VTA = M3.VTA.por.UMM)
suma <- sum(cuadro6$m3VTA) /5 #Jaime me dijo que quiere cinco niveles
cuadro6$Tratamiento.Silvicola <- 'Seleccion'
cuadro6$Volumen.por.infraestructura <- 0
cuadro6$Posibilidad.mas.columen.por.infraestructura <- cuadro6$m3VTA+cuadro6$Volumen.por.infraestructura
corte <- cuadro6%>% mutate(corte =cumsum(m3VTA))
corte$comparacion <- cut(corte$corte, breaks = suma*0:5)
    #los niveles los tuve que escoger a mano 
corte <- corte%>% mutate(Area.de.corta = '1')
corte$Area.de.corta[corte$No<=14] <- '1 (2018-2019)'
corte$Area.de.corta[corte$No>14 & corte$No<=28] <- '2 (2020-2021)'
corte$Area.de.corta[corte$No>28 & corte$No<=36] <- '3 (2022-2023)'
corte$Area.de.corta[corte$No>36 & corte$No<=59] <- '4 (2024-2025)'
corte$Area.de.corta[corte$No>59] <- '5 (2026-2027)'
corte%>%group_by(Area.de.corta)%>%summarise(sum(m3VTA)) # medio estan balanceados
cuadro6 <- corte %>% select(Area.de.corta, No, Superficie,
                 Tratamiento.Silvicola, m3VTA,
                 Volumen.por.infraestructura, Posibilidad.mas.columen.por.infraestructura)
###se termina cuadro 6
###escritura a excel
library(xlsx) # paquete de excel que puede crear problemas porque depende de java
write.xlsx(cuadro3, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro3')# tendras que cambiar la ruta de escritura
write.xlsx(cuadro4.1, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro4.1 NOM', append = TRUE)# tendras que cambiar la ruta de escritura
write.xlsx(cuadro4.2, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro4.2 NOM', append = TRUE)# tendras que cambiar la ruta de escritura
write.xlsx(cuadro4, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro4', append = TRUE)# tendras que cambiar la ruta de escritura
write.xlsx(cuadro5, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro5', append = TRUE)# tendras que cambiar la ruta de escritura
write.xlsx(cuadro6, file ='C:/Users/fou-f/Desktop/ta Garcia/cuadros.xlsx',
           sheetName = 'cuadro6', append = TRUE)# tendras que cambiar la ruta de escritura
