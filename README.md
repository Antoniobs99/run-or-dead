# run-or-dead
Improving code for data analysis from the Allen Brain Atlas database

```{r}
#Librerias que voy a ocupar
#HK= HouseKeeping
#DT= data.table
library(QuickAndDirtyABAAnalisis)
library(ABAData)
library(BiocManager)
library(ABAEnrichment)
library(tidyverse)
library(data.table)
```
```{r}
#Primer filtrado de genes que vamos a ocupar
Genes_IPC<- GetAreasGenes(structure.id = 10225,selected.dataset = "dataset_5_stages")
Genes_VFC<- GetAreasGenes(structure.id = 10185,selected.dataset = "dataset_5_stages")
```
```{r}
#Segundo filtrado paras solo quedarnos con el promedio de las 5 etapas para cada gen, por lo que nos da solo un arrow por gen y dos variables, su nombre y el promedio
#Fue "senal" Porque me marcaba error con ñ jsjs
primer_filtrado_IPC<-Genes_IPC%>%
  group_by(hgnc_symbol)%>%
  summarise(promedio_senal<-mean(signal))
primer_filtrado_VFC<-Genes_VFC%>%
  group_by(hgnc_symbol)%>%
  summarise(promedio_senal<-mean(signal))
```
```{r}
#Los volvemos DT
primer_filtrado_IPC_DT<-setDT(primer_filtrado_IPC)
primer_filtrado_VFC_DT<-setDT(primer_filtrado_VFC)
```
```{r}
#Tercer filtrado, promedio mañor a 7, algunos bocenos de los calculos estadisticos utilizados para elegir arriba (o igual) a 7, los puse hasta abajo
Expresion_Mayor7_IPC<-primer_filtrado_IPC_DT%>%
  filter(`promedio_senal <- mean(signal)`>=7)
Expresion_Mayor7_VFC<-primer_filtrado_VFC_DT%>%
  filter(`promedio_senal <- mean(signal)`>=7)
```
```{r}
HK_genes <- unique(HK_exons [ ,"Gene Name"] )#Aqui Importé la base de datos de genes HK sacados de una fuente de wikipedia (Espero eso no le quite rigurosidad) y como algunos estaban repetidos hicimos esto porque por tidyverse no me salio :p
HK_genesDT<-setDT(HK_genes) #Y lo volvemos DT
```
```{r}
#Convertimos nuestras variables a DT
Expresion_Mayor7_IPCDT<-setDT(Expresion_Mayor7_IPC)
Expresion_Mayor7_VFCDT<-setDT(Expresion_Mayor7_VFC)
#Por ultimo sacamos los genes que son HK y listo! Cada variable es por las zonas que vamos a utilizar, promedio del gen en las 5 etapas, con expresión mayor a 7 y que no sean HouseKeeping
IPC_sinHK <-Expresion_Mayor7_IPCDT[!hgnc_symbol %in% HK_genesDT$`Gene Name`]
VFC_sinHK <-Expresion_Mayor7_VFCDT[!hgnc_symbol %in% HK_genesDT$`Gene Name`]
```
#Algunos calculos estadisticos irrelevantes para el codigo pero aun así considero necesario ponerlos
```{r}
mean_1filtrado<-mean(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)
rango_1filtrado<-range(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)
max_1filtrado<-max(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)
min_1filtrado<-min(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)#it´s 0 bc are number very small
var_1filtado<-var(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)
sd_1filtrado<-sd(primer_filtrado_IPC$`promedio_senal <- mean(signal)`)
```
```{r}
valor_MECP2<-primer_filtrado_IPC_DT%>%
  filter(hgnc_symbol=="MECP2")
```
```{r}
ZConValOfMECP2<- -0.11
#Valor acomulado al -Infinito= 0.45 
#0.45-0.5= 0.05 Por lo tanto se lo agregamos a .5 y nos da el acumulado apartir de 7 para adelante
#Esto es, el 55% de los datos.
```
