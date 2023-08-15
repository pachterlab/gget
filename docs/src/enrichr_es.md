> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget enrichr 💰
Realice un análisis de enriquecimiento de conjuntos de genes en una lista de genes utilizando [Enrichr](https://maayanlab.cloud/Enrichr/).  
Regresa: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  
  
**Parámetro posicional**  
`genes`  
Lista de nombres cortos (símbolos) de los genes de interés con que realizar la análisis de enriquecimiento, p. PHF14 RBM3 MSL1 PHF21A.  
Alternativamente: usa la bandera `--ensembl` para ingresar IDs del tipo Ensembl, p. ENSG00000106443 ENSG00000102317 ENSG00000188895.  

**Otros parámetros requeridos**  
`-db` `--database`  
Base de datos a utilizar como referencia para el análisis de enriquecimiento.  
Admite cualquier base de datos enumerada [aquí](https://maayanlab.cloud/Enrichr/#libraries) o uno de los siguientes accesos directos:  
'pathway'       (KEGG_2021_Human)  
'transcription'     (ChEA_2016)  
'ontology'      (GO_Biological_Process_2021)  
'diseases_drugs'   (GWAS_Catalog_2019)   
'celltypes'      (PanglaoDB_Augmented_2021)  
'kinase_interactions'   (KEA_2015)  
  
**Parámetros optionales**  
`-bkg_l` `--background_list`  
Lista de nombres cortos (símbolos) de genes de 'background' (de fondo/control), p. NSUN3 POLRMT NLRX1.  
Alternativamente: usa la bandera `--ensembl_background` para ingresar IDs del tipo Ensembl.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

`figsize`  
Solo para Python. (ancho, alto) de la visualización en pulgadas. (Por defecto: (10,10))

`ax`   
Solo para Python. Ingresa un objeto de ejes matplotlib para personalizar la visualización.(Por defecto: None)  

  
**Banderas**  
`-e` `--ensembl`   
Usa esta bandera si `genes` se ingresa como una lista de IDs del tipo Ensembl.     

`-e_b` `--ensembl_background`  
Usa esta bandera si `background_list` se ingresa como una lista de IDs del tipo Ensembl.  

`-bkg` `--background`  
Use un conjunto de 20,625 genes 'background' de [https://maayanlab.cloud/Enrichr/](https://maayanlab.cloud/Enrichr/).
 
`-csv` `--csv`  
Solo para la Terminal. Regresa los resultados en formato CSV.    
Para Python, usa `json=True` para regresar los resultados en formato JSON.   

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  
  
`plot`  
Solo para Python. `plot=True` provee una visualización de los primeros 15 resultados (por defecto: False).  
  
  
### Por ejemplo
```bash
gget enrichr -db ontology ACE2 AGT AGTR1
```
```python
# Python
gget.enrichr(["ACE2", "AGT", "AGTR1"], database="ontology", plot=True)
```
&rarr; Regresa vías/funciones celulares relacionadas con los genes ACE2, AGT y AGTR1 de la base de datos *GO Biological Process 2021*. En Python, `plot=True` provee una visualización de los resultados:

![alt text](https://github.com/pachterlab/gget/blob/main/figures/gget_enrichr_results.png?raw=true)

<br/><br/>

**Use `gget enrichr` con una lista de genes 'background':**  
```bash
# Aquí, ingresamos los genes de interés (parámetro posicional 'genes') primero, para que no se agregan a la lista de genes 'background' detrás del parámetro '-bkgr_l'
gget enrichr \
	PHF14 RBM3 MSL1 PHF21A ARL10 INSR JADE2 P2RX7 LINC00662 CCDC101 PPM1B KANSL1L CRYZL1 ANAPC16 TMCC1 CDH8 RBM11 CNPY2 HSPA1L CUL2 PLBD2 LARP7 TECPR2 ZNF302 CUX1 MOB2 CYTH2 SEC22C EIF4E3 ROBO2 ADAMTS9-AS2 CXXC1 LINC01314 ATF7 ATP5F1 \
	-db ChEA_2022 \
	-bkg_l NSUN3 POLRMT NLRX1 SFXN5 ZC3H12C SLC25A39 ARSG DEFB29 PCMTD2 ACAA1A LRRC1 2810432D09RIK SEPHS2 SAC3D1 TMLHE LOC623451 TSR2 PLEKHA7 GYS2 ARHGEF12 HIBCH LYRM2 ZBTB44 ENTPD5 RAB11FIP2 LIPT1 INTU ANXA13 KLF12 SAT2 GAL3ST2 VAMP8 FKBPL AQP11 TRAP1 PMPCB TM7SF3 RBM39 BRI3 KDR ZFP748 NAP1L1 DHRS1 LRRC56 WDR20A STXBP2 KLF1 UFC1 CCDC16 9230114K14RIK RWDD3 2610528K11RIK ACO1 CABLES1 LOC100047214 YARS2 LYPLA1 KALRN GYK ZFP787 ZFP655 RABEPK ZFP650 4732466D17RIK EXOSC4 WDR42A GPHN 2610528J11RIK 1110003E01RIK MDH1 1200014M14RIK AW209491 MUT 1700123L14RIK 2610036D13RIK PHF14 RBM3 MSL1 PHF21A ARL10 INSR JADE2 P2RX7 LINC00662 CCDC101 PPM1B KANSL1L CRYZL1 ANAPC16 TMCC1 CDH8 RBM11 CNPY2 HSPA1L CUL2 PLBD2 LARP7 TECPR2 ZNF302 CUX1 MOB2 CYTH2 SEC22C EIF4E3 ROBO2 ADAMTS9-AS2 CXXC1 LINC01314 ATF7 ATP5F1COX15 TMEM30A NSMCE4A TM2D2 RHBDD3 ATXN2 NFS1 3110001I20RIK BC038156 C330002I19RIK ZFYVE20 POLI TOMM70A LOC100047782 2410012H22RIK RILP A230062G08RIK PTTG1IP RAB1 AFAP1L1 LYRM5 2310026E23RIK SLC7A6OS MAT2B 4932438A13RIK LRRC8A SMO NUPL2
```
```python
# Python
gget.enrichr(
	genes = [
		"PHF14", "RBM3", "MSL1", "PHF21A", "ARL10", "INSR", "JADE2", "P2RX7",
		"LINC00662", "CCDC101", "PPM1B", "KANSL1L", "CRYZL1", "ANAPC16", "TMCC1",
		"CDH8", "RBM11", "CNPY2", "HSPA1L", "CUL2", "PLBD2", "LARP7", "TECPR2", 
		"ZNF302", "CUX1", "MOB2", "CYTH2", "SEC22C", "EIF4E3", "ROBO2",
		"ADAMTS9-AS2", "CXXC1", "LINC01314", "ATF7", "ATP5F1"
	], 
	database = "ChEA_2022",
	background_list = [
		"NSUN3","POLRMT","NLRX1","SFXN5","ZC3H12C","SLC25A39","ARSG",
		"DEFB29","PCMTD2","ACAA1A","LRRC1","2810432D09RIK","SEPHS2",
		"SAC3D1","TMLHE","LOC623451","TSR2","PLEKHA7","GYS2","ARHGEF12",
		"HIBCH","LYRM2","ZBTB44","ENTPD5","RAB11FIP2","LIPT1",
		"INTU","ANXA13","KLF12","SAT2","GAL3ST2","VAMP8","FKBPL",
		"AQP11","TRAP1","PMPCB","TM7SF3","RBM39","BRI3","KDR","ZFP748",
		"NAP1L1","DHRS1","LRRC56","WDR20A","STXBP2","KLF1","UFC1",
		"CCDC16","9230114K14RIK","RWDD3","2610528K11RIK","ACO1",
		"CABLES1", "LOC100047214","YARS2","LYPLA1","KALRN","GYK",
		"ZFP787","ZFP655","RABEPK","ZFP650","4732466D17RIK","EXOSC4",
		"WDR42A","GPHN","2610528J11RIK","1110003E01RIK","MDH1","1200014M14RIK",
		"AW209491","MUT","1700123L14RIK","2610036D13RIK",
		"PHF14", "RBM3", "MSL1", "PHF21A", "ARL10", "INSR", "JADE2", 
		"P2RX7", "LINC00662", "CCDC101", "PPM1B", "KANSL1L", "CRYZL1", 
		"ANAPC16", "TMCC1","CDH8", "RBM11", "CNPY2", "HSPA1L", "CUL2", 
		"PLBD2", "LARP7", "TECPR2", "ZNF302", "CUX1", "MOB2", "CYTH2", 
		"SEC22C", "EIF4E3", "ROBO2", "ADAMTS9-AS2", "CXXC1", "LINC01314", "ATF7", 
		"ATP5F1""COX15","TMEM30A","NSMCE4A","TM2D2","RHBDD3","ATXN2","NFS1",
		"3110001I20RIK","BC038156","C330002I19RIK","ZFYVE20","POLI","TOMM70A",
		"LOC100047782","2410012H22RIK","RILP","A230062G08RIK",
		"PTTG1IP","RAB1","AFAP1L1", "LYRM5","2310026E23RIK",
		"SLC7A6OS","MAT2B","4932438A13RIK","LRRC8A","SMO","NUPL2"
	],
	plot=True
)
```
&rarr; Regresa factores de transcripción relacionados con los genes de interés y controlados con la lista de genes background de la base de datos *ChEA 2022*. En Python, `plot=True` provee una visualización de los resultados:

![alt text](https://github.com/pachterlab/gget/blob/main/figures/gget_enrichr_results_2.png?raw=true)

<br/><br/>

El siguiente ejemplo fue enviado por [Dylan Lawless](https://github.com/DylanLawless) través de un [PR](https://github.com/pachterlab/gget/pull/54) (con ajustes de [Laura Luebbert](https://github.com/lauraluebbert)):  
**Use `gget enrichr` en R y cree unq visualización similar usando [ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html).**  
TENGA EN CUENTA el cambio de ejes en comparación con la visualización en Python. 
```r
system("pip install gget")
install.packages("reticulate")
library(reticulate)
gget <- import("gget")

# Perform enrichment analysis on a list of genes
df <- gget$enrichr(list("ACE2", "AGT", "AGTR1"), database = "ontology")

# Count number of overlapping genes
df$overlapping_genes_count <- lapply(df$overlapping_genes, length) |> as.numeric()

# Only keep the top 15 results
df <- df[1:15, ]

# Plot
library(ggplot2)

df |>
	ggplot() +
	geom_bar(aes(
		x = -log10(adj_p_val),
		y = reorder(path_name, -adj_p_val)
	),
	stat = "identity",
  	fill = "lightgrey",
  	width = 0.5,
	color = "black") +
	geom_text(
		aes(
			y = path_name,
			x = (-log10(adj_p_val)),
			label = overlapping_genes_count
		),
		nudge_x = 0.75,
		show.legend = NA,
		color = "red"
	) +
  	geom_text(
		aes(
			y = Inf,
			x = Inf,
      			hjust = 1,
      			vjust = 1,
			label = "# of overlapping genes"
		),
		show.legend = NA,
		color = "red"
	) +
	geom_vline(linetype = "dotted", linewidth = 1, xintercept = -log10(0.05)) +
	ylab("Pathway name") +
	xlab("-log10(adjusted P value)")
```

#### [Más ejemplos](https://github.com/pachterlab/gget_examples)
