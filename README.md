# CARNIVAL

CARNIVAL is an R-package providing a framework to perform causal reasoning to infer a subset of signalling network from transcriptomics data. This work was originally based on [Melas et al.](https://pubs.rsc.org/en/content/articlehtml/2015/ib/c4ib00294f) with a number improved functionalities comparing to the original version.
Transcription factors’ (TFs) activities and pathway scores from gene expressions can be inferred with our in-house tools [DoRothEA](https://github.com/saezlab/DoRothEA) & [PROGENy](https://github.com/saezlab/progeny), respectively.
TFs’ activities and signed directed protein-protein interaction networks +/- drug targets and pathway scores are then used to derive a series of linear constraints to generate integer linear programming (ILP) problems. 
An ILP solver (CPLEX) is subsequently applied to identify the sub-network topology with minimised discrepancies on fitting error and model size.

More detailed descriptions of CARNIVAL, benchmarking and applicational studies can be found on it's dedicated [web-page](https://saezlab.github.io/CARNIVAL/) and in [Liu, Trairatphisan, Gjerga et al.](https://www.nature.com/articles/s41540-019-0118-z):

> Liu A*, Trairatphisan P*, Gjerga E*, Didangelos A, Barratt J, Saez-Rodriguez J. (2019). From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL. *npj Systems Biology and Applications*, https://doi.org/10.1038/s41540-019-0118-z (*equal contributions).


## Getting Started

A tutorial for preparing CARNIVAL input files starting from differentially gene expression (DEG) and for running the CARNIVAL pipeline are provided as vignettes in R-Markdown, R-script and HTML formats. The wrapper script "runCARNIVAL" was introduced to take input arguments, pre-process input descriptions, run optimisation and export results as network files and figures. Three built-in CARNIVAL examples are also supplied as case studies for users.

### Prerequisites

CARNIVAL requires the interactive version of IBM Cplex or CBC-COIN solver as the network optimiser. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio?S_PKG=CoG&cm_mmc=Search_Google-_-Data+Science_Data+Science-_-WW_IDA-_-+IBM++CPLEX_Broad_CoG&cm_mmca1=000000RE&cm_mmca2=10000668&cm_mmca7=9041989&cm_mmca8=kwd-412296208719&cm_mmca9=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_&cm_mmca10=267798126431&cm_mmca11=b&mkwid=_k_Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB_k_|470|135655&cvosrc=ppc.google.%2Bibm%20%2Bcplex&cvo_campaign=000000RE&cvo_crid=267798126431&Matchtype=b&gclid=Cj0KCQiAr93gBRDSARIsADvHiOpDUEHgUuzu8fJvf3vmO5rI0axgtaleqdmwk6JRPIDeNcIjgIHMhZIaAiwWEALw_wcB). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available for any user. Alternatively for smaller cases, users can rely on the freely available [lpSolve R-package](https://cran.r-project.org/web/packages/lpSolve/index.html).

### Installing

CARNIVAL is currently available for the installation as an R-package from our GitHub page

```R
# Install CARNIVAL from Github using devtools
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/CARNIVAL', build_vignettes = TRUE)
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CARNIVAL_directory', repos = NULL, type="source")
```

## Running CARNIVAL

To obtain the list of tutorials about how to run CARNIVAL over a few toy examples, user can start with typing the following commmand on R-console:

```R
vignette("CARNIVAL-vignette")
```

### Real case example
Now we show an example about how we can run CARNIVAL on a real case-study. We will show step-by-step how we can build a prior knowledge of signalling and a DoRothEA regulon list from [OmniPath](https://github.com/saezlab/OmnipathR). Next we show how we can estimate pathway activities via [PROGENy](https://github.com/saezlab/progeny). Next we will demonstrate how we can use the [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html for the R-package). Finally we will show how we can use CARNIVAL to combine all this information in order to infer the causal regulatory network.

Through OmnipathR package, we build the network object needed for contextualizing the signalling network with CARNIVAL.

```R
library(OmnipathR) # load OmnipathR library
##importing interactions from SignaLink3, PhosphoSite and Signor
interactions <- import_Omnipath_Interactions(filter_databases=c("SignaLink3","PhosphoSite","Signor"))
## keeping the directed interactions
interactions <- interactions[which(interactions$is_directed==1), ]
# keeping signed interactions only
interactions <- interactions[which((interactions$is_stimulation+interactions$is_inhibition)==1), ]
## now building the network data-frame and keeping it as netObj
network <- matrix(data = , nrow = nrow(interactions), ncol = 3)
network[, 1] <- interactions$source_genesymbol
network[which(interactions$is_stimulation==1), 2] <- "1"
network[which(interactions$is_inhibition==1), 2] <- "-1"
network[, 3] <- interactions$target_genesymbol
netObj <- as.data.frame(network)
```

For the analysis we consider the expression data example from the [progeny](https://github.com/saezlab/progeny). We estimate normalized (from -1 to 1) pathway activity scores via the progeny function for the first sample and then we assign the inferred activities as node weights to the progeny members from the progenyMembers object loaded from CARNIVAL package. Users can estimate the pathway activities for any sample they wish (see documentation of _assignPROGENyScores()_ function)

```R
library(progeny)
library(CARNIVAL)
expr <- as.matrix(read.csv(system.file("extdata", "human_input.csv", package = "progeny"), row.names = 1))
human_def_act <- progeny(expr, scale = TRUE, organism = "Human", top = 100, perm = 10000, z_scores = FALSE)
## loading the progeny members to assign the weights
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))
## now assigning the PROGENy weights to pathway members only for the first
## sample which we can consider for the CARNIVAL analysis
weightObj <- assignPROGENyScores(progeny = human_def_act, progenyMembers = progenyMembers, id = "gene", access_idx = 1)
```

Next we can retrieve regulons from OmniPath and create the regulon table. From there we obtain the viper regulon list through the createRegulonList function (see documentation for more details).

```R
regulon_df <- import_TFregulons_Interactions(select_organism = 9606)
regulon_df <- regulon_df[which((regulon_df$is_stimulation+regulon_df$is_inhibition)==1), ]
regulon_table <- matrix(data = , nrow = nrow(regulon_df), ncol = 3)
regulon_table[, 1] <- regulon_df$source_genesymbol
regulon_table[which(regulon_df$is_stimulation==1), 2] = "1"
regulon_table[which(regulon_df$is_inhibition==1), 2] = "-1"
regulon_table[, 3] <- regulon_df$target_genesymbol
regulons <- createRegulonList(regulon_table = regulon_table)
```

From here we can estimate the TF activities and generate the continuos measurement inputs for CARNIVAL.

```R
library(viper)
TF_activities = as.data.frame(viper::viper(eset = expr, regulon = regulons, nes = TRUE, method = 'none', minsize = 4, eset.filter = FALSE))
tfList <- generateTFList(df = TF_activities, top = "all", access_idx = 1)
```

So far we have generated the prior knowledge network which will be used to contextualize the regulatory signalling network (_netObj_) as well as a list object containing the progeny weights (_weightObj_) and the TF activities (tfList to be assigned to the _measObj_). For running CARNIVAL we have to access one element at the time from the weightObj and tfList. We use all this as inputs to run the CARNIVAL analysis through the _runCARNIVAL()_ function. Since the input targets are not known, we perform the invCARNIVAL analysis.

```R
##users need to define the path to the cplex solver
res <- runCARNIVAL(measObj = tfList[[1]], netObj = netObj, weightObj = weightObj[[1]], solverPath = solverPath, solver = "cplex", DOTfig = TRUE)
```

The solution network of this small real-case application will be also saved as a DOT figure in the working directory.

## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/CARNIVAL/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).

## References

[Melas et al.](https://pubs.rsc.org/en/content/articlehtml/2015/ib/c4ib00294f):

> Melas IN, Sakellaropoulos T, Iorio F, Alexopoulos L, Loh WY, Lauffenburger DA, Saez-Rodriguez J, Bai JPF. (2015). Identification of drug-specific pathways based on gene expression data: application to drug induced lung injury. *Integrative Biology*, Issue 7, Pages 904-920, https://doi.org/10.1039/C4IB00294F.

[DoRothEA v2 - Garcia-Alonso et al.](https://www.biorxiv.org/content/early/2018/06/03/337915):

> Garcia-Alonso L, Ibrahim MM, Turei D, Saez-Rodriguez J. (2018). Benchmark and integration of resources for the estimation of human transcription factor activities. *bioRXiv*, https://doi.org/10.1101/337915.

[PROGENy - Schubert et al.](https://www.nature.com/articles/s41467-017-02391-6):

> Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Blüthgen N, Saez-Rodriguez J. (2018). Perturbation-response genes reveal signaling footprints in cancer gene expression. *Nature Communication*, Issue 9, Nr. 20. https://doi.org/10.1038/s41467-017-02391-6.


## Acknowledgement

CARNIVAL has been developed as a computational tool to analyse -omics data within the [TransQST Consortium](https://transqst.org) and [H2020 Symbiosys ITN Training Network](https://www.h2020symbiosys.eu/).

"This project has received funding by the European Union’s H2020 program (675585 Marie-Curie ITN ‘‘SymBioSys’’) and the Innovative Medicines Initiative 2 Joint Undertaking under grant agreement No 116030. The Joint Undertaking receives support from the European Union's Horizon 2020 research and innovation programme and EFPIA."
