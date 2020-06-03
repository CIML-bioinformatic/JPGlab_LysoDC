# Differentiation paths of Peyer’s patch LysoDC are linked to sampling site positioning, migration and T cell priming 


## Article information

**Title:** Differentiation paths of Peyer’s patch LysoDC are linked to sampling site positioning, migration and T cell priming 

**Authors:** Camille Wagner 1,&, Johnny Bonnardel 1,4,5,&, Clément Da Silva 1,6,&, Lionel Spinelli 1, Cynthia Arroyo Portilla 1,2, Julie Tomas 1, Margaux Lagier 1, Lionel Chasson 1, Marion Masse 1, Marc Dalod 1, Alexandre Chollat-Namy 3, Jean- Pierre Gorvel 1 and Hugues Lelouard 1,*

1 Aix Marseille Univ, CNRS, INSERM, CIML, Marseille, France.

2 Departamento de Análisis Clínicos, Facultad de Microbiología, Universidad de Costa Rica, San José, Costa Rica.

3 Pathological anatomy and cytology laboratory, Marseille, France.

4 Current address: Department of Biomedical molecular Biology, Ghent University, Ghent, Belgium.

5 Current address: Laboratory of Myeloid Cell Ontogeny and Functional Specialisation, VIB Center for Inflammation Research, Ghent, Belgium.

6 Current address: Immunology Section, Lund University, Lund 221 84, Sweden.

& these authors contributed equally to this work.

\* correspondence: Hugues Lelouard, Centre d'Immunologie de Marseille-Luminy, Marseille, 13288, France. Phone: + 33 491269416; Fax: + 33 491269430; e-mail: lelouard@ciml.univ- mrs.fr 

**Abstract:**
The specialized monocyte-derived phagocytes termed LysoDC are hallmarks of Peyer’s patches where their main function is to sample microrganisms and detect pathogens. Here, we studied their differentiation pathway in relation with their sampling, migratory and T cell priming abilities. Among the four identified LysoDC differentiation stages that display similar phagocytic activity, one is located in follicles in a CXCR5- independent manner whereas the others reside in subepithelial domes (SED) where they proliferate and then mature, as they get closer to the epithelium. Fully mature LysoDC but not macrophages express a set of genes in common with conventional DC, and prime naïve T cells in vitro. However, at steady state, they do not migrate in naïve T cell-enriched interfollicular regions (IFR) but, upon stimulation, they express the chemokine receptor CCR7 and migrate from the SED to the IFR periphery where they strongly interact with proliferative immune cells, including helper T cells. Finally, we show that LysoDC populate human Peyer’s patches, strengthening their interest as targets for modulating intestinal immunity. 

---
---

## Goal of the github
This github project contains the instructions and material to reproduce the analysis reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Required data and builded Docker images are available on download. Intructions to reproduce the analysis are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the docker images tar file
- Load the docker images on your system
- Download the pre-processed data (Cell Ranger results) and the processed data (analysis results)

Below you will find detailed instruction for each of these steps.

#### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder "LysoDC" with all the source code. 
You must set an environment variable called WORKING_DIR with a value set to the path to this folder.

For instance, if I have chosen to clone the Git repository in "/home/spinellil/workspace", then the WORKING_DIR variable will be set to "/home/spinellil/workspace/LysoDC"

**On linux:**

    export WORKING_DIR=/home/spinellil/workspace/LysoDC

#### Download the docker images

Docker images tar file are stored on Zenodo. Open a shell command and change dir to the root of the cloned Git repository. Then execute the following commands to download the images tar files to the right project folder:

    wget URL_jpglab_lysodc_seurat.tar -o $WORKING_DIR/docker/image/jpglab_lysodc_seurat/jpglab_lysodc_seurat.tar
    wget URL_jpglab_lysodc_rnavelocity.tar -o $WORKING_DIR/docker/image/jpglab_lysodc_rnavelocity/jpglab_lysodc_rnavelocity.tar
    wget URL_jpglab_lysodc_pseudotime.tar -o $WORKING_DIR/docker/image/jpglab_lysodc_pseudotime/jpglab_lysodc_pseudotime.tar

#### Load docker images

In order to execute analysis, you must load the provided docker images onto your Docker. Docker must be installed on your system. 
See https://docs.docker.com/install/ for details on Docker installation.
Open a shell command, change dir to the folder in which you cloned the project and type:

    docker load -i $WORKING_DIR/docker/image/jpglab_lysodc_seurat/jpglab_lysodc_seurat.tar
    docker load -i $WORKING_DIR/docker/image/jpglab_lysodc_rnavelocity/jpglab_lysodc_rnavelocity.tar
    docker load -i $WORKING_DIR/docker/image/jpglab_lysodc_pseudotime/jpglab_lysodc_pseudotime.tar

Those commands may take some time. If you encounter an issue loading some docker image layer, try again. Sometimes issue would be resolved. 

#### Download the pre-processed data (Cell Ranger results) and the processed data (analysis results)

Sequencing data are available for download on GEO but are not used here. Pre-processed data (result from Cell Ranger counts) are also available on GEO. Annotation data and processed data (result of the analysis you will reproduce here) are available on Zenodo. In order to simplify the procedure, required data to reproduce the analysis have been packaged on the Zenodo repository. To dowload all those data, use the following command:

---
---

## Run the analysis

Analysis can be directly run inside docker containers by compiling Rmarkdown files. The Rmarkdown file comilation will launch the required analysis for the step
and produce a final HTML report. 

Each step report will be generated in the <WORKING_DIR>/data/step<X>/output where <WORKING_DIR> is the folder where you clone the git repository (and set the WORKING_DIR environment variable) and <X> is the number of the analysis step you are executing (1 or 5 to 8).

---

### Run step 1

#### Goal
This step is a standard quality control and first analysis of the cell heterogeneity. 

#### Output 
This step output files that contain cellular barcodes associated to cell do not passing
one of the QC test. It also output two matrices of gene expression o,e with the raw UMI counts, and one with the normalized UMI counts. 
Finally, it output the HTML report of the analysis.

#### Execution
To run the step1, ensure you have correctly downloaded the data in the folder <WORKING_DIR>/data/raw and run the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_seurat R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step1/JPGlab_LysoDC_PrimaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step1/output"), output_file = "JPGlab_LysoDC_PrimaryAnalysis.html", quiet=FALSE)'

#### Results
Once the analysis done, here is the tree the WORKING_DIR/data/step1 folder you may obtain (with the newly created "output" folder):

    step1
    └── output
        ├── excluded_cells_contamination.txt            # The list of cellular barcode of cells consired as contamination
        ├── excluded_cells_HighMitoGenePerc.txt         # The list of cellular barcode of cells with too high percentage of Mitoconcrial genes
        ├── excluded_cells_LowGeneNb.txt                # The list of cellular barcode of cells with too low number of genes
        ├── excluded_cells_LowUMINb.txt                 # The list of cellular barcode of cells with too low number of UMI
        ├── excluded_cells_proliferation.txt            # The list of cellular barcode of cells identified as proliferating
        ├── filtered_normalized_expression_matrix.csv   # The matrix of normalized gene expression
        ├── filtered_raw_expression_matrix.csv          # The matrix of raw gene expression (UMI counts)
        └── JPGlab_LysoDC_PrimaryAnalysis.html          # The HTML report of the analysis

---

### Run step 5

#### Goal
This step aims to produce a first version of the RNA velocity analysis of the data. This step is also used to identify cell considered as contamination or "in dying" state. These cells are removed from the analysis to improve the quality of the next steps.

#### Input
This step use input from previous analysis. Here is a tree présentation of the content of the folder $WORKING_DIR/data/step5 as it should be before the execution of the analysis.

    step5
    └── input
        ├── excluded_cells_contamination.txt    # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_contamination.txt
        ├── excluded_cells_HighMitoGenePerc.txt # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_HighMitoGenePerc.txt
        ├── excluded_cells_LowGeneNb.txt        # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowGeneNb.txt
        ├── excluded_cells_LowUMINb.txt         # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowUMINb.txt
        └── excluded_cells_proliferation.txt    # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_proliferation.txt

#### Output
This step output a loom file containing the Velocyto analysis of the bam file produced by CellRanger (see above). It also produce a
file containg the result of cell clustering.

**Important note:** The clustering in this step is performed thanks to the package Pagoda2, using the _makeKnnGraph_ function (on PCA space) and the _getKnnClusters_ using the 'multilevel community' method. This method is based on a Louvain clustering with maximization of the cluster modularity. This method is not fully reproducible because it is sensitive to the initial coonditions (for instance, the order of the list of nodes that may be different on different computer/CPU/memory management). We observed that run on different computers, the clustering result may vary (different number of clusters and/or different cell to cluster assignation). To allow reproducibility of our results, we provide the cluster/cell assignation table in an output file (see below). 

Since we use these clusters to identify new series of cells that we remove from analysis in the next steps (see methods in the article), to reproduce the exact same results as the ones presented in our article, you have to use this clusters file as input in the next step.

#### Execution 
To run the step5, you need to execute two analysis. The first one produce the loom file containing the information about spliced and unspliced RNA. It is a python script provided by Velocyto. The second analysis takes this loom file as input and produces the RNA velocity analysis as an HTML report.

**_1. Launch Velocyto to produce the loom file of spliced/unspliced RNA counts_**

First ensure the script $WORKING_DIR/src/step5/execute_velocito.sh have execution rights by typing the following command:

    chmod 755 $WORKING_DIR/src/step5/execute_velocito.sh

Then launch the Velocyto tool inside the suitable Docker image with the command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity $WORKING_DIR/src/step5/execute_velocito.sh

Important note: This step is computationally intensive : the process will, at some steps, use all the available CPU and memory usage may exceed 30GB RAM.

**_2. Launch Velocyto result R analysis_**

Once the loom file has been produced by the previous command, ensure that the other input files are correctly copied into the $WORKING_DIR/data/step5/input folder (see above) and launch the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step5/JPGlab_LysoDC_QuinaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step5/output"), output_file = "JPGlab_LysoDC_QuinaryAnalysis.html", quiet=FALSE)'

#### Results:
Once the analysis done, here is the tree the $WORKING_DIR/data/step5 folder you may obtain (with the newly created "output" folder):

    step5
    ├── input
    │   ├── excluded_cells_contamination.txt    
    │   ├── excluded_cells_HighMitoGenePerc.txt 
    │   ├── excluded_cells_LowGeneNb.txt        
    │   ├── excluded_cells_LowUMINb.txt         
    │   └── excluded_cells_proliferation.txt    
    └── output
        ├── 10635173.loom                        # The loom file produced by the velocyto tool
        ├── cell_cluster_mapping.tsv             # The association of cells to clusters
        └── JPGlab_LysoDC_QuinaryAnalysis.html   # The analysis report

---

### Run step 6

#### Goal
This step aims to produce a more focus analysis of the RNA velocity. Cells analyzed in the previous step as contamination or not suitable for analysis thanks to the marker genes of the cluster they are part of are eliminated, providing a clearer view of the serached processes. In this step, we concentrate on the velocyto output analysis,
looking in details at some key genes, trying several combination of Velocyto parameters to make the mapping of relative velocyto of cells on the t-SNE embedding clearer.

#### Input
This step use input from previous analysis. Here is a tree presentation of the content of the folder $WORKING_DIR/data/step6 as it should be before the execution of the analysis.

    step6
    └── input
        ├── 10635173.loom                         # this file must be a copy from $WORKING_DIR/data/step5/output/10635173.loom
        ├── cell_cluster_mapping.tsv              # this file must be a copy from $WORKING_DIR/data/step5/output/cell_cluster_mapping.tsv
        ├── excluded_cells_contamination.txt      # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_contamination.txt
        ├── excluded_cells_HighMitoGenePerc.txt   # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_HighMitoGenePerc.txt
        ├── excluded_cells_LowGeneNb.txt          # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowGeneNb.txt
        ├── excluded_cells_LowUMINb.txt           # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowUMINb.txt
        └── excluded_cells_proliferation.txt      # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_proliferation.txt

#### Output
This step output a new cluster mapping file, done once cells have been filtered. It also ouput binary files (RDS format) containing the R objects produced by the velocyto analysis (spliced and unspliced matrix counts and relative velocity estimations) and the pagoda2 analysis.

We chose to export the complete R objects of this analysis in order to keep fixed some information that may vary from run to run (for instance t-SNE embedding). Those
objects will be resued in the next steps as input files.

**Important note:**  In this step the clustering is performed with the pagoda2 package and the _makeKnnGraph_ function (on PCA space) and the _getKnnClusters_ using the 'walktrap' method. This method using a random walk procedure can generate differences when re-running the clustering compare to our run. To allow reproducibility of our results, we provide the cluster/cell assignation table in an output file (see below). We use the original clustering to remove some cells and the clusters to remove are indicated in the Rmd code src/step6/JPGlab_LysoDC_SenaryAnalysis.Rmd at line 45. If you do not use the orginal clustering, you have to change the cluster number to suit your new clustering.

#### Execution
To execute the analysis, ensure input files are correctly copied into the $WORKING_DIR/data/step6/input folder (see above) and launch the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step6/JPGlab_LysoDC_SenaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step6/output"), output_file = "JPGlab_LysoDC_SenaryAnalysis.html", quiet=FALSE)'
    
**Important note:** This step is computationally intensive : the process will, at some steps, use all the available CPU and memory usage may exceed 30GB RAM.

#### Results
Once the analysis done, here is the tree the WORKING_DIR/data/step6 folder you may obtain (with the newly created "output" folder):

    step6
    ├── input
    │   ├── 10635173.loom
    │   ├── cell_cluster_mapping.tsv
    │   ├── excluded_cells_contamination.txt
    │   ├── excluded_cells_HighMitoGenePerc.txt
    │   ├── excluded_cells_LowGeneNb.txt
    │   ├── excluded_cells_LowUMINb.txt
    │   └── excluded_cells_proliferation.txt
    └── output
        ├── cell_cluster_mapping.tsv            # New cells clusters mapping (once the cells were filtered from undesired cells identified in step 5)
        ├── filtered_emat.rds                   # The Velocyto matrix of spliced RNA counts
        ├── filtered_nmat.rds                   # The Velocyto matrix of unspliced RNA counts
        ├── JPGlab_LysoDC_SenaryAnalysis.html   # The HTML report of the step 6 analysis
        ├── r_filtered.rds                      # The pagoda2 object with the counts, clusters and t-SNE informations
        └── rvel.cd.rds                         # the Velocyto object with the computed relative velocity estimation of cells

---

### Run step 7

#### Goal
This step aims at comparing the expression level of key genes and the flows of cells relative estimate velocities in order to decipher the differentiation steps.
This study allowed us to define sepcific zone in the t-SNE embedding, called _bassins_, corresponding to zones where the flow of velocities goes out (emission bassin) 
or goes in (attraction bassin (see article).

**Important note:** The definition of the bassins coordinates is done in the R code, more precisely in the src/step8/analysis_improve_rnavalocity_plot.R at line 114. If you do not use the original t-SNE embedding, you have to change those coordinates to suit the new embedding.

#### Input
This step use input from previous analysis. Here is a tree presentation of the content of the folder $WORKING_DIR/data/step7 as it should be before the execution of the analysis.

    step7
    └── input
        ├── filtered_emat.rds         # this file must be a copy from $WORKING_DIR/data/step6/output/filtered_emat.rds
        ├── filtered_nmat.rds         # this file must be a copy from $WORKING_DIR/data/step6/output/filtered_nmat.rds
        ├── r_filtered.rds            # this file must be a copy from $WORKING_DIR/data/step6/output/r_filtered.rds
        └── rvel.cd.rds               # this file must be a copy from $WORKING_DIR/data/step6/output/rvel.cd.rds

#### Output
This step produces files of the list of marker genes for each bassin and the list of cell barcode associated to each bassin. It also produces a report in HTML format.

#### Execution
To execute the analysis, ensure input files are correctly copied into the $WORKING_DIR/data/step7/input folder (see above) and launch the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step7/JPGlab_LysoDC_SeptenaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step7/output"), output_file = "JPGlab_LysoDC_SeptenaryAnalysis.html", quiet=FALSE)'

#### Results
Once the analysis done, here is the tree the WORKING_DIR/data/step7 folder you may obtain (with the newly created "output" folder):

    step7
    ├── input
    │   ├── filtered_emat.rds         
    │   ├── filtered_nmat.rds         
    │   ├── r_filtered.rds            
    │   └── rvel.cd.rds               
    └── output
        ├── basins_markers.csv                # Markers genes of bassins of attraction and emission
        ├── cells_in_basins.tsv               # List of cells barcodes in each bassin
        └── JPGlab_LysoDC_SeptenaryAnalysis   # the HTML report of the analysis

---

### Run step 8

#### Goal
This step produces a pseudotime analysis of the data using Monocle2. The study infer the cell type (SP, DP, TP, TN - see article) according to gene signatures and uses markers genes to infer the pseudo-time. Pseudot-time resultats are then compared to Velocyto results (e.g. content of bassins of attractions/emissions).

#### Input
This step use input from previous analysis. Here is a tree presentation of the content of the folder $WORKING_DIR/data/step8 as it should be before the execution of the analysis.

    step8
    └── input
        ├── filtered_raw_expression_matrix.csv            # this file must be a copy from $WORKING_DIR/data/step1/output/filtered_raw_expression_matrix.txt
        └── r_filtered.rds                                # this file must be a copy from $WORKING_DIR/data/step6/output/r_filtered.rds

#### Output
This step generate an HTML report of the analysis. It also generates the PDF of the figures of the BEAM analysis (Branched Expression Modeling Analysis) with three list
of genes: (i) a selection of marker genes used for the Monocle analysis (ii) the TOP50 markers genes and (iii) the top 50 best differentially expressed genes (DEG)

#### Execution
To execute the analysis, ensure input files are correctly copied into the $WORKING_DIR/data/step8/input folder (see above) and launch the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_pseudotime R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step8/JPGlab_LysoDC_OctonaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step8/output"), output_file = "JPGlab_LysoDC_OctonaryAnalysis.html", quiet=FALSE)'

#### Results
Once the analysis done, here is the tree the WORKING_DIR/data/step8 folder you may obtain (with the newly created "output" folder):

    step8
    ├── input
    │   ├── filtered_raw_expression_matrix.csv            
    │   └── r_filtered.rds                                
    └── output
        ├── BEAM_SelectedMarkers.pdf              # The BEAM result on selected marker genes
        ├── BEAM_Top50Markers.pdf                 # The BEAM result on Top 50 marker genes
        ├── BEAM_TOP50DEG.pdf                     # The BEAM result on TOP 50 DEG genes
        └── JPGlab_LysoDC_OctonaryAnalysis.html   # the HTML report of the analysis
