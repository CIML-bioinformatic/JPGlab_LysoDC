# LysoDC


THIS WORK IS IN PROGRESS

This github project contains the instructions and material to reproduce the analysis reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Requested data and builded Docker images are available on download. Intructions to reproduce the analysis are provided below.

To reproduce the analysis, we have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

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


## Run the analysis

Analysis can be directly run inside docker containers by compiling Rmarkdown files. The Rmarkdown file comilation will launch the required analysis for the step
and produce a final HTML report. 

Each step report will be generated in the <WORKING_DIR>/data/step<X>/output where <WORKING_DIR> is the folder where you clone the git repository (and set the WORKING_DIR environment variable) and <X> is the number of the analysis step you are executing (1 or 5 to 8).

### Run step 1

**__Goal:__** This step is a standard quality control and first analysis of the cell heterogeneity. 

**__Output:__** This step output files that contain cellular barcodes associated to cell do not passing
one of the QC test. It also output two matrices of gene expression o,e with the raw UMI counts, and one with the normalized UMI counts. 
Finally, it output the HTML report of the analysis.

**__Execution:__** To run the step1, ensure you have correctly downloaded the data in the folder <WORKING_DIR>/data/raw and run the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_seurat R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step1/JPGlab_LysoDC_PrimaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step1/output"), output_file = "JPGlab_LysoDC_PrimaryAnalysis.html", quiet=FALSE)'

**__Results:__** One the analysis done, here is the tree the $WORKING_DIT/data/step1 folder you may obtain:

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


### Run step 5

**__Goal:__** This step aims to produce a first version of the RNA velocity analysis of the data. Thi sstep is also used to identify cell considered as contamination or "in dying" state. These cells are removed from the analysis to improve the quality of the next steps.

**__Input:__** This step use input from previous analysis. Here is a tree présentation of the content of the folder $WORKING_DIR/data/step5 as it should be before the execution
of the analysis.

    step5
    └── input
        ├── excluded_cells_contamination.txt    #(this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_contamination.txt)
        ├── excluded_cells_HighMitoGenePerc.txt #(this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_HighMitoGenePerc.txt)
        ├── excluded_cells_LowGeneNb.txt        #(this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowGeneNb.txt)
        ├── excluded_cells_LowUMINb.txt         #(this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowUMINb.txt)
        └── excluded_cells_proliferation.txt    #(this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_proliferation.txt)

**__Execution:__** To run the step5, you need to execute two analysis. The first one produce the loom file containing the information about spliced and unspliced RNA. It is a python script provided by Velocyto. The second analysis takes this loom file as input and produces the RNA velocity analysis as an HTML report.

#### 1. Launch Velocyto to produce the loom file

First ensure the script $WORKING_DIR/src/step5/execute_velocito.sh have execution rights by typing the following command:

    chmod 755 $WORKING_DIR/src/step5/execute_velocito.sh

Then launch the Velocyto tool inside the suitable Docker image with the command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity $WORKING_DIR/src/step5/execute_velocito.sh

**Important note:** This step is computationally intensive : the process will, at some steps, use all the available CPU and memory usage may exceed 30GB RAM.

#### 2. Launch Velocyto result analysis

Once the loom file has been produced by the previous command, launch the following command to produce the analysis report on RNA velocity:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step5/JPGlab_LysoDC_QuinaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step5/output"), output_file = "JPGlab_LysoDC_QuinaryAnalysis.html", quiet=FALSE)'

**__Output:__** One the analysis done, here is the tree the $WORKING_DIT/data/step5 folder you may obtain:

    step5
    ├── input
    |   ├── excluded_cells_contamination.txt    
    |   ├── excluded_cells_HighMitoGenePerc.txt 
    |   ├── excluded_cells_LowGeneNb.txt        
    |   ├── excluded_cells_LowUMINb.txt         
    |   └── excluded_cells_proliferation.txt    
    └── output
        ├── 10635173.loom                        # The loom file produced by the velocyto tool
        ├── cell_cluster_mapping.tsv             # The association of cells to clusters
        └── JPGlab_LysoDC_QuinaryAnalysis.html   # The analysis report


#### Important note

The clustering in this step is performed thanks to the package Pagoda2, using the _makeKnnGraph_ function (on PCA space) and the _getKnnClusters_ using the 'multilevel community' method. This method is based on a Louvain clustering with maximization of the cluster modularity. This method is not fully reproducible because it is sensitive to the initial coonditions (for instance, the order of the list of nodes that may be different on different computer/CPU/memory management). We observed that run on different machine, the clustering result may vary. To allow reproducibility of our results, we provide the cluster/cell assignation table in an output file (see above). 
Since we use this cluster to identify new series of cells that we remove from analysis in the next step (see methods in the article), to reproduce the exact same result as the one presented in our article, you have to use this cluster file as input in the next step.

### Run step 6

**Goal:** This step aims to produce a more focus analysis of the RNA velocity. Cells analyzed in the previous step as contamination or not suitale for analysis thanks to the marker genes of the cluster they are part of are eliminated, providing a clearer view of the serached processes.

**__Input:__**

    step6
    ├── input
        ├── 10635173.loom                         # this file must be a copy from $WORKING_DIR/data/step5/output/10635173.loom
        ├── cell_cluster_mapping.tsv              # this file must be a copy from $WORKING_DIR/data/step5/output/10635173.loom
        ├── excluded_cells_contamination.txt      # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_contamination.txt
        ├── excluded_cells_HighMitoGenePerc.txt   # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_HighMitoGenePerc.txt
        ├── excluded_cells_LowGeneNb.txt          # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowGeneNb.txt
        ├── excluded_cells_LowUMINb.txt           # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_LowUMINb.txt
        └── excluded_cells_proliferation.txt      # this file must be a copy from $WORKING_DIR/data/step1/output/excluded_cells_proliferation.txt



    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step6/JPGlab_LysoDC_SenaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step6/output"), output_file = "JPGlab_LysoDC_SenaryAnalysis.html", quiet=FALSE)'

### Run step 7

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step7/JPGlab_LysoDC_SeptenaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step7/output"), output_file = "JPGlab_LysoDC_SeptenaryAnalysis.html", quiet=FALSE)'

### Run step 8

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step8/JPGlab_LysoDC_OctonaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step8/output"), output_file = "JPGlab_LysoDC_OctonaryAnalysis.html", quiet=FALSE)'



