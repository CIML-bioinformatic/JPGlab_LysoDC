#LysoDC


THIS WORK IS IN PROGRESS


This github project contains all you need to reproduce the analysis reported in the article (and more). The project contains data, R and shell scripts and docker images.
Source code (script and dockerfile) are available in the github repository. Data and builded Docker images are available on download. Intructions are provided below

### Prepare the environments

In order to prepare the environment for analysis execution you need:
- Clone the github repository
- Download the docker images tar file
- Load the docker images on your system
- Download the raw data (sequencing data)
- Set the WORKING_DIR environment variable

##### Download the docker images

Docker images tar file are stored on Zenodo. Open a shell command and change dir to the root of the cloned Git repository. Then execute the following commands to download the images tar files to the right project folder:

    wget URL_jpglab_lysodc_seurat.tar -o docker/image/jpglab_lysodc_seurat/jpglab_lysodc_seurat.tar
    wget jpglab_lysodc_rnavelocity.tar -o docker/image/jpglab_lysodc_rnavelocity/jpglab_lysodc_rnavelocity.tar
    wget jpglab_lysodc_pseudotime.tar -o docker/image/jpglab_lysodc_pseudotime/jpglab_lysodc_pseudotime.tar

##### Load docker images

In order to execute analysis, you must load the provided docker images onto your Docker. Docker must be installed on your system. 
See https://docs.docker.com/install/ for details on Docker installation.
Open a shell command, change dir to the folder in which you cloned the project and type:

    docker load -i docker/image/jpglab_lysodc_seurat/jpglab_lysodc_seurat.tar
    docker load -i docker/image/jpglab_lysodc_rnavelocity/jpglab_lysodc_rnavelocity.tar
    docker load -i docker/image/jpglab_lysodc_pseudotime/jpglab_lysodc_pseudotime.tar

Those commands may take some time. If you encounter an issue loading some docker image layer, try again. Sometimes issue would be resolved. 

##### Download the raw data

Raw data are availbale for download on GEO. They are also available in a more convevient format on Zenodo. To dowload tem, use the following command:



##### Set the WORKING_DIR variable
 
You must set an environment variable called WORKING_DIR with a value set to the path to the folder in which you cloned git repository

**On linux:**

    export WORKING_DIR=/home/spinellil/workspace/github-CIML-LysoDC

### Run the analysis

Analysis can be directly run inside a docker container by compiling a rmarkdown file. The rmarkdown will launch the required analysis for the step
and produce a final HTML report. This report will be generated in the <WORKING_DIR>/data/step<X>/output where <WORKING_DIR> is the folder where you clone the git repository (and set the WORKING_DIR environment variable) and <X> is the number of the analysis step you are executing (1 or 5 to 8).

##### Run step 1

To run the step1, ensure you have correctly downloaded the data in the folder <WORKING_DIR>/data/raw and run the following command:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_seurat R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step1/JPGlab_LysoDC_PrimaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step1/output"), output_file = "JPGlab_LysoDC_PrimaryAnalysis.html", quiet=FALSE)'


##### Run step 5

To run the step5, you need to execute two analysis. The first one produce the loom file containing the information about spliced and unspliced RNA. It is a python script provided by Velocyto. The second analysis takes this loom file as input and produces the RNA velocity analysis as an HTML report.

**Launch Velocyto to produce the loom file

* ensure  $WORKING_DIR/src/step5/execute_velocito.sh have execution rights by typing the following command:

    chmod 755 $WORKING_DIR/src/step5/execute_velocito.sh

* launch the docker image with the command to execute:

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity $WORKING_DIR/src/step5/execute_velocito.sh

##### Launch Velocyto result analysis

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step5/JPGlab_LysoDC_QuinaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step5/output"), output_file = "JPGlab_LysoDC_QuinaryAnalysis.html", quiet=FALSE)'

####Run step 6

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR jpglab_lysodc_rnavelocity R -e 'WORKING_DIR=Sys.getenv( "WORKING_DIR");rmarkdown::render( input=file.path( WORKING_DIR, "src/step6/JPGlab_LysoDC_SenaryAnalysis.Rmd"), output_dir = file.path( WORKING_DIR, "data/step6/output"), output_file = "JPGlab_LysoDC_SenaryAnalysis.html", quiet=FALSE)'
