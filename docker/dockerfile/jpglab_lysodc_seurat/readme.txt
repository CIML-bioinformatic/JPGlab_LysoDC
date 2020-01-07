This image contains:

 - R 3.4.4
 - Rstudio server (installation requires the userconf.sh file)
 - Packages for the Seurat analysis
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t jpglab_lysodc_seurat $WORKING_DIR/docker/jpglab_lysodc_seurat

# ######################
     RUN THE IMAGE
# ######################

 docker run --name jpglab_lysodc_seurat -d -p 8787:8787 -v $WORKING_DIR:$WORKING_DIR -e PASSWORD=LysoDC -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) jpglab_lysodc_seurat
 
 Note: you can change the value after PASSWORD= to set your own password
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <you_user_name>/LysoDC (or your own password if you changed it).
 You will obtain a RStudio interface you can use to browse the code and execute step by step analysis.
 Don't forget to set the variable WORKING_DIR in R with the value you set it in environment variable $WORKING_DIR

 
