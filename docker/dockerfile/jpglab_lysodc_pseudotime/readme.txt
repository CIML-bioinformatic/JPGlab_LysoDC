This image contains:

 - R 3.4.2
 - Rstudio server (installation requires the userconf.sh file)
 - Packages for the Monocle analysis
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t jpglab_lysodc_pseudotime ~/workspace/ciml-docker/images/project/JPGlab/LysoDC/scrnaseq_pseudotime

# ######################
     RUN THE IMAGE
# ######################

docker run --name jpglab_lysodc_pseudotime -d -p 8787:8787 -v /mnt:/mnt -v /home/$USER/workspace:/home/$USER/workspace --network=host  -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) jpglab_lysodc_pseudotime
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <user>/rstudio
 
# ######################
	 NOTES
# ######################
 
 - To use knitr PDF compilation instead of Sweave, you have to go into Rstudio menu Tools->Global Options->Sweave->Weave Rnw files with.. and select "knitr".
 
