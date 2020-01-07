This image contains:

 - R 3.4.4
 - Rstudio server (installation requires the userconf.sh file)
 - Python 3.6.6
 - Volocito (http://velocyto.org/velocyto.py/install/index.html)
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t jpglab_lysodc_rnavelocity ~/workspace/ciml-docker/images/project/JPGlab/LysoDC/scrnaseq_rnavelocity

# ######################
     RUN THE IMAGE
# ######################

docker run --name jpglab_lysodc_rnavelocity -d -p 8787:8787 -v /mnt:/mnt -v /home/$USER/workspace:/home/$USER/workspace --network=host -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  jpglab_lysodc_rnavelocity
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <user>/rstudio
 
# ######################
	 NOTES
# ######################
 
 
 
