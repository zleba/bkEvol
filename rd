#/bin/sh
image=bkevol
sudo docker run  -e ARMADIR='/arma'  -v $PWD:/bkevol --rm -it  $image   /bin/bash -c "cd /bkevol/; $*"
#sudo docker run -e  DISPLAY=$DISPLAY  -e ARMADIR='/'  -v /tmp/.X11-unix:/tmp/.X11-unix -v $PWD:/bkevol --rm -it --user $(id -u) $image   /bin/bash -c "cd /bkevol/; $*"
