pwd=$PWD
mkdir -p $pwd/arma/install
cd $pwd/arma
#wget  https://sourceforge.net/projects/arma/files/armadillo-8.500.1.tar.xz/download  && tar xf download
wget  http://sourceforge.net/projects/arma/files/armadillo-9.100.5.tar.xz  && tar xf *.tar.xz
cd $pwd/arma/armadillo-9.100.5

cmake . -DCMAKE_INSTALL_PREFIX:PATH=$pwd/arma/install
make
make install
