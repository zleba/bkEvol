pwd=$PWD
mkdir -p $pwd/arma/install
cd $pwd/arma
wget  https://sourceforge.net/projects/arma/files/armadillo-8.500.1.tar.xz/download  && tar xf download
cd $pwd/arma/armadillo*

cmake . -DCMAKE_INSTALL_PREFIX:PATH=$pwd/arma/install
make
make install
