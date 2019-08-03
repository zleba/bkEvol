#ifndef _nodes_
#define _nodes_

#include <armadillo>

//For points between 0 and 1
arma::vec GetWeights(int Size);
arma::vec GetNodes(int Size);
arma::mat GetCoefs(int oldSize, bool isInverse = false);
arma::mat GetCoefs(int Size, double a, double b);
arma::mat GetCoefsGeneric(int Size, int SizeI, double a, double b);
arma::mat GetCoefsCheb(int oldSize);

#endif
