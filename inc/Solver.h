#ifndef Solver_H_
#define Solver_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <string>

#include <armadillo>
#include <mpi.h>


#if hasGPU == true
    #include "gpuBooster.h"
#endif

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "alphaSpline.h"

#include "Settings.h"
#include "utils.h"

using namespace std;

vector<double> GetWeights(int Size);
pair<arma::mat, arma::mat> GetTransMatrices(int base, int ext, bool toTrivial);


inline pair<int,int> GetRankSize()
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return make_pair(world_rank, world_size);
}


inline pair<long long,long long> GetStartEnd(bool isMPI, long long Min, long long Max)
{
    if(!isMPI) return make_pair(Min, Max);

    int rank, nrank;
    tie(rank,nrank) = GetRankSize();
    return GetStartEnd(nrank, rank, Min, Max);
}



struct Nodes {
    Nodes(int _N, double _a, double _b) : N(_N), a(_a), b(_b) {}
    Nodes() {}
    void Init(int _N, double _a, double _b) { N=_N; a=_a; b=_b; }
    int N;
    double a, b;
    vector<double> xi, SqrtExpXi,  wi;
    void CalcNodes(bool bkSolverMode = true) {
        /*
        wi.resize(N);
        for(int i = 0; i < N; ++i) {
            double Cos = cos((2*i+1)/(2.*N) * M_PI);
            xi[i] = (a + b +  Cos * (a-b) )/2.;
            //cout << i<<" "<< xi[i] << endl;
            wi[i] = M_PI / N * sqrt(1-Cos*Cos) * (b-a)/2.;

            //double Cosp = cos((i+1.)/(N+1) * M_PI);
            //xi[i] = (a + b + Cosp * (a-b) )/2.;
            //wi[i] = M_PI / (N+1) * (1-Cosp*Cosp) / sqrt(1-Cosp*Cosp);
            //cout << (1-Cosp*Cosp) << endl;
            //assert(isfinite(wi[i]));
        }
        */

        //Classical nodes ala Clenshaw–Curtis
        if(!bkSolverMode) {
            xi.resize(N);
            for(int i = 0; i < N; ++i) {
              double Cos = cos(i /(N-1.) * M_PI);
              xi[i] = (a + b +  Cos * (a-b) )/2.;
              //cout << "R " << exp(0.5*xi[i]) << endl;
            }
            wi = GetWeights(N);
            for(auto &w : wi) w *= (b-a)/2.;

        }
        //Nodes used in BKsolver
        else {
            xi.resize(2*N-1);
            for(int i = 0; i < N*2-1; ++i) {
                double Cos = cos(i /(2*N-2.) * M_PI);
                double an = a - (b-a);
                double bn = b;
                xi[i] = (an + bn +  Cos * (an-bn) )/2.;
                //cout << "R " << exp(0.5*xi[i]) << endl;
            }

            xi.erase(xi.begin(), xi.begin() + N - 1);
            //cout << "Begin " << endl; 
            //for(auto x : xi)
                //cout << "ahoj " <<x <<" "<< exp(0.5*x) << endl;
            //cout << "End " << endl; 

            wi = GetWeights(2*N-1);
            wi.erase(wi.begin(), wi.begin() + N - 1);
            for(auto &w : wi) w *= 2*(b-a)/2.;

            wi.front() *= 0.5;
        }
        cout << "I am here " <<__LINE__<< endl;


        SqrtExpXi.resize(xi.size());
        for(int i = 0; i < xi.size(); ++i)
            SqrtExpXi[i] = exp(0.5*xi[i]);

        //exit(0);
    }
    double Integrate(function<double(double)> fun) {
        double sum = 0;
        for(int i = 0; i < N; ++i) {
           sum += fun(xi[i]) * wi[i];
        }
        return sum;
    }
    double Integrate(vector<double> vals) {
        assert(vals.size() == wi.size());
        reverse(vals.begin(), vals.end());
        double sum = 0;
        for(int i = 0; i < N; ++i) {
           sum += vals[i] * wi[i];
        }
        return sum;
    }



};



struct Solver {
    Settings S;
    bool withMPI = false;


    Nodes nod, nodBase;

    double alphaS(double l, double lp);
    
    /*
    Solver(int N_) : Nint(1*(N_-1)+1), N(N_),
        nod(Nint, Lmin, Lmax), nodBase(N, Lmin, Lmax) {
        nod.CalcNodes(false);
        nodBase.CalcNodes(false);
        tie(redMat,extMat) = GetTransMatrices(N, Nint, toTrivial);
    }
    */

    Solver(Settings SS) {
        SS.Recalculate();
        S = SS;

        nod.Init(S.Nint, S.Lmin, S.Lmax);
        nodBase.Init(S.N, S.Lmin, S.Lmax);
        nod.CalcNodes(S.bkSolverGrid);
        nodBase.CalcNodes(S.bkSolverGrid);
        tie(redMat,extMat) = GetTransMatrices(S.N, S.Nint, S.toTrivial);

        alphaSpline::FixMasses( 1e-8, 4.2,	1e21);
        alphaSpline::FixParameters(2, S.asMZ, 5, 91.2);
        //cout << "Radek " << alphaSpline::alphaS(2*log(91.2))<< endl;

        InitF([](double x, double kT2) {
            return kT2*exp(-kT2);
        });
    }
    Solver(std::istream &Stream) : Solver(Settings(Stream)) {}



    //vector<arma::mat> matN, matNDiag;
    arma::cube matN, matNDiag, matNInv;
    vector<arma::vec> PhiRapN, Phi0N;

    arma::cube convF2, convFL;
    vector<arma::vec> F2rap, FLrap;
    
    //vector<arma::vec> &GetF2() {return F2rap;}
    //vector<arma::vec> &GetFL() {return FLrap;}


    arma::mat extMat, redMat;

#if hasGPU == true
    gpuBooster gpu;
#endif


    void CalcEvolKernel();
    void SetSolution(function<double(double, double)> fun);



    void SaveEvolKernels(string file) {
        //string aStag = to_string(1000*asMZ);
        //file += "_as"+ to_string(lrint(1000*asMZ));
        cout << "Saving Evol Kernels to " << file << endl;
        //matN.save(file+"/kernel_base.h5", arma::hdf5_binary);
        //matNDiag.save(file+"/kernel_diag.h5", arma::hdf5_binary);
        //matNInv.save(file+"/kernel_inv.h5", arma::hdf5_binary);

        matN.save(arma::hdf5_name(file, "kernel"));
        matNDiag.save(arma::hdf5_name(file, "diag", arma::hdf5_opts::append));
        matNInv.save(arma::hdf5_name(file, "inv", arma::hdf5_opts::append));
        S.SaveToFile(file);
    }


    void LoadEvolKernels(string file) {
        //matN.load(file+"/kernel_base.h5", arma::hdf5_binary);
        //matNDiag.load(file+"/kernel_diag.h5", arma::hdf5_binary);
        //cout << "RADEK size " << matNDiag.slice(0).n_rows << endl;
        //matNInv.load(file+"/kernel_inv.h5", arma::hdf5_binary);
        matN.load(arma::hdf5_name(file, "kernel"));
        matNDiag.load(arma::hdf5_name(file, "diag"));
        matNInv.load(arma::hdf5_name(file, "inv"));
        cout << "Load sizes " << matN.n_slices <<" "<< matNDiag.n_slices << endl;
        if(matN.n_slices == 0) {
            cout << "Problems with reading EvolKernels from file " << file << endl;
            assert(0);
        }
        S.LoadFromFile(file);

    }

    void LoadConvKernels(string file) {
        assert(convF2.load(file+"/conv_F2.h5", arma::hdf5_binary));
        assert(convFL.load(file+"/conv_FL.h5", arma::hdf5_binary));

        //MPI_Bcast(convF2.memptr(), convF2.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(convFL.memptr(), convF2.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }




    arma::vec GetLinSolution(const arma::mat &Mat, const arma::vec &y) {
       return arma::solve(Mat, y); 
    }

    vector<double> GetRegSolution(const vector<vector<double>> &MatEq, const vector<double> &y, const vector<double> &yReg);
    arma::vec IterSolution(const arma::mat &Mat, double factor, const arma::vec &y);


    void EvolveAll();
    void CalcF2L();
    void DoIteration();
    void RunIterations(int Niter, bool init = true);
    arma::vec GetRHS(const arma::vec &PHI);
    void Step(double delta);


    //Function of x and kT2
    void InitF(function<double(double, double)> fun);

    double Interpolate(double y, double L);
    void PrintBaseGrid();
    void PrintGrid();
    void PrintReduce();

    static arma::mat vector2matrix(vector<arma::vec> &Vec) {
        arma::mat matPhi(Vec.size(), Vec[0].n_rows);
        for(int i = 0; i < Vec.size(); ++i)
            matPhi.row(i) =  Vec[i].t();
        return matPhi;
    }

};



#endif 
