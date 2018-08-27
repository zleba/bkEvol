#ifndef _Settings_
#define _Settings_

#include <string>
#include <cassert>
#include <dlfcn.h>
#include <iostream>
#include <fstream>
#include <armadillo>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

struct Settings {

    double asMZ = 0.118;
    double freezingScale = 1;
    double LnFreeze2 = 2*log(freezingScale);
    double eps = 1e-7;
    int Nint = 257; // kT nodes in Nintegral
    int N = 257;// = 32*16 + 1; //must be 2*n+1
    int Nrap = 513;
    bool bkSolverGrid = false;
    bool toTrivial = true;

    double kT2Min = 1e-2, kT2Max = 1e6;
    double Lmin= log(kT2Min), Lmax = log(kT2Max);


    //const double Lmin= log(1e-4), Lmax = log(1e8);
    double mu2 = 1e-2;
    double xMin = 1e-6, xMax = 1;
    double rapMax = -log(xMin), rapMin = -log(xMax);
    //bool putZero = false;
    bool withMPI = false;

    string kernelType;

    string funStr;
    vector<tuple<double,double,double>> pars;
    int nPar;
    double (*fitFun)(double kT2, double x, const double *p);

    int maxIter;

    Settings() {}

    static Settings & I() {
        static Settings instance;
        return instance;
    }
    Settings(istream &Stream) {Init(Stream); }

    void Recalculate() {
        LnFreeze2 = 2*log(freezingScale);
        Lmin = log(kT2Min);
        Lmax = log(kT2Max);
        rapMax = -log(xMin);
        rapMin = -log(xMax);
    }

    template<typename T>
    void saveVar(string file, string name,  T var) {
        arma::Col<T> v(1);
        v(0) = var;
        v.save(arma::hdf5_name(file, name, arma::hdf5_opts::append));
    }

    template<typename T>
    void getVar(string file, string name,  T &var) {
        arma::Col<T> v(1);
        v.load(arma::hdf5_name(file, name));
        var = v(0);
    }


    void SaveToFile(string file)
    {
        saveVar(file, "alphaS",  asMZ);
        saveVar(file, "freezingScale",  freezingScale);
        saveVar(file, "eps",  eps);
        saveVar(file, "mu2",  mu2);

        saveVar(file, "NkT2",  N);
        saveVar(file, "NkT2int",  Nint);
        saveVar(file, "kT2Min",   kT2Min );
        saveVar(file, "kT2Max",   kT2Max );

        saveVar(file, "Nrap",  Nrap);
        saveVar(file, "xMin",  xMin);
        saveVar(file, "xMax",  xMax);

        saveVar(file, "bkSolverGrid",  (int) bkSolverGrid );
        saveVar(file, "toTrivial",  (int) toTrivial);
    }

    void LoadFromFile(string file)
    {
        getVar(file, "alphaS",  asMZ);
        getVar(file, "freezingScale",  freezingScale);
        getVar(file, "eps",  eps);
        getVar(file, "mu2",  mu2);

        getVar(file, "NkT2",  N);
        getVar(file, "NkT2int",  Nint);
        getVar(file, "kT2Min",   kT2Min );
        getVar(file, "kT2Max",   kT2Max );

        getVar(file, "Nrap",  Nrap);
        getVar(file, "xMin",  xMin);
        getVar(file, "xMax",  xMax);

        int bkSolverGridInt, toTrivialInt; 
        getVar(file, "bkSolverGrid",   bkSolverGridInt);
        getVar(file, "toTrivial",   toTrivialInt);
        bkSolverGrid = bkSolverGridInt;
        toTrivial    = toTrivialInt;

        Recalculate();
    }

    void printInfo()
    {
        cout << "Function " << funStr << endl;
        for(auto p : pars)
            cout <<"Par " <<  get<0>(p) <<" "<< get<1>(p) <<" "<< get<2>(p) << endl;
    }



    int getNpars(string funS)
    {
        bool isDone = false;
        int nPar = 0;
        for(int i = 0; i < 9; ++i) {
            if(funS.find("p["+to_string(i)+"]") != string::npos) {
                ++nPar;
                if(isDone) {
                    cout << "There is a gap between parameters" << endl;
                    assert(0);
                }
            }
            else {
                isDone = true;
            }
        }
        return nPar;
    }



    void Init(istream &Stream) {

        boost::property_tree::ptree tree;
        boost::property_tree::ini_parser::read_ini(Stream, tree);


        try {

            asMZ = tree.get<double>("Constants.alphaS");
            freezingScale = tree.get<double>("Constants.freezingScale");

            eps = tree.get<double>("Constants.eps");
            mu2 = tree.get<double>("Constants.mu2");
            //Rapidity properties
            Nrap = tree.get<int>("RapiditySpace.Nrap");
            xMin =  tree.get<double>("RapiditySpace.xMin");
            xMax =  tree.get<double>("RapiditySpace.xMax");

            //Running 
            kernelType  = tree.get<string>("RunningMode.kernelType");


            //Transverse properties
            N = tree.get<int>("TransverseSpace.NkT2");
            Nint = tree.get<int>("TransverseSpace.NkT2int");
            kT2Min = tree.get<double>("TransverseSpace.kT2Min");
            kT2Max = tree.get<double>("TransverseSpace.kT2Max");
            bkSolverGrid = tree.get<bool>("TransverseSpace.bkSolverGrid");
            toTrivial = tree.get<bool>("TransverseSpace.toTrivial");

            Recalculate();
            
            //Fit Properties
            maxIter = tree.get<int>("Fit.maxIter");


            funStr = tree.get<string>("Fit.function");
            nPar = getNpars(funStr);
            for(int i = 0; i < nPar; ++i) {
                string par = tree.get<string>("Fit.p"+to_string(i));
                istringstream iss(par);
                vector<double> parsNow;
                double pNow;
                while(iss >> pNow) parsNow.push_back(pNow);
                assert(parsNow.size() == 1 || parsNow.size() == 3);

                double p, pmin, pmax;
                p = parsNow[0];
                pmin = pmax = 0;
                if(parsNow.size() == 3) {
                    pmin = parsNow[1];
                    pmax = parsNow[2];
                }
                pars.push_back(make_tuple(p, pmin, pmax));

                cout << "Reading parameter "<< i <<" : " << p << " "<< pmin << " "<< pmax << endl;
                //if(iss.good()) cout << "String is good " << p <<" "<< pmin <<" "<< pmax << endl;
            }
            //exit(0);

        }
        catch(const std::exception& e) {
            cout << "Some of parameters in steering not defined:" << endl;
            cout << e.what() << endl;
            exit(1);
        }

        assert((N - 1) %2 == 0);
        assert((Nint - 1) %2 == 0);


        cout << "Used evolution parameters" << endl;
        cout << "Rapidity Space" << endl;
        cout << "xMin = " << xMin << endl;
        cout << "xMax = " << xMax << endl;

        cout << "Nrap = " << Nrap << endl;
        cout << "Nkt = " << N << endl;
        cout << "bkSolverGrid = " << bkSolverGrid << endl;

        //exit(1);

        if(bkSolverGrid) assert(N == Nint);

        cout << "Compiling function: " << endl;
        cout <<  funStr << endl;
        Compile(funStr);

        //cout << "exiting " << endl;
        //exit(0);
        //cout << "f(1,0.01) = " << fitFun(1, 0.01) << endl;
        //cout << "Done " << endl;
        //exit(0);

    }

    void WriteFunction (string str)
    {
        ofstream f("fitFun.cpp");
        f << "#include <math.h>" << endl;
        f << "extern \"C\" {" << endl;
        f << "double fitFun(double kT2, double x, const double *p) {" << endl;
        f << "    return (" << str <<");"<< endl;
        f << "}" << endl;
        f << "}" << endl;
        //return 
    }
                                    

    void Compile (string str)
    {
        //double (*fun)(double x);

        string f = __FILE__;
        f.substr(0, f.size() -  10);//inc/Settings.h 

        #define STRINGIFY(x) #x
        #define TOSTRING(x) STRINGIFY(x)

        string n = TOSTRING(pwdDir);
        n.substr(0, f.size() - 14);
        n += "/obj/fitFun.so";
        cout <<"RADEK " <<  f <<" "<< n << endl;


        
        WriteFunction(str);
        int st = system(("g++ fitFun.cpp -o " + n + " -shared -fPIC").c_str());
        assert(st >= 0);
        //system("sleep 5");
        void *lib = dlopen(n.c_str(), RTLD_LAZY);
        assert(lib);

        fitFun = reinterpret_cast<double (*)(double,double,const double*)>(dlsym(lib, "fitFun"));
        //cout << dlerror() << endl;
        assert(dlerror() == NULL);
        assert(fitFun);

        //cout << "Fun " <<  fitFun(5, 0.5) << endl;

        return;
    }
                        

};

#endif
