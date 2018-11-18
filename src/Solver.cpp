#include "Solver.h"
#include "kernels.h"

//typedef double (*FUNker)(double, double, double);

typedef double (Kernel::*FUNker)(double,double,double);

double dummy(double,double,double) {return 0;}


struct KernelPtr {
    FUNker Off, Diag;

    KernelPtr(FUNker O, FUNker D) {
        Off = O;
        Diag = D;
    }
    KernelPtr() {}
};





//--------------------------------------------------------------------
void Solver::CalcEvolKernel()
{
    /*
    int provided;
    int argc = 1;
    char **argv = new char*[2];
    argv[0] = "bkEvol";
    argv[1] = "bkEvol";
    cout << "I am Before" << endl;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    cout << "I am after" << endl;
    assert(MPI_THREAD_FUNNELED == provided);
    */


    //matN.resize(Nrap, arma::mat(N,N,arma::fill::zeros));

    alphaSpline::FixMasses( 1e-8, 4.2,	1e21);
    alphaSpline::FixParameters(2, S.asMZ, 5, 91.2);

    Kernel ker;
    ker.eps = S.eps;
    ker.mu2 = S.mu2;
    ker.rapMin = S.rapMin;
    ker.rapMax = S.rapMax;
    ker.Nrap = S.Nrap;
    ker.LnFreeze2 = S.LnFreeze2;
    ker.SqrtExpXi = nod.SqrtExpXi;
    ker.myCheck = arma::mat(ker.SqrtExpXi.size(), ker.SqrtExpXi.size(), arma::fill::zeros);

    map<string,KernelPtr> kerMap;

#define initOne(name,tag)  kerMap[STRINGIFY(name##__##tag)] = KernelPtr(&Kernel::name##__##tag##_Off, &Kernel::name##__##tag##_Diag);
#define initTwo(name)      initOne(name,Eps); initOne(name,Sub);
                                                           

    initOne(BFKLplain,Sub);
    initTwo(BFKL);
    initTwo(BFKL_res);
    initTwo(BFKL_res_kc_simp);
    initTwo(BFKL_res_kc_v_r_simp);
    initTwo(BFKL_res_kc_full);
    initTwo(BFKL_res_kc_v_r_full);
    initOne(BFKL_res_DGLAP,Eps);
    initOne(BFKL_res_DGLAP,ZEps);
    initOne(BFKL_res_kc_full_DGLAP,Eps);
    initOne(BFKL_res_kc_full_DGLAP_simp_kc,Eps);
    initOne(BFKL_res_kc_full_DGLAP_simp_kc,ZEps);
    initOne(BFKL_res_kc_full_DGLAP_full_kc,Eps);


    cout << "Available kernels:" << endl;
    for(auto &el : kerMap) {
        cout << el.first << endl;
    }

    string name =  S.kernelType;

    string base      = S.kernelType.substr(0, S.kernelType.find(':'));
    string tag       = S.kernelType.substr(S.kernelType.find(':')+1, string::npos);
    string wholeName = base+"__"+tag;
    if(!kerMap.count(wholeName)) {
        cout << "Kernel "<<base<< " with tag "+tag+ " is not implemented\nExiting" << endl;
        exit(1);
    }
    cout << "Using kernel " << base <<" with tag "<< tag<< endl;

    //FUNker KernelOff, KernelDiag;

    FUNker KernelOff  = kerMap[wholeName].Off;
    FUNker KernelDiag = kerMap[wholeName].Diag;


    matN.zeros( S.N, S.N, S.Nrap);
    matNInv.zeros( S.N, S.N, S.Nrap);

    int fac = (S.Nint-1)/(S.N-1);

    int start, end;
    tie(start,end) = GetStartEnd(withMPI, 0, S.Nrap-1);




    cout << "Start+end|nrap " << start <<" "<< end <<"|" << S.Nrap<< endl;

    double stepY = (S.rapMax - S.rapMin) / (S.Nrap-1);
    //for(int y = 0; y < Nrap; ++y) { 
    for(int y = start; y <= end; ++y) { 
        double rap = S.rapMin + stepY * y;
        double z = exp(-rap);

        arma::mat mTemp(S.Nint,S.Nint,arma::fill::zeros);
        cout << "\rMatrix init y " << y << flush;// << endl;

#pragma omp parallel for
        for(int i = 0; i < S.Nint; ++i)   //loop over L
            for(int j = 0; j < S.Nint; ++j) { //loop over L' (integral)
                double L  = nod.xi[i];
                double Lp = nod.xi[j];
                double w = nod.wi[j];

                //double l  = exp(0.5*L);
                //double lp = exp(0.5*Lp);
                double l  = nod.SqrtExpXi[i];
                double lp = nod.SqrtExpXi[j];

                //mat[y][i][j] += Kernel85(l, lp, z) * w;
                //mat[y][i][i] += Kernel85Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub81(l, lp, z) * w;
                //mat[y][i][i] += KernelSub81Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub88(l, lp, z) * w;

                //matDiag[y][i][j] = Kernel85zDiag(l, lp, z) * w; //Diag-z DGLAP part

                if(!S.toTrivial || i % fac == 0) {

                    //mTemp(i,j) += KernelBFKL(l, lp, z) * w;
                    //mTemp(i,i) += KernelBFKLDiag(l, lp, z) * w;

                    mTemp(i,j)     += (ker.*KernelOff)(l, lp, z) * w;
                    mTemp(i,i)     += (ker.*KernelDiag)(l, lp, z) * w;
                
                    if(isnan(mTemp(i,j)))
                        cout <<"Hela "<< (ker.*KernelOff)(l, lp, z) <<" "<< l <<" "<<lp<<" "<< z << endl;

                    //cout << "Hela " << mTemp(i,j) << " "<<mTemp(i,i)<<" "<<  mDiagTemp(i,j) << endl;

                    //mTemp(i,j) += Kernel86(l, lp, z) * w;
                    //mTemp(i,i) += Kernel86Diag(l, lp, z) * w;

                    //mTemp(i,j) += Kernel9(l, lp, z) * w;
                    //mTemp(i,i) += Kernel9Diag(l, lp, z) * w;

                }



                //Clasicall BFKL
                /*
                   if(i != j)
                   mat[y][i][j] += Kernel(L, Lp) * w;

                   double kerDs, kerDr;
                   tie(kerDs, kerDr) = KernelDiag(L, Lp);
                   if(i != j)
                   mat[y][i][i] += (kerDs+kerDr) * w;
                   else
                   mat[y][i][i] += (kerDr) * w;
                   */

            }
        //cout << redMat << endl;
        //cout << extMat << endl;


        //cout << extMat << endl;
        //cout << redMat << endl;
        //exit(0);

        matN.slice(y)     = stepY * redMat * mTemp * extMat;




        /*
           for(int i = 0; i < Nint; ++i)  //loop over L
           for(int j = 0; j < Nint; ++j) { //loop over L' (integral)
           double m  = matDiag[y][i][j];
           assert(isfinite(mat[y][i][j]));
           assert(isfinite(matDiag[y][i][j]));
        //cout << "Ihned "<<y<<" : " <<i<<" "<<j<<" "<< m <<" "<< mt << " "<< 2*(m-mt)/(m+mt) <<endl;
        //cout << "Mdiag "<<y<<" : " <<i<<" "<<j<<" "<< matDiag[y][i][j] <<endl;
        }
        */

        //if(y == 1) exit(0);

    }

    
    cout << "Check" << endl;
    cout << ker.myCheck << endl;


    cout << "Reduce start" << endl;
    if(withMPI) {
        //Merge things together
        MPI_Allreduce(MPI_IN_PLACE, matN.memptr(), matN.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }

    for(int y = start; y <= end; ++y) { 
        if(y == 0) {
            if(!ker.putZero)
                matNInv.slice(y) = inv(arma::mat(S.N,S.N,arma::fill::eye) );
        }
        else
            matNInv.slice(y) = inv(arma::mat(S.N,S.N,arma::fill::eye) - 0.5*matN.slice(0) );
    }
    if(withMPI) MPI_Allreduce(MPI_IN_PLACE, matNInv.memptr(), matNInv.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    cout << "Reduce done" << endl;


    //putZero = ker.putZero;


    //cout << matN << endl;
    //cout << matNInv << endl;



    //matN.save("ahoj.hdf5", arma::hdf5_binary);
    //MPI_Finalize();
    //exit(0);


    //for(int y = 0; y < Nrap; ++y) { 
    //MPI_Allreduce(MPI_IN_PLACE, matN.slice(y).memptr(),N*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //}

    /*
       for(int y = start; y <=end; ++y) {
       int rank = GetRankSize().first;
       MPI_Bcast(matN[y].memptr(), N*N, MPI_DOUBLE, rank, MPI_COMM_WORLD);
       }
       */

    //exit(0);

    /*
       if(start == 0) {
       for(int y = 0; y < Nrap; ++y) {
       const double stepY = (rapMax - rapMin) / (Nrap-1);
       for(int i = 0; i < N; ++i)
       for(int j = 0; j < N; ++j) {
       cout <<"Kernel "<<y<<" : "<< i <<" "<< j <<" "<<setprecision(10)<< matN[y](i,j) <<" "<<  stepY<< endl;
    //cout <<"Kernel "<<y<<" : "<< i <<" "<< j <<" "<<setprecision(10)<< matMPI[y](i,j) <<" "<< matMPIDiag[y](i,j)<<" "<< stepY<< endl;
    }
    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
    */




    /*
       exit(0);
       */
}






//--------------------------------------------------------------------
void Solver::EvolveAll()
{
    //const double stepY = (rapMax - rapMin) / (Nrap-1);

    cout << "Sizes " << matN.n_slices << endl;
    int Nrap = matN.n_slices;
    int N = matN.slice(0).n_rows;

    F2rap.resize(Nrap);
    FLrap.resize(Nrap);

    //In case of convF2, convFL not filled
    if(convF2.n_rows == 0) {
        cout << "Dummy init of the convolution cube" << endl;
        convF2.zeros(46, N,  Nrap);
        convFL.zeros(46, N,  Nrap);
    }


    //bool doGPU = true;
#if hasGPU == true
    if(!gpu.isInited) gpu.InitAll(matN, convF2, convFL);
    gpu.ResetVector();
#endif


    //cout << Phi0N[8] << endl;
    //return;

    int start = 0;

    //Classical approach
    if(start == 0) {
        cout << "Size of matrixes " << N <<" : "<<   endl;
        /*
        arma::mat MatEq = arma::mat(N,N,arma::fill::eye) -  matNDiag.slice(0);
        cout << "Is put Zero " << putZero << endl;
        if(putZero) {
            PhiRapN[0] = arma::vec(N, arma::fill::zeros);
        }
        else {
            PhiRapN[0] = GetLinSolution(MatEq, Phi0N[0]);
        }
        */
        PhiRapN[0] =  matNInv.slice(0) * Phi0N[0];

        F2rap[0].zeros(convF2.n_rows);
        FLrap[0].zeros(convFL.n_rows);
    }
    else { //Dummy start for DGLAP
        for(int y = 0; y <= start; ++y)
            PhiRapN[y] = Phi0N[y];
    }

    cout << "Evolving start,Nrap " << start+1<<" "<< Nrap << endl;

    for(int y = start+1; y < Nrap; ++y) {
        //Starting point of evol with 0.5 (Trapezius)
        cout << "\rEvolving " << y << flush;

        arma::vec yTemp(N, arma::fill::zeros);
        arma::vec myVec(N, arma::fill::zeros);

        F2rap[y] = 0.5*convF2.slice(y) * PhiRapN[0];// + convF2.slice(0) * PhiRapN[y]);
        FLrap[y] = 0.5*convFL.slice(y) * PhiRapN[0];// + convFL.slice(0) * PhiRapN[y]);

        //openMPI treatment
        int start=1, end;

    #if hasGPU == true //with GPU
            if(y > 1) {
                gpu.ConvoluteAll(y);
                //gpu.GetResult(y, yTemp);
                gpu.GetResults(y, yTemp, F2rap[y], FLrap[y]);

            }
            /*
               cout << "Radek start "<<y << endl;
            //cout << myVec << yTemp << endl;
            for(int i = 0; i < myVec.n_elem; ++i)
            cout << i <<" " << myVec(i) <<" "<< yTemp(i) << endl;
            cout << "Radek end " <<y<< endl;
            if(y > 2) exit(0);
            */
     #else //no GPU
        //if(!doGPU) {
            tie(start,end) = GetStartEnd(withMPI, 1, y-1); //from 1 to y-1
            //cout << "Start+end|nrap " << start <<" "<< end <<"|" << y-1<< endl;
            //Remaining without mult
            for(int d = start; d <= end; ++d)
                yTemp += matN.slice(d) * PhiRapN[y-d];

            if(withMPI) MPI_Allreduce(MPI_IN_PLACE, yTemp.memptr(), yTemp.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        //}
    #endif

        yTemp += 0.5 * matN.slice(y) * PhiRapN[0];


        //Whole right hand side
        //yTemp = stepY * yTemp + Phi0N[y];
        yTemp += Phi0N[y];

        //cout << "radek "<< y << endl;
        //cout << yTemp << endl;

        if(y % 100 == 0)
            cout <<"Rap point " << y << endl;
        /*
           auto matNow = mat[0]; //Adding diagonal DGLAP term
           for(int k = 0; k < N; ++k)
           for(int l = 0; l < N; ++l)
           matNow[k][l] += matDiag[y][k][l];
           */

        //arma::mat matEq = arma::mat(N,N,arma::fill::eye) - 0.5*stepY*matN.slice(0) -matNDiag.slice(y);

        //vector<double> Pred = PhiRap[y-1];
        //if(y >= 4)
        //for(int k = 0; k < N; ++k)
        //Pred[k] += PhiRap[y-1][k] - PhiRap[y-2][k];


        //PhiRap[y] = GetRegSolution(matEq, yTemp, Pred);
        //PhiRapN[y] = GetLinSolution(matEq, yTemp);
        PhiRapN[y] = matNInv.slice(y) * yTemp;

        F2rap[y] += 0.5*convF2.slice(0) * PhiRapN[y];
        FLrap[y] += 0.5*convFL.slice(0) * PhiRapN[y];

        #if hasGPU == true
            gpu.SetPhi(y, PhiRapN[y]);
        #endif

        //PhiRapN[y] =  IterSolution(matEq, double factor, const arma::vec &y);


        bool isGood = true;

        //cout << "Id of y = " << y << endl;
        if(!isGood && 0) {

            for(int j = 0; j < 10; ++j)
                cout << y <<" " << j <<" "<< PhiRapN[y](j) << endl;

        }
    }

}

//--------------------------------------------------------------------
void Solver::CalcF2L()
{
    F2rap.resize(S.Nrap);
    FLrap.resize(S.Nrap);
    F2rap[0] = arma::vec(convF2.n_rows, arma::fill::zeros);
    FLrap[0] = arma::vec(convFL.n_rows, arma::fill::zeros);
    for(int y = 1; y < S.Nrap; ++y) {
        //kT spectrum for particular bin y
        //Starting point of evol with 0.5 (Trapezius)
        F2rap[y] = 0.5*(convF2.slice(y) * PhiRapN[0] + convF2.slice(0) * PhiRapN[y]);
        FLrap[y] = 0.5*(convFL.slice(y) * PhiRapN[0] + convFL.slice(0) * PhiRapN[y]);
        //Remaining without mult by 0.5, convF2 should contain deltaRap factor
        for(int d = 1; d < y; ++d) {
            F2rap[y] += convF2.slice(d) * PhiRapN[y-d];
            FLrap[y] += convFL.slice(d) * PhiRapN[y-d];
        }
    }
}







//--------------------------------------------------------------------
void Solver::DoIteration()
{
    //const double stepY = (rapMax - rapMin) / (Nrap-1);

    //vector<vector<double>> PhiRapNew(Nrap);
    vector<arma::vec> PhiRapNew(S.Nrap);


    PhiRapNew[0] = Phi0N[0];// + matNDiag[0] * PhiRapN[0];

    for(int y = 1; y < S.Nrap; ++y) {

        //kT spectrum for particular bin y
        //Starting point of evol with 0.5 (Trapezius)
        arma::vec yTemp = 0.5*(matN.slice(y) * PhiRapN[0] + matN.slice(0) * PhiRapN[y]);

        //Remaining without mult
        for(int d = 1; d < y; ++d) {
            yTemp += matN.slice(d) * PhiRapN[y-d];
        }

        //Whole right hand side
        yTemp = yTemp + Phi0N[y];

        //Diag part (=virtual DGLAP term)
        //yTemp += matNDiag.slice(y) * PhiRapN[y];

        PhiRapNew[y] = yTemp;
    }
    PhiRapN = PhiRapNew;
}


//--------------------------------------------------------------------
void Solver::RunIterations(int Niter, bool init)
{
    //Init PhiRap 
    if(init) {
        assert(Phi0N.size() == PhiRapN.size());
        for(int y = 0; y < PhiRapN.size(); ++y)
            PhiRapN[y] = Phi0N[y];
    }

    //Do itrerations
    for(int i = 0; i < Niter; ++i) {
        DoIteration();
        cout << "Iteration " << i <<" done." << endl;
        cout << "Phi[kT0] = " << PhiRapN[S.Nrap-1](0) << endl;
    }
}


//--------------------------------------------------------------------
arma::vec Solver::GetRHS(const arma::vec &PHI)
{
    arma::vec dPhi = matN.slice(0) * PHI;
    return dPhi;
}


//--------------------------------------------------------------------
void Solver::Step(double delta) {
    static int y = 0;
    if(y == 0) PhiRapN[0] = Phi0N[0];

    arma::vec dPhi = GetRHS(PhiRapN[0]);
    ++y;
    PhiRapN[0] += delta * dPhi;

}


//Function of x and kT2
//--------------------------------------------------------------------
void Solver::InitF(function<double(double, double)> fun) {
    const double stepY = (S.rapMax-S.rapMin) / (S.Nrap-1);

    //New version
    Phi0N.resize(S.Nrap);
    for(int y = 0; y < S.Nrap; ++y) {
        arma::vec temp(S.Nint);
        double x = exp(-y*stepY);
        for(int i = 0; i < S.Nint; ++i) {
            double kT2 = exp(nod.xi[i]);
            temp(i) = fun(x, kT2);
        }
        Phi0N[y] = redMat * temp;
    }

    PhiRapN.resize(S.Nrap, arma::vec(S.N, arma::fill::zeros));
}

//SetSolution of x and kT2
//--------------------------------------------------------------------
void Solver::SetSolution(function<double(double, double)> fun)
{
    const double stepY = (S.rapMax-S.rapMin) / (S.Nrap-1);

    //New version
    PhiRapN.resize(S.Nrap);
    for(int y = 0; y < S.Nrap; ++y) {
        arma::vec temp(S.Nint);
        double x = exp(-y*stepY);
        for(int i = 0; i < S.Nint; ++i) {
            double kT2 = exp(nod.xi[i]);
            temp(i) = fun(x, kT2);
        }
        PhiRapN[y] = redMat * temp;
        //cout << redMat << endl;
        //exit(0);
    }
}



//--------------------------------------------------------------------
double Solver::Interpolate(double y, double L)
{

    y = max(y, S.rapMin);
    y = min(y, S.rapMax);

    double stepY = (S.rapMax-S.rapMin) / (S.Nrap-1);
    int yId = (y - S.rapMin) / stepY;
    int LId = 0;
    for(LId = S.N-1; LId >= 0; --LId)
        if(nodBase.xi[LId] <= L) break;
    //--LId;
    LId = max(0, LId);
    yId = min(S.Nrap-2, yId);

    double yLeft = S.rapMin + yId*stepY;
    double LLeft = nodBase.xi[LId];
    double LRight = nodBase.xi[LId+1];


    //if (LRight - L < 0 && L - LRight < 1e-7)
    //L = LRight;

    assert(y - yLeft >= 0);
    assert(yLeft+stepY - yLeft >= 0);



    if (LRight - L < 0 || L - LLeft < 0) {
        cout <<"Error " <<setprecision(24)<< L <<" : "<< LLeft <<" "<< LRight << endl;
        //L = LRight;
        //cout << "Ahojky " << nod.xi[0] << " "<< nod.xi[Nrap-1] << endl;
    }
    assert(L - LLeft >= 0);
    assert(LRight - L >= 0);

    if (yId > S.Nrap - 2) {
        cout << "Ahoj " << y << " :  " << yLeft << " "<< yLeft + stepY << endl;
    }

    assert(yId <= S.Nrap - 2);


    assert(LId <= S.N - 2);
    assert(yId >= 0);
    assert(LId >= 0);

    double fLL, fLR, fRL, fRR;

    fLL = PhiRapN[yId](LId);
    fLR = PhiRapN[yId](LId+1);
    fRL = PhiRapN[yId+1](LId);
    fRR = PhiRapN[yId+1](LId+1);

    bool canBeLog = fLL > 0 && fLR > 0 && fRL > 0 && fRR > 0;

    if(canBeLog) {
        fLL = log(fLL);
        fLR = log(fLR);
        fRL = log(fRL);
        fRR = log(fRR);
    }

    //cout << exp(fLL) << " "<< exp(fLR) <<" : "<< exp(fRL) << " "<< exp(fRR)<<  endl;

    double sum = 
        (y - yLeft)*(L - LLeft)*fRR  + (yLeft + stepY - y)*(L - LLeft)*fLR +
        (y - yLeft)*(LRight - L)*fRL + (yLeft + stepY - y)*(LRight - L)*fLL;
    sum /= stepY * (LRight - LLeft);
    assert(stepY != 0);
    assert((LRight-LLeft) != 0);

    if(canBeLog)
        return exp(sum);
    else 
        return sum;

}



//--------------------------------------------------------------------
void Solver::PrintBaseGrid()
{
    double stepY = (S.rapMax - S.rapMin) / (S.Nrap-1);
    double stepL = (S.Lmax - S.Lmin) / (S.N-1);


    for(int y = 0; y < S.Nrap; ++y)
        for(int l = 0; l < S.N; ++l) {
            double yNow = S.rapMin + y*stepY;
            double L  = nodBase.xi[l];

            double kt2 = exp(L);
            double x = exp(-yNow);

            cout << x << " "<< kt2 <<"  " <<  PhiRapN[y](l) << endl;
        }
}

//--------------------------------------------------------------------
void Solver::PrintGrid()
{
    double stepY = (S.rapMax - S.rapMin) / (S.Nrap-1);

    //cout <<"Interpolation "<< Interpolate(0, log(1) ) << endl;
    //return;

    double yStepMich = log(1e-2/1e-8) / 99;//(Nrap -1);
    double LStepMich = log(1e6/1e-2)  / 99;// (N -1);

    cout << "Radek " << endl;
    //for(int y = 0; y < 100; ++y) 
    for(int y = 99; y >= 0; --y) {
        double x = 1e-2 * pow(1e-8/1e-2, y/99.);
        double yNow = y*yStepMich;
        for(int i = 0; i < 100; ++i) {
            double kT2 = 0.01 * pow(1e6/0.01, i/99.);
            double LNow = i*LStepMich + log(0.01);
            if(i == 99) LNow = log(1e6)-1e-9;
            double res = Interpolate(yNow, LNow);
            /*if(yNow == 0)*/ cout << x<<" "<< kT2<<" "<< res <<  endl;
        }
    }
}


//--------------------------------------------------------------------
void Solver::PrintReduce()
{
    const double stepY = (S.rapMax - S.rapMin) / (S.Nrap-1);

    vector<double> q2Arr={0.15, 0.2, 0.25, 0.35, 0.4, 0.5, 0.65, 0.85, 1.2, 1.5, 2, 2.7, 3.5, 4.5,
        6.5, 8.5, 10, 12, 15, 18, 22, 27, 35, 45, 60, 70, 90, 120, 150, 200, 250,
        300, 400, 500, 650, 800, 1000, 1200, 1500, 2000, 3000, 5000, 8000, 12000, 20000, 30000};

    double sBeam = pow(318.12,2);

    //assert(F2rap.size() == matN.n_slices);
    //assert(FLrap.size() == matN.n_slices);

    for(int rapID = 0; rapID < S.Nrap; ++rapID) {
        double rap = S.rapMin + rapID*stepY;
        double x = exp(-rap);

        for(int i = 0; i < 46; ++i) {
            double Q2 = q2Arr[i];
            double y = Q2/(x*sBeam);

            double yPlus = 1+(1-y)*(1-y);

            //cout << " rapID " << rapID << " " << F2rap[rapID].n_rows << endl;
            //assert(F2rap[rapID].n_rows == 46);
            //assert(FLrap[rapID].n_rows == 46);
            double sRed = F2rap[rapID](i) - y*y/yPlus * FLrap[rapID](i);

            //cout << x << " "<< Q2 <<" "<< sRed << endl;
            cout << x << " "<< Q2 <<" "<< F2rap[rapID](i) <<" "<< FLrap[rapID](i) << endl;

        }
    }
}


vector<double> Solver::GetRegSolution(const vector<vector<double>> &MatEq, const vector<double> &y, const vector<double> &yReg)
{
    int N = S.N;
    //cout 
    arma::mat M(N,N);
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j){
            M(i,j) = MatEq[i][j];
        }

    arma::vec v(N), vReg(N);
    for(int i = 0; i < N; ++i) {
        v(i) = y[i]; 
        vReg(i) = 0;//yReg[i]; 
    }
    arma::mat L(N-1,N);

    for(int i = 0; i < N-1; ++i) 
        for(int j = 0; j < N; ++j) {
            if(i == j)
                L(i,j) = 1.;///(pow(yReg[i],2) + 1e-10);
            if(i == j - 1)
                L(i,j) =-1.;///(pow(yReg[i],2) + 1e-10);
            else
                L(i,j) = 0.;///(pow(yReg[i],2) + 1e-10);
        }
    arma::mat L2 = trans(L)*L;

    //double tau = 3e-1;
    double tau = 5e-1;
    arma::mat E = arma::trans(M)*M + tau*L2;
    //for(int i = 0; i < N; ++i)
    //E(i,i) += tau * L(i,i);

    //cout << " " << Ep(3,3) << " "<< Ep(3,4) << endl;
    //arma::mat E = arma::trans(M)*M + tau;
    //cout << " " << E(3,3) << " "<< E(3,4) << endl;
    //assert(0);
    arma::vec r = arma::trans(M)*v + tau*L2*vReg;

    arma::vec res = arma::solve(E,r);

    vector<double> Res(N);
    for(int i = 0; i < N; ++i)
        Res[i] = res[i];
    return Res;

}



arma::vec Solver::IterSolution(const arma::mat &Mat, double factor, const arma::vec &y)
{
    arma::vec yNow = y;
    arma::vec ySum = y;

    double diff = 0;
    for(int i = 0; i < 20; ++i) {
        yNow = factor * Mat * yNow;

        diff = 0;
        for(int k = 0; k < S.N; ++k)
            diff = max(diff, abs(yNow[k]) / (1e-13+abs(yNow(k)) + abs(ySum(k))));
        //cout <<"Diff " <<  i << " "<< diff << endl;

        ySum += yNow;

        //cout << "Diff is " << diff << endl;
        if(diff < 1e-9) break;
    }
    assert(diff < 1e-5);


    return ySum;
}

