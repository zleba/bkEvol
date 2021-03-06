<<<<<<< HEAD
#define SF TString::Format

map<double, TGraph*> LoadDataRadek(TString fname)
{
    ifstream file(fname.Data());

    double x, q2, f2, fl;
    map< double, map<double, double> > points;

    map<double, TGraph*> grMap;

    while(1) {
        file >> x >> q2 >> f2 >> fl;

        points[q2][x] = f2;


        if(grMap.count(q2) == 0) {
            grMap[q2] = new TGraph();
        }

        grMap[q2]->SetPoint(grMap[q2]->GetN(), x, f2);
        //grMap[q2]->SetTitle(SF("k_{T}^{2}=%g", kT2));


        if(!file.good()) break;
    }



    for(auto &p : points)
        cout << p.first <<" "<< p.second.size() << endl;// p.second.size() << endl;

=======
map<double, TGraph*> ReadFile(const char *fName)
{
    ifstream file(fName);

    map<double, TGraph*> grMap;
    int i = 0;
    while(1) {
        string str;
        getline(file, str);
        if(!file.good()) break;
        if(str[0] == '#') continue;
        if(str.size() < 8) continue;
        stringstream  sStream(str);
        double kT2, b, y, Phi;
        sStream >> kT2 >> b >> y >> Phi;
        kT2 *= kT2;

        if(grMap.count(y) == 0) {
            i = 0;
            grMap[y] = new TGraph();
        }
        grMap[y]->SetPoint(i++, kT2, Phi);
        cout << kT2 <<" "<< b<<" "<<y << " "<< Phi << endl;

    }
>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511
    return grMap;

}

<<<<<<< HEAD

map<double, TGraph*> ReadMichalFile(const char *fName, map<double, TGraph*>  myGr)
=======
map<double, TGraph*> ReadMichalFile(const char *fName)
>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511
{
    ifstream file(fName);

    map<double, TGraph*> grMap;
<<<<<<< HEAD
    map< double, map<double, double> > points;

    while(1) {
        if(!file.good()) break;
        double q2, y, f2, x;
        file >> x >> q2 >> f2;
        y = -log(x);

        points[q2][x] = f2;

    }


    for(auto &g : myGr) {
        double q2 = g.first;

        auto First = points.begin();
        auto Last = points.rbegin();

        while( First->first < q2 && First != points.end()) ++First;
        while( Last->first > q2 && Last != points.rend()) ++Last;

        if(First == points.end()) continue;
        if(Last == points.rend()) continue;

        --First;
        --Last;

        cout << First->first <<" "<< Last->first <<" "<< q2<<  endl;

        double dist1 =  log(q2 / First->first) /  log(Last->first / First->first);
        double dist2 = -log(q2 / Last->first) /  log(Last->first / First->first);

        grMap[q2] = new TGraph();

        for(auto xVal : First->second) {
           double x = xVal.first;
           double f2_1 = xVal.second;
           double f2_2 = (Last->second)[x];

           double res = (f2_1 * dist2 + f2_2 * dist1) / (dist1 + dist2);

           grMap[q2]->SetPoint(grMap[q2]->GetN(), x, res);

           cout <<  xVal.first << xVal.second << endl;


        }

        //cout << p.first <<" "<< p.second.size() << endl;// p.second.size() << endl;

    }


    return grMap;

}

void DrawWithRat(TGraph *gr1, TGraph *gr2, double q2)
{
   TVirtualPad *pad = gPad;

   pad->Divide(2,1, 0.0001, 0.00001);

   pad->cd(1);
   gPad->SetLeftMargin(0.15);

   gr1->SetLineColor(kBlack);
   gr2->SetLineColor(kBlue);
   gr1->Draw("alp");
   gr2->Draw("lp same");

   gr1->GetXaxis()->SetTitle("x");
   gr1->GetYaxis()->SetTitle("F_{2}");

   gr1->SetMinimum(1e-5);
   gPad->SetLogx();
   gPad->SetLogy();

   TLegend *leg = new TLegend(0.4, 0.7-0.4, 0.6, 0.9-0.4);
   leg->SetBorderSize(0);
   leg->SetHeader(SF("Q^{2} = %g", q2));
   leg->AddEntry(gr1, "Radek", "l");
   leg->AddEntry(gr2, "Michal", "l");
   leg->Draw();

   gPad->Update();
   gPad->RedrawAxis();
   gPad->Update();

   pad->cd(2);
   gPad->SetLeftMargin(0.15);

   TGraph *grRat = new TGraph;

   for(int i = 0; i < gr1->GetN(); ++i) {
        double x, val1, val2;
        gr1->GetPoint(i, x, val1);
        val2 = gr2->Eval(x);

        grRat->SetPoint(i, x, val1/val2);
   }
   grRat->Draw();
   gPad->SetLogx();
   grRat->SetMaximum(5);
   grRat->SetMinimum(0);

   grRat->GetXaxis()->SetTitle("x");
   grRat->GetYaxis()->SetTitle("#frac{F_{2}^{Rad}}{F_{2}^{Mich}}");

}


void compare()
{

    auto grRad = LoadDataRadek("../vystupF2new");
    //auto grMich = ReadMichalFile("/afs/desy.de/user/z/zlebcr/temp/F2-M-091017.dat", grRad);
    auto grMich = ReadMichalFile("/afs/desy.de/user/z/zlebcr/temp/F2-191017-2.dat", grRad);

    //graphs[10]->Draw("acp");
    //temp[10]->Draw("cp same");

    TCanvas *can = new TCanvas("can");
    can->SaveAs("ahojjj.pdf[");
    for(auto gMich : grMich) {
        can->Clear();
        double q2 = gMich.first;
        DrawWithRat(grRad[q2], grMich[q2], q2);
        can->SaveAs("ahojjj.pdf");
    }
    can->SaveAs("ahojjj.pdf]");

}
=======
    int i = 0;
    while(1) {
        string str;
        getline(file, str);
        if(!file.good()) break;
        if(str[0] == '#') continue;
        if(str.size() < 8) continue;
        stringstream  sStream(str);
        double kT2, y, Phi, x;
        sStream >> x >> kT2 >> Phi;
        y = -log(x/1e-2);

        if(grMap.count(y) == 0) {
            i = 0;
            grMap[y] = new TGraph();
        }
        grMap[y]->SetPoint(i++, kT2, Phi);
        //cout << kT2 <<" "<< b<<" "<<y << " "<< Phi << endl;

    }
    return grMap;

}


map<double, TGraph*> ReadFileRadek(const char *fName)
{
    ifstream file(fName);

    map<double, TGraph*> grMap;
    int i = 0;
    while(1) {
        string str;
        getline(file, str);
        if(!file.good()) break;
        if(str[0] == '#') continue;
        if(str.size() < 8) continue;
        if(str[0] < '0' || str[0] > '9') continue;
        stringstream  sStream(str);
        double kT2, y, Phi;
        sStream >> kT2 >>  Phi;
        kT2 *= kT2;
        y = log(1e6);

        if(grMap.count(y) == 0) {
            i = 0;
            grMap[y] = new TGraph();
        }
        grMap[y]->SetPoint(i++, kT2, Phi);
        cout << kT2 <<" "<< y << " "<< Phi << endl;

    }
    return grMap;

}


void DrawRatio(vector<TGraph*> grs)
{
    vector<TGraph *> grsRat(grs.size());
    for(auto &gr : grsRat) gr = new TGraph();
    for(int i = 0; i < grs.size(); ++i)
        grsRat[i]->SetLineColor(grs[i]->GetLineColor());


    TF1 *f = new TF1("f", "exp([0] + [1]*log(x) + [2] *log(x)^2 + [3]*log(x)^3 + [4]*log(x)^4)", 0.01*1.1, 1e5/1.0);
    grs[0]->Fit(f);

    for(int j = 0; j < grs.size(); ++j) {
    for(int i = 0; i < grs[j]->GetN(); ++i) {
        double x, y;
        grs[j]->GetPoint(i, x, y);
        double yRef = f->Eval(x);
        //double y = grs[j]->Eval(x);
        assert(grsRat[j]);
        cout << "RADEK " << grsRat[j] << " "<< i << endl;
        grsRat[j]->SetPoint(i, x, y/yRef);
        }
    }


    /*
    for(int i = 0; i < grs[0]->GetN(); ++i) {
        double x, y0;
        grs[0]->GetPoint(i, x, y0);
        //y0 = f->Eval(x);

        for(int j = 0; j < grs.size(); ++j) {
            assert(grs[j]);
            double y = grs[j]->Eval(x);
            assert(grsRat[j]);
            cout << "RADEK " << grsRat[j] << " "<< i << endl;
            grsRat[j]->SetPoint(i, x, y/y0);
        }
    }
    */

    TCanvas *can = new TCanvas("canRat", "canvas");
    gPad->SetLogx();
    //gPad->SetLogy();

    grsRat[0]->Draw("*ca");
    grsRat[1]->Draw("*c same");
    //grsRat[2]->Draw("c same");

    grsRat[0]->GetXaxis()->SetTitle("k_{T}^{2}");
    grsRat[0]->GetXaxis()->SetTitleOffset(1.3);
    grsRat[0]->GetYaxis()->SetTitle("#phi/#phi_{Radek}");
    grsRat[0]->GetYaxis()->SetRangeUser(0.6, 1.3);

    TLegend *leg = new TLegend(0.3, 0.6, 0.5, 0.8);
    leg->SetBorderSize(0);
    leg->SetTextSize(grsRat[0]->GetYaxis()->GetLabelSize());
    leg->SetHeader("Evolution of k_{T}^{2} exp(-k_{T}^{2})");
    //leg->AddEntry(grsRat[0], "BKsolver solution (x=1e-8)", "l");
    //leg->AddEntry(grsRat[1], "Michal's solution (x=1e-8)", "l");
    //leg->AddEntry(grsRat[2], "Radek's solution (x=1e-8)", "l");

    leg->AddEntry(grsRat[0], "Radek's solution of 81 (x=1e-8)", "l");
    leg->AddEntry(grsRat[1], "Michal's solution of 81 (x=1e-8)", "l");

    leg->Draw();

    //can->SaveAs("eq81rat.eps");

}


void compare()
{
    
    auto graphsBK = ReadFile("../bkresult.dat");
    auto graphsML = ReadMichalFile("/home/radek/Dropbox/system-of-equations/note_tex/new-grids-tests/eq81M.dat");
    //auto graphsRA = ReadFileRadek("/home/radek/Dropbox/system-of-equations/note_tex/new-grids-tests/RADEK_kt2expmKT2/res79normal");
    auto graphsRA = ReadFileRadek("/home/radek/Dropbox/Krakow/iterate/results/subKer81high");

    cout << "Graphs " <<  graphsML.size() << endl;
    for(auto &gr : graphsML) {
        cout << gr.first << endl;
    }
    TGraph *grBK = graphsBK.rbegin()->second;
    TGraph *grML = graphsML.rbegin()->second;
    TGraph *grRA = graphsRA.rbegin()->second;

    grML->SetLineColor(kRed);
    grRA->SetLineColor(kBlue);

    //DrawRatio({grBK, grML, grRA});
    DrawRatio({grRA, grML});

    //return;

    TCanvas *can = new TCanvas("can", "canvas");
    gPad->SetLogx();
    gPad->SetLogy();

    //graphsBK[1.381551e+01]->Draw("ac");
    grBK->Draw("ac");
    grBK->GetYaxis()->SetRangeUser(1e-3, 1e3);
    grBK->GetXaxis()->SetTitle("k_{T}^{2}");
    grBK->GetYaxis()->SetTitle("#phi");
    grBK->GetXaxis()->SetTitleOffset(1.3);
    //grML->Draw("c same");
    grRA->Draw("c same");

    //TF1 *f = new TF1("f", "exp([0] + [1]*log(x) + [2] *log(x)^2 + [3]*log(x)^3 + [4]*log(x)^4)", 0.01*1.3, 1e6/1.3);
    //grRA->Fit(f);


    TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.8);
    leg->SetBorderSize(0);
    leg->SetTextSize(grBK->GetYaxis()->GetLabelSize());
    leg->SetHeader("Evolution of k_{T}^{2} exp(-k_{T}^{2})");
    leg->AddEntry(grBK, "BKsolver solution (x=1e-8)", "l");
    leg->AddEntry(grML, "Michal's solution of 81 (x=1e-8)", "l");
    leg->AddEntry(grRA, "Radek's solution of 81 (x=1e-8)", "l");
    leg->Draw();

    //can->SaveAs("eq81.eps");
}



>>>>>>> 3218ac9a0ae3dc07b896fb991fac340ae8181511
