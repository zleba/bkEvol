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
    return grMap;

}

map<double, TGraph*> ReadMichalFile(const char *fName)
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



