#include "PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(PlottingHelper/plottingHelper_C.so)
using namespace PlottingHelper;

TGraph2D* Transform(map<double, TGraph*>   grMap)
{
    TGraph2D *g2D = new TGraph2D();
    int k = 0;
    for(auto &gr : grMap) {
        double y = gr.first;
        TGraph *grNow = gr.second;
        for(int i = 0; i < grNow->GetN(); ++i) {
           double kt2, val;
           grNow->GetPoint(i, kt2, val); 
            val = abs(val);
           g2D->SetPoint(k++, y, kt2, val);
        }
    }

    return g2D;
}


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

pair<map<double, TGraph*>, map<double, TGraph*>> ReadMichalFileBoth(const char *fName)
{
    ifstream file(fName);

    map<double, TGraph*> grMapY, grMapKt;

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
        cout << "x,kt2,phi " << x <<" "<< kT2 <<" "<< Phi << endl;

        if(grMapY.count(x) == 0) {
            grMapY[x] = new TGraph();
        }
        if(grMapKt.count(kT2) == 0) {
            grMapKt[kT2] = new TGraph();
        }

        grMapY[x]->SetPoint(grMapY[x]->GetN(), kT2, Phi);
        grMapY[x]->SetTitle(SF("x=%g", x));
        grMapY[x]->SetName(SF("x=%g", x));
        grMapKt[kT2]->SetPoint(grMapKt[kT2]->GetN(), x, Phi);
        grMapKt[kT2]->SetTitle(SF("k_{T}^{2}=%g", kT2));

        //cout << kT2 <<" "<< b<<" "<<y << " "<< Phi << endl;

    }
    return make_pair(grMapY, grMapKt);

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

TGraph *getNth(map<double, TGraph*> grMap, int n) {
    int k = 0;
    for(auto &gr : grMap) {
        if(k == n){
            cout << "xVal is " << gr.first << endl;
            return gr.second;
        }
        ++k;
    }

    return nullptr;

}

void DrawRatio(vector<TGraph*> grs)
{
    vector<TGraph *> grsRat(grs.size());
    for(auto &gr : grsRat) gr = new TGraph();
    for(int i = 0; i < grs.size(); ++i)
        grsRat[i]->SetLineColor(grs[i]->GetLineColor());


    TF1 *f = new TF1("f", "exp([0] + [1]*log(x) + [2] *log(x)^2 + [3]*log(x)^3 + [4]*log(x)^4)", 0.01*1.1, 1e2/1.0);
    grs[0]->Fit(f);

    for(int j = 0; j < grs.size(); ++j) {
    for(int i = 0; i < grs[j]->GetN(); ++i) {
        double x, y;
        grs[j]->GetPoint(i, x, y);
        double yRef = f->Eval(x);
        //double y = grs[j]->Eval(x);
        assert(grsRat[j]);
        //cout << "RADEK " << grsRat[j] << " "<< i << endl;
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




int cols[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9};

void PlotOverview(vector<map<double, TGraph*>> grMap, TString mode, double Min, double Max)
{
    gPad->SetLogx();
    gPad->SetLogy();

    vector<double> vals;
    for(auto &g : grMap[0])
        vals.push_back(g.first);

    gPad->SetRightMargin(0.2);

    int start = vals.size()*0.00;
    int end = vals.size()*1.00;

    int step = (mode == "kT") ? 30 : 70;
    int iCol =0;
    for(int i = start; i < end; i+=step) {
        if(iCol == 0) {
            grMap[0][vals[i]]->Draw("al");
        }
        else {
            grMap[0][vals[i]]->Draw("l same");
        }
        grMap[0][vals[i]]->SetLineColor(cols[iCol%9]);

        if(grMap.size() == 2) {
            grMap[1][vals[i]]->Draw("l same");
            grMap[1][vals[i]]->SetLineColor(cols[iCol%9]);
            grMap[1][vals[i]]->SetLineStyle(2);
        }
        ++iCol;
    }


    //grMap[vals[0]]->SetMinimum(Min);
    //grMap[vals[0]]->SetMaximum(Max);
    GetYaxis()->SetRangeUser(Min, Max);
    GetFrame()->SetTitle("");

    GetYaxis()->SetTitle("#phi");
    if(mode == "kT")
        GetXaxis()->SetTitle("x");
    else
        GetXaxis()->SetTitle("k_{T}^{2}");

    //gPad->Update();



    TLegend *leg = new TLegend(0.8, 0.3, 1.0, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry((TObject*)0, "Effect of DGLAP term", "h");

    leg->AddEntry((TObject*)0, "Solid - w/o DGLAP", "h");
    leg->AddEntry((TObject*)0, "Dashed - with DGLAP", "h");

    for(int i = start; i < end; i+=step) {
        double val =  vals[i];
        if(mode == "kT")
            leg->AddEntry(grMap[0][vals[i]], TString::Format("k_{T}^{2} = %3.3g", val), "l");
        else
            leg->AddEntry(grMap[0][vals[i]], TString::Format("x = %3.3g", val), "l");
    }
    leg->Draw();


}

void compareNew()
{
    
    //auto graphsBK = ReadFile("../bkresult.dat");

    map<double, TGraph*> graphsMLy, graphsMLkt;
    map<double, TGraph*> graphsRA1y, graphsRA1kt;
    map<double, TGraph*> graphsRA2y, graphsRA2kt;
    //tie(graphsMLy, graphsMLkt) = ReadMichalFileBoth("/home/radek/Dropbox/system-of-equations/note_tex/new-grids-tests/eq84M.dat");
    tie(graphsRA1y, graphsRA1kt) = ReadMichalFileBoth("bfklFull2.txt");// hope85highStat");
    tie(graphsRA2y, graphsRA2kt) = ReadMichalFileBoth("bfklFullCheb.txt");


    //auto graphsRA = ReadFileRadek("/home/radek/Dropbox/system-of-equations/note_tex/new-grids-tests/RADEK_kt2expmKT2/res79normal");
    //auto graphsRA = ReadFileRadek("/home/radek/Dropbox/Krakow/iterate/results/subKer81high");

    /*
    cout << "Graphs " <<  graphsML.size() << endl;
    for(auto &gr : graphsML) {
        cout << gr.first << endl;
    }
    TGraph *grBK = graphsBK.rbegin()->second;
    TGraph *grML = graphsML.rbegin()->second;
    TGraph *grRA = graphsRA.rbegin()->second;

    grML->SetLineColor(kRed);
    grRA->SetLineColor(kBlue);
    */

    //graphsMLy.rbegin()->second->Draw("a*c");

    //TCanvas *can1 = new TCanvas("can1", "canvas");
    //PlotOverview({graphsMLkt}, "kT", 1e-7, 1e2);
    //can1->SaveAs("allMichal.pdf");
    //TCanvas *can2 = new TCanvas("can2", "canvas");
    //PlotOverview({graphsMLy}, "x", 1e-7, 1e2);


    TCanvas *can3 = new TCanvas("can3", "canvas");
    can3->SetLogx();
    can3->SetLogy();
    graphsRA1y[1e-6]->Draw("apl");
    graphsRA2y[1e-6]->Draw("same pl");

    for(int i = 0; i < graphsRA2y[1e-6]->GetN(); ++i) {
        double x, y1, y2;
        graphsRA2y[1e-6]->GetPoint(i, x, y2);
        y1 = graphsRA1y[1e-6]->Eval(x);
        cout << x << " "<< y2/y1 << endl;
    }

    GetYaxis()->SetRangeUser(1e-6, 400);

    //PlotOverview({graphsRA1y, graphsRA2y}, "y", 1e-6, 2e2);
    can3->SaveAs("comparisonYfixed.pdf");
    return 0;
    TCanvas *can4 = new TCanvas("can4", "canvas");
    PlotOverview({graphsRA1kt, graphsRA2kt}, "kT", 1e-6, 2e2);
    can4->SaveAs("comparisonKTfixed.pdf");


    //DrawRatio( { getNth(graphsRA1y,150), getNth(graphsRA2y,150) });

    /*
    TGraph2D *gHope = Transform(graphsML);

    //gHope->Draw("TRI1");
    gHope->Draw("TRI1");
    gHope->SetMinimum(1e-3);
    gHope->SetMaximum(1e2);
    //gHope->Draw("col2");
    gPad->SetLogy();
    gPad->SetLogz();
    */


    return;

    /*
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
    */

    //can->SaveAs("eq81.eps");
}



