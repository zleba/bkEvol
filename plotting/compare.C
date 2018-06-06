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

    return grMap;

}


map<double, TGraph*> ReadMichalFile(const char *fName, map<double, TGraph*>  myGr)
{
    ifstream file(fName);

    map<double, TGraph*> grMap;
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
