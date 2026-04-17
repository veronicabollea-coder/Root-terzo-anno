#include<iostream>
#include<vector>
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"




using namespace std;




void fit_monocromatore(){
    //
    //
    vector<double> chn={
        461.32,
        526.59,
        604.66,
        680.37,
        742.93,
        809,
        880.43,
        953.96,
        1023.26,
        1101.2,
        1172.6,
        1248.33,
        1320,
        1397.7,
        1472.66,
        1542.9,
        1621.9,
        1696.76,
        1768.22,
        1844.1,
        1922.86,
    };
    vector<double> ua={
        500,
        515,
        530,
        545,
        560,
        575,
        590,
        605,
        620,
        635,
        650,
        665,
        680,
        695,
        710,
        725,
        740,
        755,
        770,
        785,
        800,
    };
    vector<double> err_chn={
        116.65,
        168,
        83,
        69.85,
        101.6,
        148.21,
        164.9,
        186.7,
        173,
        169.9,
        127.9,
        133.3,
        103.5,
        82.5,
        85.66,
        93.38,
        94.35,
        91.2,
        72.64,
        71.38,
        57.6,
    };
    vector<double> err_ua={
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
    };
    //
    TCanvas* c1= new TCanvas("c1", "curva di calibrazione monocromatore", 800, 600);
    c1->SetGrid();
    //
    //lunghezza degli array
    int len= ua.size();
    //
    //NB: se qualcosa qua non ha errori, al posto di err_chn o err_lambda, mettere nullptr al posto
    TGraphErrors* gr1= new TGraphErrors(len, chn.data(), ua.data(), err_chn.data(), err_ua.data());
    //
    gr1->SetTitle("Fit picchi monocromatore; canali ; u.a.");
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.7);
    gr1->SetLineColor(1);
    gr1->SetMarkerColor(1);
    //eventualmente modifica gli estremi in cui è calcolata la TF1
    TF1* fit= new TF1("fit", "[0]+x*[1]", 0, 3000);
    //
    //si renderà necessario fittare i parametri iniziale; fare sistema di 4 equazioni in 4 incognite per trovarli (usare matlab)
    fit->SetParameters(350, 0.2, 0, 0);
    gr1->Fit("fit", "R");
    gr1->Draw();  
    c1->Print("FC_fit_monocromatore.png", "png");
cout<<"p-value é: " << fit->GetProb()<<endl;
}
