#include<iostream>
#include<vector>
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"




using namespace std;




void calibrazione_spettrometro(){
    //
    //
    vector<double> lambda={
        365.015,
        404.655,
        435.833,
        546.074,
        579.066,
        696.543,
        763.511,
        811.531,
        842.465,
        912.297,
    };
    vector<double> chn={
        116.72,
        299.46,
        445.3,
        972.2,
        1128.8,
        1725.7,
        2075.8,
        2331.9,
        2499.57,
        2895.5,
    };
    vector<double> err_lambda={};
    vector<double> err_chn={
        0.304166666666667,
        6.5,
        6.3,
        7,
        19.24,
        7.7,
        8.4,
        10.5,
        14.33,
        8.9,
    };
    //
    TCanvas* c1= new TCanvas("c1", "curva di calibrazione spettrometro", 800, 600);
    c1->SetGrid();
    //
    //lunghezza degli array
    int len= lambda.size();
    //
    //NB: se qualcosa qua non ha errori, al posto di err_chn o err_lambda, mettere nullptr al posto
    TGraphErrors* gr1= new TGraphErrors(len, chn.data(), lambda.data(), err_chn.data(), nullptr);

    gr1->SetTitle("Fit spettrometro; canali ; lambda [nm]");
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.3);
    gr1->SetLineColor(1);
    gr1->SetMarkerColor(1);
    //
    //eventualmente modifica gli estremi in cui è calcolata la TF1
    TF1* fit= new TF1("fit", "pol3(0)", 0, 3000);
    //
    //si renderà necessario fittare i parametri iniziale; fare sistema di 4 equazioni in 4 incognite per trovarli (usare matlab)
    fit->SetParameters(350, 0.2, 0, 0);
    gr1->Fit("fit", "R");
    gr1->Draw();  
    c1->Print("FC_calibrazione_spettrometro.png", "png");
cout<<"p-value é: " << fit->GetProb()<<endl;
}










/*MATLAB: il seguente codice lo si può usare per risolvere un sistema di 4 equazioni in 4 incognite in modo veloce
%definizione dei parametri: la nostra funz. di fit è a+bx+cx^2+dx^3; le x saranno prese a punti diversi, quinid useremo x1, x2, x3, x4 per le 4 funioni


syms x1 x2 x3 x4 a b c d


x1=0; % mettici il valore del channel in cui trovi la lambda corrispettiva
eq1= a*x1+b*x1+c*x1^2+d*x1^3==0; % al posto di 0 mettici il valore atteso di lambda  


x2=0;
eq2=  a+b*x2+c*x2^2+d*x2^3==0;


x3=0;
eq3=  a+b*x3+c*x3^2+d*x3^3==0;


x4=0;
eq4=  a+b*x4+c*x4^2+d*x4^3==0;


sol= solve([eq1, eq2, eq3, eq4], [a, b, c, d]);


*/
