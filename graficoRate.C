#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>

void graficoRate() {

    // --- DATI DI ESEMPIO ---
    const int N = 2;
    double x[N]     = {10.1, 19.8};        // variabile indipendente
    double y[N]     = {2.7, 1.21 }; // variabile misurata
    double ex[N]    = {0, 0};        // errori in x (qui nulli)
    double ey[N]    = {0.3, 0.15 }; // errori verticali

    // --- CREAZIONE DEL GRAFICO CON ERRORI ---
    TGraphErrors* gr = new TGraphErrors(N, x, y, ex, ey);

    gr->SetTitle("Rate normalizzato in funzione della distanza;d[cm];R(d)/R(30)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetLineColor(1);
    gr->SetMarkerColor(1);

    // --- CANVAS ---
    TCanvas* c = new TCanvas("c", "Rate normalizzato in funzione di d", 800, 600);
    c->SetGrid();

    // --- DISEGNO ---
    gr->Draw("AP"); // A = assi, P = punti (+ barre)

    c->Update();
    c->Print("gabibbo.png", "png");

}