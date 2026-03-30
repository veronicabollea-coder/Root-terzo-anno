#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>

void graficoRate() {

    // --- DATI DI ESEMPIO ---
    const int N = 3;
    double x[N]     = {9.8, 19.5, 29.9};        // variabile indipendente
    double y[N]     = {0.45, 0.55, 1.24}; // variabile misurata
    double ex[N]    = {0, 0, 0};        // errori in x (qui nulli)
    double ey[N]    = {0.03, 0.03, 0.07}; // errori verticali

    // --- CREAZIONE DEL GRAFICO CON ERRORI ---
    TGraphErrors* gr = new TGraphErrors(N, x, y, ex, ey);

    gr->SetTitle("Rate in funzione della distanza;d[cm];R(d)[Hz]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetLineColor(kBlue + 1);
    gr->SetMarkerColor(kBlue + 1);

    // --- CANVAS ---
    TCanvas* c = new TCanvas("c", "Rate in funzione di d", 800, 600);
    c->SetGrid();

    // --- DISEGNO ---
    gr->Draw("AP"); // A = assi, P = punti (+ barre)

    c->Update();
}