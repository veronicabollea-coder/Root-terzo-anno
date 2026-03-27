#include <iostream>
#include <iomanip>
#include <math.h>   

#include "TH1F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPad.h"

double pi = 3.141592;

void computeAcceptanceCosmics(unsigned int events = 1000, 
			      double lengthX = 14.5,              // in cm 
			      double widthY = 14.5,               // in cm
			      double spacingZ = 10.0,             // in cm
			      double distanceFromXWall = 2.0,    // in m
			      double distanceFromYWall = 3.0,    // in m
			      double alpha = 0.10,              // in m-1 
			      double deltaAlpha = 0.01)         // in m-1
{
  // compute cosmic ray acceptance and attenuation of a two-slab scintillator detector 
  // (slabs assumed as infinitely thin)

  bool drawTracks = true;

  // histograms
  TH1F* hxUp = new TH1F("hxUp", "x on upper slab", 100, -0.55*lengthX, 0.55*lengthX);
  TH1F* hyUp = new TH1F("hyUp", "y on upper slab", 100, -0.6*widthY, 0.6*widthY);
  TH1F* hxDown = new TH1F("hxDown", "x on lower slab", 100, -0.55*lengthX, 0.55*lengthX);
  TH1F* hyDown = new TH1F("hyDown", "y on lower slab", 100, -0.6*widthY, 0.6*widthY);
  TH1F* hphi = new TH1F("hphi", "muon phi all", 100, -3.5, 3.5);
  TH1F* hcosTheta = new TH1F("hcosTheta", "muon cos(theta) all", 100, -0.05, 1.05);
  TH1F* hphiAcc = new TH1F("hphiAcc", "muon phi accepted", 100, -3.5, 3.5);
  TH1F* hcosThetaAcc = new TH1F("hcosThetaAcc", "muon cos(theta) accepted", 100, -0.05, 1.05);
  TH1F* hpath = new TH1F("hpath", "path of accepted muons", 100, 0.95*spacingZ, 2.0*spacingZ);
 
  // tracks
  TLine* longTrack[1000];
  TLine* trasvTrack[1000];

  // Building data
  static const int nFloors = 8;
  // thickness of every floor from bottom to top
  float thick[nFloors]  = {0.22,0.46,0.34,0.36,0.54,0.39,0.53,0.86};        // in m
  float thickTotal = 3.7;                                                   // in m
  //floor height starting from 2nd floor underground
  float height[nFloors] = {3.28,6.79,10.98,14.48,17.98,21.48,25.13,28.25};  // in m

  // miscellaneous initialization
  TRandom3 *rnd = new TRandom3(1234);
  unsigned int npassedAcc = 0;
  unsigned int npassedAccAtt = 0;
  unsigned int npassedAccAttPlus = 0;
  unsigned int npassedAccAttMinus = 0;
  unsigned int ndrawn = 0;
  double meanPath = 0.;
  double meanPathSq = 0.;

  for (unsigned int ievent = 0; ievent < events; ievent++) { 

    // track color: red if in acceptance, green if out
    int icolor = 0;
    
    // randomly choose hit point on upper slab
    Double_t xUp = rnd->Uniform(lengthX) - lengthX/2.;   // rnd between -l/2 and l/2
    Double_t yUp = rnd->Uniform(widthY) - widthY/2.;     // rnd between -w/2 and w/2
    hxUp->Fill(xUp);
    hyUp->Fill(yUp);
    
    // randomly choose theta / phi of muon track 
    Double_t phi = rnd->Uniform(2.*pi) - pi;
    
    // For theta use inverted-cumulative MonteCarlo method
    Double_t cosTheta = pow(rnd->Uniform(),0.33333);

    hphi->Fill(phi);
    hcosTheta->Fill(cosTheta);

    // ATTENUATION: propagate back to upper floors
    bool isAbsorbedAlpha = false;    bool isAbsorbedAlphaPlus = false;   bool isAbsorbedAlphaMinus = false; 
    for (unsigned int ifloor = 0; ifloor < nFloors; ifloor++) {
      Double_t xFloor = (xUp/100.) - height[ifloor]*cos(phi)*sqrt(1-cosTheta*cosTheta)/cosTheta;
      Double_t yFloor = (yUp/100.) - height[ifloor]*sin(phi)*sqrt(1-cosTheta*cosTheta)/cosTheta;
      if (xFloor < distanceFromXWall && yFloor < distanceFromYWall) {    // floor is traversed, try absorption
	Double_t absorb = rnd->Uniform();
	if (absorb > exp(-alpha*thick[ifloor]/cosTheta)) isAbsorbedAlpha = true;
	if (absorb > exp(-(alpha+deltaAlpha)*thick[ifloor]/cosTheta)) isAbsorbedAlphaPlus = true;
	if (absorb > exp(-(alpha-deltaAlpha)*thick[ifloor]/cosTheta)) isAbsorbedAlphaMinus = true;
      } 
    }

    // ACCEPTANCE: propagate down to lower slab
    bool isInAccept = false;
    Double_t xDown = xUp + spacingZ*cos(phi)*sqrt(1-cosTheta*cosTheta)/cosTheta;
    Double_t yDown = yUp + spacingZ*sin(phi)*sqrt(1-cosTheta*cosTheta)/cosTheta;
  
    if (fabs(xDown) < lengthX/2. && fabs(yDown) < widthY/2.) {
      isInAccept = true;    // is in acceptance
      npassedAcc++;
    }

    if (!isAbsorbedAlpha) {
      // define tracks
      if (ndrawn < 1000) {
	float xDraw = xDown;           float yDraw = yDown;  
	float zDraw1 = -spacingZ/2.;    float zDraw2 = -spacingZ/2.; 
	if (!isInAccept) {
	  if (xUp > xDown) xDraw = TMath::Max(xDown,(double)-0.85*lengthX);
	  else xDraw = TMath::Min(xDown,(double)0.85*lengthX);
	  if (yUp > yDown) yDraw = TMath::Max(yDown,(double)-0.85*widthY);
	  else yDraw = TMath::Min(yDown,(double)0.85*widthY);
	  zDraw1 = spacingZ*(xDraw-xDown)/(xUp-xDown) - spacingZ/2.;
	  zDraw2 = spacingZ*(yDraw-yDown)/(yUp-yDown) - spacingZ/2.;
	  icolor = 3;
	} else icolor = 2;
	longTrack[ndrawn] = new TLine(xUp,spacingZ/2.,xDraw,zDraw1);
	longTrack[ndrawn]->SetLineWidth(1);
	longTrack[ndrawn]->SetLineColor(icolor);
	trasvTrack[ndrawn] = new TLine(yUp,spacingZ/2.,yDraw,zDraw2);
	trasvTrack[ndrawn]->SetLineWidth(1);
	trasvTrack[ndrawn]->SetLineColor(icolor);
	ndrawn++;
      }
    }
      
    if (isInAccept && !isAbsorbedAlpha) {
      npassedAccAtt++;
      hxDown->Fill(xDown);
      hyDown->Fill(yDown);
      // compute Pof
      double pof = spacingZ/cosTheta;    
      meanPath += pof;         meanPathSq += pof*pof;
      hpath->Fill(pof); 
      hphiAcc->Fill(phi);
      hcosThetaAcc->Fill(cosTheta);
    }

    if (isInAccept && !isAbsorbedAlphaPlus) npassedAccAttPlus++;
    if (isInAccept && !isAbsorbedAlphaMinus) npassedAccAttMinus++;
  }

  // Results: acceptance/attenuation and mean path-of-flight 
  double acceptance = (double)npassedAcc/(double)events;
  double accErr = sqrt(acceptance*(1-acceptance)/(double)events);
  double accAtt = (double)npassedAccAtt/(double)events;
  double accAttPlus = (double)npassedAccAttPlus/(double)events;
  double accAttMinus = (double)npassedAccAttMinus/(double)events;
  double accAttErr = (fabs(accAtt-accAttPlus)+fabs(accAtt-accAttMinus))/2.;   
  meanPath /= (double)npassedAccAtt;   meanPathSq /= (double)npassedAccAtt;
  double meanPathRMS = sqrt(meanPathSq - meanPath*meanPath);
  double meanPathErr = sqrt((meanPathSq - meanPath*meanPath)/double(npassedAccAtt));

  std::cout << "Test: G(stat. err.) = " << std::setprecision(4) << acceptance << " +/- " << std::setprecision(2) << accErr << std::endl;
  std::cout << "Test: eta = " << std::setprecision(4) << exp(-alpha*thickTotal) << std::endl;
  std::cout << "Test: G*eta(factorized) = " << std::setprecision(4) << acceptance*exp(-alpha*thickTotal) << std::endl;
  std::cout << "FINAL: G*eta = " << std::setprecision(4) << accAtt << " +/- " << std::setprecision(3) << accAttErr << std::endl;
 
  std::cout << "FINAL: mean Path-Of-Flight = (" << std::setprecision(5) << meanPath << " +/- " << std::setprecision(2) << meanPathErr << ") cm - RMS: " << std::setprecision(4) << meanPathRMS << " cm" << endl;

  char fileString[200];

  if (drawTracks) {

    // Draw tracks and scintillators
 
    // Define axes
    double dummyGraphX[4] = {(double)-0.8*lengthX, (double)-0.8*lengthX, (double)0.8*lengthX, (double)0.8*lengthX};
    double dummyGraphY[4] = {(double)-0.8*widthY, (double)-0.8*widthY, (double)0.8*widthY, (double)0.8*widthY};
    double dummyGraphZ[4] = {(double)-0.6*spacingZ, (double)0.6*spacingZ, (double)-0.6*spacingZ, (double)0.6*spacingZ};
    TGraph* dummyGraphLong = new TGraph(4,dummyGraphX,dummyGraphZ);
    dummyGraphLong->SetMarkerColor(0);
    dummyGraphLong->SetTitle("Longitudinal view");
    TGraph* dummyGraphTrasv = new TGraph(4,dummyGraphY,dummyGraphZ);
    dummyGraphTrasv->SetMarkerColor(0);
    dummyGraphTrasv->SetTitle("Transverse view");
    
    // best draw settings (maybe?)
    unsigned int lineDraw = (ndrawn > 1000) ? 1000 : ndrawn;
    int canvHor = 800;        int canvVer = 800.;
    if (lengthX > 4.*spacingZ) {canvHor = 1000; canvVer = 200.;}
    
    // longitudinal view
    TCanvas cl("cl","Longitudinal view",10,10,canvHor,canvVer);
    cl.SetGrid();
    dummyGraphLong->Draw("AP");
    TLine scintUpL(-lengthX/2.,spacingZ/2.,lengthX/2.,spacingZ/2.);
    scintUpL.SetLineWidth(8);   scintUpL.SetLineColor(4);   
    TLine scintDownL(-lengthX/2.,-spacingZ/2.,lengthX/2.,-spacingZ/2.);
    scintDownL.SetLineWidth(8);   scintDownL.SetLineColor(4);   
    scintUpL.Draw("same"); 
    scintDownL.Draw("same");
    for (unsigned int il = 0; il < lineDraw; il++) {
      longTrack[il]->Draw("same");
    }
    sprintf(fileString,"longitView_l%d_w%d_sp%d_alpha%d.gif",(int)lengthX,(int)widthY,(int)spacingZ,(int)(alpha*10.));
    cl.SaveAs(fileString);
    cl.Print();
    
    // best draw settings (maybe?)
    canvHor = 800;            canvVer = 800.;
    if (widthY > 4.*spacingZ) {canvHor = 1000; canvVer = 200.;}
    
    // transverse view
    TCanvas ct("ct","Transverse view",20,20,canvHor,canvVer);
    ct.SetGrid();
    dummyGraphTrasv->Draw("AP");
    TLine scintUpT(-widthY/2.,spacingZ/2.,widthY/2.,spacingZ/2.);
    scintUpT.SetLineWidth(8);   scintUpT.SetLineColor(4);   
    TLine scintDownT(-widthY/2.,-spacingZ/2.,widthY/2.,-spacingZ/2.);
    scintDownT.SetLineWidth(8);   scintDownT.SetLineColor(4);   
    scintUpT.Draw("same"); 
    scintDownT.Draw("same");
    for (unsigned int il = 0; il < lineDraw; il++) {
      trasvTrack[il]->Draw("same");
    }
    sprintf(fileString,"transvView_l%d_w%d_sp%d_alpha%d.gif",(int)lengthX,(int)widthY,(int)spacingZ,(int)(alpha*10.));
    ct.SaveAs(fileString);
    ct.Print();
  }

  // output histograms
  sprintf(fileString,"outputHists_l%d_w%d_sp%d_alpha%d.root",(int)lengthX,(int)widthY,(int)spacingZ,(int)(alpha*10.));
  TFile f(fileString,"RECREATE");
  hxUp->Write();
  hyUp->Write();
  hxDown->Write();
  hyDown->Write();
  hphi->Write();
  hcosTheta->Write();
  hphiAcc->Write();
  hcosThetaAcc->Write();
  hpath->Write();
  f.Close();

}
