#include <TH1D.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TApplication.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <iostream>
#include "myTRandom.h"
#include "Geometry.h"

using namespace std;

namespace par{
  // Parametri dello scintillatore (NaI)
  double me      = 0.511;
  double rho     = 3.67; // g/cm^3
  double cm      = 0.01;
  // Volume scintillatore cilindro di raggio r con asse z, una base a z=0, l'altra a z=dz
  double r       = 0.02;
  double dz      = 0.025;
  double dist    = 0.30;
  double coslim  = dist/sqrt(dist*dist+r*r);
  // Fotone
  // Energia
  double en      = 0.6617;
  // Numero di eventi
  int    nev     = 10000;
  // Massimo numero di fotoni graficati
  int    ngr     = 100;
  // ROOT & utilities
  myTRandom   *rnd;
  Geometry    *geo;
};


TVector3 GetDirection(TVector3 d, double theta, double phi){

  TVector3 dp;
  if (d.Z()!=1){
    dp.SetX(d.X()*cos(theta)+sin(theta)/sqrt(1-pow(d.Z(),2))*(d.X()*d.Z()*cos(phi)-d.Y()*sin(phi)));
    dp.SetY(d.Y()*cos(theta)+sin(theta)/sqrt(1-pow(d.Z(),2))*(d.Y()*d.Z()*cos(phi)+d.X()*sin(phi)));
    dp.SetZ(d.Z()*cos(theta)-sqrt(1-pow(d.Z(),2))*sin(theta)*cos(phi));
  } else if (d.Z()==1) {
    dp.SetX(sin(theta)*cos(phi));
    dp.SetY(sin(theta)*sin(phi));
    dp.SetZ(cos(theta));
  } else if (d.Z()==-1) {
    dp.SetX(sin(theta)*cos(phi+TMath::Pi()));
    dp.SetY(sin(theta)*sin(phi+TMath::Pi()));
    dp.SetZ(-cos(theta));
  }
  return dp;

}

TGraph grPhoto("NaIPhoto.data");
TGraph grCompton("NaICompton.data");

double CalcEnergyDeposition(TVector3 x0, TVector3 d, double E, double& w,string modus){
  // Return the energy deposition

  double Edep = 0;
  do {
    
    //TODO compute mu
    double muCompton = grCompton.Eval(E)*par::rho/par::cm;
    double muPhoto   = grPhoto.Eval(E)*par::rho/par::cm;
    double muTot     = muCompton + muPhoto;
    double lambda    = 1/muTot;
    
    //TODO compute free path
    double s = -lambda*log(par::rnd->Rndm());
    // calcolo posizione 
    TVector3 x  = x0 + d*s;
    x0 = x;
    
    // check that the particle is still in the detector
    if (!par::geo->Contains(x))
      break;
    //Update graphics
    par::geo->Update(x);

    bool lfoto = (par::rnd->Rndm() < muPhoto/muTot); // True if Photoelectric
    bool lbias = (modus=="bias+");
    if (lbias){
      w = w*(1-muPhoto/muTot);
      double pkill=0.35;
      if (par::rnd->Rndm()<pkill){
	break;
      } else {
	w = w/(1-pkill);
      }
      lfoto = false;  
    }

    if (lfoto){
      Edep += E;
      E = 0;
    } else {
      double phi = par::rnd->Rndm()*TMath::TwoPi();
      double Eprim, theta;
      par::rnd->KleinNishina(E,Eprim,theta);
      d     = GetDirection(d,theta,phi);
      Edep += E-Eprim;
      E     = Eprim;
    }

  } while (E!=0);

  return Edep;

}

void SimSci(string modus){
  
  par::rnd = new myTRandom;
  par::rnd->SetSeed(time(NULL));
  
  par::geo = new Geometry("CYLINDER",par::ngr);
  par::geo->SetDim(par::r,par::dz);
  
  par::geo->Draw();  

  TH1D *hE = new TH1D("hE","",100,0,1.0);
  TStopwatch tstop;

  for (int iev=0;iev<par::nev;iev++){
    TVector3 x0,d;

    double w = 1;
    double theta=0,phi=0,alpha=0,costh=1;
    if (modus=="detsimple"){
      theta = 0; 
      phi   = 0;
    } else if (modus=="det"){
      phi = par::rnd->Rndm()*TMath::TwoPi();
      costh = 2*par::rnd->Rndm()-1;
      if (costh<par::coslim) continue;
    } else if (modus=="bias" || modus=="bias+"){
      costh = 1-(1-par::coslim)*par::rnd->Rndm();
      phi = par::rnd->Rndm()*TMath::TwoPi();
      w = (1-par::coslim)/2;
    }
    theta = acos(costh);
    alpha = par::dist/costh;
    d  = TVector3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
    x0 = TVector3(d.X()*alpha,d.Y()*alpha,0);
    
    par::geo->Event(x0);
    double Edep  = CalcEnergyDeposition(x0,d,par::en,w,modus);
    if (Edep!=0){
      double Emeas = par::rnd->Gaus(Edep,sqrt(0.02*0.02+0.04*0.04/Edep)*Edep);
      hE->Fill(Emeas,w);
    }
  }

  double time = tstop.CpuTime();
  double eval;
  double val = hE->IntegralAndError(hE->GetXaxis()->FindBin(0.35),hE->GetXaxis()->FindBin(0.4),eval);
  val  = val/hE->Integral();
  eval = eval/hE->Integral();
  std::cout << "FOM " << 1/(time*pow(eval/val,2)) << std::endl;
  
  TCanvas *c = new TCanvas("h","Energy",620,10,320,320);
  c->Draw();

  gStyle->SetOptStat("111111");
  hE->Draw("E");
  
  string file="output_"; file+=modus; file+=".root";
  TFile f(file.c_str(),"recreate");
  hE->Write();
  f.Close();

}

#ifndef __CINT__
int main(int argc,char* argv[]){
  if (argc!=2){
    cout << "Devi specificare un parametro" << endl;
    cout << " detsimple  - simulazione dettagliata semplificata (fotone lungo z) " << endl;
    cout << " det        - simulazione dettagliata con direzione corretta " << endl;
    cout << " bias       - simulazione con bias della sorgente " << endl;
    cout << " bias+      - simulazione con bias (peso statistico + roulette russa) " << endl;
    return 1;
  }
  TApplication app("app",0,NULL);

  SimSci(argv[1]);

  gPad->Update();
  app.Run(true);

  return 0;
}
#endif
