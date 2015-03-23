//
#include <iostream>
#include <stdio.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TMatrixT.h"

using namespace std;

const static int NFIT = 25;

// To allow access to data file from all functions
static FILE * dFile;

// Arrays of matrices used for one collection of NFIT events
static TMatrixT<double> AA[NFIT];
static TMatrixT<double> CC[NFIT];

// Fetch single event from data file
void getEvent(TLorentzVector * pp);

// Fill A & C matrices
void fillMatrices(TLorentzVector * pp, TMatrixT<double> &A, TMatrixT<double> &C);

// \chi^2 function
void fcn(int &npar, double* gin, double &f, double* par, int iflag);

int fit( ){
  
  // Reset ROOT
  gROOT->Reset();
  
  // As the man says...
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  
  // Canvas
  TCanvas * c1 = new TCanvas("c1","",400,400);
  TCanvas * c2 = new TCanvas("c2","",800,400);

  // Histograms
  TH2D * h_fit_MN_MZ = new TH2D("h_fit_MN_MZ","",300,0,300,300,400,650);
  h_fit_MN_MZ->SetMarkerStyle(4);
  h_fit_MN_MZ->SetMarkerColor(kBlue);
  h_fit_MN_MZ->SetMarkerSize(0.25);
  TH2D * h_fit_MX_MZ = new TH2D("h_fit_MX_MZ","",300,0,300,300,400,650);
  h_fit_MX_MZ->SetMarkerStyle(4);
  h_fit_MX_MZ->SetMarkerColor(kRed);
  h_fit_MX_MZ->SetMarkerSize(0.25);
  TH2D * h_fit_MY_MZ = new TH2D("h_fit_MY_MZ","",300,0,300,300,400,650);
  h_fit_MY_MZ->SetMarkerStyle(4);
  h_fit_MY_MZ->SetMarkerColor(kGreen);
  h_fit_MY_MZ->SetMarkerSize(0.25);
  TH1D * hmsq = new TH1D("hmsq","",300,0,650);
  TH1D * hmchi02 = new TH1D("hmchi02","",300,0,650);
  hmchi02->SetLineColor(kGreen);
  
  // Decorations
  TLine * lMZ = new TLine(0,565,300,565);
  //lsquark->SetLineWidth(3);
  lMZ->SetLineStyle(2);
  TLine * lMN = new TLine(97,400,97,650);
  lMN->SetLineStyle(2);
  lMN->SetLineColor(kBlue);
  TLine * lMX = new TLine(143,400,143,650);
  lMX->SetLineStyle(2);
  lMX->SetLineColor(kRed);
  TLine * lMY = new TLine(177,400,177,650);
  lMY->SetLineStyle(2);
  lMY->SetLineColor(kGreen);
  
  // Open file
  dFile = fopen("/home/jorgenem/git-repos/master/simulation/events/Pythia_cascade_events_no_ISR_or_FSR_20150120_only_opposite_flavour_leptons.dat","r");
  // dFile = fopen("/home/jorgenem/git-repos/master/simulation/events/simple_2500_events_gauss_and_exp_mass_smearing.dat","r");
  // dFile = fopen("/home/jorgenem/git-repos/master/simulation/events/Pythia_cascade_10000_events_everything_turned_on_20150210_only_opposite_flavour_leptons.dat","r");
  
  TLorentzVector pp[9];
  for(int iFits = 0; iFits < 100; iFits++){

    // Get events
    for(int iEvent = 0; iEvent < NFIT; iEvent++){
      getEvent(pp);
      AA[iEvent].ResizeTo(8,8);
      CC[iEvent].ResizeTo(8,1);
      fillMatrices(pp, AA[iEvent], CC[iEvent]);
    
      // Store interesting information for events
      TLorentzVector chi021 = pp[2]+pp[3]+pp[4];
      TLorentzVector chi022 = pp[6]+pp[7]+pp[8];
      TLorentzVector sq1 = chi021+pp[1];
      TLorentzVector sq2 = chi022+pp[5];
      hmchi02->Fill(chi021.M());
      hmchi02->Fill(chi022.M());
      hmsq->Fill(sq1.M());
      hmsq->Fill(sq2.M());
    }
  
    // Initialize TMinuit with 4 parameters
    TMinuit * gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);
  
    double arglist[10];
    int ierflg = 0;
  
    arglist[0] = 0.00000001;
    // arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg); // Set error definition
    arglist[0] = 2;
    gMinuit->mnexcm("SET PRI", arglist, 1, ierflg); // Set verbosity
  
    // Set starting values and step sizes for parameters
    double vstart[4] = {97.0, 143., 177., 565.};
    double step[4] = {0.1,0.1,0.1,0.1};
    gMinuit->mnparm(0, "MN", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "MX", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "MY", vstart[2], step[2], 0, 0, ierflg);
    gMinuit->mnparm(3, "MZ", vstart[3], step[3], 0, 0, ierflg);
  
    // Now ready for minimization step
    // arglist[0] = 5000;
    // arglist[1] = 1.;
    arglist[0] = 0;
    arglist[1] = 0;
    gMinuit->mnexcm("SIMPLEX", arglist ,1, ierflg);
    //gMinuit->mnexcm("MIGRAD", arglist ,2, ierflg);
    cout << "Fit status: " << ierflg << endl;
    
    // Extract parameter values at minima if we have convergence
    if(ierflg == 0){
      double vfinish[4], efinish[4];
      for(int i = 0; i < 4; i++){
        gMinuit->GetParameter(i, vfinish[i], efinish[i]);
      }
      h_fit_MN_MZ->Fill(vfinish[0],vfinish[3]);
      h_fit_MX_MZ->Fill(vfinish[1],vfinish[3]);
      h_fit_MY_MZ->Fill(vfinish[2],vfinish[3]);
    }
  }

  // // Test of function value
  // double f;
  // int n;
  // double par[4] = {97,144,180,565};
  // fcn(n, 0, f, par, 0);
  // cout << "xi^2 = " << f << endl;
  
  
  // Plot
  c1->cd();
  h_fit_MN_MZ->Draw();
  h_fit_MX_MZ->Draw("same");
  h_fit_MY_MZ->Draw("same");
  lMN->Draw("same");
  lMX->Draw("same");
  lMY->Draw("same");
  lMZ->Draw("same");
  
  
  c2->cd();
  hmsq->Draw();
  hmchi02->Draw("same");
  
  // Close file
  fclose(dFile);
  
  return 0;
}

void fcn(int &npar, double* gin, double &f, double* par, int iflag){
  
  // Define masses
  double MN2 = pow(par[0],2);
  double MX2 = pow(par[1],2);
  double MY2 = pow(par[2],2);
  double MZ2 = pow(par[3],2);
  
  // Define matrices
  TMatrixT<double> A(8,8), B(8,8), D(8,8), E(8,1), M(8,1), P(8,1);
  
  // Fill sparse B
  B(0,0) = -1; B(0,1) = 1; B(1,1) = -1; B(1,2) = 1; B(2,2) = -1; B(2,3) = 1;
  B(4,4) = -1; B(4,5) = 1; B(5,5) = -1; B(5,6) = 1; B(6,6) = -1; B(6,7) = 1;
  //B.Print();
  
  // chi2
  double chi2 = 0;
  
  //cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << endl;

  // Get chi2
  // Loop over events
  for(int i = 0; i < NFIT; i++){
    
    // Fill M
    M(0,0) = MZ2; M(1,0) = MY2; M(2,0) = MX2; M(3,0) = MN2;
    M(4,0) = MZ2; M(5,0) = MY2; M(6,0) = MX2; M(7,0) = MN2;
    //M.Print();
    
    // Invert A
    A = AA[i];
    A.Invert();
    //A.Print();
    
    D = A*B;
    //D.Print();
    
    E = A*CC[i];
    //E.Print();
    
    P = D*M+E;
    //P.Print();
    
    P.Sqr();

    chi2 += pow(P(3,0)-P(0,0)-P(1,0)-P(2,0)-MN2,2);
    chi2 += pow(P(7,0)-P(4,0)-P(5,0)-P(6,0)-MN2,2);
    
  }

  // chi2 to return
  // cout << chi2/pow(100.,4) << endl;
  f = chi2/pow(100.,4);
}

void getEvent(TLorentzVector * pp){
  
  // Temp storage for particles when reading file
  int pdg;
  float px, py, pz, EE;
  
  // Read file if not at end
  if( !feof(dFile) ){
    // Read one event
    fscanf(dFile, "%*[^\n]\n", NULL);    // Trick to skip line
    for(int iParticle = 1; iParticle < 9; iParticle++) {
      // Read particle information
      fscanf (dFile, "%i %f %f %f %f %*f", &pdg, &px, &py, &pz, &EE, NULL);
      //cout << px << " " << py << " " << pz << " " << EE << endl;
      pp[iParticle].SetPxPyPzE(px, py, pz, EE);
      
    }
    fscanf(dFile, "%*[^\n]\n", NULL);    // Need to get to end of line *sigh*
    
  }
  else { cout << "Reached end of file!" << endl; }
  
}

void fillMatrices(TLorentzVector * pp, TMatrixT<double> &A, TMatrixT<double> &C){
  // Fill A matrix
  A(0,0) = pp[1].Px(); A(0,1) = pp[1].Py(); A(0,2) = pp[1].Pz(); A(0,3) = -pp[1].E();
  A(1,0) = pp[2].Px(); A(1,1) = pp[2].Py(); A(1,2) = pp[2].Pz(); A(1,3) = -pp[2].E();
  A(2,0) = pp[3].Px(); A(2,1) = pp[3].Py(); A(2,2) = pp[3].Pz(); A(2,3) = -pp[3].E();
  A(3,0) = 0.5;        A(3,4) = 0.5;
  A(4,4) = pp[5].Px(); A(4,5) = pp[5].Py(); A(4,6) = pp[5].Pz(); A(4,7) = -pp[5].E();
  A(5,4) = pp[6].Px(); A(5,5) = pp[6].Py(); A(5,6) = pp[6].Pz(); A(5,7) = -pp[6].E();
  A(6,4) = pp[7].Px(); A(6,5) = pp[7].Py(); A(6,6) = pp[7].Pz(); A(6,7) = -pp[7].E();
  A(7,1) = 0.5;        A(7,5) = 0.5;
  //A.Print();
  A = A+A;
  //A.Print();
  
  // Fill C vector
  C(0,0) = 2*pp[1]*pp[2]+2*pp[1]*pp[3]+pp[1].M2();
  C(1,0) = 2*pp[2]*pp[3]+pp[2].M2();
  C(2,0) = pp[3].M2();
  C(3,0) = pp[4].Px()+pp[8].Px(); // pxmiss
  C(4,0) = 2*pp[5]*pp[6]+2*pp[5]*pp[7]+pp[5].M2();
  C(5,0) = 2*pp[6]*pp[7]+pp[6].M2();
  C(6,0) = pp[7].M2();
  C(7,0) = pp[4].Py()+pp[8].Py(); // pymiss
  //C.Print();
}