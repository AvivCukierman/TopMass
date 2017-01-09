#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TGenPhaseSpace.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
//#include "TRoot.h"
#include "TMarker.h"
#include "TH2F.h"

void example(double number=10000, double mtop = 175.){

    double pTmax = 0;
    //double mtop = 175.; // top quark mass in GeV
    double mW = 80.;// W boson mass in GeV
    double mB = 5.; // b quark mass in GeV.
    double had = 0.0; //hadronization parameter
    double loss = 0.0; //hadronization parameter
    double swap = 0.0; //combinatoric noise

    TH1F *mjjj = new TH1F ("","",40,0,400);

    double pi = 4.*atan(1.);
    TLorentzVector jet = TLorentzVector(0,0,0,mtop);
    vector<double> masses;
    masses.push_back(mW);
    masses.push_back(mB);
    TGenPhaseSpace event;
    event.SetDecay(jet, masses.size(), &masses[0]); 
    TGenPhaseSpace eventb;
    eventb.SetDecay(jet, masses.size(), &masses[0]); 
    vector<double> masses2;
    masses2.push_back(0);
    masses2.push_back(0);

    TRandom3 rand = TRandom3(0);

    char filename[50];
    sprintf(filename,"out%d.root",mtop);
    TFile *out = new TFile(filename,"RECREATE");
    std::vector<float> mjjj_vf;
    TTree *tree = new TTree("oTree","");
    float mjjj_f;
    tree->Branch("mjjj",&mjjj_f);
    float j1_pt; tree->Branch("j1_pt",&j1_pt);
    float j1_eta; tree->Branch("j1_eta",&j1_eta);
    float j1_phi; tree->Branch("j1_phi",&j1_phi);
    float j1_m; tree->Branch("j1_m",&j1_m);
    float j2_pt; tree->Branch("j2_pt",&j2_pt);
    float j2_eta; tree->Branch("j2_eta",&j2_eta);
    float j2_phi; tree->Branch("j2_phi",&j2_phi);
    float j2_m; tree->Branch("j2_m",&j2_m);
    float b_pt; tree->Branch("b_pt",&b_pt);
    float b_eta; tree->Branch("b_eta",&b_eta);
    float b_phi; tree->Branch("b_phi",&b_phi);
    float b_m; tree->Branch("b_m",&b_m);
    float W_pt; tree->Branch("W_pt",&W_pt);
    float W_eta; tree->Branch("W_eta",&W_eta);
    float W_phi; tree->Branch("W_phi",&W_phi);
    float W_m; tree->Branch("W_m",&W_m);
    float t_pt; tree->Branch("t_pt",&t_pt);
    float t_eta; tree->Branch("t_eta",&t_eta);
    float t_phi; tree->Branch("t_phi",&t_phi);
    float t_m; tree->Branch("t_m",&t_m);

    for (int i = 1; i<= number; i++){
        if(!(i%1000)) std::cout << i << std::endl;

        double phi = rand.Uniform(0.,2*pi);
        double theta = rand.Uniform(0.,pi);
        double sinTheta = sin(theta);
        double cosTheta = cos(theta);
           
        double q=pTmax/(2*mtop*sin(theta));
        double beta = q/sqrt(1+q*q);
        TVector3 betavec = TVector3(beta*sinTheta*cos(phi), beta*sinTheta*sin(phi), beta*cosTheta);
        event.Generate();
        TLorentzVector *W1=event.GetDecay(0); 
        TLorentzVector *b1=event.GetDecay(1); 
        W1->Boost(betavec);
        b1->Boost(betavec);
        TLorentzVector top1 = *W1 + *b1;

        int aaa = __LINE__; //for some reason it crashes if this isn't here

        eventb.Generate();
        TLorentzVector *W2=eventb.GetDecay(0); 
        TLorentzVector *b2=eventb.GetDecay(1); 
        W2->Boost(betavec);
        b2->Boost(betavec);
        TLorentzVector top2 = *W2 + *b2;

        TGenPhaseSpace event2;
        event2.SetDecay(*W1, masses2.size(), &masses2[0]); 
        event2.Generate();
        TLorentzVector *ele=event2.GetDecay(0); 
        TLorentzVector *eleneut=event2.GetDecay(1); 

        TGenPhaseSpace event2b;
        event2b.SetDecay(*W2, masses2.size(), &masses2[0]); 
        event2b.Generate();
        TLorentzVector *j1=event2b.GetDecay(0); 
        TLorentzVector *j2=event2b.GetDecay(1); 
        //Parton shower, hadronization, jet clustering all in one.
        j1->SetPtEtaPhiM(j1->Pt()*rand.Gaus(1-loss,had),j1->Eta(),j1->Phi(),0.);
        j2->SetPtEtaPhiM(j2->Pt()*rand.Gaus(1-loss,had),j2->Eta(),j2->Phi(),0.);
        b1->SetPtEtaPhiM(b1->Pt()*rand.Gaus(1-loss,had),b1->Eta(),b1->Phi(),mB);
        b2->SetPtEtaPhiM(b2->Pt()*rand.Gaus(1-loss,had),b2->Eta(),b2->Phi(),mB);

        TLorentzVector metT = TLorentzVector((*eleneut).Px(),(*eleneut).Py(),0.,(*eleneut).Pt());
        TLorentzVector lepT = TLorentzVector(ele->Px(),ele->Py(),0.,ele->Pt()); 

        if (rand.Uniform(0,1) > swap) mjjj_f = ((*j1)+(*j2)+(*b2)).M();
        else mjjj_f = ((*j1)+(*j2)+(*b1)).M();

        j1_pt = j1->Pt();
        j1_eta = j1->Eta();
        j1_phi = j1->Phi();
        j1_m = j1->M();
        j2_pt = j2->Pt();
        j2_eta = j2->Eta();
        j2_phi = j2->Phi();
        j2_m = j2->M();
        b_pt = b2->Pt();
        b_eta = b2->Eta();
        b_phi = b2->Phi();
        b_m = b2->M();
        W_pt = W2->Pt();
        W_eta = W2->Eta();
        W_phi = W2->Phi();
        W_m = W2->M();
        t_pt = top1.Pt();
        t_eta = top1.Eta();
        t_phi = top1.Phi();
        t_m = top1.M();

        tree->Fill();
    }
    out->Write();

    /*
    TCanvas *c1 = new TCanvas("","",750,750);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);

    mjjj->GetXaxis()->SetTitleOffset(1.2);
    mjjj->GetXaxis()->SetTitle("m_{jjj} [GeV]");
    mjjj->GetYaxis()->SetTitleOffset(1.8);
    mjjj->GetYaxis()->SetNdivisions(505);
    mjjj->GetYaxis()->SetTitle("Normalized to Unity");
    mjjj->SetLineColor(2);

    mjjj->Draw();

    TLatex l  = TLatex();
    l.SetNDC();
    l.SetTextColor(1);
    l.DrawLatex(0.15,0.92,"#bf{#scale[0.6]{Scalar Phase Space Only}}");
    c1->Print("mjjj.pdf");*/

}


