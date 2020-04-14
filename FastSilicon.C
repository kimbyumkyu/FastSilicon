#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TMath.h"
#include <iostream>
#include "TRandom.h"
#include "TGraph.h"
#include "TFile.h"



const Double_t VESAT = 1.1e5;
const Double_t VHSAT = 9.5e4;
const Double_t MU_E0 = 0.1350;
const Double_t MU_H0 = 0.0480; // electron mobility in m2/Vs
const Double_t ec = 1.60217662e-19;
const Double_t kB = 1.38065e-23;  // Boltzmann constant (J/K)
const Double_t m0 = 9.109383e-31; //electron mass (kg)
const Double_t me = 0.26 * m0;    // mass of electron see: DOI: 10.1109/TNS.2009.2021426
const Double_t mh = 0.36 * m0;
Double_t Temp = 300;
const Double_t vacuumpermitttivity = 8.854187e-12; // F/m
const Double_t siliconrelativepermittivity = 11.68;
//const Double_t pi4per = TMath::Pi() * 4 * vacuumpermitttivity * siliconrelativepermittivity; //F/m
const Double_t pi4per = vacuumpermitttivity * siliconrelativepermittivity; //F/m
Double_t spacedopingdensity = 1.9e12 * 1e6;                         // e.a/m^3
Double_t spacedopingdensityp = 1e15 * 1e6;                         // e.a/m^3
//typedef std::vector<Double_t> Double1D;
//typedef std::vector<Double1D> Double2D;
//typedef std::vector<Double3D> Double3D;
template <typename H>
void setstyle(TH3F *h);
void setstyle(TH2F *h);
void setstyle(TH1F *h);
void setpad(TVirtualPad *pad);
Double_t xdim = 300;
Double_t ydim = 300;
Double_t zdim = 300;
Double_t reversebiasvoltage = 140;
const Int_t scale = Int_t(xdim / zdim);
TH3F *base = nullptr;
TH3F *pbase = nullptr;

TH3F *ndensity = nullptr;
TH3F *hdensity = nullptr;
TH3F *EFX = nullptr;
TH3F *EFY = nullptr;
TH3F *EFZ = nullptr;

TGraph *gcurrent = nullptr;
TGraph *gncurrent = nullptr;
TGraph *ghcurrent = nullptr;
//TCanvas *c1 = nullptr;
//TCanvas *c2 = nullptr;
//TCanvas *c3 = nullptr;
Int_t globalmeshcounter = 0;


Double_t theta = 45;
vector<Double_t> nx, ny, nz;
vector<Double_t> hx, hy, hz;
Double_t MIP = 75; // /micrometer
Double_t accurate = 1e-3;
Double_t ptime = 0;

void setstyle(TGraph *h);

class MyMainFrame {
   RQ_OBJECT("MyMainFrame")
private:
   TGMainFrame         *fMainFrame;
   TGCompositeFrame    *fSubFrame;
   TRootEmbeddedCanvas *fMainCanvas;
   TGNumberEntry *fDetectorWidthValue;
   TGNumberEntry *fDetectorLengthValue;
   TGNumberEntry *fDetectorDepthValue;
   TGNumberEntry *fNtypeDopingDensityValue;
   TGNumberEntry *fPtypeDoingDensityValue;

   TRootEmbeddedCanvas *fCanvas1;
   TRootEmbeddedCanvas *fCanvas2;
   TRootEmbeddedCanvas *fCanvas3;
   TRootEmbeddedCanvas *fCanvas4;
   TRootEmbeddedCanvas *fCanvas5;
   TRootEmbeddedCanvas *fCanvas6;
   TRootEmbeddedCanvas *fCanvas7;
   TRootEmbeddedCanvas *fCanvas8;
   TRootEmbeddedCanvas *fCanvas9;


   TGNumberEntry *fReverseBiasVoltageValue;
   TGNumberEntry *fTemperatureValue;
   TGNumberEntry *fBFieldValue;

   TGNumberEntry *fMIPValue;
   TGNumberEntry *fIncedentAngleValue;

   public:
   MyMainFrame();
   virtual ~MyMainFrame();
   void StartAnalysis();
   void SolvePoisson(const int nmesh, const Double_t time, std::vector<TCanvas *> C);
   void InjectParticle(Double_t theta);
   void SolveCurrent(Double_t time, Double_t nmesh, std::vector<TCanvas *> C);
   void DoAnalysis();
};
MyMainFrame::~MyMainFrame() {
   // Clean up used widgets: frames, buttons, layout hints
   fMainFrame->Cleanup();
   delete fMainFrame;
}

MyMainFrame::MyMainFrame(){
   // main frame
   fMainFrame = new TGMainFrame(gClient->GetRoot(), 10, 10, kMainFrame | kVerticalFrame);
   fMainFrame->SetName("Fast Silicon");
   fMainFrame->SetLayoutBroken(kTRUE);

   // composite frame
   fSubFrame = new TGCompositeFrame(fMainFrame,2186,1400,kVerticalFrame);
   fSubFrame->SetName("fSubFrame");
   fSubFrame->SetLayoutBroken(kTRUE);

   // embedded canvas
   fCanvas1 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas1->SetName("fCanvas1");
   Int_t wfCanvas1 = fCanvas1->GetCanvasWindowId();
   TCanvas *c177 = new TCanvas("c177", 10, 10, wfCanvas1);
   fCanvas1->AdoptCanvas(c177);
   fSubFrame->AddFrame(fCanvas1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas1->MoveResize(8,8,536,424);

   // embedded canvas
   fCanvas2 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas2->SetName("fCanvas2");
   Int_t wfCanvas2 = fCanvas2->GetCanvasWindowId();
   TCanvas *c178 = new TCanvas("c178", 10, 10, wfCanvas2);
   fCanvas2->AdoptCanvas(c178);
   fSubFrame->AddFrame(fCanvas2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas2->MoveResize(552,8,536,424);

   // embedded canvas
   fCanvas3 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas3->SetName("fCanvas3");
   Int_t wfCanvas3 = fCanvas3->GetCanvasWindowId();
   TCanvas *c179 = new TCanvas("c179", 10, 10, wfCanvas3);
   fCanvas3->AdoptCanvas(c179);
   fSubFrame->AddFrame(fCanvas3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas3->MoveResize(1096,8,536,424);

   // embedded canvas
   fCanvas4 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas4->SetName("fCanvas4");
   Int_t wfCanvas4 = fCanvas4->GetCanvasWindowId();
   TCanvas *c180 = new TCanvas("c180", 10, 10, wfCanvas4);
   fCanvas4->AdoptCanvas(c180);
   fSubFrame->AddFrame(fCanvas4, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas4->MoveResize(8,440,536,424);

   // embedded canvas
   fCanvas5 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas5->SetName("fCanvas5");
   Int_t wfCanvas5 = fCanvas5->GetCanvasWindowId();
   TCanvas *c181 = new TCanvas("c181", 10, 10, wfCanvas5);
   fCanvas5->AdoptCanvas(c181);
   fSubFrame->AddFrame(fCanvas5, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas5->MoveResize(552,440,536,424);

   // embedded canvas
   fCanvas6 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas6->SetName("fCanvas6");
   Int_t wfCanvas6 = fCanvas6->GetCanvasWindowId();
   TCanvas *c182 = new TCanvas("c182", 10, 10, wfCanvas6);
   fCanvas6->AdoptCanvas(c182);
   fSubFrame->AddFrame(fCanvas6, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas6->MoveResize(1096,440,536,424);

   // embedded canvas
   fCanvas7 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas7->SetName("fCanvas7");
   Int_t wfCanvas7 = fCanvas7->GetCanvasWindowId();
   TCanvas *c183 = new TCanvas("c183", 10, 10, wfCanvas7);
   fCanvas7->AdoptCanvas(c183);
   fSubFrame->AddFrame(fCanvas7, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas7->MoveResize(8,872,536,424);

   // embedded canvas
   fCanvas8 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas8->SetName("fCanvas8");
   Int_t wfCanvas8 = fCanvas8->GetCanvasWindowId();
   TCanvas *c184 = new TCanvas("c184", 10, 10, wfCanvas8);
   fCanvas8->AdoptCanvas(c184);
   fSubFrame->AddFrame(fCanvas8, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas8->MoveResize(552,872,536,424);

   // embedded canvas
   fCanvas9 = new TRootEmbeddedCanvas(0,fSubFrame,536,424,kSunkenFrame);
   fCanvas9->SetName("fCanvas9");
   Int_t wfCanvas9 = fCanvas9->GetCanvasWindowId();
   TCanvas *c185 = new TCanvas("c185", 10, 10, wfCanvas9);
   fCanvas9->AdoptCanvas(c185);
   fSubFrame->AddFrame(fCanvas9, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fCanvas9->MoveResize(1096,872,536,424);
   TGTextButton *fTextButton951 = new TGTextButton(fSubFrame,"Start Analysis",-1,TGTextButton::GetDefaultGC()(),TGTextButton::GetDefaultFontStruct(),kRaisedFrame);
   fTextButton951->SetTextJustify(36);
   fTextButton951->SetMargins(0,0,0,0);
   fTextButton951->SetWrapLength(-1);
   fTextButton951->Resize(312,56);
   fSubFrame->AddFrame(fTextButton951, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButton951->MoveResize(400,1320,312,56);
   fTextButton951->Connect("Clicked()", "MyMainFrame", this, "StartAnalysis()");


   TGTextButton *fTextButton952 = new TGTextButton(fSubFrame,"Exit","gApplication->Terminate(0)");
   fTextButton952->SetTextJustify(36);
   fTextButton952->SetMargins(0,0,0,0);
   fTextButton952->SetWrapLength(-1);
   fTextButton952->Resize(312,56);
   fSubFrame->AddFrame(fTextButton952, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTextButton952->MoveResize(936,1320,312,56);

   ULong_t ucolor;        // will reflect user color changes
   gClient->GetColorByName("#cccccc",ucolor);
   TGLabel *fLabel953 = new TGLabel(fSubFrame,"Detector Property",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fLabel953->SetTextJustify(36);
   fLabel953->SetMargins(0,0,0,0);
   fLabel953->SetWrapLength(-1);
   fSubFrame->AddFrame(fLabel953, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabel953->MoveResize(1656,24,512,56);

   gClient->GetColorByName("#85c2a3",ucolor);
   TGLabel *fLabel954 = new TGLabel(fSubFrame,"Detector Width (um)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fLabel954->SetTextJustify(36);
   fLabel954->SetMargins(0,0,0,0);
   fLabel954->SetWrapLength(-1);
   fSubFrame->AddFrame(fLabel954, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabel954->MoveResize(1656,104,344,56);
   fDetectorWidthValue = new TGNumberEntry(fSubFrame, (Double_t) 300,9,-1,(TGNumberFormat::EStyle) 5);
   fDetectorWidthValue->SetName("fDetectorWidthValue");
   fSubFrame->AddFrame(fDetectorWidthValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fDetectorWidthValue->MoveResize(2016,112,144,35);

   TGLabel *fLabel959 = new TGLabel(fSubFrame,"Detector Length (um)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fLabel959->SetTextJustify(36);
   fLabel959->SetMargins(0,0,0,0);
   fLabel959->SetWrapLength(-1);
   fSubFrame->AddFrame(fLabel959, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fLabel959->MoveResize(1656,168,344,56);
   fDetectorLengthValue = new TGNumberEntry(fSubFrame, (Double_t) 300,9,-1,(TGNumberFormat::EStyle) 5);
   fDetectorLengthValue->SetName("fDetectorLengthValue");
   fSubFrame->AddFrame(fDetectorLengthValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fDetectorLengthValue->MoveResize(2016,176,145,35);

   TGLabel *fDetectorDepth = new TGLabel(fSubFrame,"Detector Depth (um)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fDetectorDepth->SetTextJustify(36);
   fDetectorDepth->SetMargins(0,0,0,0);
   fDetectorDepth->SetWrapLength(-1);
   fSubFrame->AddFrame(fDetectorDepth, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fDetectorDepth->MoveResize(1656,232,344,56);
   fDetectorDepthValue = new TGNumberEntry(fSubFrame, (Double_t) 300,9,-1,(TGNumberFormat::EStyle) 5);
   fDetectorDepthValue->SetName("fDetectorDepthValue");
   fSubFrame->AddFrame(fDetectorDepthValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fDetectorDepthValue->MoveResize(2016,240,144,35);

   TGLabel *fNtypeDopingDensity = new TGLabel(fSubFrame,"N-type Doping Density (m-3)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fNtypeDopingDensity->SetTextJustify(36);
   fNtypeDopingDensity->SetMargins(0,0,0,0);
   fNtypeDopingDensity->SetWrapLength(-1);
   fSubFrame->AddFrame(fNtypeDopingDensity, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNtypeDopingDensity->MoveResize(1656,296,344,56);
   fNtypeDopingDensityValue = new TGNumberEntry(fSubFrame, (Double_t) 1.9e18,9,-1,(TGNumberFormat::EStyle) 5);
   fNtypeDopingDensityValue->SetName("fNtypeDopingDensityValue");
   fSubFrame->AddFrame(fNtypeDopingDensityValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fNtypeDopingDensityValue->MoveResize(2016,304,144,35);

   TGLabel *fPtypeDopingDensity = new TGLabel(fSubFrame,"P-type Doping Density (m-3)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fPtypeDopingDensity->SetTextJustify(36);
   fPtypeDopingDensity->SetMargins(0,0,0,0);
   fPtypeDopingDensity->SetWrapLength(-1);
   fSubFrame->AddFrame(fPtypeDopingDensity, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fPtypeDopingDensity->MoveResize(1656,360,344,56);
   fPtypeDoingDensityValue = new TGNumberEntry(fSubFrame, (Double_t) 1e21,9,-1,(TGNumberFormat::EStyle) 5);
   fPtypeDoingDensityValue->SetName("fPtypeDoingDensityValue");
   fSubFrame->AddFrame(fPtypeDoingDensityValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fPtypeDoingDensityValue->MoveResize(2016,368,144,35);

   gClient->GetColorByName("#cccccc",ucolor);
   TGLabel *fBoundaryCondition = new TGLabel(fSubFrame,"Boundary Condition",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fBoundaryCondition->SetTextJustify(36);
   fBoundaryCondition->SetMargins(0,0,0,0);
   fBoundaryCondition->SetWrapLength(-1);
   fSubFrame->AddFrame(fBoundaryCondition, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fBoundaryCondition->MoveResize(1656,448,512,56);

   gClient->GetColorByName("#85c2a3",ucolor);
   TGLabel *fReverseBiasVoltage = new TGLabel(fSubFrame,"Reverse-bias Voltage (V)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fReverseBiasVoltage->SetTextJustify(36);
   fReverseBiasVoltage->SetMargins(0,0,0,0);
   fReverseBiasVoltage->SetWrapLength(-1);
   fSubFrame->AddFrame(fReverseBiasVoltage, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fReverseBiasVoltage->MoveResize(1656,528,344,56);
   fReverseBiasVoltageValue = new TGNumberEntry(fSubFrame, (Double_t) 140,9,-1,(TGNumberFormat::EStyle) 5);
   fReverseBiasVoltageValue->SetName("fReverseBiasVoltageValue");
   fSubFrame->AddFrame(fReverseBiasVoltageValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fReverseBiasVoltageValue->MoveResize(2016,536,144,35);

   TGLabel *fTemperature = new TGLabel(fSubFrame,"Temperature (K)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fTemperature->SetTextJustify(36);
   fTemperature->SetMargins(0,0,0,0);
   fTemperature->SetWrapLength(-1);
   fSubFrame->AddFrame(fTemperature, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTemperature->MoveResize(1656,592,344,56);
   fTemperatureValue = new TGNumberEntry(fSubFrame, (Double_t) 300,9,-1,(TGNumberFormat::EStyle) 5);
   fTemperatureValue->SetName("fTemperatureValue");
   fSubFrame->AddFrame(fTemperatureValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fTemperatureValue->MoveResize(2016,600,144,35);

   TGLabel *fBField = new TGLabel(fSubFrame,"B Field (T)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fBField->SetTextJustify(36);
   fBField->SetMargins(0,0,0,0);
   fBField->SetWrapLength(-1);
   fSubFrame->AddFrame(fBField, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fBField->MoveResize(1656,656,344,56);
   fBFieldValue = new TGNumberEntry(fSubFrame, (Double_t) 0,9,-1,(TGNumberFormat::EStyle) 5);
   fBFieldValue->SetName("fBFieldValue");
   fSubFrame->AddFrame(fBFieldValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fBFieldValue->MoveResize(2016,664,144,35);

   gClient->GetColorByName("#cccccc",ucolor);
   TGLabel *fBeamCondition = new TGLabel(fSubFrame,"Beam Property",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fBeamCondition->SetTextJustify(36);
   fBeamCondition->SetMargins(0,0,0,0);
   fBeamCondition->SetWrapLength(-1);
   fSubFrame->AddFrame(fBeamCondition, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fBeamCondition->MoveResize(1656,744,512,56);

   gClient->GetColorByName("#85c2a3",ucolor);
   TGLabel *fMIP = new TGLabel(fSubFrame,"MIP (e.a/um)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fMIP->SetTextJustify(36);
   fMIP->SetMargins(0,0,0,0);
   fMIP->SetWrapLength(-1);
   fSubFrame->AddFrame(fMIP, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMIP->MoveResize(1656,824,344,56);
   fMIPValue = new TGNumberEntry(fSubFrame, (Double_t) 75,9,-1,(TGNumberFormat::EStyle) 5);
   fMIPValue->SetName("fMIPValue");
   fSubFrame->AddFrame(fMIPValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fMIPValue->MoveResize(2016,832,144,35);

   TGLabel *fIncedentAngle = new TGLabel(fSubFrame,"Angle of Incedence (degree)",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
   fIncedentAngle->SetTextJustify(36);
   fIncedentAngle->SetMargins(0,0,0,0);
   fIncedentAngle->SetWrapLength(-1);
   fSubFrame->AddFrame(fIncedentAngle, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fIncedentAngle->MoveResize(1656,888,344,56);
   fIncedentAngleValue = new TGNumberEntry(fSubFrame, (Double_t) 0,5,-1,(TGNumberFormat::EStyle) 5);
   fIncedentAngleValue->SetName("fIncedentAngleValue");
   fSubFrame->AddFrame(fIncedentAngleValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   fIncedentAngleValue->MoveResize(2016,896,144,35);

   fMainFrame->AddFrame(fSubFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   fSubFrame->MoveResize(0,0,2186,1400);

   fMainFrame->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fMainFrame->MapSubwindows();

   fMainFrame->Resize(fMainFrame->GetDefaultSize());
   fMainFrame->MapWindow();
   fMainFrame->Resize(2218,1453);
  


}

void FastSilicon(){
    new MyMainFrame();

}

void MyMainFrame::DoAnalysis()
{

    vector<Int_t> NMESH = {4,10,20};

    std::vector<TString> CNAME = {"potential3d", "potential2D", "potential1D", "Efield2D", "Efield1D", "ElecDensity", "HoleDensty", "Current"};
    std::vector<TRootEmbeddedCanvas*> RootCanvas = {fCanvas1,fCanvas2,fCanvas3,fCanvas4, fCanvas5, fCanvas6, fCanvas7, fCanvas8, fCanvas9};
    std::vector<TCanvas *> C;
    int n = 0;
    for (auto name : CNAME)
    {
        auto c = RootCanvas[n] -> GetCanvas() ;
        c->SetTitle(CNAME.at(n).Data());
        c->cd();
        setpad(gPad);
        C.push_back(c);
        if (n ==0) {
            gPad->SetTopMargin(0.05);
            gPad->SetRightMargin(0.18);
        }
        if (n == 1 || n == 3)
            gPad->SetRightMargin(0.2);
        if (n == 5 || n == 6)
        {
            gPad->SetRightMargin(0.18);
            gPad->SetBottomMargin(0.18);
            gPad->SetTheta(6.882412);
            gPad->SetPhi(-34.07409);
        }
        n++;
    }

    cout << "back : " << NMESH.back() << endl;
    for (auto nmesh : NMESH)
        SolvePoisson(nmesh, 0, C);

    InjectParticle(theta);

    Int_t nbins = 200;
    Double_t logbins[nbins + 1];
    Double_t low = 1e-12;
    //Double_t high = 100e-9;
    Double_t high = 4e-9;
    Double_t logbw = (log(high) - log(low)) / nbins;
    std::vector<Double_t> time;
    for (int ij = 0; ij <= nbins; ij++)
        logbins[ij] = low * exp(ij * logbw);

    for (auto i = 0; i < nbins; i++)
    {
        time.push_back(logbins[i]);
        SolvePoisson(NMESH.back(), logbins[i], C);
        SolveCurrent(logbins[i], NMESH.back(), C);
    }


    auto fr = new TFile("NewSimulationResult.root", "recreate");
    gcurrent->SetName("Simulation");
    gncurrent->SetName("ecurrent");
    ghcurrent->SetName("hcurrent");
    gcurrent->Write();
    gncurrent->Write();
    ghcurrent->Write();
    fr->Close();

    std::cout << "finished" << std::endl;
}

void MyMainFrame::StartAnalysis(){

   xdim = fDetectorWidthValue->GetNumber();
   ydim = fDetectorLengthValue->GetNumber();
   zdim = fDetectorDepthValue->GetNumber();
   spacedopingdensity = fNtypeDopingDensityValue->GetNumber();
   spacedopingdensityp = fPtypeDoingDensityValue->GetNumber();

   reversebiasvoltage = fReverseBiasVoltageValue->GetNumber();
   Temp = fTemperatureValue->GetNumber();

   theta = fIncedentAngleValue->GetNumber();
   MIP = fMIPValue->GetNumber();

   DoAnalysis();

}

void MyMainFrame::SolveCurrent(Double_t time,Double_t nmesh , std::vector<TCanvas *> C){

    ndensity->Reset();
    hdensity->Reset();
    Double_t dt = time - ptime;
    Double_t ncurrent = 0;
    Double_t hcurrent = 0;
    if (ptime == 0){
        
        gcurrent = new TGraph;
        gncurrent = new TGraph;
        ghcurrent = new TGraph;
        gcurrent->SetPoint(gcurrent->GetN(), 0, 0);
        gncurrent->SetPoint(gncurrent->GetN(), 0, 0);
        ghcurrent->SetPoint(ghcurrent->GetN(), 0, 0);
        setstyle(gcurrent);
        gcurrent->GetXaxis()->SetTitle("time (s)");
        gcurrent->GetYaxis()->SetTitle("Current (A)");
    }

    for (auto i = 0; i < nx.size(); i++)
    {

        Int_t binx = EFX->GetXaxis()->FindBin(nx[i]);
        Int_t biny = EFY->GetYaxis()->FindBin(ny[i]);
        Int_t binz = EFZ->GetZaxis()->FindBin(nz[i]);

        if (i == 100)
            cout << "nx ny nz " << nx[i] << " " << ny[i] << " " << nz[i] << endl;
        Double_t efx = EFX->GetBinContent(binx, biny, binz);
        Double_t efy = EFY->GetBinContent(binx, biny, binz);
        Double_t efz = EFZ->GetBinContent(binx, biny, binz);
        //Double_t MU_E = MU_E0 * TMath::Power(1./(1+TMath::Power(MU_E0*TMath::Abs(ey)/VESAT,2)),1./2);
        Double_t vx = -1 * MU_E0 * efx / TMath::Sqrt(1 + TMath::Power(MU_E0 * TMath::Abs(efx) / VESAT, 2));
        Double_t vy = -1 * MU_E0 * efy / TMath::Sqrt(1 + TMath::Power(MU_E0 * TMath::Abs(efy) / VESAT, 2));
        Double_t vz = -1 * MU_E0 * efz / TMath::Sqrt(1 + TMath::Power(MU_E0 * TMath::Abs(efz) / VESAT, 2));
        Double_t dx = vx * dt * 1e6;
        Double_t dy = vy * dt * 1e6;
        Double_t dz = vz * dt * 1e6;
        if (i == 100)
            cout << "dx dy dz dt efz " << dx << " " << dy << " " << dz << " " << dt << " " << efz << endl;
        nx[i] += dx;
        ny[i] += dy;
        nz[i] += dz;
        ndensity->Fill(nx[i], ny[i], nz[i]);
        //if (nz[i]<zdim)
            //ncurrent += ec * (vz * efz);
    }

    for (auto i = 0; i < hx.size(); i++){

        Int_t binx = EFX->GetXaxis()->FindBin(hx[i]);
        Int_t biny = EFY->GetYaxis()->FindBin(hy[i]);
        Int_t binz = EFZ->GetZaxis()->FindBin(hz[i]);

        Double_t efx = EFX->GetBinContent(binx, biny, binz);
        Double_t efy = EFY->GetBinContent(binx, biny, binz);
        Double_t efz = EFZ->GetBinContent(binx, biny, binz);

        Double_t vx = MU_H0 * efx / (1 + TMath::Power(MU_H0 * TMath::Abs(efx) / VHSAT, 1));
        Double_t vy = MU_H0 * efy / (1 + TMath::Power(MU_H0 * TMath::Abs(efy) / VHSAT, 1));
        Double_t vz = MU_H0 * efz / (1 + TMath::Power(MU_H0 * TMath::Abs(efz) / VHSAT, 1));
        Double_t dx = vx * dt * 1e6;
        Double_t dy = vy * dt * 1e6;
        Double_t dz = vz * dt * 1e6;
        if (i == 100)
            cout << "hole dx dy dz dt efz " << dx << " " << dy << " " << dz << " " << dt << " " << efz << endl;
        hx[i] += dx;
        hy[i] += dy;
        hz[i] += dz;
        hdensity->Fill(hx[i], hy[i], hz[i]);
        //if (hz[i]>0)
            // hcurrent += ec * (vz * efz);
    }


    for (auto i = 1; i <= EFX -> GetNbinsX(); i++)
    {
        for (auto j = 1; j <= EFX -> GetNbinsY(); j++)
        {
            for (auto k = 1; k <= EFX -> GetNbinsZ(); k++)
            {
                Double_t efx = EFX->GetBinContent(i, j, k);
                Double_t efy = EFY->GetBinContent(i, j, k);
                Double_t efz = EFZ->GetBinContent(i, j, k);
                Double_t binwidthx = EFX->GetXaxis()->GetBinWidth(i);
                Double_t binwidthy = EFX->GetYaxis()->GetBinWidth(j);
                Double_t binwidthz = EFX->GetZaxis()->GetBinWidth(k);
                Double_t nd = ndensity -> GetBinContent(i, j, k);
                Double_t ndm = ndensity -> GetBinContent(i, j, k-1);
                Double_t hd = hdensity -> GetBinContent(i, j, k);
                Double_t hdm = hdensity -> GetBinContent(i, j, k-1);
                Double_t MU_EZ = MU_E0 * TMath::Sqrt(1. / (1 + (MU_E0 * efz / VESAT)*(MU_E0 * efz / VESAT)));
                Double_t MU_HZ = MU_H0 * 1. / (1 + MU_H0 * TMath::Abs(efz) / VHSAT);
                Double_t volume = (binwidthx * binwidthy * binwidthz) * 1e-18; //micrometer -> meter
                Double_t area = (binwidthx * binwidthy) * 1e-12; // micrometer -> meter

                ncurrent += 1 * ec / volume * nd * efz * area * MU_EZ;
                //Double_t ediffcurr = kB * Temp * MU_EZ * (nd - ndm) / volume / (binwidthz * 1e-6)*area;
                //if (TMath::Abs(ediffcurr)>0)
                //    cout<<"ncurrent : diffcurent " << ncurrent << " " << ediffcurr << endl;
                //ncurrent -= ediffcurr;


                hcurrent += 1 * ec / volume * hd * efz * area * MU_HZ;
                //Double_t hdiffcurr = kB * Temp * MU_HZ * (hd - hdm) / volume / (binwidthz * 1e-6)*area;
                //hcurrent -= hdiffcurr;
                
            }
        }
    }


    Double_t totalcurrent = fabs(ncurrent) + fabs(hcurrent);
    gcurrent->SetPoint(gcurrent->GetN(), time, totalcurrent / nmesh);
    //gcurrent->SetPoint(gcurrent->GetN(), time, totalcurrent / 140);
    gncurrent->SetPoint(gncurrent->GetN(), time,  fabs(ncurrent) / nmesh);
    ghcurrent -> SetPoint(ghcurrent->GetN(), time,  fabs(hcurrent) / nmesh);
    Int_t timepower = TMath::Log10(time) - 1;
    Double_t timea = time * TMath::Power(10, -1 * timepower);
    ndensity->SetTitle(Form ("Free electron density at #it{t} = %2.1f #times 10^{%d} s",timea, timepower ));
    hdensity->SetTitle(Form("Free hole density at #it{t} = %2.1f #times 10^{%d} s",timea, timepower ));
    C.at(5)->cd();
    ndensity->Draw("box2 ");
    C.at(6)->cd();
    hdensity->Draw("box2 ");
    C.at(7)->cd();
    gPad->SetLogx(1);

    gcurrent->GetXaxis();
    gncurrent->SetLineColor(4);
    ghcurrent->SetLineColor(2);
    gcurrent->Draw();
    gncurrent->Draw("same");
    ghcurrent->Draw("same");
    C.at(5)->Update();
    C.at(6)->Update();
    C.at(7)->Update();
    ptime = time;


    auto legend = new TLegend(0.607553, 0.68265, 0.90, 0.859187);
    legend->AddEntry(gcurrent, "Total current", "l");
    legend->AddEntry(gncurrent, "Electron current", "l");
    legend->AddEntry(ghcurrent, "Hole current", "l");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.0473934);
    legend->Draw();


}

void MyMainFrame::InjectParticle(Double_t theta){
    theta = theta * TMath::Pi() / 180.; //transform from degree to radian
    const Double_t R = zdim / TMath::Cos(theta);
    cout << "R = " << R <<" um" << endl;
    const Int_t totalmip = R * MIP;
    cout << "total MIP in the device = " << totalmip<<" ea" << endl;
    TRandom *rand = new TRandom(0);

    for (auto i = 0; i < totalmip; i++){
        Double_t r = rand->Uniform(0, R);
        Double_t x = 0;
        Double_t y = r * TMath::Sin(theta);
        Double_t z = r * TMath::Cos(theta);
        nx.push_back(x);
        ny.push_back(y);
        nz.push_back(z);

        r = rand->Uniform(0, R);
        x = 0;
        y = r * TMath::Sin(theta);
        z = r * TMath::Cos(theta);

        hx.push_back(x);
        hy.push_back(y);
        hz.push_back(z);
    }
}


void setstyle(TH3F *h)
{

    h->GetXaxis()->SetTitle("pixel width (#it{x}) [#mum]");
    h->GetYaxis()->SetTitle("pixel length (#it{y}) [#mum]");
    h->GetZaxis()->SetTitle("pixel depth (#it{z}) [#mum]");
    h->GetXaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetZaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->GetZaxis()->CenterTitle(1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetZaxis()->SetTitleSize(0.06);
    h->SetStats(0);
    h->SetTitle(0);
}
void setstyle(TH2F *h)
{

    h->GetXaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->GetZaxis()->CenterTitle(1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetZaxis()->SetTitleSize(0.06);
    h->SetStats(0);
}
void setstyle(TH1F *h)
{

    h->GetXaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->SetStats(0);
}
void setstyle(TGraph *h)
{

    h->GetXaxis()->SetTitleOffset(0.8);
    h->GetYaxis()->SetTitleOffset(0.8);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
}

void setpad(TVirtualPad *pad)
{
    pad->SetTopMargin(0.1);
    pad->SetBottomMargin(0.13);
    pad->SetRightMargin(0.05);
    pad->SetLeftMargin(0.13);
    pad->SetTheta(63.33279);
    pad->SetPhi(-24.80808);
}

void MyMainFrame::SolvePoisson(const int nmesh,const Double_t time, std::vector <TCanvas*> C)
{

    //base = new TH3F("base","base",mxy.at(m),-1*xdim/2.,xdim/2., mxy.at(m),-1*xdim/2.,xdim/2., mz.at(m), 0, zdim );
    //base = new TH3F("base", "base", nmesh, -1 * xdim / 2., xdim / 2., nmesh, -1 * xdim / 2., xdim / 2., nmesh, 0, zdim);
    base = new TH3F("base", "base", nmesh, -1 * xdim / 2., xdim / 2., nmesh, -1 * xdim / 2., xdim / 2., nmesh, 0, zdim);
    setstyle(base);

    Int_t sizex = nmesh + 2; // under and overflown taken into account
    Int_t sizez = nmesh + 2; // under and overflown taken into account
    float potential[sizex][sizex][sizez];
    Double_t deltah = 300e-6 / nmesh;
    // boundary condition

    for (auto i = 0; i < sizex; i++)
    {
        for (auto j = 0; j < sizex; j++)
        {
            for (auto k = 0; k < sizez; k++)
            {
                potential[i][j][k] = 0;
            }
        }
    }

    for (auto i = 0; i < sizex; i++)
    {
        for (auto j = 0; j < sizex; j++)
        {
            for (auto k = 0; k < sizez; k++)
            {
                if (globalmeshcounter > 0)
                {
                    Double_t xcenter = base->GetXaxis()->GetBinCenter(i);
                    Double_t ycenter = base->GetYaxis()->GetBinCenter(j);
                    Double_t zcenter = base->GetZaxis()->GetBinCenter(k);
                    Int_t binx = pbase->GetXaxis()->FindBin(xcenter);
                    Int_t biny = pbase->GetYaxis()->FindBin(ycenter);
                    Int_t binz = pbase->GetZaxis()->FindBin(zcenter);
                    potential[i][j][k] = pbase->GetBinContent(binx, biny, binz);
                    //potential[i][j][k] = pbase->Interpolate(xcenter, ycenter, zcenter);
                }
            }
        }
    }
    if (globalmeshcounter > 0)
        pbase->Delete();

    for (auto i = 0; i < sizex; i++)
    {
        for (auto j = 0; j < sizex; j++)
        {
            potential[i][j][0] = 0;
            potential[i][j][sizez - 1] = reversebiasvoltage -1. ;  // drift potential taken into accout as 1V. 
        }
    }

    Double_t value = 0;
    Double_t sumpotential = 1e15, presumpotential = 1e10;
    std::cout << "meshsize = " << nmesh << std::endl;
    //while (fabs((presumpotential - sumpotential) / sumpotential*mz.at(m)*mxy.at(m)*mxy.at(m)) > 1e-4/mz.at(m))

    if (nx.size() > 0)
    {
        ndensity = (TH3F *)base->Clone();
        hdensity = (TH3F *)base->Clone();
        ndensity->Reset();
        hdensity->Reset();
        for (auto i = 0; i < nx.size(); i++)
        {
            ndensity->Fill(nx.at(i), ny.at(i), nz.at(i));
            hdensity->Fill(hx.at(i), hy.at(i), hz.at(i));
        }
    }

    while (fabs((presumpotential - sumpotential) / sumpotential) > accurate / nmesh)
    {
        std::cout <<"t = "<<time<< " delta : " << abs((presumpotential - sumpotential) / sumpotential) << " accurate set : "<<accurate/nmesh << endl;
        presumpotential = sumpotential;
        sumpotential = 0;


        for (auto i = 1; i < sizex - 1; i++)
        {
            for (auto j = 1; j < sizex - 1; j++)
            {
                for (auto k = 1; k < sizez - 1; k++)
                {
                    Int_t newi = i;
                    Int_t newj = j;
                    if (i == 1 || i == (sizex - 2)) // size boundary condition applied
                        newi = sizex / 2;
                    if (j == 1 || j == (sizex - 2))
                        newj = sizex / 2;
                    value = potential[newi - 1][newj][k] + potential[newi + 1][newj][k] + potential[newi][newj - 1][k] + potential[newi][newj + 1][k] + potential[newi][newj][k - 1] + potential[newi][newj][k + 1];
                    //if (i == sizex / 2 && j == sizex / 2 && k == sizez / 2)
                    //    cout << "potential = " << deltah * deltah * ec * spacedopingdensity / pi4per / 6 << endl;
                    //    potential[i][j][k] = value / 6 + deltah * deltah * ec * (spacedopingdensity-spacedopingdensityp*2/base->GetZaxis()->GetBinWidth(k)) / pi4per / 6;
                    potential[i][j][k] = value / 6 ;
                    Double_t npspacedopingdensity = spacedopingdensity;
                    //if (k ==1)
                    //    npspacedopingdensity = -1*spacedopingdensityp * 2 / base->GetZaxis()->GetBinWidth(k);
                    Double_t spacedopingcontribution = deltah * deltah * ec * npspacedopingdensity / pi4per / 6;
                    Double_t depletionlength = TMath::Sqrt(2*siliconrelativepermittivity*vacuumpermitttivity*reversebiasvoltage/ec/spacedopingdensity)*1e6;
                    if (base->GetZaxis()->GetBinCenter(k)<depletionlength) potential[i][j][k] += spacedopingcontribution;
                    if (time > 0)
                    {
                        if (hdensity->GetBinContent(i,j,k)>0 || ndensity->GetBinContent(i,j,k)>0){
                            Double_t carrierdensity = (hdensity->GetBinContent(i, j, k) - ndensity->GetBinContent(i, j, k)) * ndensity->GetXaxis()->GetBinWidth(i) / ndensity->GetYaxis()->GetBinWidth(j) / ndensity->GetZaxis()->GetBinWidth(k) * 1e18;
                            Double_t carriervoltage = deltah * deltah * ec * carrierdensity / pi4per / 6 ;
                            //cout << "carrier charge voltage = " << carriervoltage << endl;
                            //cout << "carrier charge density = " << carrierdensity << endl;
                            //potential[i][j][k] += carriervoltage;
                        }
                    }
                    sumpotential += potential[i][j][k];
                    base->SetBinContent(i, j, k, potential[i][j][k]);
                }
            }
        }

        // potential 3d
        C.at(0)->cd();
        base->SetTitle("Potential 3-D");
        base->Draw("box2 ");

        // potential 2d
        C.at(1)->cd();
        auto temp = (TH3F *)base->Clone();
        temp->GetXaxis()->SetRange(temp->GetNbinsX() / 2, temp->GetNbinsX() / 2);
        auto baseprojection = (TH2F *)temp->Project3D("zy");
        setstyle(baseprojection);
        //baseprojection->SetTitle("Potential at #it{x} = 0");
        baseprojection->SetTitle(0);
        baseprojection->GetXaxis()->SetTitle("pixel length (#it{y}) [#mum]");
        baseprojection->GetYaxis()->SetTitle("pixel depth (#it{z}) [#mum]");
        baseprojection->GetZaxis()->SetTitle("potential [V]");
        baseprojection->SetTitle("Potential 2-D");
        baseprojection->Draw("colz");

        // potential 1d
        C.at(2)->cd();
        TH1F *base1d = (TH1F *)baseprojection->ProjectionY("", baseprojection->GetNbinsX() / 2, baseprojection->GetNbinsX() / 2, "");
        base1d->SetStats(0);
        //base1d->SetTitle("Potential at #it{x}, #it{y} = 0");
        base1d->SetTitle(0);
        setstyle(base1d);
        base1d->GetYaxis()->SetTitle("Potential [V]");
        base1d->GetXaxis()->SetTitle("pixel depth (#it{z}) [#mum]");
        base1d->SetLineColor(1);
        base1d->SetLineWidth(2);
        base1d->SetTitle("Potential 1-D");
        base1d->Draw("c");

        // Draw electric field in z direction
        C.at(3)->cd();
        auto efz = (TH2F *)baseprojection->Clone();
        setstyle(efz);
        efz->GetZaxis()->SetTitle("Electric field [V/m]");
        for (auto j = 1; j < sizex - 1; j++)
        {
            for (auto k = 1; k < sizez - 1; k++)
            {
                Double_t binwidth = efz->GetYaxis()->GetBinWidth(k) * 1e-6;
                efz->SetBinContent(j, k, (potential[sizex / 2][j][k - 1] - potential[sizex / 2][j][k]) / binwidth);
            }
        }
        efz->SetTitle("Electric Field in #it{z}-direction (#it{E}_{#it{z}})");
        efz->Draw("colz ARROW");

        EFX = (TH3F*)base->Clone();
        EFX->Reset();
        EFY = (TH3F*)base->Clone();
        EFY->Reset();
        EFZ = (TH3F*)base->Clone();
        EFZ->Reset();

        for (auto i = 1; i < sizex - 1; i++)
        {
            for (auto j = 1; j < sizex - 1; j++)
            {
                for (auto k = 1; k < sizez - 1; k++)
                {
                    Double_t binwidthx = EFX->GetXaxis()->GetBinWidth(i) * 1e-6;
                    Double_t binwidthy = EFY->GetYaxis()->GetBinWidth(j) * 1e-6;
                    Double_t binwidthz = EFZ->GetZaxis()->GetBinWidth(k) * 1e-6;
                    EFX->SetBinContent(i,j, k, (potential[i-1][j][k] - potential[i][j][k]) / binwidthx);
                    EFY->SetBinContent(i,j, k, (potential[i][j-1][k] - potential[i][j][k]) / binwidthy);
                    EFZ->SetBinContent(i,j, k, (potential[i][j][k-1] - potential[i][j][k]) / binwidthz);
                }
            }
        }

        C.at(4)->cd();
        TH1F *base1de = (TH1F *)efz->ProjectionY("elecz", efz->GetNbinsX() / 2, efz->GetNbinsX() / 2, "");
        base1de->SetStats(0);
        //base1d->SetTitle("Potential at #it{x}, #it{y} = 0");
        base1de->SetTitle(0);
        setstyle(base1de);
        base1de->GetYaxis()->SetTitle("Electric field [V/m]");
        base1de->GetXaxis()->SetTitle("pixel depth (#it{z}) [#mum]");
        base1de->SetLineColor(1);
        base1de->SetLineWidth(2);
        base1de->SetTitle("Electric field at #it{x,y} = 0");
        base1de->Draw("c");
        C.at(5)->cd();
        if (ndensity) ndensity->Draw("box2 ");
        C.at(6)->cd();
        if (hdensity) hdensity->Draw("box2 ");

        // update canvases
        for (auto c : C)
            c->Update();
        temp->Delete();
        //base1de->Delete();
        //base1d->Delete();
        //baseprojection->Delete();
        //base1d->Delete();
    }
    pbase = (TH3F *)base->Clone();
    globalmeshcounter++;
    std::cout << "check t = " << time << " delta : " << abs((presumpotential - sumpotential) / sumpotential) << " accurate set : " << accurate / nmesh << endl;
}
     