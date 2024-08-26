#include "include/Filipad.h"
#include "include/rootcommon.h"

double integralv2(TF1 *fflow, TF1 *fspectra, double lpt, double hpt); // via integral
double SampingMethod(TF1 *fflow, TF1 *fspectra, double lpt, double hpt); // via simple samping
double FitVN(double *x,double *par); // Vn(pt) fit function
double LevyTsallisF0(double *x, double *par){ // Spectra Tsallis fit function, Eq1 : https://arxiv.org/pdf/1510.05449.pdf https://arxiv.org/abs/1911.04878
double mass = par[3];
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / (par[1]*par[2]) ), -par[1] );
}
// Centrality and information for v2 tables from HEPData
const int NC = 6;
double cent[NC+1] =           { 0, 5,10,20,30,40,50}; //https://doi.org/10.17182/hepdata.8621
int v2Hepdatatable[NC] = {15,23,33,43,53,63}; // https://doi.org/10.17182/hepdata.83737
TGraphAsymmErrors *gr_v2pt[NC];
TGraphErrors *gr_v2ptRun3[NC];
TGraphAsymmErrors *gr_dndpt[NC];
TH1F *h_dndpt[NC]; // gr -> hist
TGraphAsymmErrors *gr_v2inte0230pub; 


// Actual calculations 
double GetIntegratedv2(int iC,double lpt, double hpt); // for one centrality bin for a given pt range.

// calculation and drawing the ratios
void GetScaleForPtDiff();// Run over various pt bins and draw centrality dependence integrated v2
void Runs();// Run one pt bins for all centrality bins and compare to the publised results.
void Run3(int ic); // 
void Run3allcent(); // Run over all centrality bins for Run3

// Test purpose..
// Test fit with miuit2 for one centrality
void FitPtMinuit(int iC);


// Loading HEPData into ROOT objects.
void LoadHEPData() {
	TFile *fflow = TFile::Open("published/HEPData-ins1666817-v1-5TeVPbPb_flow.root");  //https://doi.org/10.17182/hepdata.83737
	TFile *fspectra = TFile::Open("published/HEPData-ins1657384-v1-5TeVPbPb_spectra.root"); //https://doi.org/10.17182/hepdata.86210
    TFile *fRun3 = TFile::Open("data/fout_FT0C_v2_cos_v2.root");  // Junlee
    for(int ic=0;ic<NC;ic++) {
	   gr_dndpt[ic] = (TGraphAsymmErrors*)fspectra->Get(Form("Table 2/Graph1D_y%d",ic+1));
       TH1F *hy1 = (TH1F*)fspectra->Get(Form("Table 2/Hist1D_y%d",ic+1));
       TH1F *hy1_e1 = (TH1F*)fspectra->Get(Form("Table 2/Hist1D_y%d_e1",ic+1));
       TH1F *hy1_e2 = (TH1F*)fspectra->Get(Form("Table 2/Hist1D_y%d_e2",ic+1));
       h_dndpt[ic] = (TH1F*)hy1->Clone();
       makeHistHEPDATA(hy1,hy1_e1,hy1_e2,h_dndpt[ic]);
    }
    for(int ic=0;ic<NC;ic++) {
		gr_v2pt[ic] = (TGraphAsymmErrors*)fflow->Get(Form("Table %d/Graph1D_y1",v2Hepdatatable[ic]));
          //run 3
        gr_v2ptRun3[ic] = (TGraphErrors*)fRun3->Get(Form("gPt_%d",ic));
	}
	gr_v2inte0230pub = (TGraphAsymmErrors*)fflow->Get("Table 1/Graph1D_y1");
}

double SampingMethod(TF1 *fflow, TF1 *fspectra, double lpt, double hpt){
    for(int j=0;j<3;j++) cout << fspectra->GetParameter(j) << endl;
    TH1D *hv2 = new TH1D("hv2","",500,0.,1.0);
    Int_t Nevt  = 1e3;
    for(Int_t i = 0; i < Nevt; i++){
        double pt = fspectra->GetRandom(lpt,hpt);
        double v2 = fflow->Eval(pt);
        hv2->Fill(v2);
    }
    cout << "Simpling Done" << endl;
    cout << lpt <<"-"<<hpt <<"="<< hv2->GetMean() << endl;
    return hv2->GetMean();
}

double integralv2(TF1 *fflow, TF1 *fspectra, double lpt, double hpt){

    int nstep=1000;
    double bw = (hpt-lpt)/nstep;
    double norm = fspectra->Integral(lpt,hpt);
    double sum = 0;

    for(int i=0;i<nstep;i++) {
        double lx = lpt+bw*i;
        double hx = lpt+bw*(i+1);
        double v2 = fflow->Eval(lx+bw/2.); // need bin by bin?
        double sub = fspectra->Integral(lx,hx);
        double weight = sub/norm;
        sum = sum + v2*weight;
        //cout << i <<"\t"<< lx <<"\t"<< hx <<"\t"<< v2 <<"\t"<< norm <<"\t"<< sub <<"\t"<< weight <<"\t" <<sum<< endl;
    }       
    return sum;
}


// Main function to calculate integrated v2
double GetIntegratedv2(int iC, double lpt, double hpt){
    TString label = Form("%2.0f-%2.0f%% %.1f<p_{T}<%.1f",cent[iC],cent[iC+1],lpt,hpt);
    const char* fitter = "Minuit2";
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fitter);
    TString fitterType(fitter);
    gStyle->SetOptFit();
    
    gStyle->SetStatX(0.6);gStyle->SetStatY(0.5);
    TF1 *fLevy = new TF1("fLevy",LevyTsallisF0,0,5.,4); // function to describe the pt spectrum
    fLevy->SetParNames("N","n","C","mass");
    fLevy->SetParLimits( 0, 1e1, 5e03);
    fLevy->SetParLimits( 1, 4.5, 15);
    fLevy->SetParLimits( 2, 0.05, 0.22);
    fLevy->SetParLimits( 3, 0.05, 0.25);
    //fLevy->SetParameters(1.2e+03,7,0.02);
    fLevy->Update();

    gr_dndpt[iC]->Fit(fLevy,"R");
    // Fit vn
    gStyle->SetStatX(0.6);gStyle->SetStatY(0.75);
    TF1 *f_vn = new TF1("f_vn",FitVN,0.,20.,6);
    f_vn->SetParameters(0,1.2,2.2,2,0.1,1.2);
    f_vn->SetParName(0,"a");
    f_vn->SetParName(1,"n");
    f_vn->SetParName(2,"lamda");
    f_vn->SetParName(3,"m");
    f_vn->SetParName(4,"c1");
    f_vn->SetParName(5,"c2");
    f_vn->SetLineColor(2);
    gr_v2pt[iC]->Fit(f_vn,"R");

    double lowx=0., highx=20.;
    double ly=0., hy=0.20;
    int logx = 0;

    // v2
    Filipad *fpad;
    fpad = new Filipad(1, 0.9, 0.4, 100, 100, 1.1,2);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(logx); p->SetLogy(0); p->cd();
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "p_{T}(GeV/c)", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.65,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.06);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    gr_v2pt[iC]->SetMarkerStyle(20);
    gr_v2pt[iC]->Draw("psame");
    leg->AddEntry(gr_v2pt[iC],label,"lp");
    leg->Draw();
  

    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1.1,1.1);
    hset( *hfr1, "p_{T}(GeV/c)","#frac{Fit-Data}{Fit}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    
    TGraphAsymmErrors *gr_ratio1;
    gr_ratio1 = GetDataOverTheory(gr_v2pt[iC],f_vn);
    gr_ratio1->Draw("same");
    gPad->GetCanvas()->SaveAs(Form("figs/v2fit_c%d.pdf",iC));
    // pt spectra
    logx=1;
    Filipad *fpad1;
    fpad1 = new Filipad(2, 0.9, 0.4, 100, 100, 1.1,2);
    fpad1->Draw();
    //==== Upper pad
    TPad *p1 = fpad1->GetPad(1); //upper pad
    p1->SetTickx(); p1->SetLogx(logx); p1->SetLogy(1); p1->cd();
    ly=1.5e-5,hy=9e3;
    lowx=0.08,highx=70.;
    TH2F *hfr01 = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr01, "p_{T}(GeV/c)", "#frac{1}{N_{evt}} #frac{d^{2} N}{dp_{T}d#eta}(Gev^{-1}c)",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr01->Draw();
    //Legend definition
    TLegend *leg1 = new TLegend(0.65,0.5,0.85,0.78,"","brNDC");
    leg1->SetTextSize(0.06);leg1->SetBorderSize(0);leg1->SetFillStyle(0);//legend settings;
    gr_dndpt[iC]->SetMarkerStyle(20);
    gr_dndpt[iC]->Draw("psame");
    leg1->AddEntry(gr_dndpt[iC],label,"lp");
    leg1->Draw();

    //==== Lower pad
    p1 = fpad1->GetPad(2);
    p1->SetTickx(); p1->SetGridy(1); p1->SetLogx(logx), p1->SetLogy(0); p1->cd();
    TH2F *hfr02 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1.1,1.1);
    hset( *hfr02, "p_{T}(GeV/c)","#frac{Fit-Data}{Fit}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr02->Draw();
    
    TGraphAsymmErrors *gr_ratio2;
    gr_ratio2 = GetDataOverTheory(gr_dndpt[iC],fLevy);
    gr_ratio2->Draw("same");

    //return integralv2(f_vn,fLevy, double lpt, double hpt);
    gPad->GetCanvas()->SaveAs(Form("figs/dndptfit_c%d.pdf",iC));
    
    return SampingMethod(f_vn,fLevy,lpt,hpt);
}

// Run over various pt bins
void GetScaleForPtDiff(){
	LoadHEPData();
	const int Npt=5;
	double startPtbins[Npt] = {0.,0.2,0.28,0.5,0.8};
	double endPt = 5.0;
	double intv2[Npt][NC];
	double centmean[NC];
	TGraphAsymmErrors *gr_v2intForPtbins[Npt];
	for(int i=0;i<Npt;i++){
		for(int ic=0;ic<NC;ic++) {
			intv2[i][ic] = GetIntegratedv2(ic,startPtbins[i],endPt);
			centmean[ic] = cent[ic]+(cent[ic+1]-cent[ic])/2.;
		}
	    gr_v2intForPtbins[i] = new TGraphAsymmErrors(NC,centmean,intv2[i],0,0,0,0);
	}
	Filipad *fpad;
	fpad = new Filipad(1, 0.9, 0.4, 100, 100, 1.1,2);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    double lowx = 0,highx=60.;
    double ly=0,hy=0.2;
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "centrality[%]", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.3,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.06);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    for(int i=0;i<Npt;i++){
    	gr_v2intForPtbins[i]->SetMarkerStyle(20+i);
    	gr_v2intForPtbins[i]->SetMarkerColor(1+i);
    	gr_v2intForPtbins[i]->Draw("psame");
    	leg->AddEntry(gr_v2intForPtbins[i],Form("%.1f<p_{T}<%.1f",startPtbins[i],endPt),"lp");
    }
    leg->AddEntry(gr_v2inte0230pub,"published 0.2<p_{T}<3.0","lp");
    gr_v2inte0230pub->SetMarkerStyle(30);
    gr_v2inte0230pub->Draw("psame");

    leg->Draw();
  
    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1.0,0.5);
    hset( *hfr1, "centrality[%]","#frac{Ref-Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphAsymmErrors *gr_ratio[Npt];
    int iref=1;
    for(int i=0;i<Npt;i++){
      gr_ratio[i] = GetDataOverTheory(gr_v2intForPtbins[i],gr_v2intForPtbins[iref]);
      gr_ratio[i]->SetMarkerStyle(20+i);
      gr_ratio[i]->SetMarkerColor(1+i);
      gr_ratio[i]->Draw("psame");
	}

    TFile *fout = new TFile("vnPtintegrated.root","recreate");
    fout->cd();
    for(int i=0;i<Npt;i++){
        gr_v2intForPtbins[i]->SetName(Form("gr_v2intForPtbins%02d",i));
        gr_v2intForPtbins[i]->SetTitle(Form("%.1f<p_{T}<%.1f",startPtbins[i],endPt));
        gr_v2intForPtbins[i]->Write();
    }
    fout->Close();
    gPad->GetCanvas()->SaveAs("figs/integratedv2_ptdiff.pdf");
}

double FitVN(double *x,double *par){
    double a = par[0];
    double n = par[1];
    double lamda = par[2];
    double m = par[3];
    double c1 = par[4];
    double c2 = par[5];
    if(x[0]<4.3) {
        return  (a + 1/TMath::Power(x[0],n))*(TMath::Power(x[0]/lamda,m)/(1+TMath::Power(x[0]/lamda,m)));            }
    else {
        return  (c1*TMath::Power(x[0],c2));
    }
}

// End of the main routines the flowing things are for test purpose and plotting.

// Test fit with miuit2
void FitPtMinuit(int iC){
    LoadHEPData();
    TString label = Form("%2.0f-%2.0f%%",cent[iC],cent[iC+1]);
    const char* fitter = "Minuit2";
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fitter);
    TString fitterType(fitter);


    gStyle->SetOptFit();

    gStyle->SetStatX(0.6);gStyle->SetStatY(0.5);
    TF1 *fLevy = new TF1("fLevy",LevyTsallisF0,0,20.,4); // function to describe the pt spectrum
    fLevy->SetParNames("N","n","C","mass");
    /*
    fLevy->SetParLimits( 0, 1e1, 5e6);
    fLevy->SetParLimits( 1, 1, 25);
    fLevy->SetParLimits( 2, 0.05, 5.);
    fLevy->SetParLimits( 3, 0.05, 0.25);
    fLevy->SetParameters(1.2e+03,7,0.02);
    */
    fLevy->Update();

    //h_dndpt[iC]->Fit(fLevy,"R");

    bool ok = true;
    Int_t npass = 20;
    for (Int_t pass=0;pass<npass;pass++) {
      if (pass%100 == 0) printf("pass : %d\n",pass);
      else printf(".");
      if (pass == 0) fLevy->SetParameters(2e03,9,0.02,0.05);
      for (Int_t i=0;i<1e6;i++) {
         h_dndpt[iC]->Fill(fLevy->GetRandom());
      }
      int iret = h_dndpt[iC]->Fit(fLevy,"Q0");
      ok &= (iret == 0);
      if (iret!=0) Error("DoFit","Fit pass %d failed !",pass);
    }

      // pt spectra
    double lowx=0., highx=20.;
    double ly=0., hy=0.20;
    int logx = 0;
    logx=1;
    Filipad *fpad1;
    fpad1 = new Filipad(iC+1, 0.9, 0.4, 100, 100, 1.1,6);
    fpad1->Draw();
    //==== Upper pad
    TPad *p1 = fpad1->GetPad(1); //upper pad
    p1->SetTickx(); p1->SetLogx(logx); p1->SetLogy(1); p1->cd();
    ly=1.5e-5,hy=9e3;
    lowx=0.08,highx=70.;
    TH2F *hfr01 = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr01, "p_{T}(GeV/c)", "#frac{1}{N_{evt}} #frac{d^{2} N}{dp_{T}d#eta}(Gev^{-1}c)",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr01->Draw();
    //Legend definition
    TLegend *leg1 = new TLegend(0.65,0.5,0.85,0.78,"","brNDC");
    leg1->SetTextSize(0.06);leg1->SetBorderSize(0);leg1->SetFillStyle(0);//legend settings;
    gr_dndpt[iC]->SetMarkerStyle(20);
    gr_dndpt[iC]->Draw("psame");
    fLevy->Draw("same");
    leg1->AddEntry(gr_dndpt[iC],label,"lp");
    leg1->Draw();

    //==== Lower pad
    p1 = fpad1->GetPad(2);
    p1->SetTickx(); p1->SetGridy(1); p1->SetLogx(logx), p1->SetLogy(0); p1->cd();
    TH2F *hfr02 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1.1,1.1);
    hset( *hfr02, "p_{T}(GeV/c)","#frac{Fit-Data}{Fit}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr02->Draw();
    
    TGraphAsymmErrors *gr_ratio2;
    gr_ratio2 = GetDataOverTheory(gr_dndpt[iC],fLevy);
    gr_ratio2->Draw("same");
}
// this was for first run to check with one published pt bin.
// Run one pt bins for all centrality bins and compare to the publised results.
void Runs(){
    LoadHEPData();
    TGraphAsymmErrors *gr_v2inte;
    double v2inte[NC];
    double centmean[NC];
    double lpt=0.2,hpt=3.0;
    TString label = Form("Pb-Pb 5.02TeV, %.1f<p_{T}<%.1f",lpt,hpt);
    for(int ic=0;ic<NC;ic++) {
        v2inte[ic] = GetIntegratedv2(ic,0.2,3.0);
        centmean[ic] = cent[ic]+(cent[ic+1]-cent[ic])/2.;
    }
    gr_v2inte = new TGraphAsymmErrors(NC,centmean,v2inte,0,0,0,0);

// v2
    Filipad *fpad3;
    fpad3 = new Filipad(1, 0.9, 0.4, 100, 100, 1.1,2);
    fpad3->Draw();
    //==== Upper pad
    TPad *p = fpad3->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    double lowx = 0,highx=60.;
    double ly=0,hy=0.2;
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "centrality[%]", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.3,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.06);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    gr_v2inte0230pub->SetMarkerStyle(20);
    gr_v2inte0230pub->Draw("psame");
    gr_v2inte->SetMarkerStyle(25);
    gr_v2inte->SetMarkerColor(2);
    gr_v2inte->Draw("psame");
    leg->AddEntry(gr_v2inte0230pub,Form("%s",label.Data()),"");
    leg->AddEntry(gr_v2inte0230pub,"published","lp");
    leg->AddEntry(gr_v2inte,"calculated from dN/dp_{T} and v_{2}(p_{T})","lp");
    leg->Draw();
  
    //==== Lower pad
    p = fpad3->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -0.2,0.2);
    hset( *hfr1, "centrality[%]","Ratio",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphAsymmErrors *gr_ratio;
    gr_ratio = GetDataOverTheory(gr_v2inte,gr_v2inte0230pub);
    gr_ratio->SetMarkerStyle(25);
    gr_ratio->SetMarkerColor(2);
    gr_ratio->Draw("psame");

    gPad->GetCanvas()->SaveAs(Form("figs/integratedv2_%.1f_%.1f.pdf",lpt,hpt));
}

void Run3allcent(){

    LoadHEPData();
    for(int ic=0;ic<NC;ic++) {
        Run3(ic);
    }
}

void Run3(int iC){
    
    TString label = Form("%2.0f-%2.0f%%",cent[iC],cent[iC+1]);
    double lowx=0., highx=10.;
    double ly=0., hy=0.30;
    int logx = 0;

    // v2
    Filipad *fpad;
    fpad = new Filipad(1, 0.9, 0.4, 100, 100, 1.1,2);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(logx); p->SetLogy(0); p->cd();
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "p_{T}(GeV/c)", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.65,0.5,0.85,0.78,label,"brNDC");
    leg->SetTextSize(0.06);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    gr_v2pt[iC]->SetMarkerStyle(20);
    gr_v2pt[iC]->Draw("psame");
    leg->AddEntry(gr_v2pt[iC],"Run2","lp");
    gr_v2ptRun3[iC]->SetMarkerStyle(24);
    gr_v2ptRun3[iC]->SetMarkerColor(2);
    gr_v2ptRun3[iC]->SetLineColor(2);
    gr_v2ptRun3[iC]->Draw("psame");
    leg->AddEntry(gr_v2ptRun3[iC],"Run3","lp");
    leg->Draw();
  

    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.5,1.5);
    hset( *hfr1, "p_{T}(GeV/c)","#frac{Run3}{Run2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    
    TGraphErrors *gr_ratio1;
    gr_ratio1 = GetRatio(gr_v2ptRun3[iC],gr_v2pt[iC]);
    gr_ratio1->Draw("same");
    gPad->GetCanvas()->SaveAs(Form("figs/v2Run3_C%d.pdf",iC));
}

