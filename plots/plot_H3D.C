#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3D()
{
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Yield[11][MAXBIN],H3_Yerr[11][MAXBIN];
     Double_t H3R_x[11][MAXBIN],H3R_Q2[11][MAXBIN],H3_Born1[11][MAXBIN],H3_Rad1[11][MAXBIN];
     Double_t H3_Born2[11][MAXBIN],H3_Rad2[11][MAXBIN];
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN];
     Double_t D2R_x[11][MAXBIN],D2R_Q2[11][MAXBIN],D2_Born1[11][MAXBIN],D2_Rad1[11][MAXBIN];
     Double_t D2_Born2[11][MAXBIN],D2_Rad2[11][MAXBIN];
     Double_t H3D_ratio[11][MAXBIN],H3D_err[11][MAXBIN],D2_RC1[11][MAXBIN],H3_RC1[11][MAXBIN],D2_RC2[11][MAXBIN],H3_RC2[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             H3R_x[ii][jj]=0.0; H3R_Q2[ii][jj]=0.0; H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
             H3_Born2[ii][jj]=0.0; H3_Rad2[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2R_x[ii][jj]=0.0; D2R_Q2[ii][jj]=0.0; D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0;
             H3D_ratio[ii][jj]=0.0;H3D_err[ii][jj]=0.0;
             D2_RC1[ii][jj]=0.0;H3_RC1[ii][jj]=0.0;
             D2_RC2[ii][jj]=0.0;H3_RC2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=0;
   while(KIN<16){
       Yfile=Form("H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_Q2,H3_Yield,H3_Yerr); 
       RCfile=Form("gsmearing_newbin/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born1,H3_Rad1); 
       RCfile=Form("Bodek_newbin/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born2,H3_Rad2);  
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_Q2,D2_Yield,D2_Yerr); 
       RCfile=Form("gsmearing_newbin/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born1,D2_Rad1);  
       RCfile=Form("Bodek_newbin/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born2,D2_Rad2);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3DRaw[11];
  for(int ii=0;ii<11;ii++){
      gH3DRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          H3D_ratio[ii][jj]=H3_Yield[ii][jj]/D2_Yield[ii][jj];
	  H3D_err[ii][jj]=H3D_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2));
          if(H3D_err[ii][jj]>0.1)continue;
          gH3DRaw[ii]->SetPoint(nn,D2R_x[ii][jj],H3D_ratio[ii][jj]);
          gH3DRaw[ii]->SetPointError(nn,0.0,H3D_err[ii][jj]);
          nn++;
      }
  }
    
  TGraphErrors *gH3DRC1[11];
  TGraphErrors *gH3DRC2[11];
  TGraphErrors *gRatio[11];
  TGraph *gRCfactor1[11];
  TGraph *gRCfactor2[11];
  for(int ii=0;ii<11;ii++){
      gH3DRC1[ii]=new TGraphErrors();
      gH3DRC2[ii]=new TGraphErrors();
      gRatio[ii]=new TGraphErrors();
      gRCfactor1[ii]=new TGraph();
      gRCfactor2[ii]=new TGraph();
      int nn=0;
      int nn1=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          D2_RC1[ii][nn]=D2_Born1[ii][nn]/D2_Rad1[ii][nn];
          H3_RC1[ii][nn]=H3_Born1[ii][nn]/H3_Rad1[ii][nn];
          Double_t ratio1=H3D_ratio[ii][jj]*H3_RC1[ii][nn]/D2_RC1[ii][nn];
          Double_t ratio1_err=H3D_err[ii][jj]*H3_RC1[ii][nn]/D2_RC1[ii][nn];

          D2_RC2[ii][nn]=D2_Born2[ii][nn]/D2_Rad2[ii][nn];
          H3_RC2[ii][nn]=H3_Born2[ii][nn]/H3_Rad2[ii][nn];
          Double_t ratio2=H3D_ratio[ii][jj]*H3_RC2[ii][nn]/D2_RC2[ii][nn];
          Double_t ratio2_err=H3D_err[ii][jj]*H3_RC2[ii][nn]/D2_RC2[ii][nn];
          //if(H3D_err[ii][jj]>0.1){nn++;continue;}
          gH3DRC1[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio1);
          gH3DRC1[ii]->SetPointError(nn1,0.0,ratio1_err);
          gH3DRC2[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio2);
          gH3DRC2[ii]->SetPointError(nn1,0.0,ratio2_err);
          gRatio[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio2/ratio1);
          gRatio[ii]->SetPointError(nn1,0.0,0.0);
          gRCfactor1[ii]->SetPoint(nn1,D2R_x[ii][jj],H3_RC1[ii][nn]/D2_RC1[ii][nn]);
          gRCfactor2[ii]->SetPoint(nn1,D2R_x[ii][jj],H3_RC2[ii][nn]/D2_RC2[ii][nn]);
          nn1++;
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3DRaw[ii]->SetMarkerStyle(8);
      if(ii<9) gH3DRaw[ii]->SetMarkerColor(ii+1);
      else gH3DRaw[ii]->SetMarkerColor(ii+2);
      mg1->Add(gH3DRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3/D Raw Yield ratio;x;H3/D");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg1->AddEntry(gH3DRaw[ii],Form("H3/D kin%d",ii),"P");
      else {leg1->AddEntry(gH3DRaw[ii],Form("H3/D kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3DRC1[ii]->SetMarkerStyle(8);
      gH3DRC2[ii]->SetMarkerStyle(22);
      if(ii<9){
         gH3DRC1[ii]->SetMarkerColor(ii+1);
         gH3DRC2[ii]->SetMarkerColor(ii+1);
      }
      else {
         gH3DRC1[ii]->SetMarkerColor(ii+2);
         gH3DRC2[ii]->SetMarkerColor(ii+2);
      }
      mg2->Add(gH3DRC1[ii]);
      mg2->Add(gH3DRC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H3/D with Radiative correction;x;H3/D");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6){
         leg2->AddEntry(gH3DRC1[ii],Form("gsmearing_newbin kin%d",ii),"P");
         leg2->AddEntry(gH3DRC2[ii],Form("Bodek_newbin kin%d",ii),"P");
      }
      else{
         leg2->AddEntry(gH3DRC1[ii],Form("gsmearing_newbin kin%d",ii+nn),"P");
         leg2->AddEntry(gH3DRC2[ii],Form("Bodek_newbin kin%d",ii+nn),"P");
	 nn++;
      }
  }
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gRatio[ii]->SetMarkerStyle(8);
      if(ii<9)gRatio[ii]->SetMarkerColor(ii+1);
      else gRatio[ii]->SetMarkerColor(ii+2);
      mg3->Add(gRatio[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("Ratio of two H3/D results;x;");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg3->AddEntry(gRatio[ii],Form("Bodek_newbin/gsmearing_newbin kin%d",ii),"P");
      else {leg3->AddEntry(gRatio[ii],Form("Bodek_newbin/gsmearing_newbin kin%d",ii+nn),"P");
            nn++;
           }
  }
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gRCfactor1[ii]->SetMarkerStyle(8);
      gRCfactor2[ii]->SetMarkerStyle(22);
      if(ii<9){
         gRCfactor1[ii]->SetMarkerColor(ii+1);
         gRCfactor2[ii]->SetMarkerColor(ii+1);
      }
      else {
         gRCfactor1[ii]->SetMarkerColor(ii+2);
         gRCfactor2[ii]->SetMarkerColor(ii+2);
      }
      mg4->Add(gRCfactor1[ii]);
      mg4->Add(gRCfactor2[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("H3/D with Radiative correction;x;H3/D");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6){
         leg4->AddEntry(gRCfactor1[ii],Form("gsmearing_newbin kin%d",ii),"P");
         leg4->AddEntry(gRCfactor2[ii],Form("Bodek_newbin kin%d",ii),"P");
      }
      else{
         leg4->AddEntry(gRCfactor1[ii],Form("gsmearing_newbin kin%d",ii+nn),"P");
         leg4->AddEntry(gRCfactor2[ii],Form("Bodek_newbin kin%d",ii+nn),"P");
	 nn++;
      }
  }
  leg4->Draw();




}
