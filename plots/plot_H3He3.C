#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He3()
{
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Yield[11][MAXBIN],H3_Yerr[11][MAXBIN];
     Double_t H3R_x[11][MAXBIN],H3R_Q2[11][MAXBIN],H3_Born1[11][MAXBIN],H3_Rad1[11][MAXBIN];
     Double_t H3_Born2[11][MAXBIN],H3_Rad2[11][MAXBIN];
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN];
     Double_t He3R_x[11][MAXBIN],He3R_Q2[11][MAXBIN],He3_Born1[11][MAXBIN],He3_Rad1[11][MAXBIN];
     Double_t He3_Born2[11][MAXBIN],He3_Rad2[11][MAXBIN];
     Double_t H3He_ratio[11][MAXBIN],H3He_err[11][MAXBIN],He3_RC1[11][MAXBIN],H3_RC1[11][MAXBIN],He3_RC2[11][MAXBIN],H3_RC2[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             H3R_x[ii][jj]=0.0; H3R_Q2[ii][jj]=0.0; H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
             H3_Born2[ii][jj]=0.0; H3_Rad2[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             He3R_x[ii][jj]=0.0; He3R_Q2[ii][jj]=0.0; He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             He3_Born2[ii][jj]=0.0; He3_Rad2[ii][jj]=0.0;
             H3He_ratio[ii][jj]=0.0;H3He_err[ii][jj]=0.0;
             He3_RC1[ii][jj]=0.0;H3_RC1[ii][jj]=0.0;
             He3_RC2[ii][jj]=0.0;H3_RC2[ii][jj]=0.0;
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
       Yfile=Form("He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_Q2,He3_Yield,He3_Yerr); 
       RCfile=Form("gsmearing_newbin/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born1,He3_Rad1);  
       RCfile=Form("Bodek_newbin/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born2,He3_Rad2);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3HeRaw[11];
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          H3He_ratio[ii][jj]=H3_Yield[ii][jj]/He3_Yield[ii][jj];
	  H3He_err[ii][jj]=H3He_ratio[ii][jj]*sqrt(pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2)+pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2));
          if(H3He_err[ii][jj]>0.1)continue;
          gH3HeRaw[ii]->SetPoint(nn,He3R_x[ii][jj],H3He_ratio[ii][jj]);
          gH3HeRaw[ii]->SetPointError(nn,0.0,H3He_err[ii][jj]);
          nn++;
      }
  }
    
  TGraphErrors *gH3HeRC1[11];
  TGraphErrors *gH3HeRC2[11];
  TGraphErrors *gRatio[11];
  TGraph *gRCfactor1[11];
  TGraph *gRCfactor2[11];
  for(int ii=0;ii<11;ii++){
      gH3HeRC1[ii]=new TGraphErrors();
      gH3HeRC2[ii]=new TGraphErrors();
      gRatio[ii]=new TGraphErrors();
      gRCfactor1[ii]=new TGraph();
      gRCfactor2[ii]=new TGraph();
      int nn=0;
      int nn1=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          He3_RC1[ii][nn]=He3_Born1[ii][nn]/He3_Rad1[ii][nn];
          H3_RC1[ii][nn]=H3_Born1[ii][nn]/H3_Rad1[ii][nn];
          Double_t ratio1=H3He_ratio[ii][jj]*H3_RC1[ii][nn]/He3_RC1[ii][nn];
          Double_t ratio1_err=H3He_err[ii][jj]*H3_RC1[ii][nn]/He3_RC1[ii][nn];

          He3_RC2[ii][nn]=He3_Born2[ii][nn]/He3_Rad2[ii][nn];
          H3_RC2[ii][nn]=H3_Born2[ii][nn]/H3_Rad2[ii][nn];
          Double_t ratio2=H3He_ratio[ii][jj]*H3_RC2[ii][nn]/He3_RC2[ii][nn];
          Double_t ratio2_err=H3He_err[ii][jj]*H3_RC2[ii][nn]/He3_RC2[ii][nn];
          //if(H3He_err[ii][jj]>0.1){nn++;continue;}
          gH3HeRC1[ii]->SetPoint(nn1,He3R_x[ii][jj],ratio1);
          gH3HeRC1[ii]->SetPointError(nn1,0.0,ratio1_err);
          gH3HeRC2[ii]->SetPoint(nn1,He3R_x[ii][jj],ratio2);
          gH3HeRC2[ii]->SetPointError(nn1,0.0,ratio2_err);
          gRatio[ii]->SetPoint(nn1,He3R_x[ii][jj],ratio2/ratio1);
          gRatio[ii]->SetPointError(nn1,0.0,0.0);
          gRCfactor1[ii]->SetPoint(nn1,He3R_x[ii][jj],H3_RC1[ii][nn]/He3_RC1[ii][nn]);
          gRCfactor2[ii]->SetPoint(nn1,He3R_x[ii][jj],H3_RC2[ii][nn]/He3_RC2[ii][nn]);
          nn1++;
          nn++;
      }
  }

  int color[11]={1,2,3,4,6,7,8,9,12,30,46};

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3HeRaw[ii]->SetMarkerStyle(8);
      gH3HeRaw[ii]->SetMarkerColor(color[ii]);
      mg1->Add(gH3HeRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3/He3 Raw Yield ratio;x;H3/He3");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=1;
  leg1->SetNColumns(5);
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg1->AddEntry(gH3HeRaw[ii],Form("H3/He3 kin%d",ii),"P");
      else {leg1->AddEntry(gH3HeRaw[ii],Form("H3/He3 kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gH3HeRC1[ii]->SetMarkerStyle(8);
      gH3HeRC2[ii]->SetMarkerStyle(22);
      gH3HeRC1[ii]->SetMarkerColor(color[ii]);
      gH3HeRC2[ii]->SetMarkerColor(color[ii]);
      mg2->Add(gH3HeRC1[ii]);
      mg2->Add(gH3HeRC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H3/He3 with Radiative correction;x;H3/He3");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  leg2->SetNColumns(5);
  for(int ii=0;ii<11;ii++){
      if(ii<6){
         leg2->AddEntry(gH3HeRC1[ii],Form("gsmearing kin%d",ii),"P");
         leg2->AddEntry(gH3HeRC2[ii],Form("Bodek kin%d",ii),"P");
      }
      else{
         leg2->AddEntry(gH3HeRC1[ii],Form("gsmearing kin%d",ii+nn),"P");
         leg2->AddEntry(gH3HeRC2[ii],Form("Bodek kin%d",ii+nn),"P");
	 nn++;
      }
  }
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gRatio[ii]->SetMarkerStyle(8);
      gRatio[ii]->SetMarkerColor(color[ii]);
      mg3->Add(gRatio[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("Ratio of two H3/He3 results;x;");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  nn=1;
  leg3->SetNColumns(5);
  for(int ii=0;ii<11;ii++){
      if(ii<6)leg3->AddEntry(gRatio[ii],Form("Bodek/gsmearing kin%d",ii),"P");
      else {leg3->AddEntry(gRatio[ii],Form("Bodek/gsmearing kin%d",ii+nn),"P");
            nn++;
           }
  }
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<11;ii++){
      gRCfactor1[ii]->SetMarkerStyle(8);
      gRCfactor2[ii]->SetMarkerStyle(22);
      gRCfactor1[ii]->SetMarkerColor(color[ii]);
      gRCfactor2[ii]->SetMarkerColor(color[ii]);
      mg4->Add(gRCfactor1[ii]);
      mg4->Add(gRCfactor2[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("H3/He3 RC factor ratio;x;RC_H3/RC_He3");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  leg4->SetNColumns(5);
  nn=1;
  for(int ii=0;ii<11;ii++){
      if(ii<6){
         leg4->AddEntry(gRCfactor1[ii],Form("gsmearing kin%d",ii),"P");
         leg4->AddEntry(gRCfactor2[ii],Form("Bodek kin%d",ii),"P");
      }
      else{
         leg4->AddEntry(gRCfactor1[ii],Form("gsmearing kin%d",ii+nn),"P");
         leg4->AddEntry(gRCfactor2[ii],Form("Bodek kin%d",ii+nn),"P");
	 nn++;
      }
  }
  leg4->Draw();




}
