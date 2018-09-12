#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He3()
{
     Double_t H3_x[10][MAXBIN],H3_Q2[10][MAXBIN],H3_Yield[10][MAXBIN],H3_Yerr[10][MAXBIN];
     Double_t H3R_x[10][MAXBIN],H3R_Q2[10][MAXBIN],H3_Born1[10][MAXBIN],H3_Rad1[10][MAXBIN];
     Double_t H3_Born2[10][MAXBIN],H3_Rad2[10][MAXBIN];
     Double_t He3_x[10][MAXBIN],He3_Q2[10][MAXBIN],He3_Yield[10][MAXBIN],He3_Yerr[10][MAXBIN];
     Double_t He3R_x[10][MAXBIN],He3R_Q2[10][MAXBIN],He3_Born1[10][MAXBIN],He3_Rad1[10][MAXBIN];
     Double_t He3_Born2[10][MAXBIN],He3_Rad2[10][MAXBIN];
     Double_t H3He3_ratio[10][MAXBIN],H3He3_err[10][MAXBIN],He3_RC1[10][MAXBIN],H3_RC1[10][MAXBIN],He3_RC2[10][MAXBIN],H3_RC2[10][MAXBIN];

     for(int ii=0;ii<10;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Yield[ii][jj]=0.0; H3_Yerr[ii][jj]=0.0;
             H3R_x[ii][jj]=0.0; H3R_Q2[ii][jj]=0.0; H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
             H3_Born2[ii][jj]=0.0; H3_Rad2[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             He3R_x[ii][jj]=0.0; He3R_Q2[ii][jj]=0.0; He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             He3_Born2[ii][jj]=0.0; He3_Rad2[ii][jj]=0.0;
             H3He3_ratio[ii][jj]=0.0;H3He3_err[ii][jj]=0.0;
             He3_RC1[ii][jj]=0.0;H3_RC1[ii][jj]=0.0;
             He3_RC2[ii][jj]=0.0;H3_RC2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<16){
       Yfile=Form("H3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H3_x,H3_Q2,H3_Yield,H3_Yerr); 
       RCfile=Form("gsmearing/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born1,H3_Rad1); 
       RCfile=Form("Bodek/H3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H3R_x,H3R_Q2,H3_Born2,H3_Rad2);  
       Yfile=Form("He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_Q2,He3_Yield,He3_Yerr); 
       RCfile=Form("gsmearing/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born1,He3_Rad1);  
       RCfile=Form("Bodek/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born2,He3_Rad2);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gH3He3Raw[10];
  for(int ii=0;ii<10;ii++){
      gH3He3Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          H3He3_ratio[ii][jj]=H3_Yield[ii][jj]/He3_Yield[ii][jj];
	  H3He3_err[ii][jj]=H3He3_ratio[ii][jj]*sqrt(pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2)+pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2));
          if(H3He3_err[ii][jj]>0.1)continue;
          gH3He3Raw[ii]->SetPoint(nn,He3_x[ii][jj],H3He3_ratio[ii][jj]);
          gH3He3Raw[ii]->SetPointError(nn,0.0,H3He3_err[ii][jj]);
          nn++;
      }
  }
    
  TGraphErrors *gH3He3RC1[10];
  TGraphErrors *gH3He3RC2[10];
  TGraphErrors *gRatio[10];
  for(int ii=0;ii<10;ii++){
      gH3He3RC1[ii]=new TGraphErrors();
      gH3He3RC2[ii]=new TGraphErrors();
      gRatio[ii]=new TGraphErrors();
      int nn=0;
      int nn1=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||H3_x[ii][jj]==0)continue;
          double_t deltax1=He3_x[ii][jj]-He3R_x[ii][nn];
          Double_t deltax2=H3_x[ii][jj]-H3R_x[ii][nn];
          if(abs(deltax1)>0.0001&&abs(deltax2)>0.0001){nn++;continue;}
          He3_RC1[ii][nn]=He3_Born1[ii][nn]/He3_Rad1[ii][nn];
          H3_RC1[ii][nn]=H3_Born1[ii][nn]/H3_Rad1[ii][nn];
          Double_t ratio1=H3He3_ratio[ii][jj]*H3_RC1[ii][nn]/He3_RC1[ii][nn];
          Double_t ratio1_err=H3He3_err[ii][jj]*H3_RC1[ii][nn]/He3_RC1[ii][nn];

          He3_RC2[ii][nn]=He3_Born2[ii][nn]/He3_Rad2[ii][nn];
          H3_RC2[ii][nn]=H3_Born2[ii][nn]/H3_Rad2[ii][nn];
          Double_t ratio2=H3He3_ratio[ii][jj]*H3_RC2[ii][nn]/He3_RC2[ii][nn];
          Double_t ratio2_err=H3He3_err[ii][jj]*H3_RC2[ii][nn]/He3_RC2[ii][nn];
          if(H3He3_err[ii][jj]>0.1){nn++;continue;}
          gH3He3RC1[ii]->SetPoint(nn1,He3_x[ii][jj],ratio1);
          gH3He3RC1[ii]->SetPointError(nn1,0.0,ratio1_err);
          gH3He3RC2[ii]->SetPoint(nn1,He3_x[ii][jj],ratio2);
          gH3He3RC2[ii]->SetPointError(nn1,0.0,ratio2_err);
          gRatio[ii]->SetPoint(nn1,He3_x[ii][jj],ratio2/ratio1);
          gRatio[ii]->SetPointError(nn1,0.0,0.0);
          nn1++;
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3He3Raw[ii]->SetMarkerStyle(8);
      if(ii<9) gH3He3Raw[ii]->SetMarkerColor(ii+1);
      else gH3He3Raw[ii]->SetMarkerColor(46);
      mg1->Add(gH3He3Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("H3/He3 Raw Yield ratio;x;H3/He3");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5)leg1->AddEntry(gH3He3Raw[ii],Form("H3/He3 kin%d",ii+1),"P");
      else {leg1->AddEntry(gH3He3Raw[ii],Form("H3/He3 kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gH3He3RC1[ii]->SetMarkerStyle(8);
      gH3He3RC2[ii]->SetMarkerStyle(22);
      if(ii<9){
         gH3He3RC1[ii]->SetMarkerColor(ii+1);
         gH3He3RC2[ii]->SetMarkerColor(ii+1);
      }
      else {
         gH3He3RC1[ii]->SetMarkerColor(46);
         gH3He3RC2[ii]->SetMarkerColor(46);
      }
      mg2->Add(gH3He3RC1[ii]);
      mg2->Add(gH3He3RC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H3/He3 with Radiative correction;x;H3/He3");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5){
         leg2->AddEntry(gH3He3RC1[ii],Form("gsmearing kin%d",ii+1),"P");
         leg2->AddEntry(gH3He3RC2[ii],Form("Bodek kin%d",ii+1),"P");
      }
      else{
         leg2->AddEntry(gH3He3RC1[ii],Form("gsmearing kin%d",ii+nn),"P");
         leg2->AddEntry(gH3He3RC2[ii],Form("Bodek kin%d",ii+nn),"P");
	 nn++;
      }
  }
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gRatio[ii]->SetMarkerStyle(8);
      if(ii<9)gRatio[ii]->SetMarkerColor(ii+1);
      else gRatio[ii]->SetMarkerColor(46);
      mg3->Add(gRatio[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("Ratio of two H3/He3 results;x;");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5)leg3->AddEntry(gRatio[ii],Form("Bodek/gsmearing kin%d",ii+1),"P");
      else {leg3->AddEntry(gRatio[ii],Form("Bodek/gsmearing kin%d",ii+nn),"P");
            nn++;
           }
  }
  leg3->Draw();


}
