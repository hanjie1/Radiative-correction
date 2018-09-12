#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t He3_x[10][MAXBIN],He3_Q2[10][MAXBIN],He3_Yield[10][MAXBIN],He3_Yerr[10][MAXBIN];
     Double_t He3R_x[10][MAXBIN],He3R_Q2[10][MAXBIN],He3_Born1[10][MAXBIN],He3_Rad1[10][MAXBIN];
     Double_t He3_Born2[10][MAXBIN],He3_Rad2[10][MAXBIN];
     Double_t D2_x[10][MAXBIN],D2_Q2[10][MAXBIN],D2_Yield[10][MAXBIN],D2_Yerr[10][MAXBIN];
     Double_t D2R_x[10][MAXBIN],D2R_Q2[10][MAXBIN],D2_Born1[10][MAXBIN],D2_Rad1[10][MAXBIN];
     Double_t D2_Born2[10][MAXBIN],D2_Rad2[10][MAXBIN];
     Double_t He3D_ratio[10][MAXBIN],He3D_err[10][MAXBIN],D2_RC1[10][MAXBIN],He3_RC1[10][MAXBIN],D2_RC2[10][MAXBIN],He3_RC2[10][MAXBIN];

     for(int ii=0;ii<10;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Yield[ii][jj]=0.0; He3_Yerr[ii][jj]=0.0;
             He3R_x[ii][jj]=0.0; He3R_Q2[ii][jj]=0.0; He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             He3_Born2[ii][jj]=0.0; He3_Rad2[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2R_x[ii][jj]=0.0; D2R_Q2[ii][jj]=0.0; D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0;
             He3D_ratio[ii][jj]=0.0;He3D_err[ii][jj]=0.0;
             D2_RC1[ii][jj]=0.0;He3_RC1[ii][jj]=0.0;
             D2_RC2[ii][jj]=0.0;He3_RC2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<16){
       Yfile=Form("He3_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,He3_x,He3_Q2,He3_Yield,He3_Yerr); 
       RCfile=Form("gsmearing/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born1,He3_Rad1); 
       RCfile=Form("Bodek/He3_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,He3R_x,He3R_Q2,He3_Born2,He3_Rad2);  
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_Q2,D2_Yield,D2_Yerr); 
       RCfile=Form("gsmearing/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born1,D2_Rad1);  
       RCfile=Form("Bodek/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born2,D2_Rad2);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gHe3DRaw[10];
  for(int ii=0;ii<10;ii++){
      gHe3DRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||He3_x[ii][jj]==0)continue;
          He3D_ratio[ii][jj]=He3_Yield[ii][jj]/D2_Yield[ii][jj];
	  He3D_err[ii][jj]=He3D_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2));
          if(He3D_err[ii][jj]>0.1)continue;
          gHe3DRaw[ii]->SetPoint(nn,D2_x[ii][jj],He3D_ratio[ii][jj]);
          gHe3DRaw[ii]->SetPointError(nn,0.0,He3D_err[ii][jj]);
          nn++;
      }
  }
    
  TGraphErrors *gHe3DRC1[10];
  TGraphErrors *gHe3DRC2[10];
  TGraphErrors *gRatio[10];
  for(int ii=0;ii<10;ii++){
      gHe3DRC1[ii]=new TGraphErrors();
      gHe3DRC2[ii]=new TGraphErrors();
      gRatio[ii]=new TGraphErrors();
      int nn=0;
      int nn1=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||He3_x[ii][jj]==0)continue;
          double_t deltax1=D2_x[ii][jj]-D2R_x[ii][nn];
          Double_t deltax2=He3_x[ii][jj]-He3R_x[ii][nn];
          if(abs(deltax1)>0.0001&&abs(deltax2)>0.0001){nn++;continue;}
          D2_RC1[ii][nn]=D2_Born1[ii][nn]/D2_Rad1[ii][nn];
          He3_RC1[ii][nn]=He3_Born1[ii][nn]/He3_Rad1[ii][nn];
          Double_t ratio1=He3D_ratio[ii][jj]*He3_RC1[ii][nn]/D2_RC1[ii][nn];
          Double_t ratio1_err=He3D_err[ii][jj]*He3_RC1[ii][nn]/D2_RC1[ii][nn];

          D2_RC2[ii][nn]=D2_Born2[ii][nn]/D2_Rad2[ii][nn];
          He3_RC2[ii][nn]=He3_Born2[ii][nn]/He3_Rad2[ii][nn];
          Double_t ratio2=He3D_ratio[ii][jj]*He3_RC2[ii][nn]/D2_RC2[ii][nn];
          Double_t ratio2_err=He3D_err[ii][jj]*He3_RC2[ii][nn]/D2_RC2[ii][nn];
          if(He3D_err[ii][jj]>0.1){nn++;continue;}
          gHe3DRC1[ii]->SetPoint(nn1,D2_x[ii][jj],ratio1);
          gHe3DRC1[ii]->SetPointError(nn1,0.0,ratio1_err);
          gHe3DRC2[ii]->SetPoint(nn1,D2_x[ii][jj],ratio2);
          gHe3DRC2[ii]->SetPointError(nn1,0.0,ratio2_err);
          gRatio[ii]->SetPoint(nn1,D2_x[ii][jj],ratio2/ratio1);
          gRatio[ii]->SetPointError(nn1,0.0,0.0);
          nn1++;
          nn++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gHe3DRaw[ii]->SetMarkerStyle(8);
      if(ii<9) gHe3DRaw[ii]->SetMarkerColor(ii+1);
      else gHe3DRaw[ii]->SetMarkerColor(46);
      mg1->Add(gHe3DRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("He3/D Raw Yield ratio;x;He3/D");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5)leg1->AddEntry(gHe3DRaw[ii],Form("He3/D kin%d",ii+1),"P");
      else {leg1->AddEntry(gHe3DRaw[ii],Form("He3/D kin%d",ii+nn),"P");
             nn++;
           }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gHe3DRC1[ii]->SetMarkerStyle(8);
      gHe3DRC2[ii]->SetMarkerStyle(22);
      if(ii<9){
         gHe3DRC1[ii]->SetMarkerColor(ii+1);
         gHe3DRC2[ii]->SetMarkerColor(ii+1);
      }
      else {
         gHe3DRC1[ii]->SetMarkerColor(46);
         gHe3DRC2[ii]->SetMarkerColor(46);
      }
      mg2->Add(gHe3DRC1[ii]);
      mg2->Add(gHe3DRC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("He3/D with Radiative correction;x;He3/D");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5){
         leg2->AddEntry(gHe3DRC1[ii],Form("gsmearing kin%d",ii+1),"P");
         leg2->AddEntry(gHe3DRC2[ii],Form("Bodek kin%d",ii+1),"P");
      }
      else{
         leg2->AddEntry(gHe3DRC1[ii],Form("gsmearing kin%d",ii+nn),"P");
         leg2->AddEntry(gHe3DRC2[ii],Form("Bodek kin%d",ii+nn),"P");
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
  mg3->SetTitle("Ratio of two He3/D results;x;");

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
