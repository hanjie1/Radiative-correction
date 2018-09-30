#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Yield[5][MAXBIN],H1_Yerr[5][MAXBIN];
     Double_t H1R_x[5][MAXBIN],H1R_Q2[5][MAXBIN],H1_Born1[5][MAXBIN],H1_Rad1[5][MAXBIN];
     Double_t H1_Born2[5][MAXBIN],H1_Rad2[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_Yield[5][MAXBIN],D2_Yerr[5][MAXBIN];
     Double_t D2R_x[5][MAXBIN],D2R_Q2[5][MAXBIN],D2_Born1[5][MAXBIN],D2_Rad1[5][MAXBIN];
     Double_t D2_Born2[5][MAXBIN],D2_Rad2[5][MAXBIN];
     Double_t Dp_ratio[5][MAXBIN],Dp_err[5][MAXBIN],D2_RC1[5][MAXBIN],H1_RC1[5][MAXBIN],D2_RC2[5][MAXBIN],H1_RC2[5][MAXBIN];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Yield[ii][jj]=0.0; H1_Yerr[ii][jj]=0.0;
             H1R_x[ii][jj]=0.0; H1R_Q2[ii][jj]=0.0; H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
             H1_Born2[ii][jj]=0.0; H1_Rad2[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2R_x[ii][jj]=0.0; D2R_Q2[ii][jj]=0.0; D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0;
             Dp_ratio[ii][jj]=0.0;Dp_err[ii][jj]=0.0;
             D2_RC1[ii][jj]=0.0;H1_RC1[ii][jj]=0.0;
             D2_RC2[ii][jj]=0.0;H1_RC2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=0;
   while(KIN<5){
       Yfile=Form("H1_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,H1_x,H1_Q2,H1_Yield,H1_Yerr); 
       RCfile=Form("gsmearing_newbin/H1_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H1R_x,H1R_Q2,H1_Born1,H1_Rad1); 
       RCfile=Form("Bodek_newbin/H1_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,H1R_x,H1R_Q2,H1_Born2,H1_Rad2);  
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_Q2,D2_Yield,D2_Yerr); 
       RCfile=Form("gsmearing_newbin/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born1,D2_Rad1);  
       RCfile=Form("Bodek_newbin/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born2,D2_Rad2);  
       KIN+=1;
   }

  TGraphErrors *gDpRaw[5];
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H1_x[ii][jj]==0)continue;
          Dp_ratio[ii][jj]=D2_Yield[ii][jj]/H1_Yield[ii][jj];
	  Dp_err[ii][jj]=Dp_ratio[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(H1_Yerr[ii][jj]/H1_Yield[ii][jj],2));
//          if(Dp_err[ii][jj]>0.1)continue;
          gDpRaw[ii]->SetPoint(nn,D2R_x[ii][jj],Dp_ratio[ii][jj]);
          gDpRaw[ii]->SetPointError(nn,0.0,Dp_err[ii][jj]);
          nn++;
      }
  }
    
  TGraphErrors *gDpRC1[5];
  TGraphErrors *gDpRC2[5];
  TGraphErrors *gRatio[5];
  TGraph *gRCfactor1[5];
  TGraph *gRCfactor2[5];
  for(int ii=0;ii<5;ii++){
      gDpRC1[ii]=new TGraphErrors();
      gDpRC2[ii]=new TGraphErrors();
      gRatio[ii]=new TGraphErrors();
      gRCfactor1[ii]=new TGraph();
      gRCfactor2[ii]=new TGraph();
      int nn=0;
      int nn1=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H1_x[ii][jj]==0)continue;
         // Double_t deltax=D2_x[ii][jj]-D2R_x[ii][nn];
       //   if(abs(deltax)>0.0001){nn++;continue;}
          D2_RC1[ii][nn]=D2_Born1[ii][nn]/D2_Rad1[ii][nn];
          H1_RC1[ii][nn]=H1_Born1[ii][nn]/H1_Rad1[ii][nn];
          Double_t ratio1=Dp_ratio[ii][jj]*D2_RC1[ii][nn]/H1_RC1[ii][nn];
          Double_t ratio1_err=Dp_err[ii][jj]*D2_RC1[ii][nn]/H1_RC1[ii][nn];

          D2_RC2[ii][nn]=D2_Born2[ii][nn]/D2_Rad2[ii][nn];
          H1_RC2[ii][nn]=H1_Born2[ii][nn]/H1_Rad2[ii][nn];
          Double_t ratio2=Dp_ratio[ii][jj]*D2_RC2[ii][nn]/H1_RC2[ii][nn];
          Double_t ratio2_err=Dp_err[ii][jj]*D2_RC2[ii][nn]/H1_RC2[ii][nn];
          //if(Dp_err[ii][jj]>0.1){nn++;continue;}
          gDpRC1[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio1);
          gDpRC1[ii]->SetPointError(nn1,0.0,ratio1_err);
          gDpRC2[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio2);
          gDpRC2[ii]->SetPointError(nn1,0.0,ratio2_err);
          gRatio[ii]->SetPoint(nn1,D2R_x[ii][jj],ratio2/ratio1);
          gRatio[ii]->SetPointError(nn1,0.0,0.0);
          gRCfactor1[ii]->SetPoint(nn1,D2R_x[ii][jj],D2_RC1[ii][nn]/H1_RC1[ii][nn]);
          gRCfactor2[ii]->SetPoint(nn1,D2R_x[ii][jj],D2_RC2[ii][nn]/H1_RC2[ii][nn]);
          nn++;
          nn1++;
      }
  }

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gDpRaw[ii]->SetMarkerStyle(8);
      gDpRaw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gDpRaw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("D/p Raw Yield ratio;x;D/p");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg1->AddEntry(gDpRaw[ii],Form("D/p kin%d",ii),"P");
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gDpRC1[ii]->SetMarkerStyle(8);
      gDpRC1[ii]->SetMarkerColor(ii+1);
      gDpRC2[ii]->SetMarkerStyle(22);
      gDpRC2[ii]->SetMarkerColor(ii+1);
      mg2->Add(gDpRC1[ii]);
      mg2->Add(gDpRC2[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("D/p with Radiative correction;x;D/p");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg2->AddEntry(gDpRC1[ii],Form("gsmearing_newbin kin%d",ii),"P");
      leg2->AddEntry(gDpRC2[ii],Form("Bodek_newbin kin%d",ii),"P");
  }
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gRatio[ii]->SetMarkerStyle(8);
      gRatio[ii]->SetMarkerColor(ii+1);
      mg3->Add(gRatio[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("Ratio of two D/p results;x;");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg3->AddEntry(gRatio[ii],Form("Bodek_newbin/gsmearing_newbin kin%d",ii),"P");
  }
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<5;ii++){
      gRCfactor1[ii]->SetMarkerStyle(8);
      gRCfactor1[ii]->SetMarkerColor(ii+1);
      gRCfactor2[ii]->SetMarkerStyle(22);
      gRCfactor2[ii]->SetMarkerColor(ii+1);
      mg4->Add(gRCfactor1[ii]);
      mg4->Add(gRCfactor2[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("D/p radiative correction factor ratio;x;RC_D/RC_p");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<5;ii++){
      leg4->AddEntry(gRCfactor1[ii],Form("gsmearing_newbin kin%d",ii),"P");
      leg4->AddEntry(gRCfactor2[ii],Form("Bodek_newbin kin%d",ii),"P");
  }
  leg4->Draw();


}
