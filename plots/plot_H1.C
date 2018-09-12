#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H1()
{
     Double_t H1_x[4][MAXBIN],H1_Q2[4][MAXBIN],H1_Yield[4][MAXBIN],H1_Yerr[4][MAXBIN];
     Double_t H1R_x[4][MAXBIN],H1R_Q2[4][MAXBIN],H1_Born1[4][MAXBIN],H1_Rad1[4][MAXBIN];
     Double_t H1_Born2[4][MAXBIN],H1_Rad2[4][MAXBIN],H1_Born3[4][MAXBIN],H1_Rad3[4][MAXBIN];

     for(int ii=0;ii<4;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
		H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Yield[ii][jj]=0.0; H1_Yerr[ii][jj]=0.0;
                H1R_x[ii][jj]=0.0; H1R_Q2[ii][jj]=0.0; H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
                H1_Born2[ii][jj]=0.0; H1_Rad2[ii][jj]=0.0;H1_Born3[ii][jj]=0.0; H1_Rad3[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<5){
          Yfile=Form("H1_kin%d.txt",KIN);
          ReadYield(Yfile,KIN,H1_x,H1_Q2,H1_Yield,H1_Yerr); 
          RCfile=Form("F2GLOB/H1_kin%d_xs.out",KIN);
          ReadRC(RCfile,KIN,H1R_x,H1R_Q2,H1_Born1,H1_Rad1); 
          RCfile=Form("Bodek/H1_kin%d_xs.out",KIN);
          ReadRC(RCfile,KIN,H1R_x,H1R_Q2,H1_Born2,H1_Rad2); 
          RCfile=Form("gsmearing/H1_kin%d_xs.out",KIN);
          ReadRC(RCfile,KIN,H1R_x,H1R_Q2,H1_Born3,H1_Rad3); 
          KIN+=1;
   }

  TGraphErrors *gH1Raw[4];
  for(int ii=0;ii<4;ii++){
      gH1Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H1_x[ii][jj]==0)continue;
          gH1Raw[ii]->SetPoint(nn,H1_x[ii][jj],H1_Yield[ii][jj]);
          gH1Raw[ii]->SetPointError(nn,0.0,H1_Yerr[ii][jj]);
          nn++;
      }
  }
  TGraph *gH1Rad1[4];
  TGraph *gH1Born1[4];
  TGraph *gH1RC1[4];
  TGraph *gH1Rad2[4];
  TGraph *gH1Born2[4];
  TGraph *gH1RC2[4];
  TGraph *gH1Rad3[4];
  TGraph *gH1Born3[4];
  TGraph *gH1RC3[4];
  TGraph *gRatio_1[4];
  TGraph *gRatio_2[4];
  for(int ii=0;ii<4;ii++){
      gH1Rad1[ii]=new TGraph();
      gH1Born1[ii]=new TGraph();
      gH1RC1[ii]=new TGraph();
      gH1Rad2[ii]=new TGraph();
      gH1Born2[ii]=new TGraph();
      gH1RC2[ii]=new TGraph();
      gH1Rad3[ii]=new TGraph();
      gH1Born3[ii]=new TGraph();
      gH1RC3[ii]=new TGraph();
      gRatio_1[ii]=new TGraph();
      gRatio_2[ii]=new TGraph();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(H1R_x[ii][jj]==0)continue;
          gH1Rad1[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Rad1[ii][jj]);
          gH1Born1[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born1[ii][jj]);
          gH1RC1[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born1[ii][jj]/H1_Rad1[ii][jj]);
          gH1Rad2[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Rad2[ii][jj]);
          gH1Born2[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born2[ii][jj]);
          gH1RC2[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born2[ii][jj]/H1_Rad2[ii][jj]);
          gH1Rad3[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Rad3[ii][jj]);
          gH1Born3[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born3[ii][jj]);
          gH1RC3[ii]->SetPoint(nn,H1R_x[ii][jj],H1_Born3[ii][jj]/H1_Rad3[ii][jj]);
          gRatio_1[ii]->SetPoint(nn,H1R_x[ii][jj],(H1_Born1[ii][jj]/H1_Rad1[ii][jj])/(H1_Born2[ii][jj]/H1_Rad2[ii][jj]));
          gRatio_2[ii]->SetPoint(nn,H1R_x[ii][jj],(H1_Born2[ii][jj]/H1_Rad2[ii][jj])/(H1_Born3[ii][jj]/H1_Rad3[ii][jj]));
          nn++;
      }
  }
    

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gH1Raw[ii]->SetMarkerStyle(8);
      gH1Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gH1Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  for(int ii=0;ii<4;ii++){
      leg1->AddEntry(gH1Raw[ii],Form("H1 kin%d",ii+1),"P");
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gH1RC1[ii]->SetMarkerStyle(8);
      gH1RC1[ii]->SetMarkerColor(1);
      mg2->Add(gH1RC1[ii]);
      gH1RC2[ii]->SetMarkerStyle(22);
      gH1RC2[ii]->SetMarkerColor(3);
      mg2->Add(gH1RC2[ii]);
      gH1RC3[ii]->SetMarkerStyle(33);
      gH1RC3[ii]->SetMarkerColor(4);
      mg2->Add(gH1RC3[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("H1 RC factor;x;Born/Rad");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  leg2->AddEntry(gH1RC1[0],"F2GLOB","P");
  leg2->AddEntry(gH1RC2[0],"INEFT","P");
  leg2->AddEntry(gH1RC3[0],"gsmearing","P");
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gH1Born1[ii]->SetMarkerStyle(8);
      gH1Born1[ii]->SetMarkerColor(1);
      mg3->Add(gH1Born1[ii]);
      gH1Born2[ii]->SetMarkerStyle(22);
      gH1Born2[ii]->SetMarkerColor(3);
      mg3->Add(gH1Born2[ii]);
      gH1Born3[ii]->SetMarkerStyle(33);
      gH1Born3[ii]->SetMarkerColor(4);
      mg3->Add(gH1Born3[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("H1 Born cross section;x;Born XS");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  leg3->AddEntry(gH1Born1[0],"F2GLOB","P");
  leg3->AddEntry(gH1Born2[0],"INEFT","P");
  leg3->AddEntry(gH1Born3[0],"gsmearing","P");
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gH1Rad1[ii]->SetMarkerStyle(8);
      gH1Rad1[ii]->SetMarkerColor(1);
      mg4->Add(gH1Rad1[ii]);
      gH1Rad2[ii]->SetMarkerStyle(22);
      gH1Rad2[ii]->SetMarkerColor(3);
      mg4->Add(gH1Rad2[ii]);
      gH1Rad3[ii]->SetMarkerStyle(33);
      gH1Rad3[ii]->SetMarkerColor(4);
      mg4->Add(gH1Rad3[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("H1 Radiative cross section;x;Rad XS");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  leg4->AddEntry(gH1Rad1[0],"F2GLOB","P");
  leg4->AddEntry(gH1Rad2[0],"INEFT","P");
  leg4->AddEntry(gH1Rad3[0],"gsmearing","P");
  leg4->Draw();

  TCanvas *c5=new TCanvas("c5");
  c5->Divide(2,1);
  c5->cd(1);
  TMultiGraph *mg5=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gRatio_1[ii]->SetMarkerStyle(8);
      gRatio_1[ii]->SetMarkerColor(1);
      mg5->Add(gRatio_1[ii]);
  }
  mg5->Draw("AP");
  mg5->SetTitle("H1 RC factor ratio;x;");

  auto leg5=new TLegend(0.65,0.7,0.85,0.85);
  leg5->AddEntry(gRatio_1[0],"F2GLOB/INEFT","P");
  leg5->Draw();
  c5->cd(2);
  TMultiGraph *mg6=new TMultiGraph();
  for(int ii=0;ii<4;ii++){
      gRatio_2[ii]->SetMarkerStyle(8);
      gRatio_2[ii]->SetMarkerColor(1);
      mg6->Add(gRatio_2[ii]);
  }
  mg6->Draw("AP");
  mg6->SetTitle("RC factor ratio;x;");

  auto leg6=new TLegend(0.65,0.7,0.85,0.85);
  leg6->AddEntry(gRatio_2[0],"INEFT/gsmearing","P");
  leg6->Draw();




/*
  for(int ii=0;ii<10;ii++){
     cout<<"-------  "<<ii<<"  -------"<<endl;
     for(int jj=0;jj<MAXBIN;jj++){
         if(H1_x[ii][jj]==0)continue;
         //cout<<H1R_x[ii][jj]<<", "<<H1R_Q2[ii][jj]<<", "<<H1_Born[ii][jj]<<", "<<H1_Rad[ii][jj]<<endl;
         cout<<H1_x[ii][jj]<<", "<<H1_Q2[ii][jj]<<", "<<H1_Yield[ii][jj]<<", "<<H1_Yerr[ii][jj]<<endl;
     }
  }
*/


}
