#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_D2()
{
     Double_t D2_x[10][MAXBIN],D2_Q2[10][MAXBIN],D2_Yield[10][MAXBIN],D2_Yerr[10][MAXBIN];
     Double_t D2R_x[10][MAXBIN],D2R_Q2[10][MAXBIN],D2_Born1[10][MAXBIN],D2_Rad1[10][MAXBIN];
     Double_t D2_Born2[10][MAXBIN],D2_Rad2[10][MAXBIN],D2_Born3[10][MAXBIN],D2_Rad3[10][MAXBIN];

     for(int ii=0;ii<10;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2R_x[ii][jj]=0.0; D2R_Q2[ii][jj]=0.0; D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0;D2_Born3[ii][jj]=0.0; D2_Rad3[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<16){
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_Q2,D2_Yield,D2_Yerr); 
       RCfile=Form("F2GLOB/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born1,D2_Rad1);  
       RCfile=Form("Bodek/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born2,D2_Rad2);  
       RCfile=Form("gsmearing/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born3,D2_Rad3);  
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

  TGraphErrors *gD2Raw[10];
  for(int ii=0;ii<10;ii++){
      gD2Raw[ii]=new TGraphErrors();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0)continue;
          gD2Raw[ii]->SetPoint(nn,D2_x[ii][jj],D2_Yield[ii][jj]);
          gD2Raw[ii]->SetPointError(nn,0.0,D2_Yerr[ii][jj]);
          nn++;
      }
  }
  TGraph *gD2Rad1[10];
  TGraph *gD2Born1[10];
  TGraph *gD2RC1[10];
  TGraph *gD2Rad2[10];
  TGraph *gD2Born2[10];
  TGraph *gD2RC2[10];
  TGraph *gD2Rad3[10];
  TGraph *gD2Born3[10];
  TGraph *gD2RC3[10];
  TGraph *gRatio_1[10];
  TGraph *gRatio_2[10];
  for(int ii=0;ii<10;ii++){
      gD2Rad1[ii]=new TGraph();
      gD2Born1[ii]=new TGraph();
      gD2RC1[ii]=new TGraph();
      gD2Rad2[ii]=new TGraph();
      gD2Born2[ii]=new TGraph();
      gD2RC2[ii]=new TGraph();
      gD2Rad3[ii]=new TGraph();
      gD2Born3[ii]=new TGraph();
      gD2RC3[ii]=new TGraph();
      gRatio_1[ii]=new TGraph();
      gRatio_2[ii]=new TGraph();
      int nn=0;
      for(int jj=0;jj<MAXBIN;jj++){
          if(D2R_x[ii][jj]==0)continue;
          gD2Rad1[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Rad1[ii][jj]);
          gD2Born1[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born1[ii][jj]);
          gD2RC1[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born1[ii][jj]/D2_Rad1[ii][jj]);
          gD2Rad2[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Rad2[ii][jj]);
          gD2Born2[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born2[ii][jj]);
          gD2RC2[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born2[ii][jj]/D2_Rad2[ii][jj]);
          gD2Rad3[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Rad3[ii][jj]);
          gD2Born3[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born3[ii][jj]);
          gD2RC3[ii]->SetPoint(nn,D2R_x[ii][jj],D2_Born3[ii][jj]/D2_Rad3[ii][jj]);
          gRatio_1[ii]->SetPoint(nn,D2R_x[ii][jj],(D2_Born1[ii][jj]/D2_Rad1[ii][jj])/(D2_Born2[ii][jj]/D2_Rad2[ii][jj]));
          gRatio_2[ii]->SetPoint(nn,D2R_x[ii][jj],(D2_Born2[ii][jj]/D2_Rad2[ii][jj])/(D2_Born3[ii][jj]/D2_Rad3[ii][jj]));
          nn++;
      }
  }
    

  TCanvas *c1=new TCanvas("c1");
  TMultiGraph *mg1=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gD2Raw[ii]->SetMarkerStyle(8);
      if(ii==9)gD2Raw[ii]->SetMarkerColor(30);
      else gD2Raw[ii]->SetMarkerColor(ii+1);
      mg1->Add(gD2Raw[ii]);
  }
  mg1->Draw("AP");
  mg1->SetTitle("Data Yield;x;Yield");

  auto leg1=new TLegend(0.7,0.6,0.85,0.85);
  int nn=2;
  for(int ii=0;ii<10;ii++){
      if(ii<5)leg1->AddEntry(gD2Raw[ii],Form("D2 kin%d",ii+1),"P");
      else {leg1->AddEntry(gD2Raw[ii],Form("D2 kin%d",ii+nn),"P");
	    nn++;	    
      }
  }
  leg1->Draw();

  TCanvas *c2=new TCanvas("c2");
  TMultiGraph *mg2=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gD2RC1[ii]->SetMarkerStyle(8);
      gD2RC1[ii]->SetMarkerColor(1);
      mg2->Add(gD2RC1[ii]);
      gD2RC2[ii]->SetMarkerStyle(22);
      gD2RC2[ii]->SetMarkerColor(3);
      mg2->Add(gD2RC2[ii]);
      gD2RC3[ii]->SetMarkerStyle(33);
      gD2RC3[ii]->SetMarkerColor(4);
      mg2->Add(gD2RC3[ii]);
  }
  mg2->Draw("AP");
  mg2->SetTitle("RC factor;x;Born/Rad");

  auto leg2=new TLegend(0.7,0.6,0.85,0.85);
  leg2->AddEntry(gD2RC1[0],"F2GLOB","P");
  leg2->AddEntry(gD2RC2[0],"INEFT","P");
  leg2->AddEntry(gD2RC3[0],"gsmearing","P");
  leg2->Draw();

  TCanvas *c3=new TCanvas("c3");
  TMultiGraph *mg3=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gD2Born1[ii]->SetMarkerStyle(8);
      gD2Born1[ii]->SetMarkerColor(1);
      mg3->Add(gD2Born1[ii]);
      gD2Born2[ii]->SetMarkerStyle(22);
      gD2Born2[ii]->SetMarkerColor(3);
      mg3->Add(gD2Born2[ii]);
      gD2Born3[ii]->SetMarkerStyle(33);
      gD2Born3[ii]->SetMarkerColor(4);
      mg3->Add(gD2Born3[ii]);
  }
  mg3->Draw("AP");
  mg3->SetTitle("Born cross section;x;Born XS");

  auto leg3=new TLegend(0.7,0.6,0.85,0.85);
  leg3->AddEntry(gD2Born1[0],"F2GLOB","P");
  leg3->AddEntry(gD2Born2[0],"INEFT","P");
  leg3->AddEntry(gD2Born3[0],"gsmearing","P");
  leg3->Draw();

  TCanvas *c4=new TCanvas("c4");
  TMultiGraph *mg4=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gD2Rad1[ii]->SetMarkerStyle(8);
      gD2Rad1[ii]->SetMarkerColor(1);
      mg4->Add(gD2Rad1[ii]);
      gD2Rad2[ii]->SetMarkerStyle(22);
      gD2Rad2[ii]->SetMarkerColor(3);
      mg4->Add(gD2Rad2[ii]);
      gD2Rad3[ii]->SetMarkerStyle(33);
      gD2Rad3[ii]->SetMarkerColor(4);
      mg4->Add(gD2Rad3[ii]);
  }
  mg4->Draw("AP");
  mg4->SetTitle("Radiative cross section;x;Rad XS");

  auto leg4=new TLegend(0.7,0.6,0.85,0.85);
  leg4->AddEntry(gD2Rad1[0],"F2GLOB","P");
  leg4->AddEntry(gD2Rad2[0],"INEFT","P");
  leg4->AddEntry(gD2Rad3[0],"gsmearing","P");
  leg4->Draw();

  TCanvas *c5=new TCanvas("c5");
  c5->Divide(2,1);
  c5->cd(1);
  TMultiGraph *mg5=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
      gRatio_1[ii]->SetMarkerStyle(8);
      gRatio_1[ii]->SetMarkerColor(1);
      mg5->Add(gRatio_1[ii]);
  }
  mg5->Draw("AP");
  mg5->SetTitle("RC factor ratio;x;");

  auto leg5=new TLegend(0.65,0.7,0.85,0.85);
  leg5->AddEntry(gRatio_1[0],"F2GLOB/INEFT","P");
  leg5->Draw();
  c5->cd(2);
  TMultiGraph *mg6=new TMultiGraph();
  for(int ii=0;ii<10;ii++){
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
         if(D2_x[ii][jj]==0)continue;
         //cout<<D2R_x[ii][jj]<<", "<<D2R_Q2[ii][jj]<<", "<<D2_Born[ii][jj]<<", "<<D2_Rad[ii][jj]<<endl;
         cout<<D2_x[ii][jj]<<", "<<D2_Q2[ii][jj]<<", "<<D2_Yield[ii][jj]<<", "<<D2_Yerr[ii][jj]<<endl;
     }
  }
*/


}
