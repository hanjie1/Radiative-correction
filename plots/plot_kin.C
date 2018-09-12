#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_kin()
{
     Double_t D2_x[10][MAXBIN],D2_Q2[10][MAXBIN],D2_Yield[10][MAXBIN],D2_Yerr[10][MAXBIN];
     Double_t D2R_x[10][MAXBIN],D2R_Q2[10][MAXBIN],D2_Born1[10][MAXBIN],D2_Rad1[10][MAXBIN];
     Double_t D2_Born2[10][MAXBIN],D2_Rad2[10][MAXBIN];

     for(int ii=0;ii<10;ii++){
         for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
             D2R_x[ii][jj]=0.0; D2R_Q2[ii][jj]=0.0; D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0;
     }}

   TString Yfile;
   TString RCfile;
   Int_t KIN=1;
   while(KIN<16){
       Yfile=Form("D2_kin%d.txt",KIN);
       ReadYield(Yfile,KIN,D2_x,D2_Q2,D2_Yield,D2_Yerr);
       RCfile=Form("gsmearing/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born1,D2_Rad1);
       RCfile=Form("Bodek/D2_kin%d_xs.out",KIN);
       ReadRC(RCfile,KIN,D2R_x,D2R_Q2,D2_Born2,D2_Rad2);
       if(KIN<5)KIN+=1;
       else KIN+=2;
   }

   TGraphErrors *gRatio[10];
   TCanvas *myCanvas[10];
   for(int ii=0;ii<10;ii++){
	myCanvas[ii]=new TCanvas(Form("c%d",ii+1));
        gRatio[ii]=new TGraphErrors();
        int nn=0;
        for(int jj=0;jj<MAXBIN;jj++){
            if(D2_x[ii][jj]==0)continue;
            gRatio[ii]->SetPoint(nn,D2_x[ii][jj],D2_Yield[ii][jj]);//D2_Rad1[ii][nn]);
            gRatio[ii]->SetPointError(nn,0.0,D2_Yerr[ii][jj]);///D2_Rad1[ii][nn]);
            nn++;
        }
        gRatio[ii]->SetMarkerStyle(22);
        gRatio[ii]->SetMarkerColor(4);
        gRatio[ii]->Draw("AP");
   }


 


}
