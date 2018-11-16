#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_C12D()
{
     Double_t D2_x[1][MAXBIN],D2_Q2[1][MAXBIN],D2_Born[1][MAXBIN],D2_Rad[1][MAXBIN],D2_Born1[1][MAXBIN],D2_Rad1[1][MAXBIN];
     Double_t Carbon_x[1][MAXBIN],Carbon_Q2[1][MAXBIN],Carbon_Born[1][MAXBIN],Carbon_Rad[1][MAXBIN],Carbon_Born1[1][MAXBIN],Carbon_Rad1[1][MAXBIN];
     Double_t Dp_RC[1][MAXBIN],Dp_RC1[1][MAXBIN];
     Double_t Dp_ratio[1][MAXBIN];

     for(int ii=0;ii<1;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             Carbon_x[ii][jj]=0.0; Carbon_Q2[ii][jj]=0.0; Carbon_Born[ii][jj]=0.0; Carbon_Rad[ii][jj]=0.0;Carbon_Born1[ii][jj]=0.0; Carbon_Rad1[ii][jj]=0.0;
             Dp_RC[ii][jj]=0.0;Dp_RC1[ii][jj]=0.0;
	     Dp_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[1]={0};
   for(int ii=0;ii<1;ii++){
       Yfile="Bodek_final/Carbon_all_xs.out";
       ReadYield(Yfile,kin[ii],Carbon_x,Carbon_Q2,Carbon_Born,Carbon_Rad); 
       Yfile="f1f217/Carbon_all_xs.out";
       ReadYield(Yfile,kin[ii],Carbon_x,Carbon_Q2,Carbon_Born1,Carbon_Rad1); 
       Yfile="Bodek_final/D2_all_xs.out";
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad);
       Yfile="f1f217/D2_all_xs.out";
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1);
   }

   TGraph *hborn=new TGraph();
   TGraph *hrad=new TGraph();
   TGraph *hRC=new TGraph();
   TGraph *hborn1=new TGraph();
   TGraph *hrad1=new TGraph();
   TGraph *hRC1=new TGraph();
   TGraph *hratio=new TGraph();
    
   int nn=0;
   for(int ii=0;ii<1;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(Carbon_x[ii][jj]==0)continue;
           if(abs(Carbon_x[ii][jj]-D2_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,Carbon_x[ii][jj],Carbon_Born[ii][jj]/D2_Born[ii][jj]/6.0);
           hrad->SetPoint(nn,Carbon_x[ii][jj],Carbon_Rad[ii][jj]/D2_Rad[ii][jj]/6.0);
           if(Carbon_Rad[ii][jj]>0&&D2_Rad[ii][jj]>0)Dp_RC[ii][jj]=(Carbon_Born[ii][jj]/Carbon_Rad[ii][jj])/(D2_Born[ii][jj]/D2_Rad[ii][jj]);
           hRC->SetPoint(nn,Carbon_x[ii][jj],Dp_RC[ii][jj]);

           hborn1->SetPoint(nn,Carbon_x[ii][jj],Carbon_Born1[ii][jj]/D2_Born1[ii][jj]/6.0);
           hrad1->SetPoint(nn,Carbon_x[ii][jj],Carbon_Rad1[ii][jj]/D2_Rad1[ii][jj]/6.0);
           if(Carbon_Rad1[ii][jj]>0&&D2_Rad1[ii][jj]>0)Dp_RC1[ii][jj]=(Carbon_Born1[ii][jj]/Carbon_Rad1[ii][jj])/(D2_Born1[ii][jj]/D2_Rad1[ii][jj]);
           hRC1->SetPoint(nn,Carbon_x[ii][jj],Dp_RC1[ii][jj]);

           if(Dp_RC1[ii][jj]>0)Dp_ratio[ii][jj]=Dp_RC[ii][jj]/Dp_RC1[ii][jj];
           hratio->SetPoint(nn,Carbon_x[ii][jj],Dp_ratio[ii][jj]);
           nn++;
       }
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   c1->Divide(2,1);
   c1->cd(1);
   TMultiGraph *mg1=new TMultiGraph();
   hborn->SetMarkerStyle(8);
   hborn->SetMarkerColor(4);
   hborn1->SetMarkerStyle(8);
   hborn1->SetMarkerColor(2);
   mg1->Add(hborn);
   mg1->Add(hborn1);
   mg1->Draw("AP");
   mg1->SetTitle("Carbon/D born cross section ratio;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hborn,"Bodek","P");
   leg1->AddEntry(hborn1,"f1f217","P");
   leg1->Draw();

   c1->cd(2);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(4);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(2);
   mg2->Add(hrad);
   mg2->Add(hrad1);
   mg2->Draw("AP");
   mg2->SetTitle("Carbon/D rad cross section ratio;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->AddEntry(hrad,"Bodek","P");
   leg2->AddEntry(hrad1,"f1f217","P");
   leg2->Draw();

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   c2->Divide(2,1);
   c2->cd(1);
   TMultiGraph *mg3=new TMultiGraph();
   hRC->SetMarkerStyle(8);
   hRC->SetMarkerColor(4);
   hRC1->SetMarkerStyle(8);
   hRC1->SetMarkerColor(2);
   mg3->Add(hRC);
   mg3->Add(hRC1);
   mg3->Draw("AP");
   mg3->SetTitle("Carbon/D RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"Bodek","P");
   leg3->AddEntry(hRC1,"f1f217","P");
   leg3->Draw();

   c2->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("Carbon/D Bodek/f1f217 ratio;xbj;"); 


}
