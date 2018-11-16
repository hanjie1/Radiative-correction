#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Born1[11][MAXBIN];
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Born[11][MAXBIN],He3_Born1[11][MAXBIN];
     Double_t Dp_ratio[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Born1[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Born[ii][jj]=0.0; He3_Born1[ii][jj]=0.0;
	     Dp_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("INEFT/He3_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born); 
       Yfile=Form("f1f217/He3_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born1); 
       Yfile=Form("INEFT/D2_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born);
       Yfile=Form("f1f217/D2_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1);
   }

   TGraph *hborn=new TGraph();
   TGraph *hborn1=new TGraph();
   TGraph *hratio=new TGraph();
    
   int nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(He3_x[ii][jj]==0)continue;
           if(abs(He3_x[ii][jj]-D2_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,He3_x[ii][jj],He3_Born[ii][jj]/D2_Born[ii][jj]);
           hborn1->SetPoint(nn,He3_x[ii][jj],He3_Born1[ii][jj]/D2_Born1[ii][jj]);

           nn++;
       }
   } 

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hborn->SetMarkerStyle(8);
   hborn->SetMarkerColor(4);
   hborn1->SetMarkerStyle(8);
   hborn1->SetMarkerColor(2);
   mg1->Add(hborn);
   mg1->Add(hborn1);
   mg1->Draw("AP");
   mg1->SetTitle("He3/D born cross section ratio;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hborn,"Bodek","P");
   leg1->AddEntry(hborn1,"f1f217","P");
   leg1->Draw();


}
