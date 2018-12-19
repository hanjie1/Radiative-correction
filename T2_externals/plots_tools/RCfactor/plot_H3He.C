#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Born[11][MAXBIN],He3_Rad[11][MAXBIN],He3_Born1[11][MAXBIN],He3_Rad1[11][MAXBIN];
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Born[11][MAXBIN],H3_Rad[11][MAXBIN],H3_Born1[11][MAXBIN],H3_Rad1[11][MAXBIN];
     Double_t Dp_RC[11][MAXBIN],Dp_RC1[11][MAXBIN];
     Double_t Dp_ratio[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Born[ii][jj]=0.0; He3_Rad[ii][jj]=0.0;He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Born[ii][jj]=0.0; H3_Rad[ii][jj]=0.0;H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
             Dp_RC[ii][jj]=0.0;Dp_RC1[ii][jj]=0.0;
	     Dp_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("Bodek_final/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born,H3_Rad); 
       Yfile=Form("f1f217_new/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born1,H3_Rad1); 
       Yfile=Form("Bodek_final/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born,He3_Rad);
       Yfile=Form("f1f217_new/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born1,He3_Rad1);
   }

   TGraph *hborn=new TGraph();
   TGraph *hrad=new TGraph();
   TGraph *hRC=new TGraph();
   TGraph *hborn1=new TGraph();
   TGraph *hrad1=new TGraph();
   TGraph *hRC1=new TGraph();
   TGraph *hratio=new TGraph();
    
   int nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_x[ii][jj]-He3_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,H3_x[ii][jj],H3_Born[ii][jj]/He3_Born[ii][jj]);
           hrad->SetPoint(nn,H3_x[ii][jj],H3_Rad[ii][jj]/He3_Rad[ii][jj]);
           if(H3_Rad[ii][jj]>0&&He3_Rad[ii][jj]>0)Dp_RC[ii][jj]=(H3_Born[ii][jj]/H3_Rad[ii][jj])/(He3_Born[ii][jj]/He3_Rad[ii][jj]);
           hRC->SetPoint(nn,H3_x[ii][jj],Dp_RC[ii][jj]);

           hborn1->SetPoint(nn,H3_x[ii][jj],H3_Born1[ii][jj]/He3_Born1[ii][jj]);
           hrad1->SetPoint(nn,H3_x[ii][jj],H3_Rad1[ii][jj]/He3_Rad1[ii][jj]);
           if(H3_Rad1[ii][jj]>0&&He3_Rad1[ii][jj]>0)Dp_RC1[ii][jj]=(H3_Born1[ii][jj]/H3_Rad1[ii][jj])/(He3_Born1[ii][jj]/He3_Rad1[ii][jj]);
           hRC1->SetPoint(nn,H3_x[ii][jj],Dp_RC1[ii][jj]);

           if(Dp_RC1[ii][jj]>0)Dp_ratio[ii][jj]=Dp_RC[ii][jj]/Dp_RC1[ii][jj];
           hratio->SetPoint(nn,H3_x[ii][jj],Dp_ratio[ii][jj]);
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
   mg1->SetTitle("H3/He3 born cross section ratio;xbj;born");

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
   mg2->SetTitle("H3/He3 rad cross section ratio;xbj;rad");

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
   mg3->SetTitle("H3/He3 RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"Bodek","P");
   leg3->AddEntry(hRC1,"f1f217","P");
   leg3->Draw();

   c2->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("H3/He3 Bodek/f1f217 ratio;xbj;"); 


}
