#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Born[5][MAXBIN],H1_Rad[5][MAXBIN],H1_Born1[5][MAXBIN],H1_Rad1[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_Born[5][MAXBIN],D2_Rad[5][MAXBIN],D2_Born1[5][MAXBIN],D2_Rad1[5][MAXBIN];
     Double_t Dp_RC[5][MAXBIN],Dp_RC1[5][MAXBIN];
     Double_t Dp_ratio[5][MAXBIN];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Born[ii][jj]=0.0; H1_Rad[ii][jj]=0.0;H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             Dp_RC[ii][jj]=0.0;Dp_RC1[ii][jj]=0.0;
	     Dp_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[5]={0,1,2,3,4};
   for(int ii=0;ii<5;ii++){
       Yfile=Form("Bodek_final/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad); 
       Yfile=Form("f1f217_new/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1); 
       Yfile=Form("Bodek_final/H1_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born,H1_Rad);
       Yfile=Form("f1f217_new/H1_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born1,H1_Rad1);
   }

   TGraph *hborn=new TGraph();
   TGraph *hrad=new TGraph();
   TGraph *hRC=new TGraph();
   TGraph *hborn1=new TGraph();
   TGraph *hrad1=new TGraph();
   TGraph *hRC1=new TGraph();
   TGraph *hratio=new TGraph();
    
   int nn=0;
   for(int ii=0;ii<5;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0)continue;
           if(abs(D2_x[ii][jj]-H1_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,D2_x[ii][jj],D2_Born[ii][jj]/H1_Born[ii][jj]);
           hrad->SetPoint(nn,D2_x[ii][jj],D2_Rad[ii][jj]/H1_Rad[ii][jj]);
           if(D2_Rad[ii][jj]>0&&H1_Rad[ii][jj]>0)Dp_RC[ii][jj]=(D2_Born[ii][jj]/D2_Rad[ii][jj])/(H1_Born[ii][jj]/H1_Rad[ii][jj]);
           hRC->SetPoint(nn,D2_x[ii][jj],Dp_RC[ii][jj]);

           hborn1->SetPoint(nn,D2_x[ii][jj],D2_Born1[ii][jj]/H1_Born1[ii][jj]);
           hrad1->SetPoint(nn,D2_x[ii][jj],D2_Rad1[ii][jj]/H1_Rad1[ii][jj]);
           if(D2_Rad1[ii][jj]>0&&H1_Rad1[ii][jj]>0)Dp_RC1[ii][jj]=(D2_Born1[ii][jj]/D2_Rad1[ii][jj])/(H1_Born1[ii][jj]/H1_Rad1[ii][jj]);
           hRC1->SetPoint(nn,D2_x[ii][jj],Dp_RC1[ii][jj]);

           if(Dp_RC1[ii][jj]>0)Dp_ratio[ii][jj]=Dp_RC[ii][jj]/Dp_RC1[ii][jj];
           hratio->SetPoint(nn,D2_x[ii][jj],Dp_ratio[ii][jj]);
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
   mg1->SetTitle("D/p born cross section ratio;xbj;born");

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
   mg2->SetTitle("D/p rad cross section ratio;xbj;rad");

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
   mg3->SetTitle("Dp RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"Bodek","P");
   leg3->AddEntry(hRC1,"f1f217","P");
   leg3->Draw();

   c2->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("Dp Bodek/f1f217 ratio;xbj;"); 


}
