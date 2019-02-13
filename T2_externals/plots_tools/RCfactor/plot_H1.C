#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H1()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Born[5][MAXBIN],H1_Rad[5][MAXBIN],H1_Born1[5][MAXBIN],H1_Rad1[5][MAXBIN];
     Double_t H1_RC[5][MAXBIN],H1_RC1[5][MAXBIN];
     Double_t H1_ratio[5][MAXBIN];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Born[ii][jj]=0.0; H1_Rad[ii][jj]=0.0;H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
             H1_RC[ii][jj]=0.0;H1_RC1[ii][jj]=0.0;
	     H1_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   for(int ii=0;ii<5;ii++){
       Yfile=Form("Bodek/H1_kin%d_xs.out",ii);
       ReadYield(Yfile,ii,H1_x,H1_Q2,H1_Born,H1_Rad); 
       Yfile=Form("model211_1/H1_kin%d_xs.out",ii);
       ReadYield(Yfile,ii,H1_x,H1_Q2,H1_Born1,H1_Rad1); 
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
	   if(H1_x[ii][jj]==0)continue;
           hborn->SetPoint(nn,H1_x[ii][jj],H1_Born[ii][jj]);
           hrad->SetPoint(nn,H1_x[ii][jj],H1_Rad[ii][jj]);
           if(H1_Rad[ii][jj]>0)H1_RC[ii][jj]=H1_Born[ii][jj]/H1_Rad[ii][jj];
           //hRC->SetPoint(nn,H1_x[ii][jj],H1_RC[ii][jj]);
           hRC->SetPoint(nn,H1_x[ii][jj],H1_Born[ii][jj]/H1_Born1[ii][jj]);

           hborn1->SetPoint(nn,H1_x[ii][jj],H1_Born1[ii][jj]);
           hrad1->SetPoint(nn,H1_x[ii][jj],H1_Rad1[ii][jj]);
           if(H1_Rad1[ii][jj]>0)H1_RC1[ii][jj]=H1_Born1[ii][jj]/H1_Rad1[ii][jj];
          // hRC1->SetPoint(nn,H1_x[ii][jj],H1_RC1[ii][jj]);
           hRC1->SetPoint(nn,H1_x[ii][jj],H1_Rad[ii][jj]/H1_Rad1[ii][jj]);

           if(H1_RC1[ii][jj]>0)H1_ratio[ii][jj]=H1_RC[ii][jj]/H1_RC1[ii][jj];
           hratio->SetPoint(nn,H1_x[ii][jj],H1_ratio[ii][jj]);
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
   mg1->SetTitle("H1 born cross section;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.811,0.811);
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
   mg2->SetTitle("H1 rad cross section;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.811,0.811);
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
   mg3->SetTitle("H1 RC=born/rad;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.811,0.811);
   leg3->AddEntry(hRC,"Bodek","P");
   leg3->AddEntry(hRC1,"f1f217","P");
   leg3->Draw();

   c2->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("H1 Bodek/f1f217;xbj;"); 


}
