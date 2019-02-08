#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Rad[11][MAXBIN],D2_Born1[11][MAXBIN],D2_Rad1[11][MAXBIN];
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Born[11][MAXBIN],He3_Rad[11][MAXBIN],He3_Born1[11][MAXBIN],He3_Rad1[11][MAXBIN];
     Double_t HeD_RC[11][MAXBIN],HeD_RC1[11][MAXBIN];
     Double_t HeD_ratio[11][MAXBIN];
     Double_t D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN]; 
     Double_t He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN]; 
     Double_t HeD_data[11][MAXBIN],HeD_dataerr[11][MAXBIN]; 

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Born[ii][jj]=0.0; He3_Rad[ii][jj]=0.0;He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             HeD_RC[ii][jj]=0.0;HeD_RC1[ii][jj]=0.0;
	     HeD_ratio[ii][jj]=0.0;
             D2_Yield[ii][jj]=0.0;  D2_Yerr[ii][jj]=0.0; 
     	     He3_Yield[ii][jj]=0.0;  He3_Yerr[ii][jj]=0.0; 
             HeD_data[ii][jj]=0.0;   HeD_dataerr[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("model211/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born,He3_Rad); 
       Yfile=Form("model111/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born1,He3_Rad1); 
       Yfile=Form("model211/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad);
       Yfile=Form("model111/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1);

       Yfile=Form("D2_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],D2_Yield,D2_Yerr);
       Yfile=Form("He3_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],He3_Yield,He3_Yerr);
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
	   if(He3_x[ii][jj]==0)continue;
           if(abs(He3_x[ii][jj]-D2_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,He3_x[ii][jj],He3_Born[ii][jj]/D2_Born[ii][jj]);
           hrad->SetPoint(nn,He3_x[ii][jj],He3_Rad[ii][jj]/D2_Rad[ii][jj]);
           if(He3_Rad[ii][jj]>0&&D2_Rad[ii][jj]>0)HeD_RC[ii][jj]=(He3_Born[ii][jj]/He3_Rad[ii][jj])/(D2_Born[ii][jj]/D2_Rad[ii][jj]);
           hRC->SetPoint(nn,He3_x[ii][jj],HeD_RC[ii][jj]);

           hborn1->SetPoint(nn,He3_x[ii][jj],He3_Born1[ii][jj]/D2_Born1[ii][jj]);
           hrad1->SetPoint(nn,He3_x[ii][jj],He3_Rad1[ii][jj]/D2_Rad1[ii][jj]);
           if(He3_Rad1[ii][jj]>0&&D2_Rad1[ii][jj]>0)HeD_RC1[ii][jj]=(He3_Born1[ii][jj]/He3_Rad1[ii][jj])/(D2_Born1[ii][jj]/D2_Rad1[ii][jj]);
           hRC1->SetPoint(nn,He3_x[ii][jj],HeD_RC1[ii][jj]);

           if(HeD_RC1[ii][jj]>0)HeD_ratio[ii][jj]=HeD_RC[ii][jj]/HeD_RC1[ii][jj];
           hratio->SetPoint(nn,He3_x[ii][jj],HeD_ratio[ii][jj]);
	
           nn++;
       }
   } 

//   ofstream outfile;
//   outfile.open("HeD_ratio.csv");
   TGraphErrors *gHeDRaw[11];
   TGraphErrors *gHeD[11];
   for(int ii=0;ii<11;ii++){
       gHeDRaw[ii]=new TGraphErrors();
       gHeD[ii]=new TGraphErrors();
       int nn=0;
       for(int jj=0;jj<MAXBIN;jj++){
          if(He3_x[ii][jj]==0||D2_x[ii][jj]==0)continue;
          HeD_data[ii][jj]=He3_Yield[ii][jj]/D2_Yield[ii][jj];
          HeD_dataerr[ii][jj]=HeD_data[ii][jj]*sqrt(pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2)+pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2));
          if(HeD_dataerr[ii][jj]>0.1)continue;
          gHeDRaw[ii]->SetPoint(nn,He3_x[ii][jj],HeD_data[ii][jj]);
          gHeDRaw[ii]->SetPointError(nn,0.0,HeD_dataerr[ii][jj]);

	  Double_t tmp_ratio=HeD_data[ii][jj]*HeD_RC1[ii][jj];
	  Double_t tmp_err=HeD_dataerr[ii][jj]*HeD_RC1[ii][jj];
 	  gHeD[ii]->SetPoint(nn,He3_x[ii][jj],tmp_ratio);
 	  gHeD[ii]->SetPointError(nn,0.0,tmp_err);

//	  outfile<<He3_x[ii][jj]<<","<<He3_Q2[ii][jj]<<","<<tmp_ratio<<","<<tmp_err<<endl;
          nn++;
      }
   }
//   outfile.close();

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hborn->SetMarkerStyle(8);
   hborn->SetMarkerColor(4);
   hborn1->SetMarkerStyle(8);
   hborn1->SetMarkerColor(9);
   mg1->Add(hborn);
   mg1->Add(hborn1);
   mg1->Draw("AP");
   mg1->SetTitle("D/p born cross section ratio;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hborn,"Bodek","P");
   leg1->AddEntry(hborn1,"model111","P");
   leg1->Draw();

   TCanvas *c3=new TCanvas("c3","c3",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(4);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(8);
   int color[11]={1,2,3,4,5,6,7,8,9,46,30};
   for(int ii=0;ii<11;ii++){
       gHeDRaw[ii]->SetMarkerStyle(22);
       gHeDRaw[ii]->SetMarkerColor(color[ii]);
       gHeD[ii]->SetMarkerStyle(8);
       gHeD[ii]->SetMarkerColor(color[ii]);
       mg2->Add(gHeDRaw[ii]);
       mg2->Add(gHeD[ii]);
   }
//   mg2->Add(hrad);
//   mg2->Add(hrad1);
   mg2->Draw("AP");
   mg2->SetTitle("D/p rad cross section ratio;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->AddEntry(hrad,"Bodek","P");
   leg2->AddEntry(hrad1,"model111","P");
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
   mg3->SetTitle("HeD RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"Bodek","P");
   leg3->AddEntry(hRC1,"model111","P");
   leg3->Draw();

   c2->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("HeD Bodek/model111 ratio;xbj;"); 


}
