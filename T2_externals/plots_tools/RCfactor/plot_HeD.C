#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_HeD()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Rad[11][MAXBIN],D2_Born1[11][MAXBIN],D2_Rad1[11][MAXBIN];
     Double_t He3_x[11][MAXBIN],He3_Q2[11][MAXBIN],He3_Born[11][MAXBIN],He3_Rad[11][MAXBIN],He3_Born1[11][MAXBIN],He3_Rad1[11][MAXBIN];
     Double_t D2_Born2[11][MAXBIN],D2_Rad2[11][MAXBIN],He3_Born2[11][MAXBIN],He3_Rad2[11][MAXBIN];
     Double_t HeD_RC[11][MAXBIN],HeD_RC1[11][MAXBIN],HeD_RC2[11][MAXBIN];
     Double_t HeD_ratio[11][MAXBIN],HeD_ratio1[11][MAXBIN];
     Double_t D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN]; 
     Double_t He3_Yield[11][MAXBIN],He3_Yerr[11][MAXBIN]; 
     Double_t HeD_data[11][MAXBIN],HeD_dataerr[11][MAXBIN]; 

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Born[ii][jj]=0.0; He3_Rad[ii][jj]=0.0;He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
	     D2_Born2[ii][jj]=0.0;  D2_Rad2[ii][jj]=0.0; 
	     He3_Born2[ii][jj]=0.0; He3_Rad2[ii][jj]=0.0; 
             HeD_RC[ii][jj]=0.0;HeD_RC1[ii][jj]=0.0; HeD_RC2[ii][jj]=0.0;
	     HeD_ratio[ii][jj]=0.0; HeD_ratio1[ii][jj]=0.0;
             D2_Yield[ii][jj]=0.0;  D2_Yerr[ii][jj]=0.0; 
     	     He3_Yield[ii][jj]=0.0;  He3_Yerr[ii][jj]=0.0; 
             HeD_data[ii][jj]=0.0;   HeD_dataerr[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("model111/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born,He3_Rad); 
       Yfile=Form("model311/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born1,He3_Rad1); 
       Yfile=Form("model211/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born2,He3_Rad2); 

       Yfile=Form("model111/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad);
       Yfile=Form("model311/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1);
       Yfile=Form("model211/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born2,D2_Rad2);

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
   TGraph *hborn2=new TGraph();
   TGraph *hrad2=new TGraph();
   TGraph *hRC2=new TGraph();
   TGraph *hratio=new TGraph();
   TGraph *hratio1=new TGraph();
    
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

           hborn2->SetPoint(nn,He3_x[ii][jj],He3_Born2[ii][jj]/D2_Born2[ii][jj]);
           hrad2->SetPoint(nn,He3_x[ii][jj],He3_Rad2[ii][jj]/D2_Rad2[ii][jj]);
           if(He3_Rad2[ii][jj]>0&&D2_Rad2[ii][jj]>0)HeD_RC2[ii][jj]=(He3_Born2[ii][jj]/He3_Rad2[ii][jj])/(D2_Born2[ii][jj]/D2_Rad2[ii][jj]);
           hRC2->SetPoint(nn,He3_x[ii][jj],HeD_RC2[ii][jj]);

           if(HeD_RC1[ii][jj]>0)HeD_ratio[ii][jj]=HeD_RC[ii][jj]/HeD_RC1[ii][jj];
           hratio->SetPoint(nn,He3_x[ii][jj],HeD_ratio[ii][jj]);
	
           if(HeD_RC2[ii][jj]>0)HeD_ratio1[ii][jj]=HeD_RC[ii][jj]/HeD_RC2[ii][jj];
           hratio1->SetPoint(nn,He3_x[ii][jj],HeD_ratio1[ii][jj]);

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
   hborn->SetMarkerColor(2);
   hborn1->SetMarkerStyle(8);
   hborn1->SetMarkerColor(4);
   hborn2->SetMarkerStyle(8);
   hborn2->SetMarkerColor(1);
   mg1->Add(hborn);
   mg1->Add(hborn1);
   mg1->Add(hborn2);
   mg1->Draw("AP");
   mg1->SetTitle("He/D born cross section ratio;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hborn,"model111","P");
   leg1->AddEntry(hborn1,"model311","P");
   leg1->AddEntry(hborn2,"model211","P");
   leg1->Draw();

   TCanvas *c3=new TCanvas("c3","c3",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(4);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(8);
   hrad2->SetMarkerStyle(8);
   hrad2->SetMarkerColor(1);
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
   mg2->SetTitle("He/D rad cross section ratio;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->AddEntry(hrad,"model211","P");
   leg2->AddEntry(hrad1,"model211_noResAll","P");
   leg2->AddEntry(hrad2,"model211_ResOnlyD2H1","P");
   leg2->Draw();

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   c2->Divide(2,1);
   c2->cd(1);
   TMultiGraph *mg3=new TMultiGraph();
   hRC->SetMarkerStyle(8);
   hRC->SetMarkerColor(2);
   hRC1->SetMarkerStyle(8);
   hRC1->SetMarkerColor(4);
   hRC2->SetMarkerStyle(8);
   hRC2->SetMarkerColor(1);
   mg3->Add(hRC);
   mg3->Add(hRC1);
   mg3->Add(hRC2);
   mg3->Draw("AP");
   mg3->SetTitle("HeD RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"model111","P");
   leg3->AddEntry(hRC1,"model311","P");
   leg3->AddEntry(hRC2,"model211","P");
   leg3->Draw();

   c2->cd(2);
   TMultiGraph *mg4=new TMultiGraph();
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio1->SetMarkerStyle(8);
   hratio1->SetMarkerColor(2);
   mg4->Add(hratio);
   mg4->Add(hratio1);
   mg4->Draw("AP");
   mg4->SetTitle("He3/D RC ratio between models;xbj;");

   auto leg4=new TLegend(0.7,0.6,0.811,0.811);
   leg4->AddEntry(hratio,"model111/model311","P");
   leg4->AddEntry(hratio1,"model111/model211","P");
   leg4->Draw();



}
