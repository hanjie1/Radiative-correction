#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_H3He()
{
     Double_t He3_x[12][MAXBIN],He3_Q2[12][MAXBIN],He3_Born[12][MAXBIN],He3_Rad[12][MAXBIN],He3_Born1[12][MAXBIN],He3_Rad1[12][MAXBIN];
     Double_t H3_x[12][MAXBIN],H3_Q2[12][MAXBIN],H3_Born[12][MAXBIN],H3_Rad[12][MAXBIN],H3_Born1[12][MAXBIN],H3_Rad1[12][MAXBIN];
     Double_t He3_Born2[12][MAXBIN],He3_Rad2[12][MAXBIN],H3_Born2[12][MAXBIN],H3_Rad2[12][MAXBIN];
     Double_t H3He_RC[12][MAXBIN],H3He_RC1[12][MAXBIN],H3He_RC2[12][MAXBIN];
     Double_t H3He_ratio[12][MAXBIN],H3He_ratio1[12][MAXBIN];
     Double_t He3_Yield[12][MAXBIN],He3_Yerr[12][MAXBIN]; 
     Double_t H3_Yield[12][MAXBIN],H3_Yerr[12][MAXBIN]; 
     Double_t H3He_data[12][MAXBIN],H3He_dataerr[12][MAXBIN]; 

     for(int ii=0;ii<12;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             He3_x[ii][jj]=0.0; He3_Q2[ii][jj]=0.0; He3_Born[ii][jj]=0.0; He3_Rad[ii][jj]=0.0;He3_Born1[ii][jj]=0.0; He3_Rad1[ii][jj]=0.0;
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Born[ii][jj]=0.0; H3_Rad[ii][jj]=0.0;H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
	     He3_Born2[ii][jj]=0.0;  He3_Rad2[ii][jj]=0.0; 
	     H3_Born2[ii][jj]=0.0; H3_Rad2[ii][jj]=0.0; 
             H3He_RC[ii][jj]=0.0;H3He_RC1[ii][jj]=0.0; H3He_RC2[ii][jj]=0.0;
	     H3He_ratio[ii][jj]=0.0; H3He_ratio1[ii][jj]=0.0;
             He3_Yield[ii][jj]=0.0;  He3_Yerr[ii][jj]=0.0; 
     	     H3_Yield[ii][jj]=0.0;  H3_Yerr[ii][jj]=0.0; 
             H3He_data[ii][jj]=0.0;   H3He_dataerr[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[12]={0,1,2,3,4,5,7,9,11,13,15,16};
   for(int ii=0;ii<12;ii++){
       Yfile=Form("model121/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born,H3_Rad); 
       Yfile=Form("model122/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born1,H3_Rad1); 
       Yfile=Form("model123/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born2,H3_Rad2); 

       Yfile=Form("model121/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born,He3_Rad);
       Yfile=Form("model122/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born1,He3_Rad1);
       Yfile=Form("model123/He3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He3_x,He3_Q2,He3_Born2,He3_Rad2);

       Yfile=Form("He3_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],He3_Yield,He3_Yerr);
       Yfile=Form("H3_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],H3_Yield,H3_Yerr);
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
   for(int ii=0;ii<12;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(H3_x[ii][jj]==0)continue;
           if(abs(H3_x[ii][jj]-He3_x[ii][jj])>0.001)continue;
           hborn->SetPoint(nn,H3_x[ii][jj],H3_Born[ii][jj]/He3_Born[ii][jj]);
           hrad->SetPoint(nn,H3_x[ii][jj],H3_Rad[ii][jj]/He3_Rad[ii][jj]);
           if(H3_Rad[ii][jj]>0&&He3_Rad[ii][jj]>0)H3He_RC[ii][jj]=(H3_Born[ii][jj]/H3_Rad[ii][jj])/(He3_Born[ii][jj]/He3_Rad[ii][jj]);
           hRC->SetPoint(nn,H3_x[ii][jj],H3He_RC[ii][jj]);

           hborn1->SetPoint(nn,H3_x[ii][jj],H3_Born1[ii][jj]/He3_Born1[ii][jj]);
           hrad1->SetPoint(nn,H3_x[ii][jj],H3_Rad1[ii][jj]/He3_Rad1[ii][jj]);
           if(H3_Rad1[ii][jj]>0&&He3_Rad1[ii][jj]>0)H3He_RC1[ii][jj]=(H3_Born1[ii][jj]/H3_Rad1[ii][jj])/(He3_Born1[ii][jj]/He3_Rad1[ii][jj]);
           hRC1->SetPoint(nn,H3_x[ii][jj],H3He_RC1[ii][jj]);

           hborn2->SetPoint(nn,H3_x[ii][jj],H3_Born2[ii][jj]/He3_Born2[ii][jj]);
           hrad2->SetPoint(nn,H3_x[ii][jj],H3_Rad2[ii][jj]/He3_Rad2[ii][jj]);
           if(H3_Rad2[ii][jj]>0&&He3_Rad2[ii][jj]>0)H3He_RC2[ii][jj]=(H3_Born2[ii][jj]/H3_Rad2[ii][jj])/(He3_Born2[ii][jj]/He3_Rad2[ii][jj]);
           hRC2->SetPoint(nn,H3_x[ii][jj],H3He_RC2[ii][jj]);

           if(H3He_RC1[ii][jj]>0)H3He_ratio[ii][jj]=H3He_RC[ii][jj]/H3He_RC1[ii][jj];
           hratio->SetPoint(nn,H3_x[ii][jj],H3He_ratio[ii][jj]);
	
           if(H3He_RC2[ii][jj]>0)H3He_ratio1[ii][jj]=H3He_RC[ii][jj]/H3He_RC2[ii][jj];
           hratio1->SetPoint(nn,H3_x[ii][jj],H3He_ratio1[ii][jj]);

           nn++;
       }
   } 

//   ofstream outfile;
//   outfile.open("H3He_ratio.csv");
   TGraphErrors *gH3HeRaw[12];
   TGraphErrors *gH3He[12];
   for(int ii=0;ii<12;ii++){
       gH3HeRaw[ii]=new TGraphErrors();
       gH3He[ii]=new TGraphErrors();
       int nn=0;
       for(int jj=0;jj<MAXBIN;jj++){
          if(H3_x[ii][jj]==0||He3_x[ii][jj]==0)continue;
          H3He_data[ii][jj]=H3_Yield[ii][jj]/He3_Yield[ii][jj];
          H3He_dataerr[ii][jj]=H3He_data[ii][jj]*sqrt(pow(H3_Yerr[ii][jj]/H3_Yield[ii][jj],2)+pow(He3_Yerr[ii][jj]/He3_Yield[ii][jj],2));
          if(H3He_dataerr[ii][jj]>0.1)continue;
          gH3HeRaw[ii]->SetPoint(nn,H3_x[ii][jj],H3He_data[ii][jj]);
          gH3HeRaw[ii]->SetPointError(nn,0.0,H3He_dataerr[ii][jj]);

	  Double_t tmp_ratio=H3He_data[ii][jj]*H3He_RC1[ii][jj];
	  Double_t tmp_err=H3He_dataerr[ii][jj]*H3He_RC1[ii][jj];
 	  gH3He[ii]->SetPoint(nn,H3_x[ii][jj],tmp_ratio);
 	  gH3He[ii]->SetPointError(nn,0.0,tmp_err);

//	  outfile<<H3_x[ii][jj]<<","<<H3_Q2[ii][jj]<<","<<tmp_ratio<<","<<tmp_err<<endl;
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
   mg1->SetTitle("H3/He3 born cross section ratio;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hborn,"model121","P");
   leg1->AddEntry(hborn1,"model122","P");
   leg1->AddEntry(hborn2,"model123","P");
   leg1->Draw();

   TCanvas *c3=new TCanvas("c3","c3",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(4);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(8);
   hrad2->SetMarkerStyle(8);
   hrad2->SetMarkerColor(1);
   int color[12]={1,2,3,4,5,6,7,8,9,46,30,38};
   for(int ii=0;ii<12;ii++){
       gH3HeRaw[ii]->SetMarkerStyle(22);
       gH3HeRaw[ii]->SetMarkerColor(color[ii]);
       gH3He[ii]->SetMarkerStyle(8);
       gH3He[ii]->SetMarkerColor(color[ii]);
       mg2->Add(gH3HeRaw[ii]);
       mg2->Add(gH3He[ii]);
   }
//   mg2->Add(hrad);
//   mg2->Add(hrad1);
   mg2->Draw("AP");
   mg2->SetTitle("H3/He3 rad cross section ratio;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->AddEntry(hrad,"model121","P");
   leg2->AddEntry(hrad1,"model122","P");
   leg2->AddEntry(hrad2,"model123","P");
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
   mg3->SetTitle("H3/He3 RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"model121","P");
   leg3->AddEntry(hRC1,"model122","P");
   leg3->AddEntry(hRC2,"model123","P");
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
   mg4->SetTitle("H3/He3 RC ratio between models;xbj;");

   auto leg4=new TLegend(0.7,0.6,0.8,0.8);
   leg4->AddEntry(hratio,"model121/model122","P");
   leg4->AddEntry(hratio1,"model121/model123","P");
   leg4->Draw();



}
