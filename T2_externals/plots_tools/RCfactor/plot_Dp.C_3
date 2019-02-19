#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_Dp()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Born[5][MAXBIN],H1_Rad[5][MAXBIN],H1_Born1[5][MAXBIN],H1_Rad1[5][MAXBIN];
     Double_t D2_x[5][MAXBIN],D2_Q2[5][MAXBIN],D2_Born[5][MAXBIN],D2_Rad[5][MAXBIN],D2_Born1[5][MAXBIN],D2_Rad1[5][MAXBIN];
     Double_t H1_Born2[5][MAXBIN],H1_Rad2[5][MAXBIN],D2_Born2[5][MAXBIN],D2_Rad2[5][MAXBIN];
     Double_t Dp_RC[5][MAXBIN],Dp_RC1[5][MAXBIN],Dp_RC2[5][MAXBIN];
     Double_t Dp_ratio[5][MAXBIN],Dp_ratio1[5][MAXBIN];
     Double_t H1_Yield[5][MAXBIN],H1_Yerr[5][MAXBIN]; 
     Double_t D2_Yield[5][MAXBIN],D2_Yerr[5][MAXBIN]; 
     Double_t Dp_data[5][MAXBIN],Dp_dataerr[5][MAXBIN]; 

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Born[ii][jj]=0.0; H1_Rad[ii][jj]=0.0;H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
	     H1_Born2[ii][jj]=0.0;  H1_Rad2[ii][jj]=0.0; 
	     D2_Born2[ii][jj]=0.0; D2_Rad2[ii][jj]=0.0; 
             Dp_RC[ii][jj]=0.0;Dp_RC1[ii][jj]=0.0; Dp_RC2[ii][jj]=0.0;
	     Dp_ratio[ii][jj]=0.0; Dp_ratio1[ii][jj]=0.0;
             H1_Yield[ii][jj]=0.0;  H1_Yerr[ii][jj]=0.0; 
     	     D2_Yield[ii][jj]=0.0;  D2_Yerr[ii][jj]=0.0; 
             Dp_data[ii][jj]=0.0;   Dp_dataerr[ii][jj]=0.0; 
     }}

   TString Yfile;
   int kin[5]={0,1,2,3,4};
   for(int ii=0;ii<5;ii++){
       Yfile=Form("model111/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad); 
       Yfile=Form("model211/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1); 
       Yfile=Form("model311/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born2,D2_Rad2); 

       Yfile=Form("model111/H1_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born,H1_Rad);
       Yfile=Form("model211/H1_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born1,H1_Rad1);
       Yfile=Form("model311/H1_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born2,H1_Rad2);

       Yfile=Form("H1_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],H1_Yield,H1_Yerr);
       Yfile=Form("D2_kin%d.txt",kin[ii]);
       ReadData(Yfile,kin[ii],D2_Yield,D2_Yerr);
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

           hborn2->SetPoint(nn,D2_x[ii][jj],D2_Born2[ii][jj]/H1_Born2[ii][jj]);
           hrad2->SetPoint(nn,D2_x[ii][jj],D2_Rad2[ii][jj]/H1_Rad2[ii][jj]);
           if(D2_Rad2[ii][jj]>0&&H1_Rad2[ii][jj]>0)Dp_RC2[ii][jj]=(D2_Born2[ii][jj]/D2_Rad2[ii][jj])/(H1_Born2[ii][jj]/H1_Rad2[ii][jj]);
           hRC2->SetPoint(nn,D2_x[ii][jj],Dp_RC2[ii][jj]);

           if(Dp_RC1[ii][jj]>0)Dp_ratio[ii][jj]=Dp_RC[ii][jj]/Dp_RC1[ii][jj];
           hratio->SetPoint(nn,D2_x[ii][jj],Dp_ratio[ii][jj]);
	
           if(Dp_RC2[ii][jj]>0)Dp_ratio1[ii][jj]=Dp_RC[ii][jj]/Dp_RC2[ii][jj];
           hratio1->SetPoint(nn,D2_x[ii][jj],Dp_ratio1[ii][jj]);

           nn++;
       }
   } 

//   ofstream outfile;
//   outfile.open("Dp_ratio.csv");
   TGraphErrors *gDpRaw[5];
   TGraphErrors *gDp[5];
   for(int ii=0;ii<5;ii++){
       gDpRaw[ii]=new TGraphErrors();
       gDp[ii]=new TGraphErrors();
       int nn=0;
       for(int jj=0;jj<MAXBIN;jj++){
          if(D2_x[ii][jj]==0||H1_x[ii][jj]==0)continue;
          Dp_data[ii][jj]=D2_Yield[ii][jj]/H1_Yield[ii][jj];
          Dp_dataerr[ii][jj]=Dp_data[ii][jj]*sqrt(pow(D2_Yerr[ii][jj]/D2_Yield[ii][jj],2)+pow(H1_Yerr[ii][jj]/H1_Yield[ii][jj],2));
          if(Dp_dataerr[ii][jj]>0.1)continue;
          gDpRaw[ii]->SetPoint(nn,D2_x[ii][jj],Dp_data[ii][jj]);
          gDpRaw[ii]->SetPointError(nn,0.0,Dp_dataerr[ii][jj]);

	  Double_t tmp_ratio=Dp_data[ii][jj]*Dp_RC1[ii][jj];
	  Double_t tmp_err=Dp_dataerr[ii][jj]*Dp_RC1[ii][jj];
 	  gDp[ii]->SetPoint(nn,D2_x[ii][jj],tmp_ratio);
 	  gDp[ii]->SetPointError(nn,0.0,tmp_err);

//	  outfile<<D2_x[ii][jj]<<","<<D2_Q2[ii][jj]<<","<<tmp_ratio<<","<<tmp_err<<endl;
          nn++;
      }
   }
//   outfile.close();

   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   TMultiGraph *mg1=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(2);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(4);
   hrad2->SetMarkerStyle(8);
   hrad2->SetMarkerColor(1);
   mg1->Add(hrad);
   mg1->Add(hrad1);
   mg1->Add(hrad2);
   mg1->Draw("AP");
   mg1->SetTitle("He/D rad cross section ratio;xbj;rad");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(hrad,"model111","P");
   leg1->AddEntry(hrad1,"model211","P");
   leg1->AddEntry(hrad2,"model311","P");
   leg1->Draw();

   TCanvas *c3=new TCanvas("c3","c3",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(4);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(8);
   hrad2->SetMarkerStyle(8);
   hrad2->SetMarkerColor(1);
   int color[5]={1,2,3,4,5};
   for(int ii=0;ii<5;ii++){
       gDpRaw[ii]->SetMarkerStyle(22);
       gDpRaw[ii]->SetMarkerColor(color[ii]);
       gDp[ii]->SetMarkerStyle(8);
       gDp[ii]->SetMarkerColor(color[ii]);
       mg2->Add(gDpRaw[ii]);
       mg2->Add(gDp[ii]);
   }
//   mg2->Add(hrad);
//   mg2->Add(hrad1);
   mg2->Draw("AP");
   mg2->SetTitle("He/D rad cross section ratio;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.85,0.85);
   leg2->AddEntry(hrad,"model111","P");
   leg2->AddEntry(hrad1,"model211","P");
   leg2->AddEntry(hrad2,"model311","P");
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
   mg3->SetTitle("Dp RC=born/rad ratio;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.85,0.85);
   leg3->AddEntry(hRC,"model111","P");
   leg3->AddEntry(hRC1,"model211","P");
   leg3->AddEntry(hRC2,"model311","P");
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
   mg4->SetTitle("D/p RC ratio between models;xbj;");

   auto leg4=new TLegend(0.7,0.6,0.85,0.85);
   leg4->AddEntry(hratio,"model111/model211","P");
   leg4->AddEntry(hratio1,"model111/model311","P");
   leg4->Draw();



}
