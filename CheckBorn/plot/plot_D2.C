#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_D2()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Born1[11][MAXBIN];
     Int_t MAXBIN1=11*MAXBIN;
     Double_t D2_Dx[MAXBIN1],D2_DQ2[MAXBIN1],D2_Dborn[MAXBIN1],D2_Dborn1[MAXBIN1],D2_Derr[MAXBIN1],D2_Derr1[MAXBIN1];
     int mark[MAXBIN1];
     Double_t D2_ratio[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0;D2_Born1[ii][jj]=0.0;
	     D2_ratio[ii][jj]=0.0;
     }}

   for(int ii=0;ii<MAXBIN1;ii++){
	D2_Dx[ii]=0.0; D2_DQ2[ii]=0.0; D2_Dborn[ii]=0.0; D2_Dborn1[ii]=0.0;
	D2_Derr[ii]=0.0; D2_Derr1[ii]=0.0;
        mark[ii]=0;
   }

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("INEFT/D2_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born); 
       Yfile=Form("f1f217_new/D2_kin%d.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1); 
   }

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    ifstream Datafile;
    TString filename="../OUT/DataBorn/D2.out";
    Datafile.open(filename);
    while(tmp.ReadLine(Datafile)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          D2_Dx[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          D2_DQ2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  D2_Dborn[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  D2_Derr[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  D2_Dborn1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  D2_Derr1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  mark[nn-1]=atoi(content.Data());
          
          from=0;
          nn++;
     }
    Datafile.close();


   TGraph *hborn=new TGraph();
   TGraph *hborn1=new TGraph();
   TGraph *hratio=new TGraph();
    
   nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0)continue;
           hborn->SetPoint(nn,D2_x[ii][jj],D2_Born[ii][jj]);
           hborn1->SetPoint(nn,D2_x[ii][jj],D2_Born1[ii][jj]);

           if(D2_Born1[ii][jj]>0)D2_ratio[ii][jj]=D2_Born[ii][jj]/D2_Born1[ii][jj];
           hratio->SetPoint(nn,D2_x[ii][jj],D2_ratio[ii][jj]);
           nn++;
       }
   } 
   
   TGraphErrors *hdata=new TGraphErrors();
   TGraphErrors *hdata1=new TGraphErrors();
   
   nn=0;
   Double_t tmpNorm=1.0;
   if(abs(D2_Dx[0]-D2_x[0][0])<0.001)tmpNorm=D2_Born[0][0]/D2_Dborn[0];
   else cout<<"The sequence of data and RC file are not the same !"<<endl;
   cout<<D2_Dx[0]<<"  "<<D2_x[0][0]<<endl;
   for(int ii=0;ii<MAXBIN1;ii++){
	if(D2_Dx[ii]==0||mark[ii]==0)continue;
	hdata->SetPoint(nn,D2_Dx[ii],D2_Dborn[ii]*tmpNorm);
        hdata->SetPointError(nn,0.0,D2_Derr[ii]*tmpNorm);
	hdata1->SetPoint(nn,D2_Dx[ii],D2_Dborn1[ii]*tmpNorm);
        hdata1->SetPointError(nn,0.0,D2_Derr1[ii]*tmpNorm);
        nn++;
   }

   TGraph *hDRratio=new TGraph();
   TGraph *hDRratio1=new TGraph();

   nn=0;int mm=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0.0)continue;
           Double_t tmp1=D2_Dborn[nn]*tmpNorm/D2_Born[ii][jj];
           Double_t tmp2=D2_Dborn[nn]*tmpNorm/D2_Born1[ii][jj];
           if(mark[nn]==1){
	      hDRratio->SetPoint(mm,D2_Dx[nn],tmp1);
	      hDRratio1->SetPoint(mm,D2_Dx[nn],tmp2);
              mm++;
           }
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
   hdata->SetMarkerColor(1);
   hdata->SetMarkerStyle(8);
   mg1->Add(hborn);
   mg1->Add(hborn1);
   mg1->Add(hdata);
   mg1->Draw("AP");
   mg1->SetTitle("D2 born cross section;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.811,0.811);
   leg1->AddEntry(hborn,"Bodek","P");
   leg1->AddEntry(hborn1,"f1f217","P");
   leg1->AddEntry(hdata,"Data","P");
   leg1->Draw();

   c1->cd(2);
   hratio->SetMarkerStyle(8);
   hratio->SetMarkerColor(4);
   hratio->Draw("AP");
   hratio->SetTitle("D2 Bodek/f1f217;xbj;"); 

   TCanvas *c2=new TCanvas("c2","c2",1500,1500);
   TMultiGraph *mg2=new TMultiGraph();
   hDRratio->SetMarkerStyle(8);
   hDRratio->SetMarkerColor(4);
   hDRratio1->SetMarkerStyle(8);
   hDRratio1->SetMarkerColor(2);
   mg2->Add(hDRratio);
   mg2->Add(hDRratio1);
   mg2->Draw("AP");
//   mg2->SetTitle("");

   auto leg2=new TLegend(0.7,0.6,0.811,0.811);
   leg2->AddEntry(hDRratio,"Bodek/Data","P");
   leg2->AddEntry(hDRratio1,"f1f217/Data","P");
   leg2->Draw();

}
