#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_D2()
{
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Rad[11][MAXBIN],D2_Born1[11][MAXBIN],D2_Rad1[11][MAXBIN];
     Double_t D2_Born2[11][MAXBIN],D2_Rad2[11][MAXBIN];
     Double_t D2_RC[11][MAXBIN],D2_RC1[11][MAXBIN],D2_RC2[11][MAXBIN];
     Double_t D2_Yield[11][MAXBIN],D2_Yerr[11][MAXBIN];
     Double_t D2_Dborn[11][MAXBIN],D2_Derr[11][MAXBIN],D2_Dborn1[11][MAXBIN],D2_Derr1[11][MAXBIN];
     Double_t D2_ratio[11][MAXBIN];

     for(int ii=0;ii<11;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
	     D2_Rad2[ii][jj]=0.0; D2_Born2[ii][jj]=0.0;D2_RC2[ii][jj]=0.0;
             D2_RC[ii][jj]=0.0;D2_RC1[ii][jj]=0.0;
	     D2_Yield[ii][jj]=0.0; D2_Yerr[ii][jj]=0.0;
	     D2_ratio[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       Yfile=Form("model111/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad); 
       Yfile=Form("model211/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1); 
       Yfile=Form("model211_1/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born2,D2_Rad2); 
   }

   ifstream Datafile;
   Ssiz_t from=0;
   TString content,tmp;
   for(int ii=0;ii<11;ii++){
       int nn=0;
       from=0;
       TString filename=Form("/home/hanjie/work/MARATHON/RadCor/T2_externals/RUNPLAN/bin_avg/D2_kin%d.txt",kin[ii]); 
       Datafile.open(filename);
       while(tmp.ReadLine(Datafile)){
             if(nn==0){nn++;continue;}
             tmp.Tokenize(content,from,", ");
             tmp.Tokenize(content,from,", ");
             tmp.Tokenize(content,from,", ");
             tmp.Tokenize(content,from,", ");
	     D2_Yield[ii][nn-1]=atof(content.Data());
             tmp.Tokenize(content,from,", ");
	     D2_Yerr[ii][nn-1]=atof(content.Data());
             nn++;
             from=0;
       }
       Datafile.close();
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
    
   ofstream myfile;
   myfile.open("DataBorn/D2.out");
   myfile<<"x   Q2   Born_ineft   err   Born_f1f2   err  mark "<<endl;
   int nn=0;
   for(int ii=0;ii<11;ii++){
       for(int jj=0;jj<MAXBIN;jj++){
	   if(D2_x[ii][jj]==0)continue;
           hborn->SetPoint(nn,D2_x[ii][jj],D2_Born[ii][jj]);
           hrad->SetPoint(nn,D2_x[ii][jj],D2_Rad[ii][jj]);
           if(D2_Rad[ii][jj]>0)D2_RC[ii][jj]=D2_Born[ii][jj]/D2_Rad[ii][jj];
           hRC->SetPoint(nn,D2_x[ii][jj],D2_RC[ii][jj]);
           //hRC->SetPoint(nn,D2_x[ii][jj],D2_Born[ii][jj]/D2_Born1[ii][jj]);

           hborn1->SetPoint(nn,D2_x[ii][jj],D2_Born1[ii][jj]);
           hrad1->SetPoint(nn,D2_x[ii][jj],D2_Rad1[ii][jj]);
           if(D2_Rad1[ii][jj]>0)D2_RC1[ii][jj]=D2_Born1[ii][jj]/D2_Rad1[ii][jj];
           hRC1->SetPoint(nn,D2_x[ii][jj],D2_RC1[ii][jj]);
         //  hRC1->SetPoint(nn,D2_x[ii][jj],D2_Rad[ii][jj]/D2_Rad1[ii][jj]);

           if(D2_RC1[ii][jj]>0)D2_ratio[ii][jj]=D2_RC[ii][jj]/D2_RC1[ii][jj];
           hratio->SetPoint(nn,D2_x[ii][jj],D2_ratio[ii][jj]);

           hborn2->SetPoint(nn,D2_x[ii][jj],D2_Born2[ii][jj]);
           hrad2->SetPoint(nn,D2_x[ii][jj],D2_Rad2[ii][jj]);
           if(D2_Rad2[ii][jj]>0)D2_RC2[ii][jj]=D2_Born2[ii][jj]/D2_Rad2[ii][jj];
           hRC2->SetPoint(nn,D2_x[ii][jj],D2_RC2[ii][jj]);

           if(D2_RC2[ii][jj]>0)
               hratio1->SetPoint(nn,D2_x[ii][jj],D2_RC[ii][jj]/D2_RC2[ii][jj]);

	   D2_Dborn[ii][jj]=D2_Yield[ii][jj]*D2_RC[ii][jj];
	   D2_Derr[ii][jj]=D2_Yerr[ii][jj]*D2_RC[ii][jj];
	   D2_Dborn1[ii][jj]=D2_Yield[ii][jj]*D2_RC1[ii][jj];
	   D2_Derr1[ii][jj]=D2_Yerr[ii][jj]*D2_RC1[ii][jj];
           nn++;
       }
   } 
  
   int mark[11][MAXBIN];
   for(int ii=0;ii<11;ii++){
       Double_t tmpmax=0.0;
       int maxii=0,maxjj=0;
       for(int jj=0;jj<MAXBIN;jj++){
           mark[ii][jj]=0;
	   if(D2_x[ii][jj]==0)continue;
           if(D2_Dborn[ii][jj]>tmpmax){
              tmpmax=D2_Dborn[ii][jj];
	      maxii=ii;maxjj=jj;
           }
       }
       mark[maxii][maxjj]=1;
   }

   for(int ii=0;ii<11;ii++){
	for(int jj=0;jj<MAXBIN;jj++){
            if(D2_x[ii][jj]==0)continue;
            myfile<<D2_x[ii][jj]<<"  "<<D2_Q2[ii][jj]<<"  "<<D2_Dborn[ii][jj]<<"  "<<D2_Derr[ii][jj]<<"  "<<D2_Dborn1[ii][jj]<<"  "<<D2_Derr1[ii][jj]<<"  "<<mark[ii][jj]<<endl;

	}
   }
   myfile.close();


   TCanvas *c1=new TCanvas("c1","c1",1500,1500);
   c1->Divide(2,1);
   c1->cd(1);
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
   mg1->SetTitle("D2 born cross section;xbj;born");

   auto leg1=new TLegend(0.7,0.6,0.811,0.811);
   leg1->AddEntry(hborn,"model111","P");
   leg1->AddEntry(hborn1,"model211","P");
   leg1->AddEntry(hborn2,"model211_1","P");
   leg1->Draw();

   c1->cd(2);
   TMultiGraph *mg2=new TMultiGraph();
   hrad->SetMarkerStyle(8);
   hrad->SetMarkerColor(2);
   hrad1->SetMarkerStyle(8);
   hrad1->SetMarkerColor(4);
   hrad2->SetMarkerStyle(8);
   hrad2->SetMarkerColor(1);
   mg2->Add(hrad);
   mg2->Add(hrad1);
   mg2->Add(hrad2);
   mg2->Draw("AP");
   mg2->SetTitle("D2 rad cross section;xbj;rad");

   auto leg2=new TLegend(0.7,0.6,0.811,0.811);
   leg2->AddEntry(hrad,"model111","P");
   leg2->AddEntry(hrad1,"model211","P");
   leg2->AddEntry(hrad2,"model211_1","P");
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
   mg3->SetTitle("D2 RC=born/rad;xbj;RC");

   auto leg3=new TLegend(0.7,0.6,0.811,0.811);
   leg3->AddEntry(hRC,"model111","P");
   leg3->AddEntry(hRC1,"model211","P");
   leg3->AddEntry(hRC2,"model211_1","P");
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
   mg4->SetTitle("D2 model111/model211;xbj;"); 

   auto leg4=new TLegend(0.7,0.6,0.811,0.811);
   leg4->AddEntry(hratio,"model111/model211","P");
   leg4->AddEntry(hratio1,"model111/model211_1","P");
   leg4->Draw();

}
