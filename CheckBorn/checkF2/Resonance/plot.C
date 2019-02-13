#define MAXBIN 41
void plot()
{
     Double_t xx[MAXBIN],Q2[MAXBIN],WSQ[MAXBIN],BKG[MAXBIN],RES[MAXBIN],UNIV[MAXBIN];
     Double_t F2_dis[MAXBIN],F2_res[MAXBIN];

     ifstream file1;
     file1.open("F2D.dat");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           tmp.Tokenize(content,from," ");
           BKG[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           RES[nn]=atof(content.Data());
	
           from=0;
           tmp.ReadLine(file1);
           tmp.Tokenize(content,from," ");
           xx[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           WSQ[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           UNIV[nn]=atof(content.Data());

	   F2_dis[nn]=BKG[nn]*UNIV[nn];
	   F2_res[nn]=BKG[nn]+RES[nn];


           nn++;
           from=0;
         }
     file1.close();

      TGraph *gUNIV=new TGraph(nn,WSQ,UNIV);
      TGraph *gF2_dis=new TGraph(nn,WSQ,F2_dis);
      TGraph *gRES=new TGraph(nn,WSQ,RES);
      TGraph *gF2_res=new TGraph(nn,WSQ,F2_res);

      TCanvas *c1=new TCanvas("c1");
      c1->Divide(2,1);
      c1->cd(1);
      gUNIV->Draw("AP*");
      c1->cd(2);
      gF2_dis->Draw("AP*");

      TCanvas *c2=new TCanvas("c2");
      c2->Divide(2,1);
      c2->cd(1);
      gRES->Draw("AP*");
      c2->cd(2);
      gF2_res->Draw("AP*");
}
