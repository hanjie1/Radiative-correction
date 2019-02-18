#define MAXBIN 41
void plot_NMC()
{
     Double_t xx[MAXBIN],Q2[MAXBIN],WSQ[MAXBIN];
     Double_t F2_dis[MAXBIN],F2_res[MAXBIN];

     ifstream file1;
     file1.open("F2D_NMC.dat");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           tmp.Tokenize(content,from," ");
           xx[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           WSQ[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2_dis[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2_res[nn]=atof(content.Data());
           nn++;
           from=0;
         }
     file1.close();

      TGraph *gF2_dis=new TGraph(nn,WSQ,F2_dis);
      TGraph *gF2_res=new TGraph(nn,WSQ,F2_res);

      TCanvas *c1=new TCanvas("c1");
      c1->Divide(2,1);
      c1->cd(1);
      gF2_res->Draw("AP*");
      c1->cd(2);
      gF2_dis->Draw("AP*");

}
