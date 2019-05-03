#define MAXBIN 81
void plot_NMC()
{
     Double_t xx[MAXBIN],Q2[MAXBIN],WSQ[MAXBIN];
     Double_t F2_NMC[MAXBIN],err_lo[MAXBIN],err_hi[MAXBIN];
     Double_t F2_B[MAXBIN],R[MAXBIN],DR[MAXBIN];

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
           F2_NMC[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           err_lo[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           err_hi[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2_B[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           R[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           DR[nn]=atof(content.Data());
           nn++;
           from=0;
         }
     file1.close();

      TGraphAsymmErrors *gF2_NMC=new TGraphAsymmErrors(nn-1,xx,F2_NMC,0,0,err_lo,err_hi);
      TGraph *gF2_B=new TGraph(nn,xx,F2_B);
      TGraphErrors *gR=new TGraphErrors(nn,xx,R,0,DR);

      TCanvas *c1=new TCanvas("c1");
      TMultiGraph *mg=new TMultiGraph();
      gF2_NMC->SetMarkerColor(2);
      gF2_NMC->SetMarkerStyle(8);
      gF2_B->SetMarkerColor(4);
      gF2_B->SetMarkerStyle(8);
      mg->Add(gF2_NMC);
      mg->Add(gF2_B);
      mg->Draw("AP");
      mg->SetTitle("F2d;x");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(gF2_NMC,"NMC","P");
   leg1->AddEntry(gF2_B,"Bodek","P");
   leg1->Draw();

	TCanvas *c2=new TCanvas("c2");
	gR->Draw("AP*");
}
