#define MAXBIN 98
void plot_F2NP()
{
     Double_t x[MAXBIN],F2NP_SLAC[MAXBIN],F2NP_CJ[MAXBIN],F2NP_NMC[MAXBIN];

     ifstream file1;
     file1.open("OUT/F2NP.out");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           x[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2NP_SLAC[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2NP_CJ[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2NP_NMC[nn-1]=atof(content.Data());

           nn++;
           from=0;
         }
     file1.close();

     TGraph *gF2NP_CJ=new TGraph(nn-1,x,F2NP_CJ);
     TGraph *gF2NP_SLAC=new TGraph(nn-1,x,F2NP_SLAC);
     TGraph *gF2NP_NMC=new TGraph(nn-1,x,F2NP_NMC);

     TMultiGraph *mg=new TMultiGraph();
     gF2NP_CJ->SetMarkerColor(2);
     gF2NP_CJ->SetMarkerStyle(8);
     gF2NP_SLAC->SetMarkerColor(4);
     gF2NP_SLAC->SetMarkerStyle(8);
     gF2NP_NMC->SetMarkerColor(1);
     gF2NP_NMC->SetMarkerStyle(8);
     mg->Add(gF2NP_CJ);
     mg->Add(gF2NP_SLAC);
     mg->Add(gF2NP_NMC);
     mg->Draw("AP");


}
