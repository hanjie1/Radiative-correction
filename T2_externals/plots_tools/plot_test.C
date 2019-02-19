#define MAXBIN 27001
void plot_test()
{
     Double_t xx[MAXBIN],Q2[MAXBIN],WSQ[MAXBIN];
     Double_t F2[MAXBIN],EMC[MAXBIN],F2_final[MAXBIN];

     ifstream file1;
     file1.open("../HE3_out.dat");

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
           F2[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC[nn]=atof(content.Data());

	   F2_final[nn]=F2[nn]*EMC[nn];
           nn++;
           from=0;
         }
     file1.close();

     TGraph *gEMC=new TGraph(nn,xx,EMC);
     TGraph *gF2=new TGraph(nn,xx,F2);
     TGraph *gF2_final=new TGraph(nn,xx,F2_final);

     TCanvas *c1=new TCanvas("c1");
     c1->Divide(2,2);
     c1->cd(1);
     gEMC->Draw("AP*");
     c1->cd(2);
     gF2->Draw("AP*");
     c1->cd(3);
     gF2_final->Draw("AP*");


}
