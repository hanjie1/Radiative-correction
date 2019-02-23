#define MAXBIN 98
void plot_EMC()
{
     Double_t x[MAXBIN],EMC[MAXBIN];

     ifstream file1;
     file1.open("OUT/EMC_He3.out");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           x[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC[nn-1]=atof(content.Data());

           nn++;
           from=0;
         }
     file1.close();

     TGraph *gEMC=new TGraph(nn,x,EMC);
     gEMC->Draw("AP*");


}
