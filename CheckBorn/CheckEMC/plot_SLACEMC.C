#define MAXBIN 99
void plot_SLACEMC()
{
     Double_t x[MAXBIN],EMC_ISO[MAXBIN],EMC[MAXBIN];

     ifstream file1;
     file1.open("OUT/SLAC_EMC_H3.out");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           x[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC_ISO[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC[nn-1]=atof(content.Data());

	   EMC[nn-1]=EMC[nn-1]*2.0/3.0;
           nn++;
           from=0;
         }
     file1.close();

     TGraph *gEMC=new TGraph(nn-1,x,EMC);
     TGraph *gEMC_ISO=new TGraph(nn-1,x,EMC_ISO);

     TMultiGraph *mg=new TMultiGraph();
     gEMC->SetMarkerColor(2);
     gEMC->SetMarkerStyle(8);
     gEMC_ISO->SetMarkerColor(4);
     gEMC_ISO->SetMarkerStyle(8);
     mg->Add(gEMC);
     mg->Add(gEMC_ISO);
     mg->Draw("AP");


}
