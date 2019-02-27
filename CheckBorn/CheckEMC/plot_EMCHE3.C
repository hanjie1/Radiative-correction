#define MAXBIN 28
void plot_EMCHE3()
{
     Double_t x[MAXBIN],EMC_HALLC[MAXBIN],Stat[MAXBIN],Sys[MAXBIN],Fis[MAXBIN];
     Double_t EMC_model[MAXBIN],EMC_KP[MAXBIN],Derr[MAXBIN];

     ifstream file1;
     file1.open("OUT/CheckA3_He3_18deg.out");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           x[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           EMC_HALLC[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Stat[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Sys[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Fis[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC_model[nn-1]=atof(content.Data());

	   Derr[nn-1]=sqrt(Sys[nn-1]*Sys[nn-1]+Stat[nn-1]*Stat[nn-1]);
	   EMC_HALLC[nn-1]=EMC_HALLC[nn-1]/Fis[nn-1];
           EMC_KP[nn-1]=1.02967-0.135929*x[nn-1]+3.92009*pow(x[nn-1],2)-21.2861*pow(x[nn-1],3)
                +64.7762*pow(x[nn-1],4)-129.928*pow(x[nn-1],5)+169.609*pow(x[nn-1],6)
                -127.386*pow(x[nn-1],7)+41.0723*pow(x[nn-1],8);

           nn++;
           from=0;
         }
     file1.close();

     TGraphErrors *gEMC=new TGraphErrors(nn-1,x,EMC_HALLC,0,Derr);
     TGraph *gEMC_model=new TGraph(nn-1,x,EMC_model);
     TGraph *gEMC_KP=new TGraph(nn-1,x,EMC_KP);

     TMultiGraph *mg=new TMultiGraph();
     gEMC->SetMarkerColor(2);
     gEMC->SetMarkerStyle(8);
     gEMC_model->SetMarkerColor(4);
     gEMC_model->SetMarkerStyle(8);
     gEMC_KP->SetMarkerColor(1);
     gEMC_KP->SetMarkerStyle(8);
     mg->Add(gEMC);
     mg->Add(gEMC_model);
     mg->Add(gEMC_KP);
     mg->Draw("AP");
     mg->SetTitle("18deg");

   auto leg=new TLegend(0.7,0.6,0.85,0.85);
   leg->AddEntry(gEMC,"DATA","P");
   leg->AddEntry(gEMC_model,"NO RES","P");
   leg->AddEntry(gEMC_KP,"K&P","P");
   leg->Draw();



}
