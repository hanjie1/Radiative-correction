#define MAXBIN 80
void plot_SLACEMC()
{
     Double_t x[MAXBIN],EMC_ISO[MAXBIN],EMC_CJ[MAXBIN];
     Double_t KP_EMC[MAXBIN],EMC_SLAC[MAXBIN],EMC_NMC[MAXBIN];

     int Z=2;

     ifstream file1;
     file1.open("OUT/Compare_EMC_He3.out");

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
           EMC_CJ[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC_NMC[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           EMC_SLAC[nn-1]=atof(content.Data());

//	   EMC[nn-1]=EMC[nn-1]*2.0/3.0;

	   if(Z==1)
	      KP_EMC[nn-1]= 1.07251-1.61648*x[nn-1]+12.3626*pow(x[nn-1],2)-65.5932*pow(x[nn-1],3)+213.311*pow(x[nn-1],4)
                -423.943*pow(x[nn-1],5)+499.994*pow(x[nn-1],6)-321.304*pow(x[nn-1],7)+87.0596*pow(x[nn-1],8);
	   if(Z==2)
 	      KP_EMC[nn-1]= 1.02967-0.135929*x[nn-1]+3.92009*pow(x[nn-1],2)-21.2861*pow(x[nn-1],3)
                   +64.7762*pow(x[nn-1],4)-129.928*pow(x[nn-1],5)+169.609*pow(x[nn-1],6)
                   -127.386*pow(x[nn-1],7)+41.0723*pow(x[nn-1],8);

	   KP_EMC[nn-1]=KP_EMC[nn-1]*3.0/2.0;

           nn++;
           from=0;
         }
     file1.close();

     TGraph *gEMC_CJ=new TGraph(nn-1,x,EMC_CJ);
     TGraph *gEMC_NMC=new TGraph(nn-1,x,EMC_NMC);
     TGraph *gEMC_SLAC=new TGraph(nn-1,x,EMC_SLAC);
     TGraph *gEMC_KP=new TGraph(nn-1,x,KP_EMC);
     TGraph *gEMC_ISO=new TGraph(nn-1,x,EMC_ISO);

     TMultiGraph *mg=new TMultiGraph();
     gEMC_CJ->SetMarkerColor(4);
     gEMC_CJ->SetMarkerStyle(8);
     gEMC_CJ->SetLineWidth(3);

     gEMC_SLAC->SetMarkerColor(8);
     gEMC_SLAC->SetMarkerStyle(8);
     gEMC_SLAC->SetLineWidth(3);

     gEMC_NMC->SetMarkerColor(kYellow+1);
     gEMC_NMC->SetMarkerStyle(8);
     gEMC_NMC->SetLineWidth(3);

     gEMC_KP->SetMarkerColor(2);
     gEMC_KP->SetMarkerStyle(8);
     gEMC_KP->SetLineWidth(3);
     mg->Add(gEMC_CJ);
     mg->Add(gEMC_NMC);
     mg->Add(gEMC_SLAC);
     mg->Add(gEMC_KP);
     mg->Draw("AP");
     if(Z==1)mg->SetTitle("H3 EMC;x;");
     if(Z==2)mg->SetTitle("He3 EMC;x;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(gEMC_KP,"K&P","P");
   leg1->AddEntry(gEMC_CJ,"SLAC EMC+CJ15","P");
   leg1->AddEntry(gEMC_NMC,"SLAC EMC+NMC","P");
   leg1->AddEntry(gEMC_SLAC,"SLAC EMC+SLAC","P");
   leg1->Draw();


}
