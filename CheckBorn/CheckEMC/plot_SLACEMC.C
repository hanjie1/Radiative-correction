#define MAXBIN 80
void plot_SLACEMC()
{
     Double_t x[MAXBIN],EMC_ISO[MAXBIN],EMC[MAXBIN];
     Double_t KP_EMC[MAXBIN];

     int Z=1;

     ifstream file1;
     file1.open("OUT/SLAC_EMC_H3_CJ15.out");

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

     TGraph *gEMC=new TGraph(nn-1,x,EMC);
     TGraph *gEMC_KP=new TGraph(nn-1,x,KP_EMC);
     TGraph *gEMC_ISO=new TGraph(nn-1,x,EMC_ISO);

     TMultiGraph *mg=new TMultiGraph();
     gEMC->SetMarkerColor(2);
     gEMC->SetMarkerStyle(8);
     gEMC_KP->SetMarkerColor(4);
     gEMC_KP->SetMarkerStyle(8);
     mg->Add(gEMC);
     mg->Add(gEMC_KP);
     mg->Draw("AP");
     if(Z==1)mg->SetTitle("H3 EMC;x;");
     if(Z==2)mg->SetTitle("He3 EMC;x;");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(gEMC,"SLAC EMC+CJ15","P");
   leg1->AddEntry(gEMC_KP,"K&P","P");
   leg1->Draw();


}
