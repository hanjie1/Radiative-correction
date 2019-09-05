#define MAXBIN 75
void plot_Compare()
{
     TString filename;
     ifstream file1;
     cout<<"Input file name: ";
     cin>>filename;
     file1.open(Form("../OUT/%s",filename.Data()));

     Double_t xbj[MAXBIN]={0.0},Q2[MAXBIN]={0.0},W2[MAXBIN]={0.0};
     Double_t F2D_B[MAXBIN]={0.0},F2D_W[MAXBIN]={0.0};
     Double_t F2D_NMC[MAXBIN]={0.0},F2D_CJ[MAXBIN]={0.0};
     Double_t Ratio[MAXBIN]={0.0};
     Double_t Ratio1[MAXBIN]={0.0};
     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           xbj[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           W2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_W[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_B[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_NMC[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_CJ[nn-1]=atof(content.Data());

	   Ratio[nn-1]=F2D_NMC[nn-1]/F2D_B[nn-1];
	   Ratio1[nn-1]=F2D_B[nn-1]/F2D_CJ[nn-1];
           nn++;
           from=0;
         }
     file1.close();

     TGraphErrors *gW=new TGraphErrors(MAXBIN,xbj,F2D_W,0,0);
     TGraphErrors *gB=new TGraphErrors(MAXBIN,xbj,F2D_B,0,0);
     TGraphErrors *gNMC=new TGraphErrors(MAXBIN,xbj,F2D_NMC,0,0);
     TGraphErrors *gCJ=new TGraphErrors(MAXBIN,xbj,F2D_CJ,0,0);

     TCanvas *c1=new TCanvas("c1");
     TMultiGraph *mg1=new TMultiGraph();
     gW->SetMarkerColor(2);
     gW->SetMarkerStyle(8);
     gB->SetMarkerColor(2);
     gB->SetMarkerStyle(8);
     gNMC->SetMarkerColor(4);
     gNMC->SetMarkerStyle(8);
     gCJ->SetMarkerColor(1);
     gCJ->SetMarkerStyle(22);
//     mg1->Add(gW);
     mg1->Add(gB);
     mg1->Add(gNMC);
//     mg1->Add(gCJ);
     mg1->Draw("AL");

     auto leg1=new TLegend(0.2,0.6,0.38,0.78);
     leg1->AddEntry(gB,"Bodek et al.","P");
     leg1->AddEntry(gNMC,"NMC","P");
     leg1->Draw();

     TCanvas *c2=new TCanvas("c2");
     TGraph *gRatio=new TGraph(MAXBIN,xbj,Ratio);
     gRatio->Draw("AP*");
}
