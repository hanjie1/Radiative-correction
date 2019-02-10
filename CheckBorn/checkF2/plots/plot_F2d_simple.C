#define MAXBIN 610
void plot_F2d_simple()
{
     TString filename;
     ifstream file1;
     cout<<"Input file name: ";
     cin>>filename;
     file1.open(Form("../OUT/%s",filename.Data()));

     Double_t xbj[MAXBIN]={0.0},Q2[MAXBIN]={0.0},W2[MAXBIN]={0.0};
     Double_t F2D[MAXBIN]={0.0},F2D_stat[MAXBIN]={0.0},F2D_sys[MAXBIN]={0.0};
     Double_t F2E[MAXBIN]={0.0},F2B[MAXBIN]={0.0};
     Double_t F2D_err[MAXBIN]={0.0};
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
           F2D[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_stat[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_sys[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2E[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2B[nn-1]=atof(content.Data());

           F2D_err[nn-1]=F2D[nn-1]*sqrt(F2D_stat[nn-1]*F2D_stat[nn-1]+F2D_sys[nn-1]*F2D_sys[nn-1]);
           nn++;
           from=0;
         }
     file1.close();

     TGraphErrors *gData=new TGraphErrors(MAXBIN,xbj,F2D,0,F2D_err);
     TGraphErrors *gNMC=new TGraphErrors(MAXBIN,xbj,F2E,0,0);
     TGraphErrors *gINEFT=new TGraphErrors(MAXBIN,xbj,F2B,0,0);

     TCanvas *c1=new TCanvas("c1");
     TMultiGraph *mg1=new TMultiGraph();
     gNMC->SetMarkerColor(2);
     gNMC->SetMarkerStyle(8);
     gINEFT->SetMarkerColor(4);
     gINEFT->SetMarkerStyle(22);
     gData->SetMarkerColor(1);
     gData->SetMarkerStyle(22);
     mg1->Add(gNMC);
     mg1->Add(gINEFT);
     mg1->Add(gData);
     mg1->Draw("AP");

}
