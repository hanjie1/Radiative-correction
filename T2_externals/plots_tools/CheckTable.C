#define MAXBIN 3000
void CheckTable()
{
     TString file1;
     file1="../OUT/XStable/Bodek/Carbon_kin1_xs.out";

     ifstream infile1;
     infile1.open(file1); 

     Double_t xbj[MAXBIN],Q2[MAXBIN],Theta[MAXBIN],Ep[MAXBIN],XS_Born[MAXBIN],XS_Rad[MAXBIN];
     int nn=0;
     Ssiz_t from=0;
     TString content,tmp;
     while(tmp.ReadLine(infile1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           xbj[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Theta[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Ep[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           XS_Born[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           XS_Rad[nn-1]=atof(content.Data());
           nn++;
           from=0;
         }
     infile1.close();

     TCanvas *c1=new TCanvas("c1");
     TGraph *gBorn=new TGraph(nn-1,Ep,XS_Born);
     gBorn->SetMarkerStyle(8);
     gBorn->Draw("AP");
     TCanvas *c2=new TCanvas("c2");
     TGraph *gRad=new TGraph(nn-1,Ep,XS_Rad);
     gRad->SetMarkerStyle(8);
     gRad->Draw("AP");


}
