int MAXBIN=100;
void plot_QE()
{
     TString file1;
     file1="../OUT/RESCSD/D2_all.out";

     ifstream infile1;
     infile1.open(file1);

     Double_t xbj[MAXBIN],Q2[MAXBIN],XS_QE[MAXBIN],XS_DIS[MAXBIN];
     Double_t ratio[MAXBIN];
     int nn=0;
     Ssiz_t from=0;
     TString content,tmp;
     while(tmp.ReadLine(infile1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           xbj[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           XS_QE[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           XS_DIS[nn-1]=atof(content.Data());
	   ratio[nn-1]=XS_QE[nn-1]/XS_DIS[nn-1];
           nn++;
           from=0;
         }
     infile1.close();

     TGraph *h1=new TGraph(nn-1,xbj,ratio);
     h1->SetMarkerStyle(8);
     h1->SetMarkerColor(4);
     h1->Draw("AP");
     h1->SetTitle("SIG_RAD_QE/SIG_RAD_DIS;xbj;");


}
