void plot_F2(){
     TString f1="gsmearing_test.dat";
     TString f2="Bodek_test.dat";
     ifstream infile1;
     infile1.open(f1);
     ifstream infile2;
     infile2.open(f2);

     Double_t xbj[22]={0.0};
     Double_t F2_Bodek[22]={0.0}; 
     Double_t F2_gsmearing[22]={0.0}; 

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(infile1)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          xbj[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2_gsmearing[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
     infile1.close();

     from=0;
     nn=0;
     while(tmp.ReadLine(infile2)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          xbj[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2_Bodek[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
     infile2.close();
 
     TGraph *gGsmearing=new TGraph(22,xbj,F2_gsmearing); 
     TGraph *gBodek=new TGraph(22,xbj,F2_Bodek);

     TCanvas *c1=new TCanvas("c1");
     TMultiGraph *mg1=new TMultiGraph();
     gGsmearing->SetMarkerColor(2);
     gGsmearing->SetMarkerStyle(8);
     gBodek->SetMarkerColor(1);
     gBodek->SetMarkerStyle(8);
     mg1->Add(gGsmearing);
     mg1->Add(gBodek);
     mg1->Draw("AP"); 
     mg1->SetTitle("F2p;xbj;F2p");

     auto leg1=new TLegend(0.7,0.6,0.85,0.85);
     leg1->AddEntry(gGsmearing,"gsmearing","P");
     leg1->AddEntry(gBodek,"Bodek","P");
     leg1->Draw();

}
