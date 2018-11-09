void plot_chi2()
{
     ifstream infile1;
     infile1.open("chi2_4He0.out");
     ifstream infile2;
     infile2.open("chi2_4He1.out");
     ifstream infile3;
     infile3.open("chi2_4He2.out");
     ifstream infile4;
     infile4.open("chi2_4He3.out");
     ifstream infile5;
     infile5.open("chi2_4He4.out");
     ifstream infile6;
     infile6.open("chi2_4He5.out");
     ifstream infile7;
     infile7.open("chi2_4He6.out");
     ifstream infile8;
     infile8.open("chi2_4He7.out");
     ifstream infile9;
     infile9.open("chi2_4He8.out");
     ifstream infile10;
     infile10.open("chi2_4He9.out");

     Double_t kf[21],es[21],xvalt[21];
     Double_t chi0[21],chi1[21],chi2[21],chi3[21],chi4[21],chi5[21],chi6[21],chi7[21],chi8[21],chi9[21];

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(infile1)){
           tmp.Tokenize(content,from," ");
           kf[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           es[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           xvalt[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           chi0[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile1.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile2)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi1[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile2.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile3)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi2[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile3.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile4)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi3[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile4.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile5)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi4[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile5.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile6)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi5[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile6.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile7)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi6[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile7.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile8)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi7[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile8.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile9)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi8[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile9.close();

     nn=0;from=0;
     while(tmp.ReadLine(infile10)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           chi9[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile10.close();

     Double_t total_chi2[21]={0.0};
     for(int ii=0;ii<21;ii++){
         total_chi2[ii]=chi0[ii]+chi1[ii]+chi2[ii]+chi3[ii]+chi4[ii]+chi5[ii]+chi6[ii]+chi7[ii]+chi8[ii]+chi9[ii];
     }

     TGraph *htotal=new TGraph(21,kf,total_chi2);
     TCanvas *c0=new TCanvas("c0","c0",1500,1500);
     htotal->SetMarkerStyle(8);
     htotal->SetMarkerColor(4);
     htotal->Draw("AP");
     htotal->SetTitle("total chi2;kf;");


     TGraph *h0=new TGraph(21,kf,chi0);
     TGraph *h1=new TGraph(21,kf,chi1);
     TGraph *h2=new TGraph(21,kf,chi2);
     TGraph *h3=new TGraph(21,kf,chi3);
     TGraph *h4=new TGraph(21,kf,chi4);
     TGraph *h5=new TGraph(21,kf,chi5);
     TGraph *h6=new TGraph(21,kf,chi6);
     TGraph *h7=new TGraph(21,kf,chi7);
     TGraph *h8=new TGraph(21,kf,chi8);
     TGraph *h9=new TGraph(21,kf,chi9);

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     c1->Divide(5,2);
     c1->cd(1);
     h0->SetMarkerStyle(8);
     h0->SetMarkerColor(4);
     h0->Draw("AP");
     h0->SetTitle("data set 0;kf;");

     c1->cd(2);
     h1->SetMarkerStyle(8);
     h1->SetMarkerColor(4);
     h1->Draw("AP");
     h1->SetTitle("data set 1;kf;");

     c1->cd(3);
     h2->SetMarkerStyle(8);
     h2->SetMarkerColor(4);
     h2->Draw("AP");
     h2->SetTitle("data set 2;kf;");

     c1->cd(4);
     h3->SetMarkerStyle(8);
     h3->SetMarkerColor(4);
     h3->Draw("AP");
     h3->SetTitle("data set 3;kf;");

     c1->cd(5);
     h4->SetMarkerStyle(8);
     h4->SetMarkerColor(4);
     h4->Draw("AP");
     h4->SetTitle("data set 4;kf;");

     c1->cd(6);
     h5->SetMarkerStyle(8);
     h5->SetMarkerColor(4);
     h5->Draw("AP");
     h5->SetTitle("data set 5;kf;");

     c1->cd(7);
     h6->SetMarkerStyle(8);
     h6->SetMarkerColor(4);
     h6->Draw("AP");
     h6->SetTitle("data set 6;kf;");

     c1->cd(8);
     h7->SetMarkerStyle(8);
     h7->SetMarkerColor(4);
     h7->Draw("AP");
     h7->SetTitle("data set 7;kf;");

     c1->cd(9);
     h8->SetMarkerStyle(8);
     h8->SetMarkerColor(4);
     h8->Draw("AP");
     h8->SetTitle("data set 8;kf;");

     c1->cd(10);
     h9->SetMarkerStyle(8);
     h9->SetMarkerColor(4);
     h9->Draw("AP");
     h9->SetTitle("data set 9;kf;");






}
