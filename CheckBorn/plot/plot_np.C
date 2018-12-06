#define MAXNUM 50
void plot_np()
{
    ifstream file1,file2;
    TString myfile1="../OUT/NP/marathon_ineft.out";
    file1.open(myfile1);
    if(!file1.is_open())return 0;

    TString myfile2="../OUT/NP/marathon_f1f217.out";
    file2.open(myfile2);
    if(!file2.is_open())return 0;

    Double_t xbj[MAXNUM],Q2[MAXNUM];
    Double_t f1n[MAXNUM],f2n[MAXNUM],f1p[MAXNUM],f2p[MAXNUM],f1d[MAXNUM],f2d[MAXNUM];
    Double_t f1n1[MAXNUM],f2n1[MAXNUM],f1p1[MAXNUM],f2p1[MAXNUM],f1d1[MAXNUM],f2d1[MAXNUM];

    for(int ii=0;ii<MAXNUM;ii++){
        xbj[ii]=0.0; Q2[ii]=0.0;
        f1n[ii]=0.0;f2n[ii]=0.0;f1p[ii]=0.0;f2p[ii]=0.0;f1d[ii]=0.0;f2d[ii]=0.0;
        f1n1[ii]=0.0;f2n1[ii]=0.0;f1p1[ii]=0.0;f2p1[ii]=0.0;f1d1[ii]=0.0;f2d1[ii]=0.0;
    }


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
	  f1p[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2p[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f1n[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2n[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f1d[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2d[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    file1.close();

    nn=0;from=0;
    while(tmp.ReadLine(file2)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
	  f1p1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2p1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f1n1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2n1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f1d1[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
	  f2d1[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    file2.close();

    TGraph *hnp=new TGraph();
    TGraph *hnp1=new TGraph();
    TGraph *hdp=new TGraph();
    TGraph *hdp1=new TGraph();
    TGraph *hp=new TGraph();
    TGraph *hp1=new TGraph();
    TGraph *hn=new TGraph();
    TGraph *hn1=new TGraph();
    TGraph *hd=new TGraph();
    TGraph *hd1=new TGraph();


    for(int ii=0;ii<nn-1;ii++){
        hnp->SetPoint(ii,xbj[ii],f2n[ii]/f2p[ii]);
        hnp1->SetPoint(ii,xbj[ii],f2n1[ii]/f2p1[ii]);

        hdp->SetPoint(ii,xbj[ii],f2d[ii]/f2p[ii]/2.0);
        hdp1->SetPoint(ii,xbj[ii],f2d1[ii]/f2p1[ii]/2.0);

        hp->SetPoint(ii,xbj[ii],f2p[ii]);
        hp1->SetPoint(ii,xbj[ii],f2p1[ii]);

        hn->SetPoint(ii,xbj[ii],f2n[ii]);
        hn1->SetPoint(ii,xbj[ii],f2n1[ii]);
    }

    TCanvas *c1=new TCanvas("c1","c1",1500,1500);
    c1->Divide(2,1);
    c1->cd(1);
    TMultiGraph *mg1=new TMultiGraph();
    hnp->SetMarkerStyle(8);
    hnp->SetMarkerColor(4);
    hnp1->SetMarkerStyle(8);
    hnp1->SetMarkerColor(2);
    mg1->Add(hnp);
    mg1->Add(hnp1);
    mg1->Draw("AP");
    mg1->SetTitle("F2n/F2p;xbj;");

    auto leg1=new TLegend(0.7,0.6,0.811,0.811);
    leg1->AddEntry(hnp,"Bodek","P");
    leg1->AddEntry(hnp1,"f1f217","P");
    leg1->Draw();

    c1->cd(2);
    TMultiGraph *mg2=new TMultiGraph();
    hdp->SetMarkerStyle(8);
    hdp->SetMarkerColor(4);
    hdp1->SetMarkerStyle(8);
    hdp1->SetMarkerColor(2);
    mg2->Add(hdp);
    mg2->Add(hdp1);
    mg2->Draw("AP");
    mg2->SetTitle("F2d/F2p;xbj;");

    auto leg2=new TLegend(0.7,0.6,0.811,0.811);
    leg2->AddEntry(hdp,"Bodek","P");
    leg2->AddEntry(hdp1,"f1f217","P");
    leg2->Draw();

    TCanvas *c2=new TCanvas("c2","c2",1500,1500);
    c2->Divide(2,1);
    c2->cd(1);
    TMultiGraph *mg3=new TMultiGraph();
    hp->SetMarkerStyle(8);
    hp->SetMarkerColor(4);
    hp1->SetMarkerStyle(8);
    hp1->SetMarkerColor(2);
    mg3->Add(hp);
    mg3->Add(hp1);
    mg3->Draw("AP");
    mg3->SetTitle("F2p;xbj;");

    auto leg3=new TLegend(0.7,0.6,0.811,0.811);
    leg3->AddEntry(hp,"Bodek","P");
    leg3->AddEntry(hp1,"f1f217","P");
    leg3->Draw();

    c2->cd(2);
    TMultiGraph *mg4=new TMultiGraph();
    hn->SetMarkerStyle(8);
    hn->SetMarkerColor(4);
    hn1->SetMarkerStyle(8);
    hn1->SetMarkerColor(2);
    mg4->Add(hn);
    mg4->Add(hn1);
    mg4->Draw("AP");
    mg4->SetTitle("F2n;xbj;");

    auto leg4=new TLegend(0.7,0.6,0.811,0.811);
    leg4->AddEntry(hn,"Bodek","P");
    leg4->AddEntry(hn1,"f1f217","P");
    leg4->Draw();







}
