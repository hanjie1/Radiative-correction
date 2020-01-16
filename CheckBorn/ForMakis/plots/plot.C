void plot(){
     TString filename;
     ifstream file1;
     filename="outfile.out";
     file1.open(Form("../%s",filename.Data()));

     Double_t xbj[18]={0.0},Q2[18]={0.0},W2[18]={0.0};
     Double_t F2Dp_f1f217[18]={0.0},F2Dp_ineft[18]={0.0};
     Double_t F2d_f1f217[18]={0.0},F2p_f1f217[18]={0.0};
     Double_t F2d_ineft[18]={0.0},F2p_ineft[18]={0.0};
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
           F2Dp_f1f217[nn-1]=atof(content.Data())/2.0;
           tmp.Tokenize(content,from," ");
           F2Dp_ineft[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2d_f1f217[nn-1]=atof(content.Data())/2.0;
           tmp.Tokenize(content,from," ");
           F2p_f1f217[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2d_ineft[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2p_ineft[nn-1]=atof(content.Data());
	   nn++;
	   from=0;
      }
      file1.close();	

      TGraph *gf1f217=new TGraph(18,xbj,F2Dp_f1f217);
      TGraph *gineft=new TGraph(18,xbj,F2Dp_ineft);
      TGraph *gf1f217_d=new TGraph(18,xbj,F2d_f1f217);
      TGraph *gf1f217_p=new TGraph(18,xbj,F2p_f1f217);
      TGraph *gineft_d=new TGraph(18,xbj,F2d_ineft);
      TGraph *gineft_p=new TGraph(18,xbj,F2p_ineft);

      TCanvas *c1=new TCanvas("c1","c1",1500,1500);
      TMultiGraph *mg=new TMultiGraph();
      gf1f217->SetMarkerStyle(8);
      gf1f217->SetMarkerColor(1);
      gf1f217->SetMarkerSize(1.5);
      gineft->SetMarkerStyle(22);
      gineft->SetMarkerColor(2);
      gineft->SetMarkerSize(1.5);

      mg->Add(gf1f217);
      mg->Add(gineft);
      mg->Draw("AP");
      mg->SetTitle(";xbj;D/p");
      mg->GetYaxis()->SetRangeUser(0.8,0.9);

      auto leg1=new TLegend(0.68,0.6,0.78,0.78);
      leg1->AddEntry(gf1f217,"f1f217","P");
      leg1->AddEntry(gineft,"Bodek","P");
      leg1->Draw();

      TCanvas *c2=new TCanvas("c2","c2",1500,1500);
      TMultiGraph *mg1=new TMultiGraph();
      gf1f217_d->SetMarkerStyle(8);
      gf1f217_d->SetMarkerColor(1);
      gf1f217_d->SetMarkerSize(1.5);
      gineft_d->SetMarkerStyle(22);
      gineft_d->SetMarkerColor(2);
      gineft_d->SetMarkerSize(1.5);

      mg1->Add(gf1f217_d);
      mg1->Add(gineft_d);
      mg1->Draw("AP");
      mg1->SetTitle("deuteron;xbj;F2d");

      auto leg2=new TLegend(0.68,0.6,0.78,0.78);
      leg2->AddEntry(gf1f217_d,"f1f217","P");
      leg2->AddEntry(gineft_d,"Bodek","P");
      leg2->Draw();

      TCanvas *c3=new TCanvas("c3","c3",1500,1500);
      TMultiGraph *mg2=new TMultiGraph();
      gf1f217_p->SetMarkerStyle(8);
      gf1f217_p->SetMarkerColor(1);
      gf1f217_p->SetMarkerSize(1.5);
      gineft_p->SetMarkerStyle(22);
      gineft_p->SetMarkerColor(2);
      gineft_p->SetMarkerSize(1.5);

      mg2->Add(gf1f217_p);
      mg2->Add(gineft_p);
      mg2->Draw("AP");
      mg2->SetTitle("proton;xbj;F2p");

      auto leg3=new TLegend(0.68,0.6,0.78,0.78);
      leg3->AddEntry(gf1f217_p,"f1f217","P");
      leg3->AddEntry(gineft_p,"Bodek","P");
      leg3->Draw();
}
