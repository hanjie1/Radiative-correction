void R_plot()
{
          TString filename="F2dis_os1tm1ht1mec1_Dav18_He3Salme";
	  ifstream file1;
          file1.open(filename);

          Double_t xbj[101]={0.0},Q2[101]={0.0};
	  Double_t F2D[101]={0.0},F2H3[101]={0.0},F2HE[101]={0.0};
	  Double_t F2n[101]={0.0},F2p[101]={0.0};
          Ssiz_t from=0;
          TString content,tmp;
          int nn=0;
          while(tmp.ReadLine(file1)){
                tmp.Tokenize(content,from,"  ");
                xbj[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                Q2[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                F2p[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                F2n[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                F2D[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                F2H3[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                F2HE[nn]=atof(content.Data());
                nn++;
                from=0;
              }
          file1.close();

	  Double_t RH3[101]={0.0},RHe[101]={0.0},RD[101]={0.0};
	  Double_t RH3D[101]={0.0},RHeD[101]={0.0},RH3He[101]={0.0};

	  TGraph *gRH3=new TGraph();
	  TGraph *gRHE=new TGraph();
	  TGraph *gRD=new TGraph();
	  TGraph *gH3D=new TGraph();
	  TGraph *gHeD=new TGraph();
	  TGraph *gH3He=new TGraph();
	  for(int ii=0;ii<nn;ii++){
		if(xbj[ii]>0.85)continue;
		RH3[ii]=F2H3[ii]*3.0/(F2p[ii]+2.0*F2n[ii]);
		RHe[ii]=F2HE[ii]*3.0/(2.0*F2p[ii]+F2n[ii]);
		RD[ii]=F2D[ii]*2.0/(F2n[ii]+F2p[ii]);

		RH3D[ii]=RH3[ii]/RD[ii];
		RHeD[ii]=RHe[ii]/RD[ii];
		RH3He[ii]=RH3[ii]/RHe[ii];

		gRH3->SetPoint(ii,xbj[ii],RH3[ii]);
		gRHE->SetPoint(ii,xbj[ii],RHe[ii]);
		gRD->SetPoint(ii,xbj[ii],RD[ii]);
		gH3D->SetPoint(ii,xbj[ii],RH3D[ii]);
		gHeD->SetPoint(ii,xbj[ii],RHeD[ii]);
		gH3He->SetPoint(ii,xbj[ii],RH3He[ii]);
	  }

	  TCanvas *c1=new TCanvas("c1","c1",1500,1000);
	  gRH3->SetLineStyle(1);
	  gRH3->SetLineColor(2);
	  gRH3->SetLineWidth(3);

	  gRHE->SetLineStyle(5);
	  gRHE->SetLineColor(4);
	  gRHE->SetLineWidth(3);

	  gRD->SetLineStyle(2);
	  gRD->SetLineColor(1);
	  gRD->SetLineWidth(3);

	  TMultiGraph *mg=new TMultiGraph();
	  mg->Add(gRH3);
	  mg->Add(gRHE);
	  mg->Add(gRD);
 	  mg->Draw("AL");
	  mg->SetTitle(";Bjorken x;EMC-type ratio;");

          auto leg1=new TLegend(0.2,0.6,0.35,0.85);
          leg1->AddEntry(gRD,"#scale[0.8]{R_{21}}","L");
   	  leg1->AddEntry(gRHE,"#scale[0.8]{R_{32}}","L");
   	  leg1->AddEntry(gRH3,"#scale[0.8]{R_{31}}","L");
	  leg1->SetMargin(0.5);
  	  leg1->Draw();

   	  c1->Print("R_EMC.pdf");

	  TCanvas *c2=new TCanvas("c2","c2",1500,1000);
	  gHeD->SetLineStyle(2);
	  gHeD->SetLineColor(8);
	  gHeD->SetLineWidth(3);

	  gH3D->SetLineStyle(5);
	  gH3D->SetLineColor(4);
	  gH3D->SetLineWidth(3);

	  gH3He->SetLineStyle(1);
	  gH3He->SetLineColor(2);
	  gH3He->SetLineWidth(3);

	  TMultiGraph *mg1=new TMultiGraph();
	  mg1->Add(gHeD);
	  mg1->Add(gH3D);
	  mg1->Add(gH3He);
 	  mg1->Draw("AL");
	  mg1->SetTitle(";Bjorken x;ratio;");

          auto leg2=new TLegend(0.6,0.6,0.75,0.85);
   	  leg2->AddEntry(gHeD,"#scale[0.7]{R_{32}/R_{21}}","L");
   	  leg2->AddEntry(gH3D,"#scale[0.7]{R_{31}/R_{21}}","L");
   	  leg2->AddEntry(gH3He,"#scale[0.7]{R_{31}/R_{32}}","L");
	  leg2->SetMargin(0.5);
  	  leg2->Draw();

   	  c2->Print("ratioR.pdf");

}
