void fitplot()
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

	  Double_t H3D[101]={0.0},HED[101]={0.0};
	  Double_t F2np[101]={0.0},H3HE[101]={0.0};
	  Double_t REMCD[101]={0.0},superR[101]={0.0};
	  Double_t REMCHe[101]={0.0},REMCH3[101]={0.0};
	  Double_t HEDISO[101]={0.0},H3DISO[101]={0.0};
	  for(int ii=0;ii<nn;ii++){
		H3D[ii]=F2H3[ii]/F2D[ii];
		HED[ii]=F2HE[ii]/F2D[ii];
		F2np[ii]=F2n[ii]/F2p[ii];
		H3HE[ii]=F2H3[ii]/F2HE[ii];
		REMCD[ii]=F2D[ii]*2.0/(F2n[ii]+F2p[ii]);
		REMCHe[ii]=F2HE[ii]*3.0/(F2n[ii]+2.0*F2p[ii]);
		REMCH3[ii]=F2H3[ii]*3.0/(2.0*F2n[ii]+F2p[ii]);
		superR[ii]=(F2HE[ii]/(2.0*F2p[ii]+F2n[ii]))/(F2H3[ii]/(F2p[ii]+2.0*F2n[ii]));

		HEDISO[ii]=F2HE[ii]*3.0/(F2D[ii]*2.0)*(F2n[ii]+F2p[ii])/(2.0*F2p[ii]+F2n[ii]);
		H3DISO[ii]=F2H3[ii]*3.0/(F2D[ii]*2.0)*(F2n[ii]+F2p[ii])/(F2p[ii]+2.0*F2n[ii]);
	  }

	  TGraph *gH3=new TGraph(nn,xbj,H3D);
	  TGraph *gHE=new TGraph(nn,xbj,HED);
	  TGraph *gH3ISO=new TGraph(nn,xbj,H3DISO);
	  TGraph *gHEISO=new TGraph(nn,xbj,HEDISO);
	  TGraph *gH3HE=new TGraph(nn,xbj,H3HE);
	  TGraph *gNP=new TGraph(nn,xbj,F2np);
	  TGraph *gREMCD=new TGraph(nn,xbj,REMCD);
	  TGraph *gREMCH3=new TGraph(nn,xbj,REMCH3);
	  TGraph *gREMCHe=new TGraph(nn,xbj,REMCHe);
	  TGraph *gsuperR=new TGraph(nn,xbj,superR);

	  TCanvas *c1=new TCanvas("c1");
 	  gH3->Draw("AP*");
	  gH3->SetTitle("H3/D2");

	  TCanvas *c2=new TCanvas("c2");
 	  gHE->Draw("AP*");
	  gHE->SetTitle("HE/D2");

	  TCanvas *c3=new TCanvas("c3");
 	  gH3HE->Draw("AP*");
	  gH3HE->SetTitle("H3/HE");

	  TCanvas *c4=new TCanvas("c4");
 	  gNP->Draw("AP*");
	  gNP->SetTitle("F2n/F2p");

	  TCanvas *c5=new TCanvas("c5");
 	  gREMCD->Draw("AP*");
	  gREMCD->SetTitle("F2D/(F2n+F2p)");

	  TCanvas *c6=new TCanvas("c6");
 	  gsuperR->Draw("AP*");
	  gsuperR->SetTitle("super ratio");

	  TCanvas *c7=new TCanvas("c7");
 	  gREMCH3->Draw("AP*");
	  gREMCH3->SetTitle("F2H3/(2*F2n+F2p)");

	  TCanvas *c8=new TCanvas("c8");
 	  gREMCHe->Draw("AP*");
	  gREMCHe->SetTitle("F2He/(F2n+2*F2p)");

	  TCanvas *c9=new TCanvas("c9");
 	  gHEISO->Draw("AP*");
	  gHEISO->SetTitle("He3 EMC isoscalar");

	  TCanvas *c10=new TCanvas("c10");
 	  gH3ISO->Draw("AP*");
	  gH3ISO->SetTitle("H3 EMC isoscalar");
}
