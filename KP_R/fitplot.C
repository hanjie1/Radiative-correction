void fitplot()
{
          TString filename="F2dis_os1tm1ht1mec1_Dav18_He3Salme";
	  ifstream file1;
          file1.open(filename);

          Double_t xbj[101]={0.0},Q2[101]={0.0};
	  Double_t F2D[101]={0.0},F2H3[101]={0.0},F2HE[101]={0.0};
          Ssiz_t from=0;
          TString content,tmp;
          int nn=0;
          while(tmp.ReadLine(file1)){
                tmp.Tokenize(content,from,"  ");
                xbj[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                Q2[nn]=atof(content.Data());
                tmp.Tokenize(content,from,"  ");
                tmp.Tokenize(content,from,"  ");
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
	  for(int ii=0;ii<nn;ii++){
		H3D[ii]=F2H3[ii]/F2D[ii];
		HED[ii]=F2HE[ii]/F2D[ii];
	  }

	  TGraph *gH3=new TGraph(nn,xbj,H3D);
	  TGraph *gHE=new TGraph(nn,xbj,HED);

	  TCanvas *c1=new TCanvas("c1");
 	  gH3->Draw("AP*");
	  gH3->SetTitle("H3/D2");

	  TCanvas *c2=new TCanvas("c2");
 	  gHE->Draw("AP*");
	  gHE->SetTitle("HE/D2");
}
