Double_t F2np(Double_t xbj=0.0){
     Double_t ratio=0.0;
     ratio=0.579363+7.49797*xbj-81.0879*xbj*xbj+395.865*pow(xbj,3)-1109.55*pow(xbj,4)+1893.08*pow(xbj,5)-1951.56*pow(xbj,6)
           +1119.91*pow(xbj,7)-274.497*pow(xbj,8);
     return ratio;
}

Double_t F2H3HE(Double_t xbj=0.0){
     Double_t ratio=0.0;
     ratio=0.982849-0.19572*xbj-2.64732*xbj*xbj+9.12714*pow(xbj,3)-2.8969*pow(xbj,4)-36.4586*pow(xbj,5)+73.4516*pow(xbj,6)
           -56.3114*pow(xbj,7)+15.7526*pow(xbj,8);
     return ratio;
}

int Xshift()
{
	ofstream outfile;
	outfile.open("Ratio_error.dat");
	outfile<<"x   F2n/F2p   F2H3/F2HE"<<endl;
	TGraph *NP=new TGraph();
	TGraph *H3HE=new TGraph();
	for(int ii=0;ii<24;ii++){
	   Double_t x1=0.15+ii*0.03;
	   Double_t yF2np1=F2np(x1);
	   Double_t x2=x1*(1.0+2.0e-4);
	   Double_t yF2np2=F2np(x2);
	   Double_t delta=(yF2np2-yF2np1)/yF2np1;
	
	   Double_t yH3HE1=F2H3HE(x1);
	   Double_t yH3HE2=F2H3HE(x2);
	   Double_t delta_H3HE=(yH3HE2-yH3HE1)/yH3HE1;

	   NP->SetPoint(ii,x1,delta);
	   H3HE->SetPoint(ii,x1,delta_H3HE);

	   outfile<<x1<<"  "<<delta<<"  "<<delta_H3HE<<endl; 
	   cout<<x1<<"  "<<x2<<"  "<<delta<<"  "<<delta_H3HE<<endl; 
	}
	outfile.close();

	TCanvas *c1=new TCanvas("c1","c1",1500,1500);
	c1->Divide(2,1);
	c1->cd(1);
        NP->SetMarkerStyle(8);
	NP->Draw("AP");

	c1->cd(2);
        H3HE->SetMarkerStyle(8);
	H3HE->Draw("AP");

	

	return 0;
}
