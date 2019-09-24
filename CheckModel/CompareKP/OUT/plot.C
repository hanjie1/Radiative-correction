void plot(){
    Double_t x[1100]={0.0},F2p[1100]={0.0},F2p_err[1100]={0.0},F2d[1100]={0.0},F2d_err[1100]={0.0};
    Double_t F2p_KP[1100]={0.0},F2d_KP[1100]={0.0};
    ifstream file;
    file.open("F2dp_Whitlow.out");
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;

    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          F2p[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2p_err[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2d[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2d_err[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2p_KP[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          F2d_KP[nn-1]=atof(content.Data());

          from=0;
          nn++;
     }
    file.close();

    TGraphErrors *gDiff_p=new TGraphErrors();
    TGraphErrors *gDiff_d=new TGraphErrors();

    for(int ii=0;ii<nn-1;ii++){
	Double_t tmpDiff=(F2p_KP[ii]-F2p[ii])/F2p_KP[ii];
	Double_t tmperr=F2p_err[ii]/F2p_KP[ii];
	gDiff_p->SetPoint(ii,x[ii],tmpDiff);
	gDiff_p->SetPointError(ii,0,tmperr);
	tmpDiff=(F2d_KP[ii]-F2d[ii])/F2d_KP[ii];
	tmperr=F2d_err[ii]/F2d_KP[ii];
	gDiff_d->SetPoint(ii,x[ii],tmpDiff);
	gDiff_d->SetPointError(ii,0,tmperr);
    }

    TCanvas *c1=new TCanvas("c1","c1",1000,1000);
    c1->Divide(1,2);
    c1->cd(1);
    gDiff_p->SetMarkerStyle(8);
    gDiff_p->SetFillStyle(3001);
    gDiff_p->SetFillColor(kCyan-3);
    gDiff_p->SetTitle(";Bjorken x;(F_{2}^{p}(KP)-F_{2}^{p}(Whitlow))/F_{2}^{p}(KP)");
    gDiff_p->Draw("APE3");
 
    c1->cd(2);
    gDiff_d->SetMarkerStyle(8);
    gDiff_d->SetFillStyle(3001);
    gDiff_d->SetFillColor(kCyan-3);
    gDiff_d->SetTitle(";Bjorken x;(F_{2}^{d}(KP)-F_{2}^{d}(Whitlow))/F_{2}^{d}(KP)");
    gDiff_d->Draw("APE3");
    c1->Print("compare.pdf");
}
