#define MAXNUM 3000
#define PI TMath::Pi()
#define Mp 0.938272
void plotData_Color()
{
    int Z,A;
    cout<<"Input Z: ";
    cin>>Z;
    cout<<"Input A: ";
    cin>>A;

    TString filename;
    if(Z==1&&A==2)filename="2H.dat"; 
    if(Z==2&&A==3)filename="4He.dat"; 
    if(Z==2&&A==4)filename="4He.dat"; 

    ifstream infile1;
    infile1.open(Form("data/%s",filename.Data()));

    Double_t E0[MAXNUM]={0.0},Nu[MAXNUM]={0.0},Theta[MAXNUM]={0.0},XS[MAXNUM]={0.0},XS_err[MAXNUM]={0.0};
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    while(tmp.ReadLine(infile1)){
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          E0[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Theta[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Nu[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          XS[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          XS_err[nn]=atof(content.Data());
          nn++;
          from=0;
    }
    infile1.close();
    
    nn=nn-1;
    Double_t Q2[MAXNUM]={0.0},xbj[MAXNUM]={0.0},W2[MAXNUM]={0.0};
    Double_t Eprime[MAXNUM]={0.0};
    for(int ii=0;ii<nn;ii++){
	Q2[ii]=4.0*E0[ii]*(E0[ii]-Nu[ii])*sin(Theta[ii]/2.0*PI/180.0)*sin(Theta[ii]/2.0*PI/180.0);
        xbj[ii]=Q2[ii]/(2.0*Mp*Nu[ii]);
	W2[ii]=Mp*Mp+2*Mp*Nu[ii]-Q2[ii];
        Eprime[ii]=E0[ii]-Nu[ii];
    }
    TCanvas *c1=new TCanvas("c1");
    TGraphErrors *gXbj=new TGraphErrors();
    int npoint=0;
    for(int ii=0;ii<nn;ii++){
	if(E0[ii]!=2.020)continue;
	if(Theta[ii]!=15.022)continue;
	gXbj->SetPoint(npoint,xbj[ii],Q2[ii]);
        gXbj->SetPointError(npoint,0,0);
        npoint++;
    }
    gXbj->SetMarkerStyle(8);
    gXbj->SetMarkerColor(4);
    gXbj->Draw("AP");
    gXbj->SetTitle("Q2 vs. xbj(E0=2.020,Theta=15.022);xbj;Q2");

    TCanvas *c2=new TCanvas("c2");
    TGraphErrors *gNu=new TGraphErrors();
    npoint=0;
    for(int ii=0;ii<nn;ii++){
	if(E0[ii]!=2.020)continue;
	if(Theta[ii]!=15.022)continue;
	gNu->SetPoint(npoint,Nu[ii],XS[ii]);
        gNu->SetPointError(npoint,0,XS_err[ii]);
        npoint++;
    }
    gNu->SetMarkerStyle(8);
    gNu->SetMarkerColor(4);
    gNu->Draw("AP");
    gNu->SetTitle("XS vs. nu;nu;XS");

    TCanvas *c3=new TCanvas("c3");
    TGraphErrors *gXS=new TGraphErrors();
    npoint=0;
    for(int ii=0;ii<nn;ii++){
	if(E0[ii]!=2.020)continue;
	if(Theta[ii]!=15.022)continue;
	gXS->SetPoint(npoint,xbj[ii],XS[ii]);
        gXS->SetPointError(npoint,0,XS_err[ii]);
        npoint++;
    }
    gXS->SetMarkerStyle(8);
    gXS->SetMarkerColor(4);
    gXS->Draw("AP");
    gXS->SetTitle("XS vs. xbj;xbj;XS");





}
