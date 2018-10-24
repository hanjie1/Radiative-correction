#define MAXNUM 3000
#define PI TMath::Pi()
#define Mp 0.938272
void plotData()
{
    int Z,A;
    cout<<"Input Z: ";
    cin>>Z;
    cout<<"Input A: ";
    cin>>A;

    TString filename;
    if(Z==1&&A==2)filename="2H.dat"; 
    if(Z==2&&A==3)filename="3He.dat"; 
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
    
    ofstream outfile;
    outfile.open("4Hei_Q2_xbj.dat");

    ofstream outfile1;
    outfile1.open("4He_50.dat");
    ofstream outfile2;
    outfile2.open("data/4He_50.dat");

    nn=nn-1;
    Double_t Q2[MAXNUM]={0.0},xbj[MAXNUM]={0.0},W2[MAXNUM]={0.0};
    Double_t Eprime[MAXNUM]={0.0};
    for(int ii=0;ii<nn;ii++){
	Q2[ii]=4.0*E0[ii]*(E0[ii]-Nu[ii])*sin(Theta[ii]/2.0*PI/500.0)*sin(Theta[ii]/2.0*PI/500.0);
        xbj[ii]=Q2[ii]/(2.0*Mp*Nu[ii]);
	W2[ii]=Mp*Mp+2*Mp*Nu[ii]-Q2[ii];
        Eprime[ii]=E0[ii]-Nu[ii];
        outfile<<E0[ii]<<"  "<<Theta[ii]<<"  "<<Q2[ii]<<"  "<<xbj[ii]<<endl; 
        if(E0[ii]==5.766&&Theta[ii]==50)outfile1<<left<<fixed<<setprecision(3)<<setw(6)<<E0[ii]<<" "<<setw(6)<<Eprime[ii]<<" "<<setw(6)<<Theta[ii]<<endl;
        if(E0[ii]==5.766&&Theta[ii]==50)outfile2<<E0[ii]<<" "<<Eprime[ii]<<" "<<Theta[ii]<<" "<<XS[ii]<<" "<<XS_err[ii]<<endl;
    }
    outfile.close();

    outfile1.close();
    outfile2.close();

    TCanvas *c1=new TCanvas("c1");
    TGraphErrors *gXbj=new TGraphErrors();
    int npoint=0;
    for(int ii=0;ii<nn;ii++){
	if(E0[ii]!=5.766)continue;
	if(Theta[ii]!=50)continue;
	gXbj->SetPoint(npoint,xbj[ii],Q2[ii]);
        gXbj->SetPointError(npoint,0,0);
        npoint++;
    }
    gXbj->SetMarkerStyle(8);
    gXbj->SetMarkerColor(4);
    gXbj->Draw("AP");
    gXbj->SetTitle("Q2 vs. xbj(E0=5.766,Theta=50);xbj;Q2");

    TCanvas *c2=new TCanvas("c2");
    TGraphErrors *gNu=new TGraphErrors();
    npoint=0;
    for(int ii=0;ii<nn;ii++){
	if(E0[ii]!=5.766)continue;
	if(Theta[ii]!=50)continue;
	gNu->SetPoint(npoint,Nu[ii],XS[ii]);
        gNu->SetPointError(npoint,0,XS_err[ii]);
        npoint++;
    }
    gNu->SetMarkerStyle(8);
    gNu->SetMarkerColor(4);
    gNu->Draw("AP");
    gNu->SetTitle("XS vs. nu;nu;XS");

    TCanvas *c3=new TCanvas("c3");
    TGraph *gQ2E0=new TGraph(nn,xbj,Q2);
    gQ2E0->SetMarkerStyle(8);
    gQ2E0->Draw("AP");
    gQ2E0->SetTitle("Q2 vs. xbj;xbj;Q2");

}
