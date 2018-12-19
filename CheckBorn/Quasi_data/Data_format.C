#include <iomanip>
#define MAXNUM 2135
using namespace std;
void Data_format()
{
    ifstream infile;
    infile.open("2H.dat");

    Double_t E0[MAXNUM]={0.0},Nu[MAXNUM]={0.0},Theta[MAXNUM]={0.0},XS[MAXNUM]={0.0},XS_err[MAXNUM]={0.0};
    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    while(tmp.ReadLine(infile)){
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
    infile.close();

    ofstream outfile;
    outfile.open("D2_QSI.dat");
    outfile<<"E0    Ep     Theta     x      Q2     w2     sigma     sigma_err"<<endl;
    for(int ii=0;ii<MAXNUM;ii++){
	Double_t Ep=E0[ii]-Nu[ii];
	Double_t Q2=4.0*E0[ii]*Ep*sin(Theta[ii]*TMath::Pi()/2.0/180.0)*sin(Theta[ii]*TMath::Pi()/2.0/180.0);
	Double_t x=Q2/(2.0*0.938272*Nu[ii]);
	Double_t W2=0.938272*0.938272+2*0.938272*Nu[ii]-Q2;
	if(x<1.0){
	  outfile<<right<<fixed;
	  outfile<<setprecision(4)<<setw(9)<<E0[ii]<<setw(9)<<Ep<<setw(9)<<setprecision(2)<<Theta[ii]
                 <<setw(9)<<setprecision(3)<<x<<setw(9)<<Q2<<setw(9)<<W2<<setw(15)<<scientific<<XS[ii]<<setw(15)<<XS_err[ii]<<endl; 
	}
    } 
    outfile.close();


}
