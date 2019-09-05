#include <TMath.h>
#define MAXBIN 70
void genPlot(){
     Double_t xx[MAXBIN]={0.0};
     Double_t Q2[MAXBIN]={0.0};

     for(int ii=0;ii<MAXBIN;ii++){
	xx[ii]=0.15+0.01*ii;
        Q2[ii]=14.0*xx[ii];
     }

     Double_t Ep[MAXBIN]={0.0},theta[MAXBIN]={0.0};
     Double_t E0=10.59;
     Double_t M=0.938272;

     ofstream file;
     TString filename="ForPlot.inp";
     file.open(filename);
     file<<"Marathon"<<endl;
     file<<filename.Data()<<endl;
     file<<endl;
     file<<endl;
     file<<endl;
     file<<"E    Ep    theta"<<endl;
     for(int ii=0;ii<MAXBIN;ii++){
         Ep[ii]=E0-Q2[ii]/(2.0*M*xx[ii]);
         theta[ii]=asin(sqrt(Q2[ii]/(4.0*E0*Ep[ii])))*2.0*180.0/TMath::Pi();
         file<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep[ii]<<" "<<theta[ii]<<endl;
     }
     file.close();



}
