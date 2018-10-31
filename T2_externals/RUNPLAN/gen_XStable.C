#include <TMath.h>
#define MAXBIN 3000  
void gen_XStable(){
     Double_t xbj[MAXBIN]={0.0},Q2[MAXBIN]={0.0};
     Double_t Ep[MAXBIN]={0.0},theta[MAXBIN]={0.0};
     Double_t E0=10.5901;
     Double_t Theta0=17.5717;
     Double_t Ep0=3.1;
     Double_t M=0.938272;

     Double_t delp=4.0;   //dp/p=+/-4.0%
     Double_t delTh=2.5;  //delta_theta=+/-2.5 deg;

     Int_t MaxBinP=(2*delp/100.0*Ep0)/0.005+1+1;
     Int_t MaxBinTh=delTh/0.1*2+1;

     ofstream file;
     file.open("Carbon_kin1.inp");
     file<<"Marathon"<<endl;    
     file<<"Carbon kin1"<<endl;    
     file<<endl;    
     file<<endl;    
     file<<endl;    
     file<<"E    Ep    theta"<<endl;    
     for(int ii=0;ii<MaxBinP;ii++){
          Ep[ii]=Ep0*(1-delp/100.0)+0.005*ii;
      for(int jj=0;jj<MaxBinTh;jj++){
	  theta[jj]=Theta0-delTh+jj*0.1;
          file<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep[ii]<<" "<<theta[jj]<<endl;
      }
     }
     file.close();

}
