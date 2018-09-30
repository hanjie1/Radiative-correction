#include <TMath.h>
#define MAXBIN 35  
void All_genEtheta(){
   TString target[4]={"H1","D2","He3","H3"};
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   Double_t Ep_central=3.1;
   Double_t Th_central[11]={16.8075,17.5718,19.1125,20.5751,21.9402,23.2064,25.5858,27.7642,29.8089,31.7275,33.5551};
   const int nBin[11]={4,4,4,5,5,5,6,6,8,7,6};
   const double xmin[11]={0.159,0.176,0.210,0.247,0.283,0.318,0.391,0.462,0.535,0.606,0.678};
   const double xmax[11]={0.237,0.258,0.300,0.343,0.385,0.427,0.511,0.591,0.672,0.753,0.828};


   for(int jj=0;jj<4;jj++){
    for(int kk=0;kk<11;kk++){
     Double_t xbj[MAXBIN]={0.0};
     Double_t dBin=(xmax[kk]-xmin[kk])/(nBin[kk]*1.0);
     for(int ii=0;ii<nBin[kk];ii++){
         xbj[ii]=xmin[kk]+ii*dBin+dBin/2.0;
cout<<xbj[ii]<<", ";
         //xbj[ii]=0.25+ii*0.01;
     }
     cout<<endl;

     TString filename=Form("%s_kin%d.inp",target[jj].Data(),kin[kk]);

     Double_t Ep[MAXBIN]={0.0},theta[MAXBIN]={0.0};
     Double_t E0=10.59;
     Double_t M=0.938272;
     Double_t Q2_central=4.0*E0*Ep_central*pow(sin(Th_central[kk]*TMath::Pi()/(2.0*180)),2);
cout<<Q2_central<<endl;
     ofstream file;
     file.open(filename);
     file<<"Marathon"<<endl;    
     file<<filename.Data()<<endl;    
     file<<endl;    
     file<<endl;    
     file<<endl;    
     file<<"E    Ep    theta"<<endl;    
     for(int ii=0;ii<MAXBIN;ii++){
         if(xbj[ii]==0.0)continue;
         Ep[ii]=E0-Q2_central/(2.0*M*xbj[ii]);
         theta[ii]=asin(sqrt(Q2_central/(4.0*E0*Ep[ii])))*2*180/TMath::Pi();
         file<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep[ii]<<" "<<theta[ii]<<endl;
     }
     file.close();
    }
   }

}
