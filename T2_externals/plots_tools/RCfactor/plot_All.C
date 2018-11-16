#include <fstream>
#include "ReadFile.h"
using namespace std;

void plot_All()
{
     Double_t H1_x[5][MAXBIN],H1_Q2[5][MAXBIN],H1_Born[5][MAXBIN],H1_Rad[5][MAXBIN],H1_Born1[5][MAXBIN],H1_Rad1[5][MAXBIN];
     Double_t D2_x[11][MAXBIN],D2_Q2[11][MAXBIN],D2_Born[11][MAXBIN],D2_Rad[11][MAXBIN],D2_Born1[11][MAXBIN],D2_Rad1[11][MAXBIN];
     Double_t He_x[11][MAXBIN],He_Q2[11][MAXBIN],He_Born[11][MAXBIN],He_Rad[11][MAXBIN],He_Born1[11][MAXBIN],He_Rad1[11][MAXBIN];
     Double_t H3_x[11][MAXBIN],H3_Q2[11][MAXBIN],H3_Born[11][MAXBIN],H3_Rad[11][MAXBIN],H3_Born1[11][MAXBIN],H3_Rad1[11][MAXBIN];

     for(int ii=0;ii<5;ii++){
	 for(int jj=0;jj<MAXBIN;jj++){
             if(ii==0&&jj<5)
                H1_x[ii][jj]=0.0; H1_Q2[ii][jj]=0.0; H1_Born[ii][jj]=0.0; H1_Rad[ii][jj]=0.0;H1_Born1[ii][jj]=0.0; H1_Rad1[ii][jj]=0.0;
             D2_x[ii][jj]=0.0; D2_Q2[ii][jj]=0.0; D2_Born[ii][jj]=0.0; D2_Rad[ii][jj]=0.0;D2_Born1[ii][jj]=0.0; D2_Rad1[ii][jj]=0.0;
             He_x[ii][jj]=0.0; He_Q2[ii][jj]=0.0; He_Born[ii][jj]=0.0; He_Rad[ii][jj]=0.0;He_Born1[ii][jj]=0.0; He_Rad1[ii][jj]=0.0;
             H3_x[ii][jj]=0.0; H3_Q2[ii][jj]=0.0; H3_Born[ii][jj]=0.0; H3_Rad[ii][jj]=0.0;H3_Born1[ii][jj]=0.0; H3_Rad1[ii][jj]=0.0;
     }}

   TString Yfile;
   int kin[11]={0,1,2,3,4,5,7,9,11,13,15};
   for(int ii=0;ii<11;ii++){
       if(ii<5){
          Yfile=Form("Bodek_final/H1_kin%d_xs.out",kin[ii]);
          ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born,H1_Rad); 
          Yfile=Form("f1f217/H1_kin%d_xs.out",kin[ii]);
          ReadYield(Yfile,kin[ii],H1_x,H1_Q2,H1_Born1,H1_Rad1); 
       }
       Yfile=Form("Bodek_final/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born,D2_Rad); 
       Yfile=Form("f1f217/D2_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],D2_x,D2_Q2,D2_Born1,D2_Rad1); 
       Yfile=Form("Bodek_final/He_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x,He_Q2,He_Born,He_Rad); 
       Yfile=Form("f1f217/He_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],He_x,He_Q2,He_Born1,He_Rad1); 
       Yfile=Form("Bodek_final/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born,H3_Rad); 
       Yfile=Form("f1f217/H3_kin%d_xs.out",kin[ii]);
       ReadYield(Yfile,kin[ii],H3_x,H3_Q2,H3_Born1,H3_Rad1); 
   }

   

}
