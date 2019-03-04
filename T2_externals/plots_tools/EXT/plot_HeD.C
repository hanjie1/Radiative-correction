#include <fstream>
#include "ReadFile.h"
#define MAXBIN 8
using namespace std;

void plot_HeD(){
	Double_t He3_x[MAXBIN],He3_Q2[MAXBIN],D2_x[MAXBIN],D2_Q2[MAXBIN];
        Double_t He3_RC[MAXBIN],He3_RC_gas[MAXBIN],He3_RC_pre[MAXBIN],He3_RC_wall[MAXBIN],He3_RC_pos[MAXBIN];
        Double_t D2_RC[MAXBIN],D2_RC_gas[MAXBIN],D2_RC_pre[MAXBIN],D2_RC_wall[MAXBIN],D2_RC_pos[MAXBIN];

	for(int ii=0;ii<MAXBIN;ii++){
	    He3_x[ii]=0; He3_Q2[ii]=0; D2_x[ii]=0;  D2_Q2[ii]=0;
            He3_RC[ii]=0;  He3_RC_gas[ii]=0;  He3_RC_pre[ii]=0;
            He3_RC_wall[ii]=0;  He3_RC_pos[ii]=0;
            D2_RC[ii]=0;  D2_RC_gas[ii]=0;  D2_RC_pre[ii]=0;
            D2_RC_wall[ii]=0;  D2_RC_pos[ii]=0;
	}	

        int kin=16;
	TString Yfile;

        Yfile=Form("all/He3_kin%d_xs.out",kin);
        ReadYield(Yfile,He3_x,He3_Q2,He3_RC);
        Yfile=Form("gas/He3_kin%d_xs.out",kin);
        ReadYield(Yfile,He3_x,He3_Q2,He3_RC_gas);
        Yfile=Form("pre_tar/He3_kin%d_xs.out",kin);
        ReadYield(Yfile,He3_x,He3_Q2,He3_RC_pre);
        Yfile=Form("post_tar/He3_kin%d_xs.out",kin);
        ReadYield(Yfile,He3_x,He3_Q2,He3_RC_pos);
        Yfile=Form("walls/He3_kin%d_xs.out",kin);
        ReadYield(Yfile,He3_x,He3_Q2,He3_RC_wall);

        Yfile=Form("all/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,D2_x,D2_Q2,D2_RC);
        Yfile=Form("gas/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,D2_x,D2_Q2,D2_RC_gas);
        Yfile=Form("pre_tar/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,D2_x,D2_Q2,D2_RC_pre);
        Yfile=Form("post_tar/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,D2_x,D2_Q2,D2_RC_pos);
        Yfile=Form("walls/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,D2_x,D2_Q2,D2_RC_wall);

        ofstream outfile;
        outfile.open("OUT/HeD_RC.csv",ios::app);
        outfile<<"kin"<<kin<<endl;
        outfile<<"x,Q2,RC(all) ratio,RC(gas) ratio,RC(pre-target) ratio,RC(wall) ratio,RC(post-target) ratio"<<endl;
        for(int ii=0;ii<MAXBIN;ii++){
            if(He3_x[ii]==0||D2_x[ii]==0)continue;
	    if(abs(He3_x[ii]-D2_x[ii])>0.01)continue;
	    Double_t tmpr_RC=He3_RC[ii]/D2_RC[ii];
	    Double_t tmpr_RC_gas=He3_RC_gas[ii]/D2_RC_gas[ii];
	    Double_t tmpr_RC_pre=He3_RC_pre[ii]/D2_RC_pre[ii];
	    Double_t tmpr_RC_wall=He3_RC_wall[ii]/D2_RC_wall[ii];
	    Double_t tmpr_RC_pos=He3_RC_pos[ii]/D2_RC_pos[ii];
            outfile<<fixed;
            outfile<<setprecision(3)<<D2_x[ii]<<","<<D2_Q2[ii]<<",";
            outfile<<setprecision(3)<<tmpr_RC<<","<<tmpr_RC_gas<<","<<tmpr_RC_pre<<","<<tmpr_RC_wall<<","<<tmpr_RC_pos<<endl;;
        }
        outfile.close();


}
