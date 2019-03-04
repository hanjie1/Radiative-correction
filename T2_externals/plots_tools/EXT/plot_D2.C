#include <fstream>
#include "ReadFile.h"
#define MAXBIN 8
using namespace std;

void plot_D2()
{
	Double_t x[MAXBIN],Q2[MAXBIN];
	Double_t RC[MAXBIN],RC_gas[MAXBIN],RC_pre[MAXBIN],RC_wall[MAXBIN],RC_pos[MAXBIN];
	TString Yfile;
	for(int ii=0;ii<MAXBIN;ii++){
	    x[ii]=0; Q2[ii]=0;
	    RC[ii]=0;  RC_gas[ii]=0;  RC_pre[ii]=0;
	    RC_wall[ii]=0;  RC_pos[ii]=0;
	}

        int kin=16;

        Yfile=Form("all/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,x,Q2,RC);
        Yfile=Form("gas/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,x,Q2,RC_gas);
        Yfile=Form("pre_tar/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,x,Q2,RC_pre);
        Yfile=Form("post_tar/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,x,Q2,RC_pos);
        Yfile=Form("walls/D2_kin%d_xs.out",kin);
        ReadYield(Yfile,x,Q2,RC_wall);

	TGraph *gRC=new TGraph();
	TGraph *gRC_gas=new TGraph();
	TGraph *gRC_pre=new TGraph();
	TGraph *gRC_wall=new TGraph();
	TGraph *gRC_pos=new TGraph();

	Double_t d_gas[MAXBIN],d_pre[MAXBIN],d_wall[MAXBIN],d_pos[MAXBIN];

	for(int ii=0;ii<MAXBIN;ii++){
	    if(x[ii]==0)continue;
	    d_gas[ii]=(RC_gas[ii]-RC[ii])/RC[ii];
	    d_pre[ii]=(RC_pre[ii]-RC[ii])/RC[ii];
	    d_wall[ii]=(RC_wall[ii]-RC[ii])/RC[ii];
	    d_pos[ii]=(RC_pos[ii]-RC[ii])/RC[ii];


	    gRC->SetPoint(ii,x[ii],0);
	    gRC_gas->SetPoint(ii,x[ii],d_gas[ii]);
	    gRC_pre->SetPoint(ii,x[ii],d_pre[ii]);
	    gRC_wall->SetPoint(ii,x[ii],d_wall[ii]);
	    gRC_pos->SetPoint(ii,x[ii],d_pos[ii]);
	}

	ofstream outfile;
	outfile.open("OUT/D2_RC.csv",ios::app);
	outfile<<"kin"<<kin<<endl;
	outfile<<"x,Q2,RC(all),RC(gas),RC(pre-target),RC(wall),RC(post-target)"<<endl;
	for(int ii=0;ii<MAXBIN;ii++){
	    if(x[ii]==0)continue;
	    outfile<<fixed;
	    outfile<<setprecision(3)<<x[ii]<<","<<Q2[ii]<<",";
	    outfile<<setprecision(3)<<RC[ii]<<","<<RC_gas[ii]<<","<<RC_pre[ii]<<","<<RC_wall[ii]<<","<<RC_pos[ii]<<endl;;
	}
	outfile.close();



	TCanvas *c1=new TCanvas("c1");
	TMultiGraph *mg=new TMultiGraph();
	int color[5]={2,1,3,4,5};
        gRC->SetMarkerStyle(8);
        gRC->SetMarkerColor(color[0]);
        gRC_gas->SetMarkerStyle(8);
        gRC_gas->SetMarkerColor(color[1]);
        gRC_pre->SetMarkerStyle(8);
        gRC_pre->SetMarkerColor(color[2]);
        gRC_wall->SetMarkerStyle(8);
        gRC_wall->SetMarkerColor(color[3]);
        gRC_pos->SetMarkerStyle(8);
        gRC_pos->SetMarkerColor(color[4]);
	mg->Add(gRC);
	mg->Add(gRC_gas);	
	mg->Add(gRC_pre);	
	mg->Add(gRC_wall);	
	mg->Add(gRC_pos);	
	mg->Draw("AP");	

        auto leg1=new TLegend(0.7,0.6,0.811,0.811);
        leg1->AddEntry(gRC,"All","P");
        leg1->AddEntry(gRC_gas,"gas","P");
        leg1->AddEntry(gRC_pre,"pre","P");
        leg1->AddEntry(gRC_wall,"wall","P");
        leg1->AddEntry(gRC_pos,"pos","P");
        leg1->Draw();

}

