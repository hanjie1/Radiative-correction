#define MAXBIN 610
void plot_D2()
{
     TString filename;
     ifstream file1;
     cout<<"Input file name: ";
     cin>>filename;
     file1.open(Form("../OUT/%s",filename.Data()));

     Double_t xbj[MAXBIN]={0.0},Q2[MAXBIN]={0.0},W2[MAXBIN]={0.0};
     Double_t F2D[MAXBIN]={0.0},F2D_stat[MAXBIN]={0.0},F2D_sys[MAXBIN]={0.0};
     Double_t F2E[MAXBIN]={0.0},F2B[MAXBIN]={0.0};
     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           xbj[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           Q2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           W2[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D[nn-1]=atof(content.Data())/1000.0;
           tmp.Tokenize(content,from," ");
           F2D_stat[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_sys[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2E[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2B[nn-1]=atof(content.Data());
           nn++;
           from=0;
         }
     file1.close();

     TGraphErrors *gData[4];
     TGraphErrors *gE[4];
     TGraphErrors *gB[4];
     TGraphErrors *gDE[4];
     TGraphErrors *gDB[4];

     for(int ii=0;ii<4;ii++){
	 gData[ii]=new TGraphErrors();
	 gE[ii]=new TGraphErrors();
	 gB[ii]=new TGraphErrors();
	 gDE[ii]=new TGraphErrors();
	 gDB[ii]=new TGraphErrors();
     }

     ofstream myfile;
     myfile.open("Bad_D2_whitlow.out");
     myfile<<"x     Q2     W2    F2_data    F2_data_err    F2_f1f217    F2_ineft"<<endl;
     int np=0;
     for(int ii=0;ii<MAXBIN;ii++){
         Double_t tmp_err=sqrt(F2D_stat[ii]*F2D_stat[ii]+F2D_sys[ii]*F2D_sys[ii]);
	 Double_t tmp_r=F2D[ii]/F2E[ii];
         Double_t tmp_r1=F2D[ii]/F2B[ii];
         if(Q2[ii]<8){
            gData[0]->SetPoint(np,xbj[ii],F2D[ii]);
            gData[0]->SetPointError(np,0,tmp_err);
            gE[0]->SetPoint(np,xbj[ii],F2E[ii]);
            gB[0]->SetPoint(np,xbj[ii],F2B[ii]);
            gDE[0]->SetPoint(np,xbj[ii],tmp_r);
            gDE[0]->SetPointError(np,0,tmp_err/F2E[ii]);
            gDB[0]->SetPoint(np,xbj[ii],tmp_r1);
            gDB[0]->SetPointError(np,0,tmp_err/F2B[ii]);
         }
         if(Q2[ii]<10 && Q2[ii]>=8){
            gData[1]->SetPoint(np,xbj[ii],F2D[ii]);
            gData[1]->SetPointError(np,0,tmp_err);
            gE[1]->SetPoint(np,xbj[ii],F2E[ii]);
            gB[1]->SetPoint(np,xbj[ii],F2B[ii]);
            gDE[1]->SetPoint(np,xbj[ii],tmp_r);
            gDE[1]->SetPointError(np,0,tmp_err/F2E[ii]);
            gDB[1]->SetPoint(np,xbj[ii],tmp_r1);
            gDB[1]->SetPointError(np,0,tmp_err/F2B[ii]);
         }
         if(Q2[ii]<12 && Q2[ii]>=10){
            gData[2]->SetPoint(np,xbj[ii],F2D[ii]);
            gData[2]->SetPointError(np,0,tmp_err);
            gE[2]->SetPoint(np,xbj[ii],F2E[ii]);
            gB[2]->SetPoint(np,xbj[ii],F2B[ii]);
            gDE[2]->SetPoint(np,xbj[ii],tmp_r);
            gDE[2]->SetPointError(np,0,tmp_err/F2E[ii]);
            gDB[2]->SetPoint(np,xbj[ii],tmp_r1);
            gDB[2]->SetPointError(np,0,tmp_err/F2B[ii]);
         }
         if(Q2[ii]<14 && Q2[ii]>=12){
            gData[3]->SetPoint(np,xbj[ii],F2D[ii]);
            gData[3]->SetPointError(np,0,tmp_err);
            gE[3]->SetPoint(np,xbj[ii],F2E[ii]);
            gB[3]->SetPoint(np,xbj[ii],F2B[ii]);
            gDE[3]->SetPoint(np,xbj[ii],tmp_r);
            gDE[3]->SetPointError(np,0,tmp_err/F2E[ii]);
            gDB[3]->SetPoint(np,xbj[ii],tmp_r1);
            gDB[3]->SetPointError(np,0,tmp_err/F2B[ii]);
         }
         if(abs(tmp_r-1.0)>=0.2||abs(tmp_r1-1.0)>=0.2)
           myfile<<xbj[ii]<<"   "<<Q2[ii]<<"   "<<W2[ii]<<"   "<<F2D[ii]<<"   "<<tmp_err<<"   "<<F2E[ii]<<"   "<<F2B[ii]<<endl;
         np++;
     }
     myfile.close();

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     c1->cd();
     TPad *pad1=new TPad("pad1","",0,0.,1,1);
     pad1->Draw();
     pad1->Divide(1,4,0.0,0.0);
     pad1->cd(1);
     TMultiGraph *mg1=new TMultiGraph();
     gDE[0]->SetMarkerColor(4);
     gDE[0]->SetMarkerStyle(8);
     gDB[0]->SetMarkerColor(1);
     gDB[0]->SetMarkerStyle(8);
     mg1->Add(gDE[0]);
     mg1->Add(gDB[0]);
     mg1->Draw("AP");
     mg1->GetXaxis()->SetLimits(0,1);
     mg1->GetXaxis()->SetLabelSize(0.1);
     mg1->GetYaxis()->SetRangeUser(0.75,1.25);
     mg1->GetYaxis()->SetLabelSize(0.1);
     mg1->GetYaxis()->SetNdivisions(6);
     TLine *l1=new TLine(0,1,1,1);
     l1->SetLineColor(2);
     l1->Draw();
     TLatex *t1 = new TLatex();
     t1->SetNDC();
     t1->SetTextFont(32);
     t1->SetTextSize(0.1);
     t1->SetTextColor(1);
     t1->DrawLatex(0.8,0.8,"Q_{2}<8");

     pad1->cd(2);
     TMultiGraph *mg2=new TMultiGraph();
     gDE[1]->SetMarkerColor(4);
     gDE[1]->SetMarkerStyle(8);
     gDB[1]->SetMarkerColor(1);
     gDB[1]->SetMarkerStyle(8);
     mg2->Add(gDE[1]);
     mg2->Add(gDB[1]);
     mg2->Draw("AP");
     mg2->GetXaxis()->SetLimits(0,1);
     mg2->GetXaxis()->SetLabelSize(0.1);
     mg2->GetYaxis()->SetRangeUser(0.75,1.25);
     mg2->GetYaxis()->SetLabelSize(0.1);
     mg2->GetYaxis()->SetNdivisions(6);
     TLine *l2=new TLine(0,1,1,1);
     l2->SetLineColor(2);
     l2->Draw();
     TLatex *t2 = new TLatex();
     t2->SetNDC();
     t2->SetTextFont(32);
     t2->SetTextSize(0.1);
     t2->SetTextColor(1);
     t2->DrawLatex(0.8,0.8,"8=<Q_{2}<10");

     pad1->cd(3);
     TMultiGraph *mg3=new TMultiGraph();
     gDE[2]->SetMarkerColor(4);
     gDE[2]->SetMarkerStyle(8);
     gDB[2]->SetMarkerColor(1);
     gDB[2]->SetMarkerStyle(8);
     mg3->Add(gDE[2]);
     mg3->Add(gDB[2]);
     mg3->Draw("AP");
     mg3->GetXaxis()->SetLimits(0,1);
     mg3->GetXaxis()->SetLabelSize(0.1);
     mg3->GetYaxis()->SetRangeUser(0.75,1.25);
     mg3->GetYaxis()->SetLabelSize(0.1);
     mg3->GetYaxis()->SetNdivisions(6);
     TLine *l3=new TLine(0,1,1,1);
     l3->SetLineColor(2);
     l3->Draw();
     TLatex *t3 = new TLatex();
     t3->SetNDC();
     t3->SetTextFont(32);
     t3->SetTextSize(0.1);
     t3->SetTextColor(1);
     t3->DrawLatex(0.8,0.8,"10<=Q_{2}<12");

     pad1->cd(4);
     TMultiGraph *mg4=new TMultiGraph();
     gDE[3]->SetMarkerColor(4);
     gDE[3]->SetMarkerStyle(8);
     gDB[3]->SetMarkerColor(1);
     gDB[3]->SetMarkerStyle(8);
     mg4->Add(gDE[3]);
     mg4->Add(gDB[3]);
     mg4->Draw("AP");
     mg4->GetXaxis()->SetLimits(0,1);
     mg4->GetXaxis()->SetLabelSize(0.1);
     mg4->GetYaxis()->SetRangeUser(0.75,1.25);
     mg4->GetYaxis()->SetLabelSize(0.1);
     mg4->GetYaxis()->SetNdivisions(6);
     TLine *l4=new TLine(0,1,1,1);
     l4->SetLineColor(2);
     l4->Draw();
     TLatex *t4 = new TLatex();
     t4->SetNDC();
     t4->SetTextFont(32);
     t4->SetTextSize(0.1);
     t4->SetTextColor(1);
     t4->DrawLatex(0.8,0.8,"12<=Q_{2}<14");
}
