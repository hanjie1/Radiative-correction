#define MAXBIN 610
void plot_F2d()
{
     TString filename;
     ifstream file1;
     cout<<"Input file name: ";
     cin>>filename;
     file1.open(Form("../OUT/%s",filename.Data()));

     Double_t xbj[MAXBIN]={0.0},Q2[MAXBIN]={0.0},W2[MAXBIN]={0.0};
     Double_t F2D[MAXBIN]={0.0},F2D_stat[MAXBIN]={0.0},F2D_sys[MAXBIN]={0.0};
     Double_t F2E[MAXBIN]={0.0},F2B[MAXBIN]={0.0};
     Double_t F2D_err[MAXBIN]={0.0};
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
           F2D[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_stat[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2D_sys[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2E[nn-1]=atof(content.Data())/2.0;
           tmp.Tokenize(content,from," ");
           F2B[nn-1]=atof(content.Data());

	   F2D_err[nn-1]=F2D[nn-1]*sqrt(F2D_stat[nn-1]*F2D_stat[nn-1]+F2D_sys[nn-1]*F2D_sys[nn-1]);
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

     int np1=0,np2=0,np3=0,np4=0;
     for(int ii=0;ii<MAXBIN;ii++){
	 if(Q2[ii]>=14||W2[ii]>=14)continue;
	 Double_t tmp_r=F2D[ii]/F2E[ii];
         Double_t tmp_r1=F2D[ii]/F2B[ii];
         if(Q2[ii]<8){
            gData[0]->SetPoint(np1,xbj[ii],F2D[ii]);
            gData[0]->SetPointError(np1,0,F2D_err[ii]);
            gE[0]->SetPoint(np1,xbj[ii],F2E[ii]);
            gB[0]->SetPoint(np1,xbj[ii],F2B[ii]);
            gDE[0]->SetPoint(np1,xbj[ii],tmp_r);
            gDE[0]->SetPointError(np1,0,F2D_err[ii]/F2E[ii]);
            gDB[0]->SetPoint(np1,xbj[ii],tmp_r1);
            gDB[0]->SetPointError(np1,0,F2D_err[ii]/F2B[ii]);
	    np1++;
         }
         if(Q2[ii]<10 && Q2[ii]>=8){
            gData[1]->SetPoint(np2,xbj[ii],F2D[ii]);
            gData[1]->SetPointError(np2,0,F2D_err[ii]);
            gE[1]->SetPoint(np2,xbj[ii],F2E[ii]);
            gB[1]->SetPoint(np2,xbj[ii],F2B[ii]);
            gDE[1]->SetPoint(np2,xbj[ii],tmp_r);
            gDE[1]->SetPointError(np2,0,F2D_err[ii]/F2E[ii]);
            gDB[1]->SetPoint(np2,xbj[ii],tmp_r1);
            gDB[1]->SetPointError(np2,0,F2D_err[ii]/F2B[ii]);
	    np2++;
         }
         if(Q2[ii]<12 && Q2[ii]>=10){
            gData[2]->SetPoint(np3,xbj[ii],F2D[ii]);
            gData[2]->SetPointError(np3,0,F2D_err[ii]);
            gE[2]->SetPoint(np3,xbj[ii],F2E[ii]);
            gB[2]->SetPoint(np3,xbj[ii],F2B[ii]);
            gDE[2]->SetPoint(np3,xbj[ii],tmp_r);
            gDE[2]->SetPointError(np3,0,F2D_err[ii]/F2E[ii]);
            gDB[2]->SetPoint(np3,xbj[ii],tmp_r1);
            gDB[2]->SetPointError(np3,0,F2D_err[ii]/F2B[ii]);
	    np3++;
         }
         if(Q2[ii]<14 && Q2[ii]>=12){
            gData[3]->SetPoint(np4,xbj[ii],F2D[ii]);
            gData[3]->SetPointError(np4,0,F2D_err[ii]);
            gE[3]->SetPoint(np4,xbj[ii],F2E[ii]);
            gB[3]->SetPoint(np4,xbj[ii],F2B[ii]);
            gDE[3]->SetPoint(np4,xbj[ii],tmp_r);
            gDE[3]->SetPointError(np4,0,F2D_err[ii]/F2E[ii]);
            gDB[3]->SetPoint(np4,xbj[ii],tmp_r1);
            gDB[3]->SetPointError(np4,0,F2D_err[ii]/F2B[ii]);
 	    np4++;
         }
     }

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     c1->cd();
     TPad *pad1=new TPad("pad1","",0,0.,1,1);
     pad1->Draw();
     pad1->Divide(1,4,0.0,0.0);
     pad1->cd(1);
     TMultiGraph *mg1=new TMultiGraph();
     gDE[0]->SetMarkerColor(2);
     gDE[0]->SetMarkerStyle(8);
     gDB[0]->SetMarkerColor(4);
     gDB[0]->SetMarkerStyle(22);
     mg1->Add(gDE[0]);
     mg1->Add(gDB[0]);
     mg1->Draw("AP");
     mg1->GetXaxis()->SetLimits(0,1);
     mg1->GetXaxis()->SetLabelSize(0.1);
     mg1->GetYaxis()->SetRangeUser(0.78,1.22);
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
     t1->DrawLatex(0.9,0.8,"Q_{2}<8");
     auto leg1=new TLegend(0.2,0.6,0.38,0.78);
     leg1->AddEntry(gDE[0],"f1f217","P");
     leg1->AddEntry(gDB[0],"Bodek","P");
     leg1->Draw();


     pad1->cd(2);
     TMultiGraph *mg2=new TMultiGraph();
     gDE[1]->SetMarkerColor(2);
     gDE[1]->SetMarkerStyle(8);
     gDB[1]->SetMarkerColor(4);
     gDB[1]->SetMarkerStyle(22);
     mg2->Add(gDE[1]);
     mg2->Add(gDB[1]);
     mg2->Draw("AP");
     mg2->GetXaxis()->SetLimits(0,1);
     mg2->GetXaxis()->SetLabelSize(0.1);
     mg2->GetYaxis()->SetRangeUser(0.78,1.22);
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
     t2->DrawLatex(0.9,0.8,"8=<Q_{2}<10");

     pad1->cd(3);
     TMultiGraph *mg3=new TMultiGraph();
     gDE[2]->SetMarkerColor(2);
     gDE[2]->SetMarkerStyle(8);
     gDB[2]->SetMarkerColor(4);
     gDB[2]->SetMarkerStyle(22);
     mg3->Add(gDE[2]);
     mg3->Add(gDB[2]);
     mg3->Draw("AP");
     mg3->GetXaxis()->SetLimits(0,1);
     mg3->GetXaxis()->SetLabelSize(0.1);
     mg3->GetYaxis()->SetRangeUser(0.78,1.22);
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
     t3->DrawLatex(0.9,0.8,"10<=Q_{2}<12");

     pad1->cd(4);
     TMultiGraph *mg4=new TMultiGraph();
     gDE[3]->SetMarkerColor(2);
     gDE[3]->SetMarkerStyle(8);
     gDB[3]->SetMarkerColor(4);
     gDB[3]->SetMarkerStyle(22);
     mg4->Add(gDE[3]);
     mg4->Add(gDB[3]);
     mg4->Draw("AP");
     mg4->GetXaxis()->SetLimits(0,1);
     mg4->GetXaxis()->SetLabelSize(0.1);
     mg4->GetYaxis()->SetRangeUser(0.78,1.22);
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
     t4->DrawLatex(0.9,0.8,"12<=Q_{2}<14");

     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     c2->cd();
     TPad *pad2=new TPad("pad2","",0,0.,1,1);
     pad2->Draw();
     pad2->Divide(1,4,0.0,0.0);
     pad2->cd(1);
     TMultiGraph *mg5=new TMultiGraph();
     gData[0]->SetMarkerColor(2);
     gData[0]->SetMarkerStyle(8);
     gE[0]->SetMarkerColor(4);
     gE[0]->SetMarkerStyle(8);
     gB[0]->SetMarkerColor(1);
     gB[0]->SetMarkerStyle(8);
     mg5->Add(gData[0]);
     mg5->Add(gE[0]);
     mg5->Add(gB[0]);
     mg5->Draw("AP");
     mg5->GetXaxis()->SetLimits(0,1);
     mg5->GetXaxis()->SetLabelSize(0.1);
//     mg5->GetYaxis()->SetRangeUser(0.5,1.5);
//     mg5->GetYaxis()->SetLabelSize(0.1);
//     mg5->GetYaxis()->SetNdivisions(6);
     TLatex *t5 = new TLatex();
     t5->SetNDC();
     t5->SetTextFont(32);
     t5->SetTextSize(0.1);
     t5->SetTextColor(1);
     t5->DrawLatex(0.8,0.8,"Q_{2}<8");

     pad2->cd(2);
     TMultiGraph *mg6=new TMultiGraph();
     gData[1]->SetMarkerColor(2);
     gData[1]->SetMarkerStyle(8);
     gE[1]->SetMarkerColor(4);
     gE[1]->SetMarkerStyle(8);
     gB[1]->SetMarkerColor(1);
     gB[1]->SetMarkerStyle(8);
     mg6->Add(gData[1]);
     mg6->Add(gE[1]);
     mg6->Add(gB[1]);
     mg6->Draw("AP");
     mg6->GetXaxis()->SetLimits(0,1);
     mg6->GetXaxis()->SetLabelSize(0.1);
//     mg5->GetYaxis()->SetRangeUser(0.5,1.5);
//     mg5->GetYaxis()->SetLabelSize(0.1);
//     mg5->GetYaxis()->SetNdivisions(6);
     TLatex *t6 = new TLatex();
     t6->SetNDC();
     t6->SetTextFont(32);
     t6->SetTextSize(0.1);
     t6->SetTextColor(1);
     t6->DrawLatex(0.8,0.8,"8<=Q_{2}<10");

     pad2->cd(3);
     TMultiGraph *mg7=new TMultiGraph();
     gData[2]->SetMarkerColor(2);
     gData[2]->SetMarkerStyle(8);
     gE[2]->SetMarkerColor(4);
     gE[2]->SetMarkerStyle(8);
     gB[2]->SetMarkerColor(1);
     gB[2]->SetMarkerStyle(8);
     mg7->Add(gData[2]);
     mg7->Add(gE[2]);
     mg7->Add(gB[2]);
     mg7->Draw("AP");
     mg7->GetXaxis()->SetLimits(0,1);
     mg7->GetXaxis()->SetLabelSize(0.1);
//     mg7->GetYaxis()->SetRangeUser(0.5,1.5);
//     mg7->GetYaxis()->SetLabelSize(0.1);
//     mg7->GetYaxis()->SetNdivisions(6);
     TLatex *t7 = new TLatex();
     t7->SetNDC();
     t7->SetTextFont(32);
     t7->SetTextSize(0.1);
     t7->SetTextColor(1);
     t7->DrawLatex(0.8,0.8,"10<=Q_{2}<12");

     pad2->cd(4);
     TMultiGraph *mg8=new TMultiGraph();
     gData[3]->SetMarkerColor(2);
     gData[3]->SetMarkerStyle(8);
     gE[3]->SetMarkerColor(4);
     gE[3]->SetMarkerStyle(8);
     gB[3]->SetMarkerColor(1);
     gB[3]->SetMarkerStyle(8);
     mg8->Add(gData[3]);
     mg8->Add(gE[3]);
     mg8->Add(gB[3]);
     mg8->Draw("AP");
     mg8->GetXaxis()->SetLimits(0,1);
     mg8->GetXaxis()->SetLabelSize(0.1);
//     mg8->GetYaxis()->SetRangeUser(0.5,1.5);
//     mg8->GetYaxis()->SetLabelSize(0.1);
//     mg8->GetYaxis()->SetNdivisions(6);
     TLatex *t8 = new TLatex();
     t8->SetNDC();
     t8->SetTextFont(32);
     t8->SetTextSize(0.1);
     t8->SetTextColor(1);
     t8->DrawLatex(0.8,0.8,"12<=Q_{2}<14");

}
