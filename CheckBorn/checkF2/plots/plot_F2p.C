#define MAXBIN 187
void plot_F2p()
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
           F2D[nn-1]=atof(content.Data());
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

     Double_t xx[10]={0.08,0.125,0.175,0.25,0.35,0.45,0.55,0.65,0.75,0.85};
     TGraphErrors *gData[10];
     TGraphErrors *gE[10];
     TGraphErrors *gB[10];
     TGraphErrors *gDE[10];
     TGraphErrors *gDB[10];

     for(int ii=0;ii<10;ii++){
	 gData[ii]=new TGraphErrors();
	 gE[ii]=new TGraphErrors();
	 gB[ii]=new TGraphErrors();
	 gDE[ii]=new TGraphErrors();
	 gDB[ii]=new TGraphErrors();
     }

     ofstream myfile;
     myfile.open("Bad_f2p_SLAC.out");
     myfile<<"x     Q2     W2    F2_data    F2_data_err    F2_f1f217    F2_ineft"<<endl;

     int np[10]={0};
     for(int ii=0;ii<MAXBIN;ii++){
         if(Q2[ii]>14||W2[ii]>14)continue;
         for(int jj=0;jj<10;jj++){
	     if(xbj[ii]!=xx[jj])continue;
             gData[jj]->SetPoint(np[jj],Q2[ii],F2D[ii]);
             Double_t tmp_err=sqrt(F2D_stat[ii]*F2D_stat[ii]+F2D_sys[ii]*F2D_sys[ii]);
             gData[jj]->SetPointError(np[jj],0,tmp_err);
             gE[jj]->SetPoint(np[jj],Q2[ii],F2E[ii]);
             gB[jj]->SetPoint(np[jj],Q2[ii],F2B[ii]);
             gDE[jj]->SetPoint(np[jj],Q2[ii],F2D[ii]/F2E[ii]);
             gDE[jj]->SetPointError(np[jj],0,tmp_err/F2E[ii]);
             gDB[jj]->SetPoint(np[jj],Q2[ii],F2D[ii]/F2B[ii]);
             gDB[jj]->SetPointError(np[jj],0,tmp_err/F2B[ii]);

             if(abs(F2D[ii]/F2E[ii]-1.0)>=0.2||abs(F2D[ii]/F2B[ii]-1.0)>=0.2)
                myfile<<xbj[ii]<<"   "<<Q2[ii]<<"   "<<W2[ii]<<"   "<<F2D[ii]<<"   "<<tmp_err<<"   "<<F2E[ii]<<"   "<<F2B[ii]<<endl;
             np[jj]=np[jj]+1;
             break;
         }
     }
     myfile.close();

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     TMultiGraph *mg=new TMultiGraph();
     for(int ii=0;ii<10;ii++){
         gData[ii]->SetMarkerStyle(8);
         gData[ii]->SetMarkerColor(1);
         gE[ii]->SetMarkerStyle(8);
         gE[ii]->SetMarkerColor(2);
         gB[ii]->SetMarkerStyle(8);
         gB[ii]->SetMarkerColor(4);
         mg->Add(gData[ii]);
         mg->Add(gE[ii]);
         mg->Add(gB[ii]);
     }
     mg->Draw("AP"); 

     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     c2->cd();
     TPad *pad1=new TPad("pad1","",0,0.,1.0,1.0);
     pad1->Draw();
     pad1->Divide(1,5,0.0,0.0);

     for(int ii=0;ii<5;ii++){
         pad1->cd(ii+1);
         TMultiGraph *mg1=new TMultiGraph();
         gDE[ii]->SetMarkerStyle(8);
         gDE[ii]->SetMarkerColor(2);
         gDB[ii]->SetMarkerStyle(22);
         gDB[ii]->SetMarkerColor(4);
         mg1->Add(gDE[ii]);
         mg1->Add(gDB[ii]);
         mg1->Draw("AP"); 
         TLine *l1=new TLine(0,1,7.5,1);
         l1->SetLineColor(2);
         l1->Draw();
         mg1->GetXaxis()->SetLimits(0,7.5);
         mg1->GetXaxis()->SetLabelSize(0.1);
         mg1->GetYaxis()->SetRangeUser(0.747,1.253);
         mg1->GetYaxis()->SetLabelSize(0.09);
         mg1->GetYaxis()->SetNdivisions(6);
         TLatex *t1 = new TLatex();
         t1->SetNDC();
         t1->SetTextFont(32);
         t1->SetTextSize(0.1);
         t1->SetTextColor(1);
         t1->DrawLatex(0.8,0.8,Form("x=%.3f",xx[ii]));
         if(ii==0){
            auto leg1=new TLegend(0.68,0.6,0.78,0.78);
            leg1->AddEntry(gDE[ii],"f1f217","P");
            leg1->AddEntry(gDB[ii],"Bodek","P");
            leg1->Draw();
         }
     }
     c2->Update();

     TCanvas *c3=new TCanvas("c3","c3",1500,1500);
     c3->cd();
     TPad *pad2=new TPad("pad2","",0,0.,1,1);
     pad2->Draw();
     pad2->Divide(1,5,0.0,0.0);

     for(int ii=0;ii<5;ii++){
         pad2->cd(ii+1);
         TMultiGraph *mg1=new TMultiGraph();
         gDE[ii+5]->SetMarkerStyle(8);
         gDE[ii+5]->SetMarkerColor(2);
         gDB[ii+5]->SetMarkerStyle(22);
         gDB[ii+5]->SetMarkerColor(4);
         mg1->Add(gDE[ii+5]);
         mg1->Add(gDB[ii+5]);
         mg1->Draw("AP"); 
         TLine *l1=new TLine(3,1,14.5,1);
         l1->SetLineColor(2);
         l1->Draw();
         mg1->GetXaxis()->SetLimits(3,14.5);
         mg1->GetXaxis()->SetLabelSize(0.1);
         mg1->GetYaxis()->SetRangeUser(0.68,1.32);
         mg1->GetYaxis()->SetLabelSize(0.09);
         mg1->GetYaxis()->SetNdivisions(8);
         TLatex *t1 = new TLatex();
         t1->SetNDC();
         t1->SetTextFont(32);
         t1->SetTextSize(0.1);
         t1->SetTextColor(1);
         t1->DrawLatex(0.8,0.8,Form("x=%.3f",xx[ii+5]));
         if(ii==0){
            auto leg1=new TLegend(0.68,0.6,0.78,0.78);
            leg1->AddEntry(gDE[ii+5],"f1f217","P");
            leg1->AddEntry(gDB[ii+5],"Bodek","P");
            leg1->Draw();
         }

     }
}
