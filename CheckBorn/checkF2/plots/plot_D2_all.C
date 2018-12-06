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

     TGraphErrors *gData=new TGraphErrors(MAXBIN,xbj,F2D,0,F2D_stat);
     TGraphErrors *gE=new TGraphErrors(MAXBIN,xbj,F2E,0,0);
     TGraphErrors *gB=new TGraphErrors(MAXBIN,xbj,F2B,0,0);
     TGraphErrors *gDE=new TGraphErrors(MAXBIN);
     TGraphErrors *gDB=new TGraphErrors(MAXBIN);

     ofstream myfile;
     myfile.open("Bad_D2_whitlow.out");
     myfile<<"x     Q2     W2    F2_data    F2_data_err    F2_f1f217    F2_ineft"<<endl;
     int np=0;
     for(int ii=0;ii<MAXBIN;ii++){
         Double_t tmp_err=sqrt(F2D_stat[ii]*F2D_stat[ii]+F2D_sys[ii]*F2D_sys[ii]);
	 Double_t tmp_r=F2D[ii]/F2E[ii];
         gDE->SetPoint(np,xbj[ii],tmp_r);
         gDE->SetPointError(np,0,tmp_err/F2E[ii]);
         Double_t tmp_r1=F2D[ii]/F2B[ii];
         gDB->SetPoint(np,xbj[ii],tmp_r1);
         gDB->SetPointError(np,0,tmp_err/F2B[ii]);
         if(abs(tmp_r-1.0)>=0.2||abs(tmp_r1-1.0)>=0.2)
           myfile<<xbj[ii]<<"   "<<Q2[ii]<<"   "<<W2[ii]<<"   "<<F2D[ii]<<"   "<<tmp_err<<"   "<<F2E[ii]<<"   "<<F2B[ii]<<endl;
         np++;
     }
     myfile.close();

     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     TMultiGraph *mg=new TMultiGraph();
     gData->SetMarkerStyle(8);
     gData->SetMarkerColor(1);
     gE->SetMarkerStyle(8);
     gE->SetMarkerColor(2);
     gB->SetMarkerStyle(8);
     gB->SetMarkerColor(4);
     mg->Add(gData);
     mg->Add(gE);
     mg->Add(gB);
     mg->Draw("AP"); 

     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     TMultiGraph *mg1=new TMultiGraph();
     gDE->SetMarkerStyle(8);
     gDE->SetMarkerColor(2);
     gDB->SetMarkerStyle(8);
     gDB->SetMarkerColor(4);
     mg1->Add(gDE);
     mg1->Add(gDB);
     mg1->Draw("AP");

}
