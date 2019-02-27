#define MAXBIN 28
void plot_HE3()
{
     Double_t x[MAXBIN],F2A[MAXBIN],eF2A[MAXBIN];
     Double_t F2_model[MAXBIN],F2_KP[MAXBIN],ratio1[MAXBIN],ratio2[MAXBIN];

     ifstream file1;
     file1.open("OUT/XS_18deg.out");

     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(file1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
	   if(atof(content.Data())>0.95){from=0;continue;}
           x[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           F2A[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           eF2A[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2_model[nn-1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           F2_KP[nn-1]=atof(content.Data());

	   ratio1[nn-1]=F2_model[nn-1]/F2A[nn-1];
	   ratio2[nn-1]=F2_KP[nn-1]/F2A[nn-1];

           nn++;
           from=0;
         }
     file1.close();

     TGraphErrors *gF2A=new TGraphErrors(nn-1,x,F2A,0,eF2A);
     TGraph *gF2A_model=new TGraph(nn-1,x,F2_model);
     TGraph *gF2A_KP=new TGraph(nn-1,x,F2_KP);
     TGraph *gratio_model=new TGraph(nn-1,x,ratio1);
     TGraph *gratio_KP=new TGraph(nn-1,x,ratio2);

     TCanvas *c1=new TCanvas("c1");
     c1->Divide(2,1);
     c1->cd(1);
     TMultiGraph *mg=new TMultiGraph();
     gF2A->SetMarkerColor(2);
     gF2A->SetMarkerStyle(8);
     gF2A_model->SetMarkerColor(4);
     gF2A_model->SetMarkerStyle(8);
     gF2A_KP->SetMarkerColor(1);
     gF2A_KP->SetMarkerStyle(8);
     mg->Add(gF2A);
     mg->Add(gF2A_model);
     mg->Add(gF2A_KP);
     mg->Draw("AP");
     mg->SetTitle("18deg");

   auto leg=new TLegend(0.7,0.6,0.85,0.85);
   leg->AddEntry(gF2A,"F2_He3(DATA)","P");
   leg->AddEntry(gF2A_model,"F2_He3(no res)","P");
   leg->AddEntry(gF2A_KP,"F2_He3(w/ res)","P");
   leg->Draw();

     c1->cd(2);
     TMultiGraph *mg1=new TMultiGraph();
     gratio_KP->SetMarkerColor(2);
     gratio_KP->SetMarkerStyle(8);
     gratio_model->SetMarkerColor(4);
     gratio_model->SetMarkerStyle(8);
     mg1->Add(gratio_KP);
     mg1->Add(gratio_model);
     mg1->Draw("AP");
     mg1->SetTitle("18deg");

   auto leg1=new TLegend(0.7,0.6,0.85,0.85);
   leg1->AddEntry(gratio_model,"F2_He3(no res)/Data","P");
   leg1->AddEntry(gratio_KP,"F2_He3(w/ res)/Data","P");
   leg1->Draw();




}
