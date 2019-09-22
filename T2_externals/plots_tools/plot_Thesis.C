#define MAXBIN 70
void plot_Thesis(){
     TString file1;
     file1="../OUT/ForPlot_H1.out";

     ifstream infile1;
     infile1.open(file1);

     Double_t xbj[MAXBIN],Q2[MAXBIN],XS_EL[MAXBIN],XS_QE[MAXBIN],XS_DIS[MAXBIN];
     Double_t B_QE[MAXBIN],B_DIS[MAXBIN];
     Double_t Rad[MAXBIN],Born[MAXBIN];
     int nn=0;
     Ssiz_t from=0;
     TString content,tmp;
     while(tmp.ReadLine(infile1)){
           if(nn==0){nn++;continue;}
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           xbj[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           Q2[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           Born[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           B_DIS[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           B_QE[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           Rad[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           XS_EL[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           XS_QE[nn-1]=content.Atof();
           tmp.Tokenize(content,from," ");
           XS_DIS[nn-1]=content.Atof();

           nn++;
           from=0;
         }
     infile1.close();
     nn=nn-1;
     TGraph *gR_EL=new TGraph(nn,xbj,XS_EL);
     TGraph *gR_QE=new TGraph(nn,xbj,XS_QE);
     TGraph *gR_DIS=new TGraph(nn,xbj,XS_DIS);
     TGraph *gR=new TGraph(nn,xbj,Rad);
     TGraph *gB=new TGraph(nn,xbj,Born);
     TGraph *gB_QE=new TGraph(nn,xbj,B_QE);
     TGraph *gB_DIS=new TGraph(nn,xbj,B_DIS);

     TCanvas *c1=new TCanvas("c1","c1",1500,1100);
     TMultiGraph *mg=new TMultiGraph();
     gR_EL->SetLineStyle(8);
     gR_EL->SetLineWidth(3);
     gR_EL->SetLineColor(1);

     gR_QE->SetLineStyle(3);
     gR_QE->SetLineWidth(3);
     gR_QE->SetLineColor(1);

     gR_DIS->SetLineStyle(10);
     gR_DIS->SetLineWidth(3);
     gR_DIS->SetLineColor(kRed-3);

     gR->SetLineStyle(1);
     gR->SetLineWidth(3);
     gR->SetLineColor(8);

     gB->SetLineStyle(7);
     gB->SetLineWidth(3);
     gB->SetLineColor(4);
//     gB->SetMarkerStyle(22);
//     gB->SetMarkerColor(2);
  
     mg->Add(gR_EL);
     mg->Add(gR_QE);
     mg->Add(gR_DIS);
     mg->Add(gR);
     mg->Add(gB);

     mg->Draw("AL");
     mg->SetTitle(";Bjorken x;d#sigma/d#Omega/dE_{p}(nB/Sr/GeV)");

     auto leg = new TLegend(0.65,0.5,0.85,0.85);
     leg->AddEntry(gR_EL,"#sigma^{#scale[0.5]{rad}}_{#scale[0.5]{EL}}","L");
     leg->AddEntry(gR_QE,"#sigma^{#scale[0.5]{rad}}_{#scale[0.5]{QE}}","L");
     leg->AddEntry(gR_DIS,"#sigma^{#scale[0.5]{rad}}_{#scale[0.5]{DIS}}","L");
     leg->AddEntry(gR,"#sigma^{#scale[0.5]{rad}}","L");
     leg->AddEntry(gB,"#sigma^{#scale[0.5]{born}}","L");
     leg->Draw();
     leg->SetMargin(0.7);
     //leg->SetBorderSize(1);

    c1->Print("XS_RC.pdf"); 
}
