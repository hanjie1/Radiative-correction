#define MAXNUM 100
void Compare()
{
    TString file1,file2,file3,file4;
    file1="data/4He_18.dat";
    file2="model_result/He4_18_newfit.out";
    file3="model_result/He4_18_newgsmearing.out";
    file4="model_result/He4_18_newgsmearing.out";

    ifstream infile1,infile2,infile3,infile4;
    infile1.open(file1);
    infile2.open(file2);
    infile3.open(file3);
    infile4.open(file4);

    Double_t xbj[MAXNUM]={0.0},Q2[MAXNUM]={0.0},XS[MAXNUM]={0.0},XS_err[MAXNUM]={0.0};
    Double_t XS_Bodek[MAXNUM]={0.0};
    Double_t XS_gsmearing[MAXNUM]={0.0};
    Double_t XS_new[MAXNUM]={0.0};
    Double_t W2[MAXNUM]={0.0};

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    while(tmp.ReadLine(infile1)){
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          XS[nn]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          XS_err[nn]=atof(content.Data());
          nn++;
          from=0;
    }
    infile1.close();

    nn=0;from=0;
    while(tmp.ReadLine(infile2)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          xbj[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          XS_Bodek[nn-1]=atof(content.Data());
	  W2[nn-1]=0.938272*0.938272+Q2[nn-1]*(1.0/xbj[nn-1]-1);
          nn++;
          from=0;
    }
    infile2.close();

    nn=0;from=0;
    while(tmp.ReadLine(infile3)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          XS_gsmearing[nn-1]=atof(content.Data());
          nn++;
          from=0;
    }
    infile3.close();

    nn=0;from=0;
    while(tmp.ReadLine(infile4)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          XS_new[nn-1]=atof(content.Data());
          nn++;
          from=0;
    }
    infile4.close();

    nn=nn-1;
    TGraphErrors *gData=new TGraphErrors(nn,xbj,XS,0,XS_err);
    TGraphErrors *gModel1=new TGraphErrors(nn,xbj,XS_Bodek,0,0);
    TGraphErrors *gModel2=new TGraphErrors(nn,xbj,XS_gsmearing,0,0);
    TGraphErrors *gModel3=new TGraphErrors(nn,xbj,XS_new,0,0);

    TMultiGraph *mg=new TMultiGraph();
    gData->SetMarkerStyle(8);
    gData->SetMarkerColor(1);
    gData->SetLineColor(1);
    gModel1->SetMarkerStyle(8);
    gModel1->SetMarkerColor(4);
    gModel1->SetLineColor(4);
    gModel2->SetMarkerStyle(8);
    gModel2->SetMarkerColor(4);
    gModel2->SetLineColor(4);
    gModel3->SetMarkerStyle(8);
    gModel3->SetMarkerColor(6);
    gModel3->SetLineColor(6);
    mg->Add(gData);
    mg->Add(gModel1);
//    mg->Add(gModel2);
//    mg->Add(gModel3);
    mg->Draw("APL");
    mg->SetTitle("He4 cross section;xbj;XS");

    auto leg1=new TLegend(0.7,0.6,0.85,0.85);
    leg1->AddEntry(gData,"Data E0=5.766,Theta=18","P");
    leg1->AddEntry(gModel1,"F1F217","P");
//    leg1->AddEntry(gModel2,"RESCSD","P");
//    leg1->AddEntry(gModel3,"original2","P");
    leg1->Draw();


}
