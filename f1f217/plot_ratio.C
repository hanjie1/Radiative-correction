#define MAXNUM 500
void plot_ratio()
{
     ifstream infile1;
     ifstream infile2;
     ifstream infile3;
     ifstream infile4;
     ifstream infile5;
     ifstream infile11;
     ifstream infile12;
     ifstream infile13;
     ifstream infile14;
     ifstream infile15;

//     infile1.open("D2.dat");
     infile1.open("rat_q0.5_leps.dat");
     infile2.open("rat_q1_leps.dat");
     infile3.open("rat_q2_leps.dat");
     infile4.open("rat_q3_leps.dat");
     infile5.open("rat_q4_leps.dat");
     infile11.open("rat_q0.5_heps.dat");
     infile12.open("rat_q1_heps.dat");
     infile13.open("rat_q2_heps.dat");
     infile14.open("rat_q3_heps.dat");
     infile15.open("rat_q4_heps.dat");

     Double_t rat05[MAXNUM],rat1[MAXNUM],rat2[MAXNUM],rat3[MAXNUM],rat4[MAXNUM];
     Double_t rat05_err[MAXNUM],rat1_err[MAXNUM],rat2_err[MAXNUM],rat3_err[MAXNUM],rat4_err[MAXNUM];
     Double_t rat105[MAXNUM],rat11[MAXNUM],rat12[MAXNUM],rat13[MAXNUM],rat14[MAXNUM];
     Double_t rat105_err[MAXNUM],rat11_err[MAXNUM],rat12_err[MAXNUM],rat13_err[MAXNUM],rat14_err[MAXNUM];
     Double_t W2_05[MAXNUM],W2_1[MAXNUM],W2_2[MAXNUM],W2_3[MAXNUM],W2_4[MAXNUM];
     Double_t W21_05[MAXNUM],W21_1[MAXNUM],W21_2[MAXNUM],W21_3[MAXNUM],W21_4[MAXNUM];
     Ssiz_t from=0;
     TString content,tmp;
     int nn=0;
     while(tmp.ReadLine(infile1)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W2_05[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           //W2_05[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat05[nn]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat05_err[nn]=atof(content.Data());
           nn++;
           from=0;
     }
     infile1.close();

     int nn1=0;
     while(tmp.ReadLine(infile2)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W2_1[nn1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat1[nn1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat1_err[nn1]=atof(content.Data());
           nn1++;
           from=0;
     }
     infile2.close();

     int nn2=0;
     while(tmp.ReadLine(infile3)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W2_2[nn2]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat2[nn2]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat2_err[nn2]=atof(content.Data());
           nn2++;
           from=0;
     }
     infile3.close();

     int nn3=0;
     while(tmp.ReadLine(infile4)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W2_3[nn3]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat3[nn3]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat3_err[nn3]=atof(content.Data());
           nn3++;
           from=0;
     }
     infile4.close();

     int nn4=0;
     while(tmp.ReadLine(infile5)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W2_4[nn4]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat4[nn4]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat4_err[nn4]=atof(content.Data());
           nn4++;
           from=0;
     }
     infile5.close();

     int nnp=0;
     while(tmp.ReadLine(infile11)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W21_05[nnp]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat105[nnp]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat105_err[nnp]=atof(content.Data());
           nnp++;
           from=0;
     }
     infile11.close();

     int nnp1=0;
     while(tmp.ReadLine(infile12)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W21_1[nnp1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat11[nnp1]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat11_err[nnp1]=atof(content.Data());
           nnp1++;
           from=0;
     }
     infile12.close();

     int nnp2=0;
     while(tmp.ReadLine(infile13)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W21_2[nnp2]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat12[nnp2]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat12_err[nnp2]=atof(content.Data());
           nnp2++;
           from=0;
     }
     infile13.close();

     int nnp3=0;
     while(tmp.ReadLine(infile14)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W21_3[nnp3]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat13[nnp3]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat13_err[nnp3]=atof(content.Data());
           nnp3++;
           from=0;
     }
     infile14.close();

     int nnp4=0;
     while(tmp.ReadLine(infile15)){
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           W21_4[nnp4]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           tmp.Tokenize(content,from," ");
           rat14[nnp4]=atof(content.Data());
           tmp.Tokenize(content,from," ");
           rat14_err[nnp4]=atof(content.Data());
           nnp4++;
           from=0;
     }
     infile15.close();

     TGraphErrors *gl05=new TGraphErrors(nn,W2_05,rat05,0,rat05_err);
     TGraphErrors *gl1=new TGraphErrors(nn1,W2_1,rat1,0,rat1_err);
     TGraphErrors *gl2=new TGraphErrors(nn2,W2_2,rat2,0,rat2_err);
     TGraphErrors *gl3=new TGraphErrors(nn3,W2_3,rat3,0,rat3_err);
     TGraphErrors *gl4=new TGraphErrors(nn4,W2_4,rat4,0,rat4_err);
     TGraphErrors *gh05=new TGraphErrors(nnp,W21_05,rat105,0,rat105_err);
     TGraphErrors *gh1=new TGraphErrors(nnp1,W21_1,rat11,0,rat11_err);
     TGraphErrors *gh2=new TGraphErrors(nnp2,W21_2,rat12,0,rat12_err);
     TGraphErrors *gh3=new TGraphErrors(nnp3,W21_3,rat13,0,rat13_err);
     TGraphErrors *gh4=new TGraphErrors(nnp4,W21_4,rat14,0,rat14_err);

     
     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     c1->cd();
     TPad *pad1=new TPad("pad1","",0,0.,1,1);
     pad1->Draw();
     pad1->Divide(1,5,0.0,0.0);
     pad1->cd(1);
     TMultiGraph *mg1=new TMultiGraph();
     gl05->SetMarkerStyle(8);
     gl05->SetMarkerColor(2);
     gh05->SetMarkerStyle(8);
     gh05->SetMarkerColor(4);
     mg1->Add(gl05);
     mg1->Add(gh05);
     mg1->Draw("AP");
     mg1->GetXaxis()->SetLimits(0,10);
     mg1->GetXaxis()->SetLabelSize(0.1);
     mg1->GetYaxis()->SetRangeUser(0.75,1.25);
     mg1->GetYaxis()->SetLabelSize(0.1);
     mg1->GetYaxis()->SetNdivisions(6);
     TLine *l1=new TLine(0,1,10,1);
     l1->SetLineColor(2);
     l1->Draw();
     TLatex *t1 = new TLatex();
     t1->SetNDC();
     t1->SetTextFont(32);
     t1->SetTextSize(0.1);
     t1->SetTextColor(1);
     t1->DrawLatex(0.8,0.8,"0.1<Q_{2}<0.5");

     pad1->cd(2);
     TMultiGraph *mg2=new TMultiGraph();
     gl1->SetMarkerStyle(8);
     gl1->SetMarkerColor(2);
     gh1->SetMarkerStyle(8);
     gh1->SetMarkerColor(4);
     mg2->Add(gl1);
     mg2->Add(gh1);
     mg2->Draw("AP");
     mg2->GetXaxis()->SetLimits(0,10);
     mg2->GetXaxis()->SetLabelSize(0.1);
     mg2->GetYaxis()->SetRangeUser(0.75,1.25);
     mg2->GetYaxis()->SetLabelSize(0.1);
     mg2->GetYaxis()->SetNdivisions(6);
     TLine *l2=new TLine(0,1,10,1);
     l2->SetLineColor(2);
     l2->Draw();
     TLatex *t2 = new TLatex();
     t2->SetNDC();
     t2->SetTextFont(32);
     t2->SetTextSize(0.1);
     t2->SetTextColor(1);
     t2->DrawLatex(0.8,0.8,"0.5<Q_{2}<1.5");

     pad1->cd(3);
     TMultiGraph *mg3=new TMultiGraph();
     gl2->SetMarkerStyle(8);
     gl2->SetMarkerColor(2);
     gh2->SetMarkerStyle(8);
     gh2->SetMarkerColor(4);
     mg3->Add(gl2);
     mg3->Add(gh2);
     mg3->Draw("AP");
     mg3->GetXaxis()->SetLimits(0,10);
     mg3->GetXaxis()->SetLabelSize(0.1);
     mg3->GetYaxis()->SetRangeUser(0.75,1.25);
     mg3->GetYaxis()->SetLabelSize(0.1);
     mg3->GetYaxis()->SetNdivisions(6);
     TLine *l3=new TLine(0,1,10,1);
     l3->SetLineColor(2);
     l3->Draw();
     TLatex *t3 = new TLatex();
     t3->SetNDC();
     t3->SetTextFont(32);
     t3->SetTextSize(0.1);
     t3->SetTextColor(1);
     t3->DrawLatex(0.8,0.8,"1.5<Q_{2}<2.5");

     pad1->cd(4);
     TMultiGraph *mg4=new TMultiGraph();
     gl3->SetMarkerStyle(8);
     gl3->SetMarkerColor(2);
     gh3->SetMarkerStyle(8);
     gh3->SetMarkerColor(4);
     mg4->Add(gl3);
     mg4->Add(gh3);
     mg4->Draw("AP");
     mg4->GetXaxis()->SetLimits(0,10);
     mg4->GetXaxis()->SetLabelSize(0.1);
     mg4->GetYaxis()->SetRangeUser(0.75,1.25);
     mg4->GetYaxis()->SetLabelSize(0.1);
     mg4->GetYaxis()->SetNdivisions(6);
     TLine *l4=new TLine(0,1,10,1);
     l4->SetLineColor(2);
     l4->Draw();
     TLatex *t4 = new TLatex();
     t4->SetNDC();
     t4->SetTextFont(32);
     t4->SetTextSize(0.1);
     t4->SetTextColor(1);
     t4->DrawLatex(0.8,0.8,"2.5<Q_{2}<3.5");

     pad1->cd(5);
     TMultiGraph *mg5=new TMultiGraph();
     gl4->SetMarkerStyle(8);
     gl4->SetMarkerColor(2);
     gh4->SetMarkerStyle(8);
     gh4->SetMarkerColor(4);
     mg5->Add(gl4);
     mg5->Add(gh4);
     mg5->Draw("AP");
     mg5->GetXaxis()->SetLimits(0,10);
     mg5->GetXaxis()->SetLabelSize(0.1);
     mg5->GetYaxis()->SetRangeUser(0.75,1.25);
     mg5->GetYaxis()->SetLabelSize(0.1);
     mg5->GetYaxis()->SetNdivisions(6);
     TLine *l5=new TLine(0,1,10,1);
     l5->SetLineColor(2);
     l5->Draw();
     TLatex *t5 = new TLatex();
     t5->SetNDC();
     t5->SetTextFont(32);
     t5->SetTextSize(0.1);
     t5->SetTextColor(1);
     t5->DrawLatex(0.8,0.8,"3.5<Q_{2}");
}
