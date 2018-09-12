void xs_ratio()
{
     TString file1,file2,file3,file4;
     file1="../OUT/marathon_example_HE3.out";
     file2="../OUT/marathon_example_T2.out";
     file3="../OUT/marathon_D2.out";
     file4="../OUT/marathon_H1.out";

     string word[13],line;
     ifstream infile1,infile2,infile3,infile4;
     infile1.open(file1);
     infile2.open(file2);
     infile3.open(file3);
     infile4.open(file4);
   
     Double_t He3_xbj[11],He3_Q2[11],He3_born[11],He3_Bin[11],He3_Bqe[11],He3_rad[11];
     Double_t He3_Rel[11],He3_Rqe[11],He3_Rin[11],He3_Ccor[11];
     int nn=0;
     getline(infile1,line);
     while(getline(infile1,line)){
           istringstream str(line);
           str>>word[0]>>word[1]>>word[2];
	   str>>He3_xbj[nn]>>He3_Q2[nn]>>He3_born[nn]>>He3_Bin[nn]>>He3_Bqe[nn];
           str>>He3_rad[nn]>>He3_Rel[nn]>>He3_Rqe[nn]>>He3_Rin[nn]>>He3_Ccor[nn];
           cout<<He3_xbj[nn]<<"  "<<He3_born[nn]<<endl;
	   nn++;
     }
     infile1.close();

     nn=0;
     Double_t H3_xbj[11],H3_Q2[11],H3_born[11],H3_Bin[11],H3_Bqe[11],H3_rad[11];
     Double_t H3_Rel[11],H3_Rqe[11],H3_Rin[11],H3_Ccor[11];
     getline(infile2,line);
     while(getline(infile2,line)){
           istringstream str(line);
           str>>word[0]>>word[1]>>word[2];
           str>>H3_xbj[nn]>>H3_Q2[nn]>>H3_born[nn]>>H3_Bin[nn]>>H3_Bqe[nn];
           str>>H3_rad[nn]>>H3_Rel[nn]>>H3_Rqe[nn]>>H3_Rin[nn]>>H3_Ccor[nn];
           nn++;
     }
     infile2.close();

     nn=0;
     Double_t D2_xbj[11],D2_Q2[11],D2_born[11],D2_Bin[11],D2_Bqe[11],D2_rad[11];
     Double_t D2_Rel[11],D2_Rqe[11],D2_Rin[11],D2_Ccor[11];
     getline(infile3,line);
     while(getline(infile3,line)){
           istringstream str(line);
           str>>word[0]>>word[1]>>word[2];
           str>>D2_xbj[nn]>>D2_Q2[nn]>>D2_born[nn]>>D2_Bin[nn]>>D2_Bqe[nn];
           str>>D2_rad[nn]>>D2_Rel[nn]>>D2_Rqe[nn]>>D2_Rin[nn]>>D2_Ccor[nn];
           nn++;
     }
     infile3.close();

     nn=0;
     Double_t H1_xbj[11],H1_Q2[11],H1_born[11],H1_Bin[11],H1_Bqe[11],H1_rad[11];
     Double_t H1_Rel[11],H1_Rqe[11],H1_Rin[11],H1_Ccor[11];
     getline(infile4,line);
     while(getline(infile4,line)){
           istringstream str(line);
           str>>word[0]>>word[1]>>word[2];
           str>>H1_xbj[nn]>>H1_Q2[nn]>>H1_born[nn]>>H1_Bin[nn]>>H1_Bqe[nn];
           str>>H1_rad[nn]>>H1_Rel[nn]>>H1_Rqe[nn]>>H1_Rin[nn]>>H1_Ccor[nn];
           nn++;
     }
     infile4.close();

     Double_t He3_cor[11],H3_cor[11],D2_cor[11],H1_cor[11],r1_cor[11],r2_cor[11];
     for(int ii=0;ii<11;ii++){
         He3_cor[ii]=He3_rad[ii]/He3_born[ii];
         H3_cor[ii]=H3_rad[ii]/H3_born[ii];
         D2_cor[ii]=D2_rad[ii]/D2_born[ii];
         H1_cor[ii]=H1_rad[ii]/H1_born[ii];
         r1_cor[ii]=D2_cor[ii]/H1_cor[ii];
         r2_cor[ii]=He3_cor[ii]/H3_cor[ii];
     }

     TGraph *hHe3=new TGraph(11,He3_xbj,He3_rad);
     TGraph *hH3=new TGraph(11,H3_xbj,H3_rad);
     TGraph *hD2=new TGraph(11,D2_xbj,D2_rad);
     TGraph *hH1=new TGraph(11,H1_xbj,H1_rad);

     TGraph *hHe3_b=new TGraph(11,He3_xbj,He3_born);
     TGraph *hH3_b=new TGraph(11,H3_xbj,H3_born);
     TGraph *hD2_b=new TGraph(11,D2_xbj,D2_born);
     TGraph *hH1_b=new TGraph(11,H1_xbj,H1_born);

     TGraph *hHe3_c=new TGraph(11,He3_xbj,He3_cor);
     TGraph *hH3_c=new TGraph(11,H3_xbj,H3_cor);
     TGraph *hD2_c=new TGraph(11,D2_xbj,D2_cor);
     TGraph *hH1_c=new TGraph(11,H1_xbj,H1_cor);

     TGraph *hcor1=new TGraph(11,H1_xbj,r1_cor);
     TGraph *hcor2=new TGraph(11,H3_xbj,r2_cor);
  
     TCanvas *c1=new TCanvas("c1","c1",1500,1500);
     c1->SetLogy();
     TMultiGraph *mg=new TMultiGraph();
     hHe3->SetMarkerStyle(8);
     hHe3->SetMarkerColor(4);
     hH3->SetMarkerStyle(8);
     hH3->SetMarkerColor(2);     
     hD2->SetMarkerStyle(8);
     hD2->SetMarkerColor(6);     
     hH1->SetMarkerStyle(8);
     hH1->SetMarkerColor(1);

     hHe3_b->SetMarkerStyle(26);
     hHe3_b->SetMarkerColor(4);
     hH3_b->SetMarkerStyle(26);
     hH3_b->SetMarkerColor(2);
     hD2_b->SetMarkerStyle(26);
     hD2_b->SetMarkerColor(6);
     hH1_b->SetMarkerStyle(26);
     hH1_b->SetMarkerColor(1);
     
     mg->Add(hHe3);
     mg->Add(hH3);
     mg->Add(hD2);
     mg->Add(hH1);
     mg->Add(hHe3_b);
     mg->Add(hH3_b);
     mg->Add(hD2_b);
     mg->Add(hH1_b);
     mg->Draw("AP");
     mg->SetTitle(";xbj;");

     auto leg = new TLegend(0.7,0.7,0.85,0.85);
     leg->AddEntry(hHe3,"He3 sig_rad","P");
     leg->AddEntry(hH3,"H3 sig_rad","P");
     leg->AddEntry(hD2,"D2 sig_rad","P");
     leg->AddEntry(hH1,"H1 sig_rad","P");
     leg->AddEntry(hHe3_b,"He3 sig_born","P");
     leg->AddEntry(hH3_b,"H3 sig_born","P");
     leg->AddEntry(hD2_b,"D2 sig_born","P");
     leg->AddEntry(hH1_b,"H1 sig_born","P");
     leg->Draw();

     TCanvas *c2=new TCanvas("c2","c2",1500,1500);
     TMultiGraph *mg1=new TMultiGraph();
     hHe3_c->SetMarkerStyle(8);
     hHe3_c->SetMarkerColor(4);
     hH3_c->SetMarkerStyle(8);
     hH3_c->SetMarkerColor(2);
     hD2_c->SetMarkerStyle(8);
     hD2_c->SetMarkerColor(6);
     hH1_c->SetMarkerStyle(8);
     hH1_c->SetMarkerColor(1);

     mg1->Add(hHe3_c);
     mg1->Add(hH3_c);
     mg1->Add(hD2_c);
     mg1->Add(hH1_c);
     mg1->Draw("AP");
     mg1->SetTitle("RadCor=sig_rad/sig_born;xbj;");

     auto leg1 = new TLegend(0.7,0.7,0.85,0.85);
     leg1->AddEntry(hHe3_c,"He3","P");
     leg1->AddEntry(hH3_c,"H3","P");
     leg1->AddEntry(hD2_c,"D2","P");
     leg1->AddEntry(hH1_c,"H1","P");
     leg1->Draw();


     TCanvas *c3=new TCanvas("c3","c3",1500,1500);
     c3->Divide(2,1);
     c3->cd(1);
     hcor1->SetMarkerStyle(8);
     hcor1->SetMarkerColor(4);
     hcor1->Draw("AP");
     hcor1->SetTitle("RadCor(D2)/RadCor(H1);xbj;");

     c3->cd(2);
     hcor2->SetMarkerStyle(8);
     hcor2->SetMarkerColor(4);
     hcor2->Draw("AP");
     hcor2->SetTitle("RadCor(He3)/RadCor(H3);xbj;");
}
