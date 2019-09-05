void Bin_center_Dp()
{
     Double_t xbj[6]={0.19,0.22,0.25,0.29,0.32,0.34};
     
     ofstream outfile;
     outfile.open("CenterBin_Dp.inp");

     outfile<<"Marathon"<<endl;
     outfile<<"Bin center"<<endl;
     outfile<<endl;
     outfile<<endl;
     outfile<<endl;
     outfile<<"E    Ep    theta"<<endl;

     Double_t E0=10.589;
     Double_t M=0.938272;

     for(int ii=0;ii<6;ii++){
         Double_t Q2=14.0*xbj[ii]; 
         Double_t Ep=E0-Q2/(2.0*M*xbj[ii]);
         Double_t theta=asin(sqrt(Q2/(4.0*E0*Ep)))*2*180/TMath::Pi();
         outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
     }
     outfile.close();
}
