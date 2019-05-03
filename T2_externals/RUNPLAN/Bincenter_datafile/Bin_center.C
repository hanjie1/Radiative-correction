void Bin_center()
{
     Double_t xbj[18]={0.19,0.22,0.25,0.29,0.33,0.36,0.385,0.43,0.48,0.51,0.55,0.59,0.63,0.67,0.7,0.74,0.78,0.82};
     
     ofstream outfile;
     outfile.open("CenterBin.inp");

     outfile<<"Marathon"<<endl;
     outfile<<"Bin center"<<endl;
     outfile<<endl;
     outfile<<endl;
     outfile<<endl;
     outfile<<"E    Ep    theta"<<endl;

     Double_t E0=10.589;
     Double_t M=0.938272;

     for(int ii=0;ii<18;ii++){
         Double_t Q2=14.0*xbj[ii]; 
         Double_t Ep=E0-Q2/(2.0*M*xbj[ii]);
         Double_t theta=asin(sqrt(Q2/(4.0*E0*Ep)))*2*180/TMath::Pi();
         outfile<<fixed<<setprecision(3)<<E0<<" "<<fixed<<setprecision(4)<<Ep<<" "<<theta<<endl;
     }
     outfile.close();
}
