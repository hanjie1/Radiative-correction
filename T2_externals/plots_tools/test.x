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

