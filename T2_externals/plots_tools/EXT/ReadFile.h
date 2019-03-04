TString Yieldpath="/home/hanjie/work/MARATHON/RadCor/T2_externals/OUT/EXT_RC/";
int ReadYield(TString filename,Double_t x[],Double_t Q2[],Double_t RC[]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          x[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          RC[nn-1]=atof(content.Data());
          from=0;
          nn++;
     }
    file.close();
    return 1;

}


