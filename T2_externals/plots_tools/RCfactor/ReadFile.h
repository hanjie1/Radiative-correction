#define MAXBIN 35
TString Yieldpath="/home/hanjie/work/MARATHON/RadCor/T2_externals/OUT/";
int ReadYield(TString filename,int kin,Double_t x[][MAXBIN],Double_t Q2[][MAXBIN],Double_t Born[][MAXBIN],Double_t Rad[][MAXBIN]){
    ifstream file;
    TString myfile=Yieldpath+filename;
    file.open(myfile);
    if(!file.is_open())return 0;

    Ssiz_t from=0;
    TString content,tmp;
    int nn=0;
    int KKin=0;
    if(kin<=5)KKin=kin;
    if(kin==7)KKin=kin-1;
    if(kin==9)KKin=kin-2;
    if(kin==11)KKin=kin-3;
    if(kin==13)KKin=kin-4;
    if(kin==15)KKin=kin-5;
    
    while(tmp.ReadLine(file)){
          if(nn==0){nn++;continue;}
          tmp.Tokenize(content,from," ");
          x[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Q2[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          tmp.Tokenize(content,from," ");
          Born[KKin][nn-1]=atof(content.Data());
          tmp.Tokenize(content,from," ");
          Rad[KKin][nn-1]=atof(content.Data());
//          cout<<x[KKin][nn-1]<<"  "<<Q2[KKin][nn-1]<<"  "<<Born[KKin][nn-1]<<"  "<<Rad[KKin][nn-1]<<endl;
          from=0;
          nn++;
     }
    file.close();
    return 1;

}

