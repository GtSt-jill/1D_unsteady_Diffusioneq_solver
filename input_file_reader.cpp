void input_file_reader()
{
    int i;
//     double sum(0);

//    srand((unsigned int)time(NULL));
//    for(int i=1;i<n_p;i++){
//        x[i]=(double)rand();
//        sum+=x[i];
//    }

//    for(i=1; i<n_p; i++){
//        x[i]=x[i]/sum*L;
//    }

//    x[0]=0.0;
//    for(i=1; i<n_p; i++){
//        x[i]=x[i-1]+x[i];
//    }
//    for(i=0;i<n_p;i++){
//        cout<<x[i]<<endl;
//    }
    // 等間隔メッシュ
    for(i=0;i<n_p;i++){
        x[i]=double(i)/double(n_e)*L;
    }
    
    for(i=0;i<n_e;i++){
        el[2*i]  =i;
        el[2*i+1]=i+1;
    }
}
