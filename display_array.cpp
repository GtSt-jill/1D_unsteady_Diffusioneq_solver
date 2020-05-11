void display_array_double(vector <double> S, int N){
    int i;
    for(i=0;i<N;i++){
        cout<<S[i]<<endl;
    }
}

void confirm_array_double(vector <double> S, int array_size, string file_name){
    int i;
    char fname[20];
    sprintf(fname,"array/%s.dat",file_name.c_str());
    ofstream out(fname);
    for(i=0;i<array_size;i++){
        out<<S[i]<<endl;
    }
    out.close();
}

void confirm_array_int(vector <int> S, int array_size, string file_name){
    int i;
    char fname[20];
    sprintf(fname,"array/%s.dat",file_name.c_str());
    ofstream out(fname);
    for(i=0;i<array_size;i++){
        out<<S[i]<<endl;
    }
    out.close();
}

// 行列確認用関数
void confirm_Matrix(vector <double> S, string file_name){
    int i,j;
    int node;
    vector <double> S_all;
    S_all.assign(n_p*n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            S_all[n_p*i+node]=S[j];
        }
    }
    
    char fname[20];
    sprintf(fname,"array/%s.dat",file_name.c_str());
    ofstream out(fname);
    for(i=0;i<n_p;i++){
        for(j=0;j<n_p;j++){
            out<<S_all[n_p*i+j]<<" ";
        }
        out<<endl;
    }
    out.close();
}
