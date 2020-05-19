// 要素番号iの要素の面積を求める関数
double element_area(int i){
    int n0,n1,n2;
    int itmp=2*i;
    double x0,x1;
    double area;
    
    n0=el[itmp];
    n1=el[itmp+1];
    
    x0=x[n0];
    x1=x[n1];
    
    area=x1-x0;
    
    return area;
}

//質量行列の構成
void Calculate_Mmatrix(){
    int i,j,jstart,jend,k;
    int j0,j1,j2;
    int n0,n1;
    int itmp,node_tmp;
    double Area;
    
    M.assign(n_link.size(),0.0);
    
    for(i=0;i<n_e;i++){
        Area=element_area(i);
        itmp=(d_o_f+1)*i;
        n0=el[itmp  ];
        n1=el[itmp+1];
        
        //n0に関して
        for(j=n_link_start[n0];j<n_link_start[n0+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0) j0=j;
            else if(node_tmp==n1) j1=j;
        }
        M[j0]+=Area/3;
        M[j1]+=Area/6;
        
        //n1に関して
        for(j=n_link_start[n1];j<n_link_start[n1+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0) j0=j;
            else if(node_tmp==n1) j1=j;
        }
        M[j0]+=Area/6;
        M[j1]+=Area/3;
    }
    cout<<"M行列構成したよ"<<endl;
    // confirm_Matrix(M,"M");
}
/*
// 剛性行列の構成
void Calculate_Kmatrix(){
    int i,j,jstart,jend,k;
    int j0,j1;
    int n0,n1;

    int itmp,node_tmp;
    double Area;
    double x0,x1;
    
    K.assign(n_link.size(),0.0);
    
    for(i=0;i<n_e;i++){
        Area=element_area(i);
        itmp=2*i;
        n0=el[itmp  ];
        n1=el[itmp+1];
        
        x0=x[n0];
        x1=x[n1];
        
        //n0について
        for(j=n_link_start[n0];j<n_link_start[n0+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0)      j0=j;
            else if(node_tmp==n1) j1=j;
        }
        K[j0]+=Di/(x1-x0);
        K[j1]+=-Di/(x1-x0);
        
        //n1について
        for(j=n_link_start[n1];j<n_link_start[n1+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0)      j0=j;
            else if(node_tmp==n1) j1=j;
        }
        K[j0]+=-Di/(x1-x0);
        K[j1]+=Di/(x1-x0);
    }
    cout<<"K行列構成したよ"<<endl;
    // confirm_Matrix(K,"K");
}
*/

void Calculate_Amatrix(){
    int i,j,jstart,jend,k;
    int j0,j1;
    int n0,n1;

    int itmp,node_tmp;
    double Area;
    double x0,x1;
    
    A.assign(n_link.size(),0.0);
    
    for(i=0; i<n_e; i++){
        Area=element_area(i);
        itmp=2*i;
        n0=el[itmp  ];
        n1=el[itmp+1];
        
        x0=x[n0];
        x1=x[n1];
        
        //n0について
        for(j=n_link_start[n0];j<n_link_start[n0+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0)      j0=j;
            else if(node_tmp==n1) j1=j;
        }
        A[j0]+=-ad*0.5;
        A[j1]+=ad*0.5;
        
        //n1について
        for(j=n_link_start[n1];j<n_link_start[n1+1];j++){
            node_tmp=n_link[j];
            if(node_tmp==n0)      j0=j;
            else if(node_tmp==n1) j1=j;
        }
        A[j0]+=-ad*0.5;
        A[j1]+=ad*0.5;
    }
    cout<<"A行列構成したよ"<<endl;
    // confirm_Matrix(A,"A");
}

// Calculate Correction Matrix generated from SUPG method.
void Calculate_Cormatrix(){
    int i;
    double dx = L / double(n_e);  // element length
    
    A_c.assign(n_link.size(),0.0); // Correction of Advection Matrix
    M_c.assign(n_link.size(),0.0); // Correction of Mass Matrix
    
    if(Di==0.0){
        
    }
    else{
        double Di_opt = ad*dx*(1/tanh(ad*dx/2/Di)-2*Di/ad/dx)/2; // optimized diffusion coefficient
        // 実際は陽に持たなくて良いが，プログラムの可読性を高めるため構成．大規模化したいときはメモリを食うので消したほうが良い．
        A_c.assign(n_link.size(),0.0); // Correction of Advection Matrix
        M_c.assign(n_link.size(),0.0); // Correction of Mass Matrix
        
        for(i=0;i<n_link.size();i++){
            A_c[i] = K[i]*Di_opt/Di;
        }
        
        for(i=0;i<n_link.size();i++){
            M_c[i] = A[i]*Di_opt/ad/ad;
        }
    }
    
    cout<<"A_c行列構成したよ"<<endl;
    cout<<"M_c行列構成したよ"<<endl;
}

//解析領域全体にかかる外力項 体積力にあたる
double Body_flux_function(double X){
//    return 100*(x*x+y+y);//とりあえず定数
    return X*X;//とりあえず定数
}

void Calculate_F_body(){
    int i,j;
    int n0,n1;
    int node_tmp;
    double x0,x1;
    double xg;
    double Area;
    double f;
    
    F_body.assign(n_p,0.0);

    for(i=0;i<n_e;i++){
        Area=element_area(i);
        int itmp=(d_o_f+1)*i;
        n0=el[itmp];
        n1=el[itmp+1];

        x0=x[n0];
        x1=x[n1];

        xg=(x0+x1)/2;
        
        f=Body_flux_function(xg);
        
        F_body[n0]+=f*Area/2;
        F_body[n1]+=f*Area/2;
    }

    cout<<"外力ベクトル構成したよ"<<endl;
}

// 境界での法線方向の熱流速関数 Fの構成に用いる．
double Boundary_flux_function(double x, double y){
    return x*y;
}

// 一旦保留　境界条件はあとで発展させる．
// void Calculate_F_boundary(){
//    F_boundary.assign(n_p,0.0);
// }

void Calculate_matrix(){
    Calculate_Mmatrix();
    // Calculate_Kmatrix();
    // Calculate_Amatrix();
    // Calculate_Cormatrix();
    // Calculate_F_body();
}
