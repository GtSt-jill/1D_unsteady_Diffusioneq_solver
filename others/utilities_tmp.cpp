// * 数値誤差を回避できる内積計算のための関数
double dot_product(vector <double> &v, vector <double> &w, vector <double> &tmp, int n, int lg2n)
{
    int i,ii;
    for(i=0; i<n_p; i++)
    {
        tmp[i]=v[i]*w[i];            //vとwの内積
    }
    for(ii=0; ii<lg2n; ii++)                    //lg2n回繰り返す
    {
        for(i=0; i<n/2; i++){tmp[i]+=tmp[(n-1)-i];}
        n/=2;
    }
    return tmp[0];
}

// In the case that advection term is considerd, it is liable to forget that the A of a equation AX=b is asymmetric, so CG method can't be applied.
// Jacobi method for static analysis
void Jacobi_static(){
    int i,j;
    int node;
    int iter;
    double error;
    
    vector <double> KA; //K_matとA_matの和
    KA.assign(n_link.size(),0.0);
    for(i=0 ;i<n_link.size(); i++){
        KA[i]=K[i]+A[i]+A_c[i];
    }
    
    // 境界条件の代入
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0) u[i]=x[i];
        else u[i]=0.0;
    }
    
    vector <double> KAu_b;
    KAu_b.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            KAu_b[i]+=KA[j]*u[node];
        }
    }
    
    // post-processing considering B.C.
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(node==i){KA[j]=1.0;}
                else{KA[j]=0.0;}
            }
        }
        else{
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(set_bc_list[node]==0){
                    KA[j]=0.0;
                }
            }
        }
    }
    // For confirming KA matrix
    // int node_tmp;
    // vector <double> KA_all;
    // KA_all.assign(n_p*n_p,0.0);
    // for(i=0;i<n_p;i++){
    //     for(j=n_link_start[i];j<n_link_start[i+1];j++){
    //        node_tmp=n_link[j];
    //        KA_all[n_p*i+node_tmp]=KA[j];
    //    }
    // }
    // ofstream out1("array/KA.dat");
    // for(i=0;i<n_p;i++){
    //     for(j=0;j<n_p;j++){
    //        out1<<setw(4)<<KA_all[n_p*i+j]<<" ";
    //     }
    //     out1<<endl;
    // }
    // out1.close();
    
    // KAの対角項と非対角項の生成　要確認
    vector <double> KA_d; // 対角項
    vector <double> KA_t; // 非対角項
    KA_d.assign(n_p,0.0);
    KA_t.assign(n_link.size(),0.0);
    
    for(i=0;i<n_link.size();i++){
        KA_t[i]=KA[i];
    }
    
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            if(i==node){
                KA_d[i]=KA[j];
                KA_t[j]=KA[j]-KA_d[i];
            }
        }
    }
    
    // Arrays for iterative calcuration
    vector <double> b_tmp;
    vector <double> u_tmp;
    b_tmp.assign(n_p,0.0);
    u_tmp.assign(n_p,0.0);
    
    //-------------- iterative calculation starts ---------------
    for(iter=1;iter<Iter_Max;iter++){
        for(i=0;i<n_p;i++){
            u_tmp[i]=u[i];
        }
        
        for(i=0;i<n_p;i++){
            b_tmp[i]=0.0; // 足し込みしてる場合は初期化忘れず！
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                b_tmp[i]-=KA_t[j]*u[node];
            }
        }
        
        for(i=0;i<n_p;i++){
            b_tmp[i]-=KAu_b[i];
        }
        
        error=0.0; // ! ATTENTION : Initializing an error
        for(i=0;i<n_p;i++){
            u[i]=b_tmp[i]/KA_d[i]; // update a solution
            error+=sqrt((u[i]-u_tmp[i])*(u[i]-u_tmp[i])); // calculate an error
        }
        
        if(iter==1) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        if(iter%500==0) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        if(error<EPS){break;} // Convergence Judgement
    }
    // --------------- iterative calculation ends ----------------
    
    //--------- 境界条件の再代入 ----------
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=x[i];
        }
    }
    
    //------------- file-export --------------
    ofstream out("array/u.csv");
    out<<"x"<<","<<"SUPG method"<<","<<"Analytical Solution"<<endl;
    for(i=0;i<n_p;i++){
        out<<x[i]<<","<<u[i]<<","<<(exp(ad/Di*x[i])-1)/(exp(ad/Di*L)-1)<<endl;
    }
    out.close();
}

// Jacobi method for dynamical analysis
void Jacobi_dynamics(vector <double> b){
    int i,j;
    int node;
    int iter;
    double error;
    
    // 境界条件の適用
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0) u[i]=0.0;
        else u[i]=0.0;
    }
    
    //
    vector <double> M_all;
    M_all.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        M_all[i]=M[i]+M_c[i];
    }
    
    // 境界条件による生まれる外力項の計算
    vector <double> M_allu_b;
    M_allu_b.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            M_allu_b[i]+=M_all[j]*u[node];
        }
    }
    
    // post-processing considering B.C.
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(node==i){M_all[j]=1.0;}
                else{M_all[j]=0.0;}
            }
        }
        else{
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(set_bc_list[node]==0){
                    M_all[j]=0.0;
                }
            }
        }
    }
    
    // KAの対角項と非対角項の生成　要確認
    vector <double> M_all_d; // 対角項
    vector <double> M_all_t; // 非対角項
    M_all_d.assign(n_p,0.0);
    M_all_t.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        M_all_t[i]=M_all[i];
    }
    
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            if(i==node){
                M_all_d[i]=M_all[j];
                M_all_t[j]=M_all[j]-M_all_d[i];
            }
        }
    }
    
    // Arrays for iterative calcuration
    vector <double> b_tmp;
    vector <double> u_tmp;
    b_tmp.assign(n_p,0.0);
    u_tmp.assign(n_p,0.0);
    
    //-------------- iterative calculation starts ---------------
    for(iter=1;iter<Iter_Max;iter++){
        // memorize vector u of the previous step
        for(i=0;i<n_p;i++){
            u_tmp[i]=u[i];
        }
        
        for(i=0;i<n_p;i++){
            b_tmp[i]=b[i]; // ! 足し込みしてる場合は初期化忘れず！
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                b_tmp[i]-=M_all_t[j]*u[node];
            }
        }
        
        for(i=0;i<n_p;i++){
            b_tmp[i]-=M_allu_b[i];
        }
        
        error=0.0; // ! ATTENTION : Initializing an error
        for(i=0;i<n_p;i++){
            u[i]=b_tmp[i]/M_all_d[i]; // update a solution
            error+=sqrt((u[i]-u_tmp[i])*(u[i]-u_tmp[i])); // calculate an error
        }
        
        // if(iter==1) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        // if(iter%500==0) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        if(error<EPS){break;} // Convergence Judgement
    }
    // --------------- iterative calculation ends ----------------
    
    //--------- 境界条件の再代入 ----------
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=0.0;
        }
    }
}

// Dynamical analysis function.
void Dynamics(){
    int i,j;
    int node;
    int n;
    
    // intitial condtions of u
    for(i=0;i<n_p;i++){
        if(x[i]<=L/2) u[i]=x[i];
        else if(x[i]>L/2) u[i]=-x[i]+L;
    }
    // 初期条件確認用
    // ofstream out("test.dat");
    // for(i=0;i<n_p;i++){
    //    out<<x[2*i]<<" "<<x[2*i+1]<<" "<<u[i]<<endl;
    // }
    // out.close();
    
    // 毎ステップで解ベクトルに乗ずる行列．時間によらない定数のためnのfor-loop文の外で定義．
    vector <double> MdtS;
    MdtS.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        MdtS[i]=M[i]+M_c[i]-dt*(A[i]+A_c[i]+K[i]);
    }
    
    vector <double> MdtSu;
    MdtSu.assign(n_p,0.0);
    
    // file exporting
    
    //--------------- dynamical analysis start ----------------
    for(n=1;n<1001;n++){
        if(n%50==0) printf("step number : %d\n",n);
        
        // 前進差分法 θ=0
        for(i=0;i<n_p;i++){
            MdtSu[i]=0.0;
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                MdtSu[i]+=MdtS[j]*u[node];
            }
        }
        
        // perform Jacobi method
        Jacobi_dynamics(MdtSu);
        
        // File-initializing
        if(n==1){
            ofstream outtmp("array/u_dynamics.dat");
            outtmp.close();
        }
        // File-exporting
        ofstream out("array/u_dynamics.dat",ios::app);
        if(n%1==0){
            for(i=0;i<n_p;i++){
                out<<x[i]<<" "<<u[i]<<endl;
            }
        }
        out.close();
    } // end of the for-loop of n
    // -------------- dynamical analysis end -------------------
}

//Ku=Fを解いている共役勾配法 静的問題に用いる．※Kが対称行列のときのみ適用可能．
void Conjugate_gradient_method(){
    int i,j;
    int node;
    vector <double> KA; //K_matとA_matの和
    KA.assign(n_link.size(),0.0);
    for(i=0 ;i<n_link.size(); i++){
        KA[i]=K[i]+A[i];
    }
    
    int node_tmp;
    vector <double> KA_all;
    KA_all.assign(n_p*n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node_tmp=n_link[j];
            KA_all[n_p*i+node_tmp]=KA[j];
        }
    }
    ofstream out1("array/KA.dat");
    for(i=0;i<n_p;i++){
        for(j=0;j<n_p;j++){
            out1<<setw(4)<<KA_all[n_p*i+j]<<" ";
        }
        out1<<endl;
    }
    out1.close();
    
    //--------- 内積計算用 -------------
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのtemporary vector
    n_p_Lg2 = (int)ceil(log((double)n_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);
    //--------------------------------
    
    //境界条件の設定
    for(i=0;i<n_p;i++){u[i]=0.0;}
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=x[i];
        }
    }
    
    vector <double> KAu; //KAとuの積
    KAu.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            KAu[i]+=KA[j]*u[node];
        }
    }
    
    vector <double> KAtKAu;
    KAtKAu.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            KAtKAu[node]+=KA[j]*KAu[i];
        }
    }
    
    // conjugate gradient method
    int iter;
    double alpha,beta;
    double c1,c2,c3;
    double err;
    
    vector <double> r;
    vector <double> p;
    vector <double> KAp;//KA*p
    vector <double> KAtKAp;
    
    r.assign(n_p,0.0);
    p.assign(n_p,0.0);
    KAp.assign(n_p,0.0);
    KAtKAp.assign(n_p,0.0);
    
    for(i=0;i<n_p;i++){u[i]=0.0;}
    
    for(i=0;i<n_p;i++){
        r[i]=-KAtKAu[i];
        p[i]=-KAtKAu[i];
    }
    //------------------- Iterative calculation starts ------------------
    for(iter=1;iter<Iter_Max;iter++){
        //境界条件のケア
        for(i=0;i<n_p;i++){
            if(set_bc_list[i]==0){
                p[i]=0.0;
            }
        }
        
        for(i=0;i<n_p;i++){
            KAp[i]=0.0;
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                KAp[i]+=KA[j]*p[node];
            }
        }
        for(i=0;i<n_p;i++){KAtKAp[i]=0.0;}
        for(i=0;i<n_p;i++){
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                KAtKAp[node]+=KA[j]*KAp[i];
            }
        }
        
        c1=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        c2=dot_product(p,KAtKAp,tmp_vec,d_p_terms,n_p_Lg2);
        
        alpha=c1/c2;
        
        for(i=0;i<n_p;i++){
            u[i]+=alpha*p[i];
            r[i]-=alpha*KAtKAp[i];
        }
        
        // 境界条件のケア
        for(i=0;i<n_p;i++){if(set_bc_list[i]==0){r[i]=0.0;}}
        
        err=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        err/=n_p;
        
        if(iter==1) printf("LOOP : %d    %.10e\n", iter, err);
        if(iter%50==0) printf("LOOP : %d    %.10e\n", iter, err);
        if(EPS > err) break;
        
        c3=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        beta=c3/c1;
        
        for(i=0; i<n_p; i++){p[i]=r[i]+beta*p[i];}
    }
    //----------------- Iterative calculation ends ------------------
    
    //--------- 境界条件の再代入 ----------
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=x[i];
        }
    }
    
    //------------- file-export --------------
    ofstream out("array/u.csv");
    for(i=0;i<n_p;i++){
        out<<x[i]<<","<<u[i]<<endl;
    }
    out.close();
}

// Mx=bを解く共役勾配法 動的解析に用いる
void cg_solver_dynamics(vector <double> b){
    int i,j;
    int node;
    
    //--------- 内積計算用 -------------
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのテンポラリーベクター
    n_p_Lg2 = (int)ceil(log((double)n_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);
    //--------------------------------
    
    //境界条件の設定
    for(i=0;i<n_p;i++){u[i]=0.0;}
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=0.0;
        }
    }
    
    // 境界条件を考慮するための外力ベクトルMuの構成
    vector <double> Mu;
    Mu.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            Mu[i]+=M[j]*u[node];
        }
    }
    
    // -------------- conjugate gradient method -----------------
    int iter;
    double alpha,beta;
    double c1,c2,c3;
    double err;
    
    vector <double> r;
    vector <double> p;
    vector <double> Mp;

    r.assign(n_p,0.0); // 残差ベクトル
    p.assign(n_p,0.0); // 修正ベクトル
    Mp.assign(n_p,0.0);
    
    for(i=0;i<n_p;i++){u[i]=0.0;} // 初期化
    for(i=0;i<n_p;i++){
        r[i]=b[i]-Mu[i];
        p[i]=r[i];
    }
    
    // ---------- iterative calculation start -------------
    for(iter=1;iter<Iter_Max;iter++){
        //境界条件のケア
        for(i=0;i<n_p;i++){
            if(set_bc_list[i]==0){
                p[i]=0.0;
            }
        }
        
        for(i=0;i<n_p;i++){
            Mp[i]=0.0;
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                Mp[i]+=M[j]*p[node];
            }
        }

        c1=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        c2=dot_product(p,Mp,tmp_vec,d_p_terms,n_p_Lg2);
        alpha=c1/c2;
        
        for(i=0;i<n_p;i++){
            u[i]+=alpha*p[i];
            r[i]-=alpha*Mp[i];
        }
        
        //　境界条件のケア
        for(i=0;i<n_p;i++){if(set_bc_list[i]==0){r[i]=0.0;}}
        
        err=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        err/=n_p;
        
        if(iter%50==0&&iter<1000){cout<<"error:"<<err<<endl;} // confirm an error.
        if(EPS > err) break;

        c3=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2);
        beta=c3/c1;
        // if(iter%50==0&&iter<1000){cout<<"beta:"<<beta<<endl;}

        for(i=0; i<n_p; i++){p[i]=r[i]+beta*p[i];}
    } // end of the for-loop of iter
    // ---------- iterative calculation end -------------
    
    //境界条件の再代入
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=0.0;
        }
    }
}