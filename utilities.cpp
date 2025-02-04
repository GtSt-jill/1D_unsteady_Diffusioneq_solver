//数値誤差を回避できる内積計算のための関数
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

// * 疎行列処理した行列とベクトルの掛け算のための関数
void matmul(vector <double> X, vector <double> v, vector <double> &ans){
    int i, j;
    int node;

    #pragma omp parallel for private(j, node)
    for(i=0; i<n_p; i++){
        ans[i] = 0.0;
        for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
            node = n_link[j];
            ans[i] += X[j] * v[node];
        }
    }
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
    
    // KAの対角項と非対角項の生成 要確認
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
        
        error=0.0; // ATTENTION : Initializing an error
        for(i=0;i<n_p;i++){
            u[i]=b_tmp[i]/KA_d[i]; // update a solution
            error+=sqrt((u[i]-u_tmp[i])*(u[i]-u_tmp[i])); // calculate an error
        }
        error=error/n_p;
        
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
    vector <double> MS;
    MS.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        MS[i]=(M[i]+M_c[i])*2+dt*(A[i]+A_c[i]+K[i]);
    }
    
    // 境界条件による生まれる外力項の計算
    vector <double> MSu_b;
    MSu_b.assign(n_p,0.0);
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            MSu_b[i]+=MS[j]*u[node];
        }
    }
    
    // post-processing considering B.C.
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(node==i){MS[j]=1.0;}
                else{MS[j]=0.0;}
            }
        }
        else{
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                if(set_bc_list[node]==0){
                    MS[j]=0.0;
                }
            }
        }
    }
    
    // KAの対角項と非対角項の生成 要確認
    vector <double> MS_d; // 対角項
    vector <double> MS_t; // 非対角項
    MS_d.assign(n_p,0.0);
    MS_t.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        MS_t[i]=MS[i];
    }
    
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            if(i==node){
                MS_d[i]=MS[j];
                MS_t[j]=MS[j]-MS_d[i];
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
            b_tmp[i]=b[i]; // 足し込みしてる場合は初期化忘れず！
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                b_tmp[i]-=MS_t[j]*u[node];
            }
        }
        
        for(i=0;i<n_p;i++){
            b_tmp[i]-=MSu_b[i];
        }
        
        error=0.0; // ATTENTION : Initializing an error
        for(i=0;i<n_p;i++){
            u[i]=b_tmp[i]/MS_d[i]; // update a solution
            error+=(u[i]-u_tmp[i])*(u[i]-u_tmp[i]); // calculate an error
        }
        error = sqrt(error);
        error = error/n_p;
        // if(iter==1) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        // if(iter%500==0) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        if(error<EPS){break;} // Convergence Judgement
        if(iter==Iter_Max-1) cout<<error<<endl;
    }
    // printf("iter : %d\n", iter);
    // --------------- iterative calculation ends ----------------
    
    //--------- 境界条件の再代入 ----------
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=0.0;
        }
    }
}

// Mx=bを解く共役勾配法 動的解析に用いる
void BiCG_dynamics(vector <double> b){
    int i,j;
    int node;
    
    vector <double> MS;
    MS.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        MS[i]=(M[i]+M_c[i])*2+dt*(A[i]+A_c[i]+K[i]);
    }

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
    vector <double> MSu;
    MSu.assign(n_p,0.0);

    matmul(MS, u, MSu);
    
    // -------------- conjugate gradient method -----------------
    int iter;
    double alpha,beta;
    double c1,c2,c3;
    double err;
    
    vector <double> r;
    vector <double> p;
    vector <double> MSp;
    vector <double> r_bi;
    vector <double> p_bi;
    vector <double> MSp_bi;

    r.assign(n_p,0.0); // 残差ベクトル
    p.assign(n_p,0.0); // 修正ベクトル
    MSp.assign(n_p,0.0);
    r_bi.assign(n_p,0.0);
    p_bi.assign(n_p,0.0);
    MSp_bi.assign(n_p,0.0);
    
    for(i=0;i<n_p;i++){u[i]=0.0;} // 初期化

    for(i=0;i<n_p;i++){
        r[i] = b[i]-MSu[i];
        p[i] = b[i]-MSu[i];
        r_bi[i] = r[i];
        p_bi[i] = p[i];
    }
    
    // ---------- iterative calculation start -------------
    for(iter=1;iter<Iter_Max;iter++){
        //境界条件のケア
        for(i=0;i<n_p;i++){
            if(set_bc_list[i]==0){
                p[i]=0.0;
                p_bi[i] = 0.0;
            }
        }

        matmul(MS, p, MSp);

        for(i=0; i<n_p; i++) MSp_bi[i] = 0.0;

        for(i=0;i<n_p;i++){
            for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
                node = n_link[j];
                MSp_bi[node] += MS[j]*p_bi[i];
            }
        }

        // 境界条件のケア
        for(i=0;i<n_p;i++){
            if(set_bc_list[i]==0){
                r[i] = 0.0;
                r_bi[i] = 0.0;
            }
        }

        c1=dot_product(r,r_bi,tmp_vec,d_p_terms,n_p_Lg2);
        c2=dot_product(p_bi,MSp,tmp_vec,d_p_terms,n_p_Lg2);

        alpha = c1/c2;
        
        for(i=0;i<n_p;i++){
            u[i] += alpha*p[i];
            r[i] -= alpha*MSp[i];
            r_bi[i] -= alpha*MSp_bi[i];
        }
        
        // 境界条件のケア
        for(i=0;i<n_p;i++){
            if(set_bc_list[i]==0){
                r[i] = 0.0;
                r_bi[i] = 0.0;
            }
        }
        
        err = dot_product(r,r_bi,tmp_vec,d_p_terms,n_p_Lg2);
        err /= n_p;
        
        // if(iter%50==0&&iter<1000){cout<<"error:"<<err<<endl;} // confirm an error.
        // if(iter%50 == 0) printf("iter : %d\terror : %f\n", iter, err);
        if(EPS > err) break;

        c3 = dot_product(r,r_bi,tmp_vec,d_p_terms,n_p_Lg2);
        beta = c3/c1;

        for(i=0; i<n_p; i++){
            p[i] = r[i] + beta*p[i];
            p_bi[i] = r_bi[i] + beta*p_bi[i];
        }
    } // end of the for-loop of iter
    // printf("iter : %d\n", iter);
    // ---------- iterative calculation end -------------
    
    //境界条件の再代入
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
    const int n_o_f = 100; // ? number of frames
    
    // Intitial condtions of u
    for(i=0;i<n_p;i++){
        if(x[i]<=L/4) u[i]=4*x[i];
        else if(x[i]>L/4&&x[i]<=L/2) u[i]=-4*x[i]+2*L;
        else if(x[i]>L/2) u[i]=0.0;
    }
    
    // 毎ステップで解ベクトルに乗ずる行列．時間によらない定数のためnのfor-loop文の外で定義．
    vector <double> M2dtS;
    M2dtS.assign(n_link.size(),0.0);
    for(i=0;i<n_link.size();i++){
        M2dtS[i]=(M[i]+M_c[i])*2-dt*(A[i]+A_c[i]+K[i]);
    }
    
    vector <double> M2dtSu;
    M2dtSu.assign(n_p,0.0);
    
    // file exporting
    
    char fname[100];// 配列の大きさに注意！！ zsh: abort というエラー出る．
    sprintf(fname,"animation/u_dynamics_ad=%.2f_Di=%.2f.dat",ad,Di);
    ofstream out(fname,ios::app); // 上書き設定
    
    // File-initializing ios::appによるファイル上書きのため初期化が必要．
    ofstream out_tmp(fname);
    out_tmp.close();
    
    //!--------------------- dynamical analysis start -------------------
    for(n = 1; n < max_step; n++){
        if(n%50==0) printf("step number : %d\n",n);
        
        // File-exporting
        if(n % int(max_step / n_o_f) == 1){
            for(i = 0; i < n_p; i++){
                out << x[i] << " " << u[i] << endl;
            }
            // paraview_visualize(n); // paraviewへの出力
        }

        // 中央差分法  差分法の結果と等価
        for(i=0;i<n_p;i++){
            M2dtSu[i]=0.0;
            for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
                node = n_link[j];
                M2dtSu[i] += M2dtS[j]*u[node];
            }
        }
        
        // perform Jacobi scheme
        Jacobi_dynamics(M2dtSu);
        // BiCG_dynamics(M2dtSu);
        
        
    }
    out.close();
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