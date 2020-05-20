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
    
    // KAの対角項と非対角項の生成　要確認
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
            error+=sqrt((u[i]-u_tmp[i])*(u[i]-u_tmp[i])); // calculate an error
        }
        
        error=error/n_p;
        // if(iter==1) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        // if(iter%500==0) cout<<"iter:"<<iter<<" error:"<<error<<endl;
        if(error<EPS){break;} // Convergence Judgement
        if(iter==Iter_Max-1) cout<<error<<endl;
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
    const int n_o_f = 100;
    
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
    
    //--------------- dynamical analysis start ----------------
    char fname[100];// 配列の大きさに注意！！ zsh: abort というエラー出る．
    sprintf(fname,"array/u_dynamics_ad=%.2f_Di=%.2f.dat",ad,Di);
    ofstream out(fname,ios::app);
    
    // File-initializing ios::appによるファイル上書きのため初期化が必要．
    ofstream out_tmp(fname);
    out_tmp.close();
    
    for(n = 1; n < max_step; n++){
        if(n%50==0) printf("step number : %d\n",n);
        
        // 中央差分法  差分法の結果と等価
        for(i=0;i<n_p;i++){
            M2dtSu[i]=0.0;
            for(j=n_link_start[i];j<n_link_start[i+1];j++){
                node=n_link[j];
                M2dtSu[i]+=M2dtS[j]*u[node];
            }
        }
        
        // perform Jacobi method
        Jacobi_dynamics(M2dtSu);
        
        // File-exporting
        if(n % int(max_step / n_o_f) == 1){
            for(i = 0; i < n_p; i++){
                out << x[i] << " " << u[i] << endl;
            }
            // paraview_visualize(n); // paraviewへの出力
        }
    }// end of the for-loop of n
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
        
        c1=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        c2=dot_product(p,KAtKAp,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        
        alpha=c1/c2;
        
        for(i=0;i<n_p;i++){
            u[i]+=alpha*p[i];
            r[i]-=alpha*KAtKAp[i];
        }
        
        // 境界条件のケア
        for(i=0;i<n_p;i++){if(set_bc_list[i]==0){r[i]=0.0;}}
        
        err=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        err/=n_p;
        
        if(iter==1) printf("LOOP : %d    %.10e\n", iter, err);
        if(iter%50==0) printf("LOOP : %d    %.10e\n", iter, err);
        if(EPS > err) break;
        
        c3=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, n_p);
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

        c1=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2,n_p);
        c2=dot_product(p,Mp,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        alpha=c1/c2;
        
        for(i=0;i<n_p;i++){
            u[i]+=alpha*p[i];
            r[i]-=alpha*Mp[i];
        }
        
        // 境界条件のケア
        for(i=0;i<n_p;i++){if(set_bc_list[i]==0){r[i]=0.0;}}
        
        err=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        err/=n_p;
        
        if(iter%50==0&&iter<1000){cout<<"error:"<<err<<endl;} // confirm an error.
        if(EPS > err) break;

        c3=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, n_p);
        beta=c3/c1;
        // if(iter%50==0&&iter<1000){cout<<"beta:"<<beta<<endl;}

        for(i=0; i<n_p; i++){p[i]=r[i]+beta*p[i];}
    } // end of the for-loop of iter
    // ---------- iterative calculation end -------------
    
    // 境界条件の再代入
    for(i=0;i<n_p;i++){
        if(set_bc_list[i]==0){
            u[i]=0.0;
        }
    }
}

// For SSFEM

void cg_method_SSFEM(vector <double> b){
    int i, j, k, l, m;
    int node;
    int NL = n_link.size();
    int a = d_o_s*(KL+1);
    
    //--------- 内積計算用 ------------- 
    // ! n_pをN_pにするのをわすれずに！！
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのテンポラリーベクター
    n_p_Lg2 = (int)ceil(log((double)N_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2 
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);
    //--------------------------------
    
    //境界条件の設定 境界では確率的な揺らぎは0に設定している．
    for(i=0;i<N_p;i++){u[i]=0.0;}
    for(i=0;i<d_o_s;i++){
        for(j=0; j<n_p; j++){
            if(set_bc_list[j]==0) {u[n_p*i+j]=0.0;}
        }
    }

    vector <double> MU;
    MU.assign(N_p, 0.0);
    // ----------------------------MUの計算 start--------------------------
    for(i=0; i<d_o_s; i++){
            for(j=0; j<d_o_s; j++){
                for(k=0; k<n_p; k++){
                    MU[n_p*i+k] = 0.0;
                    for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                        node = n_link[l];
                        MU[n_p*i+k] += PPx[a*i + d_o_s*j + 0] * M[l] * u[n_p*j+node];
                    }
                    
                }
            }
        }

    for(i=0; i<d_o_s; i++){
        for(j=0; j<d_o_s; j++){
            for(k=0; k<n_p; k++){
                for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                    node = n_link[l];
                    for(m=0; m<KL+1; m++){
                        MU[n_p*i+k] += PPx[a*i + d_o_s*j + m] * K[NL*m+l] * u[n_p*j+node] * dt / 2;
                    }
                }
            }
        }
    }
    // ----------------------------MUの計算 end--------------------------

    // Setting for CG method.
    int iter;
    double alpha,beta;
    double c1,c2,c3;
    double err;
    
    vector <double> r;
    vector <double> p;
    vector <double> Mp;

    r.assign(N_p,0.0); // 残差ベクトル
    p.assign(N_p,0.0); // 修正ベクトル
    Mp.assign(N_p,0.0);
    
    for(i=0;i<N_p;i++){u[i]=0.0;} // 初期化
    for(i=0;i<N_p;i++){
        r[i]=b[i]-MU[i];
        p[i]=r[i];
    }

    // ---------- iterative calculation of CG start -------------
    for(iter=1;iter<Iter_Max;iter++){
        //境界条件のケア
        for(i=0; i<d_o_s; i++){
            for(j=0; j<n_p; j++){
                if(set_bc_list[j]==0){
                    p[n_p*i+j]=0.0;
                }
            }
        }
        
        // ----------------------------Mpの計算 start--------------------------
        for(i=0; i<d_o_s; i++){
            for(j=0; j<d_o_s; j++){
                for(k=0; k<n_p; k++){
                    Mp[n_p*i+k] = 0.0;
                    for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                        node = n_link[l];
                        Mp[n_p*i+k] += PPx[a*i + d_o_s*j + 0] * M[l] * p[n_p*j+node];
                    }
                    
                }
            }
        }

        for(i=0; i<d_o_s; i++){
            for(j=0; j<d_o_s; j++){
                for(k=0; k<n_p; k++){
                    for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                        node = n_link[l];
                        for(m=0; m<KL+1; m++){
                            Mp[n_p*i+k] += PPx[a*i + d_o_s*j + m] * K[NL*m+l] * p[n_p*j+node] * dt / 2;
                        }
                    }
                }
            }
        }
        // ----------------------------MUの計算 end--------------------------

        c1 = dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2,N_p);
        c2 = dot_product(p,Mp,tmp_vec,d_p_terms,n_p_Lg2, N_p);
        alpha = c1/c2;
        
        for(i=0;i<N_p;i++){
            u[i] += alpha*p[i];
            r[i] -= alpha*Mp[i];
        }
        
        // 境界条件のケア
        for(i=0; i<d_o_s; i++){
            for(j=0; j<n_p; j++){
                if(set_bc_list[j]==0){r[n_p*i+j]=0.0;}
            }
        }
        
        err = dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, N_p);
        err /= N_p;
        
        if(iter%50==0&&iter<1000){cout<<"error:"<<err<<endl;} // confirm an error.
        if(EPS > err) break;

        c3 = dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2, N_p);
        beta = c3/c1;
        // if(iter%50==0&&iter<1000){cout<<"beta:"<<beta<<endl;}

        for(i=0; i<N_p; i++){p[i]=r[i]+beta*p[i];}
    } // end of the for-loop of iter
    // ---------- iterative calculation end -------------
    
    // 境界条件の再代入
    for(i=0; i<d_o_s; i++){
        for(j=0; j<n_p; j++){
            if(set_bc_list[j]==0) u[n_p*i+j]=0.0;
        }
    }
}

void Dynamics_SSFEM(){
    int i, j, k, l, m;
    int node;
    int n;
    int a = (KL+1)*d_o_s;
    int NL = n_link.size();

    // Setting an initial condition.
    for(i=0; i<n_p; i++){
        double x_tmp,y_tmp;
        x_tmp=x[i];
        if(x_tmp>0.5&&x_tmp<0.8){
            u[i]=x_tmp*x_tmp;
        }else{
            u[i]=0.0;
        }
    }
    vector <double> b;
    b.assign(N_p, 0.0);

    for(n = 1; n<3; n++){
        printf("step_number : %d\n", n);
        // MB*Uの計算
        for(i=0; i<d_o_s; i++){
            for(j=0; j<d_o_s; j++){
                for(k=0; k<n_p; k++){
                    b[n_p*i+k] = 0.0;
                    for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                        node = n_link[l];
                        b[n_p*i+k] += PPx[a*i + d_o_s*j + 0] * M[l] * u[n_p*j+node];
                    }
                }
            }
        }

        // K*dt/2 * Uの計算
        for(i=0; i<d_o_s; i++){
            for(j=0; j<d_o_s; j++){
                for(k=0; k<n_p; k++){
                    for(l=n_link_start[k]; l<n_link_start[k+1]; l++){
                        node = n_link[l];
                        for(m=0; m<KL+1; m++){
                            b[n_p*i+k] -= PPx[a*i + d_o_s*j + m] * K[NL*m+l] * u[n_p*j+node] * dt / 2;
                        }
                    }
                }
            }
        }

        // conjugate gradient method
        cg_method_SSFEM(b);
    }
}