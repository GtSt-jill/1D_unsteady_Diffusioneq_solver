//数値誤差を回避できる内積計算のための関数
double dot_product(vector <double> &v, vector <double> &w, vector <double> &tmp, int n, int lg2n, int N)
{
    int i,ii;
    for(i=0; i<N; i++)
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

double Covariance_Function(double x1, double x2){
    double d = abs(x1 - x2);
    return c_o_v * c_o_v * exp(- d / stan_len);
}

void Calculate_Covmatrix(){
    int     i, j, k, l;
    int     e0, e1; // element
    double  x0g, x1g; // gravity point
    double  l0, l1; // element length

    Cov.assign(n_p*n_p, 0.0);

    for(i=0; i<n_p; i++){
        for(j=0; j<n_p; j++){
            for(k=e_link_start[i]; k<e_link_start[i+1]; k++){
               e0  = e_link[k];
               x0g = (x[el[2*e0]] + x[el[2*e0+1]]) / 2;
               l0  = x[el[2*e0+1]] - x[el[2*e0]];
               for(l=e_link_start[j]; l<e_link_start[j+1]; l++){
                   e1  = e_link[l];
                   x1g = (x[el[2*e1]] + x[el[2*e1+1]]) / 2;
                   l1  = x[el[2*e1+1]] - x[el[2*e1]];
                   Cov[n_p*i+j] += Covariance_Function(x0g, x1g)*l0*l1;
               }
            }
        }
    }
    printf("Covを構成したよ．\n");

    confirm_array_double(Cov, n_p*n_p, "Cov_mat");
}

// AAq = bを解くプログラム
void eigen_cgmethod(vector <double> AA, vector <double> &q, vector <double> b){
    int i, j;
    int iter;
    int node;
    double c1, c2, c3;
    double alpha, beta;
    double error;

    vector <double> r;
    vector <double> p;
    vector <double> AAq;
    vector <double> AAp;

    r.assign(n_p,0.0);
    p.assign(n_p,0.0);
    AAq.assign(n_p,0.0);
    AAp.assign(n_p,0.0);

    // 内積計算用
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのテンポラリーベクター
    n_p_Lg2 = (int)ceil(log((double)n_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);

    // 初期ベクトルの設定 
    for(i=0; i<n_p; i++) {q[i] = 0;}

    // AA*qの計算
    for(i=0; i<n_p; i++){
        for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
            node = n_link[j];
            AAq[i] = AA[j] * q[node];
        }
    }

    // initializeing
    for(i=0;i<n_p;i++){
        r[i]=b[i]-AAq[i];
        p[i]=b[i]-AAq[i];
    }

    //  ------------------iteration start-------------------
    for(iter=1; iter<Iter_Max; iter++){
        // AA*pの計算
        for(i=0; i<n_p; i++){
            AAp[i] = 0.0;
            for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
                node = n_link[j];
                AAp[i] += AA[j] * p[node];
            }
        }
        
        c1    = dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2,n_p);
        c2    = dot_product(p,AAp,tmp_vec,d_p_terms,n_p_Lg2,n_p);
        alpha = c1/c2;
        
        for(i=0;i<n_p;i++){
            q[i] = q[i]+alpha*p[i];
            r[i] = r[i]-alpha*AAp[i];
        }
        
        error = sqrt(dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2,n_p));
        if(error < EPS / n_p) break;
        // cout << error << endl;
        c3=dot_product(r,r,tmp_vec,d_p_terms,n_p_Lg2,n_p);
        beta = c3/c1;

        for(i=0;i<n_p;i++) {p[i] = r[i] + beta*p[i];}
    }
}

void Power_method(){
    int i,j,k,l;
    int node;
    double err;
    int iter; // 冪乗法の繰り返し回数
    int iter_max = 1000;
    double eigen_norm;
    double eigenvalue_tmp;
    double ev1ev2;

    // 内積計算用
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのテンポラリーベクター
    n_p_Lg2 = (int)ceil(log((double)n_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);
    
    vector <double> Cov_ev2(n_p);//Cov*eigenvector2
    vector <double> M_ev1(n_p);//B*eigenvector1
    vector <double> Cov_ev1(n_p);//y=Cov*eigen_vec
    vector <double> M_inv_Cov_ev(n_p);//z=B^(-1)*Cov*eigen_vec
    vector <double> eigenvector_tmp;

    eigenvector_tmp.assign(n_p,0.0);

    // KL展開の次数分だけ固有値を求める．
    for(k = 0; k < KL; k++){
        for(i=0;i<n_p;i++){eigenvector_tmp[i]=eigenvectors[k*n_p+i];}
        eigenvector_tmp[0]=1.0;

        for(iter=1; iter<iter_max; iter++){
            if(k == 0){
                for(i=0; i<n_p; i++){
                    Cov_ev1[i] = 0.0; // ! 発散してしまうので注意
                    for(j=0; j<n_p; j++){
                        Cov_ev1[i] += Cov[n_p*i+j] * eigenvector_tmp[j];
                    }
                }
            }else{
                // Cov*ev2の計算
                for(i=0;i<n_p;i++){
                    Cov_ev2[i]=0;
                    for(j=0; j<n_p; j++){
                        Cov_ev2[i] += Cov[n_p*i+j]*eigenvector_tmp[j];
                    }
                }
                
                vector <double> pre_eigenvector(n_p);
                for(l=0;l<k;l++){
                    for(i=0;i<n_p;i++){pre_eigenvector[i]=eigenvectors[l*n_p+i];}
                        
                        ev1ev2=dot_product(pre_eigenvector,eigenvector_tmp,tmp_vec,d_p_terms,n_p_Lg2,n_p);
                        
                        for(i=0; i<n_p; i++){
                            M_ev1[i] = 0.0;
                            for(j=n_link_start[i]; j<n_link_start[i+1]; j++){
                                node = n_link[j];
                                M_ev1[i] += M[j] * pre_eigenvector[node];
                            }
                        }

                        for(i=0;i<n_p;i++) {Cov_ev2[i]-=eigenvalues[l]*ev1ev2*M_ev1[i];}
                    }
                for(i=0;i<n_p;i++){Cov_ev1[i]=Cov_ev2[i];}
            }

            eigen_cgmethod(M, M_inv_Cov_ev, Cov_ev1);

            eigenvalue_tmp = eigenvalues[k];

            eigen_norm = sqrt(dot_product(M_inv_Cov_ev, M_inv_Cov_ev, tmp_vec, d_p_terms, n_p_Lg2, n_p));

            eigenvalues[k] = eigen_norm*eigen_norm / dot_product(eigenvector_tmp,M_inv_Cov_ev,tmp_vec,d_p_terms,n_p_Lg2,n_p);

            err = abs(eigenvalues[k]-eigenvalue_tmp); // / abs(eigenvalues[k]); // abs(eigen_value1);

            if(err<EPS) break;
            for(i=0;i<n_p;i++){eigenvector_tmp[i]=M_inv_Cov_ev[i]/eigen_norm;}
        }
        printf("Iteration : %d \teigenvalue %d : %f\n", iter, k+1, eigenvalues[k]);
        // cout<<"iteration:"<<iter<<" "<<"eigenvalue"<<k+1<<":"<<eigenvalues[k]<<endl;
        for(i=0;i<n_p;i++){eigenvectors[k*n_p+i]=eigenvector_tmp[i];}
        confirm_array_double(eigenvectors, n_p*KL, "eigenvectors");
    }

}