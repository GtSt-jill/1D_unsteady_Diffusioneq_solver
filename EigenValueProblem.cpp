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

double Covariance_Function(double x1, double x2){
    double d = abs(x1 - x2);
    return c_o_v * exp(d / stan_len);
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

void eigen_cgmethod(){

}

void Solve_EigenvalueProblem(){
    int i,j,k,l;
    double err;
    double z_norm;
    double eigenvalue_tmp;
    int ed;
    int e0,e1;
    int iter; // 冪乗法の繰り返し回数
    int iter_max = 100;

    // 内積計算用
    int n_p_Lg2, d_p_terms;
    vector <double> tmp_vec; // 内積計算での情報落ちを防ぐためのテンポラリーベクター
    n_p_Lg2 = (int)ceil(log((double)n_p)/M_LN2);    // log((double)n_p) = ln(n_p), M_LN2:ln2
    d_p_terms = 1;
    for(i=0; i<n_p_Lg2; i++){d_p_terms *= 2;}
    tmp_vec.assign(d_p_terms,0.0);
    
    double ev1ev2;
    vector <double> Cov_ev2(n_p);//Cov*eigenvector2
    vector <double> B_ev1(n_p);//B*eigenvector1
    vector <double> y(n_p);//y=Cov*eigen_vec
    vector <double> z(n_p);//z=B^(-1)*Cov*eigen_vec
    vector <double> eigenvector_tmp;
    eigenvector_tmp.assign(n_p,0.0);

    // KL展開の次数分だけ固有値を求める．
    for(k = 0; k < KL; k++){
        
    }

}