void Assign_array(){
    x.assign(n_p,0.0);
    el.assign(n_e*2,0);
    set_bc_list.assign(n_p*d_o_s,1);
    u.assign(n_p*d_o_s,0.0);
}

// e_linkとn_linkの構成
void InitialAdjacency(){
    int i, j, jstart, jend;
    int itmp;
    int e0, e1; //線要素
    
    vector <int>::iterator p;
    vector <int>::iterator p_start, p_end;
    vector <int> n_link_tmp;
    
//    --- Link list for element -------- start ------------
    e_link_start.assign(n_p+1,0);
    
    for(i=0; i<n_e; i++){
        itmp=2*i;
        e_link_start[el[itmp]]++;
        e_link_start[el[itmp+1]]++;
    } // Number of elements linked to each particle
    
    for(i=1; i<=n_p; i++) e_link_start[i]+=e_link_start[i-1];
    
    for(i=n_p; i>0; i--) e_link_start[i]=e_link_start[i-1];
    
    e_link_start[0]=0;
    
    e_link.assign(e_link_start[n_p], -1);
    
    for(i=0; i<n_e; i++){
        itmp = 2*i;
        e0 = el[itmp];
        e1 = el[itmp+1];

        jstart = e_link_start[e0];
        jend = e_link_start[e0+1];
        for(j=jstart; j<jend; j++){
            if(e_link[j]==-1){
                e_link[j]=i;
                j=jend; // break;
            }
        }
        
        jstart=e_link_start[e1];
        jend=e_link_start[e1+1];
        for(j=jstart; j<jend; j++){
            if(e_link[j]==-1){
                e_link[j]=i;
                j=jend;
            }
        }
    }
    
//    --- Link list for node -------- start ------------
    n_link_start.assign(n_p+1, 0);  // Initialize  n_link_start[]
    n_link.reserve(e_link.size());
    
    for(i=0; i<n_p; i++){
        jstart=e_link_start[i];
        jend=e_link_start[i+1];
        for(j=jstart; j<jend; j++){
            itmp=2*e_link[j];
            n_link_tmp.push_back(el[itmp]);
            n_link_tmp.push_back(el[itmp+1]);
        }
        
        sort(n_link_tmp.begin(),n_link_tmp.end()); // sort．
        p_end=unique(n_link_tmp.begin(),n_link_tmp.end()); // 連続する要素をまとめる．

        p_start=n_link_tmp.begin();// この3行を使えば，1. 自分自身，2. 対称性によるredundantな節点，のリンクも含まれる．
        for(p=p_start; p<p_end; p++) n_link.push_back(*p);
        n_link_start[i+1]=n_link_start[i]+(p_end-p_start);



        n_link_tmp.clear(); //毎回，n_link_tmpはクリアしないといけない
    }
    
    n_link.shrink_to_fit();
}
