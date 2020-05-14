void set_bc(){
    // int i;
    set_bc_list[0]=0;
    set_bc_list[n_p-1]=0;
    printf("境界条件設定したよ\n");
}
/*
void set_bc(){
    int i,j,k,l;
    int node;
    int e0,e1;
    int count=0;
    
    for(i=0;i<n_p;i++){
        for(j=n_link_start[i];j<n_link_start[i+1];j++){
            node=n_link[j];
            
            for(k=e_link_start[i];k<e_link_start[i+1];k++){
                e0=e_link[k];
                
                for(l=e_link_start[node];l<e_link_start[node+1];l++){
                    e1=e_link[l];
                    if(e0==e1) count++;
                }
            }
            if(count==1){
                set_bc_list[i]=0;
                set_bc_list[node]=0;
            }
            count=0;
        }
    }
    //確認用
//    for(i=0;i<n_p;i++){
//        cout<<set_bc_list[i]<<endl;
//    }
    cout<<"境界節点のリスト作ったよ"<<endl;
}
*/
