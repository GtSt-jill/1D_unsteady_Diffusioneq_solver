// create data sets for paraview to visualize results
void paraview_visualize(int step_number){
    int i;
    char fname[50];
    // FILE *data;
    sprintf(fname,"u_dynamics_data/tmp_%d.vtk",step_number);
    // data=fopen(fname,"w");
    // fprintf(data,"# vtk DataFile Version 3.0\ntest\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS\t%d\tfloat\n",n_p);
    // for(i=0;i<n_p;i++){
    //     fprintf(data,"%f\t%f\t%f\n",x[i],0.0,0.0);// \tは水平タブ
    // }
    // fclose(data);
    
    ofstream out(fname);

    out<<"# vtk DataFile Version 3.0"<<endl;
    out<<"test"<<endl<<"ASCII"<<endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<endl;
    out<<"POINTS  "<<n_p<<" float"<<endl;
    for(i=0;i<n_p;i++){
        out<<x[i]<<" "<<0.0<<" "<<0.0<<endl;
    }
    out<<"CELLS "<<n_e<<" "<<n_e*(2+1)<<endl;
    for(i=0;i<n_e;i++){
        out<<2<<" "<<el[2*i]<<" "<<el[2*i+1]<<endl;
    }

    out<<"CELL_TYPES "<<n_e<<endl;
    for(i=0;i<n_e;i++){
        out<<2<<endl;
    }

    out<<"POINT_DATA "<<n_p<<endl<<"SCALARS velocity float 1"<<endl;
    out<<"LOOKUP_TABLE default"<<endl;
    for(i=0;i<n_p;i++){
        out<<u[i]<<endl;
    }
    out.close();
}

void export_data(){
    FILE *data;
    char fname[20];
    sprintf(fname,"test444.dat");
    int i;
    data=fopen(fname,"w");
    fprintf(data,"this is a test\n");
    fclose(data);
}
