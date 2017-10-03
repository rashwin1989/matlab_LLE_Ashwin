
  // input data
  MY_PRINTF("\n--------------------------------------------------\n");
  MY_PRINTF(  "-----Input data: input.dat+-----------------------\n");
  MY_PRINTF("n_species: %d\n", n);
  MY_PRINTF("species:   ");
  for (i=0;i<n;i++) MY_PRINTF("%d ", idx_species[i]); MY_PRINTF("\n");
  MY_PRINTF("groups:   ");
  for (i=0;i<n;i++) MY_PRINTF("%d ", idx_groups[i]);  MY_PRINTF("\n");

  MY_PRINTF("y1:       ");
  for (i=0;i<n;i++) MY_PRINTF("%6.4lf ", Y01[i]);  MY_PRINTF("\n");
  MY_PRINTF("y2:       ");
  for (i=0;i<n;i++) MY_PRINTF("%6.4lf ", Y02[i]);  MY_PRINTF("\n");
  MY_PRINTF("x1:       ");
  for (i=0;i<n;i++) MY_PRINTF("%6.4lf ", x01[i]);  MY_PRINTF("\n");
  MY_PRINTF("x2:       ");
  for (i=0;i<n;i++) MY_PRINTF("%6.4lf ", x02[i]);  MY_PRINTF("\n");
  
  MY_PRINTF("T0:     %lf\n", T0);
  MY_PRINTF("T1:     %lf\n", T1);
  MY_PRINTF("P:      %le\n", P );
  MY_PRINTF("m12:    %lf\n", m012);
  MY_PRINTF("c1:    %lf\n", c01);
  
  MY_PRINTF( "---------------------------------------------------\n");

  // species data
  MY_PRINTF("\n--------------------------------------------------\n");
  MY_PRINTF(  "-----Species' properties--------------------------\n");
  for (i=0; i<n; i++){
    MY_PRINTF("%5.1f,%6.1f,%4.1f,%6.1f,%6.1f,%4.1f,%5.1f,%5.3f,%4.2f,%6.1f\n",
        Pc[i]*1e-5,Tc[i],w[i],MW[i],Tb[i],SG[i],H8[i],k[i],dm[i],Vc[i]);
  }
  MY_PRINTF( "---------------------------------------------------\n");

