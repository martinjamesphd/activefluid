#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<hdf5.h>
#include<math.h>
#include"ns2d.h"

void store_spectrum(double *e_spectra, hid_t spectra_file_id, int k)
{
 int y_size_f = sys_size/2+1;
 char dset[10];
 
 hid_t       dataset_id, dataspace_id;
 
 sprintf(dset,"%d",k);
 
 hsize_t     dims[1];

 /* Create the data space for the dataset. */
 dims[0] = y_size_f;
 dataspace_id = H5Screate_simple(1, dims, NULL);

 /* Create the dataset. */
 dataset_id = H5Dcreate2(spectra_file_id, dset, H5T_IEEE_F64BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                          
 	    
 H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, e_spectra);
		      
 H5Dclose(dataset_id);

  
 H5Sclose(dataspace_id);
}

hid_t init_spectra_storage()
{

 int y_size_f = sys_size/2+1,i;
 char dset[]="k";
 
 double *k_vect=(double*)malloc(y_size_f*sizeof(double));
 
 hid_t file_id;
 
 hid_t       dataset_id, dataspace_id;
 
 file_id = H5Fcreate("spectra.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
  
 hsize_t     dims[1];

 /* Create the data space for the dataset. */
 dims[0] = y_size_f;
 dataspace_id = H5Screate_simple(1, dims, NULL);

 /* Create the dataset. */
 dataset_id = H5Dcreate2(file_id, dset, H5T_IEEE_F64BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
 
 for(i=0;i<y_size_f;i++)
  k_vect[i]=scale*i;                        
 	    
 H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, k_vect);
		      
 H5Dclose(dataset_id);

  
 H5Sclose(dataspace_id);
 
 return file_id;
}

void store_field(double *omega, hid_t field_file_id, int k)
{
 int size = sys_size*(sys_size+2);
 char dset[10];
 
 hid_t       dataset_id, dataspace_id;
 
 sprintf(dset,"%d",k);
 
 hsize_t     dims[1];

 /* Create the data space for the dataset. */
 dims[0] = size;
 dataspace_id = H5Screate_simple(1, dims, NULL);

 /* Create the dataset. */
 dataset_id = H5Dcreate2(field_file_id, dset, H5T_IEEE_F64BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                          
 	    
 H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, omega);
		      
 H5Dclose(dataset_id);

  
 H5Sclose(dataspace_id);
}

hid_t init_field_storage()
{

 char dsetx[]="parameters";
 char dsety[]="parameter_values";
 
 char parameters[6]={"NslabL"};

 double parameter_values[6]={max_iter,sys_size,nlin,mu-1,mu2,2*M_PI/scale};

 parameter_values[0]=max_iter;
 
 hid_t file_id;
 
 hid_t       dataset_id, dataspace_id;
 
 file_id = H5Fcreate("field.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
  
 hsize_t     dims[1];

 /* Create the data space for the dataset. */
 dims[0] = 6;
 dataspace_id = H5Screate_simple(1, dims, NULL);

 /* Create the dataset. */
 dataset_id = H5Dcreate2(file_id, dsetx, H5T_C_S1, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
 
 H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, parameters);
 
 H5Dclose(dataset_id);

  
 H5Sclose(dataspace_id);
 
 dataspace_id = H5Screate_simple(1, dims, NULL);

 /* Create the dataset. */
 dataset_id = H5Dcreate2(file_id, dsety, H5T_IEEE_F64BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
 
 H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, parameter_values);
 
 H5Dclose(dataset_id);

  
 H5Sclose(dataspace_id);
 
 return file_id;
}
