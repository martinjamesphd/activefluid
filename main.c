/*This program solves active fluids model in 2 dimensions
* for a periodic velocity using a pseudo-
* spectral algorithm*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fftw3.h>
#include<hdf5.h>
#include"ns2d.h"

/* Following are the global variables.
  These variables do not change value
  during the execution of program*/
fftw_plan p_for_xvel, p_inv_xvel, p_for_yvel, p_inv_yvel, p_for_w, p_inv_w;
double *exp_dt;			//variable to store time increments and aliasing (refer function init_den_state())
int sys_size, nthreads, nthreads_omp, tau_num,max_iter; 
double dt, vis, vis2, nlin, mu, mu2, scale=1.0; 

//starts function MAIN() 
int main(int argc, char *argv[])
{
	
	int i, j, k, jk, k1, pflag, s_flag;
	double energy, file_return;
	char variable[32];
	FILE *fp, *fp_e;
	
	hid_t spectra_file_id, field_file_id;
	
	fp = fopen("input","r");
	
	file_return=fscanf(fp,"%s",variable);
	
	while((strcmp(variable,"size")==0||strcmp(variable,"N")==0||strcmp(variable,"dt")==0||strcmp(variable,"lambda")==0||strcmp(variable,"alpha")==0||strcmp(variable,"beta")==0||strcmp(variable,"threads")==0||strcmp(variable,"threads_omp")==0||strcmp(variable,"flag")==0||strcmp(variable,"data_points")==0||strcmp(variable,"L")==0)&&file_return!=EOF)
	{
	 
	printf("\n%s\n",variable);
		
	if(strcmp(variable,"size")==0)
	 fscanf(fp,"%d", &sys_size);
	
	if(strcmp(variable,"N")==0)
	 fscanf(fp,"%d", &max_iter);
	
	if(strcmp(variable,"dt")==0)
	 fscanf(fp,"%lE", &dt);
	
	if(strcmp(variable,"lambda")==0)
	 fscanf(fp,"%lE", &nlin);
	
	if(strcmp(variable,"alpha")==0)
	 fscanf(fp,"%lE", &mu);
	
	if(strcmp(variable,"beta")==0)
	 fscanf(fp,"%lE", &mu2);
	
	if(strcmp(variable,"threads")==0)
	 fscanf(fp,"%d", &nthreads);
	
	if(strcmp(variable,"threads_omp")==0)
	 fscanf(fp,"%d", &nthreads_omp);
	
	if(strcmp(variable,"flag")==0)
	 fscanf(fp,"%d",  &pflag);
	
	if(strcmp(variable,"data_points")==0)
	 fscanf(fp,"%d", &s_flag);
	
	if(strcmp(variable,"L")==0)
	 fscanf(fp,"%lE", &scale);
	
	file_return=fscanf(fp,"%s",variable);
	}
	
	fclose(fp);
	
	scale=2*M_PI/scale;
	vis=-2;
	vis2=-1;
	
	mu=mu+1;
	
	int y_size_f = sys_size/2+1;
	int y_size_r = sys_size+2;
	
	double scale2 = 1.0/(sys_size*sys_size);
	
	double enstrophy;

	//Memory allocation
	double *omega 			= (double*)malloc(sys_size*(sys_size+2)*sizeof(double));
	fftw_complex *omega_ft		= (fftw_complex*)omega;
	fftw_complex *omega_t		= (fftw_complex*)malloc((sys_size*y_size_f)*sizeof(fftw_complex));
	double *x_vel 			= (double*)malloc(sys_size*(sys_size+2)*sizeof(double));
	fftw_complex *x_vel_ft		= (fftw_complex*)x_vel;		
	double *y_vel 			= (double*)malloc(sys_size*(sys_size+2)*sizeof(double));
	fftw_complex *y_vel_ft		= (fftw_complex*)y_vel;	
	double *e_spectra		= (double*)malloc(y_size_f*sizeof(double));
	
	exp_dt				= (double*)malloc(sys_size*y_size_f*sizeof(double));

	//fftw variables for forward and inverse transform
	fftw_init_threads();			//initializes threads for fftw
	fftw_plan_with_nthreads(nthreads);
	
	p_for_xvel = fftw_plan_dft_r2c_2d(sys_size, sys_size, x_vel, x_vel_ft, FFTW_ESTIMATE);
        p_inv_xvel = fftw_plan_dft_c2r_2d(sys_size, sys_size, x_vel_ft, x_vel, FFTW_ESTIMATE);
        p_for_yvel = fftw_plan_dft_r2c_2d(sys_size, sys_size, y_vel, y_vel_ft, FFTW_ESTIMATE);
        p_inv_yvel = fftw_plan_dft_c2r_2d(sys_size, sys_size, y_vel_ft, y_vel, FFTW_ESTIMATE);
        p_for_w = fftw_plan_dft_r2c_2d(sys_size, sys_size, omega, omega_ft, FFTW_ESTIMATE);
        p_inv_w = fftw_plan_dft_c2r_2d(sys_size, sys_size, omega_ft, omega, FFTW_ESTIMATE);
	
	double time=0.0;
	int t_flag=max_iter/s_flag;		//flags to store data
		

	
	if(argv[2][0]=='0')				//if no input omega file is given
	{						//initialize omega at t=0 by a random
		init_den_state();
		init_omega(omega_ft);
	}
	
	else						//else read from the file
	{
		fp = fopen(argv[2],"r");
		fread(omega_ft,sys_size*(sys_size+2)*sizeof(double),1,fp);
		init_den_state();
		fclose(fp);
	}
	
	fftw_execute(p_inv_w);	//inverse fourier transforms and fourier transforms
	for(j=0;j<sys_size;j++)							//the initialized omega
		for(k=0;k<y_size_r;k++)
		{	
			jk=j*y_size_r+k;
			omega[jk]*=scale2;
		}
	fftw_execute(p_for_w);
	find_vel_ft( omega_ft, x_vel_ft, y_vel_ft);
	
	fp_e=fopen("energy","w");	//file to store energy vs time data
	fprintf(fp_e,"#time-Energy-Enstropy\n");
	
	spectra_file_id=init_spectra_storage();
	
	field_file_id=init_field_storage();
	
	//time marching
	for(k=0;k<s_flag;k++)
	{
	
		enstrophy=0.0;
		fftw_execute(p_inv_w);

		for(j=0;j<sys_size;j++)						
			for(k1=0;k1<sys_size;k1++)
			{	
				jk=j*y_size_r+k1;
				omega[jk]*=scale2;
				enstrophy+=omega[jk]*omega[jk];
				
			}
		
		if(pflag==0)
		{
		store_field(omega,field_file_id,k);
		find_e_spectra(x_vel_ft, y_vel_ft, e_spectra);	//calculates energy spectra
		store_spectrum(e_spectra,spectra_file_id,k);
		}
		fftw_execute(p_for_w);

		
		
		enstrophy*=scale2;
		energy=find_energy(e_spectra);		//calculate total energy
		fprintf(fp_e,"%f %E %E\n", time, energy, enstrophy);
		fflush(fp_e);

		fp=fopen("init","w");		//stores current vorticity profile (in fourier space)
		fwrite(&omega_ft,sizeof(omega_ft[0]),sys_size*y_size_f,fp);//for future runs
		fclose(fp);
				
		for(i=0;i<t_flag;i++)	
		{	
			/*Solves the differential in fourier space using second order Runka Kutta scheme
			takes current omega_ft, x_velocity_ft and y_velocity_ft and gives the updated omega_ft
			omega_t is a temporary variable.*/
			time=time+dt;
			solve_rk2(omega_ft, omega_t, x_vel_ft, y_vel_ft);    
		}
			
	}
	fclose(fp_e);		

	fp=fopen("init","w");		//stores current vorticity profile (in fourier space)
	fwrite(omega_ft,sys_size*(sys_size+2)*sizeof(double),1,fp);//for future runs
	fclose(fp);
		
	fftw_destroy_plan(p_for_w);	//free heap variables
	fftw_destroy_plan(p_inv_w);
	fftw_destroy_plan(p_for_xvel);
	fftw_destroy_plan(p_inv_xvel);
	fftw_destroy_plan(p_for_yvel);
	fftw_destroy_plan(p_inv_yvel);
	
	H5Fclose (field_file_id);
	H5Fclose (spectra_file_id);
	
	free(omega);
	free(omega_t);
	free(x_vel);
	free(y_vel);
	free(exp_dt);
	free(e_spectra);
	
	return 0;
}//end of MAIN()
