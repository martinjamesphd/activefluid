#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<time.h>
#include"ns2d.h"

//function to initialize time increment and density of states
void init_den_state()
{
	int i, j, ij;
	int y_size_f=sys_size/2+1;
	int y_alias = sys_size/4;
	
	for(i=0;i<y_size_f;i++)
		for(j=0;j<y_size_f;j++)
		{
			ij=i*y_size_f+j;
			exp_dt[ij]=exp(-(vis*scale*scale*(i*i+j*j)-vis2*scale*scale*(i*i+j*j)*scale*scale*(i*i+j*j)+mu)*dt/2);
			
			if((i*i+j*j)>y_alias*y_alias)	//accounts for aliasing (vorticity evolution uses this factor)
				exp_dt[ij]=0.0;
		}
	 
	for(i=y_size_f;i<sys_size;i++)	//i greater than N/2 corresponds to kx = i-N (where N is sys_size)
		for(j=0;j<y_size_f;j++)
		{
		 	ij=i*y_size_f+j;
		 	exp_dt[ij]=exp(-(vis*scale*scale*((sys_size-i)*(sys_size-i)+j*j)-vis2\
		 	  		*scale*scale*((sys_size-i)*(sys_size-i)+j*j)*scale*\
					scale*((sys_size-i)*(sys_size-i)+j*j)+mu)*dt/2);
		 	
		 	if(((sys_size-i)*(sys_size-i)+j*j)>y_alias*y_alias)
		 		exp_dt[ij]=0.0;
		} 
		exp_dt[0]=0.0;
}

//initializes omega in fourier space
void init_omega(fftw_complex *omega_ft)
{
	int i, j, ij;
	int y_size_r = sys_size+2;
	
	double *omega=(double*)omega_ft;
	
	srand(time(NULL));	//seeds random function using current time
	
	for(i=0;i<sys_size;i++)
		for(j=0;j<sys_size;j++)
		{
			ij=i*y_size_r+j;
			
			omega[ij]=(double)((double)rand()/RAND_MAX)*10-5.0;
		}
	
	fftw_execute(p_for_w);
}
