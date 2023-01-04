#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<math.h>
#include<omp.h>
#include"ns2d.h"

//function to calculate stream function from vorticity in fourier space	
void find_vel_ft(fftw_complex *omega_ft, fftw_complex *x_vel_ft, fftw_complex *y_vel_ft)
{
	int i, j, ij, y_size=sys_size/2+1;
	double temp_re, temp_im;
	
	y_vel_ft[0][0]=0;
	y_vel_ft[0][1]=0;	

	//openmp_s
	# pragma omp parallel private(i, j, ij, temp_re, temp_im)
	{
		
	# pragma omp for schedule(static) nowait
	for(j=1;j<y_size;j++)//takes care when kx=0 and ky=0
	{
		y_vel_ft[j][0]=-omega_ft[j][0]/(scale*scale*j*j);
		y_vel_ft[j][1]=-omega_ft[j][1]/(scale*scale*j*j);
	}

	# pragma omp for schedule(static) nowait
	for(i=1;i<y_size;i++)		//calculates and temporarily stores stream function in y_vel_ft
	{
		for(j=0;j<y_size;j++)
		{
			ij=i*y_size+j;
			y_vel_ft[ij][0]=-omega_ft[ij][0]/(scale*scale*(i*i+j*j));
			y_vel_ft[ij][1]=-omega_ft[ij][1]/(scale*scale*(i*i+j*j));
		}
	}

	# pragma omp for schedule(static)
	for(i=y_size;i<sys_size;i++)//i greater than N/2 corresponds to kx = i-N (where N is sys_size)
	{
		for(j=0;j<y_size;j++)
		{
			ij=i*y_size+j;
			y_vel_ft[ij][0]=-omega_ft[ij][0]/(scale*scale*((sys_size-i)*(sys_size-i)+j*j));
			y_vel_ft[ij][1]=-omega_ft[ij][1]/(scale*scale*((sys_size-i)*(sys_size-i)+j*j));
		}
	}
	
	# pragma omp for schedule(static)
	for(i=0;i<sys_size;i++)		//calculates fourier of x velocity from stream function (currently stored in y_vel_ft)
	{
		for(j=0;j<y_size;j++)
		{
			ij=i*y_size+j;
			x_vel_ft[ij][0]=scale*y_vel_ft[ij][1]*(j);
			x_vel_ft[ij][1]=scale*y_vel_ft[ij][0]*(-j);
			
		}
	}
	
	# pragma omp for schedule(static) nowait
	for(i=0;i<y_size;i++)		//calculates fourier of x velocity from stream function (currently stored in y_vel_ft)
	{
		for(j=0;j<y_size;j++)
		{
			ij=i*y_size+j;
			temp_re=scale*y_vel_ft[ij][1]*(-i);
			temp_im=scale*y_vel_ft[ij][0]*(i);
			y_vel_ft[ij][0]=temp_re;
			y_vel_ft[ij][1]=temp_im;
		}
	}
		
	# pragma omp for schedule(static)
	for(i=y_size;i<sys_size;i++)//i greater than N/2 corresponds to kx = i-N (where N is sys_size)
	{
		for(j=0;j<y_size;j++)
		{
			ij=i*y_size+j;
			temp_re=scale*y_vel_ft[ij][1]*(-(i-sys_size));
			temp_im=scale*y_vel_ft[ij][0]*((i-sys_size));
			y_vel_ft[ij][0]=temp_re;
			y_vel_ft[ij][1]=temp_im;
		}
	}
	}
	//openmp_e
}

//function to calculate energy spectrum
void find_e_spectra(fftw_complex *x_vel_ft, fftw_complex *y_vel_ft, double *e_spectra)
{
	extern int sys_size;
	int i, j, ij, k, y_size = sys_size/2+1;
	
	double energyNorm = 1.0/((double)sys_size*sys_size*sys_size*sys_size);
	
	for(k=0;k<y_size;k++)
		e_spectra[k]=0;
		
	for(i=0;i<y_size;i++)
		for(j=1;j<y_size;j++)
		{
			ij = i*y_size+j;
			k = (int)sqrt(i*i+j*j);
			if(k<y_size)
			{
				e_spectra[k]+=2*x_vel_ft[ij][0]*x_vel_ft[ij][0];
				e_spectra[k]+=2*x_vel_ft[ij][1]*x_vel_ft[ij][1];
				e_spectra[k]+=2*y_vel_ft[ij][0]*y_vel_ft[ij][0];
				e_spectra[k]+=2*y_vel_ft[ij][1]*y_vel_ft[ij][1];
			}
		}
		
	for(i=0;i<y_size;i++)
		{
			ij = i*y_size;
			k = i;
			if(k<y_size)
			{
				e_spectra[k]+=x_vel_ft[ij][0]*x_vel_ft[ij][0];
				e_spectra[k]+=x_vel_ft[ij][1]*x_vel_ft[ij][1];
				e_spectra[k]+=y_vel_ft[ij][0]*y_vel_ft[ij][0];
				e_spectra[k]+=y_vel_ft[ij][1]*y_vel_ft[ij][1];
			}
		}
		
	for(i=y_size;i<sys_size;i++)//i greater than N/2 corresponds to kx = i-N (where N is sys_size)
		for(j=1;j<y_size;j++)
		{
			ij=i*y_size+j;
			k = (int)sqrt((sys_size-i)*(sys_size-i)+j*j);
			if(k<y_size)
			{
				e_spectra[k]+=2*x_vel_ft[ij][0]*x_vel_ft[ij][0];
				e_spectra[k]+=2*x_vel_ft[ij][1]*x_vel_ft[ij][1];
				e_spectra[k]+=2*y_vel_ft[ij][0]*y_vel_ft[ij][0];
				e_spectra[k]+=2*y_vel_ft[ij][1]*y_vel_ft[ij][1];
			}
		}
		
	for(i=y_size;i<sys_size;i++)
		{
			ij=i*y_size;
			k = (sys_size-i);
			if(k<y_size)
			{
				e_spectra[k]+=x_vel_ft[ij][0]*x_vel_ft[ij][0];
				e_spectra[k]+=x_vel_ft[ij][1]*x_vel_ft[ij][1];
				e_spectra[k]+=y_vel_ft[ij][0]*y_vel_ft[ij][0];
				e_spectra[k]+=y_vel_ft[ij][1]*y_vel_ft[ij][1];
			}
		}
	
	for(k=0;k<y_size;k++)
		e_spectra[k]*=0.5*energyNorm/scale;
}

//function to calculate total energy
double find_energy(double *e_spectra)
{
	extern int sys_size;
	int k, y_size=sys_size/2+1;
	
	double energy=0.0;
	
	for(k=0;k<y_size;k++)
		energy+=e_spectra[k]*scale;
		
	return energy;
}

//function to calculate epsilon
double find_epsilon(fftw_complex *omega_ft)
{
	double epsilon=0, y_size=sys_size/2+1;
	int i, j, ij;
	
	for(i=0;i<sys_size;i++)
	{
		ij=i*y_size;
		epsilon+=vis*(omega_ft[ij][0]*omega_ft[ij][0]+omega_ft[ij][1]*omega_ft[ij][1])/(sys_size*sys_size);
		ij=ij+y_size-1;
		epsilon+=vis*(omega_ft[ij][0]*omega_ft[ij][0]+omega_ft[ij][1]*omega_ft[ij][1])/(sys_size*sys_size);
	}
		
	for(i=0;i<sys_size;i++)
		for(j=1;j<y_size-1;j++)
		{
			ij=i*y_size+j;
			epsilon+=2*vis*(omega_ft[ij][0]*omega_ft[ij][0]+omega_ft[ij][1]*omega_ft[ij][1])/(sys_size*sys_size);
		}
		
	epsilon/=(sys_size*sys_size);
	
	return epsilon;
}
