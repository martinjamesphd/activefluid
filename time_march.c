#include<stdio.h>
#include<stdlib.h>
#include<fftw3.h>
#include<math.h>
#include"ns2d.h"

//function to evaluate jacobian (stores jacobian_ft in variable 'omega')
void find_jacobian_ft(double *omega, double *x_vel, double *y_vel, double *vx_vsq_ft_0, double *vy_vsq_ft_0)
{
	int i,j,ij;
	int y_size_f = sys_size/2+1;
	int y_size_r = sys_size+2;
	
	double *vx_vsq 	= (double*)malloc(sys_size*(sys_size+2)*sizeof(double));
	double *vy_vsq 	= (double*)malloc(sys_size*(sys_size+2)*sizeof(double));
	
	fftw_complex *vx_omega_ft = (fftw_complex*)x_vel;
	fftw_complex *vy_omega_ft = (fftw_complex*)y_vel;
	fftw_complex *jacobian_ft = (fftw_complex*)omega;
	fftw_complex *vx_vsq_ft = (fftw_complex*)vx_vsq;
	fftw_complex *vy_vsq_ft = (fftw_complex*)vy_vsq;
	
	fftw_plan p_for_vx_vsq = fftw_plan_dft_r2c_2d(sys_size, sys_size, vx_vsq, vx_vsq_ft, FFTW_ESTIMATE);
	fftw_plan p_for_vy_vsq = fftw_plan_dft_r2c_2d(sys_size, sys_size, vy_vsq, vy_vsq_ft, FFTW_ESTIMATE);
	
	//openmp_s
	# pragma omp parallel for schedule(static) private(i, j, ij)
	for(i=0;i<sys_size;i++)
		for (j=0;j<y_size_r;j++)
		{
			ij = i*y_size_r + j;
			vx_vsq[ij] = mu2*x_vel[ij]*(x_vel[ij]*x_vel[ij]+y_vel[ij]*y_vel[ij]); 
			vy_vsq[ij] = mu2*y_vel[ij]*(x_vel[ij]*x_vel[ij]+y_vel[ij]*y_vel[ij]);
			x_vel[ij] = nlin*x_vel[ij]*omega[ij]; 
			y_vel[ij] = nlin*y_vel[ij]*omega[ij];
		}
	//openmp_e
	
	fftw_execute(p_for_yvel);
	fftw_execute(p_for_xvel);
	
	fftw_execute(p_for_vx_vsq);
	fftw_execute(p_for_vy_vsq);
	
	//openmp_s
	# pragma omp parallel
	{
	# pragma omp for schedule(static) private(i, j, ij) 
	for(i=0;i<y_size_f;i++)
		for (j=0; j<y_size_f; j++)
	   	{
			ij = i*y_size_f + j;
			jacobian_ft[ij][0] = scale*((-(i*vx_omega_ft[ij][1] + j*vy_omega_ft[ij][1]))\
						+(-(i*vy_vsq_ft[ij][1] - j*vx_vsq_ft[ij][1])));
			jacobian_ft[ij][1] = scale*(i*vx_omega_ft[ij][0] + j*vy_omega_ft[ij][0]\
						+i*vy_vsq_ft[ij][0] - j*vx_vsq_ft[ij][0]);
	   	}
	   	
	# pragma omp for schedule(static) private(i, j, ij)	
	for(i=y_size_f;i<sys_size;i++)//i greater than N/2 corresponds to kx = i-N (where N is sys_size)
	   	for(j=0;j<y_size_f;j++)
	   	{
	   		ij=i*y_size_f+j;
	   		jacobian_ft[ij][0] = scale*((-((i-sys_size)*vx_omega_ft[ij][1]+j*vy_omega_ft[ij][1]))\
	   					+(-((i-sys_size)*vy_vsq_ft[ij][1]-j*vx_vsq_ft[ij][1])));
	   		jacobian_ft[ij][1] = scale*((i-sys_size)*vx_omega_ft[ij][0]+j*vy_omega_ft[ij][0]\
	   					+(i-sys_size)*vy_vsq_ft[ij][0]-j*vx_vsq_ft[ij][0]);
	   	}
	}
	//openmp_e
	
	free(vx_vsq);
	free(vy_vsq);
	
	fftw_destroy_plan(p_for_vx_vsq);	//free heap variables
	fftw_destroy_plan(p_for_vy_vsq);
}

//function to solve using Runke Kutta 2
void solve_rk2 (fftw_complex *omega_ft, fftw_complex *omega_t, fftw_complex *x_vel_ft, fftw_complex *y_vel_ft)
{
	double *omega = (double*)omega_ft;
	double *x_vel = (double*)x_vel_ft;
	double *y_vel = (double*)y_vel_ft;
	double max_x_vel=0, max_y_vel=0;
	fftw_complex *jacobian_ft = omega_ft;	
	
	int i, j, ij;
	int y_size_f=sys_size/2+1;
	int y_size_r=sys_size+2;
	
	double vx_vsq_ft_0, vy_vsq_ft_0;
	
	//openmp_s
	# pragma omp parallel for schedule(static) private(i, j, ij)
	for(i=0;i<sys_size;i++)
	{
		for (j=0; j<y_size_f; j++)
	   	{
	   		ij = i*y_size_f + j;
	   		omega_t[ij][0] = omega_ft[ij][0];
	   		omega_t[ij][1] = omega_ft[ij][1];
	   	}
	}
	//openmp_e
	
	fftw_execute(p_inv_yvel);
	fftw_execute(p_inv_xvel);
	fftw_execute(p_inv_w);    
	
	//openmp_s
	//# pragma omp parallel for schedule(static) private(i, j, ij)
	for(i=0;i<sys_size;i++)
	{
		for (j=0; j<sys_size; j++)
		{
			ij = i*y_size_r + j;
			omega[ij] = omega[ij]/(sys_size*sys_size);	//normalising for FT inverse
			x_vel[ij] = x_vel[ij]/(sys_size*sys_size);
			y_vel[ij] = y_vel[ij]/(sys_size*sys_size);
			
			if(x_vel[ij]>max_x_vel) max_x_vel = x_vel[ij];
			if(-x_vel[ij]>max_x_vel) max_x_vel = -x_vel[ij];
			if(y_vel[ij]>max_y_vel) max_y_vel = y_vel[ij];
			if(-y_vel[ij]>max_y_vel) max_y_vel = -y_vel[ij];
		}
	}
	//openmp_e
	dt=(max_x_vel>max_y_vel)?(0.5*M_PI/(sys_size*scale*max_x_vel)):(0.5*M_PI/(sys_size*scale*max_y_vel));
		
	find_jacobian_ft(omega, x_vel, y_vel, &vx_vsq_ft_0, &vy_vsq_ft_0);         //update jacobian_ft 
	
	//RK step1
	//openmp_s
	# pragma omp parallel for schedule(static) private(i, j, ij)
	for(i=0;i<sys_size;i++)
	{
		for (j=0; j<y_size_f; j++)
	   	{	
	   		ij = i*y_size_f + j;
	   		omega_ft[ij][0] = exp_dt[ij]*omega_t[ij][0] + exp_dt[ij]*0.5*dt*(-jacobian_ft[ij][0]);
			omega_ft[ij][1] = exp_dt[ij]*omega_t[ij][1] + exp_dt[ij]*0.5*dt*(-jacobian_ft[ij][1]);
		}
	}
	//openmp_e		

	find_vel_ft( omega_ft, x_vel_ft, y_vel_ft);
	
	fftw_execute(p_inv_yvel);
	fftw_execute(p_inv_xvel);
	fftw_execute(p_inv_w);
	
	//openmp_s
	# pragma omp parallel for schedule(static) private(i, j, ij)
	for(i=0;i<sys_size;i++)
	{
		for (j=0; j<y_size_r; j++)
		{
			ij = i*y_size_r + j;
			omega[ij] = omega[ij]/(sys_size*sys_size);	//normalising for FT inverse
			x_vel[ij] = x_vel[ij]/(sys_size*sys_size);
			y_vel[ij] = y_vel[ij]/(sys_size*sys_size);
		}
	}
	//openmp_e
	
	find_jacobian_ft(omega, x_vel, y_vel, &vx_vsq_ft_0, &vy_vsq_ft_0);		//update jacobian_ft wrt new omega
	
	//RK step2
	//openmp_s
	# pragma omp parallel for schedule(static) private(i, j, ij)	
	for(i=0;i<sys_size;i++)
	{
		for (j=0; j<y_size_f; j++)
		{
			ij = i*y_size_f + j;
			omega_ft[ij][0] = exp_dt[ij]*exp_dt[ij]*omega_t[ij][0]+exp_dt[ij]*dt*(-jacobian_ft[ij][0]);
			omega_ft[ij][1] = exp_dt[ij]*exp_dt[ij]*omega_t[ij][1]+exp_dt[ij]*dt*(-jacobian_ft[ij][1]);
		}
	}
	//openmp_e

	find_vel_ft( omega_ft, x_vel_ft, y_vel_ft);
	init_den_state();
}
