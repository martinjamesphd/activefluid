#ifndef ns2d_h
#define ns2d_h

#include<hdf5.h>

extern fftw_plan p_for_xvel, p_inv_xvel, p_for_yvel, p_inv_yvel, p_for_w, p_inv_w;
extern double *exp_dt;			//variable to store time increments and aliasing (refer function init_den_state())
extern int sys_size, nthreads, nthreads_omp, tau_num,max_iter; 
extern double dt, vis, vis2, nlin, mu, mu2, scale; 
	
void init_den_state();

void init_omega(fftw_complex *omega_ft);

void find_vel_ft(fftw_complex *omega_ft, fftw_complex *x_vel_ft, fftw_complex *y_vel_ft);

void find_e_spectra(fftw_complex *x_vel, fftw_complex *y_vel, double *e_spectra);

void find_jacobian_ft(double *omega, double *x_vel, double *y_vel, double *vx_vsq_ft_0, double *vy_vsq_ft_0);

double find_energy(double *e_spectra);

double find_epsilon(fftw_complex *omega_ft);

void solve_rk2 (fftw_complex *omega_ft, fftw_complex *omega_t, fftw_complex *x_vel_ft, fftw_complex *y_vel_ft);

void store_spectrum(double *e_spectra, hid_t spectra_file_id, int k);

hid_t init_spectra_storage();

void store_field(double *omega, hid_t field_file_id, int k);

hid_t init_field_storage();

#endif
