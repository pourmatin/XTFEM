#include <Kokkos_Core.hpp>
#include "kokkos_kernel.h"
#include <stdio.h>

/* Thu July 21 2016  */
/* Developed by Hossein Pourmatin @ UTD for XTFEM 3D */

// Type of a one-dimensional length-N array of int.

typedef Kokkos::View<int*> view_1D_int;
typedef view_1D_int::HostMirror host_view_1D_int;

// Type of a one-dimensional length-N array of double.
typedef Kokkos::View<double*> view_1D_double;
typedef view_1D_double::HostMirror host_view_1D_double;

extern "C" /* IMPORTANT! */

void kokkos_host_(int rank, int ne, int nn, int ndofn, int nnse, int nnte,
    int ngpe, int ncmp, int nip, int *nodes, double *Nx, double *Ny,   
    double *Nz, double *Nt, double *sol,double *eps_old, double *eps_pl, 
    double *sig_eff, double *sig_back, double *p, double *D, double *ws,  
    double *del_po, double *sig_d)
{

/*------------------------------------------------------------------------------------------------------
    Initialize memory and copy data to device
------------------------------------------------------------------------------------------------------*/
    /* Calculate data size */
    size_t size_nodes, size_Nx, size_Nt, size_sol;
    size_t size_eps, size_p;
    size_nodes = nnse*ngpe*ne;
    size_Nx = nnse*ngpe*ne;
    size_Nt = 2*nnte*nip;
    size_sol = 2*nnte*nn*ndofn;
    size_eps = ncmp*ngpe*ne;
    size_p = ngpe*ne;
    
    /*  Declare Device Views    */
    view_1D_int d_nodes ("Node", size_nodes);
    view_1D_double d_Nx ("Shape_Func_x", size_Nx);
    view_1D_double d_Ny ("Shape_Func_y", size_Nx);
    view_1D_double d_Nz ("Shape_Func_z", size_Nx);
    view_1D_double d_Nt ("Shape_Func_t", size_Nt);
    view_1D_double d_sol ("Displacement", size_sol);
    view_1D_double d_eps_old ("Strain_Old", size_eps);
    view_1D_double d_eps_pl ("Strain_Plastic", size_eps);
    view_1D_double d_sig_eff ("Stress_Eff", size_eps);
    view_1D_double d_sig_back ("Stress_Back", size_eps);
    view_1D_double d_p ("Plasticity", size_p);
    view_1D_double d_D ("Damage", size_p);
    view_1D_double d_ws ("Energy", size_p);
    view_1D_double d_del_po ("delta", size_p);
    view_1D_double d_sig_d ("Stress_D", size_p);
    
    /*  Declare Host Views    */
    host_view_1D_int h_nodes = Kokkos::create_mirror_view (d_nodes);
    host_view_1D_double h_Nx = Kokkos::create_mirror_view (d_Nx);
    host_view_1D_double h_Ny = Kokkos::create_mirror_view (d_Ny);
    host_view_1D_double h_Nz = Kokkos::create_mirror_view (d_Nz);
    host_view_1D_double h_Nt = Kokkos::create_mirror_view (d_Nt);
    host_view_1D_double h_sol = Kokkos::create_mirror_view (d_sol);
    host_view_1D_double h_eps_old = Kokkos::create_mirror_view (d_eps_old);
    host_view_1D_double h_eps_pl = Kokkos::create_mirror_view (d_eps_pl);
    host_view_1D_double h_sig_eff = Kokkos::create_mirror_view (d_sig_eff);
    host_view_1D_double h_sig_back = Kokkos::create_mirror_view (d_sig_back);
    host_view_1D_double h_p = Kokkos::create_mirror_view (d_p);
    host_view_1D_double h_D = Kokkos::create_mirror_view (d_D);
    host_view_1D_double h_ws = Kokkos::create_mirror_view (d_ws);
    host_view_1D_double h_del_po = Kokkos::create_mirror_view (d_del_po);
    host_view_1D_double h_sig_d = Kokkos::create_mirror_view (d_sig_d);
    
    /*  Initialize host Views   */
    for (int i=0; i<size_nodes; i++) {
        h_nodes(i) = nodes[i];
    }
    for (int i=0; i<size_Nx; i++) {
        h_Nx(i) = Nx[i];
        h_Ny(i) = Ny[i];
        h_Nz(i) = Nz[i];
    }
    for (int i=0; i<size_Nt; i++) {
        h_Nt(i) = Nt[i];
    }
    for (int i=0; i<size_sol; i++) {
        h_sol(i) = sol[i];
    }
    for (int i=0; i<size_eps; i++) {
        h_eps_old(i) = eps_old[i];
        h_eps_pl(i) = eps_pl[i];
        h_sig_eff(i) = sig_eff[i];
        h_sig_back(i) = sig_back[i];
    }
    for (int i=0; i<size_p; i++) {
        h_sig_d(i) = sig_d[i];
        h_del_po(i) = del_po[i];
	//printf("del_po( %i ) = %f\n",i,del_po[i]);
        h_p(i) = p[i];
        h_D(i) = D[i];
        h_ws(i) = ws[i];
    }
    
    /* Copy data to device */
    Kokkos::deep_copy (d_nodes, h_nodes);
    Kokkos::deep_copy (d_Nx, h_Nx);
    Kokkos::deep_copy (d_Ny, h_Ny);
    Kokkos::deep_copy (d_Nz, h_Nz);
    Kokkos::deep_copy (d_Nt, h_Nt);
    Kokkos::deep_copy (d_sol, h_sol);
    Kokkos::deep_copy (d_eps_old, h_eps_old);
    Kokkos::deep_copy (d_eps_pl, h_eps_pl);
    Kokkos::deep_copy (d_sig_eff, h_sig_eff);
    Kokkos::deep_copy (d_sig_back, h_sig_back);
    Kokkos::deep_copy (d_sig_d, h_sig_d);
    Kokkos::deep_copy (d_del_po, h_del_po);
    Kokkos::deep_copy (d_p, h_p);
    Kokkos::deep_copy (d_D, h_D);
    Kokkos::deep_copy (d_ws, h_ws);
    /*if (rank==0) {
	for (int i=0; i<size_sol; i++) {
		if (sol[i] > 0.000001) printf("sol[%i] = %e\n",i, sol[i]);
	}
    }*/
/*------------------------------------------------------------------------------
            Invoke kokkos kernel
------------------------------------------------------------------------------*/
    typedef Kokkos::TeamPolicy<>    team_policy ;
    int threadsPerTeam = 256;
    //unsigned int team_size = team_policy::team_size_recommended();
    int perTeamSize = 2048;
    int perThreadSize = 32;
    int League = (ngpe*ne+threadsPerTeam-1)/threadsPerTeam;
    /* Launch "League" number of teams of the maximum number of threads per team */
    const team_policy policy( League , threadsPerTeam );
    Kokkos::parallel_for(ngpe*ne, Kokkos_Kernel(rank, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, d_nodes, d_Nx, d_Ny, d_Nz, d_Nt, d_sol, d_eps_old, d_eps_pl, d_sig_eff, d_sig_back, d_p, d_D, d_ws, d_del_po, d_sig_d) );
/*------------------------------------------------------------------------------
    Copy data from device and cleanup memory
------------------------------------------------------------------------------*/
    /* Copy data from device */
    Kokkos::deep_copy (h_eps_old, d_eps_old);
    Kokkos::deep_copy (h_eps_pl, d_eps_pl);
    Kokkos::deep_copy (h_sig_eff, d_sig_eff);
    Kokkos::deep_copy (h_sig_back, d_sig_back);
    Kokkos::deep_copy (h_sig_d, d_sig_d);
    Kokkos::deep_copy (h_del_po, d_del_po);
    Kokkos::deep_copy (h_p, d_p);
    Kokkos::deep_copy (h_D, d_D);
    Kokkos::deep_copy (h_ws, d_ws);
    
    /*  Copy the results from View to the original arrays for MPI   */
    for (int i=0; i<size_eps; i++) {
        eps_old[i] = h_eps_old(i);
        eps_pl[i] = h_eps_pl(i);
        sig_eff[i] = h_sig_eff(i);
        sig_back[i] = h_sig_back(i);
    }
    for (int i=0; i<size_p; i++) {
        sig_d[i] = h_sig_d(i);
        del_po[i] = h_del_po(i);
        p[i] = h_p(i);
        D[i] = h_D(i);
        ws[i] = h_ws(i);
	//std::cout<<"ws[ "<<i<<" ] = "<<ws[i]<<endl;
    }
}
