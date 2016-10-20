//
//  kokkos_kernel.h
//  
//
//  Created by Hossein Pourmatin on 7/28/16.
//
//

#ifndef kokkos_kernel_h
#define kokkos_kernel_h

#include <Kokkos_Core.hpp>
#include "damage.h"

typedef Kokkos::TeamPolicy<>              team_policy ;
typedef team_policy::member_type team_member ;
// Type of a one-dimensional length-N array of int.
typedef Kokkos::View<int*> view_1D_int;

// Type of a one-dimensional length-N array of double.
typedef Kokkos::View<double*> view_1D_double;

class Kokkos_Kernel {
    
    int ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip;
    view_1D_int nodes;
    view_1D_double Nx, Ny, Nz, Nt, sol, eps_old, eps_pl, sig_eff, sig_back, p, d, ws, del_po, sig_d;
    
    /* THIS CODE SHOULD BE ONLY USED FOR C3D8 ELEMENT TYPE */
    void dispinterp(int, int, view_1D_double) const;
    void calstrain( int, view_1D_double, view_1D_double ) const;
    void caldamage( int, view_1D_double ) const;
    
public:
    Kokkos_Kernel(int _ne, int _nn, int _ndofn, int _nnse, int _nnte, int _ngpe, int _ncmp, int _nip, view_1D_int _nodes, view_1D_double _Nx, view_1D_double _Ny, view_1D_double _Nz,     view_1D_double _Nt, view_1D_double _sol, view_1D_double _eps_old, view_1D_double _eps_pl,view_1D_double _sig_eff, view_1D_double _sig_back, view_1D_double _p, view_1D_double _d, view_1D_double _ws, view_1D_double _del_po, view_1D_double _sig_d):
    ne(_ne), nn(_nn), ndofn(_ndofn), nnse(_nnse), nnte(_nnte), ngpe(_ngpe), ncmp(_ncmp), nip(_nip),
    nodes(_nodes), Nx(_Nx), Ny(_Ny), Nz(_Nz), Nt(_Nt), sol(_sol), eps_old(_eps_old),
    eps_pl(_eps_pl), sig_eff(_sig_eff), sig_back(_sig_back), p(_p), d(_d), ws(_ws), del_po(_del_po),
    sig_d(_sig_d) {};
    
    KOKKOS_INLINE_FUNCTION
    void operator() ( const int ii/*team_member & thread*/ ) const {
        //int ii = thread.league_rank() * thread.team_size() + thread.team_rank();
        if (ii < ngpe*ne)
        {
            /* loop over temporal interpolation points */
            for (int jj = 0; jj < nip; jj++)
            {
                /* Interpolate displacement solution */
                double disp[24] = {0.0};
                //dispinterp(ii, jj, disp);
                for (int i = 0; i < nnse; i++)
                {
                    int ndid = nodes(ii*nnse+i) - 1;
                    for (int j = 0; j < ndofn; j++)
                    {
                        for (int k = 0; k < 2*nnte; k++)
                        {
                            int dofid = k*nn*ndofn + ndid*ndofn + j;
                            disp[i*ndofn+j] += Nt(jj*2*nnte+k)*sol(dofid);
                        }
                    }
                }
                
                /* Calculate strain */
                double strain[6] = {0.0};
                //calstrain( ii, disp, strain );
                for (int i = 0; i < nnse; i++)
                {
                    strain[0] += Nx(ii*nnse+i)*disp[i*ndofn+0];
                    strain[1] += Ny(ii*nnse+i)*disp[i*ndofn+1];
                    strain[2] += Nz(ii*nnse+i)*disp[i*ndofn+2];
                    strain[3] += Ny(ii*nnse+i)*disp[i*ndofn+0] + Nx(ii*nnse+i)*disp[i*ndofn+1];
                    strain[4] += Nz(ii*nnse+i)*disp[i*ndofn+1] + Ny(ii*nnse+i)*disp[i*ndofn+2];
                    strain[5] += Nz(ii*nnse+i)*disp[i*ndofn+0] + Nx(ii*nnse+i)*disp[i*ndofn+2];
                }
                /* Calculate damage */
                //caldamage( ii, strain );
                double del_p = 0.0;
                double dmat[6][6] = {0.0};
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        dmat[i][j] = COEF*nu;
                    }
                }
                for (int i = 0; i < 3; i++)
                {
                    dmat[i][i] = COEF*(1.0-nu);
                    dmat[i+3][i+3] = COEF*(1.0-2.0*nu)/2.0;
                }
                /*------------------------------------------------------------------------------
                 Elastic prediction phase
                 ------------------------------------------------------------------------------*/
                double c1 = 1.0/(1.0-b*d(ii));
                double c2 = (a-b)*d(ii)/((1.0-b*d(ii))*(1.0-a*d(ii)));
                double c3 = b*(1.0-d(ii))/(1.0-b*d(ii));
                double test_eps_l[6] = {0.0}; /* total test local strain */
                double test_eps_el[6] = {0.0}; /* test local elastic strain*/
                for (int i = 0; i < 6; i++)
                {
                    test_eps_l[i] = c1*strain[i] + c3*eps_pl(ii*ncmp+i);
                }
                double tr_eps_new = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    tr_eps_new += strain[i];
                }
                tr_eps_new = c2*tr_eps_new/3.0;
                for (int i = 0; i < 3; i++)
                {
                    test_eps_l[i] += tr_eps_new;
                }
                for (int i = 0; i < 6; i++)
                {
                    test_eps_el[i] = test_eps_l[i] - eps_pl(ii*ncmp+i);
                }
                double test_sig_eff[6] = {0.0};
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        test_sig_eff[i] += dmat[i][j]*test_eps_el[j];
                    }
                }
                double tr_test_sig_eff = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    tr_test_sig_eff += test_sig_eff[i];
                }
                tr_test_sig_eff = tr_test_sig_eff/3.0;
                double test_dev_sig_eff[6] = {0.0};
                test_dev_sig_eff[0] = test_sig_eff[0] - tr_test_sig_eff;
                test_dev_sig_eff[1] = test_sig_eff[1] - tr_test_sig_eff;
                test_dev_sig_eff[2] = test_sig_eff[2] - tr_test_sig_eff;
                test_dev_sig_eff[3] = test_sig_eff[3];
                test_dev_sig_eff[4] = test_sig_eff[4];
                test_dev_sig_eff[5] = test_sig_eff[5];
                double test_sig[6] = {0.0};
                for (int i = 0; i < 6; i++)
                {
                    test_sig[i] = test_dev_sig_eff[i] - sig_back(ii*ncmp+i);
                }
                double test_sig_eq = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    test_sig_eq += 1.5*test_sig[i]*test_sig[i];
                }
                for (int i = 3; i < 6; i++)
                {
                    test_sig_eq += 3.0*test_sig[i]*test_sig[i];
                }
                test_sig_eq = sqrt(test_sig_eq);
                if (test_sig_eq <= sig_f)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        sig_eff(ii*ncmp+i) = test_sig_eff[i];
                        eps_old(ii*ncmp+i) = strain[i];
                    }
                    sig_d(ii) = test_sig_eq;
                }
                /*------------------------------------------------------------------------------
                 Plastic correction phase
                 ------------------------------------------------------------------------------*/
                else
                {
                    double dstrain[6] = {0.0};
                    for (int i = 0; i < 6; i++)
                    {
                        dstrain[i] = strain[i] - eps_old(ii*ncmp+i);
                    }
                    double dstress[6] = {0.0};
                    for (int i = 0; i < 6; i++)
                    {
                        for (int j = 0; j < 6; j++)
                        {
                            dstress[i] += dmat[i][j]*dstrain[j];
                        }
                        dstress[i] = dstress[i]/(1.0-b*d(ii));
                    }
                    double c4 = K*(a-b)*d(ii)/((1.0-b*d(ii))*(1.0-a*d(ii)));
                    double tr_dstrain = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        tr_dstrain += dstrain[i];
                    }
                    for (int i = 0; i < 3; i++)
                    {
                        dstress[i] += c4*tr_dstrain;
                    }
                    double Qs[6] = {0.0};
                    for (int i = 0; i < 6; i++)
                    {
                        Qs[i] = sig_back(ii*ncmp+i) - sig_eff(ii*ncmp+i) - dstress[i];
                    }
                    double gamma = Cy*(1.0-d(ii)) + 3.0*G*(1.0-b)/(1.0-b*d(ii));
                    double tr_Qs = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        tr_Qs += Qs[i];
                    }
                    tr_Qs = tr_Qs/3.0;
                    double dev_Qs[6] = {0.0};
                    dev_Qs[0] = Qs[0] - tr_Qs;
                    dev_Qs[1] = Qs[1] - tr_Qs;
                    dev_Qs[2] = Qs[2] - tr_Qs;
                    dev_Qs[3] = Qs[3];
                    dev_Qs[4] = Qs[4];
                    dev_Qs[5] = Qs[5];
                    double Qs_eq = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        Qs_eq += 1.5*dev_Qs[i]*dev_Qs[i];
                    }
                    for (int i = 3; i < 6; i++)
                    {
                        Qs_eq += 3.0*dev_Qs[i]*dev_Qs[i];
                    }
                    Qs_eq = sqrt(Qs_eq);
                    double tr_sigs = -tr_Qs;
                    del_p = (Qs_eq - sig_f)/gamma;
                    p(ii) += del_p;
                    double dev_sigs[6] = {0.0};
                    double c5 = 1.0 + gamma*del_p/sig_f;
                    for (int i = 0; i < 6; i++)
                    {
                        dev_sigs[i] = -dev_Qs[i]/c5;
                    }
                    double sigs[6] = {0.0};
                    sigs[0] = dev_sigs[0] + tr_sigs;
                    sigs[1] = dev_sigs[1] + tr_sigs;
                    sigs[2] = dev_sigs[2] + tr_sigs;
                    sigs[3] = dev_sigs[3];
                    sigs[4] = dev_sigs[4];
                    sigs[5] = dev_sigs[5];
                    /*------------------------------------------------------------------------------
                     Update internal variables
                     ------------------------------------------------------------------------------*/
                    double dstrain_pl[6] = {0.0};
                    for (int i = 0; i < 6; i++)
                    {
                        dstrain_pl[i] = 3.0*dev_sigs[i]*del_p/(2.0*sig_f);
                        eps_pl(ii*ncmp+i) += dstrain_pl[i];
                    }
                    for (int i = 0; i < 6; i++)
                    {
                        sig_back(ii*ncmp+i) += 2.0*Cy*(1.0-d(ii))*dstrain_pl[i]/3.0;
                    }
                    for (int i = 0; i < 6; i++)
                    {
                        sig_eff(ii*ncmp+i) = sigs[i] + sig_back(ii*ncmp+i);
                    }
                    double tr_sig_eff = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        tr_sig_eff += sig_eff(ii*ncmp+i);
                    }
                    tr_sig_eff = tr_sig_eff/3.0;
                    double dev_sig_eff[6] = {0.0};
                    dev_sig_eff[0] = sig_eff(ii*ncmp+0) - tr_sig_eff;
                    dev_sig_eff[1] = sig_eff(ii*ncmp+1) - tr_sig_eff;
                    dev_sig_eff[2] = sig_eff(ii*ncmp+2) - tr_sig_eff;
                    dev_sig_eff[3] = sig_eff(ii*ncmp+3);
                    dev_sig_eff[4] = sig_eff(ii*ncmp+4);
                    dev_sig_eff[5] = sig_eff(ii*ncmp+5);
                    double sig_eff_eq = 0.0;
                    for (int i = 0; i < 3; i++)
                    {
                        sig_eff_eq += 1.5*dev_sig_eff[i]*dev_sig_eff[i];
                    }
                    for (int i = 3; i < 6; i++)
                    {
                        sig_eff_eq += 3.0*dev_sig_eff[i]*dev_sig_eff[i];
                    }
                    sig_eff_eq = sqrt(sig_eff_eq);
                    double sig_new = fabs(sig_eff_eq - sig_f);
                    double sig_old = fabs(sig_d(ii) - sig_f);
                    /* Note: int abs(int n), double fabs(double x) */
                    ws(ii) += (sig_new*del_p + sig_old*del_po(ii))/2.0;
                    if (ws(ii) >= wd)
                    {
                        double sig_eff_p[3] = {0.0};
                        double sig_eff_n[3] = {0.0};
                        for (int i = 0; i < 3; i++)
                        {
                            if (sig_eff(ii*ncmp+i) < 0.0)
                            {
                                sig_eff_n[i] = sig_eff(ii*ncmp+i);
                            }
                            else
                            {
                                sig_eff_p[i] = sig_eff(ii*ncmp+i);
                            }
                        }
                        double Y_p = (1.0+nu)/(2.0*E)*(sig_eff_p[0]*sig_eff_p[0] + 
                                                       sig_eff_p[1]*sig_eff_p[1] + sig_eff_p[2]*sig_eff_p[2] + 
                                                       h*((1.0-d(ii))/(1.0-h*d(ii)))*((1.0-d(ii))/(1.0-h*d(ii)))*
                                                       (sig_eff_n[0]*sig_eff_n[0] + sig_eff_n[1]*sig_eff_n[1] + 
                                                        sig_eff_n[2]*sig_eff_n[2]));
                        tr_sig_eff = 0.0;
                        for (int i = 0; i < 3; i++)
                        {
                            tr_sig_eff += sig_eff(ii*ncmp+i);
                        }
                        double Y_n = 0.0;
                        if (tr_sig_eff >= 0.0)
                        {
                            Y_n = nu/(2.0*E)*tr_sig_eff*tr_sig_eff;
                        }
                        else
                        {
                            Y_n = nu/(2.0*E)*h*((1.0-d(ii))/(1.0-h*d(ii)))*
                            ((1.0-d(ii))/(1.0-h*d(ii)))*tr_sig_eff*tr_sig_eff;
                        }
                        double Y = Y_p - Y_n;
                        d(ii) += pow(Y/S, s)*del_p;
                    }
                    else
                    {
                        d(ii) = 0.0;
                    }
                    for (int i = 0; i < 6; i++)
                    {
                        eps_old(ii*ncmp+i) = strain[i];
                    }
                    sig_d(ii) = sig_eff_eq;
                    del_po(ii) = del_p;
                }
            }
        }
    }
};

#endif /* kokkos_kernel_h */
