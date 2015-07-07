/* Matrix Method Finite Difference Schrodinger/Poisson Solver 3.1*/
/* C.R.J. Fisher;
07_2012 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

void Normalise(gsl_vector *Single_eigvector_ptr, int Dim, double dx, FILE *fp_2)
{
    int n;
    double temp;
    double norm_factor, vector_sum=0;
    for (n=0; n<Dim; n++){
        temp = gsl_vector_get(Single_eigvector_ptr, n);
        vector_sum = vector_sum + temp*temp*dx;
    }
    norm_factor = sqrt(vector_sum);
    for (n=0; n<Dim; n++){
        temp = (gsl_vector_get(Single_eigvector_ptr, n))/norm_factor;
        gsl_vector_set(Single_eigvector_ptr, n, temp);
        fprintf(fp_2, "n%E", gsl_vector_get(Single_eigvector_ptr, n));
        /*output mth value of nth eigenvector*/
    }
}

int Schrodinger_Solve(double dx, int Dim, int barrier_L, int barrier_R, double E_well, double *Meff_ptr, double *subband_ptr,
gsl_eigen_symmv_workspace *Workspace, gsl_matrix *Hamil_ptr, gsl_vector *Eigvalue_ptr, gsl_matrix *Eigvector_ptr,
gsl_vector *s_eigvector_ptr, gsl_vector *Potential_ptr, gsl_matrix *Norm_Eigvector_ptr, FILE *fp_1, FILE *fp_2 ){
    int n, m, N=0;
    double Hbar=1.054571E-34, Element;
    double t_0=(Hbar*Hbar)/(2*dx*dx);
    for (n = 0; n < Dim; n++)
    {
        for (m = 0; m < Dim; m++)
        {
            if (n == m) /*main diagonal*/
            {
                Element = t_0/Meff_ptr[n-1] + t_0/Meff_ptr[n+1] + gsl_vector_get(Potential_ptr, n);
            }
            else if (n == m-1) /*subdiagonal*/
            {
                Element = -1*t_0/Meff_ptr[n-1];
            }
            else if (n == m+1) /*superdiagonal*/
            {
                Element = -1*t_0/Meff_ptr[n+1];
            }
            else
            {
                Element = 0;
                /*all other matrix elements*/
            }
            gsl_matrix_set(Hamil_ptr, n, m, Element);
            /*transfer value to matrix*/
        }
    }
    gsl_eigen_symmv(Hamil_ptr, Eigvalue_ptr, Eigvector_ptr, Workspace);
    /*solve matrix*/
    gsl_eigen_symmv_sort(Eigvalue_ptr, Eigvector_ptr, GSL_EIGEN_SORT_VAL_ASC);
    /*re-order the results in ascending order*/
    /*output eigenvalues*/
    for (n = 0; n < Dim; n++){
        double eigenvalue = gsl_vector_get(Eigvalue_ptr, n);
        if (eigenvalue>E_well && eigenvalue<0) /*select confined states*/
        {
            fprintf(fp_1, "n%E", (eigenvalue/1.6E-22));
            /*output in meV*/
            subband_ptr[n] = eigenvalue; /*For Fermi-Dirac */
            N = N+1;
            /* sum the number of eigenstates*/
        }
    }
    /*output eigenvectors*/
    for (n = 0; n < N; n++){  /*retrieve eigenvectors for confined states only*/
        fprintf(fp_2, "nnnn");
        gsl_matrix_get_col(s_eigvector_ptr, Eigvector_ptr, n);
        /*retrieve nth eigenvector*/
        Normalise(s_eigvector_ptr, Dim, dx, fp_2);
        /*normalise eigenvector*/
        gsl_matrix_set_col(Norm_Eigvector_ptr, n, s_eigvector_ptr);
    }
    return N;
}

double FD_Bisect(double *Subband_f, double E_fermi_super, double E_fermi_diff, double c_doping_yz, double Meff, int N, int T)
{
    double n_sp = 0, n_s = 0;
    double Hbar=1.054571E-34, Kb = 1.38065E-23;
    int n, m;
    /*temperature set to 300K*/
    double DOS = Meff/(Hbar*Hbar*3.14159);
    double E_fermi_sub = E_fermi_super + E_fermi_diff;
    double E_fermi_mid = (E_fermi_sub + E_fermi_super)/2;
    for (n=0; n<50; n++){
        for (m = 0; m < N; m++){ /*add contributions from all subbands*/
            n_sp = DOS*Kb*T*log(exp((E_fermi_mid - Subband_f[m])/(Kb*T)) + 1);
            n_s = n_sp + n_s;
        }
        if (n_s > c_doping_yz){
            E_fermi_super = E_fermi_mid;
            E_fermi_mid = (E_fermi_super + E_fermi_sub)/2;
        }
        else {
            E_fermi_sub = E_fermi_mid;
            E_fermi_mid = (E_fermi_sub + E_fermi_super)/2;
        }
        n_s = 0;
    }
    printf("nPrecise Fermi energy is %E eV", E_fermi_mid);
    return E_fermi_mid;
}

double FermiDiracPopulate(double *Subband_f, double c_doping_yz, double Meff, double E_fermi, int N, int T)
{
    int m;
    double FL_precise, n_sp = 0, n_s = 0;
    /* initial fermi-energy roughly at intrinsic level (@300K)*/
    double Hbar=1.054571E-34, Kb = 1.38065E-23;
    double d_E_fermi = E_fermi/100;
    double DOS = Meff/(Hbar*Hbar*3.14159);
    for (n_s = 0; E_fermi < 0; E_fermi = E_fermi - d_E_fermi){
        printf("n%E", E_fermi);
        for (m = 0; m < N; m++){ /*add contributions from all subbands*/
            n_sp = DOS*Kb*T*log(exp((E_fermi - Subband_f[m])/(Kb*T)) + 1);
            /*NB. subband energies were calculated to be positive (should be negative)*/
            n_s = n_sp + n_s;
        }
        //printf("ncharge density for E_fermi of %E is %E per meter squared", E_fermi, n_s);
        if (n_s > c_doping_yz){
            printf("nRaw Fermi energy is %E eV", E_fermi);
            FL_precise = FD_Bisect(Subband_f, E_fermi, d_E_fermi, c_doping_yz, Meff, N, T);
            break;
        }
        n_s = 0;
    }
    return FL_precise;
}

double subband_contrib(double E_fermi, double E, double Meff, int T)
{
    double n_sp = 0;
    double Hbar=1.054571E-34, Kb = 1.38065E-23;
    double DOS = Meff/(Hbar*Hbar*3.14159);
    n_sp = DOS*Kb*T*log(exp((E_fermi + E)/(Kb*T)) + 1);
    return n_sp;
}

void Charge_Density(gsl_matrix *Norm_eigenvector_ptr, gsl_vector *s_eigvector_ptr, gsl_vector *Charge_Density_ptr, double *subband_ptr, double E_fermi, double dx, double c_doping_yz, double Meff, int Dim, int N, int T, FILE *fp_3)
{
    int n, m;
    double phi, n_sp, temp;
    double c_doping_xyz, q = 1.602176E-19;
    temp = sqrt(c_doping_yz);
    c_doping_xyz = temp*temp*temp;
    for (n=0; n<N; n++)
    {
        gsl_matrix_get_col(s_eigvector_ptr, Norm_eigenvector_ptr, n);
        n_sp = subband_contrib(E_fermi, subband_ptr[n], Meff, T);
        printf("ncharge contribution from subband %d is %E coulombs per meter squared", n, n_sp);
        for (m=0; m<Dim; m++)
        {
            phi = gsl_vector_get(s_eigvector_ptr, m);
            temp = gsl_vector_get(Charge_Density_ptr, m) - q*phi*phi*n_sp;
            gsl_vector_set(Charge_Density_ptr, m, temp);
        }
    }
    /*Output*/
    for(m=0; m<Dim; m++){
        //fprintf(fp_3, "n%E", gsl_vector_get(Charge_Density_ptr, m));
        temp = gsl_vector_get(Charge_Density_ptr, m) + q*c_doping_xyz;
        //temp = temp*dx*dx; /*for Poisson Solver*/
        gsl_vector_set(Charge_Density_ptr, m, temp);
        fprintf(fp_3, "n%E", gsl_vector_get(Charge_Density_ptr, m));
    }
}

double charge_density_check(gsl_vector *charge_density_ptr, int Dim)
{
    double charge_tot=0;
    int n;
    for(n=0; n<Dim; n++)
    {
        charge_tot = charge_tot + gsl_vector_get(charge_density_ptr, n);
    }
    return charge_tot;
}

//
//double Poisson_Solve(int Dim, int barrier_R, int barrier_L, double dx, gsl_vector *Main_Diagonal_ptr, gsl_vector *Super_Diagonal_ptr, gsl_vector *Sub_Diagonal_ptr,
//  gsl_vector *Charge_Density_ptr, gsl_vector *Elec_Potential_ptr, FILE *fp_3, FILE *fp_4){
    //
    //int n;
    //
    ///*constants*/
    //double dielectric_Ge = 16, dielectric_Si = 11.9, permitt_0 = 8.854187E-12;
    //double abs_permitt_Ge = dielectric_Ge*permitt_0, abs_permitt_Si = dielectric_Si*permitt_0;
    //
    //
    ///*assign vectors*/
    ////for(n=1; n<=Dim; n++){
        ////gsl_vector_set(Main_Diagonal_ptr, (n-1), -2);
    ////}
    ////for(n=1; n<=(Dim-1); n++){
        ////gsl_vector_set(Super_Sub_Diagonal_ptr, (n-1), 1);
    ////}
    ////for(n=1; n<=Dim; n++){
        ////double m=n/100;
        //////double temp = (1-(2*m*m));
        ////double temp = 1;
        ////gsl_vector_set(Charge_Density_ptr, (n-1), temp );
        /*/************************************************************\
        *
        \************************************************************/
        function 1-2x^2*/
        ////printf("%E", m);
        ////printf("ncharge density at point %E is %E", n, temp);
        ////
    ////}
    //
    //for (n=0; n<(Dim); n++) {
        // if(n>=barrier_R){
            // gsl_vector_set(Main_Diagonal_ptr, n, -2*abs_permitt_Si);
        // }
        // if(n<barrier_R && n>=barrier_L){
            // gsl_vector_set(Main_Diagonal_ptr, n, -2*abs_permitt_Ge);
        // }
        // if (n<barrier_L){
            // gsl_vector_set(Main_Diagonal_ptr, n, -2*abs_permitt_Si);
        // }
    //}
    //for (n=0; n<(Dim-1); n++) {
        // if(n>=(barrier_R-1)){
            // gsl_vector_set(Sub_Diagonal_ptr, n, 1*abs_permitt_Si);
        // }
        // if(n<(barrier_R-1) && n>(barrier_L-2)){
            // gsl_vector_set(Sub_Diagonal_ptr, n, 1*abs_permitt_Ge);
        // }
        // if (n<(barrier_L-1)){
            // gsl_vector_set(Sub_Diagonal_ptr, n, 1*abs_permitt_Si);
        // }
    //}
    //for (n=0; n<(Dim-1); n++) {
        // if(n>=barrier_R){
            // gsl_vector_set(Super_Diagonal_ptr, n, 1*abs_permitt_Si);
        // }
        //  if(n<barrier_R && n>=barrier_L){
            // gsl_vector_set(Super_Diagonal_ptr, n, 1*abs_permitt_Ge);
        // }
        //  if (n<barrier_L){
            // gsl_vector_set(Super_Diagonal_ptr, n, 1*abs_permitt_Si);
        // }
    //}
    //
    //
    /////*solve matrix*///
    //gsl_linalg_solve_tridiag(Main_Diagonal_ptr, Super_Diagonal_ptr, Sub_Diagonal_ptr, Charge_Density_ptr, Elec_Potential_ptr);
    //
    ////gsl_linalg_solve_symm_tridiag(Main_Diagonal_ptr, Super_Sub_Diagonal_ptr, Charge_Density_ptr, Potential_ptr);
    //
    //
    ///*output potential & charge distribution*/
    //for(n=0; n<Dim; n++){
        //fprintf(fp_4, "n%E", gsl_vector_get(Elec_Potential_ptr, n));
    //}
    //for(n=0; n<Dim; n++){
        //fprintf(fp_3, "n%E", gsl_vector_get(Charge_Density_ptr, n));
    //}
    //return 0;
//}
          
int main()
{
    int n=0, N=0;
    int Dim =800;
    int barrier_L=200, barrier_R=600;
    double dx=1E-10, E_well = -1.6022E-20, Meff_barrier = 1.73078256E-31, Meff_well = 6.1032E-32; /* Meff for Ge is [3.7348E-32]; Meff for GaAs is [6.1032E-32]; Meff for Si is [1.73078256E-31]*/
    double E_fermi, T=100;
    double c_doping_yz = 1E16;
    FILE *fp_1, *fp_2;
    FILE *fp_3, *fp_4;
    fp_1=fopen("//home//chuck/Documents//Project Code//MM Schrodinger Poisson Solver 3.1//eigenvalues.txt", "w+");
    fp_2=fopen("//home//chuck/Documents//Project Code//MM Schrodinger Poisson Solver 3.1//eigenvectors.txt", "w+");
    fp_3=fopen("//home//chuck/Documents//Project Code//MM Schrodinger Poisson Solver 3.1//charge_distribution.txt", "w+");
    fp_4=fopen("//home//chuck/Documents//Project Code//MM Schrodinger Poisson Solver 3.1//elec_potential.txt", "w+");
    double *Subband_ptr;
    Subband_ptr = calloc(Dim, sizeof(double));
    //double *Subband_prev_ptr;
    //Subband_prev_ptr = calloc(Dim, sizeof(double));
    double *Meff_ptr;
    Meff_ptr = calloc((Dim+1), sizeof(double));
    gsl_eigen_symmv_workspace *Workspace = gsl_eigen_symmv_alloc(Dim);
    gsl_matrix *Hamil_ptr = gsl_matrix_alloc(Dim, Dim);
    gsl_vector *Eigvalue_ptr = gsl_vector_alloc(Dim);
    gsl_matrix *Eigvector_ptr = gsl_matrix_alloc(Dim, Dim);
    gsl_vector *Potential_ptr = gsl_vector_alloc(Dim);
    gsl_matrix *Norm_Eigvector_ptr = gsl_matrix_alloc(Dim, Dim);
    gsl_vector *S_Eigvector_ptr = gsl_vector_alloc(Dim);
    gsl_vector *Charge_Density_ptr = gsl_vector_alloc(Dim);
    //gsl_vector *Main_Diagonal_ptr = gsl_vector_alloc(Dim);
    //gsl_vector *Super_Diagonal_ptr = gsl_vector_alloc(Dim-1);
    //gsl_vector *Sub_Diagonal_ptr = gsl_vector_alloc(Dim-1);
    ////gsl_vector *Super_Sub_Diagonal_ptr = gsl_vector_alloc(Dim-1);
    //gsl_vector *Elec_Potential_ptr = gsl_vector_alloc(Dim);
    /*initialise potential well*/
    for (n=0; n<(Dim); n++) {
        if(n>=barrier_R){
            gsl_vector_set(Potential_ptr, n, 0);
        }
        if(n<barrier_R && n>=barrier_L){
            gsl_vector_set(Potential_ptr, n, E_well);
        }
        if (n<barrier_L){
            gsl_vector_set(Potential_ptr, n, 0);
        }
    }
    /*initialise effective mass*/
    for (n=0; n<(Dim+1); n++) {
        if(n>barrier_R && n<Dim){
            Meff_ptr[n]=Meff_barrier;
        }
        if(n<(barrier_R-1) && n>barrier_L){
            Meff_ptr[n]=Meff_well;
        }
        if (n<(barrier_L) && n>0){
            Meff_ptr[n] = Meff_barrier;
        }
        if (n==0 || n==Dim){
            Meff_ptr[n]= Meff_barrier/2;
        }
        if (n==barrier_R || n==(barrier_R-1) || n==barrier_L || n==(barrier_L-1)){
            Meff_ptr[n]= (Meff_barrier+Meff_well)/2;
        }
    }
    N = Schrodinger_Solve(dx, Dim, barrier_L, barrier_R, E_well, Meff_ptr, Subband_ptr, Workspace, Hamil_ptr, Eigvalue_ptr, Eigvector_ptr, S_Eigvector_ptr, Potential_ptr, Norm_Eigvector_ptr, fp_1, fp_2);
    E_fermi = FermiDiracPopulate(Subband_ptr, c_doping_yz, Meff_well, E_well, N, T);
    Charge_Density(Norm_Eigvector_ptr, S_Eigvector_ptr, Charge_Density_ptr, Subband_ptr, E_fermi, dx, c_doping_yz, Meff_well, Dim, N, T, fp_3);
    double Charge_tot = charge_density_check(Charge_Density_ptr, Dim);
    printf("n%E", Charge_tot);
    //Poisson_Solve(Dim, barrier_R, barrier_L, dx, Main_Diagonal_ptr, Super_Diagonal_ptr, Sub_Diagonal_ptr, Charge_Density_ptr, Elec_Potential_ptr, fp_3, fp_4);
    //Convergence_Check()
    //Subband_store()
    /*free memory*/
    gsl_matrix_free(Eigvector_ptr);
    gsl_vector_free(Eigvalue_ptr);
    gsl_matrix_free(Hamil_ptr);
    gsl_eigen_symmv_free(Workspace);
    gsl_matrix_free(Norm_Eigvector_ptr);
    gsl_vector_free (S_Eigvector_ptr);
    free(Subband_ptr);
    free(Meff_ptr);
    gsl_vector_free(Charge_Density_ptr);
    //gsl_vector_free (Main_Diagonal_ptr);
    //gsl_vector_free (Super_Diagonal_ptr);
    //gsl_vector_free(Sub_Diagonal_ptr);
    //gsl_vector_free(Elec_Potential_ptr);
    //free(Subband_prev_ptr)
    return 0;
}