#include<stdio.h>
#include<math.h>
#include<iostream>
using namespace std;

//==============================================================//

#define G 9.8
#define E 0.5

//==============================================================//

double maxx(double a, double b){
    if(a > b){
        return a;
    }else{
        return b;
    }
}

//==============================================================//

double minn(double a, double b){
    if(a < b){
        return a;
    }else{
        return b;
    }
}

//==============================================================//

int main()
{

    //==============================================================//

    /* variable yang dibutuhkan */
    double L            = 100;
    double dx           = 0.1;
    int    Nx           = L/dx;

    double x_nol        = 19.85;
    double H            = 0.3;
    double h_nol        = 1;
    double t_nol        = sqrt(G/h_nol);
    double x_satu       = x_nol + ( (1/sqrt(3*H/(4*h_nol))) * acosh(sqrt(20)) );

    double Tfinal       = 30/t_nol;
    //double Tfinal       = 25/t_nol;
    //double Tfinal       = 20/t_nol;
    //double Tfinal       = 15/t_nol;
    double dt           = 0.01;
    int    Nt           = L/dt;

    double energi       = 0;
    double energy_head  = 0;

    double* x           = new double[Nx+1];

    double* h_new       = new double[Nx+1];
    double* h_old       = new double[Nx+1];
    double* h_bin_old   = new double[Nx+1];
    double* h_bin_new   = new double[Nx+1];
    double* h_bar       = new double[Nx+1];

    double* q           = new double[Nx+1];
    double* q_setengah  = new double[Nx+1];
    double* q_bar       = new double[Nx+1];

    double* eta_new     = new double[Nx+1];
    double* eta_old     = new double[Nx+1];
    double* eta_head    = new double[Nx+1];

    double* d           = new double[Nx+1];

    double* unew        = new double[Nx+1];
    double* uold        = new double[Nx+1];
    double* ubin        = new double[Nx+1];

    double* shock       = new double[Nx+1];

    //==============================================================//

    /* initial state */
    for(int i= 0; i<=Nx; i++){

        x[i] = ((i*dx) - 10) ;

        if(x[i] <= x_nol){
            d[i] = -x[i]/x_nol;
        }else {
            d[i] = -1;
        }

        eta_old[i]  = maxx( (H/h_nol* ( pow( (1/(cosh( (sqrt(3*H/(4*h_nol))) * (x[i]-x_satu) ))) , 2 ) ) ) , d[i] );
        eta_new[i]  = eta_old[i];

        h_old[i]    = eta_old[i] - d[i];
        h_new[i]    = h_old[i];

        uold[i]     = -1*(( (sqrt(G*(H+h_nol))) * eta_old[i] ) / ( 1+eta_old[i] ));
        unew[i]     = uold[i];

    }

    for(int k=0; k<=Nx; k++)
    {
        shock[k] = ( (h_new[k+1]*unew[k+1]) + (h_new[k-1]*unew[k-1]) ) / (h_new[k-1]+h_new[k+1]);
    }

    //==============================================================//

    /* buat file untuk menyimpan initial output */
    FILE* fp = fopen("Initial-Output1D.dat","w");
    for(int k=0; k<=Nx; k++)
    {
        fprintf(fp,"%g %g %g\n",x[k],eta_old[k],d[k]);
    }
    fclose(fp);

    //==============================================================//

    /* buat file untuk menyimpan energi */
    FILE* fp2 = fopen("output_energi_runup_stelling.dat","w");

    //==============================================================//

    /* buat file untuk menyimpan energi */
    FILE* fp3 = fopen("output_energi_head.dat","w");

    //==============================================================//

    /* buat file gambar untuk simulasi */
    int    nama_shock = 1;
    double time       = 0;

    FILE* gnuplot1=popen("C://gnuplot/bin/gnuplot.exe","w");
    fprintf(gnuplot1,"reset\n");
    fprintf(gnuplot1,"set title 'time = %g'\n",time);
    fprintf(gnuplot1,"set yrange[-2:0.3]\n");
    fprintf(gnuplot1,"set xrange[0:20]\n");
    fprintf(gnuplot1,"set term png\n");
    fprintf(gnuplot1,"set output \"hasil_shock/graph%d.png\"\n",nama_shock);
    fprintf(gnuplot1,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot1,"%g %g\n",x[i],shock[i]);
    }
    fprintf(gnuplot1,"end\n");

    nama_shock++;

    //==============================================================//

    /* buat file gambar untuk simulasi */
    int    nama = 1;

    FILE* gnuplot = popen("C://gnuplot/bin/gnuplot.exe","w");
    fprintf(gnuplot,"reset\n");
    fprintf(gnuplot,"set title 'time = %g'\n",time);
    fprintf(gnuplot,"set yrange[-0.2:0.5]\n");
    fprintf(gnuplot,"set xrange[-10:20]\n");
    fprintf(gnuplot,"set term png\n");
    fprintf(gnuplot,"set output \"hasil/graph%d.png\"\n",nama);
    fprintf(gnuplot,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot,"%g %g\n",x[i],eta_new[i]);
    }
    fprintf(gnuplot,"end\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot,"%g %g\n",x[i],d[i]);
    }
    fprintf(gnuplot,"end\n");

    nama++;
    cout << "time : " << time << "\n";

    //==============================================================//

    /* jalankan simulasi */
    while(time < Tfinal){

        //==============================================================//

        time += dt;

        //==============================================================//

        /* langkah 1 : loop untuk h_bin_old */
        for(int k=1; k<=Nx-1; k++)
        {
            if(uold[k] >= 0){
                h_bin_old[k] = h_old[k];
            }else if(uold[k] < 0){
                h_bin_old[k] = h_old[k+1];
            }else{
                h_bin_old[k] = maxx(eta_old[k],eta_old[k+1]) + minn(d[k], d[k+1]);
            }
        }

        //==============================================================//

        /* langkah 2 : loop untuk q */
        for(int k=1; k<=Nx-1; k++)
        {
            q[k] = h_bin_old[k] * uold[k];
        }

        /* Boundary condition */
        q[0]    = 0;
        q[Nx]   = q[Nx-1];

        //==============================================================//

        /* langkah 3 : loop untuk h_new */
        for(int k=1; k<=Nx; k++)
        {
            h_new[k] = h_old[k] - ((dt/dx)*(q[k] - q[k-1]));
        }

        /* update h_old */
        for(int k=1; k<=Nx; k++)
        {
            h_old[k] = h_new[k];
        }

        //==============================================================//

        /* langkah 4 : loop untuk h_bar */
        for(int k=1; k<=Nx-1; k++)
        {
            h_bar[k] = 0.5*(h_new[k] + h_new[k+1]);
        }

        //==============================================================//

        /* langkah 5 : loop untuk eta_head sebelum menuju safiro filter */
        for(int k=1; k<=Nx; k++)
        {
            //eta_head[k] = h_new[k] + d[k];
            eta_new[k] = h_new[k] + d[k];
        }

        /* Boundary condition menggunakan safiro atau tanpa safiro*/
        eta_head[0]     = eta_head[1];
        eta_head[Nx]    = eta_head[Nx-1];

        //==============================================================//

        /* langkah 6 : loop untuk eta_new safiro filter */
        for(int k=1; k<=Nx-1; k++)
        {
            //eta_new[k] = ((1-E)*eta_head[k]) + (0.5*E*(eta_head[k-1]+eta_head[k+1]));
        }

        /* Boundary condition menggunakan safiro maupun tanpa safiro*/
        eta_new[0]   = eta_new[1];
        //eta_new[Nx]  = eta_new[Nx-1];

        /* update eta_old */
        for(int k=1; k<=Nx; k++)
        {
            eta_old[k] = eta_new[k];
        }

        //==============================================================//

        /* langkah 7 : loop untuk q_bar */
        for(int k=1; k<=Nx; k++)
        {
            q_bar[k] = 0.5 * (q[k-1] + q[k]);
        }

        //==============================================================//

        /* langkah 8 : loop untuk ubin */
        for(int k=1; k<=Nx; k++)
        {
            if(q_bar[k] >= 0){
                ubin[k] = uold[k-1];
            }else{
                ubin[k] = uold[k];
            }
        }

        //==============================================================//

        /* langkah 9 : loop untuk unew */
        for(int k=1; k<=Nx-1; k++)
        {
            if(h_bar[k] == 0){
                unew[k] = 0;
            }else{
                unew[k] = uold[k] - (dt/dx)*(1/h_bar[k])*( (q_bar[k+1]*ubin[k+1]) - (q_bar[k]*ubin[k]) - (uold[k]*(q_bar[k+1] - q_bar[k])) ) - (G*(dt/dx)*(eta_new[k+1]-eta_new[k]));
            }
        }

        /* Boundary condition */
        unew[0]     = 0;
        unew[Nx]    = 0;

        /* Update u_old */
        for(int k=0; k<=Nx; k++)
        {
            uold[k] = unew[k];
        }

        //==============================================================//

        for(int k=0; k<=Nx; k++)
        {
            shock[k] = ( (h_new[k+1]*unew[k+1]) + (h_new[k-1]*unew[k-1]) ) / (h_new[k-1]+h_new[k+1]);
        }

        //==============================================================//

        /* buat file gambar untuk simulasi */
        if(nama_shock%100==0){

            FILE* gnuplot1=popen("C://gnuplot/bin/gnuplot.exe","w");
            fprintf(gnuplot1,"reset\n");
            fprintf(gnuplot1,"set title 'time = %g'\n",time);
            fprintf(gnuplot1,"set yrange[-2:0.3]\n");
            fprintf(gnuplot1,"set xrange[0:20]\n");
            fprintf(gnuplot1,"set term png\n");
            fprintf(gnuplot1,"set output \"hasil_shock/graph%d.png\"\n",nama);
            fprintf(gnuplot1,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
            for(int i=0; i<=Nx; i++){
                fprintf(gnuplot1,"%g %g\n",x[i],shock[i]);
            }
            fprintf(gnuplot,"end\n");
        }

        //==============================================================//

        /* buat file gambar untuk simulasi */
        if(nama%100==0){

            FILE* gnuplot=popen("C://gnuplot/bin/gnuplot.exe","w");
            fprintf(gnuplot,"reset\n");
            fprintf(gnuplot,"set title 'time = %g'\n",time);
            fprintf(gnuplot,"set yrange[-0.2:0.5]\n");
            fprintf(gnuplot,"set xrange[-10:20]\n");
            fprintf(gnuplot,"set term png\n");
            fprintf(gnuplot,"set output \"hasil/graph%d.png\"\n",nama);
            fprintf(gnuplot,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
            for(int i=0; i<=Nx; i++){
                fprintf(gnuplot,"%g %g\n",x[i],eta_new[i]);
            }
            fprintf(gnuplot,"end\n");
            for(int i=0; i<=Nx; i++){
                fprintf(gnuplot,"%g %g\n",x[i],d[i]);
            }
            fprintf(gnuplot,"end\n");
        }

        //==============================================================//

        nama++;
        nama_shock++;
        cout << "time : " << time << "\n";

        //==============================================================//


        /* print energi ke dalam file */
        for(int k=0; k<=Nx; k++)
        {
            energi = energi + ( (0.5*h_new[k]*unew[k]*unew[k]) + (0.5*G*h_new[k]*h_new[k]) + (G*h_new[k]*d[k]) );
        }
        for(int k=1; k<=Nx-1; k++)
        {
            energy_head = energy_head + ( (unew[k]*unew[k]/(2*G))+eta_new[k] );
        }
        fprintf(fp2,"%g %g\n",time,energi);
        fprintf(fp3,"%g %g\n",time,energy_head);

        //==============================================================//

    }

    //==============================================================//

    gnuplot1=popen("C://gnuplot/bin/gnuplot.exe","w");
    fprintf(gnuplot1,"reset\n");
    fprintf(gnuplot1,"set title 'time = %g'\n",time);
    fprintf(gnuplot1,"set yrange[-2:0.3]\n");
    fprintf(gnuplot1,"set xrange[0:20]\n");
    fprintf(gnuplot1,"set term png\n");
    fprintf(gnuplot1,"set output \"hasil_shock/graph%d.png\"\n",nama_shock);
    fprintf(gnuplot1,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot1,"%g %g\n",x[i],shock[i]);
    }
    fprintf(gnuplot1,"end\n");

    nama_shock++;

    gnuplot=popen("C://gnuplot/bin/gnuplot.exe","w");
    fprintf(gnuplot,"reset\n");
    fprintf(gnuplot,"set title 'time = %g'\n",time);
    fprintf(gnuplot,"set yrange[-0.2:0.5]\n");
    fprintf(gnuplot,"set xrange[-10:20]\n");
    fprintf(gnuplot,"set term png\n");
    fprintf(gnuplot,"set output \"hasil/graph%d.png\"\n",nama);
    fprintf(gnuplot,"plot '-' with lines lw 3 lc rgb 'blue', '-' with lines lw 3 lc rgb 'black'\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot,"%g %g\n",x[i],eta_new[i]);
        cout << "time : " << time << "\n";
    }
    fprintf(gnuplot,"end\n");
    for(int i=0; i<=Nx; i++){
        fprintf(gnuplot,"%g %g\n",x[i],d[i]);
    }
    fprintf(gnuplot,"end\n");

    //==============================================================//


    //==============================================================//

    /* buat file untuk menyimpan hasil output */
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t15.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t20.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t25.dat","w");
    FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t30.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t15_safiro.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t20_safiro.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t25_safiro.dat","w");
    //FILE* fp1=fopen("Output1D_Runup_SkemaStelling_t30_safiro.dat","w");
    for(int k=0; k<=Nx; k++)
    {
        fprintf(fp1,"%g %g %g\n",x[k],eta_new[k],d[k]);
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    //==============================================================//

    return 0;
}

