//###################################################################################################################################//
//###################################################################################################################################//
//*************************************   NUMERICAL SOLUTION FOR BLASIUS EQUATION    ************************************************//
//*************************************   AUTHOR: PRAKASH (AE15D015)                 ************************************************//
//*************************************   DATE  : OCTOBER 22, 2015                   ************************************************//
//*************************************   LAST MODIFIED: October 24,2015             ************************************************//
//###################################################################################################################################//
//###################################################################################################################################//

// In this program f = f; f'= g1; f'' = g2; f''' = g3;

#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;

#define N 32

int i, j;
long double z, z_older, z_old, z_new, del_x;
long double eta = 10.0, g1_old_inf, g1_older_inf;
long double x[N+1], f[N+1], g1[N+1], g2[N+1];
long double k1_f[N+1], k2_f[N+1], k3_f[N+1], k4_f[N+1];
long double k1_g1[N+1], k2_g1[N+1], k3_g1[N+1], k4_g1[N+1];
long double k1_g2[N+1], k2_g2[N+1], k3_g2[N+1], k4_g2[N+1];
long double* Runge_Kutta(long double);
void plt();

//************ Function calculating f, g1, g2 using fourth order RK ***************************//
long double* Runge_Kutta(long double z){
	for(i=0; i<N; i++){
		f[0] = 0.0; g1[0] = 0.0; g2[0] = z;
		k1_f[i] = del_x*g1[i];
		k1_g1[i] = del_x*g2[i];
		k1_g2[i] = del_x*(-0.5*f[i]*g2[i]);

		k2_f[i] = del_x*(g1[i] + (0.5*k1_g1[i]));
		k2_g1[i] = del_x*(g2[i] + (0.5*k1_g2[i]));
		k2_g2[i] = del_x*(-0.5*(f[i] + (0.5*k1_f[i]))*(g2[i] + (0.5*k1_g2[i])));

		k3_f[i] = del_x*(g1[i] + (0.5*k2_g1[i]));
		k3_g1[i] = del_x*(g2[i] + (0.5*k2_g2[i]));
		k3_g2[i] = del_x*(-0.5*(f[i] + (0.5*k2_f[i]))*(g2[i] + (0.5*k2_g2[i])));
		
		k4_f[i] = del_x*(g1[i] + k3_g1[i]);
		k4_g1[i] = del_x*(g2[i] + k3_g2[i]);
		k4_g2[i] = del_x*(-0.5*(f[i] + k3_f[i])*(g2[i] + k3_g2[i]));

		f[i+1] = f[i] + ((1.0/6.0)*(k1_f[i] + 2.0*k2_f[i] + 2.0*k3_f[i] + k4_f[i]));
		g1[i+1] = g1[i] + ((1.0/6.0)*(k1_g1[i] + 2.0*k2_g1[i] + 2.0*k3_g1[i] + k4_g1[i]));
		g2[i+1] = g2[i] + ((1.0/6.0)*(k1_g2[i] + 2.0*k2_g2[i] + 2.0*k3_g2[i] + k4_g2[i]));
	}
	
	return f;
	return g1;
	return g2;
}

int main(){

	ofstream myfile;
	myfile.open ("Blasius.dat",ios::out);
	myfile.flags( ios::dec | ios::scientific);
	myfile.precision(7);

	/*******  Grid generation ***********/
	del_x = eta/N;
	for(i=0; i<=N; i++) x[i] = i*del_x;

	/******** Initial calculation with guess values *******/
	z_older = 0.2;  
	Runge_Kutta(z_older);
	g1_older_inf = g1[N];

	z_old = 0.25;
	Runge_Kutta(z_old);
	g1_old_inf = g1[N];

	z_new = z_old - (((g1_old_inf - 1.0)*(z_old - z_older))/(g1_old_inf - g1_older_inf)); //NR-Method

	/***************** Loop for z_new convergence **************************/
	while((z_new - z_old)>1E-15){
		z_older = z_old;
		z_old = z_new;		
	
		f[0] = 0.0; g1[0] = 0.0; g2[0] = z_older; 
		Runge_Kutta(z_older);
		g1_older_inf = g1[N];
	
		f[0] = 0.0; g1[0] = 0.0; g2[0] = z_old;
		Runge_Kutta(z_old);
		g1_old_inf = g1[N];

		z_new = z_old - (((g1_old_inf - 1.0)*(z_old - z_older))/(g1_old_inf - g1_older_inf)); //Newton-Raphson Method
				
	}
	
	cout<<"f'''(0)="<<z_new<<endl;
	cout.flags( ios::dec | ios::scientific);
	cout.precision(7);
	
	for(i=0; i<=N; i++){
		myfile<<x[i]<<"\t"<<f[i]<<"\t"<<g1[i]<<"\t"<<g2[i]<<endl;
	}
	myfile.close();
	plt();

}

/*********** Function for plotting using gnuplot ****************/
void plt(){
	FILE *plot = popen("gnuplot -persist","w");  
	fprintf(plot, "set title 'Blasius solution'\n");
	fprintf(plot,"set xlabel 'eta'\n");
	fprintf(plot,"set ylabel 'u/Uinf'\n");
	fprintf(plot,"set xrange [0:7] \n");
	fprintf(plot,"set yrange [0:1.0] \n");
	fprintf(plot, "set grid xtics \n");
	fprintf(plot, "set grid ytics \n");
	fprintf(plot, "plot 'Blasius.dat' u 1:2 title 'f' w lp, 'Blasius.dat' u 1:3 title 'g1' w lp, 'Blasius.dat' u 1:4 title 'g2' w lp \n");
	fclose(plot);
}

	








