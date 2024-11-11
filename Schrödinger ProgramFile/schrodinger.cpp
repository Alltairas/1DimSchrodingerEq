#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <filesystem>
#include <cstring>
#include "time.h"

using namespace std;


#include "schrodinger.h"

/***********************************************************************************************
			Fonction de la class schrodinger
***********************************************************************************************/

void schrodinger::resultat()//écrit les résultats des coefficient dans un fichier
{

	string nomfichier;
	valeurdonne();
	nomfichier = cheminfichier+"/resultat.txt";

	ofstream fichier2;
	fichier2.open(nomfichier.c_str());
	
	fichier2 <<"E= "<<E<<endl;
	fichier2 <<"m= "<<m<<endl;
	fichier2<<"a= "<<a<<endl;
	fichier2<<"N= "<<N<<endl;
	fichier2 <<"V0= "<<V0<<endl;
	
	
	
	
	//barière de potentiel
	if( choixdef==1){
		fichier2<<" 1: barriere de potentiel " <<endl;
		fichier2<<"écart théorique numérique =  "<<abs(rtheo-Rmeth1)<<endl;
		fichier2<<"écart théorique numérique par euler est =  "<<abs(rtheo-Rmeth2)<<endl;
	}
	//gaussienne
	else if(choixdef==2){
		fichier2<<" 2: gaussienne " <<endl;
		fichier2<<"sigma= "<<sigma<<endl;
		fichier2 <<"x0gauss= "<<x0gauss<<endl;
		
	}
	//2 barière de potentiel
	else if(choixdef==3){	
		fichier2<<" 3: 2 barrères de potentiel " <<endl;
		
	}
	//3 barières de potentiel	
	else if(choixdef==4){
		fichier2<<"  4: 3 barrière de potentiel " <<endl;
		
	}
	//barrière triangulaire 
	else if(choixdef==5){
		fichier2<<"  5: barrière triangulaire " <<endl;
	
		
	}
	//barriere alpha decay
	else if(choixdef==6){
		fichier2 <<"λ= "<<lambda<<endl;
		fichier2<<" 6:barriere décroissance Alpha"<<endl;
		
		
	}
	//barriere Electron trap avec une infinité de barrière potentiel 
	else{
		fichier2<<" 7:barriere electron trap"<<endl;
		
		
	}
	fichier2.precision(15);
	fichier2 <<"R calculé numériquement = "<<Rmeth1<<endl;
	fichier2 <<"T calculé numériquement = "<<Tmeth1<<endl;
	fichier2<<"Le coefficient de reflexion numérique par euler est R= "<<Rmeth2<<endl;
	fichier2<<"Le coefficient de transmission numérique par euler est T= "<<1-Rmeth2<<endl;
	fichier2<<"écart entre les 2 méthodes est =  "<<abs(Rmeth1-Rmeth2)<<endl;
	fichier2<<"energie resonance = "<<result[0]<<" coeff t ="<<result[1]<<"  à la taille "<<result[2]<<endl;

	fichier2.close();
}




void schrodinger :: valeurdonne()
{
	ifstream fichier;
	fichier.open("valeur.txt");
	fichier>>cheminfichier>>choixdef>>E>>m>>a>>N>>V0>>lambda;
	fichier.close();

	sigma=(1.0/(V0*sqrt(2.0*M_PI)));
        x0gauss=a/2;
	k=sqrt(2.0*m*E)/(hbar);
	pas = a/N;
	
} 






double schrodinger :: f(double x)
{    
    //barière de potenciel
    if( choixdef==1 ){

        if(x>0 && x<=a){
            return V0;
        }    
        else{
            return 0;
        }
    }
    //gaussienne
    else if(choixdef==2){
   
        return (1.0/(sigma*sqrt(2*M_PI)))*exp(-((x-x0gauss)*(x-x0gauss))/(2.0*sigma*sigma));
    
    }
    //2 barières de potentiel
    else if(choixdef==3){
        if((x>=0 && x<=(a/3))||(x>=(a*2/3) && x<=a)){
            return V0;
        }    
        else{
            return 0;
        }
    }
    //3 barières de potentiel
    else if(choixdef==4){
        if((x>=0 && x<=(a*20/100))||(x>=(a*40/100) && x<=(a*60/100)) || (x>=(a*80/100) && x<=(a)))
        {
            return V0;
        }
        else{
            return 0;
        }
    }
    //barrière triangulaire    
    else if(choixdef==5){
        return (2 * V0/ M_PI) * asin(sin(2 * (M_PI / (a*2)) * (double)x ));
    }
    //barriere de Alpha decay
    else if(choixdef==6){
	if(x<0){return -V0;
	}else{return (V0/N)*(1/x);
	}
	    	
    }
    //barriere Electron trap avec une infinité de barrière potentiel 
    else {
		int j;
	
	for(int i=0;i<N;i++){
 		if(x>=(i*400.0*(a/N)) && x<= ((i+1)*400.0*(a/N))){j=1;}
		
		i++;
	}
	if(j==1){ j=0;
		return V0;
	}	
	else{

	return 0;
	}
}
}



complex<double> schrodinger::psi(double x)
{
	return exp(1i*k*x);

}


void schrodinger::calculecoeff()
{
	complex<double> psipr=(Psiderive-Psi0)/pas;
	Ameth1=((1i*k*Psi0)+psipr)/(2.0*1i*k);
	Bmeth1=((1i*k*Psi0)-psipr)/(2.0*1i*k);
	Rmeth1=norm(Bmeth1/Ameth1);
	Tmeth1=norm(1.0/Ameth1);
	
}


void schrodinger::suitepsi()//methode 1 pour trouver les coefficient transmission et réflexion
{

	complex<double> psixmpas;
	complex<double> psix=psi(a);
	complex<double> psixppas=psi(a+pas);
	double x;	

	for(int i=0;i<=N;i++){
		x=a-i*pas;
		psixmpas=-psixppas+((2.0*m*(pas*pas))/(hbar*hbar)*(f(x)-E)+2)*psix;
		psixppas=psix;
		psix=psixmpas;
		
		if(i==(N-1)){
			Psiderive=psix;
		}
		
		if(i==N){
			Psi0=psix;
		}
	}
	calculecoeff();
}


void schrodinger::enregistre()//enregistre dans des fichiers la fonction d'ondes pour les deux méthodes
{
	vector<double> T(N+1);
	vector<vector<double>> U(5,T);
	string nomfichier;
	double x;


	for(int i=0;i<=N;i++){
		x=(-a+i*pas);
		
		U[0][i]=x;
		U[1][i]=real(Ameth1*psi(x)+Bmeth1*psi(-x));
		U[2][i]=imag(Ameth1*psi(x)+Bmeth1*psi(-x));
		U[3][i]=real(psi(x)+Bmeth2*psi(-x));
		U[4][i]=imag(psi(x)+Bmeth2*psi(-x));
		
		
	}

	nomfichier="/ComplexPsiValues1.txt";
	ecrire(nomfichier, U,1);
		
	
	complex<double> psixmpas,psix,psixppas,psixmpas2,psix2,psixppas2,C2;
	
  	psix=psi(a);
	psixppas=psi(a+pas);
	
   	psixmpas2=psi(-pas)+Bmeth2*psi(-(-pas));
	psix2=1.0+Bmeth2;

	for(int i=0;i<=N;i++){
		x=a-i*pas;
		U[0][N-i]=x;

		psixmpas=-psixppas+((2.0*m*(pas*pas))/(hbar*hbar)*(f(x)-E)+2)*psix;
		psixppas=psix;
		psix=psixmpas;
		
		U[1][N-i]=real(psix);
		U[2][N-i]=imag(psix);
		
		x=i*pas;		
		
		
		psixppas2=-psixmpas2+((2.0*m*(pas*pas))/(hbar*hbar)*(f(i*pas)-E)+2)*psix2;
		psixmpas2=psix2;
		psix2=psixppas2;
		

		U[3][i]=real(psix2);
		U[4][i]=imag(psix2);

		if(i==N){
			
			C2=psix2*psi(-(a));
			
		}
	}
	
	nomfichier = "/ComplexPsiValues2.txt";
	ecrire(nomfichier, U,1);
	
	
	for(int i=0;i<=N;i++){
		x=a+i*pas;
		U[0][i]=x;
		U[1][i]=real(psi(x));
		U[2][i]=imag(psi(x));
		U[3][i]=real(C2*psi(x));
		U[4][i]=imag(C2*psi(x));
	}


	nomfichier = "/ComplexPsiValues3.txt";
	ecrire(nomfichier, U,1);

	vector<double> T2(3*N+1);
	vector<vector<double>> U2(2,T2);
	for(int i=(-N);i<=(2*N);i++){
		U2[0][i+N]=i*pas;
		U2[1][i+N]=f(i*pas);
		
	
	}
	
	nomfichier ="/potentiel.txt";
	ecrire(nomfichier, U2,2);

}

void schrodinger::ecrire(string nomfichier,vector<vector<double>> U,int n)
{
	ofstream fichier;
	nomfichier=cheminfichier+nomfichier;
	remove(nomfichier.c_str());
	
	fichier.open(nomfichier.c_str());
	if(n==1){
		fichier<<"#Données"<<endl;
		fichier<<"#x"<<"    "<<"Psi real Meth1"<<"    "<<"Psi imag Meth1"<<"    "<<"Psi real Meth2"<<"    "<<"Psi imag Meth2"<<endl;
	
		for(int i=0;i<=N;i++){
		
			fichier<<U[0][i]<<"   		"<<U[1][i]<<"  	 	"<<U[2][i]<<"   		"<<U[3][i]<<"  	 	"<<U[4][i]<<endl;

		}
	}
	if(n==2){
		fichier<<"#Données"<<endl;
		fichier<<"#x"<<"      "<<"Energie Potentiel"<<endl;
		for(int i=0;i<=(3*N);i++){
		
			fichier<<U[0][i]<<"       "<<U[1][i]<<endl;
		
		}

	}
	if(n==3){
		fichier<<"#Données"<<endl;
		fichier<<"#Energie"<<"      "<<"Coefficient de Transmition"<<"      "<<"a"<<"      "<<"Coefficient de Transmition"<<endl;
		for(int i=0;i<=2000;i++){
		
			fichier<<U[0][i]<<"		      "<<U[1][i]<<"		      "<<U[2][i]<<"		      "<<U[3][i]<<endl;
		
		}
		

	}

	fichier.close();

}

void schrodinger::coef2()//calcule des coefficients de transmission et réflexion analytiques.
{
	double unhbar=1./hbar;
        double unhbar2=unhbar*unhbar;
	double r,t;


        //E<V0
        if (E<V0)
        {
                double q2a = sqrt((2*m*V0*a*a*unhbar2)*(1-E/V0));
                double sinh2=sinh(q2a);
                sinh2 = sinh2*sinh2;
                t=(1.0/(1+sinh2*V0*V0/(4*E*(V0-E))));
                r=1-t; 
            
        }
        //E=V0 
        else if (E==V0)
        {
                t= (1/(1+(m*a*a*V0*unhbar2/(2))));
                r=(1-t);
        }
        //E>V0
        else
        {
                double q3a = sqrt((2*m*V0*a*a*unhbar2)*(E/V0-1));
                t=(1/(1+((sin(q3a)*sin(q3a))/(4*E*(E-V0)/(V0*V0)))));
                r=(1-t);
               
        }
	
	rtheo=r;
	ttheo=t;

}

void schrodinger::methode2()//utilise la minimisation pour trouver les coefficients de transmission et de réflexion
{
	minimisation coefB;
	
	double B[2];
	B[0]=0.1;
	B[1]=0.1;
 	coefB.steepestDescent(B);
 	Bmeth2=B[0]+1i*B[1];
	Rmeth2 = norm(Bmeth2);

}




void schrodinger::coefTenfonctiondeEeta()//calcule du coefficient de transmission en fonction de E et de a et recherche du premier maximum
{
	vector<double> V(2001);
	vector<double> res(3);
	vector<vector<double>> W(4,V);
	double max(0),Tmax;
	double NRJ=E/1000;
	int n=0;
	string nomfichier;
	
	
	valeurdonne();

	
	for(int i=0;i<=2000;i++){
		
		E=i*NRJ;	
		suitepsi();
		W[0][i]=E;
		W[1][i]=Tmeth1;

		if(Tmeth1>max){
			max=Tmeth1;
			}
		if(n==0){
			if(Tmeth1<max){
				n=i;
				Tmax=max;
			}
		}
	}

	res[0]=n*NRJ;
	res[1]=Tmax;
	res[2]=a;
	valeurdonne();
	double taille=a/1000;
	
	for(int i=0;i<=2000;i++){
		
		a=i*taille;	
		suitepsi();
		W[2][i]=a;
		W[3][i]=Tmeth1;	
	}
		
	nomfichier="/Coeff_transmission_fonction_Eeta.txt";
	
	ecrire(nomfichier,W,3);
	result=res;

}




/***********************************************************************************************
			Fonction de la class minimisation
***********************************************************************************************/
void minimisation::resolution (double B[])
{	
	complex<double> psixmpas,psix,psixppas,C;
	complex<double> ValB=B[0]+1i*B[1];
	double x;	
	
   	psixmpas=A.psi(-A.pas)+ValB*A.psi(-(-A.pas));
	psix=1.0+ValB;
	
	for(int i=0;i<=A.N;i++){
		x=i*A.pas;
		psixppas=-psixmpas+((2.0*A.m*(A.pas*A.pas))/(A.hbar*A.hbar)*(A.f(x)-A.E)+2)*psix;
		psixmpas=psix;
		psix=psixppas;
	
		if(i==A.N){
			
			C=psix*A.psi(-A.a);
			
		}
	}
	
	psix=C*A.psi(A.a);
	psixppas=C*A.psi(A.a+A.pas);
	
	for(int i=0;i<=A.N;i++){
		x=A.a-i*A.pas;
		psixmpas=-psixppas+((2.0*A.m*(A.pas*A.pas))/(A.hbar*A.hbar)*(A.f(x)-A.E)+2)*psix;
		psixppas=psix;
		psix=psixmpas;
	
		if(i==A.N){
			ValB=psix-1.0;	
		}
	}
	B[0]=real(ValB);
	B[1]=imag(ValB);
	
	
}
double minimisation::fmethode2 (double B[])
{
	
	int i;
	double BN[2] = {B[0],B[1]};
	double n2 = 0;
	resolution(BN);
	for (i=0; i<n; i++)
	{
		n2 += (B[i]-BN[i])*(B[i]-BN[i]);
	}
	return n2;
}

// Fonction qui calcule le gradient, numeriquement, de f en tout point x. En fait pour notre cas n = 2 et x represente B
void minimisation::calculGradient(double x[],  double h, double grad[])
{
	
	double y[n];
	double f0 = fmethode2(x);
	int i;

	for (i=0; i<n; i++) y[i]=x[i];
	for (i=0; i<n; i++)
	{
		y[i] += h;
		grad[i] = (fmethode2(y)-f0)/h;
		y[i] -= h;
	}
}
void minimisation::steepestDescent (double x[])
{	
	A.valeurdonne();
	double h = del;
	double f0 = fmethode2(x);
	double gradf[n];
	double deriv = 0;
	double b;
	double f1;
	int i;

	calculGradient(x, h, gradf);
	

	for (i=0; i<n; ++i) deriv += gradf[i]*gradf[i];

	deriv = sqrt(deriv);
	b = d/deriv;

	while (deriv > del)
	{
		
		for (i=0; i<n; ++i) x[i] -= b*gradf[i];

		h /= 1.5;// Un parametre (1.5) a ajuster en cas de divergence
		calculGradient(x, h, gradf);
		deriv = 0;
		

		for (i=0; i<n; ++i) deriv += gradf[i]*gradf[i];

		deriv = sqrt(deriv);
		b = d/deriv;

		f1 = fmethode2(x);

		if (f1 > f0) d /= 2;
		else f0 = f1;
	}
}







