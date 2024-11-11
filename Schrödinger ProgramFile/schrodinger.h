#ifndef SCHRODINGER_H_INCLUDED
#define SCHRODINGER_H_INCLUDED

class schrodinger
{
	private :
		
		complex<double> Psi0,Psiderive,Bmeth2,Ameth1,Bmeth1;
		double choixdef,sigma,x0gauss,rtheo,ttheo,Rmeth2,lambda,Rmeth1,Tmeth1;
		void ecrire(string nomfichier,vector<vector<double>> U,int n);
		vector<double> result;

	public :
		double hbar=1.0; 
		double E,m,a,N,V0,k,pas;
		string cheminfichier;

		void valeurdonne();
		complex<double> psi(double x);	
		void calculecoeff();
		void suitepsi();
		double f(double x);
		void coef2();
		void enregistre();
		void coefTenfonctiondeEeta();
		void methode2();
		void resultat();
};

class minimisation
{
//minimisation de 
	private:
		double del=0.0001;
		double d=1000;
		int n=2;
		schrodinger A;
		
		void resolution (double B[]);
		double fmethode2 (double B[]);
		void calculGradient(double x[],  double h, double grad[]);
	public:
		
		void steepestDescent (double x[]);
	


};



#endif

