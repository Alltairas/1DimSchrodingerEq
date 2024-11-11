#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>
#include <limits>
#include <filesystem>
#include <cstring>
#include "time.h"

using namespace std;
#include "schrodinger.h"
int main()
{

	schrodinger A;
	A.valeurdonne();
	
	//création d'un dossier pour mettre les résultats
	string rep="mkdir "+A.cheminfichier;
	const char *cmd = rep.c_str();
   	system(cmd); 
	
	A.suitepsi();
	A.coef2();
	A.methode2();
	
	A.enregistre();
	A.coefTenfonctiondeEeta();
	A.resultat();
	return 0;

}
