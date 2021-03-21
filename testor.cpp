#include "Itmean.h"

int main() { //Test the IterativeMean family!
	double foo[] = {1.1,0.2,0.3,1};	
	vector<double> bar(foo, foo + sizeof(foo) / sizeof(double) );
	ClocksAtSea boo(&bar);
	AGM loo(&bar);
	Geothdian xoo(&bar);
	for(int i=0; i<bar.size(); i++) {
		cout << boo[i] << ' ';
	}
	cout << '\n' << boo.calculate() << ' ' << loo.calculate() << ' ' << xoo.calculate() << '\n';
	cout << boo.iterations() << '\n';
	//Confirm the result of https://xkcd.com/2435/
	double moo[] = {1.0,1.0,2.0,3.0,5.0};
	vector<double> lar(moo, moo + sizeof(moo) / sizeof(double) );
	Geothdian xconfirm(&lar);
	cout << xconfirm.calculate() << '\n';
	return 0;
}
