#include "Itmean.h"

int main() { //Test the IterativeMean family!
	double foo[] = {1.1,0.2,0.3,1};	
	vector<double> bar(foo, foo + sizeof(foo) / sizeof(double) );
	ClocksAtSea boo(bar);
	AGM loo(bar);
	for(int i=0; i<bar.size(); i++) {
		cout << boo[i] << ' ';
	}
	cout << '\n' << boo.calculate() << ' ' << loo.calculate() << '\n';
	cout << boo.iterations() << '\n';
	return 0;
}