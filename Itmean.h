#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <numeric>
#include <valarray>
using namespace std;

class IterativeMean { //Creates a class framework for building child means using recursive operations
	protected:
		IterativeMean() {} //Iterative mean is not a complete class and should not be instantiated.
		vector<double> raw; //Data to be averaged is stored here.
		unsigned int iteration=0; //Number of iterations gone through
		double precision=0.000001; //Precision to which mean should be calculated.
		double mean; //The working mean, updated each iteration.
		virtual double step_average() {} //This function should be defined to represent a single averaging step in the derived class.
		virtual void reset() { return; } //This function should be defined to prepare the derived class for averaging the data from the first iteration.
	public:
		void set_precision(double precis)
			{ precision=precis; }
		double calculate();
		double calculate(unsigned int iteration_depth);
		double iterate()
			{ iteration++; double next_mean=step_average(); double diff=abs(mean-next_mean); mean=next_mean; return diff; }
		unsigned int size()
			{ return raw.size(); }
		double operator[] (unsigned int i)
			{ return raw[i]; }
		unsigned int iterations()
			{ return iteration; }
};

double IterativeMean::calculate() {
	reset();
	double x=2*precision;
	while(x>precision) {
		x=iterate();
	}
	return mean;
}

double IterativeMean::calculate(unsigned int iteration_depth) {
	reset();
	while(iteration<iteration_depth) {
		iterate();
	}
	return mean;
}


class AGM:public IterativeMean { //The Arithmetic-Geometric Mean is an iterative classic!
	private:
		vector<double> a;
	protected:
		virtual double step_average() //A single step returns the arithmetic mean of the two means from the last step (arithmetic and geometric).
			{ double a1_next=arithmetic_mean(a); a[0]=geometric_mean(a); a[1]=a1_next; return a1_next; }
			//This could be made faster for cases with a zero by returning the geometric mean.
		virtual void reset() { a[0]=geometric_mean(raw); a[1]=arithmetic_mean(raw); } //The initial means are taken from the data.
	public:
		AGM(vector<double> data):a(2,0) { raw=data; }
		double arithmetic_mean(vector<double> data);
		double geometric_mean(vector<double> data);
};

double AGM::arithmetic_mean(vector<double> data) {
	double sum_all=accumulate(data.begin(),data.end(),0.0);
	double mean=sum_all/data.size();
	return mean;
}

double AGM::geometric_mean(vector<double> data) {
	double product_all=accumulate(data.begin(),data.end(),1.0,multiplies<double>());
	double mean=pow(product_all,(1.0/data.size()));
	return mean;
}


class ClocksAtSea:public IterativeMean { //An experimental mean based on the idea of collective indication, e.g. time by clocks at sea.
	private:
		vector<double> set_times;
	protected:
		virtual double step_average();
		virtual void reset() { set_times.assign(raw.size(),0); }
	public:
		ClocksAtSea(vector<double> data) { raw=data; }
};

double ClocksAtSea::step_average() { //Calcluate the arithmetic mean of all but the furthest outlying set_times element, and set that element to the new average.
	//Since step_times advances at each iteration, this mean is divided by total iterations before being passed back.
	for(unsigned int i=0;i<set_times.size();i++) {
		set_times[i]+=raw[i];
	}
	double sum_all=accumulate(set_times.begin(),set_times.end(),0.0);
	double naivmean=sum_all/set_times.size(); //Calculate the arithmetic mean of all set_times; this is the naive mean.
	valarray<double> off_times(set_times.size());
	for(unsigned int i=0;i<set_times.size();i++) { //Calculate the difference of the elements of set_times from the naive mean
		off_times[i]=set_times[i]-naivmean;
	}
	//Find the greatest difference from the naive mean, remove it from the average, and set the corresponding element to the new average.
	double new_mean;
	if(abs(*max_element(set_times.begin(),set_times.end())-naivmean)>abs(*min_element(set_times.begin(),set_times.end())-naivmean)) {
		new_mean=(naivmean*set_times.size()-*max_element(set_times.begin(),set_times.end()))/(set_times.size()-1);
		*max_element(set_times.begin(),set_times.end())=new_mean;
	}
	else if(abs(*max_element(set_times.begin(),set_times.end())-naivmean)<abs(*min_element(set_times.begin(),set_times.end())-naivmean)) {
		new_mean=(naivmean*set_times.size()-*min_element(set_times.begin(),set_times.end()))/(set_times.size()-1);
		*min_element(set_times.begin(),set_times.end())=new_mean;
	}
	else { //If there are two furthest elements equidistant in opposite directions, reset them both.
		new_mean=naivmean;
		*max_element(set_times.begin(),set_times.end())=new_mean;
		*min_element(set_times.begin(),set_times.end())=new_mean;
	}
	return new_mean/iteration;
}