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
		vector<double> raw; //Data to be averaged.
		unsigned int iteration=0; //Number of iterations gone through
		double precision=0.000001; //Precision to which mean should be calculated.
		double mean; //The working mean, updated each iteration.
		virtual double step_average()=0; //This function should be defined to represent a single averaging step in the derived class.
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
		double arithmetic_mean(vector<double> data); //Define these utility means for derived classes
		double geometric_mean(vector<double> data);
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

double IterativeMean::arithmetic_mean(vector<double> data) {
	double sum_all=accumulate(data.begin(),data.end(),0.0);
	double mean=sum_all/data.size();
	return mean;
}

double IterativeMean::geometric_mean(vector<double> data) {
	double product_all=accumulate(data.begin(),data.end(),1.0,multiplies<double>());
	double mean=pow(product_all,(1.0/data.size()));
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
		AGM() {}
		AGM(vector<double>* input):a(2,0) { raw=*input; }
};

class Geothdian:public IterativeMean { //https://xkcd.com/2435/
	//Pretty similar to the AGM, honestly!
	private:
		vector<double> a;
	protected:
		virtual double step_average()
			{ sort(a.begin(),a.end()); double med_next=a[1]; double a1_next=arithmetic_mean(a); a[0]=geometric_mean(a); a[1]=a1_next; a[2]=med_next; return med_next; }
		virtual void reset() { sort(raw.begin(),raw.end()); a[0]=geometric_mean(raw); a[1]=arithmetic_mean(raw); if(raw.size()%2){ a[2]=raw[raw.size()/2]; }else{ a[2]=double(raw[raw.size()/2]+raw[(raw.size()/2)+1])/2.0; } } //Almost long enough to consider expanding.
	public:
		Geothdian() {}
		Geothdian(vector<double>* input):a(3,0) { raw=*input; } //Geothdian will sort the passed in vector by size.
};


class ClocksAtSea:public IterativeMean { //An experimental mean based on the idea of collective indication, e.g. time by clocks at sea.
	private:
		vector<double> set_times;
	protected:
		virtual double step_average();
		virtual void reset() { set_times.assign(raw.size(),0); }
	public:
		ClocksAtSea(vector<double>* input) { raw=*input; }
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
	auto set_max=max_element(set_times.begin(),set_times.end());
	auto set_min=min_element(set_times.begin(),set_times.end());
	if(abs(*set_max-naivmean)>abs(*set_min-naivmean)) {
		new_mean=(naivmean*set_times.size()-*set_max)/(set_times.size()-1);
		*set_max=new_mean;
	}
	else if(abs(*set_max-naivmean)<abs(*set_min-naivmean)) {
		new_mean=(naivmean*set_times.size()-*set_min)/(set_times.size()-1);
		*set_min=new_mean;
	}
	else { //If there are two furthest elements equidistant in opposite directions, reset them both.
		new_mean=naivmean;
		*set_max=new_mean;
		*set_min=new_mean;
	}
	return new_mean/iteration;
}

class ExtendedAGM { //An AGM-style mean - extended to accept negative values!
	protected:
		vector<double> raw; //Data to be averaged.
		double precision=0.000001; //Precision to which mean should be calculated.
		double mean;
		AGM negAGM;
		AGM posAGM;
		void initialize(); //Prep for computation
	public:
		ExtendedAGM(vector<double>* input) { raw=*input; }
		void set_precision(double precis)
			{ precision=precis; }
		double calculate();
		double calculate(unsigned int iteration_depth);
		unsigned int size()
			{ return raw.size(); }
		double operator[] (unsigned int i)
			{ return raw[i]; }
};

void ExtendedAGM::initialize() { //AGM iteration can't handle negative values, so split the Real line.
	vector<double> negs;
	vector<double> poss;
	for(unsigned int i=0;i<raw.size();i++) {
		if(raw[i]>0) {
			poss.push_back(raw[i]);
		}
		else if(raw[i]<0) {
			negs.push_back(-raw[i]);
		}
		else if(raw[i]==0) {
			poss.push_back(raw[i]);
			negs.push_back(raw[i]);
		}
	}
	AGM negAGM(&negs);
	AGM posAGM(&poss);
	this->negAGM=negAGM;
	this->posAGM=posAGM;
}

double ExtendedAGM::calculate() {
	initialize();
	negAGM.set_precision(precision);
	posAGM.set_precision(precision);
	double negmean=negAGM.calculate();
	double posmean=posAGM.calculate();
	this->mean=posmean-negmean;
	return this->mean;
}

double ExtendedAGM::calculate(unsigned int iteration_depth) {
	initialize();
	double negmean=negAGM.calculate(iteration_depth);
	double posmean=posAGM.calculate(iteration_depth);
	this->mean=posmean-negmean;
	return this->mean;
}
