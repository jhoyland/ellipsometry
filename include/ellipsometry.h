#ifndef __ELLIPSOMETRY_H
#define __ELLIPSOMETRY_H

#include <fresnel.h>
#include <vector>
#include <iostream>

using namespace std;

#define LAMBDA_HeNe 633

class single_film_model
{
	complex<double> n1;
	complex<double> n2;
	complex<double> n3;

	double d;

	double a1;
	double a2;
	double a3;

	complex<double> Rp;
	complex<double> Rs;

	complex<double> r12s;
	complex<double> r23s;
	complex<double> r12p;
	complex<double> r23p;

	complex<double> beta;

	double wavelen;

	double _delta;
	double _psi;

public:

	void set_ambient(double,double);
	void set_film(double,double);
	void set_substrate(double,double);
	void set_thickness(double);
	void set_aoi(double);
	void set_aoi_degrees(double);
	void set_wavelen(double);

	void calculate();

	double delta();
	double psi();


};

class ellipsometry_model
{
	vector<layer_data> layers;

	complex<double> Rp;
	complex<double> Rs;

	double _delta;
	double _psi;

	double lambda;

public:

	ellipsometry_model();

	void add_layer(double,double,double);
	void set_lambda(double);
	double get_lambda();

    void set_nk_top(double,double);
    void set_nk_bottom(double,double);

    void set_n_top(double);
    void set_k_top(double);

    void set_n_bottom(double);
    void set_k_bottom(double);

    layer_data& top_ambient();
    layer_data& bottom_ambient();

	void set_aoi(double);
	void set_aoi_degrees(double);
	void set_aoi();
	void calculate();

	double psi();
	double delta();
        
    double psi_degrees();
	double delta_degrees();
        
    void output(ostream&);

    layer_data& operator[](int);
    layer_data& get_layer(int i);
};

#endif

