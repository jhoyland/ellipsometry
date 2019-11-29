#ifndef __FRESNEL_H
#define __FRESNEL_H

#include <string>
#include <complex>
// layer_data data structure constructors



using namespace std;

struct layer_data
{
	double t; // Thickness
	double angle; // Angle
	complex<double> n; // Complex refractive index

	layer_data();
	layer_data(const layer_data& l);
	layer_data(complex<double> _n, double _t, double _a);
	layer_data(double _n, double _k, double _t, double _a);
	layer_data(complex<double> _n, double _t);
	layer_data(double _n, double _k, double _t);

};

class element_matrix
{
	complex<double> a;
	complex<double> b;
	complex<double> c;
	complex<double> d;

public:
	// Element matrix constructors, destructor and assignment

	element_matrix();
	element_matrix(const element_matrix& orig);
	~element_matrix();
	element_matrix& operator =(const element_matrix& orig);

	// layer_data matrix calculation functions - using indices and aoi

	void interface_p(const complex<double>& n1, const complex<double>& n2, const double& theta_i);
	void interface_s(const complex<double>& n1, const complex<double>& n2, const double& theta_i);
	void propagation(const complex<double>& n, const double& t, const double& lambda, const double& theta);

	// layer_data matrix calculation functions - using layer_data structures
	// WARNING: explicit transmission angles are used meaning it is assumed the layer_datas have already had the appropreate angles calculated

	void interface_p(const layer_data& L1, const layer_data& L2);
	void interface_s(const layer_data& L1, const layer_data& L2);
	void propagation(const layer_data& L1, const double& lambda);
	void identity();

	complex<double> r();

	// Matrix algebra

	element_matrix& operator*=(const element_matrix& _m2);



	friend element_matrix operator*(const element_matrix& _m1, const element_matrix& _m2);
	friend ostream& operator<<(ostream& s, const element_matrix& m);

};


// Fresnel coefficients - transmission angle calculated using refractive indices and aoi
// Remember: for p, the n and cos have opposite subscript in each term. For s they have the same

complex<double> fresnel_rp(const complex<double>& n1, const complex<double>& n2, const double& theta_i);
complex<double> fresnel_rs(const complex<double>& n1, const complex<double>& n2, const double& theta_i);
complex<double> fresnel_tp(const complex<double>& n1, const complex<double>& n2, const double& theta_i);
complex<double> fresnel_ts(const complex<double>& n1, const complex<double>& n2, const double& theta_i);

// Fresnel coefficients - transmission angle supplied as parameter (make sure it's right!)

complex<double> fresnel_rp(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t);
complex<double> fresnel_rs(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t);
complex<double> fresnel_tp(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t);
complex<double> fresnel_ts(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t);

#endif
