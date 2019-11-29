
#include <fresnel.h>
#include <string>
// layer_data data structure constructors

layer_data::layer_data()
{
	t = angle = 0;
	n = complex<double>(1,0);
}

layer_data::layer_data(const layer_data& l):t(l.t),angle(l.angle),n(l.n) {}


layer_data::layer_data(complex<double> _n, double _t, double _a)
{
	t = _t;
	angle = _a;
	n = _n;
}

layer_data::layer_data(double _n, double _k, double _t, double _a)
{
	t = _t;
	angle = _a;
	n = complex<double>(_n,_k);
}


layer_data::layer_data(complex<double> _n, double _t)
{
	t = _t;
	angle = 0;
	n = _n;
}

layer_data::layer_data(double _n, double _k, double _t)
{
	t = _t;
	angle = 0;
	n = complex<double>(_n,_k);
}

// Element matrix constructors, destructor and assignment

element_matrix::element_matrix() {

    a = b = c = d = complex<double>(0,0);

}

element_matrix::element_matrix(const element_matrix& orig) {

    a = orig.a;     b = orig.b;
    c = orig.c;     d = orig.d;

}

element_matrix::~element_matrix() {}

element_matrix& element_matrix::operator =(const element_matrix& orig) {

    a = orig.a;     b = orig.b;
    c = orig.c;     d = orig.d;

    return *this;
}

// layer_data matrix calculation functions - using indices and aoi

void element_matrix::interface_p(const complex<double>& n1, const complex<double>& n2, const double& theta_i)
{
	complex<double> r = fresnel_rp(n1,n2,theta_i);
	complex<double> t_rec = complex<double>(1,0) / fresnel_tp(n1,n2,theta_i);

	a = d = t_rec;
	b = c = r * t_rec;
}

void element_matrix::interface_s(const complex<double>& n1, const complex<double>& n2, const double& theta_i)
{
	complex<double> r = fresnel_rs(n1,n2,theta_i);
	complex<double> t_rec = complex<double>(1,0) / fresnel_ts(n1,n2,theta_i);

	a = d = t_rec;
	b = c = r * t_rec;
}

void element_matrix::propagation(const complex<double>& n, const double& t, const double& lambda, const double& theta)
{
	complex<double> beta = 2.0 * M_PI * (t / lambda) * n * complex<double>(0,1) * cos(theta);

	a = exp(-beta);
	b = c = complex<double>(0,0);
	d = exp(beta);
}

// layer_data matrix calculation functions - using layer_data structures
// WARNING: explicit transmission angles are used meaning it is assumed the layer_datas have already had the appropreate angles calculated

void element_matrix::interface_p(const layer_data& L1, const layer_data& L2)
{
	complex<double> r = fresnel_rp(L1.n,L2.n,L1.angle);//,L2.angle);
	complex<double> t_rec = complex<double>(1,0) / fresnel_tp(L1.n,L2.n,L1.angle);//,L2.angle);

	a = d = t_rec;
	b = c = r * t_rec;
}

void element_matrix::interface_s(const layer_data& L1, const layer_data& L2)
{
	complex<double> r = fresnel_rs(L1.n,L2.n,L1.angle);//,L2.angle);
	complex<double> t_rec = complex<double>(1,0) / fresnel_ts(L1.n,L2.n,L1.angle);//,L2.angle);

	a = d = t_rec;
	b = c = r * t_rec;
}

void element_matrix::propagation(const layer_data& L1, const double& lambda)
{
	complex<double> beta = 2.0 * M_PI * (L1.t / lambda) * L1.n * complex<double>(0,1) * cos(L1.angle);

	a = exp(-beta);
	b = c = complex<double>(0,0);
	d = exp(beta);
}

// Identity matrix

void element_matrix::identity()
{
	a = d = complex<double>(1,0);
	b = c = complex<double>(0,0);
}


complex<double> element_matrix::r()
{
	return c / a;
}



// Matrix algebra

element_matrix& element_matrix::operator*=(const element_matrix& _m2)
{

    complex<double> _a,_b,_c,_d;

    _a = a * _m2.a + b * _m2.c;
    _b = a * _m2.b + b * _m2.d;
    _c = c * _m2.a + d * _m2.c;
    _d = c * _m2.b + d * _m2.d;

    a = _a;
    b = _b;
    c = _c;
    d = _d;

    return *this;

}


element_matrix operator*(const element_matrix& _m1, const element_matrix& _m2)
{

    element_matrix r;

    r.a = _m1.a * _m2.a + _m1.b * _m2.c;
    r.b = _m1.a * _m2.b + _m1.b * _m2.d;
    r.c = _m1.c * _m2.a + _m1.d * _m2.c;
    r.d = _m1.c * _m2.b + _m1.d * _m2.d;

    return r;

}


ostream& operator<<(ostream& s, const element_matrix& m)
{
    s << "\n" << m.a << "\t" << m.b;
    s << "\n" << m.c << "\t" << m.d;
    return s;
}

// Fresnel coefficients - transmission angle calculated using refractive indices and aoi
// Remember: for p, the n and cos have opposite subscript in each term. For s they have the same

complex<double> fresnel_rp(const complex<double>& n1, const complex<double>& n2, const double& theta_i)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = sqrt(1 - pow((n1.real() / n2.real())*sin(theta_i),2)); // From snell's law and trig identity: cos^2(a) + sin^2(a) = 1

	complex<double> term1 = n1*cos_theta_t;
	complex<double> term2 = n2*cos_theta_i;

	return (term2 - term1)/(term1 + term2);
}

complex<double> fresnel_rs(const complex<double>& n1, const complex<double>& n2, const double&  theta_i)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = sqrt(1 - pow((n1.real() / n2.real())*sin(theta_i),2)); // From snell's law and trig identity: cos^2(a) + sin^2(a) = 1

	complex<double> term1 = n1*cos_theta_i;
	complex<double> term2 = n2*cos_theta_t;

	return (term1 - term2)/(term1 + term2);
}

complex<double> fresnel_tp(const complex<double>& n1, const complex<double>& n2, const double&  theta_i)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = sqrt(1 - pow((n1.real() / n2.real())*sin(theta_i),2)); // From snell's law and trig identity: cos^2(a) + sin^2(a) = 1

	complex<double> term1 = n1*cos_theta_t;
	complex<double> term2 = n2*cos_theta_i;

	return (2.0 * n1 * cos_theta_i)/(term1 + term2);
}

complex<double> fresnel_ts(const complex<double>& n1, const complex<double>& n2, const double&  theta_i)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = sqrt(1 - pow((n1.real() / n2.real())*sin(theta_i),2)); // From snell's law and trig identity: cos^2(a) + sin^2(a) = 1

	complex<double> term1 = n1*cos_theta_i;
	complex<double> term2 = n2*cos_theta_t;

	return (2.0 * n1 * cos_theta_i)/(term1 + term2);
}

// Fresnel coefficients - transmission angle supplied as parameter (make sure it's right!)

complex<double> fresnel_rp(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = cos(theta_t);

	complex<double> term1 = n1*cos_theta_t;
	complex<double> term2 = n2*cos_theta_i;

	return (term1 - term2)/(term1 + term2);
}

complex<double> fresnel_rs(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = cos(theta_t);

	complex<double> term1 = n1*cos_theta_i;
	complex<double> term2 = n2*cos_theta_t;

	return (term1 - term2)/(term1 + term2);
}

complex<double> fresnel_tp(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = cos(theta_t);

	complex<double> term1 = n1*cos_theta_t;
	complex<double> term2 = n2*cos_theta_i;

	return (2.0 * n1 * cos_theta_i)/(term1 + term2);
}

complex<double> fresnel_ts(const complex<double>& n1, const complex<double>& n2, const double& theta_i, const double& theta_t)
{
	double cos_theta_i = cos(theta_i);
	double cos_theta_t = cos(theta_t);

	complex<double> term1 = n1*cos_theta_i;
	complex<double> term2 = n2*cos_theta_t;

	return (2.0 * n1 * cos_theta_i)/(term1 + term2);
}

