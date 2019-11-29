#include <ellipsometry.h>

void single_film_model::set_ambient(double n,double k)
{
	n1 = complex<double>(n,k);
}
void single_film_model::set_film(double n,double k)
{
	n2 = complex<double>(n,k);
}
void single_film_model::set_substrate(double n,double k)
{
	n3 = complex<double>(n,k);
}
void single_film_model::set_thickness(double t)
{
	d = t;
}
void single_film_model::set_aoi(double a)
{
	a1 = a;

}
void single_film_model::set_aoi_degrees(double a)
{
	set_aoi(M_PI * a / 180.0);

}
void single_film_model::set_wavelen(double w)
{
	wavelen = w;
}

void single_film_model::calculate()
{
	//cout << "\nCalculation: d= " << d;

	a2 = asin( (n1.real() / n2.real()) * sin(a1) );
	a3 = asin( (n2.real() / n3.real()) * sin(a2) );

	r12s = fresnel_rs(n1,n2,a1,a2);
	r12p = fresnel_rp(n1,n2,a1,a2);
	r23s = fresnel_rs(n2,n3,a2,a3);
	r23p = fresnel_rp(n2,n3,a2,a3);

	beta = (2.0*M_PI*(d/wavelen)*n2*cos(a2));

	complex<double> prop_term = exp(-2.0 * complex<double>(0,1) * beta);

	Rs = (r12s + r23s * prop_term) / (1.0 + r12s * r23s * prop_term);
	Rp = (r12p + r23p * prop_term) / (1.0 + r12p * r23p * prop_term);

	complex<double> result = Rp / Rs;

	_psi = atan(abs(result));
	_delta = arg(result);
	if(_delta < 0) _delta += 2.0*M_PI;
}

double single_film_model::delta() {return _delta;}
double single_film_model::psi() {return _psi;}


ellipsometry_model::ellipsometry_model()
{
	layer_data top_ambient;
	layer_data bottom_ambient;

	layers.push_back(top_ambient);
	layers.push_back(bottom_ambient);

	lambda = LAMBDA_HeNe * 1e-9;
}

void ellipsometry_model::output(ostream& os)
{
	vector<layer_data>::iterator it;
        int l=0;
        os << "\nLayers:";
    	for(it=layers.begin();it<layers.end();it++)
	{
            os<<"\n"<<l;
            os<<"\t t="<<it->t;
            os<<"\t angle="<<it->angle * 360/M_2_PI;
            os<<"\t n="<<it->n;
            l++;
	}
}

// Sets the angles of incidence for each of the interfaces in the stack using snell's law
// Uses the angle of the initial (ambient) layer as the starting point

void ellipsometry_model::set_aoi()
{
	vector<layer_data>::iterator it;

	double input_angle;
	double input_n;

	bool first = true;

	for(it=layers.begin();it<layers.end();it++)
	{
		if(first)
		{
			first = false;
		}
		else
		{
			it->angle = asin( input_n * sin(input_angle) / it->n.real() );
		}

                
		input_n = it->n.real();
		input_angle = it->angle;
	}
}

// Sets the initial (ambient) angle of incidence and then calls set_aoi() to calculate the rest.

void ellipsometry_model::set_aoi(double a)
{
	layers[0].angle = a;
	this->set_aoi();
}

// Convenience function to allow input of angle in degrees

void ellipsometry_model::set_aoi_degrees(double a)
{
    this->set_aoi(a * M_PI / 180.0);
}

void ellipsometry_model::add_layer(double n, double k, double t)
{
	layer_data bottom_ambient = layers.back();
	layers.pop_back();
	layers.push_back(layer_data(n,k,t));
	layers.push_back(bottom_ambient);
}


void ellipsometry_model::calculate()
{
        set_aoi();

	// Iterators for stepping over layer structure

	vector<layer_data>::iterator layer_0;
	vector<layer_data>::iterator layer_1;

	// Propagation, interface and intermediate transfer matrices

	element_matrix M_p;
	element_matrix I_p;
	element_matrix M_s;
	element_matrix I_s;
	element_matrix P;

	// All matrices initialized as identity matrices

	M_p.identity();
	I_p.identity();
	M_s.identity();
	I_s.identity();
	P.identity();

	// Layer counter

	int layer_no = 0;

	// Begin layer_1 iterating over layers

	for(layer_1=layers.begin();layer_1<layers.end();layer_1++)
	{
		if(layer_no>=1)
		{
                      /*  cout << "\n--------------" << layer_no;*/
                        
			if(layer_no==1)
			{
				// Start layer_0 iterating when layer_1 has already moved on to first layer after ambient

				layer_0=layers.begin();
			}

			// Once layer_0 is pointing to a real layer (not ambient) calculate propagation matrix

			if(layer_no>1)
			{
				P.propagation(*layer_0,lambda);
                     /*           
                                cout << "\nP[" << layer_no-1 << "]";
                                cout << P;
                                */
			}

			// Calculate interface matrices

			I_s.interface_s(*layer_0,*layer_1);
			I_p.interface_p(*layer_0,*layer_1);
                        
                 /*       cout << "\nIs[" << layer_no-1 << "," << layer_no << "]";
                        cout << I_s;
                        
                        cout << "\nIp[" << layer_no-1 << "," << layer_no << "]";
                        cout << I_p;*/

			// Multiply matrices. For the case where layer_0 is at the ambient and layer_1 is at the first
			// real layer, P still equals identity.

			M_s *= (P * I_s);
			M_p *= (P * I_p);
                        
             /*           cout << "\n--------------";
                        
                        cout << "\nMs[" << layer_no-1 << "," << layer_no << "]";
                        cout << M_s;
                        
                        cout << "\nMp[" << layer_no-1 << "," << layer_no << "]";
                        cout << M_p;
                        
                        cout << "\n--------------";
                        */
                        layer_0++;
		}

                layer_no++;

	}

	// Calculate the reflection coefficients

	Rs = M_s.r();
	Rp = M_p.r();
        
    /*    cout << "\n\n==================\nRs=" << Rs;
        
        cout << "\nRp=" << Rp << "\n==================\n";*/

	// Calculate ellipsometric angles

	complex<double> result = Rp / Rs;
        
  //      cout << "\nresult=" << result << "\n\n==================\n";

	_psi = atan(abs(result));
	_delta = arg(result);
        if(_delta < 0) _delta += 2.0*M_PI;

}

double ellipsometry_model::delta() {return _delta;}
double ellipsometry_model::psi() {return _psi;}

double ellipsometry_model::delta_degrees() {return _delta * 180.0 / M_PI;}
double ellipsometry_model::psi_degrees() {return _psi * 180.0 / M_PI;}

void ellipsometry_model::set_lambda(double l) {lambda = l;}
double ellipsometry_model::get_lambda() {return lambda;}

layer_data& ellipsometry_model::operator[](int i) {return layers[i];}
layer_data& ellipsometry_model::get_layer(int i) {return layers[i];}

layer_data& ellipsometry_model::top_ambient() {return layers.front();}
layer_data& ellipsometry_model::bottom_ambient() {return layers.back();}

void ellipsometry_model::set_nk_top(double n, double k){
    layers.front().n = complex<double>(n,k);
}

void ellipsometry_model::set_nk_bottom(double n, double k){
    layers.back().n = complex<double>(n,k);
}

void ellipsometry_model::set_n_top(double _n) {
    layers.front().n = complex<double>(_n,imag(layers.front().n));
}

void ellipsometry_model::set_k_top(double _k) {
    layers.front().n = complex<double>(real(layers.front().n),_k);
}

void ellipsometry_model::set_n_bottom(double _n) {
    layers.back().n = complex<double>(_n,imag(layers.back().n));
}

void ellipsometry_model::set_k_bottom(double _k) {
    layers.back().n = complex<double>(real(layers.back().n),_k);
}

