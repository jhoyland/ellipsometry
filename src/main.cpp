#include <ellipsometry.h>

//#define __MATRIX_MODEL

int main()
{
	#if defined __MATRIX_MODEL

	ellipsometry_model model;

	#else

	single_film_model model;

	#endif

	double layer_n = 1.4;
	double layer_k = 1e-6;

	double substrate_n = 1.5104;
	double substrate_k = 9.3e-9;   //BK7

	double t0 = 4.5e-6;
	double t1 = 5.5e-6;
	double dt = 1e-7;

	double steps = (t1-t0)/dt;

	#if defined __MATRIX_MODEL

	model.set_nk_top(1,0);
	model.set_nk_bottom(substrate_n,substrate_k);
	model.add_layer(layer_n,layer_k,t0);
	model.set_aoi_degrees(60);

	#else

	model.set_ambient(1,0);
	model.set_film(layer_n,layer_k);
	model.set_substrate(substrate_n,substrate_k);
	model.set_thickness(t0);
	model.set_aoi_degrees(60);
	model.set_wavelen(LAMBDA_HeNe * 1e-9);

	#endif

	double t[(int)steps];
	double p[(int)steps];
	double d[(int)steps];

	cout << "\nn\tthick\tpsi\tdelta"; 

	for(int n = 0; n<(int)steps; n++)
	{
		t[n] = n * dt + t0;
		#if defined __MATRIX_MODEL
		model[1].t = t[n];
		#else
		model.set_thickness(t[n]);
		#endif
		model.calculate();
		p[n] = model.psi();
		d[n] = model.delta();

		cout << "\n" << n << "\t" << t[n] << "\t" << p[n] << "\t" << d[n]; 
	}
}




