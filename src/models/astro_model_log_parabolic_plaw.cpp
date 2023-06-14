#include "../../include/models/astro_model_log_parabolic_plaw.h"

NuFit::astro_model_log_parabolic_plaw::astro_model_log_parabolic_plaw(double emin, double emax)
{
	// model name
	model_name = std::string("log_parabolic_powerlaw");

	// specify names and ordering (indices) of model parameters
	std::string par0("astro_norm");
	std::string par1("astro_index");
	std::string par2("log_slope");
	par_names.push_back(par0);
	par_names.push_back(par1);
	par_names.push_back(par2);

	// store the ordering in a map
	for (unsigned int i=0; i<par_names.size(); ++i)
		parameters.insert(std::pair<std::string, unsigned int>(par_names[i], i));
	
	npars = parameters.size();
    this->emin = emin;
    this->emax = emax;

	return;
} 

double NuFit::astro_model_log_parabolic_plaw::get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const {
	// assumes pars is ordered according to names provided above

    if (energy<this->emin || energy>this->emax)
    {
        return 0;
    }
    double Gamma = pars[1]+pars[2]*TMath::Log(energy/1.e5);    
	return pars[0] * 1.e-18 * TMath::Power(energy / 1.e5, -1.0 * Gamma);
}
