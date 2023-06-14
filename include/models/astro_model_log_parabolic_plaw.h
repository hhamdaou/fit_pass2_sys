/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the single powerlaw astro physical flux model
 * 
 */

#ifndef _ASTRO_MODEL_LOG_PARABOLIC_PLAW_H
#define _ASTRO_MODEL_LOG_PARABOLIC_PLAW_H

#include "../helpers.h"
#include "../astro_model_base.h"

namespace NuFit
{

	class astro_model_log_parabolic_plaw: public astro_model_base {
		public:
			astro_model_log_parabolic_plaw(double emin=0, double emax=1.e9);
			~astro_model_log_parabolic_plaw() { };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;
            double emin;
            double emax;
	
	};

}

#endif

