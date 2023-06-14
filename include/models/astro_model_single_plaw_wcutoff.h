/* author: Hans Niederhausen
 * date: 03/09/2017
 *
 * this class manages the single powerlaw with cutoff astro physical flux model
 * 
 */

#ifndef _ASTRO_MODEL_SINGLE_PLAW_WCUTOFF_H
#define _ASTRO_MODEL_SINGLE_PLAW_WCUTOFF_H

#include "../helpers.h"
#include "../astro_model_base.h"

namespace NuFit
{

	class astro_model_single_plaw_wcutoff: public astro_model_base {
		public:
			astro_model_single_plaw_wcutoff();
			~astro_model_single_plaw_wcutoff() { };
			double get_flux(const double *pars, double &energy, double &coszen, double &ra, NuFit::helpers::neutrino_type &ptype) const;
	
	};

}

#endif

