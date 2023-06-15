#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/models/astro_model_broken_plaw.h"
#include "include/bootstrap/toymc.h"



int main(int argc, char **argv)
{

	std::string outdir=std::string("/data/user/hhamdaoui/fit_pass2_sys_output/");

	/** create analysis. */
	NuFit::analysis wrapper;

	/** need to specifically create the analysis. */
	/** this function contains all the analysis specific details and can be edited in ./src/analysis.cpp */
	wrapper.create(); 

	/** get occasional prints of minimization progress from inside likelihood function */
	//wrapper.set_verbosity(true);

	/** use stats class to analyze NuFit::analysis */
	NuFit::stats min(wrapper);
	/** option to change model (NuFit::analysis uses single powerlaw as default) */
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//NuFit::astro_model_log_parabolic_plaw *new_astro_model = new NuFit::astro_model_log_parabolic_plaw();
	NuFit::astro_model_broken_plaw *new_astro_model = new NuFit::astro_model_broken_plaw();
	min.change_astro_model(new_astro_model);
	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */

    //NuFit::helpers::par_options log_slope("log_slope", 0, 0.001, -0.5, 0.5);
    NuFit::helpers::par_options astro_index1("astro_index1", 2.6, 0.01, 0, 10.0);
    NuFit::helpers::par_options astro_index2("astro_index2", 2.6, 0.01, 0, 10.0);
    NuFit::helpers::par_options astro_logEbreak("astro_logEbreak", 4.3, 0.01, 0, 10.0);

	NuFit::helpers::par_options astro_norm("astro_norm", 1.58, 0.01, 0., 10.);
	//NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.01, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.5, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.1, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.1, 0.0, 100.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -0.5, 0.5);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.0, 10.0);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.0, 10.0);
    NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.01, -1.0, 2.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.01, -0.2, 0.2);
    NuFit::helpers::par_options selfveto("selfveto", 750, 50, 5.0, 3000);
    NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0, 0.01, 0, 1);

    int nsteps = 29;
    double s_min = 1;
    double s_max = 3;

    //scan
    NuFit::helpers::scan_options s("astro_norm", nsteps, s_min, s_max);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s.name, s));

    std::map<std::string, double> asimov_point;

    /*
    asimov_point.insert(std::pair<std::string, double>(muon_norm.name, 1.15));
    asimov_point.insert(std::pair<std::string, double>(conv_norm.name, 1.57));
    asimov_point.insert(std::pair<std::string, double>(prompt_norm.name, 1.06));
    asimov_point.insert(std::pair<std::string, double>(astro_norm.name, 1.81));
    asimov_point.insert(std::pair<std::string, double>(astro_index.name, 2.58));
    asimov_point.insert(std::pair<std::string, double>(dom_eff.name, 0.91));
    asimov_point.insert(std::pair<std::string, double>(scattering.name, 0.99));
    asimov_point.insert(std::pair<std::string, double>(absorption.name, 0.97));
    asimov_point.insert(std::pair<std::string, double>(holeicep0.name, -1.06));
    asimov_point.insert(std::pair<std::string, double>(holeicep1.name, 0.01));
    asimov_point.insert(std::pair<std::string, double>(cr_index.name, 0.056));
    asimov_point.insert(std::pair<std::string, double>(selfveto.name, 363));
    asimov_point.insert(std::pair<std::string, double>(hadronicinteraction.name, 1.0));
    */

    /*
    asimov_point.insert(std::pair<std::string, double>(muon_norm.name, 1.14));
    asimov_point.insert(std::pair<std::string, double>(conv_norm.name, 1.59));
    asimov_point.insert(std::pair<std::string, double>(prompt_norm.name, 0.17));
    asimov_point.insert(std::pair<std::string, double>(astro_norm.name, 2.26));
    asimov_point.insert(std::pair<std::string, double>(astro_index.name, 2.69));
    asimov_point.insert(std::pair<std::string, double>(log_slope.name, 0.16));
    asimov_point.insert(std::pair<std::string, double>(dom_eff.name, 0.91));
    asimov_point.insert(std::pair<std::string, double>(scattering.name, 1.02));
    asimov_point.insert(std::pair<std::string, double>(absorption.name, 0.97));
    asimov_point.insert(std::pair<std::string, double>(holeicep0.name, -0.57));
    asimov_point.insert(std::pair<std::string, double>(holeicep1.name, 0.01));
    asimov_point.insert(std::pair<std::string, double>(cr_index.name, 0.046));
    asimov_point.insert(std::pair<std::string, double>(selfveto.name, 918));
    asimov_point.insert(std::pair<std::string, double>(hadronicinteraction.name, 0.87));
    */

    asimov_point.insert(std::pair<std::string, double>(muon_norm.name, 1.12));
    asimov_point.insert(std::pair<std::string, double>(conv_norm.name, 1.59));
    asimov_point.insert(std::pair<std::string, double>(prompt_norm.name, 1.42));
    asimov_point.insert(std::pair<std::string, double>(astro_norm.name, 1.72));
    asimov_point.insert(std::pair<std::string, double>(astro_index1.name, 0.70));
    asimov_point.insert(std::pair<std::string, double>(astro_index2.name, 2.83));
    asimov_point.insert(std::pair<std::string, double>(astro_logEbreak.name, 4.41));
    asimov_point.insert(std::pair<std::string, double>(dom_eff.name, 0.90));
    asimov_point.insert(std::pair<std::string, double>(scattering.name, 1.02));
    asimov_point.insert(std::pair<std::string, double>(absorption.name, 0.96));
    asimov_point.insert(std::pair<std::string, double>(holeicep0.name, -0.31));
    asimov_point.insert(std::pair<std::string, double>(holeicep1.name, -0.01));
    asimov_point.insert(std::pair<std::string, double>(cr_index.name, 0.045));
    asimov_point.insert(std::pair<std::string, double>(selfveto.name, 1439));
    asimov_point.insert(std::pair<std::string, double>(hadronicinteraction.name, 0.89));
    min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d")+std::string(".txt"), scan_profile, outdir+std::string("injected_hist")+std::string(".root"), asimov_point);
    return 0;
}
