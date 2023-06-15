#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/models/astro_model_segmented.h"
#include "include/models/astro_model_broken_plaw.h"
#include "include/bootstrap/toymc.h"



int main(int argc, char **argv)
{

	std::string outdir=std::string("/data/user/hhamdaoui/fit_pass2_sys_output/");

	if (argc<5) {
		std::cout << "expect six arguments: jobid[1,2] and njobs[1,2] and astro_norm_seed, astro_index_seed. quitting!" << std::endl;
		return 1;
	}

	// get arguments relevant for profile scan
	int jobid1 = std::stoi(argv[1]);
	int njobs1 = std::stoi(argv[2]);
	int jobid2 = std::stoi(argv[3]);
	int njobs2 = std::stoi(argv[4]);
    float astro_norm_seed;
    float astro_index_seed;
    if (argc<6) {
        std::cout << "using default astro flux seed. astro_norm_seed = 1.58, astro_index_seed = 2.53." << std::endl;
        astro_norm_seed = 1.58;
        astro_index_seed = 2.53;
    }else{
        astro_norm_seed = std::stoi(argv[5]);
        astro_index_seed = std::stoi(argv[6]);
    }
    //int random_seed = jobid1;
//        int njobs = 1;

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
	//NuFit::astro_model_single_plaw *new_astro_model = new NuFit::astro_model_single_plaw(1.e3,1.e9);
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//NuFit::astro_model_log_parabolic_plaw *new_astro_model = new NuFit::astro_model_log_parabolic_plaw();
	//NuFit::astro_model_broken_plaw *new_astro_model = new NuFit::astro_model_broken_plaw();
	//NuFit::astro_model_segmented *new_astro_model = new NuFit::astro_model_segmented();
	//min.change_astro_model(new_astro_model);
	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */

	NuFit::helpers::par_options astro_norm("astro_norm", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_index("astro_index", 2.6, 0.01, 0, 10.0);
	//NuFit::helpers::par_options log_slope("log_slope", 0, 0.001, -0.5, 0.5);
	//NuFit::helpers::par_options astro_index1("astro_index1", 2.6, 0.01, 0, 10.0);
	//NuFit::helpers::par_options astro_index2("astro_index2", 2.6, 0.01, 0, 10.0);
	//NuFit::helpers::par_options astro_logEbreak("astro_logEbreak", 4.3, 0.01, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.2, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.85, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 2.0, 0.01, 0.0, 20.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -1, 1);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.0, 2.0);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.0, 2.0);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.0, 2.0);
	NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.01, -3.0, 2.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.01, -0.5, 0.5);
	NuFit::helpers::par_options selfveto("selfveto", 300, 50, 50, 3000);
	NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0.5, 0.01, 0, 1);

	/** .. package everything */
	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(log_slope.name, log_slope));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index1.name, astro_index1));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index2.name, astro_index2));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_logEbreak.name, astro_logEbreak));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm.name, conv_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm.name, muon_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm.name, prompt_norm));	
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index.name, cr_index));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff.name, dom_eff));
    options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering.name, scattering));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption.name, absorption));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep0.name, holeicep0));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep1.name, holeicep1));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(selfveto.name, selfveto));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(hadronicinteraction.name, hadronicinteraction));

	/** write seed histograms */
	std::map<std::string, double> pars;
	pars.insert(std::pair<std::string, double>(muon_norm.name, muon_norm.seed));
	pars.insert(std::pair<std::string, double>(conv_norm.name, conv_norm.seed));
	//pars.insert(std::pair<std::string, double>(log_slope.name, log_slope.seed));
	pars.insert(std::pair<std::string, double>(prompt_norm.name, prompt_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_norm.name, astro_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_index.name, astro_index.seed));
	//pars.insert(std::pair<std::string, double>(astro_index1.name, astro_index1.seed));
	//pars.insert(std::pair<std::string, double>(astro_index2.name, astro_index2.seed));
	//pars.insert(std::pair<std::string, double>(astro_logEbreak.name, astro_logEbreak.seed));
	pars.insert(std::pair<std::string, double>(cr_index.name, cr_index.seed));
    pars.insert(std::pair<std::string, double>(dom_eff.name, dom_eff.seed));
    pars.insert(std::pair<std::string, double>(scattering.name, scattering.seed));
	pars.insert(std::pair<std::string, double>(absorption.name, absorption.seed));
	pars.insert(std::pair<std::string, double>(holeicep0.name, holeicep0.seed));
	pars.insert(std::pair<std::string, double>(holeicep1.name, holeicep1.seed));
	pars.insert(std::pair<std::string, double>(selfveto.name, selfveto.seed));
	pars.insert(std::pair<std::string, double>(hadronicinteraction.name, hadronicinteraction.seed));


	/** .. and send to minimizer */
	min.set_options(options);
    //min.set_tolerance(1000);
    min.set_tolerance(1);

    /*
    int nsteps_norm = 11;
    double norm_min = 0.0;
    double norm_max = 1.0;
    double ds = (norm_max-norm_min)/(nsteps_norm-1);
    if (nsteps_norm % njobs2 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job2 = nsteps_norm / njobs2;
    double norm_min_tj = norm_min + jobid2 * nsteps_job2 * ds;
    double norm_max_tj = norm_min_tj + (nsteps_job2-1) * ds;

    //scan
    NuFit::helpers::scan_options s_astro_norm("hadronicinteraction", nsteps_job2, norm_min_tj, norm_max_tj);
    */

    int nsteps_norm = 10;
    double emin = 1.e2;
    double emax = 1.e8;
    if (jobid1<20)
    {
        //emin = std::pow(10,3.5+0.05*jobid1); // broken power law and log parabola
        emin = std::pow(10,3.0+0.05*jobid1); //Single power law
        emax = std::pow(10,8);
    }
    else
    {
        emin = std::pow(10,2);
        emax = std::pow(10,6+0.05*(jobid1-20));
    }
    std::cout<<"emin: "<<emin<<" "<<"emax: "<<emax<<std::endl;
	NuFit::astro_model_single_plaw *new_astro_model1 = new NuFit::astro_model_single_plaw(emin,emax);
	//NuFit::astro_model_log_parabolic_plaw *new_astro_model1 = new NuFit::astro_model_log_parabolic_plaw(emin,emax);
	//NuFit::astro_model_broken_plaw *new_astro_model1 = new NuFit::astro_model_broken_plaw(emin,emax);
	//NuFit::astro_model_segmented *new_astro_model = new NuFit::astro_model_segmented();
	min.change_astro_model(new_astro_model1);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;

    min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".txt"), scan_profile, outdir+std::string("injected_hist")+std::string(".root"));
    return 0;
}
