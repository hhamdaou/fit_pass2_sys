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

	std::string outdir=std::string("/data/user/zzhang1/fit_pass2_sys/output/");

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
	NuFit::astro_model_single_plaw *new_astro_model = new NuFit::astro_model_single_plaw(1.e5,1.e6);
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//NuFit::astro_model_log_parabolic_plaw *new_astro_model = new NuFit::astro_model_log_parabolic_plaw();
	//NuFit::astro_model_broken_plaw *new_astro_model = new NuFit::astro_model_broken_plaw();
	//NuFit::astro_model_segmented *new_astro_model = new NuFit::astro_model_segmented();
	min.change_astro_model(new_astro_model);
	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */
    /*
	//NuFit::helpers::par_options astro_norm("astro_norm", 1.58, 0.001, 0., 10.);
	NuFit::helpers::par_options astro_norm("astro_norm", astro_norm_seed, 0.003, 0., 10.);
	//NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.01, 0, 10.0);
	NuFit::helpers::par_options astro_index("astro_index", astro_index_seed, 0.003, 0, 10.0);
    */

    /*
	NuFit::helpers::par_options astro_norm("astro_norm", 1.80981, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_index("astro_index", 2.58259, 0.01, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.14708, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.57392, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 1.06014, 0.01, 0.0, 20.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0561507, 0.01, -1, 1);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 0.912188, 0.01, 0.0, 2.0);
    NuFit::helpers::par_options scattering("scattering", 0.987459, 0.01, 0.0, 2.0);
	NuFit::helpers::par_options absorption("absorption", 0.971261, 0.01, 0.0, 2.0);
	NuFit::helpers::par_options holeicep0("holeicep0", -1.05944, 0.01, -3.0, 2.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.00911564, 0.01, -0.5, 0.5);
	NuFit::helpers::par_options selfveto("selfveto", 363.308, 50, 50, 3000);
	//NuFit::helpers::par_options selfveto("selfveto", 3000, 50, 50, 3000);
	NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 1.0, 0.01, 0, 1);
    */

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
	NuFit::helpers::par_options cross_section("cross_section", 300, 50, 50, 3000);

	NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0.5, 0.01, 0, 1);
	//NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 1.0, 0.0, 0, 1);

    /*
	NuFit::helpers::par_options astro_norm("astro_norm", 0.1, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_index("astro_index", 1.5, 0.01, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 0.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.00, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 1.0, 0.01, 0.0, 20.0);
    */


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
    /*
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm.name, conv_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm.name, muon_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm.name, prompt_norm));	
    */

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

    /*
	NuFit::helpers::par_options astro_norm_0("astro_norm_0", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_1("astro_norm_1", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_2("astro_norm_2", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_3("astro_norm_3", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_4("astro_norm_4", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_5("astro_norm_5", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_6("astro_norm_6", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_7("astro_norm_7", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_8("astro_norm_8", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_9("astro_norm_9", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_10("astro_norm_10", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_11("astro_norm_11", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_12("astro_norm_12", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_13("astro_norm_13", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_14("astro_norm_14", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_15("astro_norm_15", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_16("astro_norm_16", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_17("astro_norm_17", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_norm_18("astro_norm_18", 1.8, 0.01, 0., 10.0);
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

	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_0.name, astro_norm_0));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_1.name, astro_norm_1));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_2.name, astro_norm_2));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_3.name, astro_norm_3));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_4.name, astro_norm_4));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_5.name, astro_norm_5));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_6.name, astro_norm_6));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_7.name, astro_norm_7));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_8.name, astro_norm_8));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_9.name, astro_norm_9));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_10.name, astro_norm_10));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_11.name, astro_norm_11));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_12.name, astro_norm_12));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_13.name, astro_norm_13));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_14.name, astro_norm_14));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_15.name, astro_norm_15));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_16.name, astro_norm_16));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_17.name, astro_norm_17));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_18.name, astro_norm_18));
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

	std::map<std::string, double> pars;
	pars.insert(std::pair<std::string, double>(astro_norm_0.name, astro_norm_0.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_1.name, astro_norm_1.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_2.name, astro_norm_2.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_3.name, astro_norm_3.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_4.name, astro_norm_4.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_5.name, astro_norm_5.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_6.name, astro_norm_6.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_7.name, astro_norm_7.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_8.name, astro_norm_8.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_9.name, astro_norm_9.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_10.name, astro_norm_10.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_11.name, astro_norm_11.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_12.name, astro_norm_12.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_13.name, astro_norm_13.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_14.name, astro_norm_14.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_15.name, astro_norm_15.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_16.name, astro_norm_16.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_17.name, astro_norm_17.seed));
	pars.insert(std::pair<std::string, double>(astro_norm_18.name, astro_norm_18.seed));
	pars.insert(std::pair<std::string, double>(muon_norm.name, muon_norm.seed));
	pars.insert(std::pair<std::string, double>(conv_norm.name, conv_norm.seed));
	pars.insert(std::pair<std::string, double>(prompt_norm.name, prompt_norm.seed));
	pars.insert(std::pair<std::string, double>(cr_index.name, cr_index.seed));
    pars.insert(std::pair<std::string, double>(dom_eff.name, dom_eff.seed));
    pars.insert(std::pair<std::string, double>(scattering.name, scattering.seed));
	pars.insert(std::pair<std::string, double>(absorption.name, absorption.seed));
	pars.insert(std::pair<std::string, double>(holeicep0.name, holeicep0.seed));
	pars.insert(std::pair<std::string, double>(holeicep1.name, holeicep1.seed));
	pars.insert(std::pair<std::string, double>(selfveto.name, selfveto.seed));
	pars.insert(std::pair<std::string, double>(hadronicinteraction.name, hadronicinteraction.seed));
    */

	/** .. and send to minimizer */
	min.set_options(options);
    //min.set_tolerance(1000);
    min.set_tolerance(1);

    min.fit(false); // ** if true -> get profile LLH errors after minimization from ROOT Minuit2
    //std::string outfile("./output/hists_seed.root");
	//wrapper.get_histograms(outfile,pars);
    //min.get_bestpars(pars);
	//outfile = std::string("./output/hists_fit.root");
	//wrapper.get_histograms(outfile, pars);

    /*
    int nsteps_selfveto = 29;
    // hadronic interaction
    //double selfveto_min = -1;
    //double selfveto_max = 2;

    //selfveto
    //double selfveto_min = 100;
    //double selfveto_max = 3000;
    
    //conv_norm
    double selfveto_min = 0.0;
    double selfveto_max = 3.0; 
    double ds = (selfveto_max-selfveto_min) / (nsteps_selfveto-1);
    if (nsteps_selfveto % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_selfveto / njobs1;
    double selfveto_min_tj = selfveto_min + jobid1 * nsteps_job1 * ds;
    double selfveto_max_tj = selfveto_min_tj + (nsteps_job1-1) * ds;

    //scan
    NuFit::helpers::scan_options s_selfveto("conv_norm", nsteps_job1, selfveto_min_tj, selfveto_max_tj);
    //NuFit::helpers::scan_options s_selfveto("hadronicinteraction", nsteps_job1, selfveto_min_tj, selfveto_max_tj);
    //NuFit::helpers::scan_options s_selfveto("selfveto", nsteps_job1, selfveto_min_tj, selfveto_max_tj);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_selfveto.name, s_selfveto));

	min.scan_profile_llh(outdir+std::string("outfile_profile_llh_part_")+std::to_string(jobid1)+std::string(".txt"), scan_profile);
    
    //min.get_bestpars(pars);
    */

    //////////////////////////////////////////////////////////
    /*
    int nsteps_index = 13;
    //double index_min = 2.2;
    //double index_max = 2.8;
    double index_min = 1.2;
    double index_max = 1.8;
    double ds = (index_max-index_min) / (nsteps_index-1);
    if (nsteps_index % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_index / njobs1;
    double index_min_tj = index_min + jobid1 * nsteps_job1 * ds;
    double index_max_tj = index_min_tj + (nsteps_job1-1) * ds;

    int nsteps_norm = 16;
    //double norm_min = 1.0;
    //double norm_max = 2.5;
    double norm_min = 0.0;
    double norm_max = 0.3;
    ds = (norm_max-norm_min)/(nsteps_norm-1);
    if (nsteps_norm % njobs2 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job2 = nsteps_norm / njobs2;
    double norm_min_tj = norm_min + jobid2 * nsteps_job2 * ds;
    double norm_max_tj = norm_min_tj + (nsteps_job2-1) * ds;

    //scan
    NuFit::helpers::scan_options s_astro_norm("astro_norm", nsteps_job2, norm_min_tj, norm_max_tj);
    NuFit::helpers::scan_options s_astro_index("astro_index", nsteps_job1, index_min_tj, index_max_tj);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index.name, s_astro_index));
    */

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

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
    */

    //min.scan_profile_llh_asimov(outdir+std::string("outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".txt"), scan_profile, outdir+std::string("injected_hist")+std::string(".root"));


    //min.get_bestpars(pars);

    /** write best fit histograms */
	//std::string outfile = outdir+std::string("hists_fit.root");
    //wrapper.get_histograms(outfile, pars);
	
	/** example of how to do a profile llh scan */
    /*
	int nsteps_index = 5;	
	double index_min = 2.3;
	double index_max = 2.6;

	double ds = (index_max-index_min) / (nsteps_index-1);

	if (nsteps_index % njobs != 0) {
		std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
		return 1;
	}

	int nsteps_job = nsteps_index / njobs;
        */
	//double index_min_tj = index_min + jobid * nsteps_job * ds;
	//double index_max_tj = index_min_tj + (nsteps_job-1) * ds;
	
        //scan
        /*
	NuFit::helpers::scan_options s_astro_norm("astro_norm", 20, 1.5, 1.7);
	//NuFit::helpers::scan_options s_astro_index("astro_index", nsteps_job, index_min_tj, index_max_tj);
	NuFit::helpers::scan_options s_astro_index("astro_index", 20, 2.3, 2.6);

	//std::cout << ds << " " << index_min_tj << " " << index_max_tj << " " << nsteps_job << std::endl;
	
	std::map<std::string, NuFit::helpers::scan_options> scan;
	scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
	scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index.name, s_astro_index));
	

	min.set_flush_rate(1); // set frequency of file dumps of scan results during scan

	//min.scan_profile_llh(outdir+std::string("outfile_profile_llh_part")+std::to_string(jobid)+std::string(".txt"), scan);
	//min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part")+std::string(".txt"), scan);
        
                        
        std::map<std::string, double> seed;
        seed.insert(std::pair<std::string, double>(muon_norm.name,1.66672));
        seed.insert(std::pair<std::string, double>(muon_norm1.name,1));
        seed.insert(std::pair<std::string, double>(conv_norm.name,1.03438));
        seed.insert(std::pair<std::string, double>(prompt_norm.name,5.02903e-5));
	min.scan_llh(outdir+std::string("outfile_llh_2d_part")+std::string(".txt"),scan,seed);
        */
	
	//scan_profile
	/*
	NuFit::helpers::scan_options s_astro_norm_profile("astro_norm", 10, 1.5, 1.7);
	NuFit::helpers::scan_options s_astro_index_profile("astro_index", 10, 2.3, 2.6);

	std::map<std::string, NuFit::helpers::scan_options> scan_profile;
	scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm_profile.name, s_astro_norm_profile));
	scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index_profile.name, s_astro_index_profile));
	
	min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part")+std::string(".txt"), scan_profile);
        */
    

    //toymc
    //std::string outfile = outdir+std::string("hists_fit.root");
    //wrapper.get_histograms(outfile, pars);

    /*
    pars[astro_norm.name]=1.66;
    pars[astro_index.name]=2.53;
    pars[muon_norm.name]=1.15;
    pars[conv_norm.name]=2.01;
    pars[prompt_norm.name]=1.0;
    pars[cr_index.name]=0.015;
    pars[absorption.name]=0.975;
    pars[dom_eff.name]=0.936;
    pars[hadronicinteraction.name]=1.66;
    pars[holeicep0.name]=-0.68;
    pars[holeicep1.name]=0.007;
    pars[scattering.name]=1.05;
    pars[selfveto.name]=296;
    */
    
    //pars[astro_norm.name]=1.66;
    //pars[astro_index.name]=2.53;
    //pars[muon_norm.name]=1.14;
    //pars[conv_norm.name]=1.50;
    //pars[prompt_norm.name]=1.0;
    //pars[cr_index.name]=0.036;
    //pars[absorption.name]=0.945;
    //pars[dom_eff.name]=0.915;
    //pars[hadronicinteraction.name]=1.0;
    //pars[holeicep0.name]=-0.86;
    //pars[holeicep1.name]=0.015;
    //pars[scattering.name]=1.03;
    //pars[selfveto.name]=348;

    /*
    pars[astro_norm.name]=1.66;
    pars[astro_index.name]=2.53;
    pars[muon_norm.name]=1.15;
    pars[conv_norm.name]=1.57;
    pars[prompt_norm.name]=1.0;
    pars[cr_index.name]=0.056;
    pars[absorption.name]=0.97;
    pars[dom_eff.name]=0.91;
    pars[hadronicinteraction.name]=1.0;
    pars[holeicep0.name]=-1.06;
    pars[holeicep1.name]=0.008;
    pars[scattering.name]=0.99;
    pars[selfveto.name]=355;
    */

    pars[astro_norm.name]=2.25518;
    pars[astro_index.name]=2.6945;
    //pars[log_slope.name]=0.156656;
    pars[muon_norm.name]=1.13886;
    pars[conv_norm.name]=1.58503;
    pars[prompt_norm.name]=0.170396;
    pars[cr_index.name]=0.0457014;
    pars[absorption.name]=0.967379;
    pars[dom_eff.name]=0.907906;
    pars[hadronicinteraction.name]=0.871643;
    pars[holeicep0.name]=-0.569008;
    pars[holeicep1.name]=0.0146574;
    pars[scattering.name]=1.02159;
    pars[selfveto.name]=917.949;

    double random_seed = 1;
    NuFit::toymc sim(wrapper, random_seed);

    //std::string outfile_toy = outdir+std::string("hists_toymc_job")+std::to_string(random_seed)+std::string(".root");
    //std::string outfile_aux = outdir+std::string("outfile_aux_job")+std::to_string(random_seed)+std::string(".txt");
    std::string outfile_toy = outdir+std::string("hists_toymc_job")+std::to_string(0)+std::string(".root");
    std::string outfile_aux = outdir+std::string("outfile_aux_job")+std::to_string(0)+std::string(".txt");

    unsigned int nsamples = 500;
    //unsigned int nsamples = 2;
    std::cout<< pars.size();
    //sim.run_simulation_waux(outfile_toy, outfile_aux, pars, nsamples);
    //sim.run_simulation(outfile_toy,pars,nsamples);

    std::map<std::string, double> point2;
    //point2.insert(std::pair<std::string, double>("astro_norm", 0.9));
    //point2.insert(std::pair<std::string, double>("astro_index", 2.5));
    min.set_flush_rate(1);
    
    
    int toy_id_low = int((nsamples/njobs1)*jobid1);
    int toy_id_high = int((nsamples/njobs1)*(jobid1+1));
    std::string outfile_toyfits = outdir+std::string("outfile_toy_")+std::to_string(toy_id_low)+std::string(".txt");
    //pars[prompt_norm.name]=0.0;
    //pars[astro_norm.name]=2.0;
    //pars[astro_index.name]=3.0;
    //std::cout<<pars.size();
    
    
    
    //min.run_toyfits(outfile_toy, outfile_toyfits, point2, pars, 0, nsamples, true, outfile_aux);
    min.run_toyfits(outfile_toy, outfile_toyfits, point2, pars, toy_id_low, toy_id_high, false);
    
    return 0;
}
