#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/bootstrap/toymc.h"



int main(int argc, char **argv)
{

	std::string outdir=std::string("/data/user/zzhang1/fit_pass2_sys/output/");

	if (argc<6) {
		std::cout << "expect six arguments: jobid[1,2] and njobs[1,2] and out_subdir. quitting!" << std::endl;
		return 1;
	}

	// get arguments relevant for profile scan
	int jobid1 = std::stoi(argv[1]);
	int njobs1 = std::stoi(argv[2]);
	int jobid2 = std::stoi(argv[3]);
	int njobs2 = std::stoi(argv[4]);
    std::string out_subdir = argv[5];

    int random_seed = jobid1;
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
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//min.change_astro_model(new_astro_model);

	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */

	NuFit::helpers::par_options astro_norm("astro_norm", 1.8, 0.01, 0., 5.0);
	NuFit::helpers::par_options astro_index("astro_index", 2.6, 0.01, 0, 5.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.2, 0.01, 0.0, 5.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.85, 0.01, 0.0, 5.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 2.0, 0.01, 0.0, 20.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -1.0, 1.0);
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

	/** .. and send to minimizer */
	min.set_options(options);
    min.set_tolerance(10);

    // run global fit to provide seed to profile llh scan when it fails.
    min.fit(false); // ** if true -> get profile LLH errors after minimization from ROOT Minuit2
    
    //////////////////////////////////////////////////////////

    //int nsteps_index = 13;
    int nsteps_index = 49;
    double index_min = 2.2;
    double index_max = 2.8;
    //double index_min = 1;
    //double index_max = 4.0;
    double ds = (index_max-index_min) / (nsteps_index-1);
    if (nsteps_index % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_index / njobs1;
    double index_min_tj = index_min + jobid1 * nsteps_job1 * ds;
    double index_max_tj = index_min_tj + (nsteps_job1-1) * ds;

    //int nsteps_norm = 16;
    int nsteps_norm = 61;
    double norm_min = 1.0;
    double norm_max = 2.5;
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

    /*
    // index = sv; norm = conv_norm
    int nsteps_1 = 13;
    double p1_min = 120;
    double p1_max = 3000;
    double ds = (p1_max-p1_min) / (nsteps_1-1);
    if (nsteps_1 % njobs1 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job1 = nsteps_1 / njobs1;
    double p1_min_tj = p1_min + jobid1 * nsteps_job1 * ds;
    double p1_max_tj = p1_min_tj + (nsteps_job1-1) * ds;

    int nsteps_2 = 16;
    double p2_min = 0.5;
    double p2_max = 2.75;
    ds = (p2_max-p2_min)/(nsteps_2-1);
    if (nsteps_2 % njobs2 != 0) {
        std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
        return 1;
    }
    int nsteps_job2 = nsteps_2 / njobs2;
    double p2_min_tj = p2_min + jobid2 * nsteps_job2 * ds;
    double p2_max_tj = p2_min_tj + (nsteps_job2-1) * ds;

    //scan
    NuFit::helpers::scan_options s_1("selfveto", nsteps_job1, p1_min_tj, p1_max_tj);
    NuFit::helpers::scan_options s_2("conv_norm", nsteps_job2, p2_min_tj, p2_max_tj);

    std::map<std::string, NuFit::helpers::scan_options> scan_profile;
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_1.name, s_1));
    scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_2.name, s_2));
    */

    min.scan_profile_llh(outdir+std::string("/")+out_subdir+std::string("/outfile_profile_llh_2d_part_")+std::to_string(jobid1)+std::string("_")+std::to_string(jobid2)+std::string(".txt"), scan_profile);

    return 0;
}
