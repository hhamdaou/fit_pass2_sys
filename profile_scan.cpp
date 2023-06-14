#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/models/astro_model_log_parabolic_plaw.h"
#include "include/models/astro_model_segmented.h"
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
    //NuFit::astro_model_log_parabolic_plaw *new_astro_model = new NuFit::astro_model_log_parabolic_plaw();

    NuFit::astro_model_segmented *new_astro_model = new NuFit::astro_model_segmented();
    min.change_astro_model(new_astro_model);

	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */

    /*
	NuFit::helpers::par_options astro_norm("astro_norm", 1.8, 0.01, 0., 10.0);
	NuFit::helpers::par_options astro_index("astro_index", 2.6, 0.01, 0, 10.0);
    //NuFit::helpers::par_options log_slope("log_slope", 0, 0.001, -0.5, 0.5);
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
    */

    // background check step1
    /*
	NuFit::helpers::par_options astro_norm("astro_norm", 1.66, 0.0, 0., 10.);
	NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.0, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.2, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.0, 0.0, 100.0);
    NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.0, -1, 1);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.0, 0.7, 1.3);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.0, 0.79, 1.21);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.0, 0.79, 1.21);
	NuFit::helpers::par_options holeicep0("holeicep0", 0.0, 0.0, -10.0, 10.0);
	NuFit::helpers::par_options holeicep1("holeicep1", 0.0, 0.0, -2, 2);
	NuFit::helpers::par_options hadronicinteraction("hadronicinteraction", 0.0, 0.0, -1, 2);
    */

    // background check step2
    /*
	NuFit::helpers::par_options astro_norm("astro_norm", 1.66, 0.0, 0., 10.);
	NuFit::helpers::par_options astro_index("astro_index", 2.53, 0.0, 0, 10.0);
    NuFit::helpers::par_options muon_norm("muon_norm", 1.2, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options conv_norm("conv_norm", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.0, 0.0, 0.0, 100.0);
    NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -1, 1);
    NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.7, 1.3);
    NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.79, 1.21);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.79, 1.21);
    */

	/** .. package everything */
    /*
	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
    //options.insert(std::pair<std::string, NuFit::helpers::par_options>(log_slope.name, log_slope));
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
    */

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

	/** .. and send to minimizer */
	min.set_options(options);
    min.set_tolerance(1);

    // run global fit to provide seed to profile llh scan when it fails.
    min.fit(false); // ** if true -> get profile LLH errors after minimization from ROOT Minuit2
    
    // only option name, low and high limit is useful for profile scan.
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */
    /*
    NuFit::helpers::par_options astro_norm_prescan("astro_norm", 0, 0.0, 1.2, 2.5);
    NuFit::helpers::par_options astro_index_prescan("astro_index", 0, 0.0, 2.3, 2.9);
    //NuFit::helpers::par_options log_slope_prescan("log_slope", 0, 0.0, -0.1, 0.4);
    NuFit::helpers::par_options prompt_norm_prescan("prompt_norm", 0, 0.0, 0.0, 4);
    NuFit::helpers::par_options muon_norm_prescan("muon_norm", 0, 0.0, 1.0, 1.3);
    NuFit::helpers::par_options conv_norm_prescan("conv_norm", 0, 0.0, 1.0, 2.0);
    NuFit::helpers::par_options cr_index_prescan("delta_cr", 0, 0.0, -0.1, 0.1);
    NuFit::helpers::par_options dom_eff_prescan("dom_efficiency", 0, 0.0, 0.7, 1.1);
    NuFit::helpers::par_options scattering_prescan("scattering", 0, 0.0, 0.9, 1.1);
    NuFit::helpers::par_options absorption_prescan("absorption", 0, 0.0, 0.9, 1.1);
    NuFit::helpers::par_options holeicep0_prescan("holeicep0", 0, 0.0, -2.0, 1.0);
    NuFit::helpers::par_options holeicep1_prescan("holeicep1", 0, 0.0, -0.2, 0.2);
    NuFit::helpers::par_options selfveto_prescan("selfveto", 0, 0.0, 50.0, 3000.0);
    NuFit::helpers::par_options hadronicinteraction_prescan("hadronicinteraction", 0, 0.0, 0, 1);

    // element in map will be sorted automatically by string value
	std::map<std::string, NuFit::helpers::par_options> prescan_options;
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption_prescan.name, absorption_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index_prescan.name, astro_index_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_prescan.name, astro_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm_prescan.name, conv_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index_prescan.name, cr_index_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff_prescan.name, dom_eff_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(hadronicinteraction_prescan.name, hadronicinteraction_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep0_prescan.name, holeicep0_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep1_prescan.name, holeicep1_prescan));
	//prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(log_slope.name, log_slope));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm_prescan.name, muon_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm_prescan.name, prompt_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering_prescan.name, scattering_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(selfveto_prescan.name, selfveto_prescan));
    */

    NuFit::helpers::par_options astro_norm_0_prescan("astro_norm_0", 0, 10.0, 7.0, 10.0);
    NuFit::helpers::par_options astro_norm_1_prescan("astro_norm_1", 0, 10.0, 7.0, 10.0);
    NuFit::helpers::par_options astro_norm_2_prescan("astro_norm_2", 0, 10.0, 7.0, 10.0);
    NuFit::helpers::par_options astro_norm_3_prescan("astro_norm_3", 0, 0.02, 0.0, 3.0);
    NuFit::helpers::par_options astro_norm_4_prescan("astro_norm_4", 0, 0.00, 0.0, 3.0);
    NuFit::helpers::par_options astro_norm_5_prescan("astro_norm_5", 0, 3.24, 0.0, 6.0);
    NuFit::helpers::par_options astro_norm_6_prescan("astro_norm_6", 0, 3.08, 0.0, 6.0);
    NuFit::helpers::par_options astro_norm_7_prescan("astro_norm_7", 0, 5.02, 3.0, 7.0);
    NuFit::helpers::par_options astro_norm_8_prescan("astro_norm_8", 0, 2.32, 1.0, 3.5);
    NuFit::helpers::par_options astro_norm_9_prescan("astro_norm_9", 0, 1.80, 1.0, 3.5);
    NuFit::helpers::par_options astro_norm_10_prescan("astro_norm_10", 0, 0.20, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_11_prescan("astro_norm_11", 0, 0.39, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_12_prescan("astro_norm_12", 0, 0.48, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_13_prescan("astro_norm_13", 0, 0.00, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_14_prescan("astro_norm_14", 0, 0.00, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_15_prescan("astro_norm_15", 0, 0.00, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_16_prescan("astro_norm_16", 0, 0.00, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_17_prescan("astro_norm_17", 0, 0.00, 0.0, 1.5);
    NuFit::helpers::par_options astro_norm_18_prescan("astro_norm_18", 0, 1.80, 1.2, 1.5);
    NuFit::helpers::par_options prompt_norm_prescan("prompt_norm", 0, 0.0, 0.0, 4);
    NuFit::helpers::par_options muon_norm_prescan("muon_norm", 0, 0.0, 1.0, 1.3);
    NuFit::helpers::par_options conv_norm_prescan("conv_norm", 0, 0.0, 1.0, 2.0);
    NuFit::helpers::par_options cr_index_prescan("delta_cr", 0, 0.0, -0.1, 0.1);
    NuFit::helpers::par_options dom_eff_prescan("dom_efficiency", 0, 0.0, 0.7, 1.1);
    NuFit::helpers::par_options scattering_prescan("scattering", 0, 0.0, 0.9, 1.1);
    NuFit::helpers::par_options absorption_prescan("absorption", 0, 0.0, 0.9, 1.1);
    NuFit::helpers::par_options holeicep0_prescan("holeicep0", 0, 0.0, -2.0, 1.0);
    NuFit::helpers::par_options holeicep1_prescan("holeicep1", 0, 0.0, -0.2, 0.2);
    NuFit::helpers::par_options selfveto_prescan("selfveto", 0, 0.0, 50.0, 3000.0);
    NuFit::helpers::par_options hadronicinteraction_prescan("hadronicinteraction", 0, 0.0, 0, 1);

    // element in map will be sorted automatically by string value
	std::map<std::string, NuFit::helpers::par_options> prescan_options;
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption_prescan.name, absorption_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_0_prescan.name, astro_norm_0_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_1_prescan.name, astro_norm_1_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_2_prescan.name, astro_norm_2_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_3_prescan.name, astro_norm_3_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_4_prescan.name, astro_norm_4_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_5_prescan.name, astro_norm_5_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_6_prescan.name, astro_norm_6_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_7_prescan.name, astro_norm_7_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_8_prescan.name, astro_norm_8_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_9_prescan.name, astro_norm_9_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_10_prescan.name, astro_norm_10_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_11_prescan.name, astro_norm_11_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_12_prescan.name, astro_norm_12_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_13_prescan.name, astro_norm_13_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_14_prescan.name, astro_norm_14_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_15_prescan.name, astro_norm_15_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_16_prescan.name, astro_norm_16_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_17_prescan.name, astro_norm_17_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm_18_prescan.name, astro_norm_18_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm_prescan.name, conv_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index_prescan.name, cr_index_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff_prescan.name, dom_eff_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(hadronicinteraction_prescan.name, hadronicinteraction_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep0_prescan.name, holeicep0_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeicep1_prescan.name, holeicep1_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm_prescan.name, muon_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm_prescan.name, prompt_norm_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering_prescan.name, scattering_prescan));
	prescan_options.insert(std::pair<std::string, NuFit::helpers::par_options>(selfveto_prescan.name, selfveto_prescan));
    int variable_index = 0;
    for (std::map<std::string, NuFit::helpers::par_options>::iterator it=prescan_options.begin(); it!=prescan_options.end(); it++)
    {
        if (variable_index == jobid2)
        {
            int nsteps = 20;
            double limit_low = it->second.limit_low;
            double limit_high = it->second.limit_high; 
            double ds = (limit_high-limit_low) / (nsteps-1);
            if (nsteps % njobs1 != 0) {
                std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
                return 1;
            }
            int nsteps_job1 = nsteps / njobs1;
            double min_tj = limit_low + jobid1 * nsteps_job1 * ds;
            double max_tj = min_tj + (nsteps_job1-1) * ds;

            //scan
            NuFit::helpers::scan_options s(it->first, nsteps_job1, min_tj, max_tj);

            std::map<std::string, NuFit::helpers::scan_options> scan_profile;
            scan_profile.insert(std::pair<std::string, NuFit::helpers::scan_options>(s.name, s));
	        min.scan_profile_llh(outdir+std::string("outfile_profile_llh_")+s.name+std::string("_part_")+std::to_string(jobid1)+std::string(".txt"), scan_profile);
        }
        variable_index += 1;
    }

    return 0;
}
