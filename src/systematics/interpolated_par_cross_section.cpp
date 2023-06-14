#include "../../include/systematics/interpolated_par_cross_section.h"

NuFit::interpolated_par_cross_section::interpolated_par_cross_section(std::string par_name_, std::vector<std::string> analysis_names_, const std::map<std::string, NuFit::hists*> &map_analyses): interpolated_par(par_name_, analysis_names_)

{
    // need to know what analyses are implemented -> check /src/analysis.cpp
    std::string dir=basedir+std::string("sys/txt/");
    std::string baselinedir=basedir+std::string("baseline/");    

    std::string dir_muon=basedir+std::string("sys/txt/");
    std::string baselinedir_muon=basedir+std::string("baseline/");

    // cascade event selection first
    std::string cascade("cascade_all");
    if (std::find(analysis_names.begin(), analysis_names.end(), cascade) != analysis_names.end())
    {    
        std::vector<double> binsx = map_analyses.at(cascade)->get_binsx();
        std::vector<double> binsy = map_analyses.at(cascade)->get_binsy();
        std::vector<double> binsz = map_analyses.at(cascade)->get_binsz();

        NuFit::interpolated_sys *nue_cross_section_cascade = new interpolated_sys(par_name, cascade, "NuE", binsx, binsy, binsz);
        nue_cross_section_cascade -> add_simulated_point(5, dir + std::string("step5/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(100, dir + std::string("step100/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(250, dir + std::string("step250/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(500, dir + std::string("step500/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(750, dir + std::string("step750/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(1000, baselinedir + std::string("nue_cascade.txt"),true); //thi is baseline hist
        nue_cross_section_cascade -> add_simulated_point(1001, baselinedir + std::string("nue_cascade.txt")); 
        nue_cross_section_cascade -> add_simulated_point(1500, dir + std::string("step1500/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(2000, dir + std::string("step2000/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(2500, dir + std::string("step2500/nue_cascade.txt"));
        nue_cross_section_cascade -> add_simulated_point(3000, dir + std::string("step3000/nue_cascade.txt"));
        nue_cross_section_cascade -> create_correction_functions(true);

        NuFit::interpolated_sys *numu_cross_section_cascade = new interpolated_sys(par_name, cascade, "NuMu", binsx, binsy, binsz); 
        numu_cross_section_cascade -> add_simulated_point(5, dir + std::string("step5/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(100, dir + std::string("step100/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(250, dir + std::string("step250/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(500, dir + std::string("step500/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(750, dir + std::string("step750/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(1000, baselinedir + std::string("numu_cascade.txt"),true); //thi is baseline hist
        numu_cross_section_cascade -> add_simulated_point(1001, baselinedir + std::string("numu_cascade.txt")); 
        numu_cross_section_cascade -> add_simulated_point(1500, dir + std::string("step1500/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(2000, dir + std::string("step2000/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(2500, dir + std::string("step2500/numu_cascade.txt"));
        numu_cross_section_cascade -> add_simulated_point(3000, dir + std::string("step3000/numu_cascade.txt"));
        numu_cross_section_cascade -> create_correction_functions(true);

        NuFit::interpolated_sys *nutau_cross_section_cascade = new interpolated_sys(par_name, cascade, "NuTau", binsx, binsy, binsz);
        nutau_cross_section_cascade -> add_simulated_point(5, dir + std::string("step5/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(100, dir + std::string("step100/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(250, dir + std::string("step250/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(500, dir + std::string("step500/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(750, dir + std::string("step750/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(1000, baselinedir + std::string("nutau_cascade.txt"),true); //thi is baseline hist
        nutau_cross_section_cascade -> add_simulated_point(1001, baselinedir + std::string("nutau_cascade.txt")); 
        nutau_cross_section_cascade -> add_simulated_point(1500, dir + std::string("step1500/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(2000, dir + std::string("step2000/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(2500, dir + std::string("step2500/nutau_cascade.txt"));
        nutau_cross_section_cascade -> add_simulated_point(3000, dir + std::string("step3000/nutau_cascade.txt"));
        nutau_cross_section_cascade -> create_correction_functions(true);

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> cascade_corrections;
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_cross_section_cascade));
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_cross_section_cascade));
        cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_cross_section_cascade));
        //cascade_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_cross_section_cascade));

        // and now add
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(cascade, cascade_corrections));
    }

    // muon event selection
    std::string muon("muon");
    if (std::find(analysis_names.begin(), analysis_names.end(), muon) != analysis_names.end())
    {    
        std::vector<double> binsx = map_analyses.at(muon)->get_binsx();
        std::vector<double> binsy = map_analyses.at(muon)->get_binsy();
        std::vector<double> binsz = map_analyses.at(muon)->get_binsz();

        NuFit::interpolated_sys *nue_cross_section_muon = new interpolated_sys(par_name, muon, "NuE", binsx, binsy, binsz);
        nue_cross_section_muon -> add_simulated_point(5, dir + std::string("step5/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(100, dir + std::string("step100/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(250, dir + std::string("step250/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(500, dir + std::string("step500/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(750, dir + std::string("step750/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(1000, baselinedir + std::string("nue_muon.txt"),true); //thi is baseline hist
        nue_cross_section_muon -> add_simulated_point(1001, baselinedir + std::string("nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(1500, dir + std::string("step1500/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(2000, dir + std::string("step2000/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(2500, dir + std::string("step2500/nue_muon.txt"));
        nue_cross_section_muon -> add_simulated_point(3000, dir + std::string("step3000/nue_muon.txt"));
        nue_cross_section_muon -> create_correction_functions(true);

        NuFit::interpolated_sys *numu_cross_section_muon = new interpolated_sys(par_name, muon, "NuMu", binsx, binsy, binsz); 
        numu_cross_section_muon -> add_simulated_point(5, dir + std::string("step5/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(100, dir + std::string("step100/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(250, dir + std::string("step250/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(500, dir + std::string("step500/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(750, dir + std::string("step750/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(1000, baselinedir + std::string("numu_muon.txt"),true); //thi is baseline hist
        numu_cross_section_muon -> add_simulated_point(1001, baselinedir + std::string("numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(1500, dir + std::string("step1500/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(2000, dir + std::string("step2000/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(2500, dir + std::string("step2500/numu_muon.txt"));
        numu_cross_section_muon -> add_simulated_point(3000, dir + std::string("step3000/numu_muon.txt"));
        numu_cross_section_muon -> create_correction_functions(true);

        NuFit::interpolated_sys *nutau_cross_section_muon = new interpolated_sys(par_name, muon, "NuTau", binsx, binsy, binsz);
        nutau_cross_section_muon -> add_simulated_point(5, dir + std::string("step5/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(100, dir + std::string("step100/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(250, dir + std::string("step250/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(500, dir + std::string("step500/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(750, dir + std::string("step750/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(1000, baselinedir + std::string("nutau_muon.txt"),true); //thi is baseline hist
        nutau_cross_section_muon -> add_simulated_point(1001, baselinedir + std::string("nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(1500, dir + std::string("step1500/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(2000, dir + std::string("step2000/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(2500, dir + std::string("step2500/nutau_muon.txt"));
        nutau_cross_section_muon -> add_simulated_point(3000, dir + std::string("step3000/nutau_muon.txt"));
        nutau_cross_section_muon -> create_correction_functions(true);

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> muon_corrections;
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_cross_section_muon));
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_cross_section_muon));
        muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_cross_section_muon));
        //muon_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_cross_section_muon));

        // store 
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(muon, muon_corrections));
    }

    // hybrid event selection
    std::string hybrid("hybrid");
    if (std::find(analysis_names.begin(), analysis_names.end(), hybrid) != analysis_names.end())
    {
        std::vector<double> binsx = map_analyses.at(hybrid)->get_binsx();
        std::vector<double> binsy = map_analyses.at(hybrid)->get_binsy();
        std::vector<double> binsz = map_analyses.at(hybrid)->get_binsz();    

        NuFit::interpolated_sys *nue_cross_section_hybrid = new interpolated_sys(par_name, hybrid, "NuE", binsx, binsy, binsz);
        nue_cross_section_hybrid -> add_simulated_point(5, dir + std::string("step5/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(100, dir + std::string("step100/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(250, dir + std::string("step250/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(500, dir + std::string("step500/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(750, dir + std::string("step750/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(1000, baselinedir + std::string("nue_hybrid.txt"),true); //thi is baseline hist
        nue_cross_section_hybrid -> add_simulated_point(1001, baselinedir + std::string("nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(1500, dir + std::string("step1500/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(2000, dir + std::string("step2000/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(2500, dir + std::string("step2500/nue_hybrid.txt"));
        nue_cross_section_hybrid -> add_simulated_point(3000, dir + std::string("step3000/nue_hybrid.txt"));
        nue_cross_section_hybrid -> create_correction_functions(true);

        NuFit::interpolated_sys *numu_cross_section_hybrid = new interpolated_sys(par_name, hybrid, "NuMu", binsx, binsy, binsz); 
        numu_cross_section_hybrid -> add_simulated_point(5, dir + std::string("step5/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(100, dir + std::string("step100/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(250, dir + std::string("step250/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(500, dir + std::string("step500/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(750, dir + std::string("step750/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(1000, baselinedir + std::string("numu_hybrid.txt"),true); //thi is baseline hist
        numu_cross_section_hybrid -> add_simulated_point(1001, baselinedir + std::string("numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(1500, dir + std::string("step1500/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(2000, dir + std::string("step2000/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(2500, dir + std::string("step2500/numu_hybrid.txt"));
        numu_cross_section_hybrid -> add_simulated_point(3000, dir + std::string("step3000/numu_hybrid.txt"));
        numu_cross_section_hybrid -> create_correction_functions(true);

        NuFit::interpolated_sys *nutau_cross_section_hybrid = new interpolated_sys(par_name, hybrid, "NuTau", binsx, binsy, binsz);
        nutau_cross_section_hybrid -> add_simulated_point(5, dir + std::string("step5/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(100, dir + std::string("step100/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(250, dir + std::string("step250/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(500, dir + std::string("step500/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(750, dir + std::string("step750/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(1000, baselinedir + std::string("nutau_hybrid.txt"),true); //thi is baseline hist
        nutau_cross_section_hybrid -> add_simulated_point(1001, baselinedir + std::string("nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(1500, dir + std::string("step1500/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(2000, dir + std::string("step2000/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(2500, dir + std::string("step2500/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> add_simulated_point(3000, dir + std::string("step3000/nutau_hybrid.txt"));
        nutau_cross_section_hybrid -> create_correction_functions(true);

        // map from flavor to correction
        std::unordered_map<std::string, NuFit::interpolated_sys *> hybrid_corrections;
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuE", nue_cross_section_hybrid));
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuMu", numu_cross_section_hybrid));
        hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("NuTau", nutau_cross_section_hybrid));
        //hybrid_corrections.insert(std::pair<std::string, NuFit::interpolated_sys *>("Muon", nue_cross_section_hybrid));

        // store
        binfits.insert(std::pair<std::string, std::unordered_map<std::string, NuFit::interpolated_sys *>>(hybrid, hybrid_corrections));
    }
    
    // specify what histograms are effected by this systematic
    components.clear();
    components.push_back(std::string("Conv"));
    components.push_back(std::string("Prompt"));
    components.push_back(std::string("Astro"));
    return;
}

NuFit::interpolated_par_cross_section::~interpolated_par_cross_section() { }

