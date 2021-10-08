// [0] Intro {{{1
// ----------------------------------------------------------------------------
/*
 Project title: "Markov health transition probabilities"

 (c) Juergen Jung, 2020
 Towson University, Towson
 jjung@towson.edu
 ----------------------------------------------------------------------------
 This version uses: SPost13 commands and NOT SPost9 commands such as prvalue etc.
 https://jslsoc.sitehost.iu.edu/spost13.htm
 https://jslsoc.sitehost.iu.edu/web_spost13/sp13_faqspost9.htm
 ----------------------------------------------------------------------------

 Uses: Rand HRS Data and MEPS data for 2000, 2002, 2004, 2006, and 2008
 HRS is only available for every other year. So we limit ourselves to years where both
 the HRS and MEPS data are available. The transition probabilities are therefore
 2-year transition probabilities.

 ----------------------------------------------------------------------------
 We use 5 states of health

 Use mgen to generate predictions with confidence intervals for certain
 characteristics. That is we 'evaluate' the predictions at certain population
 characteristics.
 ----------------------------------------------------------------------------
*/
// }}}1
capture log close
log using "$DoDir/LogFiles/200_oLogitPlot_MEPS_vs_HRS.smcl", replace
timer clear 1
// [1] Initialize {{{1
// ----------------------------------------------------------------------------
clear all
set matsize 100
set more off

// Turn off Stata colors in plots

//set scheme s1manual
//ssc install blindschemes, replace all
set scheme plotplainblind

if missing("$initOn") {
    di "--------------------------"
    di "ERROR: run 000_Init.do first"
    di "--------------------------"
}

local agel = 20
local ageh = 95

// Graph options: MEPS vs HRS
// ----------------------------
global line_meps "msymbol(s) connect(l) mcolor(orange_red*1.75) lcolor(orange_red*1.75) lpattern(dash)"
global CI_meps "color(orange_red*0.3)"

global line_hrs "msymbol(o) connect(l) mcolor(gray*1.1) lcolor(gray*1.1) lpattern(dot)"
global CI_hrs "color(gray*0.2)"
global line_hrs_adj "msymbol(sh) connect(l) mcolor(ebblue*1.4) lcolor(ebblue*1.4) lpattern(dot)"


// Graph options: Males vs. Females
// ---------------------------------
global line_f "msymbol(s) connect(l) mcolor(orange_red*1.5) lcolor(orange_red*1.5) lpattern(dash)"
global CI_f "color(orange_red*0.3)"
// HRS adjusted line
global line_f_adj "msymbol(sh) connect(l) mcolor(orange_red*1.85) lcolor(orange_red*1.85) lpattern(dash)"

global line_m "msymbol(o) connect(l) mcolor(ebblue*1.5) lcolor(ebblue*1.5) lpattern(dot)"
global CI_m "color(ebblue*0.3)"
// HRS adjusted line
global line_m_adj "msymbol(oh) connect(l) mcolor(ebblue*1.85) lcolor(ebblue*1.85) lpattern(dot)"

// Graph options: White vs. Black
// -------------------------------
global line_b0 "msymbol(s) connect(l) mcolor(orange_red*1.5) lcolor(orange_red*1.5) lpattern(dash)"
global CI_b0 "color(orange_red*0.3)"
// HRS adjusted line
global line_b0_adj "msymbol(sh) connect(l) mcolor(orange_red*1.85) lcolor(orange_red*1.85) lpattern(dash)"

global line_b1 "msymbol(o) connect(l) mcolor(ebblue*1.5) lcolor(ebblue*1.5) lpattern(dot)"
global CI_b1 "color(ebblue*0.3)"
// HRS adjusted line
global line_b1_adj "msymbol(oh) connect(l) mcolor(ebblue*1.85) lcolor(ebblue*1.85) lpattern(dot)"

// Graph options: Non-Smoke vs. Smoker
// ------------------------------------
global line_s0 "msymbol(s) connect(l) mcolor(orange_red*1.5) lcolor(orange_red*1.5) lpattern(dash)"
global CI_s0 "color(orange_red*0.3)"
// HRS adjusted line
global line_s0_adj "msymbol(sh) connect(l) mcolor(orange_red*1.85) lcolor(orange_red*1.85) lpattern(dash)"

global line_s1 "msymbol(o) connect(l) mcolor(ebblue*1.5) lcolor(ebblue*1.5) lpattern(dot)"
global CI_s1 "color(ebblue*0.3)"
// HRS adjusted line
global line_s1_adj "msymbol(oh) connect(l) mcolor(ebblue*1.85) lcolor(ebblue*1.85) lpattern(dot)"


//}}}1
// [2] Set HH or Individual level weights {{{1
// ----------------------------------------------------------------------------
// Load data
use "$DatDir/MEPS_HRS_1992_2017_1r.dta", clear
do $DoDir/f_globalLists.do

// 0. Individual level (Head + non Head)
// 1. HH level (Head of HH)

//local I_obs_level "HH"
//local I_obs_level "HIEU"
local I_obs_level "Person"

//local I_HIEU_level = 0

if "`I_obs_level'" == "HH" {
    // Household level
    keep if d_head == 1

    local weight famweight
    local int_weight int_famweight
    local levelObs "Household"
    label variable age "Age of head of household"

    sort famID year
    by famID: gen first = _n == 1
    gen obsID = sum(first)
    label variable obsID "Obs unit: Family"
    order obsID first famID persID year
}
else if "`I_obs_level'" == "HIEU" {
    // HIEU level
    keep if d_head_HIEU == 1

    local weight hieuweight
    local int_weight int_hieuweight
    local levelObs "HIEU"
    label variable age "Age of head of HIEU"

    sort hieuID year
    by hieuID: gen first = _n == 1
    gen obsID = sum(first)
    label variable obsID "Obs unit: HIEU"
    order obsID first hieuID persID year
}
else if "`I_obs_level'" == "Person" {
    // Individual level
    keep if d_head == 1 | d_head == 0
    local weight persweight
    local int_weight int_persweight
    local levelObs "Individual"
    label variable age "Age of individual"

    gen obsID =  persID
    label variable obsID "Obs unit: Individual"

}
xtset persID year
sort persID year

//}}}1
// [3] Estimate Ordered Logit - reg F_healthStatusD healthStatus Xv  {{{1
// ----------------------------------------------------------------------------
// ologit using a lag health state variable in regression.
// This allows for different intercept by lagged health state.
// Controlled for cohort effect

timer on 1

foreach dataName in "HRS" "MEPS" {

    preserve
        if "`dataName'" == "MEPS" {
            keep if d_MEPS == 1
        }
        else {
            keep if d_HRS == 1
        }

        // Run regression for MEPS and HRS separately
        if "`dataName'" == "HRS" {
            ologit F_healthStatusD i.healthStatus i.initHealth $Xv $periodInteraction d_cohort* if age>=`agel' & age<=`ageh' [pweight = `weight'], cluster(persID) difficult iterate(1000)
        }
        else {
            ologit F_healthStatusD i.healthStatus $Xv $periodInteraction d_cohort* if age>=`agel' & age<=`ageh' [pweight = `weight'], cluster(persID) difficult iterate(1000)
        }
        estimates store res_ologit_`dataName'

        local i_h = 0
        foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

            local i_h = `i_h' + 1 // counts up from 1 to 5 health status

            display"-----------------------------------------"
            display `i_h'
            display"-----------------------------------------"

            // Run regression
            // ------------------------------------------------------------------------
            // The Model: P[Excellent,Very Good, Good,Fair,Poor,Dead |  `from_status']
            estimates restore res_ologit_`dataName'

            // All
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h') stub(a_`dataName'_`i_h'_to_)

            // Females, white, with partner
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_female=1) stub(f_`dataName'_`i_h'_to_)

            // Males, white, with partner
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_female=0) stub(m_`dataName'_`i_h'_to_)

            // Non-Black
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_black=0) stub(b0_`dataName'_`i_h'_to_)

            // Black
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_black=1) stub(b1_`dataName'_`i_h'_to_)

            // non-Smoker
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_smoke=0) stub(s0_`dataName'_`i_h'_to_)

            // Smoker
            mgen, at(age=(`agel'(5)`ageh') healthStatus = `i_h' d_smoke=1) stub(s1_`dataName'_`i_h'_to_)



            rename m_`dataName'_`i_h'_to_age agev
            label var agev "Age"

            // Export the results into a .dat file, that can be read by MatLab
            if "`dataName'" == "MEPS" {
                export delimited agev a_`dataName'_`i_h'_to_* ///
                f_`dataName'_`i_h'_to_* m_`dataName'_`i_h'_to_* ///
                b0_`dataName'_`i_h'_to_* b1_`dataName'_`i_h'_to_* ///
                s0_`dataName'_`i_h'_to_* s1_`dataName'_`i_h'_to_* ///
                    in 1/16 using "$StataOutDir/`dataName'_1999-2017_Markov_Age_`i_h'_Dead.csv", delimiter(",") replace
            }
            else {
                export delimited agev a_`dataName'_`i_h'_to_* ///
                f_`dataName'_`i_h'_to_* m_`dataName'_`i_h'_to_* ///
                b0_`dataName'_`i_h'_to_* b1_`dataName'_`i_h'_to_* ///
                s0_`dataName'_`i_h'_to_* s1_`dataName'_`i_h'_to_* ///
                    in 1/16 using "$StataOutDir/`dataName'_1992-2016_Markov_Age_`i_h'_Dead.csv", delimiter(",") replace
            }
                        // Drop these variables
            drop agev
            drop a_`dataName'_`i_h'_to_*
            drop m_`dataName'_`i_h'_to_*
            drop f_`dataName'_`i_h'_to_*
            drop b0_`dataName'_`i_h'_to_*
            drop b1_`dataName'_`i_h'_to_*
            drop s0_`dataName'_`i_h'_to_*
            drop s1_`dataName'_`i_h'_to_*
        }
    restore
}
//}}}1
timer off 1
timer list 1
// [4] Merge All Markov Results Together  {{{1
// ----------------------------------------------------------------------------
foreach i_h in 1 2 3 4 5 {
    di "-----------------------------------"
    di `i_h'
    di "-----------------------------------"

    // MEPS
    clear
    import delimited using "$StataOutDir/MEPS_1999-2017_Markov_Age_`i_h'_Dead.csv", delimiter(",") varnames(1)
    sort agev
    save "$DatDir/MEPS_1999-2017_Age_`i_h'_Dead.dta", replace

    // HRS
    clear
    import delimited using "$StataOutDir/HRS_1992-2016_Markov_Age_`i_h'_Dead.csv", delimiter(",") varnames(1)
    sort agev
    save "$DatDir/HRS_1992-2016_Age_`i_h'_Dead.dta", replace
}

clear
input agev
20
25
30
35
40
45
50
55
60
65
70
75
80
85
90
95
end

foreach i_h in 1 2 3 4 5 {
    sort agev
    merge agev using "$DatDir/MEPS_1999-2017_Age_`i_h'_Dead.dta"
    tab _merge
    drop _merge
}
foreach i_h in 1 2 3 4 5 {
    sort agev
    merge agev using "$DatDir/HRS_1992-2016_Age_`i_h'_Dead.dta"
    tab _merge
    drop _merge
}

// Stata format save
save "$DatDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit.dta", replace

// .csv format save  --> for Python Adjustments and Plots
export delimited using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit.csv", delimiter(",") replace
//}}}1
di "-----------------------------------------------------------------"
di "Now Call Python Script: Python/003_mainAdjustOLogitMarkov.py"
di "-----------------------------------------------------------------"
!rm "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv"
!$AnacondaDir/python $PythonDir/003A_mainAdjustOLogitMarkov.py

while (1) {
    capture insheet using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv", clear
    if _rc == 0 continue, break
    sleep 10000  // Wait 10 seconds
}
// [5] Make Plots - Data Cleaning Program {{{1
//-----------------------------------------------------------------------------

capture program drop f_prepareData
program define f_prepareData

// Cleaning steps: Express all variables as percentages, add labels, etc.

    capture drop v1

    qui foreach v of varlist _all {
        replace `v' = `v' * 100
        label variable `v' %
    }
    replace agev = agev/100
    label variable agev "Age"

    // To get good y-scales on graphs we need to drop observations that we do
    // not use in the graphs

    // Drop obs smaller 50 for HRS data
    foreach x of varlist *_hrs_* {
        replace `x' = . if agev < 50
    }

    // Drop obs larger 70 MEPS data
    foreach x of varlist *_meps_* {
        replace `x' = . if agev > 70
    }

    gen age_mepsv = agev if agev <=70
    label variable age_mepsv "Age"

    gen age_hrsv  = agev if agev >=50
    label variable age_hrsv "Age"
end
//}}}1
// [5.A] Plot Adjustd Matrix - All: Graph MEPS and HRS vs. HRS-adjusted Markov Probabilities Separately {{{1
//-----------------------------------------------------------------------------
// You need to run the python script: Python/003_mainAdjustMarkovPeriodsDead.py
// which uses "$DatDir/MEPS_HRS_1992-2017_Age_All_Dead.dta" and adjust the HRS markov
// transition matrices from 2-year to 1-year matrices using an
// Algorithm based on: J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing
// Cycle Lengths in State-Transition Models: Challenges and Solutions,” Medical
// Decision Making, vol. 36, no. 8, pp. 952–964, 2016.

clear
import delimited using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv", delimiter(",") varnames(1)

// Call our cleaning step
f_prepareData

// Start plot
local h_from = 1
foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

    local h_to   = 1
    foreach to_status in "Excellent" "Very Good" "Good" "Fair" "Poor" "Dead" {

        if "`to_status'" != "Dead" {

            graph twoway    (rarea a_meps_`h_from'_to_ll`h_to' a_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_meps) ///
                        (rarea  a_hrs_`h_from'_to_ll`h_to' a_hrs_`h_from'_to_ul`h_to' age_hrsv, $CI_hrs) ///
                        (connected a_meps_`h_from'_to_pr`h_to' age_mepsv, $line_meps)  ///
                        (connected  a_hrs_`h_from'_to_pr`h_to' age_hrsv, $line_hrs)  ///
                        (connected  a_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_hrs_adj), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(off) name(fig`h_to', replace)
        }
        else {

            // Make last plot with legend
            graph twoway    (rarea a_meps_`h_from'_to_ll`h_to' a_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_meps) ///
                        (rarea  a_hrs_`h_from'_to_ll`h_to' a_hrs_`h_from'_to_ul`h_to' age_hrsv, $CI_hrs) ///
                        (connected a_meps_`h_from'_to_pr`h_to' age_mepsv, $line_meps)  ///
                        (connected  a_hrs_`h_from'_to_pr`h_to' age_hrsv, $line_hrs)  ///
                        (connected  a_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_hrs_adj), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(size(vsmall) ring(0) position(10) col(1) ///
                    order(3 "MEPS" 4 "HRS" 5 "HRS-Adj.")) name(fig`h_to', replace)
                    //order(3 "MEPS" 1 "95% CI" 4 "HRS" 2 "95% CI" 5 "HRS-Adj."))
        }

        local h_to = `h_to' + 1
    }

    // Make an empty plot, just to get legend out!!
    //graph twoway (rarea a_meps_1_to_ll1 a_meps_1_to_ul1 age_mepsv, $CI_meps) ///
    //             (rarea a_hrs_1_to_ll1 a_hrs_1_to_ul1 age_hrsv, $CI_hrs) ///
    //             (connected a_meps_1_to_pr1 age_mepsv, $line_meps)  ///
    //             (connected a_hrs_1_to_pr1 age_hrsv, $line_hrs)  ///
    //             (connected a_hrs_1_to_pr1_adj age_hrsv, $line_hrs_adj), ///
    //        yscale(off) xscale(off)  ///
    //        legend(ring(0) position(0) col(2) ///
    //        order(3 "MEPS" 1 "95% CI" 4 "HRS" 2 "95% CI" 5 "HRS-Adj.")) name(jleg, replace)
    //
    //gr_edit .plotregion1.draw_view.setstyle, style(no)

    // -----------------------------------------------
    //graph combine fig1 fig2 fig3 fig4 fig5 fig6 jleg, ///
    graph combine fig1 fig2 fig3 fig4 fig5 fig6, ///
        xcommon col(2) imargin(1 1 1 1) graphregion(margin(l=1 r=1)) ///
        iscale(0.8) ysize(20.0) xsize(16.0) ///
        name(ffig, replace)
    //
    graph export "$TexDir/MEPSvsHRS_TransGraph`h_from'_Dead_OLogit_initHealth.svg", replace
    //
    local h_from = `h_from' + 1
}
//}}}1
// [5.B] Plot Adjusted Matrix - MALES vs. Females: Graph MEPS and HRS Markov Probabilities Separately {{{1
//-------------------------------------------------------------------------
// You need to run the python script: Python/003_mainAdjustMarkovPeriodsDead.py
// which uses "$DatDir/MEPS_HRS_1992-2017_Age_All_Dead.dta" and adjust the HRS markov
// transition matrices from 2-year to 1-year matrices using an
// Algorithm based on: J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing
// Cycle Lengths in State-Transition Models: Challenges and Solutions,” Medical
// Decision Making, vol. 36, no. 8, pp. 952–964, 2016.

clear
import delimited using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv", delimiter(",") varnames(1)

// Call our cleaning step
f_prepareData

// Start plot
local h_from = 1
foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

    local h_to   = 1
    foreach to_status in "Excellent" "Very Good" "Good" "Fair" "Poor" "Dead" {

        if "`to_status'" != "Dead" {
            graph twoway (rarea f_meps_`h_from'_to_ll`h_to' f_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_f) ///
                     (rarea m_meps_`h_from'_to_ll`h_to' m_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_m) ///
                     (connected f_meps_`h_from'_to_pr`h_to' age_mepsv, $line_f)  ///
                     (connected m_meps_`h_from'_to_pr`h_to' age_mepsv, $line_m)  ///
                     (connected  f_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_f_adj)  ///
                     (connected  m_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_m_adj), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(off) name(fig`h_to', replace)

        }
        else {
            // Make last plot with legend
            graph twoway (rarea f_meps_`h_from'_to_ll`h_to' f_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_f) ///
                     (rarea m_meps_`h_from'_to_ll`h_to' m_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_m) ///
                     (connected f_meps_`h_from'_to_pr`h_to' age_mepsv, $line_f)  ///
                     (connected m_meps_`h_from'_to_pr`h_to' age_mepsv, $line_m)  ///
                     (connected  f_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_f_adj)  ///
                     (connected  m_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_m_adj), ///
                xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                legend(size(vsmall) ring(0) position(10) col(1) ///
                order(3 "Fem MEPS" 4 "Male MEPS" 5 "Fem HRS" 6 "Male HRS")) ///
                name(fig`h_to', replace)
        }

        local h_to = `h_to' + 1
    }

    //// Make an empty plot, just to get legend out!!
    //graph twoway (rarea f_meps_1_to_ll1 f_meps_1_to_ul1 age_mepsv, $CI_f) ///
    //             (rarea m_meps_1_to_ll1 m_meps_1_to_ul1 age_mepsv, $CI_m) ///
    //             (connected f_meps_1_to_pr1 age_mepsv, $line_f)  ///
    //             (connected m_meps_1_to_pr1 age_mepsv, $line_m)  ///
    //             (connected f_hrs_1_to_pr1_adj age_hrsv, $line_f_adj)  ///
    //             (connected m_hrs_1_to_pr1_adj age_hrsv, $line_m_adj), ///
    //        yscale(off) xscale(off)  ///
    //        legend(ring(0) position(0) col(2) ///
    //        order(3 "Females-MEPS" 1 "95% CI" 4 "Males-MEPS" 2 "95% CI" 5 "Females-HRS-Adj." 6 "Males-HRS-Adj.")) name(jleg, replace)

    //gr_edit .plotregion1.draw_view.setstyle, style(no)

    // -----------------------------------------------
    //graph combine fig1 fig2 fig3 fig4 fig5 fig6 jleg, ///
    graph combine fig1 fig2 fig3 fig4 fig5 fig6, ///
        xcommon col(2) imargin(1 1 1 1) graphregion(margin(l=1 r=1)) ///
        iscale(0.8) ysize(20.0) xsize(16.0) ///
        name(ffig,replace)
    //
    graph export "$TexDir/MEPSvsHRS_FemalesVsMales_TransGraph`h_from'_Dead_OLogit_initHealth.svg", replace
    //
    local h_from = `h_from' + 1
}
//}}}1
// [5.C] Plot Adjusted Matrix - White (Non-Black) vs. Black: Graph MEPS and HRS Probabilities Separately {{{1
//-------------------------------------------------------------------------
// You need to run the python script: Python/003_mainAdjustMarkovPeriodsDead.py
// which uses "$DatDir/MEPS_HRS_1992-2017_Age_All_Dead.dta" and adjust the HRS markov
// transition matrices from 2-year to 1-year matrices using an
// Algorithm based on: J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing
// Cycle Lengths in State-Transition Models: Challenges and Solutions,” Medical
// Decision Making, vol. 36, no. 8, pp. 952–964, 2016.

clear
import delimited using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv", delimiter(",") varnames(1)

// Call our cleaning step
f_prepareData

// Start plot
local h_from = 1
foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

    local h_to   = 1
    foreach to_status in "Excellent" "Very Good" "Good" "Fair" "Poor" "Dead" {

        if "`to_status'" != "Dead" {
            graph twoway (rarea b0_meps_`h_from'_to_ll`h_to' b0_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_b0) ///
                     (rarea b1_meps_`h_from'_to_ll`h_to' b1_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_b1) ///
                     (connected b0_meps_`h_from'_to_pr`h_to' age_mepsv, $line_b0)  ///
                     (connected b1_meps_`h_from'_to_pr`h_to' age_mepsv, $line_b1)  ///
                     (connected  b0_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_b0_adj)  ///
                     (connected  b1_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_b1_adj), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(off) name(fig`h_to', replace)

        }
        else {
            // Make last plot with legend
            graph twoway (rarea b0_meps_`h_from'_to_ll`h_to' b0_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_b0) ///
                     (rarea b1_meps_`h_from'_to_ll`h_to' b1_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_b1) ///
                     (connected b0_meps_`h_from'_to_pr`h_to' age_mepsv, $line_b0)  ///
                     (connected b1_meps_`h_from'_to_pr`h_to' age_mepsv, $line_b1)  ///
                     (connected  b0_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_b0_adj)  ///
                     (connected  b1_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_b1_adj), ///
                xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                legend(size(vsmall) ring(0) position(10) col(1) ///
                order(3 "White MEPS" 4 "Black MEPS" 5 "White HRS" 6 "Black HRS")) ///
                name(fig`h_to', replace)
        }

        local h_to = `h_to' + 1
    }

    //// Make an empty plot, just to get legend out!!
    //graph twoway (rarea b0_meps_1_to_ll1 b0_meps_1_to_ul1 age_mepsv, $CI_b0) ///
    //             (rarea b1_meps_1_to_ll1 b1_meps_1_to_ul1 age_mepsv, $CI_b1) ///
    //             (connected b0_meps_1_to_pr1 age_mepsv, $line_b0)  ///
    //             (connected b1_meps_1_to_pr1 age_mepsv, $line_b1)  ///
    //             (connected b0_hrs_1_to_pr1_adj age_hrsv, $line_b0_adj)  ///
    //             (connected b1_hrs_1_to_pr1_adj age_hrsv, $line_b1_adj), ///
    //        yscale(off) xscale(off)  ///
    //        legend(ring(0) position(0) col(2) ///
    //        order(3 "White-MEPS" 1 "95% CI" 4 "Black-MEPS" 2 "95% CI" 5 "White-HRS-Adj." 6 "Black-HRS-Adj.")) name(jleg, replace)

    //gr_edit .plotregion1.draw_view.setstyle, style(no)

    // -----------------------------------------------
    //graph combine fig1 fig2 fig3 fig4 fig5 fig6 jleg, ///
    graph combine fig1 fig2 fig3 fig4 fig5 fig6, ///
        xcommon col(2) imargin(1 1 1 1) graphregion(margin(l=1 r=1)) ///
        iscale(0.8) ysize(20.0) xsize(16.0) ///
        name(ffig,replace)
    //
    graph export "$TexDir/MEPSvsHRS_WhiteVsBlack_TransGraph`h_from'_Dead_OLogit_initHealth.svg", replace
    //
    local h_from = `h_from' + 1
}
//}}}1
// [5.D] Plot Adjusted Matrix - Smoker vs. Non-Smoker: Graph MEPS and HRS Probabilities Separately {{{1
//-------------------------------------------------------------------------
// You need to run the python script: Python/003_mainAdjustMarkovPeriodsDead.py
// which uses "$DatDir/MEPS_HRS_1992-2017_Age_All_Dead.dta" and adjust the HRS markov
// transition matrices from 2-year to 1-year matrices using an
// Algorithm based on: J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing
// Cycle Lengths in State-Transition Models: Challenges and Solutions,” Medical
// Decision Making, vol. 36, no. 8, pp. 952–964, 2016.

clear
import delimited using "$StataOutDir/MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv", delimiter(",") varnames(1)

// Call our cleaning step
f_prepareData

// Start plot
local h_from = 1
foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

    local h_to   = 1
    foreach to_status in "Excellent" "Very Good" "Good" "Fair" "Poor" "Dead" {

        if "`to_status'" != "Dead" {
            graph twoway (rarea s0_meps_`h_from'_to_ll`h_to' s0_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_s0) ///
                     (rarea s1_meps_`h_from'_to_ll`h_to' s1_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_s1) ///
                     (connected s0_meps_`h_from'_to_pr`h_to' age_mepsv, $line_s0)  ///
                     (connected s1_meps_`h_from'_to_pr`h_to' age_mepsv, $line_s1)  ///
                     (connected  s0_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_s0_adj)  ///
                     (connected  s1_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_s1_adj), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(off) name(fig`h_to', replace)

        }
        else {
            // Make last plot with legend
            graph twoway (rarea s0_meps_`h_from'_to_ll`h_to' s0_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_s0) ///
                     (rarea s1_meps_`h_from'_to_ll`h_to' s1_meps_`h_from'_to_ul`h_to' age_mepsv, $CI_s1) ///
                     (connected s0_meps_`h_from'_to_pr`h_to' age_mepsv, $line_s0)  ///
                     (connected s1_meps_`h_from'_to_pr`h_to' age_mepsv, $line_s1)  ///
                     (connected  s0_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_s0_adj)  ///
                     (connected  s1_hrs_`h_from'_to_pr`h_to'_adj age_hrsv, $line_s1_adj), ///
                xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                legend(size(vsmall) ring(0) position(10) col(1) ///
                order(3 "Non-S MEPS" 4 "Smoke MEPS" 5 "Non-S HRS" 6 "Smoke HRS")) ///
                name(fig`h_to', replace)
        }

        local h_to = `h_to' + 1
    }

    //// Make an empty plot, just to get legend out!!
    //graph twoway (rarea b0_meps_1_to_ll1 b0_meps_1_to_ul1 age_mepsv, $CI_b0) ///
    //             (rarea b1_meps_1_to_ll1 b1_meps_1_to_ul1 age_mepsv, $CI_b1) ///
    //             (connected b0_meps_1_to_pr1 age_mepsv, $line_b0)  ///
    //             (connected b1_meps_1_to_pr1 age_mepsv, $line_b1)  ///
    //             (connected b0_hrs_1_to_pr1_adj age_hrsv, $line_b0_adj)  ///
    //             (connected b1_hrs_1_to_pr1_adj age_hrsv, $line_b1_adj), ///
    //        yscale(off) xscale(off)  ///
    //        legend(ring(0) position(0) col(2) ///
    //        order(3 "White-MEPS" 1 "95% CI" 4 "Black-MEPS" 2 "95% CI" 5 "White-HRS-Adj." 6 "Black-HRS-Adj.")) name(jleg, replace)

    //gr_edit .plotregion1.draw_view.setstyle, style(no)

    // -----------------------------------------------
    //graph combine fig1 fig2 fig3 fig4 fig5 fig6 jleg, ///
    graph combine fig1 fig2 fig3 fig4 fig5 fig6, ///
        xcommon col(2) imargin(1 1 1 1) graphregion(margin(l=1 r=1)) ///
        iscale(0.8) ysize(20.0) xsize(16.0) ///
        name(ffig,replace)
    //
    graph export "$TexDir/MEPSvsHRS_NonSmokerVsSmoker_TransGraph`h_from'_Dead_OLogit_initHealth.svg", replace
    //
    local h_from = `h_from' + 1
}
//}}}1
log close
