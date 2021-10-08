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

 The default is to use the mean. The default for dummy variables is of course question-
 able. Nevertheless I used mean values of education dummies and whether people live
 in partnerships or not. Another way would be to make separate graphs for these sub-
 populations.

 This produces smoothers graphs over the age range, than the previous method, where
 we just predict and then report the mean for each age group.
 ----------------------------------------------------------------------------
*/
// }}}1
capture log close
log using "$DoDir/LogFiles/100A_Table_oLogit_2TimePeriods_initialHealth.smcl", replace
// [1] Initialize {{{1
// ----------------------------------------------------------------------------
clear all
set matsize 100
set more off

// Turn off Stata colors in plots
set scheme s1manual

if missing("$initOn") {
    di "--------------------------"
    di "ERROR: run 000_Init.do first"
    di "--------------------------"
}

local agel = 20
local ageh = 95

// Graph options
global jLine_f "msymbol(s) connect(l) mcolor(orange_red*1.75) lcolor(orange_red*1.75) lpattern(dash)"
global jCI_f "color(orange_red*0.3)"

global jLine_m "msymbol(o) connect(l) mcolor(ebblue*1.2) lcolor(ebblue*1.2) lpattern(dot)"
global jCI_m "color(ebblue*0.3)"

// Load data
use "$DatDir/MEPS_HRS_1992_2017_1r.dta", clear
do $DoDir/f_globalLists.do
//}}}1
// [2] Set HH or Individual level weights {{{1
// ----------------------------------------------------------------------------
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
// [3] Programs to save and load estimation results to/from a file {{{1
// ----------------------------------------------------------------------------
// [] Define: f_saveEstToFile {{{2
// ----------------------------------------------------------------------------
capture program drop f_saveEstToFile
program define f_saveEstToFile
    // Program save estimation results and e(sample) info a subfolder /StoredEstimates
    // under the name `in_name'

    // Input 'variables'
    // -------------------------------------------------
    //`1' - name of file
    //
    //`2' - ....
    //`3' - .....
    // -------------------------------------------------

    local in_name    = "`1'"
    //local in_...  = "`2'"
    //local in_...  = "`3'"

    // Store Results in File
    preserve
        estimates save "$DoDir/StoredEstimates/`in_name'", replace // store the results of ologit into a .ster file
        keep if e(sample)
        compress
        save "$DoDir/StoredEstimates/`in_name'_esample.dta", replace // store e(sample) with it
    restore
end
//}}}2
// [] Define: f_loadEstFromFile {{{2
// ----------------------------------------------------------------------------
capture program drop f_loadEstFromFile
program define f_loadEstFromFile
    // Program loads estimation results and e(sample) from a subfolder /StoredEstimates
    // under the name `in_name'

    // Input 'variables'
    // -------------------------------------------------
    //`1' - name of file
    //
    //`2' - ....
    //`3' - .....
    // -------------------------------------------------

    local in_name    = "`1'"
    //local in_...  = "`2'"
    //local in_...  = "`3'"

    // Now stored e(sample) .dta file
    use "$DoDir/StoredEstimates/`in_name'_esample.dta", clear
    capture drop temp_sample
    gen temp_sample = 1
    // Load the results of ologit from a .ster file
    estimates use "$DoDir/StoredEstimates/`in_name'"
    // Set e(sample) so post estimation command work as expected
    estimates esample: if temp_sample == 1
end
//}}}2
//}}}1
// [4] Run regression {{{1
// --------------------------------------------------------------------------------
// HRS
// The ologit is run on a lag health state in the regression - so it only allows for a different intercept
ologit F_healthStatusD i.healthStatus i.initHealth $Xv $periodInteraction  d_cohort* if d_MEPS==0 [pweight = `weight'], cluster(persID) difficult iterate(1000)
estimates store res_ologit_initHealth_HRS // store ologit results
f_saveEstToFile "ologit_initHealth_Table_HRS"

// MEPS
// NOTE: In MEPS we cannot control for inital health because T<3
ologit F_healthStatusD $Xv $periodInteraction i.healthStatus d_cohort* if d_MEPS==1 [pweight = `weight'], cluster(persID) difficult iterate(1000)
estimates store res_ologit_initHealth_MEPS // store ologit results
f_saveEstToFile "ologit_Table_MEPS"
//}}}1
// [4.A] Incl. Death: Calculate  ologit predictions - reg F_healthStatus healthStatus Xv {{{1
// ----------------------------------------------------------------------------
// All genders between 20 and 95
// [1] Excellent Health:
// [2] Very Good Health:
// [3] Good Health:
// [4] Fair Health:
// [5] Poor Health:
// [6] Dead:

//foreach data_name in HRS MEPS {
foreach data_name in HRS {

    forvalues i=1/5 {

        display "------------------------"
        display `i'
        display "------------------------"

        if "`data_name'" == "HRS" {
            di "Load HRS estimation results"
            f_loadEstFromFile "ologit_initHealth_Table_HRS"
            local time_lag = 2
        }
        else {
            di "Load MEPS estimation results"
            f_loadEstFromFile "ologit_Table_MEPS"
            local time_lag = 1
        }

        di "Calculating mtable ... "
        mtable, at(healthStatus = `i') statistics(est se) // marginal effects averaged across all obs

        mat Pt=r(table)
        mat Pa_pr = Pt[1,1..6]*100 // to express in percent
        mat Pa_se = Pt[2,1..6]*100

        if `i' == 1 {
            mat Pa=Pa_pr // start new matrix
            mat Pa=(Pa \ Pa_se) // add standard errors
        }
        else {
            mat Pa= (Pa \ Pa_pr) // add next line probabilities
            mat Pa= (Pa \ Pa_se) // add standard errors
        }
    }

    mat jvoid2= (.,.,.,.,.,.,.)
    mat j100 = (100\.\100\.\100\.\100\.\100\.)
    mat temp = (Pa, j100) // tag on 100 at the end to indicate probs sum to 100%
    mat jout = (jvoid2 \ temp)

    mat rownames jout = Health_(t) Excellent (s.e.) Very_Good (s.e.) Good  (s.e.) Fair (s.e.) Poor (s.e.)
    mat colnames jout = Excellent Very_Good Good Fair Poor Dead Sum
    mat list jout

    if "`data_name'" == "HRS" {
        // HRS Output Option
        estout matrix(jout, fmt("1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            )) ///
            using "$TexDir/oLogitTMatrix_`data_name'_Death_initHealth.tex", replace style(tex) ///
            prehead(\begin{tabular}{l*{8}{r}} \toprule)  posthead(\midrule) ///
            prefoot(\bottomrule) postfoot(\end{tabular}) substitute(_ \,)  ///
            mlabels(Health_(t+`time_lag'), span prefix(\multicolumn{7}{c}{) suffix(}))
    }
    else {
        // MEPS Output Option
        estout matrix(jout, fmt("1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            "1 1 3 1 3 1 3 1 3 1 3" ///
                            )) ///
            using "$TexDir/oLogitTMatrix_`data_name'_Death.tex", replace style(tex) ///
            prehead(\begin{tabular}{l*{8}{r}} \toprule)  posthead(\midrule) ///
            prefoot(\bottomrule) postfoot(\end{tabular}) substitute(_ \,)  ///
            mlabels(Health_(t+`time_lag'), span prefix(\multicolumn{7}{c}{) suffix(}))
    }

}
//}}}1
// [4.B] Incl. Death: Calculate Marginal Effects of ologit {{{1
// ----------------------------------------------------------------------------
// based on ologit F_healthStatus healthStatus Xv
//
// Change long label in table for pre-tax HH income
label variable preGovIncHH "HH Gross Income (in \\$1,000)"
label variable age "Age"

foreach data_name in HRS MEPS {

    local i = 1
    foreach o in Excellent Very_Good Good Fair Poor Dead {

        if "`data_name'" == "HRS" {
            di "Load HRS estimation results"
            f_loadEstFromFile "ologit_initHealth_Table_HRS"

            di "Calculating margins ... "
            margins, dydx(healthStatus initHealth age d_period1 d_period2 d_female d_partner d_smoke d_black d_hispanic d_hiSchool d_college preGovIncHH) predict(outcome(`i')) post
        }
        else {
            di "Load MEPS estimation results"
            f_loadEstFromFile "ologit_Table_MEPS"

            di "Calculating margins ... "
            margins, dydx(healthStatus age d_period1 d_period2 d_female d_partner d_smoke d_black d_hispanic d_hiSchool d_college preGovIncHH) predict(outcome(`i')) post
        }

        eststo, title(`o')
        local i = `i' + 1
    }

    if "`data_name'" == "HRS" {
        esttab using "$TexDir/oLogit_TableTransProbMargEff_`data_name'_PeriodInteraction_InitHealth.tex", replace noobs b(3) se(3) mtitles nonumbers label booktabs nonote
    }
    else {
        esttab using "$TexDir/oLogit_TableTransProbMargEff_`data_name'_PeriodInteraction.tex", replace noobs b(3) se(3) mtitles nonumbers label booktabs nonote
    }

    eststo clear
}

// (Does not work - cannot get results in esttab)
// Version B. Marginal effects using mtable
// ----------------------------------------
//estimate restore res_ologit
//mtable, dydx(L_healthStatus age d_female) statistics(est se)

//eststo res_mtable_dydx
//esttab res_mtable_dydx using "$TexDir/jTableTransProbMargEff_mtable.tex", replace
//}}}1
log close
