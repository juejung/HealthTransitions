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

local ageLo = 20
local ageHi = 95

// Graph options: MEPS vs HRS
global line_meps "msymbol(s) connect(l) mcolor(orange_red*1.75) lcolor(orange_red*1.75) lpattern(dash)"
global CI_meps "color(orange_red*0.3)"

global line_hrs "msymbol(o) connect(l) mcolor(gray*1.1) lcolor(gray*1.1) lpattern(dot)"
global CI_hrs "color(gray*0.2)"
global line_hrs_adj "msymbol(s) connect(l) mcolor(ebblue*1.4) lcolor(ebblue*1.4) lpattern(dot)"


// Graph options: Males vs. Females
global line_f "msymbol(s) connect(l) mcolor(orange_red*1.5) lcolor(orange_red*1.5) lpattern(dash)"
global CI_f "color(orange_red*0.3)"
// HRS adjusted line
global line_f_adj "msymbol(s) connect(l) mcolor(orange_red*1.85) lcolor(orange_red*1.85) lpattern(dash)"

global line_m "msymbol(o) connect(l) mcolor(ebblue*1.5) lcolor(ebblue*1.5) lpattern(dot)"
global CI_m "color(ebblue*0.3)"
// HRS adjusted line
global line_m_adj "msymbol(o) connect(l) mcolor(ebblue*1.85) lcolor(ebblue*1.85) lpattern(dot)"

// Load data
use "$DatDir/MEPS_HRS_1992_2017_1r.dta", clear
//}}}1
// [2] Count by 5 year age group and store results in matrices {{{1
// ----------------------------------------------------------------------------
label variable L_healthStatus "Health Status (t-1)"
label variable healthStatus "Health Status (t)"
label variable F_healthStatus "Health Status (t+1)"

label variable L_healthStatusD "Health Status (t-1)"
label variable healthStatusD "Health Status (t)"
label variable F_healthStatusD "Health Status (t+1) incl. Death"

// MEPS and HRS
forvalues i_meps = 0/1 {
    forvalues i_age = 1/15 {
        mat A = J(5,6,.)
        mat B`i_age'`i_meps' = J(5,6,.)

        if `i_age' == 1 {
            local age_low = 20
            local age_hi  = 24
        }
        else {
            local age_low = 15 + `i_age' * 5 - 2
            local age_hi  = 15 + `i_age' * 5 + 2
        }

        //tabulate  healthStatus F_healthStatusD if d_MEPS==`i_meps' & lage == `i_age', row nofreq matcell(A)
        tabulate  healthStatus F_healthStatusD if d_MEPS==`i_meps' & age >= `age_low' & age <= `age_hi', row nofreq matcell(A)


        mata : st_matrix("C", colsum(st_matrix("A")))

        forvalues i = 1/5 {
            forvalues j = 1/6 {
                mat B`i_age'`i_meps'[`i',`j'] = A[`i',`j']/C[1, `i']
            }
        }
        matlist B`i_age'`i_meps'
    }
}

//}}}1
// [3] Convert matrices to data {{{1
// ----------------------------------------------------------------------------
// Move to data and plot
// svmat by lage --> somehow
// make column vector meps_a1to1v(age) meps_a1to2v(age) .... meps_5to6v(agev) and then use svmat to store
// make column vector hrs_a1to1v(age) hrs_a1to2v(age) .... hrs_5to6v(agev) and then use svmat to store
//
clear // clears data

mat lagev = J(15,1,.)
forvalues i_age = 1/15 {
    mat lagev[`i_age',1] = 20 + (`i_age' - 1) * 5
}
svmat lagev

forvalues i_meps = 0/1 {
    forvalues i_row = 1/5 {
        forvalues i_col = 1/6 {

            mat p_`i_meps'_`i_row'_`i_col'v = J(15,1,.)

            forvalues i_age = 1/15 {

                mat p_`i_meps'_`i_row'_`i_col'v[`i_age',1] = B`i_age'`i_meps'[`i_row',`i_col']
            }
            svmat p_`i_meps'_`i_row'_`i_col'v
        }
    }
}

save "$DatDir/MEPS_HRS_1992-2017_TransitionCounts_Age.dta", replace
//}}}1
// [4] Make Plots - Data Cleaning Program {{{1
//-----------------------------------------------------------------------------
capture program drop f_prepareData
program define f_prepareData

// Cleaning steps: Express all variables as percentages, add labels, etc.

    qui foreach v of varlist _all {
        replace `v' = `v' * 100
        label variable `v' %
    }
    replace lagev = lagev/100
    label variable lagev "Age"

end
//}}}1
// [5] Plot counted relative frequencies by age group {{{1
//-----------------------------------------------------------------------------
use "$DatDir/MEPS_HRS_1992-2017_TransitionCounts_Age.dta", clear

// Call our cleaning step
f_prepareData

// Start plot
local h_from = 1
foreach from_status in "Excellent" "Very Good" "Good" "Fair" "Poor" {

    local h_to   = 1
    foreach to_status in "Excellent" "Very Good" "Good" "Fair" "Poor" "Dead" {

        if "`to_status'" != "Dead" {
            graph twoway    (connected p_1_`h_from'_`h_to'v1 lagev1, $line_meps)  ///
                        (connected p_0_`h_from'_`h_to'v1 lagev1, $line_hrs), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(off) name(fig`h_to', replace)
        }
        else {
            graph twoway    (connected p_1_`h_from'_`h_to'v1 lagev1, $line_meps)  ///
                        (connected p_0_`h_from'_`h_to'v1 lagev1, $line_hrs), ///
                    xlabel(20(10)95) title("P(`to_status' | `from_status')") ///
                    legend(ring(0) position(10) col(1) order(1 "MEPS" 2 "HRS")) ///
                    name(fig`h_to', replace)
        }
        local h_to = `h_to' + 1
    }

    // Make an empty plot, just to get legend out!!
    //graph twoway (connected p_1_1_1v1 lagev1, $line_meps)  ///
    //             (connected p_0_1_1v1 lagev1, $line_hrs), ///
    //        yscale(off) xscale(off)  ///
    //        legend(ring(0) position(0) col(2) order(1 "MEPS" 2 "HRS")) ///
    //        name(jleg, replace)

    //gr_edit .plotregion1.draw_view.setstyle, style(no)

    // -----------------------------------------------
    //graph combine fig1 fig2 fig3 fig4 fig5 fig6 jleg, ///
    graph combine fig1 fig2 fig3 fig4 fig5 fig6, ///
        xcommon col(2) imargin(1 1 1 1) graphregion(margin(l=1 r=1)) ///
        iscale(0.8) ysize(20.0) xsize(16.0) ///
        name(ffig,replace)
    //
    graph export "$TexDir/MEPSvsHRS_Count_TransGraph`h_from'_Dead.svg", replace
    //
    local h_from = `h_from' + 1
}
//}}}1
