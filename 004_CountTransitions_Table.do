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
// [3.A] Unweighted: Make count tables {{{1
// ----------------------------------------------------------------------------
// MEPS
// ----------------------------------------------------------------------------
label variable L_healthStatus "Health Status (t-1)"
label variable healthStatus "Health Status (t)"
label variable F_healthStatus "Health Status (t+1)"

label variable L_healthStatusD "Health Status (t-1)"
label variable healthStatusD "Health Status (t)"
label variable F_healthStatusD "Health Status (t+1) incl. Death"

tabout  healthStatus F_healthStatus if d_MEPS==1  ///
   using "$TexDir/TableTransProb1_MEPS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(20)

// Incl. Death
// NOTE: In MEPS healthStatus = poor & healthStatusD=6 is possible if person died in this year!!!
tabout  healthStatus F_healthStatusD if d_MEPS==1  ///
   using "$TexDir/Table_Count_TransProb_MEPS_Dead.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(20)

// HRS
// ----------------------------------------------------------------------------
label variable L_healthStatus "Health Status (t-2)"
label variable healthStatus "Health Status (t)"
label variable F_healthStatus "Health Status (t+2)"

label variable L_healthStatusD "Health Status (t-2)"
label variable healthStatusD "Health Status (t)"
label variable F_healthStatusD "Health Status (t+2) incl. Death"

tabout  healthStatus F_healthStatus if d_MEPS==0  ///
   using "$TexDir/Table_Count_TransProb_HRS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(20)

// Incl. Death
// NOTE: In HRS healthStatus = poor & healthStatusD=6 is NOT possible since person ALWAYS dies in next year with zero other observations
tabout  healthStatusD F_healthStatusD if d_MEPS==0 ///
   using "$TexDir/Table_Count_TransProb_HRS_Dead.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(20)

//}}}1
// [3.B] Weighted: Make count tables {{{1
// ----------------------------------------------------------------------------
// MEPS
// ----------------------------------------------------------------------------
label variable L_healthStatus "Health Status (t-1)"
label variable healthStatus "Health Status (t)"
label variable F_healthStatus "Health Status (t+1)"

label variable L_healthStatusD "Health Status (t-1)"
label variable healthStatusD "Health Status (t)"
label variable F_healthStatusD "Health Status (t+1) incl. Death"

xtset, clear
svyset [pweight=`weight'], strata(varstr) psu(varpsu) //vce(bootstrap) // vce(linearized) is not supported by tabout

// Unweighted
tab healthStatus F_healthStatus if d_MEPS == 1, row
// Weigthed
tab healthStatus F_healthStatus  [iweight = `weight'] if d_MEPS == 1, row
svy, sub(d_MEPS): tab healthStatus F_healthStatus, row

// This tabout command results in the same probabilities based on the weighted frequency count
tabout  healthStatus F_healthStatus [iweight = `weight'] if d_MEPS == 1 ///
   using "$TexDir/Table_Count_TransProb_MEPS_Weighted.tex", replace ///
   c(row) f(1p) style(tex) font(bold) twidth(20)

// Note: Does not work. Issues with vce(linearized)
//svy, sub(d_MEPS): tabout  healthStatus F_healthStatus  ///
//   using "$TexDir/Table_Count_TransProb_MEPS_Weighted.tex", replace ///
//   svy c(freq row) f(0c 1p) style(tex) font(bold) twidth(20)

// Unweighted
tab healthStatus F_healthStatusD if d_MEPS == 1, row
// Weigthed
svy, sub(d_MEPS): tab healthStatus F_healthStatusD, row
// This tabout command results in the same probabilities based on the weighted frequency count
// Incl. Death
// NOTE: In MEPS healthStatus = poor & healthStatusD=6 is possible if person died in this year!!!
tabout  healthStatus F_healthStatusD [iweight = `weight'] if d_MEPS==1  ///
   using "$TexDir/Table_Count_TransProb_MEPS_Dead_Weighted.tex", replace ///
   c(row) f(1p) style(tex) font(bold) twidth(20)

// HRS
// ----------------------------------------------------------------------------
label variable L_healthStatus "Health Status (t-2)"
label variable healthStatus "Health Status (t)"
label variable F_healthStatus "Health Status (t+2)"

label variable L_healthStatusD "Health Status (t-2)"
label variable healthStatusD "Health Status (t)"
label variable F_healthStatusD "Health Status (t+2) incl. Death"

tabout  healthStatus F_healthStatus [iweight = `weight'] if d_MEPS==0  ///
   using "$TexDir/Table_Count_TransProb_HRS_Weighted.tex", replace ///
   c(row) f(1p) style(tex) font(bold) twidth(20)

// Incl. Death
// NOTE: In HRS healthStatus = poor & healthStatusD=6 is NOT possible since person ALWAYS dies in next year with zero other observations
tabout  healthStatusD F_healthStatusD [iweight = `weight'] if d_MEPS==0 ///
   using "$TexDir/Table_Count_TransProb_HRS_Dead_Weighted.tex", replace ///
   c(row) f(1p) style(tex) font(bold) twidth(20)

//}}}1
// [4] Make count tables by gender {{{1
// ----------------------------------------------------------------------------
// MEPS
// ----------------------------------------------------------------------------
label variable healthStatus "Health Status (t)"
label variable F_healthStatusD "Health Status (t+1) incl. Death"

// FEMALE
tabout healthStatus F_healthStatusD if d_MEPS==1 & age>=`agel' & age<=`ageh' & d_female==1 ///
   using "$TexDir/Table_Count_TransProbFem_MEPS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(18)

// MALE
tabout healthStatus F_healthStatusD if d_MEPS==1 & age>=`agel' & age<=`ageh' & d_female==0 ///
   using "$TexDir/Table_Count_TransProbMale_MEPS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(18)

// HRS
// ----------------------------------------------------------------------------
label variable healthStatus "Health Status (t)"
label variable F_healthStatusD "Health Status (t+2) incl. Death"

// FEMALE
tabout healthStatus F_healthStatusD if d_MEPS==0 & age>=`agel' & age<=`ageh' & d_female==1 ///
   using "$TexDir/Table_Count_TransProbFem_HRS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(18)

// MALE
tabout healthStatus F_healthStatusD if d_MEPS==0 & age>=`agel' & age<=`ageh' & d_female==0 ///
   using "$TexDir/Table_Count_TransProbMale_HRS.tex", replace ///
   c(freq row) f(0c 1p) style(tex) font(bold) twidth(18)
//}}}1
