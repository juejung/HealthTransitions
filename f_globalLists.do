// ----------------------------------------------------------------------------
//                      f_globalList.do file
// ----------------------------------------------------------------------------
global MEPS_HRS_Keep_Listv ///
    year age lage d_female d_partner d_black d_hispanic numChildren region ///
    laborIncome laborIncHH preGovIncHH ///
    healthStatus L_healthStatus F_healthStatus   ///
    healthStatusD L_healthStatusD F_healthStatusD ///
    d_healthy BMI d_smoke d_dead ///
    F_d_attrit F_d_attritDeath F_d_attritNoDeath ///
    healthExpenditure healthExp_noCohort ///
    healthExpenditureHH ///
    OOPExpenditure OOPExpenditureHH OOPExp_noCohort ///
    d_noHiSchool d_hiSchool d_college ///
    d_insured d_public_insurance d_private_insurance d_head ///
    d_cohort* ///
    persID persweight int_persweight ///
    famID famweight int_famweight ///
    varstr varpsu ///
    d_NMES d_MEPS d_HRS ///
    l5cohort l10cohort period

global MEPS_Listv  healthCapital

// ----------------------------------------------------------------------------
//           Set global lists for summary stats and regressions
// ----------------------------------------------------------------------------
global MEPS_HRS_SumStat1   d_year*
global MEPS_HRS_SumStat2   age_year*

global MEPS_HRS_SumStat3   ///
    age d_female d_partner d_black ///
    d_noHiSchool d_hiSchool d_college ///
    laborIncome laborIncHH preGovIncHH ///
    d_Excellent d_VeryGood d_Good d_Fair d_Poor ///
    d_iniExcellent d_iniVeryGood d_iniGood d_iniFair d_iniPoor ///
    d_dead F_d_attritDeath ///
    d_healthy BMI d_smoke ///
    OOPExpenditure OOPExpenditureHH ///
    d_insured d_public_insurance d_private_insurance d_head

global periodInteraction i.d_period1 i.d_period2 ///
       i.d_black#i.d_period1 i.d_black#i.d_period2 ///
       i.d_female#i.d_period1 i.d_female#i.d_period2

global Xv_ageOnly ///
    c.age c.age#c.age c.age#c.age#c.age


global Xv ///
    c.age c.age#c.age c.age#c.age#c.age ///
    i.d_partner ///
    i.d_smoke ///
    i.d_female i.d_female#c.age  ///
    i.d_black i.d_black#c.age ///
    i.d_hispanic ///
    i.d_hiSchool i.d_college c.preGovIncHH i.region

global Xv_cmp ///
    c.age c.age#c.age c.age#c.age#c.age ///
    i.d_partner ///
    i.d_female i.d_female#c.age  ///
    i.d_black i.d_black#c.age ///
    i.d_hispanic ///
    i.d_hiSchool i.d_college i.region


global Xv_higherPoly ///
    c.age c.age#c.age c.age#c.age#c.age ///
    i.d_partner ///
    i.d_smoke ///
    i.d_female i.d_female#c.age i.d_female#c.age#c.age  ///
    i.d_black  i.d_black#c.age  i.d_black#c.age#c.age   ///
    i.d_hispanic ///
    i.d_hiSchool i.d_college c.preGovIncHH i.region


global Xv_esample ///
    age ///
    d_female d_partner d_smoke d_black d_hispanic ///
    d_hiSchool d_college preGovIncHH region

