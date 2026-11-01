*==============================================================================
*  01_analysis.do
*  Replication package
*  "The Dual-Pilot Policy and Urban Green Total Factor Productivity:
*   Evidence from Chinese Prefecture-Level Cities"
*
*  Purpose : reproduce Tables 3-10 of the manuscript
*  Input   : data_EN.dta  (built by 00_data_prep.do)
*  Output  : output/Table3_Baseline.rtf      ... output/Table10_Synergy.rtf
*
*  Software: Stata 16 or later; ddml + pystacked (R package of the same name
*            is not required, only a working Python installation)
*
*  Estimator: partially linear double machine learning
*             - 3-fold cross-fitting, 1:2 sample split (Table 6 relaxes this)
*             - random forest for both nuisance functions (Table 6 relaxes this)
*             - seed 42, robust standard errors
*
*  Required packages (run once):
*      ssc install ddml, replace
*      ssc install pystacked, replace
*      ssc install reghdfe, replace
*      ssc install ftools, replace
*      ssc install winsor2, replace
*      ssc install estout, replace
*==============================================================================

clear all
set more off
set varabbrev off
capture log close
capture mkdir output

*------------------------------------------------------------------------------
* 0. Globals
*------------------------------------------------------------------------------
global Y      GTFP
global D      Low_Inn
global CTRL   EAL FDI FDO UR HC IL
global CTRL2  EAL_sq FDI_sq FDO_sq UR_sq HC_sq IL_sq
global FE     i.year i.city_id
global X      $CTRL $CTRL2 $FE
global SEED   42

* esttab formatting shared by all tables
global TABOPT "b(%9.4f) se(%9.4f) star(* 0.10 ** 0.05 *** 0.01) stats(N, labels("Observations") fmt(%9.0f)) nocons compress nogaps"

*------------------------------------------------------------------------------
* 1. Wrapper programs around ddml
*
*    Wrapping the five estimation steps in a program guarantees that every
*    specification in the paper uses an identical estimator. Only the
*    arguments change across tables.
*------------------------------------------------------------------------------

*--- 1.1 Partially linear DML -------------------------------------------------
capture program drop dml_plr
program define dml_plr
    syntax varlist(min=2 max=2), X(string) [K(integer 3) LEARN(string) SEED(integer 42)]
    tokenize "`varlist'"
    local y "`1'"
    local d "`2'"
    if "`learn'" == "" local learn "rf"
    set seed `seed'
    quietly ddml init partial, kfolds(`k')
    quietly ddml E[D|X]: pystacked `d' `x', type(reg) method(`learn')
    quietly ddml E[Y|X]: pystacked `y' `x', type(reg) method(`learn')
    quietly ddml crossfit
    ddml estimate, robust
end

*--- 1.2 Interactive DML ------------------------------------------------------
capture program drop dml_int
program define dml_int
    syntax varlist(min=2 max=2), X(string) [K(integer 3) LEARN(string) SEED(integer 42)]
    tokenize "`varlist'"
    local y "`1'"
    local d "`2'"
    if "`learn'" == "" local learn "rf"
    set seed `seed'
    quietly ddml init interactive, kfolds(`k')
    quietly ddml E[Y|X,D]: pystacked `y' `x', type(reg) method(`learn')
    quietly ddml E[D|X]:   pystacked `d' `x', type(reg) method(`learn')
    quietly ddml crossfit
    ddml estimate, robust
end

*--- 1.3 IV-DML ---------------------------------------------------------------
capture program drop dml_iv
program define dml_iv
    syntax varlist(min=3 max=3), X(string) [K(integer 3) LEARN(string) SEED(integer 42)]
    tokenize "`varlist'"
    local y "`1'"
    local d "`2'"
    local z "`3'"
    if "`learn'" == "" local learn "rf"
    set seed `seed'
    quietly ddml init iv, kfolds(`k')
    quietly ddml E[Y|X]: pystacked `y' `x', type(reg) method(`learn')
    quietly ddml E[D|X]: pystacked `d' `x', type(reg) method(`learn')
    quietly ddml E[Z|X]: pystacked `z' `x', type(reg) method(`learn')
    quietly ddml crossfit
    ddml estimate, robust
end

*------------------------------------------------------------------------------
* 8. Table 2: descriptive statistics of the main variables
*------------------------------------------------------------------------------
use "data_EN.dta", clear

estpost summarize GTFP Low_Inn EAL IL FDI FDO UR HC, detail
esttab using "output/Table2_Descriptive.rtf", replace ///
    cells("count(fmt(%9.0f)) mean(fmt(%9.3f)) sd(fmt(%9.3f)) min(fmt(%9.3f)) max(fmt(%9.3f))") ///
    noobs nonumber nomtitle ///
    collabels("Obs" "Mean" "SD" "Min" "Max") ///
    title("Table 2. Descriptive statistics of main variables")

di as text _n ">>> 00_data_prep.do finished: data_EN.dta and Table2_Descriptive.rtf written"

*==============================================================================
*                        TABLE 3  BASELINE REGRESSION
*  (1) city and year fixed effects only
*  (2) linear controls added
*  (3) linear and quadratic controls added  [baseline specification]
*==============================================================================
di as text _n(2) "=== Table 3: baseline regression ==="
use "data_EN.dta", clear

dml_plr $Y $D, x(i.year i.city_id)
estimates store t3_1

dml_plr $Y $D, x($CTRL $FE)
estimates store t3_2

dml_plr $Y $D, x($X)
estimates store t3_3

esttab t3_1 t3_2 t3_3 using "output/Table3_Baseline.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("(1)" "(2)" "(3)") ///
    title("Table 3. Baseline regression results") ///
    addnotes("Notes: partially linear DML with 3-fold cross-fitting and random forest; seed 42. Robust standard errors in parentheses.")

*==============================================================================
*                     TABLE 4  ROBUSTNESS: SAMPLE AND SPECIFICATION
*==============================================================================
di as text _n(2) "=== Table 4: robustness ==="

*--- (1) exclude autonomous regions and directly administered municipalities
use "data_EN.dta", clear
keep if city_adj == 1
dml_plr $Y $D, x($X)
estimates store t4_1

*--- (2) drop 2015-2016 (supply-side structural reform)
use "data_EN.dta", clear
drop if inrange(year, 2015, 2016)
dml_plr $Y $D, x($X)
estimates store t4_2

*--- (3) winsorize at the 1st and 99th percentiles
use "data_EN.dta", clear
winsor2 $Y $CTRL $CTRL2, cuts(1 99) replace
dml_plr $Y $D, x($X)
estimates store t4_3

*--- (4) winsorize at the 10th and 90th percentiles
use "data_EN.dta", clear
winsor2 $Y $CTRL $CTRL2, cuts(10 90) replace
dml_plr $Y $D, x($X)
estimates store t4_4

*--- (5) province-by-year fixed effects
use "data_EN.dta", clear
encode province, gen(prov_id)
dml_plr $Y $D, x($CTRL $CTRL2 i.year i.city_id i.year#i.prov_id)
estimates store t4_5

esttab t4_1 t4_2 t4_3 t4_4 t4_5 using "output/Table4_Robustness.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("(1)" "(2)" "(3)" "(4)" "(5)") ///
    title("Table 4. Robustness test results") ///
    addnotes("Notes: column (1) excludes autonomous regions and municipalities; column (2) drops 2015-2016; columns (3)-(4) winsorize at the 1/99 and 10/90 percentiles; column (5) adds province-by-year fixed effects.")

*==============================================================================
*                TABLE 5  ROBUSTNESS: CONCURRENT POLICIES
*  Each column adds a different concurrent policy to the covariate set;
*  column (5) controls for all four at once.
*==============================================================================
di as text _n(2) "=== Table 5: concurrent policies ==="
use "data_EN.dta", clear

dml_plr $Y $D, x($X SMA)
estimates store t5_1

dml_plr $Y $D, x($X BIGDATA)
estimates store t5_2

dml_plr $Y $D, x($X TRA)
estimates store t5_3

dml_plr $Y $D, x($X SPE)
estimates store t5_4

dml_plr $Y $D, x($X SMA BIGDATA TRA SPE)
estimates store t5_5

esttab t5_1 t5_2 t5_3 t5_4 t5_5 using "output/Table5_Concurrent.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("(1)" "(2)" "(3)" "(4)" "(5)") ///
    title("Table 5. Robustness test excluding the interference of relevant policies") ///
    addnotes("Notes: SMA = smart city pilot; BIGDATA = national big data pilot zone; TRA = carbon emission trading pilot; SPE = special air pollutant emission limit policy.")

*==============================================================================
*          TABLE 6  ROBUSTNESS: RE-ESTIMATED DML SPECIFICATIONS
*  (1) 1:4 split (5 folds)     (2) 1:7 split (8 folds)
*  (3) GBDT learner            (4) SVM learner        (5) interactive DML
*==============================================================================
di as text _n(2) "=== Table 6: re-estimated DML ==="
use "data_EN.dta", clear

dml_plr $Y $D, x($X) k(5)
estimates store t6_1

dml_plr $Y $D, x($X) k(8)
estimates store t6_2

dml_plr $Y $D, x($X) learn(gradboost)
estimates store t6_3

dml_plr $Y $D, x($X) learn(svm)
estimates store t6_4

dml_int $Y $D, x($X)
estimates store t6_5

esttab t6_1 t6_2 t6_3 t6_4 t6_5 using "output/Table6_Reestimated.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("1:4" "1:7" "GBDT" "SVM" "Interactive") ///
    title("Table 6. Robustness checks for the re-estimated double machine learning model") ///
    addnotes("Notes: columns (1)-(2) change the cross-fitting split; columns (3)-(4) change the learner; column (5) replaces the partially linear model with an interactive model.")

*==============================================================================
*                    TABLE 7  ENDOGENEITY: IV-DML
*  The instruments are interacted with (year - 2008) so that they vary over
*  time, as required for a panel setting.
*==============================================================================
di as text _n(2) "=== Table 7: IV-DML ==="
use "data_EN.dta", clear
generate double z_rugged = IV_RUGGED * (year - 2008)
generate double z_dist   = IV_DIST   * (year - 2008)
label variable z_rugged "Instrument: terrain ruggedness x (year - 2008)"
label variable z_dist   "Instrument: distance to Hangzhou x (year - 2008)"

*--- (1) terrain ruggedness
dml_iv $Y $D z_rugged, x($X)
estimates store t7_1
capture ivreg2 Y1_pystacked_1 (D1_pystacked_1 = Z1_pystacked_1), nocons robust first

*--- (2) distance to Hangzhou
dml_iv $Y $D z_dist, x($X)
estimates store t7_2
capture ivreg2 Y1_pystacked_1 (D1_pystacked_1 = Z1_pystacked_1), nocons robust first

esttab t7_1 t7_2 using "output/Table7_IV.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("(1)" "(2)") ///
    title("Table 7. Endogeneity test results") ///
    addnotes("Notes: IV-DML with 3-fold cross-fitting and random forest. Kleibergen-Paap rk LM and Wald F statistics are reported after ddml estimate. The Stock-Yogo (2005) critical value for a 10% maximal IV size with one endogenous regressor is 16.38.")

*==============================================================================
*                      TABLE 8  MECHANISM ANALYSIS
*  (1)-(2) green innovation bias      (3)-(5) green technology diffusion
*  (6)-(7) green factor agglomeration
*==============================================================================
di as text _n(2) "=== Table 8: mechanisms ==="
use "data_EN.dta", clear

foreach m in GIB1 GIB2 GTD1 GTD2 GTD3 GFA1 GFA2 {
    dml_plr `m' $D, x($X)
    estimates store mech_`m'
}

esttab mech_GIB1 mech_GIB2 mech_GTD1 mech_GTD2 mech_GTD3 mech_GFA1 mech_GFA2 ///
    using "output/Table8_Mechanism.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("GIB1" "GIB2" "GTD1" "GTD2" "GTD3" "GFA1" "GFA2") ///
    title("Table 8. Results of mechanism tests") ///
    addnotes("Notes: GIB = green innovation bias; GTD = green technology diffusion; GFA = green factor agglomeration.")

*==============================================================================
*                     TABLE 9  HETEROGENEITY ANALYSIS
*  Splits are pre-computed median splits stored in data_EN.dta;
*  1 always denotes the "high" group.
*==============================================================================
di as text _n(2) "=== Table 9: heterogeneity ==="

use "data_EN.dta", clear
keep if INNO_FDN == 1
dml_plr $Y $D, x($X)
estimates store t9_1

use "data_EN.dta", clear
keep if INNO_FDN == 0
dml_plr $Y $D, x($X)
estimates store t9_2

use "data_EN.dta", clear
keep if ER_high == 1
dml_plr $Y $D, x($X)
estimates store t9_3

use "data_EN.dta", clear
keep if ER_high == 0
dml_plr $Y $D, x($X)
estimates store t9_4

use "data_EN.dta", clear
keep if MAR == 1
dml_plr $Y $D, x($X)
estimates store t9_5

use "data_EN.dta", clear
keep if MAR == 0
dml_plr $Y $D, x($X)
estimates store t9_6

esttab t9_1 t9_2 t9_3 t9_4 t9_5 t9_6 using "output/Table9_Heterogeneity.rtf", replace $TABOPT ///
    keep(Low_Inn) mtitles("High IF" "Low IF" "High ER" "Low ER" "High MAR" "Low MAR") ///
    title("Table 9. Results of heterogeneity test") ///
    addnotes("Notes: IF = innovation foundation; ER = environmental regulation intensity; MAR = marketization level.")

*==============================================================================
*                TABLE 10  SYNERGY BETWEEN THE TWO PILOTS
*
*  (1) net effect of the low-carbon pilot, sampled on cities that are never
*      designated as innovative cities
*  (2) net effect of the innovative city pilot, sampled on cities that are
*      never designated as low-carbon cities
*  (3) dual-pilot effect within the restricted sample of cities that
*      eventually hold both designations
*  (4) low-carbon cities, treatment = subsequent innovative city designation
*  (5) innovative cities, treatment = subsequent low-carbon designation
*==============================================================================
di as text _n(2) "=== Table 10: synergy ==="
use "data_EN.dta", clear
bysort city_id: egen ever_LCCP = max(LCCP)
bysort city_id: egen ever_ICP  = max(ICP)
label variable ever_LCCP "City ever designated as a low-carbon pilot"
label variable ever_ICP  "City ever designated as an innovative city pilot"

*--- (1) low-carbon pilot net effect
preserve
keep if ever_ICP == 0
dml_plr $Y LCCP, x($X)
estimates store t10_1
restore

*--- (2) innovative city pilot net effect
preserve
keep if ever_LCCP == 0
dml_plr $Y ICP, x($X)
estimates store t10_2
restore

*--- (3) dual-pilot effect in the restricted sample
preserve
keep if ever_LCCP == 1 & ever_ICP == 1
dml_plr $Y $D, x($X)
estimates store t10_3
restore

*--- base-period comparability before the path comparison
*    (2009 is prior to every designation, so the two groups should be
*     statistically indistinguishable in the outcome at that point)
preserve
keep if year == 2009
generate byte treat_up = (ever_LCCP == 1 & ever_ICP == 1)
ttest GTFP, by(treat_up) unequal
restore

*--- (4) low-carbon cities upgrading to the innovative city pilot
preserve
keep if ever_LCCP == 1
generate byte treat_icp = (ever_ICP == 1)
ttest GTFP if year == 2009, by(treat_icp) unequal
dml_plr $Y ICP, x($X)
estimates store t10_4
restore

*--- (5) innovative cities upgrading to the low-carbon pilot
preserve
keep if ever_ICP == 1
generate byte treat_lccp = (ever_LCCP == 1)
ttest GTFP if year == 2009, by(treat_lccp) unequal
dml_plr $Y LCCP, x($X)
estimates store t10_5
restore

esttab t10_1 t10_2 t10_3 t10_4 t10_5 using "output/Table10_Synergy.rtf", replace $TABOPT ///
    keep(Low_Inn LCCP ICP) ///
    mtitles("Low" "Inn" "Low_Inn" "Low-carbon first" "Innovation first") ///
    title("Table 10. Synergy effect analysis results") ///
    addnotes("Notes: columns (1)-(3) report the coefficient on the treatment variable named in the column header; columns (4)-(5) report the effect of acquiring the second designation. Coefficients in columns (1)-(3) are estimated on different samples and their comparison is descriptive rather than a formal test of equality.")

di as text _n(2) ">>> 01_analysis.do finished: Tables 3-10 written to output/"
