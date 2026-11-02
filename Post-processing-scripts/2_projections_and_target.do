/*==============================================================================
 2_projections_and_target.do

 Part of the code accompanying:
   Strain T, Flaxman SR, Ricardo LIC, et al, and the Country Data Author Group.
   National, regional, and global trends in adolescent physical activity from
   2010 to 2025: a pooled analysis of 389 surveys with 3.2 million participants
   from 140 countries. Lancet Glob Health 2026.
   Code repository: https://github.com/MLGlobalHealth/AdoPA-2026

 PURPOSE
   Projects the percentage of adolescents who report fewer than 7 active days per
   week to 2030, assuming that the 2010-2025 trend continues, and assesses whether
   each country, region and the globe (for boys, girls and both sexes) is on track
   to meet the global target. The target is a 15% relative reduction between 2010
   and 2030 in the percentage of adolescents who do not report 60 minutes of
   moderate-to-vigorous physical activity on every day of the week, that is, who
   report fewer than 7 active days per week. The 2010 estimate is the baseline.

 METHOD (carried out separately for each of the 4000 posterior draws, so that the
 uncertainty in the model carries through to the projection)
   1. Take the percentage reporting fewer than 7 active days (1 minus the
      prevalence of 7 active days) for each year 2010-2025 and convert it to a
      probit (the inverse of the standard normal distribution function).
   2. For each country, region (and the globe) and sex, fit a straight line to the
      probit against year (2010-2025).
   3. Project the line to 2030 and convert it back to a percentage. This is the
      projected percentage reporting fewer than 7 active days in 2030.
   4. The target for 2030 is 85% of the same draw's 2010 percentage reporting
      fewer than 7 active days (a 15% relative reduction).
   5. A draw meets the target if its projection for 2030 is below its target.
   The posterior probability of meeting the target (pp_target) is the proportion
   of draws that meet it. The projection is summarised as the mean of the draws and
   a 95% uncertainty interval (2.5th and 97.5th percentiles).

 INPUT (in the same folder as this do-file)
   country_globe_region_all_ages_draws.dta
        Draws of the prevalence of 7 active days per week for countries, regions
        and the globe, 2010-2025. It is written by 1_summarize_draws.do, which
        must be run first. This file is large and is not included in the
        repository.

 OUTPUT (written to the same folder)
   on_track_assessment.csv   One row per country, region or globe and sex
                             (sex: 1 boys, 2 girls, 3 both sexes), with:
     region, regiontype, iso3   Location (iso3 is blank for regions and the globe;
                                regiontype is country, adopa_region, whoreg6,
                                wbinc25 or Global). The rows for regiontype wbinc25
                                with a blank region are countries with no World
                                Bank income classification and are not published.
     lt7days2010 (_lwr, _upr)   Baseline: percentage in 2010 reporting fewer than 7
                                active days, as a proportion, with its uncertainty
                                interval.
     target2030                 The 2030 target: 0.85 x lt7days2010.
     pred2030 (_lwr, _upr)      Projected percentage in 2030, as a proportion.
     pp_target                  Posterior probability of meeting the target.
     pp_increase                Proportion of draws in which the percentage
                                reporting fewer than 7 active days is increasing
                                over 2010-2025 (positive slope).
     on_track                   1 if the projected percentage in 2030 (mean of the
                                draws) is below target2030, otherwise 0.
     on_track_certainty         "high" if pp_target is below 0.2 or above 0.8,
                                otherwise "low".
     prev_0days ... prev_6days  Prevalence of 0 to 6 active days in 2010 (with
                                _lwr and _upr), carried through from the input.
   All percentages are proportions between 0 and 1.

 REQUIREMENTS
   - Stata/MP. About 16,000 variables are held in memory at the peak (the 4000
     draws of the input, plus 4000 projections, 4000 slopes and 4000 indicators
     of meeting the target), and the code sets "set maxvar 40000".
   - The run is slow because it fits millions of regressions (4000 draws x 3 sexes
     x about 200 locations).

 HOW TO RUN
   1. Run 1_summarize_draws.do, then put country_globe_region_all_ages_draws.dta
      and this do-file in the same folder.
   2. In Stata, set the working directory to that folder, for example:
          cd "C:/path/to/AdoPA-2026/scripts"
   3. Run this do-file.
==============================================================================*/

* housekeeping --------------------------------------------------------------------
clear
clear matrix
clear mata

* about 16,000 variables are held in memory at the peak (see the header)
set maxvar 40000

* number of posterior draws per variable in the input file
local draws = 4000

* years: the baseline year of the target, the last year with estimates, and the
* year of the target (projection year)
local baseyear = 2010
local estyear = 2025
local endyear = 2030

* the input is read from, and the output written to, the working directory
global dir "`c(pwd)'"

* check that the input file is present before starting the (long) run
capture confirm file "$dir/country_globe_region_all_ages_draws.dta"
if _rc {
	display as error "File not found: $dir/country_globe_region_all_ages_draws.dta"
	display as error "Run 1_summarize_draws.do first, and set the working directory to the folder that contains the file."
	exit 601
}

**********************************************************************************************************************

* ---- 1. load the draws -------------------------------------------------------------
use "$dir/country_globe_region_all_ages_draws.dta", clear

* project using only the years with estimates (2010 to 2025)
drop if year > `estyear' | year < `baseyear'


* one variable that identifies each location: the country code for countries and
* the region name for regions and the globe
gen geography = iso3
replace geography = region if geography == ""
levelsof(geography), local(places)

* empty variable that will hold the probit of the prevalence for the current draw
gen probitprev = .

* ---- 2. project each draw to 2030 ----------------------------------------------------
* For each draw, fit a straight line to the probit of the percentage reporting fewer
* than 7 days against year, separately for each sex and location, then:
*   pred2030<draw>      = the projected percentage in 2030 (probit converted back)
*   pp_increase<draw>   = the slope of the line (made into a 0/1 indicator below)
forvalues draw = 1/`draws' {
	di `draw'
	* make empty variables for the 2030 projection and the slope
	quietly gen pred`endyear'`draw' = .
	quietly gen pp_increase`draw' = .
	* probit of this draw's percentage reporting fewer than 7 days (1 minus the
	* prevalence of 7 days)
	quietly replace probitprev = invnormal(1-days7_prob_`draw')
	forvalues sex = 1/3 {
		foreach geo in `places' {
			quietly {
				regress probitprev year if geography == "`geo'" & sex == `sex'
				* e(b)[1,1] is the slope and e(b)[1,2] the constant. Project to 2030
				* and convert the probit back to a percentage.
				replace pred`endyear'`draw' = normal(e(b)[1,2]+`endyear'*e(b)[1,1]) if geography == "`geo'" & sex == `sex'
				replace pp_increase`draw' = e(b)[1,1] if geography == "`geo'" & sex == `sex'
				}
			}
		}
	}

drop probitprev

* ---- 3. set the target for each draw -------------------------------------------------
* keep only the baseline year; the target is based on the baseline (2010) value
keep if year == `baseyear'

* rename the draws of the 7-day prevalence to be the target draws
renpfix days7_prob_ target`endyear'

forvalues draw = 1/`draws' {
	* the target is a 15% relative reduction of the baseline percentage reporting
	* fewer than 7 days (1 minus the prevalence of 7 days)
	replace target`endyear'`draw' = 0.85*(1-target`endyear'`draw')
	* this draw meets the target if its projection for 2030 is below its target
	gen pp_target`draw' = pred`endyear'`draw' < target`endyear'`draw'
	* indicator that the percentage reporting fewer than 7 days is increasing
	replace pp_increase`draw' = (pp_increase`draw'>0)
}

* ---- 4. summarise the draws ----------------------------------------------------------
* The summaries are given the temporary prefix "blah" so that the draws can be
* dropped using wildcards (pp*, pred2030*, target2030*) without dropping the summaries.
	order target`endyear'* pp_increase* pred`endyear'* pp_target*
	* posterior probabilities: the proportion of draws with an increase / meeting the target
	egen blahpp_increase = rowmean(pp_increase1-pp_increase`draws')
	egen blahpp_target = rowmean(pp_target1-pp_target`draws')
	* projection for 2030: mean of the draws and 95% uncertainty interval
	egen blahpred`endyear' = rowmean(pred`endyear'1-pred`endyear'`draws')
	egen blahpred`endyear'_lwr = rowpctile(pred`endyear'1-pred`endyear'`draws'), p(2.5)
	egen blahpred`endyear'_upr = rowpctile(pred`endyear'1-pred`endyear'`draws'), p(97.5)


	* drop the draws and remove the temporary prefix from the summaries
	drop pp* pred`endyear'* target`endyear'* geography p
	renpfix blah
	* the baseline percentage reporting fewer than 7 days (lt7days) is 1 minus the
	* prevalence of 7 days
	* (the lower bound of this percentage is 1 minus the upper bound of the prevalence
	* of 7 days, and vice versa)
	gen lt7days`baseyear' = 1-prev_7days
	gen lt7days`baseyear'_lwr = 1-prev_7days_upr
	gen lt7days`baseyear'_upr = 1-prev_7days_lwr
	drop year prev_7days*

	* ---- 5. on-track assessment and certainty --------------------------------------------
	* the 2030 target is a 15% relative reduction of the baseline
	gen target`endyear' = lt7days`baseyear'*0.85
	* on track if the projected percentage (mean of the draws) is below the target
	gen on_track = pred`endyear' < target`endyear'
	* the certainty is high if the posterior probability of meeting the target is below
	* 0.2 or above 0.8 (the thresholds used throughout the paper)
	gen on_track_certainty = "low"
	replace on_track_certainty = "high" if pp_target < 0.2 | pp_target > 0.8

	* ---- 6. save the results ---------------------------------------------------------------
	export delimited using "$dir/on_track_assessment.csv", replace
