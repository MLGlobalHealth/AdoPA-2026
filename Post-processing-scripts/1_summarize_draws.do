/*==============================================================================
 1_summarize_draws.do

 Part of the code accompanying:
   Strain T, Flaxman SR, Ricardo LIC, et al, and the Country Data Author Group.
   National, regional, and global trends in adolescent physical activity from
   2010 to 2025: a pooled analysis of 389 surveys with 3.2 million participants
   from 140 countries. Lancet Glob Health 2026.
   Code repository: https://github.com/MLGlobalHealth/AdoPA-2026

 PURPOSE
   Takes the posterior draws of the Bayesian cumulative ordinal model (one set of
   draws per country, sex, single year of age and year) and computes summary
   statistics with uncertainty for:
     - each country,
     - regions (the eight analysis regions, WHO regions, and World Bank income
       groups) and the globe, obtained by population-weighting the country draws.
   Statistics are produced for boys (sex = 1), girls (sex = 2) and both sexes
   combined (sex = 3). The statistics are:
     - the percentage of adolescents reporting 0, 1, ..., 7 active days per week
       (prev_0days ... prev_7days; stored as proportions between 0 and 1), and
     - the mean number of active days per week (mean_days_mean).
   An "active day" is a day with at least 60 minutes of moderate-to-vigorous
   physical activity. Each statistic has a point estimate (mean of the draws) and
   a 95% uncertainty interval (2.5th and 97.5th percentiles of the draws; _lwr and
   _upr).

 INPUTS (all in the same folder as this do-file; see "HOW TO RUN")
   mcmc_draws_<model>_all_countries.csv
                                     Posterior draws from the model, written by
                                     fit.r to results/. One row per iso3 x year x
                                     age x sex. Columns: iso3, year, age, sex, then
                                     one column per posterior draw for each of
                                     days0_prob ... days7_prob (probability of
                                     reporting 0 ... 7 active days) and mean_days
                                     (mean active days), named days0_prob_1,
                                     days0_prob_2, ... The number of draws is read
                                     from the file. This file is large and is NOT
                                     included in the repository.
   countries regions.dta             Lookup of each of the 195 countries (iso3) to its
                                     analysis region (adopa_region), World Bank
                                     income group (wbinc25, 2025 classification;
                                     blank for 2 countries with no classification)
                                     and WHO region (whoreg6; 7_NonMS is a country
                                     that is not a WHO Member State). The file also
                                     holds the country name (whoname), which is used
                                     by the results scripts.
   WPP24.dta                         UN World Population Prospects 2024 population
                                     (variable p, in thousands) by iso3, sex
                                     (1 boys, 2 girls), single year of age (11-17)
                                     and year (2000-2030), used as weights.
   countries_with_incl_surveys.dta   List of the 140 countries (iso3) that have at
                                     least one included survey; used to flag which
                                     country estimates are published.

 OUTPUTS (written to the same folder)
   stats_by_age.dta / .csv           Statistics by single year of age (11-17)
   stats_all_ages.dta / .csv         Statistics for all ages 11-17 combined
   country_globe_region_all_ages_draws.dta
                                     Draws of the 7-day prevalence (days7_prob_1 to
                                     days7_prob_<draws>) for all ages combined,
                                     2010-2025, for countries, regions and the
                                     globe, with the summary statistics. It is the
                                     input to 2_projections_and_target.do.

 OUTPUT VARIABLES (stats_*.dta / .csv)
   year, iso3, sex (1 boys, 2 girls, 3 both), age (stats_by_age only),
   region and regiontype (adopa_region, wbinc25, whoreg6, Global or country),
   p (adolescent population in thousands, used for weighting), mean_days_mean/_lwr/_upr,
   prev_0days ... prev_7days with _lwr/_upr, and publish ("yes"/"no"; see below).

 REQUIREMENTS
   - The draws file has 4 + 9 x (number of draws) columns and is loaded in full.
     With fit.r's default settings (800 draws) that is 7,204 columns, which runs in
     Stata/SE or MP. The published run used 4000 draws (36,004 columns, about
     16 GB), which needs Stata/MP and "set maxvar 40000".
   - The run takes a long time because of the loops over the draws.

 HOW TO RUN
   1. Put this do-file, the three .dta inputs and the draws .csv written by fit.r
      (results/mcmc_draws_<model>_all_countries.csv) in one folder.
   2. In Stata, set the working directory to that folder, for example:
          cd "C:/path/to/AdoPA-2026/Post-processing-scripts"
   3. Set the model name (the <model> part of the draws file name, e.g.
      model_20260224_0731) at the top of the code below, then run this do-file.
==============================================================================*/

* housekeeping --------------------------------------------------------------------
* temporary files used to hold the country and aggregate results during the run
tempfile countries globe _adopa_region _wbinc25 _whoreg6 all
clear

* the draws file has 4 + 9 x draws columns, which can exceed Stata's default limit
* on the number of variables (set maxvar 40000 needs Stata/MP; fall back to the
* Stata/SE maximum otherwise)
capture set maxvar 40000
if _rc capture set maxvar 32767

* name of the model run: the <model> part of results/mcmc_draws_<model>_all_countries.csv
* written by fit.r (the model name "model" plus the run timestamp)
local model = "model_YYYYMMDD_HHMM"

* the number of posterior draws per variable is read from the draws file after loading

* all inputs are read from, and all outputs are written to, the working directory
global dir "`c(pwd)'"

* check that every input file is present before starting the (long) run
foreach f in "countries regions.dta" "WPP24.dta" "countries_with_incl_surveys.dta" "mcmc_draws_`model'_all_countries.csv" {
	capture confirm file "$dir/`f'"
	if _rc {
		display as error "File not found: $dir/`f'"
		display as error "Set the working directory to the folder containing this do-file and its input files (see the header)."
		exit 601
	}
}

**********************************************************************************************************************

* The code below is run twice, once for each treatment of age:
*   by_age    = statistics for each single year of age (11 to 17)
*   all_ages  = statistics for all ages 11 to 17 combined (the main results)
local agetreatments = "by_age all_ages"

foreach type in `agetreatments' {

	* ---- 1. load the country-level draws and add regions and population ------------
	import delimited "$dir/mcmc_draws_`model'_all_countries.csv", clear

	* number of posterior draws per variable: 4 identifier columns, then 9 variables per draw
	local draws = (c(k) - 4) / 9
	display "Posterior draws per variable: `draws'"

	* add each country's analysis region, WHO region and World Bank income group
	merge m:1 iso3 using "$dir/countries regions.dta", keep(master match) keepusing(adopa_region wbinc25 whoreg6) nogen

	* add the adolescent population (variable p) of each country, sex, age and year
	merge 1:1 iso3 sex age year using "$dir/WPP24.dta", keep(master match) nogen

	* Population-weighted aggregation is done separately within every draw so that the
	* uncertainty of the model carries through to the regional and global estimates.
	* Step 1: multiply each draw by the population. Summing these products over
	* countries gives the numerator of a population-weighted mean; the denominator
	* (the summed population) is divided out in step 4 below.
	local d = 0
	while `d' < `draws' {
		local d = `d' +1
		foreach day of numlist 0/7 {
			replace days`day'_prob_`d' = days`day'_prob_`d'*p
		}
		replace mean_days_`d' = mean_days_`d'*p
	}


	* keep the population-weighted country draws; they are also reported as countries
	save `countries', replace


	* ---- 2. aggregate countries into regions and the globe -----------------------
	* Step 2: sum the population-weighted draws within each region. This is done for
	* each of the three region classifications:
	*   adopa_region = the eight analysis regions used in the paper
	*   wbinc25      = World Bank income groups (2025 classification)
	*   whoreg6      = WHO regions
	local regions = "adopa_region wbinc25 whoreg6"

	foreach regiontype in `regions' {
		preserve
		collapse (sum) mean_days* days* p, by(year sex `regiontype' age) fast
		ren `regiontype' region
		gen regiontype = "`regiontype'"
		save `_`regiontype'', replace
		restore
	}

	* Globe: sum over all countries
	collapse (sum) mean_days* days* p, by(year sex age) fast
	gen region = "Global"
	gen regiontype = "Global"

	* append the three sets of region values to the globe values
	foreach regiontype in `regions' {
		append using `_`regiontype''
		}

	* append the country values (rows with an iso3 code are countries)
	append using `countries'
	replace regiontype = "country" if iso3 != ""

	save `all', replace

	* ---- 3. both sexes combined -----------------------------------------------------
	* Step 3: sum boys and girls to get both sexes (sex = 3), then add them to the
	* sex-specific rows. The sums are of population-weighted draws, so the result is
	* weighted by the number of boys and girls.
	collapse (sum) mean_days* days* p, by(year age region regiontype iso3) fast
	gen sex = 3

	append using `all'

	* Treatment of age. This is the only part of the code that differs between
	* by_age and all_ages: ages are either kept separate (by_age) or summed
	* (all_ages). Summing the population-weighted draws over ages weights the ages
	* by their population.
	if "`type'" == "by_age" {
		collapse (sum) mean_days* days* p, by(iso3 year sex region regiontype age) fast
	}

	if "`type'" == "all_ages" {
		collapse (sum) mean_days* days* p, by(iso3 year sex region regiontype) fast
	}


	* ---- 4. convert the weighted sums back to population-weighted averages ----------
	* Step 4: divide each draw by the total population, so that each draw is now the
	* population-weighted probability (or mean) for the region, country, sex or age
	* group. For a single country and sex this reverses step 1.
	local d = 0
	while `d' < `draws' {
		local d = `d' +1
		foreach day of numlist 0/7 {
			replace days`day'_prob_`d' = days`day'_prob_`d'/p
		}
		replace mean_days_`d' = mean_days_`d'/p
	}


	* ---- 5. summarise the draws: point estimate and 95% uncertainty interval --------
	* mean active days per week: mean of the draws, and 2.5th / 97.5th percentiles
	egen mean_days_mean = rowmean(mean_days_1-mean_days_`draws')
	egen mean_days_lwr = rowpctile(mean_days_1-mean_days_`draws'), p(2.5)
	egen mean_days_upr = rowpctile(mean_days_1-mean_days_`draws'), p(97.5)


	* prevalence of each number of active days (0 to 7), as a proportion
	foreach day of numlist 0/7 {
		egen prev_`day'days = rowmean(days`day'_prob_1-days`day'_prob_`draws')
		egen prev_`day'days_lwr = rowpctile(days`day'_prob_1-days`day'_prob_`draws'), p(2.5)
		egen prev_`day'days_upr = rowpctile(days`day'_prob_1-days`day'_prob_`draws'), p(97.5)
	}

	* the draws for 1 to 6 active days are not needed again; drop them to save memory
	drop days1_prob_1-days6_prob_`draws'


	* ---- 6. save the draws used for the projections -----------------------------------
	* For the all-ages treatment, save the draws of the 7-day prevalence for 2010-2025
	* (the estimation period). 2_projections_and_target.do uses them to project to
	* 2030 and to compute the posterior probability of meeting the target. The draws
	* of the 0-day prevalence and of mean days, and the summary statistics of mean
	* days, are not needed there, so they are dropped from this file to save space
	* (they are still in the stats_*.dta files saved below). Note that this file is
	* still large.
	if "`type'" == "all_ages" {
		preserve
		keep if year >= 2010 & year < 2026
		drop days0_prob_1-days0_prob_`draws' mean_days_*
		save "$dir/country_globe_region_all_ages_draws.dta", replace
		restore
	}

	* the draws are no longer needed
	drop days0_prob_1-days7_prob_`draws' mean_days_1-mean_days_`draws'


	* ---- 7. flag which estimates are published --------------------------------------
	* Estimates are computed for all countries (they feed the regional and global
	* values), but country-level estimates are published only for countries with at
	* least one included survey. Regional and global rows (which have no iso3 code)
	* are always published.
	merge m:1 iso3 using "$dir/countries_with_incl_surveys.dta", keepusing(whoname)
	gen publish = "yes"
	replace publish = "no" if _merge == 1 & regiontype == "country"
	drop _merge whoname

	* ---- 8. save the results ----------------------------------------------------------
	save "$dir/stats_`type'.dta", replace
	export delimited using "$dir/stats_`type'.csv", replace
}
