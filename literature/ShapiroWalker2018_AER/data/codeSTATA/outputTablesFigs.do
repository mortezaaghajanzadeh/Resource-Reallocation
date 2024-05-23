clear all
set trace on
set tracedepth 2
set more off
set seed 123
set type double, permanently
version 14.0
set matsize 10000
set maxvar 15000
set linesize 120


sysuse auto
/*===========================================================================================*/
/*                                     Main Program                                          */
/*===========================================================================================*/
capture program drop main
program define main
    paths

    // =============== 0 Comment in/out subprograms you wish to run ================

	outputF1
	outputF2
	outputF3
	outputF6
	outputT2
	

end
// program main


/*===========================================================================================*/
/*                                    Sub Programs                                           */
/*===========================================================================================*/
 
/*---------------------------------------------------------*/
/* Define Path Macros 					                   */
/*---------------------------------------------------------*/
capture program drop paths
program define paths
	
*	adopath + /Users/reedwalker/Desktop/Dropbox/ado
*	adopath + /bulk/fac/rwalker/ado
*	global datadir "/home/faculty/rwalker/bulk/pollutionTrade/"
*	global datadir `"/Users/reedwalker/Desktop/Dropbox/pollution_productivity_trade/pollution_productivity_trade"'
	global datadir "E:/Dropbox/pollution_productivity_trade/replicationFiles"
	global inventoryYears "1990 1996 1999 2002 2005 2008"
	
end	
//paths;



capture program drop spliceNaicsSic
program define spliceNaicsSic
	* `1' here is the variable we're splicing naics/sic for

	* define variable as proportional change since 1990
	gen  `1'90  = `1' if year == 1990
	egen `1'90b = max(`1'90)
	replace `1' = `1' / `1'90b

	* sic/naics change in 1997 can make series discontinuous. 
	* asm v. census can also make series discontinuous
	* so, calculate 1995-1996 trend, compare its linear prediction for 1998
	* against actual 1998, call that ratio the splice, and multiple 1998+ by it
	gen trend`1' = l2.`1' - l3.`1'
	gen  `1'N = l2.`1' + 2 * trend`1' if year == 1998
	gen  `1'splice = `1'N / `1' if year == 1998
	egen `1'spliceM = max(`1'splice)
	la var `1'spliceM "ratio between listed 1998 value and the 1998 value fitted from 95-96 linear trend"
	replace `1'N = `1' if year < 1997
	* now adjust 1998+ naics years by this splice value
	replace `1'N = `1' * `1'spliceM if year > 1997
	* exclude 1997 since it's a census year (below will exclude all census years)
	replace `1'N = . if year == 1997
	drop `1'splice `1'spliceM trend`1' `1'90 `1'90b
	replace `1' = `1'N
end


capture program drop outputF1
program define outputF1

	u "$datadir/dataSTATA/f1.dta", clear
	
	* main graph: separately by pollutant
	tw  (line      vship year, sort lcolor(blue)) 											///
		(connected co    year, sort lcolor(red)   lpattern(dash) msymb(D)  mcolor(red)) 	///
		(connected nox   year, sort lcolor(red)   lpattern(dash) msymb(T)  mcolor(red)) 	///
		(connected pm25  year, sort lcolor(green) lpattern(dash) msymb(S)  mcolor(green)) 	///
		(connected pm10  year, sort lcolor(green) lpattern(dash) msymb(Dh) mcolor(green)) 	///
		(connected so2   year, sort lcolor(navy)  lpattern(dash) msymb(Th) mcolor(navy)) 	///
		(connected voc   year, sort lcolor(navy)  lpattern(dash) msymb(Sh) mcolor(navy)), 	///
		graphr(color(white)) yscale(noline) 												///
		legend(order(1 "Real Output (2008$)" 												///
		2 "CO" 3 "NO{subscript:x}" 4 "PM{subscript:2.5}" 5 "PM{subscript:10}" 6 "SO{subscript:2}" 7 "VOCs") rows(2)) 						///
		xtit("") 																		    ///
		ytit("1990=100") ylabel(#6)			

	graph export "$datadir/figures/f1.eps", replace


end

* tw  (line      vship year, sort lcolor(blue)) (line so2   year, sort lcolor(cranberry)  lpattern(dash) mcolor(navy)), graphr(color(white)) yscale(noline) legend(order(1 "Real Output (2008$)" 2 "Sulfur Dioxide")) xtit("") ytit("1990=100") xlab(1990 1995 2000 2005 2008)

capture program drop outputF2
program define outputF2

	insheet using "$censusdir/pollutionHeterogeneity.txt", names clear

	capture graph drop _all
	
	tw scatter noxintensityltfp ltfp || lfit noxintensityltfp ltfp || , title("") ///
		name(noxtfp) legend(off) xtitle("Log Total Factor Productivity") ///
		ytitle("Log NO{subscript:x} Emissions Per Dollar Real Output") ///
		graphregion(fcolor(white))			
	graph export "$datadir/figures/2_readinProductivityScatter_tfpnox.eps", replace
		

end




capture program drop outputF3
program define outputF3
syntax, []
	
	*bring in pollution data
	local co   = "CO"
	local nox  = "NO{subscript:x}"
	local pm25 = "PM{subscript:2.5}"
	local pm10 = "PM{subscript:10}" 
	local so2  = "SO{subscript:2}"
	local voc  = "VOC"
	
	use "$datadir/dataSTATA/1_readinNEIfacilityCollapseSIC.dta", clear
	keep if sic >= 2000 & sic <= 3999
	collapse (sum) co nox pm10 pm25 so2 voc, by(year)
	foreach var in co nox pm10 pm25 so2 voc{
		sum `var' if year == 1990
		replace `var' = `var' / r(mean)
	}	
	tempfile temp
	save "`temp.dta'", replace

	* scale effect from product totals is ...sharepv.txt
	* pv is total value of specified products per year
	* compoutput variables are composition effect using product output
	// bring in composition effect using product totals 
	capture graph drop _all
	insheet using "$censusdir/1_createProductDecompositionproductSharepv.txt", names clear
	merge 1:1 year using "`temp.dta'", keep(1 3) nogen
	save "`temp.dta'", replace
	
	// bring in composition effect using industry totals 
	insheet using "$censusdir/1_createIndustryDecomposition.txt", names clear
	foreach var of varlist comp*{
		rename `var' `var'ind
	}
	merge 1:1 year using "`temp.dta'", keep(1 3) nogen

	
	tsset year
	foreach var of varlist compoutput* pv {	
		spliceNaicsSic `var'
	}

	foreach var of varlist compoutput* pv vship{
		replace `var' = . if inlist(year,1992,1997,2002,2007)
	}

	foreach var of varlist co nox pm10 pm25 so2 voc compoutput* pv vship{
		replace `var'  = `var' * 100
	}
	
	* scale effect from nber-ces is this ...vship.txt
	capture graph drop _all
	insheet using "$censusdir/1_createProductDecompositionproductSharevship.txt", names clear
	merge 1:1 year using "`temp.dta'", keep(1 3) nogen
	
	tsset year
	foreach var of varlist compoutput* {		
		spliceNaicsSic `var'
	}

	
	* figure 3: nitrogen oxides emissions from united states manufacturing (main decomposition)
	foreach var in nox{	
		tw     line `var'            year, lw(medthick) lp(shortdash_dot) ///
			|| line  compoutput`var' year, lw(medthick) lp(dash)  		  ///
			|| line vship            year, lw(medthick)  ||,  			  ///
			legend(label(1 "Scale, Composition, & Technique (NEI)") 	  ///
			       label(2 "Scale & Composition (NEI+Census)") 			  ///
				   label(3 "Scale (Census)") col(1) order(3 2 1) 		  ///
			bplacement(sw) ring(0))  									  ///
			xtitle(Year) ytitle("1990 = 100") 							  ///
			graphr(color(white))
		graph export "$censusdir/2_readinProductDecompositionspliceexcludecmf_Sharepv_`var'.eps", replace
	}	
		

end



* Figure 6: NOx Pollution Tax Changes as a Function fo NOx Budget Trading Program Status
capture program drop outputF6
program define outputF6

	insheet using "$censusdir/1105_walker_release_req3719_20140821/eventStudy.txt", names clear
	
	tw connect estimate year, lw(medthick) || rline min95 max95 year , lp(dash) lc(maroon) || , ///
		legend(off) xtitle("Year") ///
		ytitle("Predicted Marginal Effect of NBP on Pollution Tax") ///
		yline(0, lc(black)) xline(2002.5, lc(black)) ///
		graphr(color(white))
		
	graph export "$datadir/figures/2_nbpEventStudy.eps", replace



end


* Table 2
capture program drop outputT2
program define outputT2
	
	u "$datadir/dataSTATA/t2_at3.dta", clear

	foreach p in co nox pm10 pm25 so2 voc total_poll {
		g poll_per_cost_`p' = `p' / costs
		egen T = total(poll_per_cost_`p')
		g mean_poll_per_cost_`p' = T / 17
		g alpha_`p' = 0.011 * (poll_per_cost_`p' / mean_poll_per_cost_`p')
		drop T mean_poll_per_cost_`p'
		}

	outsheet poll_per_cost_total_poll alpha_total_poll using "$datadir/tables/t2_cols12.txt", replace
end




/*---------------------------------------------------------*/
 /* Run Main Program                                       */
 /*--------------------------------------------------------*/
capture program drop addstat
program addstat, eclass
ereturn scalar `1' = `2'
end

main
//main program
exit

