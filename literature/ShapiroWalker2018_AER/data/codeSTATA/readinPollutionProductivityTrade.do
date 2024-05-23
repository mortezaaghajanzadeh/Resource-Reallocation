set trace on
set tracedepth 2
set more off
timer clear
clear
set matsize 10000


sysuse auto
/*===========================================================================================*/
/*                                     Main Program                                          */
/*===========================================================================================*/
capture program drop main
program define main
    paths

    // =============== 0 Comment in/out subprograms you wish to run ================

	*"Table 1: POLLUTION ELASTICITY IV"  	
	pollElasticityCountyIndNAICSiv2, pdf pollutionTrim polluterLinear collapseWeight weight(q) collapseMean logPlus1 abatementCapital anyNonattainment

	*"Table 2: CALCULATE PARAMETERS"
	calculateElasticitySubst, collapseWeight paceOnly
	calculateShapeParameter, collapseWeight paceOnly pdf cutoff(90)
	pollElasticityISICmanual, collapseWeight abatementCapital paceOnly

	
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

	global datadir 	[REDACTED]
	adopath + [REDACTED]
	adopath + [REDACTED]
	
end
//paths

/*---------------------------------------------------------*/
**** Regress pollution/output on abatement expenditure shares
/*---------------------------------------------------------*/
capture program drop pollElasticityCountyIndNAICSiv2
program define pollElasticityCountyIndNAICSiv2			   
syntax, [pdf polluterShare(real 1) polluterLinear balancedPlants laborShare weight(string) pollutionTrim collapseMean  ///
		 collapseWeight paceOnly logPlus1 abateNominal abatementCapital  intoNonattainment ///
		 anyNonattainment totalNonattainment wholeNonattainment partialNonattainment abatementTotal]

	tempfile file
	if "`pdf'"=="pdf"{
		writeln `file' "\documentclass{article}"
		writeln `file' "\usepackage{booktabs}"
		writeln `file' "\usepackage{pdflscape}"
		writeln `file' "\usepackage[margin=1cm]{geometry}"
		writeln `file' "\begin{document}"
	}
	
	use $datadir/[REDACTED]/1_pollutionProductivityTrade_createPACEdataNAICS.dta, clear

	keep if naics>=310000 & naics<=339999

	if "`paceOnly'"=="paceOnly" keep if paceWeight!=.
	
	if "`balancedPlants'"=="balancedPlants"{
		bysort lbdnum naics: egen cellsize=sum(1)
		keep if cellsize==2
	}
	
	if "`abateNominal'"=="abateNominal"{
		replace atot=atot*pimat
		replace paoc=paoc*pimat
		replace totalAbatementExpenditure = atot + acapstock*eqrkl
	}	

	if "`abatementTotal'"=="abatementTotal" replace atot=paoc
	
	*aggregate to naics6-digit level
	gen plants=1
	if "`collapseWeight'"=="collapseWeight"{
		foreach var in totalExpenditures sw q kstst ksteq{
			replace `var'=`var'*wt
		}
		foreach var in alab atot totalAbatementExpenditure acr adep{
			replace `var' = `var'*paceWeight
		}
	}   
	
	destring fipsst cou, replace force
	drop if fipsst==. | cou==. | fipsst==0 | cou==0
	gen fips=fipsst*1000+cou
		
	if "`collapseMean'"=="collapseMean" collapse (mean) totalExpenditures sw q alab atot totalAbatementExpenditure* acr adep plants kstst ksteq, by(fips naics year) fast
	else collapse (sum) totalExpenditures sw q alab atot totalAbatementExpenditure* acr adep plants kstst ksteq , by(fips naics year) fast
	destring naics, replace
	bysort fips naics: keep if _N==2 // keeping a balanced panel county-industry output
	tempfile temp
	save `temp', replace
	
	use $datadir/[REDACTED]/1_pollutionProductivityTrade_neiCountyNAICSaggregate, clear
	keep if year==1990 | year == 2005
	bysort fips naics: keep if _N==2 // keeping a balanced panel in county-industry emissions
	merge 1:1 fips naics year using `temp'
	keep if _m==3
	drop _m
	drop if naics==.
	tempfile temp
	save `temp', replace

	*bring in nonattainment data
	use $datadir/[REDACTED]/1_readinNonattainmentCountyYear`wholeNonattainment'`partialNonattainment'.dta
	tsset fips year
	foreach var of varlist non*{
		replace `var' = 0 if l.`var'==0 & `var'==1 & inlist(year,2004,2005)
	}

	keep if year==1989 | year==2005 
	replace year =1990 if year==1989
	merge 1:m fips year using `temp'
	drop if _m==1
	drop _m
	foreach var of varlist non*{
		replace `var'=0 if `var'==.
	}

	if "`totalNonattainment'"=="totalNonattainment"{
		egen nonattainment=rowtotal(non*)
	}
	else egen nonattainment=rowmax(non*)
	
	*bring in polluter definition
	merge m:1 naics using $datadir/[REDACTED]/1_pollutionProductivityTrade_NEIindustryPolluters
	drop if _m==2
	drop _m
	foreach var in co nox so2 voc pm10 pm25{
		gen polluter`var' = ind`var'share>=`=`polluterShare'/100' if ind`var'share!=.
		if "`polluterLinear'"=="polluterLinear" replace polluter`var' = ind`var'share if ind`var'share!=.
		replace  polluter`var'=0 if ind`var'share==.
	}

	foreach var of varlist polluter*{
		replace `var'=0 if `var'==.
	}
	egen polluter=rowmax(polluter*)
	
	egen id = group(naics fips)
	xtset id year, delta(15)
	*what is the average effect of any newly designated nonattainment status for any polluter

	foreach var in CO NO2 SO2 PM Ozone{
		if "`anyNonattainment'"=="anyNonattainment" replace nonattainment`var'= nonattainment
		if "`intoNonattainment'"=="intoNonattainment" replace nonattainment`var'=l.nonattainment`var' if nonattainment`var'==0 & l.nonattainment`var'!=0 & year==2005
	}
	if "`intoNonattainment'"=="intoNonattainment" replace nonattainment=l.nonattainment if nonattainment==0 & l.nonattainment!=0 & year==2005
	
	
	*define reduced form nonattainment treatment (i.e. being an industry polluting X when county switched into nonattainment for X)
	*and lower order interaction terms
	gen noncoPollco = nonattainmentCO*polluterco
	gen nonnoxPollnox = nonattainmentNO2*polluternox
	gen nono3Pollnox = nonattainmentOzone*polluternox
	gen nonso2Pollso2= nonattainmentSO2*polluterso2
	gen nonpm25Pollpm25 = nonattainmentPM*polluterpm25
	gen nonpm10Pollpm10 = nonattainmentPM*polluterpm10
	gen nono3Pollvoc = nonattainmentOzone*pollutervoc
	
	gen nonPollco = nonattainment*polluterco
	gen nonPollnox = nonattainment*polluternox
	gen nonPollso2= nonattainment*polluterso2
	gen nonPollpm25 = nonattainment*polluterpm25
	gen nonPollpm10 = nonattainment*polluterpm10
	gen nonPollvoc = nonattainment*pollutervoc
	gen nonPoll = nonattainment*polluter	
	
	foreach var in co nox so2 pm25 pm10 voc {
		gen polluter`var'year = polluter`var'*(year==2005)
	}
	gen polluteryear=polluter*(year==2005)
	 
	foreach var of varlist co-pm25W atot{
		replace `var'=0 if `var'==.
	}
	foreach var of varlist co-pm25W{
		if "`logPlus1'"=="logPlus1" replace `var'=1 if `var'==0
	}
	if "`pollutionTrim'"=="pollutionTrim"{
		foreach var in co nox so2 voc pm10 pm25{
			replace `var' = `var'W
		}
	}
	
	gen naics3=floor(naics/1000)
	gen laborRatio=log(1-alab/sw)
	*total abatement expenditures taken as given. alternatively could add in capital stock*rental rate 
	gen totalExpenseRatio = log(1-atot/totalExpenditures)
	if "`abatementCapital'"=="abatementCapital"	replace totalExpenseRatio = log(1-totalAbatementExpenditure/totalExpenditures)
	if "`abatementCapital'"=="abatementCapital"	& "`abatementTotal'"=="abatementTotal" replace totalExpenseRatio = log(1-totalAbatementExpenditurePaoc/totalExpenditures)
	egen totalPoll=rowtotal(co nox so2 voc pm10 pm25)
	gen kst=kstst+ksteq 
	
	foreach var in co nox so2 voc pm10 pm25 totalPoll{
		gen `var'OutputIntensity=log((`var')/q)
	}
	
	if "`weight'"!=""{
		local weightvar "[aweight=`weight'1990]"
		local weightTitle "weight`weight'"
		gen `weight'1990 = `weight' if year==1990
		bysort naics: egen m`weight'1990 =mean(`weight'1990)
		replace `weight'1990 =m`weight'1990
		drop m`weight'1990 
	}
	
	gen double naicsXyear = naics*10000+year
	gen state=floor(fips/1000)
	gen stateXyear = state*10000+year
	if "`n6y'"=="n6y" local 		
			
			
	xtset, clear
	capture est drop _all 
	xi: xtivreg2  coOutputIntensity (totalExpenseRatio=noncoPollco) nonattainmentCO pollutercoyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f1) saverfprefix(rf1)
	eststo est1
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho1 
	xi: xtivreg2  coOutputIntensity totalExpenseRatio nonattainmentCO pollutercoyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols1
	
	xi: xtivreg2  noxOutputIntensity (totalExpenseRatio=nono3Pollnox) nonattainmentOzone polluternoxyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f2) saverfprefix(rf2)
	eststo est2
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho2 
	xi: xtivreg2  noxOutputIntensity totalExpenseRatio nonattainmentOzone polluternoxyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols2
	
	xi: xtivreg2  pm10OutputIntensity (totalExpenseRatio=nonpm10Pollpm10) nonattainmentPM polluterpm10year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f3) saverfprefix(rf3)
	eststo est3
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho3
	xi: xtivreg2  pm10OutputIntensity totalExpenseRatio nonattainmentPM polluterpm10year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols3
	
	xi: xtivreg2  pm25OutputIntensity (totalExpenseRatio=nonpm25Pollpm25) nonattainmentPM polluterpm25year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f4) saverfprefix(rf4)
	eststo est4
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho4 
	xi: xtivreg2  pm25OutputIntensity totalExpenseRatio nonattainmentPM polluterpm25year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols4
	
	xi: xtivreg2  so2OutputIntensity (totalExpenseRatio=nonso2Pollso2) nonattainmentSO2 polluterso2year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f5) saverfprefix(rf5)
	eststo est5
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho5
	xi: xtivreg2  so2OutputIntensity totalExpenseRatio nonattainmentSO2 polluterso2year i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols5
	
	xi: xtivreg2  vocOutputIntensity (totalExpenseRatio=nono3Pollvoc) nonattainmentOzone pollutervocyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f6) saverfprefix(rf6)
	eststo est6
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho6
	xi: xtivreg2  vocOutputIntensity totalExpenseRatio nonattainmentOzone pollutervocyear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols6
	
	xi: xtivreg2  totalPollOutputIntensity (totalExpenseRatio=nonPoll) nonattainment polluteryear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) savefirst saverf savefprefix(f7) saverfprefix(rf7)
	eststo est7
	nlcom (rho: 1/(_b[totalExpenseRatio]+1)), post
	eststo rho7 
	xi: xtivreg2  totalPollOutputIntensity totalExpenseRatio nonattainment polluteryear i.year `weightvar', ///
		cluster(naics3 fips) fe i(id) partial(_I*) 
	eststo ols7
	
*******add in extra column in all tables below*******
 		
	writeln `file' "\begin{table}[p]"
	writeln `file' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"
	writeln `file' "\caption{Pollution Elasticity Parameter: Industry Level IV, by Pollutant}"
	writeln `file' "{ \begin{center}"
	writeln `file' "\begin{tabular}{lcccccccc} "				
	writeln `file' " \toprule  "
	writeln `file' ""
	local rsq="r2_a"
	esttab f1* f2* f3* f4* f5* f6* f7* ///
		using `file', ///
		keep(nonPoll) rename(nono3Pollvoc nonPoll nonpm25Pollpm25 nonPoll nonpm10Pollpm10 nonPoll nono3Pollnox nonPoll noncoPollco nonPoll nonso2Pollso2 nonPoll) /// 
		order(nonPoll) ///
		cells(b(star fmt(%9.3f)) se(par)) stats( ) /// 			
		legend title("") nodepvars nomtitles collabels("")  ///
		coef(nonPoll "Nonattain*Polluter") ///
		prehead(" ") posthead(" 	&	CO	&	NOx (O3) &	PM10	&	PM25	& SO2 & VOC (O3) & Total (Any) \\ \midrule  &\multicolumn{7}{c}{Panel A: First Stage} \\ \cmidrule{2-8} \addlinespace" ) ///
		postfoot("") prefoot("\cmidrule{2-8} ") dropped(---) ///
		nogaps booktabs append se nonotes label wrap  varwidth(45) ///
		substitute([htbp] [p]) star(* 0.10 ** 0.05 *** 0.01)
		
	writeln `file' " &\multicolumn{7}{c}{Panel B: Reduced Form} \\ \cmidrule{2-8}"
	local rsq="r2_a"
	esttab rf1* rf2* rf3* rf4* rf5* rf6* rf7* ///
		using `file', ///
		keep(nonPoll) rename(nono3Pollvoc nonPoll nonpm25Pollpm25 nonPoll nonpm10Pollpm10 nonPoll nono3Pollnox nonPoll noncoPollco nonPoll nonso2Pollso2 nonPoll) /// 
		order(nonPoll) ///
		cells(b(star fmt(%9.3f)) se(par)) stats( ) /// 			
		legend title("") nodepvars nomtitles nonumbers collabels("") ///
		coef(nonPoll "Nonattain*Polluter") ///
		prehead(" ") posthead(" 	") ///
		postfoot("") prefoot("\cmidrule{2-8} ") dropped(---) ///
		nogaps booktabs append se nonotes label wrap  varwidth(45) ///
		substitute([htbp] [p]) star(* 0.10 ** 0.05 *** 0.01)
		
	writeln `file' " &\multicolumn{6}{c}{Panel C: Instrumental Variables} \\ \cmidrule{2-8}"
	esttab est* ///
		using `file', ///
		keep(totalExpenseRatio ) /// 
		order(totalExpenseRatio) ///
		cells(b(star fmt(%9.3f)) se(par)) stats(N rkf, label(N "First Stage F") fmt(%9.0g %9.2g)) /// 			
		legend title("") nodepvars nomtitles nonumbers collabels("") ///
		coef(totalExpenseRatio "Abatement Expenditure Ratio") ///
		prehead(" ") posthead(" ") ///
		postfoot("") prefoot("\cmidrule{2-8} ") dropped(---) ///
		nogaps booktabs append se nonotes label wrap  varwidth(45) ///
		substitute([htbp] [p]) star(* 0.10 ** 0.05 *** 0.01)
		
	writeln `file' " \cmidrule{2-8} "
	writeln `file' " &\multicolumn{7}{c}{Panel D: Pollution Elasticity Parameter} \\ \cmidrule{2-8}"
	esttab rho* ///
		using `file', ///
		keep(rho) /// 
		order(rho) ///
		cells(b(star fmt(%9.3f)) se(par)) stats( ) /// 			
		legend title("") nodepvars nonumbers nomtitles collabels("") ///
		coef(rho "Pollution Elasticity (Rho)") ///
		prehead(" ") posthead(" " ) ///
		postfoot("") prefoot(" ") dropped(---) ///
		nogaps booktabs append se nonotes label wrap  varwidth(45) ///
		substitute([htbp] [p]) star(* 0.10 ** 0.05 *** 0.01)
		
	writeln `file' "\midrule"
	writeln `file' "County-NAICS FE	&	X	&	X 	&	X 	&	X 	&	X 	&	X  &	X \\ " 
	writeln `file' "\bottomrule "
	writeln `file' "\end{tabular}"
	writeln `file' "\end{center} }"
	writeln `file' "{\footnotesize \emph{Notes}: }"
	writeln `file' "\end{table}"												
	writeln `file' "\clearpage"  
	
	writeln `file' "\begin{table}[p]"
	writeln `file' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"
	writeln `file' "\caption{Pollution Elasticity Parameter: Industry Level OLS, by Pollutant}"
	writeln `file' "{ \begin{center}"
	writeln `file' "\begin{tabular}{lcccccccc} "				
	writeln `file' " \toprule  "
	writeln `file' ""
	local rsq="r2_a"
	esttab ols* ///
		using `file', ///
		keep(totalExpenseRatio ) /// 
		order(totalExpenseRatio) ///
		cells(b(star fmt(%9.3f)) se(par)) stats(N, label(N ) fmt(%9.0g %9.2g)) /// 			
		legend title("") nodepvars nomtitles nonumbers collabels("") ///
		coef(totalExpenseRatio "Abatement Expenditure Ratio") ///
		prehead(" ") posthead(" ") ///
		postfoot("") prefoot("\cmidrule{2-7} ") dropped(---) ///
		nogaps booktabs append se nonotes label wrap  varwidth(45) ///
		substitute([htbp] [p]) star(* 0.10 ** 0.05 *** 0.01)
	writeln `file' "\midrule"
	writeln `file' "County-NAICS FE	&	X	&	X 	&	X 	&	X 	&	X 	&	X &	X \\ " 
	writeln `file' "\bottomrule "
	writeln `file' "\end{tabular}"
	writeln `file' "\end{center} }"
	writeln `file' "{\footnotesize \emph{Notes}: }"
	writeln `file' "\end{table}"												
	writeln `file' "\clearpage"	
 
	if "`pdf'"=="pdf"{
		writeln `file' "\end{document}"
	}
	
	!cp `file' $datadir/[REDACTED]/1_pollutionProductivityTrade_pollElasticityCountyIndNAICSiv2`laborShare'`balancedPlants'`paceOnly'`logPlus1'`abateNominal'`abatementCapital'`collapseMean'`collapseWeight'`weightTitle'`anyNonattainment'`totalNonattainment'`wholeNonattainment'`partialNonattainment'`intoNonattainment'`abatementTotal'`pollutionTrim'pollShare`polluterShare'`polluterLinear'.tex
	
	if "`pdf'"=="pdf"{
		cd $datadir/[REDACTED]/
		!pdflatex $datadir/[REDACTED]/1_pollutionProductivityTrade_pollElasticityCountyIndNAICSiv2`laborShare'`balancedPlants'`paceOnly'`logPlus1'`abateNominal'`abatementCapital'`collapseMean'`collapseWeight'`weightTitle'`anyNonattainment'`totalNonattainment'`wholeNonattainment'`partialNonattainment'`intoNonattainment'`abatementTotal'`pollutionTrim'pollShare`polluterShare'`polluterLinear'.tex
		cd $datadir/[REDACTED]/
	}
	
		
end
  
 
/*---------------------------------------------------------*/
**** Total Input Costs Divided by Total Revenue
/*---------------------------------------------------------*/
capture program drop calculateElasticitySubst
program define calculateElasticitySubst
syntax, [balancedPlants collapseMean collapseWeight paceOnly ]

	use $datadir/[REDACTED]/1_pollutionProductivityTrade_createPACEdataISIC.dta, clear

	if "`paceOnly'"=="paceOnly" keep if paceWeight!=.
	
	if "`balancedPlants'"=="balancedPlants"{
		bysort lbdnum isicrev3 myindustry naics1997: keep if _N==2
	}
	
	gen plants=1
	if "`collapseWeight'"=="collapseWeight"{
		foreach var in totalExpenditures sw q kstst ksteq{
			replace `var'=`var'*wt
		}
		foreach var in alab atot totalAbatementExpenditure acr adep{
			replace `var' = `var'*paceWeight
		}
	}   
	
	destring fipsst cou, replace force
	drop if fipsst==. | cou==. | fipsst==0 | cou==0
	
		
	levelsof myindustry, local(ind)
	gen count=1
	
	
	*DROPPING P99+ EXPENDITURE/OUTPUT PLANTS
	foreach var in totalExpenditures q{
		sum `var', d
		replace `var'=. if `var'>=r(p99) & `var'!=.
	}
	drop if totalExpenditures==. | q==.

	
	if "`collapseMean'"=="collapseMean" collapse (mean) totalExpenditures q count, by(myindustry year) fast
	else collapse (sum) totalExpenditures q count, by(myindustry year) fast
	 
	keep if year==1990
		 
	gen inputShare=totalExpenditures/q
	gen sigma = 1/(1-inputShare)

	keep myindustry inputShare* sigma* count
	list

	outsheet using $datadir/[REDACTED]/1_readinPPT_calculateElasticitySubst`balancedPlants'`paceOnly'`collapseMean'`collapseWeight'`weightTitle'.csv, comma replace
		
end


/*---------------------------------------------------------*/
**** Shape Parameter: sales rank on log of sales
/*---------------------------------------------------------*/
capture program drop calculateShapeParameter
program define calculateShapeParameter
syntax, [balancedPlants collapseMean collapseWeight paceOnly pdf cutoff(integer 75)]

	use if year==1990 using $datadir/[REDACTED]/1_pollutionProductivityTrade_createPACEdataISIC.dta, clear

	if "`paceOnly'"=="paceOnly" keep if paceWeight!=.
	
	if "`balancedPlants'"=="balancedPlants"{
		bysort lbdnum isicrev3 myindustry naics1997: keep if _N==2
	}
	
	gen plants=1
	if "`collapseWeight'"=="collapseWeight"{
		foreach var in totalExpenditures sw q kstst ksteq{
			replace `var'=`var'*wt
		}
		foreach var in alab atot totalAbatementExpenditure acr adep{
			replace `var' = `var'*paceWeight
		}
	}   
	
	destring fipsst cou, replace force
	drop if fipsst==. | cou==. | fipsst==0 | cou==0
	
	levelsof myindustry, local(ind)
	gen count=1

	tempfile file
	if "`pdf'"=="pdf"{
		writeln `file' "\documentclass{article}"
		writeln `file' "\usepackage{booktabs}"
		writeln `file' "\usepackage{pdflscape}"
		writeln `file' "\usepackage[margin=1cm]{geometry}"
		writeln `file' "\begin{document}"
	}
	
	writeln `file' "\begin{table}[p]"
	writeln `file' "\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"
	writeln `file' "\caption{Shape Parameter Estimates}"
	writeln `file' "{ \begin{center}"
	writeln `file' "\begin{tabular}{cccccccc} "				
	writeln `file' " \toprule  "
	writeln `file' " \addlinespace"
	writeln `file' " &	Regression &	Regression & 	Input	&	Elasticity		& Shape &	Shape 	&		\\"
	writeln `file' " Industry	&	Estimate	&	Standard Error	& Share	&	of Substitution &  Estimate	&	Standard Error	&	N	\\"
	writeln `file' " \midrule"
	
	*DROPPING P99+ EXPENDITURE/OUTPUT PLANTS
	foreach var in totalExpenditures q{
		sum `var', d
		replace `var'=. if `var'>=r(p99) & `var'!=.
	}
	drop if totalExpenditures==. | q==.	
	preserve
	collapse (sum) totalExpenditures q , by(myindustry)
	gen inputShare=totalExpenditures/q
	gen sigma = 1/(1-inputShare)
	foreach var in `ind'{
		sum sigma if myindustry==`var'
		local `var' `r(mean)'
		sum inputShare if myindustry==`var'
		local `var'b `r(mean)'
	}
	restore
	
	// only looking at domestic shipments as per di Giovanni and Levchenko
	replace q=q-exp
	
	collapse (sum) q, by(firmid isicrev3 myindustry naics1997 year)
	gen lq=log(q)
	drop if lq==.	
	tempfile temp
	save `temp', replace
	
	// visual evidence suggests [REDACTED] is where linearity starts to happen
	sum lq, d
	local cut = (`r(max)'-`r(min)')/10
	egen group = cut(lq), at(`r(min)'(`cut')`r(max)')
	egen salesRank = rank(q), field
	replace salesRank = 1-(_N-salesRank)/_N
	gen logSalesRank=log(salesRank)
	collapse (mean) lq logSalesRank , by(group )
	tw scatter logSalesRank lq || lfit logSalesRank lq

	foreach var in `ind'{
		use if myindustry==`var' using `temp', clear
		egen salesRank = rank(q), field
		replace salesRank = 1-(_N-salesRank)/_N
		gen logSalesRank=log(salesRank)
		sum q, d
		drop if q<r(p`cutoff')
		sum lq
		local cut = (`r(max)'-`r(min)')/10
		egen group = cut(lq), at(`r(min)'(`cut')`r(max)')
		preserve
		collapse (mean) logSalesRank lq , by(group)
		tw scatter logSalesRank lq || lfit logSalesRank lq 
		restore
		reg logSalesRank lq, robust
		nlcom (theta: _b[lq]*(``var''-1)) (beta: _b[lq]), post 
		local beta`var': display %9.2f _b[theta]
		local se`var': display %4.2f _se[theta]
		local se`var' "(`se`var'')"
		local beta`var'2: display %9.2f _b[beta]
		local se`var'2: display %4.2f _se[beta]
		local se`var'2 "(`se`var'2')"
		local N`var': display %9.0f e(N)
		local elasticity: display %4.2f ``var''
		local inputShare: display %4.2f ``var'b'
		writeln `file' "`var' 	&	`beta`var'2'	&	`se`var'2'	&	`inputShare'	& `elasticity'	& `beta`var''	&	`se`var''	&	`N`var''	\\ "
	}
	
	writeln `file' "\bottomrule "
	writeln `file' "\end{tabular}"
	writeln `file' "\end{center} }"
	writeln `file' "{\footnotesize \emph{Notes}: Standard errors clustered by four digit NAICS.}"
	writeln `file' "\end{table}"
	
	if "`pdf'"=="pdf"{
		writeln `file' "\end{document}"
	}
	
	!cp `file' $datadir/[REDACTED]/1_readinPPT_calculateShapeParameter`balancedPlants'`paceOnly'`collapseMean'`collapseWeight'`weightTitle'cutoffP`cutoff'.tex
	cd $datadir/[REDACTED]/
	!pdflatex $datadir/[REDACTED]/1_readinPPT_calculateShapeParameter`balancedPlants'`paceOnly'`collapseMean'`collapseWeight'`weightTitle'cutoffP`cutoff'.tex
	cd $datadir/[REDACTED]/
	
end	




/*---------------------------------------------------------*/
**** Regress pollution/output on abatement expenditure shares
/*---------------------------------------------------------*/
capture program drop pollElasticityISICmanual
program define pollElasticityISICmanual
syntax, [balancedPlants collapseMean collapseWeight paceOnly abateNominal abatementCapital abatementTotal]

	use $datadir/[REDACTED]/1_pollutionProductivityTrade_createPACEdataISIC.dta, clear

	if "`paceOnly'"=="paceOnly" keep if paceWeight!=.
	
	if "`balancedPlants'"=="balancedPlants"{
		bysort lbdnum isicrev3 myindustry naics1997: keep if _N==2
	}
	
	if "`abateNominal'"=="abateNominal"{
		replace atot=atot*pimat
		replace paoc=paoc*pimat
		replace totalAbatementExpenditure = atot + acapstock*eqrkl
	}	

	if "`abatementTotal'"=="abatementTotal" replace atot=paoc
	
	
	gen plants=1
	if "`collapseWeight'"=="collapseWeight"{
		foreach var in totalExpenditures sw q kstst ksteq{
			replace `var'=`var'*wt
		}
		foreach var in alab atot totalAbatementExpenditure acr adep{
			replace `var' = `var'*paceWeight
		}
	}   
	
	destring fipsst cou, replace force
	drop if fipsst==. | cou==. | fipsst==0 | cou==0
	
	log using 1_readinDataDisclosure.log, replace
	levelsof myindustry, local(ind)
	foreach var in `ind'{ 
		preserve
		keep if myindustry==`var'
		display "**********************************Industry-`var'**********************************"
		disclose_nk_p, disclosureVar(tvs)
		restore
	}
	
	if "`collapseMean'"=="collapseMean" collapse (mean) totalExpenditures sw q alab atot totalAbatementExpenditure acr adep plants kstst ksteq, by(myindustry year) fast
	else collapse (sum) totalExpenditures sw q alab atot totalAbatementExpenditure acr adep plants kstst ksteq , by(myindustry year) fast
	 
	gen totalExpenseRatio = atot/totalExpenditures
	if "`abatementCapital'"=="abatementCapital"	replace totalExpenseRatio = totalAbatementExpenditure/totalExpenditures
	if "`abatementCapital'"=="abatementCapital"	& "`abatementTotal'"=="abatementTotal" replace totalExpenseRatio = log(1-totalAbatementExpenditurePaoc/totalExpenditures)
	
	keep if year==1990
	keep myindustry totalExpenseRatio 
	 
	outsheet using $datadir/[REDACTED]/1_pollutionProductivityTrade_pollElasticityISICmanual`balancedPlants'`paceOnly'`abateNominal'`abatementCapital'`collapseMean'`collapseWeight'`weightTitle'`abatementTotal'.txt, replace comma
		
end

 
 
 *---------------------------------------------------------*/
 /* Run Main Program                                       */
 /*--------------------------------------------------------*/


main
//main program
exit
