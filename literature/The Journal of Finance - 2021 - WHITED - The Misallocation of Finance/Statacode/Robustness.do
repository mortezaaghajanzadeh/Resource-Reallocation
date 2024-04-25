* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

*===================*
* Size Intersection *
*===================*

* Create intersection sample first before computing standard errors

* U.S.

gen D_si = inflator*lt
gen E_si = inflator*(lse-lt)
drop if D_si <= 0 | D_si == . | E_si <= 0 | E_si == .

gen industry3 = floor(sic/10)
gen industry2 = floor(sic/100)
gen industry1 = floor(sic/1000)

gen pay_emp = pay/emp_NBER
bysort fyear industry3: egen pay_emp_industry3 = mean(pay/emp_NBER)
replace pay_emp = pay_emp_industry3 if pay_emp == .
bysort fyear industry2: egen pay_emp_industry2 = mean(pay/emp_NBER)
replace pay_emp = pay_emp_industry2 if pay_emp == .
bysort fyear industry1: egen pay_emp_industry1 = mean(pay/emp_NBER)
replace pay_emp = pay_emp_industry1 if pay_emp == .

gen PF_si = inflator*(oibdp+pay_emp*emp)
*gen PF_si = inflator*(prcc_f*csho)+D_si
drop if PF_si <= 0 | PF_si == .

gen assets = inflator*at

scalar num_years = 9
forvalues jj=1/`=num_years'{
 winsor assets if fyear == `jj'+1998, gen(assets`jj') h(10)
}

drop if assets1 != . & (assets1 < 4.645215 | assets1 > 1508.094)
drop if assets2 != . & (assets2 < 5.442798 | assets2 > 1790.724)
drop if assets3 != . & (assets3 < 4.711926 | assets3 > 2021.394)
drop if assets4 != . & (assets4 < 5.567475 | assets4 > 2023.837)
drop if assets5 != . & (assets5 < 5.713070 | assets5 > 2000.106)
drop if assets6 != . & (assets6 < 5.307580 | assets6 > 2104.040)
drop if assets7 != . & (assets7 < 6.622001 | assets7 > 3029.363)
drop if assets8 != . & (assets8 < 7.212945 | assets8 > 3446.628)
drop if assets9 != . & (assets9 < 8.435747 | assets9 > 3767.053)

* China

gen D_si = inflator*totalliabilities/1000
gen E_si = inflator*totalequity/1000
drop if D_si <= 0 | D_si == . | E_si <= 0 | E_si == .

gen industry3 = floor(sector/10)

gen PF_si = inflator*(totalprofit+vatpayable+wagespayable+depreciationcurrentyear)/1000
drop if PF_si <= 0 | PF_si == .

gen assets = inflator*totalassets/1000/8

scalar num_years = 9
forvalues jj=1/`=num_years'{
 winsor assets if fyear == `jj'+1998, gen(assets`jj') h(10)
}

drop if assets1 != . & (assets1 < 4.645215 | assets1 > 1508.094)
drop if assets2 != . & (assets2 < 5.442798 | assets2 > 1790.724)
drop if assets3 != . & (assets3 < 4.711926 | assets3 > 2021.394)
drop if assets4 != . & (assets4 < 5.567475 | assets4 > 2023.837)
drop if assets5 != . & (assets5 < 5.713070 | assets5 > 2000.106)
drop if assets6 != . & (assets6 < 5.307580 | assets6 > 2104.040)
drop if assets7 != . & (assets7 < 6.622001 | assets7 > 3029.363)
drop if assets8 != . & (assets8 < 7.212945 | assets8 > 3446.628)
drop if assets9 != . & (assets9 < 8.435747 | assets9 > 3767.053)

*===================*
* State-Owned Firms *
*===================*

* Use China (Original)

drop if sector < 1300 | sector >= 4400
drop if inflator*sales/1.19453 < 5000

gen stateownership = (capitalstate+capitalcollective)/capitalpaidin
drop if capitalstate+capitalcollective < 0 | capitalpaidin < 0

drop if totalassets <= 0 | fixedassets < 0 | capitalpaidin < 0
drop if capitalpaidin != capitalstate+capitalcollective+capitalcorporate+capitalpersonal+capitalhkmt+capitalforeign
drop if sales < 0 | costofgoodssold < 0
drop if currentassets < 0 | inventory < 0 | accountsreceivable < 0 | currentassets < inventory+accountsreceivable
drop if depreciationcurrentyear < 0 | totalliabilities < 0
drop if wagespayable < 0
drop if fixedassets > totalassets
drop if currentassets > totalassets

rename year fyear

*drop if stateownership >= 0.5
*drop if stateownership < 0.5

* Compute sector gamma

gen industry2 = floor(sector/100)

gen gamma_s = .
replace gamma_s = 1.490976654 if industry2 == 13
replace gamma_s = 1.392491385 if industry2 == 14
replace gamma_s = 1.472384382 if industry2 == 15
replace gamma_s = 2.054080711 if industry2 == 16
replace gamma_s = 1.526464682 if industry2 == 17
replace gamma_s = 1.486347602 if industry2 == 18
replace gamma_s = 1.526755886 if industry2 == 19
replace gamma_s = 1.608595255 if industry2 == 20
replace gamma_s = 1.407367625 if industry2 == 21
replace gamma_s = 1.506641667 if industry2 == 22
replace gamma_s = 1.530490644 if industry2 == 23
replace gamma_s = 1.476333073 if industry2 == 24
replace gamma_s = 1.453136808 if industry2 == 25
replace gamma_s = 1.517054308 if industry2 == 26
replace gamma_s = 1.444064939 if industry2 == 27
replace gamma_s = 1.473580833 if industry2 == 28
replace gamma_s = 1.478408919 if industry2 == 29
replace gamma_s = 1.503288204 if industry2 == 30
replace gamma_s = 1.510329972 if industry2 == 31
replace gamma_s = 1.441427227 if industry2 == 32
replace gamma_s = 1.459015071 if industry2 == 33
replace gamma_s = 1.531359615 if industry2 == 34
replace gamma_s = 1.515730229 if industry2 == 35
replace gamma_s = 1.441313663 if industry2 == 36
replace gamma_s = 1.555507141 if industry2 == 37
replace gamma_s = 1.973600948 if industry2 == 39
replace gamma_s = 1.549345766 if industry2 == 40
replace gamma_s = 1.478156505 if industry2 == 41
replace gamma_s = 1.378515606 if industry2 == 42
replace gamma_s = 1.480411856 if industry2 == 43
