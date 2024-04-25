* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

scalar num_reps = 100
scalar num_years = 9

matrix fraction = J(num_reps,num_years,.)
matrix percent = J(num_reps,num_years,.)

forvalues ii=1/`=num_reps'{

use "US Compustat Sample Data.dta"

bsample

scalar sigma = 1.77
*replace gamma_s = 1000

gen D_si = inflator*lt
*gen D_si = inflator*(lt-act)
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
*gen PF_si = inflator*(prcc_f*csho+D_si)
drop if PF_si <= 0 | PF_si == .

bysort fyear industry3: egen D_s = sum(D_si)
bysort fyear industry3: egen E_s = sum(E_si)

gen alpha_s = D_s^(1/gamma_s)/(D_s^(1/gamma_s)+E_s^(1/gamma_s))

gen A_si = PF_si^(sigma/(sigma-1))/(alpha_s*D_si^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))
*gen A_si = PF_si^(sigma/(sigma-1))/(D_si^alpha_s*E_si^(1-alpha_s))

bysort fyear industry3: egen A_s_temp = sum(A_si^(sigma-1))
gen D_si_eff = D_s*A_si^(sigma-1)/A_s_temp
gen E_si_eff = E_s*A_si^(sigma-1)/A_s_temp

gen F_si = A_si*(alpha_s*D_si^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))
gen F_si_eff = A_si*(alpha_s*D_si_eff^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si_eff^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))
*gen F_si = A_si*D_si^alpha_s*E_si^(1-alpha_s)
*gen F_si_eff = A_si*D_si_eff^alpha_s*E_si_eff^(1-alpha_s)

bysort fyear industry3: egen F_s = sum(F_si^((sigma-1)/sigma))
replace F_s = F_s^(sigma/(sigma-1))
bysort fyear industry3: egen F_s_eff = sum(F_si_eff^((sigma-1)/sigma))
replace F_s_eff = F_s_eff^(sigma/(sigma-1))

bysort fyear: egen PF = sum(PF_si)
bysort fyear industry3: egen theta_s = sum(PF_si)
replace theta_s = theta_s/PF

gen F_temp = F_s^theta_s
gen F_eff_temp = F_s_eff^theta_s

collapse (mean) F_temp F_eff_temp, by(fyear industry3)

bysort fyear: egen F = sum(log(F_temp))
bysort fyear: egen F_eff = sum(log(F_eff_temp))

replace F = exp(F)
replace F_eff = exp(F_eff)

gen ra_gain = F/F_eff
gen ra_gain_percent = 100/ra_gain-100

forvalues jj=1/`=num_years'{
	summ ra_gain if fyear == `jj'+1998
	matrix fraction[`ii',`jj'] = r(mean)
	summ ra_gain_percent if fyear == `jj'+1998
	matrix percent[`ii',`jj'] = r(mean)
}

*summ ra_gain
*matrix fraction[`ii',1] = r(mean)
*summ ra_gain_percent
*matrix percent[`ii',1] = r(mean)

clear

}

matrix list fraction
matrix list percent
