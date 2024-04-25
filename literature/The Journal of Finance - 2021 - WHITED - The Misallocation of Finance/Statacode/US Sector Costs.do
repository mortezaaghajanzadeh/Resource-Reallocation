* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

gen industry3 = floor(sic/10)
gen industry2 = floor(sic/100)
gen industry1 = floor(sic/1000)

merge m:1 industry3 fyear using Beta

drop if _merge == 2

bysort fyear industry2: egen beta2 = mean(beta)
replace beta = beta2 if beta == .

*---------------*
* Misallocation *
*---------------*

scalar sigma = 1.77

gen D_si = inflator*lt
gen E_si = inflator*(lse-lt)
drop if D_si <= 0 | D_si == . | E_si <= 0 | E_si == .

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

*-----*
* r_s *
*-----*

bysort fyear industry3: egen interest = sum(xint)
bysort fyear industry3: egen debt = sum(dlc+dltt)

gen r_s = interest/debt
replace r_s = 0.03 if r_s == 0
replace r_s = 1 if r_s > 1

*----------*
* lambda_s *
*----------*

* Unlevered beta used
gen lambda_s = 0.03+0.06*(beta/(1+(1-0.35)*D_s/E_s))

*---------------*
* Misallocation *
*---------------*

gen alpha_s = r_s*D_s^(1/gamma_s)/(r_s*D_s^(1/gamma_s)+lambda_s*E_s^(1/gamma_s))

gen A_si = PF_si^(sigma/(sigma-1))/(alpha_s*D_si^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))

bysort fyear industry3: egen A_s_temp = sum(A_si^(sigma-1))
gen D_si_eff = D_s*A_si^(sigma-1)/A_s_temp
gen E_si_eff = E_s*A_si^(sigma-1)/A_s_temp

gen F_si = A_si*(alpha_s*D_si^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))
gen F_si_eff = A_si*(alpha_s*D_si_eff^((gamma_s-1)/gamma_s)+(1-alpha_s)*E_si_eff^((gamma_s-1)/gamma_s))^(gamma_s/(gamma_s-1))

bysort fyear industry3: egen F_s = sum(F_si^((sigma-1)/sigma))
replace F_s = F_s^(sigma/(sigma-1))
bysort fyear industry3: egen F_s_eff = sum(F_si_eff^((sigma-1)/sigma))
replace F_s_eff = F_s_eff^(sigma/(sigma-1))

bysort fyear: egen PF = sum(PF_si)
bysort fyear industry3: egen theta_s = sum(PF_si)
replace theta_s = theta_s/PF

gen F_temp = F_s^theta_s
gen F_eff_temp = F_s_eff^theta_s

gen r_si = alpha_s*(sigma-1)/sigma*PF_si/(alpha_s*D_si+(1-alpha_s)*E_si^((gamma_s-1)/gamma_s)*D_si^(1/gamma_s))
gen lambda_si = (1-alpha_s)*(sigma-1)/sigma*PF_si/(alpha_s*D_si^((gamma_s-1)/gamma_s)*E_si^(1/gamma_s)+(1-alpha_s)*E_si)

collapse (mean) F_temp F_eff_temp, by(fyear industry3)

bysort fyear: egen F = sum(log(F_temp))
bysort fyear: egen F_eff = sum(log(F_eff_temp))

replace F = exp(F)
replace F_eff = exp(F_eff)

gen ra_gain = F/F_eff
bysort fyear: summ ra_gain

collapse (mean) ra_gain, by(fyear)

summ ra_gain
