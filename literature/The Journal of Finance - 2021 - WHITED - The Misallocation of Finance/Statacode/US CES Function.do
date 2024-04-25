* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

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
*gen PF_si = inflator*(prcc_f*csho+D_si)
drop if PF_si <= 0 | PF_si == .

gen F_si = PF_si/piship

*-----*
* CES *
*-----*

gen f_i = log(F_si)

gen d_i = log(D_si)
gen e_i = log(E_si)

gen de_i = (d_i-e_i)^2

xtset gvkey fyear

xtreg f_i d_i e_i de_i, fe
test _b[d_i] + _b[e_i] = 1, coef

*gmm ((-f_i+{rho}*l.f_i)+{beta_A}*(1-{rho})+{beta_D}*(d_i-{rho}*l.d_i)+(1-{beta_D})*(e_i-{rho}*l.e_i)+{beta_DE}*(de_i-{rho}*l.de_i)), instruments(l.f_i l.d_i l.e_i l.de_i)

matrix key = J(20,1,.)
matrix ces = J(20,6,.)

local ii = 1

forvalues jj=20/39{
	quietly summ industry2 if industry2 == `jj'
	if r(min) != . {
		matrix key[`ii',1] = `jj'
		local ii = `ii'+1
	}
	display `ii'
}

forvalues jj=1/20{
	capture noisily xtreg f_i d_i e_i de_i if industry2 == key[`jj',1], fe
	capture noisily test _b[d_i] + _b[e_i] = 1, coef
	capture noisily matrix temp = J(1,5,.)
	capture noisily matrix temp = e(b)
	clear results
	capture noisily matrix ces[`jj',1] = key[`jj',1]
	capture noisily matrix ces[`jj',2] = temp[1,1]
	capture noisily matrix ces[`jj',3] = temp[1,2]
	capture noisily matrix ces[`jj',4] = temp[1,3]
	capture noisily matrix ces[`jj',5] = temp[1,4]
	capture noisily matrix ces[`jj',6] = temp[1,5]
}

matrix list ces
