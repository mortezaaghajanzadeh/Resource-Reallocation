* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

gen D_si = inflator*totalliabilities/1000
gen E_si = inflator*totalequity/1000
drop if D_si <= 0 | D_si == . | E_si <= 0 | E_si == .

gen PF_si = inflator*(totalprofit+vatpayable+wagespayable+depreciationcurrentyear)/1000
drop if PF_si <= 0 | PF_si == .

* gen F_si = PF_si
gen piship = industrialoutputcurrentprice/industrialoutputconstantprice
gen F_si = PF_si/piship

gen industry3 = floor(sector/10)
gen industry2 = floor(sector/100)
gen industry1 = floor(sector/1000)

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

*gmm ((-f_i+{rho}*l.f_i)+{beta_A}*(1-{rho})+exp({beta_D})*(d_i-{rho}*l.d_i)+(1-exp({beta_D}))*(e_i-{rho}*l.e_i)+{beta_DE}*(de_i-{rho}*l.de_i)), instruments(l.f_i l.d_i l.e_i l.de_i l2.f_i l2.d_i l2.e_i l2.de_i)

matrix key = J(31,1,.)
matrix ces = J(31,6,.)

local ii = 1

forvalues jj=13/43{
	quietly summ industry2 if industry2 == `jj'
	if r(min) != . {
		matrix key[`ii',1] = `jj'
		local ii = `ii'+1
	}
	display `ii'
}

forvalues jj=1/31{
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
