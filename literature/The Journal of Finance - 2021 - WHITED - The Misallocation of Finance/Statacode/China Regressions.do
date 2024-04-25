* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

* Get r_si and lambda_si from "China Misallocation"

gen location = 0
replace location = 1 if zipcode >= 100000 & zipcode < 109999
replace location = 1 if zipcode >= 200000 & zipcode < 209999
replace location = 1 if zipcode >= 510000 & zipcode < 511499
replace location = 1 if zipcode >= 518000 & zipcode < 518999

gen state = 0
replace state = 1 if capitalstate > 0

gen foreign = 0
replace foreign = 1 if capitalhkmt+capitalforeign > 0

gen size = log(inflator*totalassets/1000)

gen time = fyear-1998

gen age = fyear-starttimeyear
drop if age < 0
replace age = 100 if age > 100
gen young = 0
replace young = 1 if age <= 3

reg r_si location state foreign size time young
reg lambda_si location state foreign size time young
