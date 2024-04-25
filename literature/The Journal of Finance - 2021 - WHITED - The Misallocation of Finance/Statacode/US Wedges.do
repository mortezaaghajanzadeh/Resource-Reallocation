* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

* Get r_si and lambda_si from "US Misallocation"

* Year

scalar num_years = 9
matrix cost_year = J(2*num_years,4,.)

forvalues ii=1/`=scalar(num_years)'{

summ r_si if fyear == `ii'+1998
matrix cost_year[2*`ii'-1,1] = r(mean)

bootstrap r(mean), reps(100): summ r_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii',1] = e(se)

summ lambda_si if fyear == `ii'+1998
matrix cost_year[2*`ii'-1,2] = r(mean)

bootstrap r(mean), reps(100): summ lambda_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii',2] = e(se)

summ r_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii'-1,3] = r(p50)

bootstrap r(p50), reps(100): summ r_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii',3] = e(se)

summ lambda_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii'-1,4] = r(p50)

bootstrap r(p50), reps(100): summ lambda_si if fyear == `ii'+1998, detail
matrix cost_year[2*`ii',4] = e(se)

}

* Firm size

scalar num_firm_sizes = 8
matrix cost_firm_size = J(2*num_firm_sizes,4,.)

forvalues ii=1/`=scalar(num_firm_sizes)'{

summ r_si if firm_size == `ii'
matrix cost_firm_size[2*`ii'-1,1] = r(mean)

bootstrap r(mean), reps(100): summ r_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii',1] = e(se)

summ lambda_si if firm_size == `ii'
matrix cost_firm_size[2*`ii'-1,2] = r(mean)

bootstrap r(mean), reps(100): summ lambda_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii',2] = e(se)

summ r_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii'-1,3] = r(p50)

bootstrap r(p50), reps(100): summ r_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii',3] = e(se)

summ lambda_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii'-1,4] = r(p50)

bootstrap r(p50), reps(100): summ lambda_si if firm_size == `ii', detail
matrix cost_firm_size[2*`ii',4] = e(se)

}

* Print out results

matrix list cost_year
matrix list cost_firm_size
