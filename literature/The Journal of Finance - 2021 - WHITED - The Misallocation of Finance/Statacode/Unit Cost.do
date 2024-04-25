* Toni Whited and Jake Zhao
* University of Michigan and Peking University HSBC Business School
* Spring 2020

* Get r_si and lambda_si from "US/China Misallocation"

gen ratio_dev = D_si/(D_si+E_si)-D_s/(D_s+E_s)
gen unit_cost = D_si/(D_si+E_si)*r_si+E_si/(D_si+E_si)*lambda_si

scalar num_firm_sizes = 8
matrix cost_firm_size = J(2*num_firm_sizes,2,.)

forvalues ii=1/`=scalar(num_firm_sizes)'{

summ unit_cost if firm_size == `ii' & ratio_dev > 0, detail
matrix cost_firm_size[2*`ii'-1,1] = r(p50)

bootstrap r(p50), reps(100): summ unit_cost if firm_size == `ii' & ratio_dev > 0, detail
matrix cost_firm_size[2*`ii',1] = e(se)

summ unit_cost if firm_size == `ii' & ratio_dev < 0, detail
matrix cost_firm_size[2*`ii'-1,2] = r(p50)

bootstrap r(p50), reps(100): summ unit_cost if firm_size == `ii' & ratio_dev < 0, detail
matrix cost_firm_size[2*`ii',2] = e(se)

}

matrix list cost_firm_size
