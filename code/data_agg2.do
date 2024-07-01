// This code is used to aggregate the data at the 2-digit level of the SNI 2007 classification

// I need to redefine the variable agg_2_level to create a new variable that will be used to aggregate the data at the 2-digit level in order to meet the confidentiality requirements of the data

gen agg_2_level = sni_2007_2d
replace agg_2_level = 11.5 if sni_2007_2d == 11
replace agg_2_level = 11.5 if sni_2007_2d == 12
replace agg_2_level = 14.5 if sni_2007_2d == 14
replace agg_2_level = 14.5 if sni_2007_2d == 15
replace agg_2_level = 18.5 if sni_2007_2d == 18
replace agg_2_level = 18.5 if sni_2007_2d == 19
replace agg_2_level = 21.5 if sni_2007_2d == 21
replace agg_2_level = 21.5 if sni_2007_2d == 22
replace agg_2_level = 26.5 if sni_2007_2d == 26
replace agg_2_level = 26.5 if sni_2007_2d == 27
replace agg_2_level = 30.5 if sni_2007_2d == 30
replace agg_2_level = 30.5 if sni_2007_2d == 31
replace agg_2_level = 32.5 if sni_2007_2d == 32
replace agg_2_level = 32.5 if sni_2007_2d == 33

// Now we can test the number of observations in each level

tab agg_2_level


// Now We can use the agg_2_level to aggregate the data

collapse (sum) vlist , by(agg_2_level)

