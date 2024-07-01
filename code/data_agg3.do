// This code is used to aggregate the data at the 3-digit level of the SNI 2007 classification

// I need to redefine the variable agg_3_level to create a new variable that will be used to aggregate the data at the 3-digit level in order to meet the confidentiality requirements of the data

gen agg_3_level = sni_2007_3d
replace agg_3_level = 101.5 if sni_2007_3d == 101
replace agg_3_level = 101.5 if sni_2007_3d == 102
replace agg_3_level = 103.5 if sni_2007_3d == 103
replace agg_3_level = 103.5 if sni_2007_3d == 104
replace agg_3_level = 108.5 if sni_2007_3d == 108
replace agg_3_level = 108.5 if sni_2007_3d == 109
replace agg_3_level = 110.5 if sni_2007_3d == 110
replace agg_3_level = 110.5 if sni_2007_3d == 120
replace agg_3_level = 132.5 if sni_2007_3d == 131
replace agg_3_level = 132.5 if sni_2007_3d == 132
replace agg_3_level = 132.5 if sni_2007_3d == 133
replace agg_3_level = 141.5 if sni_2007_3d == 141
replace agg_3_level = 141.5 if sni_2007_3d == 143
replace agg_3_level = 151.5 if sni_2007_3d == 151
replace agg_3_level = 151.5 if sni_2007_3d == 152
replace agg_3_level = 191.5 if sni_2007_3d == 191
replace agg_3_level = 191.5 if sni_2007_3d == 192
replace agg_3_level = 202.5 if sni_2007_3d == 202
replace agg_3_level = 202.5 if sni_2007_3d == 203
replace agg_3_level = 202.5 if sni_2007_3d == 204
replace agg_3_level = 205.5 if sni_2007_3d == 205
replace agg_3_level = 205.5 if sni_2007_3d == 206
replace agg_3_level = 211.5 if sni_2007_3d == 211
replace agg_3_level = 211.5 if sni_2007_3d == 212
replace agg_3_level = 232.5 if sni_2007_3d == 232
replace agg_3_level = 232.5 if sni_2007_3d == 233
replace agg_3_level = 234.5 if sni_2007_3d == 234
replace agg_3_level = 234.5 if sni_2007_3d == 235
replace agg_3_level = 236.5 if sni_2007_3d == 236 //
replace agg_3_level = 236.5 if sni_2007_3d == 237 // This could be assigned to 239 as well
replace agg_3_level = 242.5 if sni_2007_3d == 242
replace agg_3_level = 242.5 if sni_2007_3d == 243
replace agg_3_level = 251.5 if sni_2007_3d == 251
replace agg_3_level = 251.5 if sni_2007_3d == 252
replace agg_3_level = 253.5 if sni_2007_3d == 253
replace agg_3_level = 253.5 if sni_2007_3d == 254
replace agg_3_level = 255.5 if sni_2007_3d == 255
replace agg_3_level = 255.5 if sni_2007_3d == 256
replace agg_3_level = 255.5 if sni_2007_3d == 257
replace agg_3_level = 261.5 if sni_2007_3d == 261
replace agg_3_level = 261.5 if sni_2007_3d == 262
replace agg_3_level = 263.5 if sni_2007_3d == 263
replace agg_3_level = 263.5 if sni_2007_3d == 264
replace agg_3_level = 265.5 if sni_2007_3d == 265
replace agg_3_level = 265.5 if sni_2007_3d == 266
replace agg_3_level = 267.5 if sni_2007_3d == 267
replace agg_3_level = 267.5 if sni_2007_3d == 268
replace agg_3_level = 271.5 if sni_2007_3d == 271
replace agg_3_level = 271.5 if sni_2007_3d == 272
replace agg_3_level = 273.5 if sni_2007_3d == 273
replace agg_3_level = 275.5 if sni_2007_3d == 275
replace agg_3_level = 275.5 if sni_2007_3d == 279
replace agg_3_level = 283.5 if sni_2007_3d == 283
replace agg_3_level = 283.5 if sni_2007_3d == 284
replace agg_3_level = 292.5 if sni_2007_3d == 291 //
replace agg_3_level = 292.5 if sni_2007_3d == 292 // This could be assigned to 293 as well
replace agg_3_level = 301.5 if sni_2007_3d == 301
replace agg_3_level = 301.5 if sni_2007_3d == 302
replace agg_3_level = 303.5 if sni_2007_3d == 303
replace agg_3_level = 303.5 if sni_2007_3d == 304
replace agg_3_level = 309.5 if sni_2007_3d == 309
replace agg_3_level = 309.5 if sni_2007_3d == 310
replace agg_3_level = 321.5 if sni_2007_3d == 321
replace agg_3_level = 321.5 if sni_2007_3d == 322
replace agg_3_level = 323.5 if sni_2007_3d == 323
replace agg_3_level = 323.5 if sni_2007_3d == 324
replace agg_3_level = 325.5 if sni_2007_3d == 325
replace agg_3_level = 325.5 if sni_2007_3d == 329
replace agg_3_level = 331.5 if sni_2007_3d == 331
replace agg_3_level = 331.5 if sni_2007_3d == 332

// Now we can test the number of observations in each level

tab agg_3_level


// Now We can use the agg_3_level to aggregate the data

collapse (sum) vlist , by(agg_3_level)