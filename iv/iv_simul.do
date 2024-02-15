********************************************************************************
********************************************************************************
* Consistency of 2SLS


clear


********************************************************************************
********************************************************************************
* OVB Correlation Process 1

* Generate program
capture program drop mcsimul
program mcsimul, eclass
clear 
local par = 3.2

set obs 1000

gen z=rnormal()

matrix corruv=[1,.7\.7,1]
drawnorm u v , corr(corruv)

* Data Generation Process for X includes an error term "v"
gen x = 1 + 0.1*z + v
* DGP for Y includes error term "u"
gen y = 0.1 + `par' * x + u

* u and v are correlated
corr v u
 
reg y x
matrix b=  _b[x]

reg x z
matrix b=b,_b[z]

ivregress 2sls y (x=z)
matrix b=b,_b[x]
matrix colname b= ols fs iv

ereturn post b
end


simulate, reps(1000): mcsimul

twoway ( kdensity _b_ols ) ( kdensity _b_iv ) , xline(3.2) legend(order(1 "OLS" 2 "2SLS" ))



********************************************************************************
********************************************************************************
* OVB Correlation Process 2

clear

capture program drop mcsimul
program mcsimul, eclass
clear 
local par = 3.2

set obs 1000

gen z = rnormal()

gen u = rnormal(0,1)

* Data Generation Process for X includes an error term "v"
gen x = 1 + 0.1*z + 0.8*u

* DGP for Y includes error term "u"
gen y = 0.1 + `par' * x + u

* u and v are correlated
* corr v u
 
reg y x
matrix b=  _b[x]

reg x z
matrix b=b,_b[z]

ivregress 2sls y (x=z)
matrix b=b,_b[x]
matrix colname b= ols fs iv

ereturn post b
end


simulate, reps(1000): mcsimul

twoway ( kdensity _b_ols ) ( kdensity _b_iv ) , xline(3.2)
sum _b_iv
