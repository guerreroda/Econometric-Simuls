********************************************************************************
********************************************************************************
* Difference in Difference Consistency


set seed 100

* Define number of Simulations
local reps = 100

* Define sample size
local N = 500

* Define treatment effects
local Beta = 1.5 


********************************************************************************
********************************************************************************
* MC simulations

forvalues i=1(1)`reps' {
	
	* Generate 
	clear
	set obs `N'
	
	gen i = _n
	gen t = 0
	
	* Append new time period to each individual
	preserve
	clear
	set obs `N'

	gen i = _n
	gen t = 1
	tempfile batch
	save `batch'
	restore
	append using `batch'
	
	
	* Define treatment group
	egen treat = max(round(runiform(0,1))), by(i)
	
	* DGP: Unobservable Error Term
	gen Epsilon = rnormal(0,1)
	
	* DGP: Outcome
	gen Y = 10 + `Beta' * t*treat + 1.1*t + Epsilon

	* Estimate DID ATE & save estimate
	reg Y treat##t
	local b1_r`i' = _b[1.treat#1.t]
	
	* Drop 50% of the sample.
	gen _r_ = rnormal(0,1)
	qui sum _r_ 
	drop if _r_ < `r(mean)'

	reg Y treat##t
	local b1_h`i' = _b[1.treat#1.t]
	
}

* Consistency
clear 
set obs `reps'
gen b1_r = .
gen b1_h = .
forvalues i=1(1)`reps' {

	replace b1_r = `b1_r`i'' if _n == `i'
	replace b1_h = `b1_h`i'' if _n == `i'
	
}

twoway (kdensity b1_r, color(blue%30)) ///        
       (kdensity b1_h, color(green%30)), ///   
       legend(order(1 "Full Sample" 2 "Half Sample" ) ring(0) pos(1)) xline( `Beta' , lcolor(red))
