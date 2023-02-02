clear
local Simulations 5000
local N 4000
set seed 100

	* Placebo Cutoff Range
	local placebo_min = 20
	local placebo_max = 60
	* cutoff
	local C = 40
	* bandwidth
	local B = 10
	* Effect
	local delta = 1.5
	

* Base File with placebo Tests
	clear
	gen beta=.
	gen left=.
	gen right=.
	gen cutoff=.
	tempfile Base
	save `Base' , replace

forvalues s=1(1)`Simulations' {
		
	clear
	set obs `N'

	*****************************************
	*****************************************
	* Fuzzy DGP

	cap drop X
	gen X = round(runiform(10,80))

	cap drop Assignment
	gen Assignment = rnormal(0,1)

	cap drop Treat
	gen Treat = 0
	replace Treat = 1 if Assignment>1.6 & X<`C'
	replace Treat = 1 if Assignment>0.7 & X>=`C'

	cap drop cutoff
	gen cutoff= X>=`C'
	bysort cutoff: sum Treat

	* Error term
	cap drop Mu
	gen Mu = rnormal(0,1)

	* Unobservables and Treatment are correlated.
	cap drop Unobs
	gen Unobs = rnormal(0,1) + 2.1*Treat
	corr Unobs Treat

	cap drop distance
	gen distance = X - `C'

	* Variables near the cutoff are identical except in the treatment
	cap drop y
	gen y = 4 + `delta'Treat + Mu if abs(distance)<(`B'+1)
	replace y = 4 + `delta'*Treat + 0.1*Unobs + Mu if abs(distance)>`B'
	
	*****************************************
	*****************************************
	* Placebo Test

	forvalues i=`placebo_min'(1)`placebo_max'{
		di "cutoff: `i'"
		qui rdrobust y X  , c( `i' ) p(1) fuzzy(Treat) kernel(uni) h(`B')
		
		local beta_`i' =  e(tau_cl)
		local left_`i'  = e(ci_l_cl)
		local right_`i'  = e(ci_r_cl)

	}

	clear
	set obs 600

	gen beta=.
	gen left=.
	gen right=.
	gen cutoff=.

	forvalues i = `placebo_min'(1)`placebo_max' {
		replace cutoff = `i' if (_n == `i')
		replace beta = `beta_`i'' if (_n == `i')
		replace left = `left_`i'' if (_n== `i')
		replace right = `right_`i'' if (_n== `i')
	}
	
	gen Simulation = `s'
	gen FRD = 1

	drop if cutoff==.

	tempfile S
	save `S' , replace
	
	use `Base' , replace
	append using `S'
	tempfile Base
	save `Base' , replace
	
	
	*****************************************
	*****************************************
	* Sharp DGP

	clear
	set obs `N'

	cap drop X
	gen X = round(runiform(10,80))

	cap drop Treat
	gen Treat = (X>=`C')
	cap drop cutoff
	gen cutoff= X>=`C'
	bysort cutoff: sum Treat

	* Error term
	cap drop Mu
	gen Mu = rnormal(0,1)

	* Unobservables and Treatment are correlated.
	cap drop Unobs
	gen Unobs = rnormal(0,1) + 2.1*Treat
	corr Unobs Treat

	cap drop distance
	gen distance = X - `C'

	* Variables near the cutoff are identical except in the treatment
	cap drop y
	gen y = 4 + `delta'*Treat + Mu if abs(distance)<(`B'+1)
	replace y = 4 + `delta'*Treat + 0.1*Unobs + Mu if abs(distance)>`B'

	*****************************************
	*****************************************
	* Placebo Test

	forvalues i=`placebo_min'(1)`placebo_max'{
		di "cutoff: `i'"
		qui rdrobust y X  , c( `i' ) p(1) fuzzy(Treat) kernel(uni) h(`B')
		
		local beta_`i' =  e(tau_cl)
		local left_`i'  = e(ci_l_cl)
		local right_`i'  = e(ci_r_cl)

	}

	clear
	set obs `N'

	gen beta=.
	gen left=.
	gen right=.
	gen cutoff=.

	forvalues i = `placebo_min'(1)`placebo_max' {
		replace cutoff = `i' if (_n == `i')
		replace beta = `beta_`i'' if (_n == `i')
		replace left = `left_`i'' if (_n== `i')
		replace right = `right_`i'' if (_n== `i')
	}
	
	gen Simulation = `s'
	gen FRD = 0

	drop if cutoff==.


	tempfile S
	save `S' , replace
	
	use `Base' , replace
	append using `S'
	tempfile Base
	save `Base' , replace
	
}

clear
use `Base' , replace

	
	gen significant = (left>0) | (right<0)
	sum significant


preserve
keep if FRD==1
drop if abs(right)>100 | abs(left)>100
twoway (scatter beta cutoff if significant==0, msymbol(diamond) mcolor(black) ) (scatter beta cutoff if significant==1, msymbol(diamond) mcolor(red) ) (rcap right left cutoff ), yline(0 , lpattern(solid) ) xline(`C', lpattern(dash) lcolor(black)) graphregion(color(white)) scheme(s2mono) leg(off) name(FRD, replace) title("Fuzzy RD")
restore

preserve
keep if FRD==0
drop if abs(right)>100 | abs(left)>100
twoway (scatter beta cutoff if significant==0, msymbol(diamond) mcolor(black) ) (scatter beta cutoff if significant==1, msymbol(diamond) mcolor(red) ) (rcap right left cutoff ), yline(0 , lpattern(solid) ) xline(`C', lpattern(dash) lcolor(black)) graphregion(color(white)) scheme(s2mono) leg(off) name(SRD, replace) title("Sharp RD")
restore


graph combine SRD FRD , graphregion(color(white))
