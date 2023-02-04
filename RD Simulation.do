	*****************************************
	*****************************************

	* Regression Discontinuity Design
	* This routine simulates a Data Generation Process
	* defined by a discontinuity (C) on the running variable (X)
	* the effect delta applies around a bandwidth where units are not affected by OVB.

clear
local Simulations 10000
local N 10000
set seed 1

	* cutoff
	local C = 40
	* bandwidth
	local B = 10
	* Effect
	local delta = 1.5
	

* Base File with placebo Tests
	clear
	gen beta_ols = .
	gen beta_rd = .
	gen left_rd = .
	gen left_ols = .
	gen right_rd = .
	gen right_ols = .
	gen beta_iv = .
	gen left_iv = .
	gen right_iv = .
	gen frd = .
	tempfile Base
	save `Base' , replace

forvalues s=1(1)`Simulations' {
	di "Simulation `s'"

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
	gen y = 4 + `delta'*Treat + Mu if abs(distance)<(`B'+1)
	replace y = 4 + `delta'*Treat + 0.1*Unobs + Mu if abs(distance)>`B'
	
	*****************************************
	*****************************************
	* Estimators Test
	
	qui reg y Treat
	local beta_ols =  _b[Treat]
	local left_ols  = _b[Treat] - _se[Treat]*1.96
	local right_ols  = _b[Treat] + _se[Treat]*1.96

	qui ivreg y (Treat = cutoff)
	local beta_iv =  _b[Treat]
	local left_iv = _b[Treat] - _se[Treat]*1.96
	local right_iv = _b[Treat] + _se[Treat]*1.96

	qui  rdrobust y X  , c( `C' ) p(1) fuzzy(Treat) kernel(uni) h(`B')
	local beta_rd =  e(tau_cl)
	local left_rd  = e(ci_l_cl)
	local right_rd  = e(ci_r_cl)
	
	clear
	set obs 1

	gen beta_ols = `beta_ols'
	gen left_ols = `left_ols'
	gen right_ols = `right_ols'
	
	gen beta_iv = `beta_iv'
	gen left_iv = `left_iv'
	gen right_iv = `right_iv'
	
	gen beta_rd = `beta_rd'
	gen left_rd = `left_rd'
	gen right_rd = `right_rd'

	gen simulation = `s'
	gen frd = 1
	
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
	* Estimators Test
	
	qui reg y Treat
	local beta_ols =  _b[Treat]
	local left_ols  = _b[Treat] - _se[Treat]*1.96
	local right_ols  = _b[Treat] + _se[Treat]*1.96

	qui ivreg y (Treat = cutoff)
	local beta_iv =  _b[Treat]
	local left_iv = _b[Treat] - _se[Treat]*1.96
	local right_iv = _b[Treat] + _se[Treat]*1.96

	qui  rdrobust y X  , c( `C' ) p(1) kernel(uni) h(`B')
	local beta_rd =  e(tau_cl)
	local left_rd  = e(ci_l_cl)
	local right_rd  = e(ci_r_cl)
	
	
	clear
	set obs 1

	gen beta_ols = `beta_ols'
	gen left_ols = `left_ols'
	gen right_ols = `right_ols'
	
	gen beta_iv = `beta_iv'
	gen left_iv = `left_iv'
	gen right_iv = `right_iv'

	gen beta_rd = `beta_rd'
	gen left_rd = `left_rd'
	gen right_rd = `right_rd'

	gen simulation = `s'
	gen frd = 0
	
	tempfile S
	save `S' , replace
	
	use `Base' , replace
	append using `S'
	tempfile Base
	save `Base' , replace
	
}

clear
use `Base' , replace
net install grc1leg.pkg

local delta = 1.5
local con if frd==0
local le  xtitle("Coefficient", size(small))  legend(order(1 "OLS" 2 "IV" 3 "RD")) xline(`delta', lpattern(solid) lcolor(red) ) graphregion(color(white)) scheme(s2mono)
twoway (kdensity beta_ols `con', lpattern(solid)) (kdensity beta_iv `con', lpattern(dot))  (kdensity beta_rd `con' , lpattern(dash)) , title("Shard DGP") name(srd, replace) `le'

local con if frd==1
twoway (kdensity beta_ols `con', lpattern(solid)) (kdensity beta_iv `con', lpattern(dot))  (kdensity beta_rd `con' , lpattern(dash)) , title("Fuzzy DGP") name(frd, replace) `le'

grc1leg srd frd, graphregion(color(white)) scheme(s2mono) l1("Density", size(small))
gr_edit .legend.Edit, style(cols(3)) style(rows(0)) keepstyles

