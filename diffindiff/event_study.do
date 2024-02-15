* EVENT STUDY SIMULATION

clear
set seed 10

gen id = .
gen year = .
gen treatment = .


* GENERATE `size' NUMBER OF INDIVIDUALS
local size = 100
forvalues i = 1(1)`size' {	
	preserve
	clear
	set obs 10
	gen id = `i'
	gen year = 2000 + _n 
	gen treatment = runiform(0,1)
	qui sum treatment
	replace treatment = r(mean)>0.5
	
	tempfile Individual
	save `Individual' , replace
	restore
	
	append using `Individual'
	
}

keep id year treatment 
* GENERATE OUTCOME VARIABLE + RANDOM DISTURBANCE
* Treatment Effects = 2.5

gen mu = rnormal(0,1)
gen Y = 10 + 2.5 * treatment * (year>2005) + mu
* NOTE: DYNAMIC Treatment: You can add heterogeneous treatment effects using:
* gen Y = 10 + 2.5 * treatment * (year>2005) * ((year - 2005)*0.2) + mu

drop mu 

sum Y 
bysort treatment : sum Y 
sum Y if year<2006
sum Y if year>2005


* Simple Difference in Difference
cap gen _Post = year>2005
reg Y treatment##_Post , r


* Event Study: Method 1
eststo event : reg Y ib2005.year##treat i.id , r

forvalues x = 2001(1)2010 {
	local mylabel = "`mylabel' " + "`x'.year#1.treatment = `x'"
}
coefplot event , keep( *.year#1.treatment ) vertical omitted baselevels label yline(0) name(event1, replace) coeflabels(  `mylabel' )


* Event Study: Method 2
cap drop _lag* _lead* ref 
levelsof year , local(Years)
local c = 0 
local p = 0 
foreach Y of local Years {
	if `Y' < 2005 {
		gen _lag`c' = year==`Y' * treat
		label variable _lag`c' "`Y'"
		local c = `c' + 1
	}
	if `Y' > 2005 {
		gen _lead`p' = year==`Y' * treat
		label variable _lead`p' "`Y'"
		local p = `p' + 1
	}

}

gen ref = 1
label variable ref "2005"

eststo event2 : reg Y _lag* ref _lead* i.id i.year , r
coefplot event2 , keep( _lag* ref _lead* ) vertical omitted baselevels label yline(0) name(event2, replace)

graph combine event1 event2 , ycomm xcomm
