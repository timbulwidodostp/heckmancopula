*! version 1.0.0 T.Hasebe April 2012
* Given dependence parameter and copula, compute Kendall's tau
# delimit ;
program define _tau, rclass ;
	syntax anything(everything), COPula(string) [DISplay] ;

	if "`display'"!="" {;	// treat anything as scalar
		tempname theta tau;
		scalar `theta' = `anything';
		if "`copula'"=="product" {;scalar `tau' = 0; };
		else if "`copula'"=="gaussian" {;scalar `tau'=2*asin(`theta')/_pi;};
		else if "`copula'"=="amh" {;scalar `tau'=(3*`theta'-2)/(3*`theta')-(2/3)*(1-1/`theta')^2 *ln(1-`theta'); }; 
		else if "`copula'"=="fgm" {;scalar `tau'=2*`theta'/9 ;};
		else if "`copula'"=="clayton" {; scalar `tau' = `theta'/(`theta'+2); };
		else if "`copula'"=="gumbel" {; scalar `tau' = (`theta'-1)/`theta'; };	
		else if "`copula'"=="frank" | "`copula'"=="joe" {; mata: tau_inte_scalar("`theta'", "`copula'"); scalar `tau' = r(tau);
			*if "`copula'"=="joe" scalar `tau' = 0 if `tau'<0;
		};
		dis `tau'; return scalar tau = `tau';
	};	//end of if display != ""

	else {;
		tempvar theta; tokenize `anything'; qui gen double `theta' = `1';
		macro shift; tokenize `*';
		local tau `1'; macro shift; local rest `*';
					
		tau_var `tau' `rest', theta(`theta') copula(`copula');
	};
	
end;
 
program define tau_var ;
	syntax newvarlist(max=1) [if] [in], theta(varlist max=1) copula(string);
	
	*dis "here";
	*marksample touse;
	
	tempvar tau;
	if "`copula'"=="product" {; qui gen double `tau' = 0; };
	else if "`copula'"=="gaussian" {;qui gen double `tau'=2*asin(`theta')/_pi;};
	else if "`copula'"=="amh" {;qui gen double `tau'=(3*`theta'-2)/(3*`theta')-(2/3)*(1-1/`theta')^2 *ln(1-`theta'); }; 
	else if "`copula'"=="fgm" {;qui gen double `tau'=2*`theta'/9 ;};
	else if "`copula'"=="clayton" {; qui gen double `tau' = `theta'/(`theta'+2); };
	else if "`copula'"=="gumbel" {; qui gen double `tau' = (`theta'-1)/`theta'; };	
	else if "`copula'"=="frank" | "`copula'"=="joe" {;
		qui gen double `tau' = .; local vars `theta' `tau';
		mata: tau_inte_var("`copula'", "`vars'");
		if "`copula'"=="joe" replace `tau' = 0 if `tau' < 0;
	};
	
	gen `varlist' = `tau' ; 
	
end; 

#delimit cr
mata:
version 10.1
mata set matalnum on
mata set matastrict on
void tau_inte_scalar(string scalar theta_t, string scalar copula)	
{
	real matrix s, phi	
	real scalar theta, S, i, j, integ, tau
	
	theta = st_numscalar(theta_t)
	S = 10000
	s = J(1,S-1,0)
	for (i =1 ; i < S; i++) {
		s[.,i] = i/S 
	}		
	
	if (copula == "frank") {
		phi = (-ln((exp(-theta*s):-1):/(exp(-theta):-1))):*(1:-exp(theta*s))/theta
	}
	else if (copula == "joe") {
		phi = ((1:-(1:-s):^theta):*log(1:-(1:-s):^theta)):/(theta:*(1:-s):^(theta:-1))
	}
	integ = 0 
	for (i = 1; i < S-1; i++) {
		j = i + 1
		integ = integ + (phi[1,i]+phi[1,j]) 
	}
	tau = 1+4*integ/(2*(S-1))
	st_numscalar("r(tau)",tau) 
}


void tau_inte_var(string scalar copula, string scalar vars)
{
	real matrix s 
	real colvector integ, tau, theta
	real scalar  S, i, j, N
	
	st_view(X=.,.,tokens(vars))
	theta = X[.,1]
	N = rows(theta)
	S = 10000
	s = J(N,S-1,0)
	for (i =1 ; i < S; i++) {
		s[.,i] = i/S * J(N,1,1)
	}		
	
	if (copula == "frank") {
		phi = (-ln((exp(-theta:*s):-1):/(exp(-theta):-1))):*(1:-exp(theta:*s)):/theta
	}
	else if (copula == "joe") {
		phi = ((1:-(1:-s):^theta):*log(1:-(1:-s):^theta)):/(theta:*(1:-s):^(theta:-1))
	}
	
	integ = J(N,1,0)
	for (i = 1; i < S-1; i++) {
		j = i + 1
		integ[.,1] = integ + (phi[.,i]+phi[.,j]) 
	}
	tau = 1:+4*integ/(2*(S-1))
	X[.,2] = tau
	
}

mata set matalnum off
mata set matastrict off
end
