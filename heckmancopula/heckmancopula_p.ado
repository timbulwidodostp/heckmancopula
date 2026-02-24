*! version 1.2.0 August 7, 2012, 2012 
* Author Takuya Hasebe 
*! version 1.0.1 January 2012 

*** predict command for heckmancopula ***
#delimit ;

program define heckmancopula_p;
	version 10.1;
	syntax anything [if] [in] [, psel xb xbsel y_c0 y_c1 cll ];
	
	marksample touse;
	
	syntax newvarname [if] [in] [, psel xb xbsel y_c0 y_c1 cll];
	
	if "`e(cmd)'" != "heckmancopula" error 301;
	local check "`psel' `xb' `xbsel' `y_c0' `y_c1' `cll'";
	if wordcount("`check'")>1 {;
		dis as error "Only one statistic is allowed.";
	};
	if wordcount("`check'")==0 {;
		dis as text "Option psel is assumed";
		local psel "psel";
	};
	
	local y_s: word 1 of `e(depvar)';
	local y: word 2 of `e(depvar)';
	
	local margin  "`e(margin1)'";
	local copula  "`e(copula)'";
	local margins "`e(margsel)'";
	
	tempname beta betas beta1 lnsig1 sig1 atheta1 theta1 df1;
	
	matrix `beta' = e(b);
	
	matrix `betas' = `beta'[1,"select:"];
	matrix `beta1' = `beta'[1,"`y':"];
	scalar `lnsig1' = _b[/lnsigma];
	scalar `sig1' = exp(`lnsig1');
	if "`copula'"!="product" scalar `atheta1' = _b[/atheta];
	
	if "`margin'"=="t" {; 
		if "`e(df)'" == "" scalar `df1' = exp(_b[/lndf]);
		else scalar `df1' = `e(df)';
	};
	
	
	tempvar xbs xb1 prob_sel;
	
	local sign1 = 2*`e(negative)'-1;
	

	matrix score double `xbs' = `betas' if `touse';
	
	if ~missing("`xbsel'") {;
		generate `typlist' `varlist' = `xbs' if `touse';
		exit;
	};
	
	tempvar z Fs;
	qui gen double `z' = `xbs';
	if "`margins'" == "probit" {;
		qui gen double `Fs' = normal(-`z') if `touse';
	};
	else if "`margins'" == "logit" {;
		qui gen double `Fs' = 1/(1+exp(`z')) if `touse';
	};
	
	if ~missing("`psel'") {;
		generate `typlist' `varlist' = 1-`Fs' if `touse';
		exit ;
	};
	
	tempvar xb1 ;
	
	matrix score double `xb1' = `beta1' if `touse';
		
	if ~missing("`xb'") {;
		generate `typlist' `varlist' = `xb1' if `touse';
		exit ;
	};
	
	tempvar mill z_minus;
	qui gen double `z_minus' = - `xbs' if `touse';
		
	if "`margin'" != "t" {; scalar `df1' = 1;};
	
	if ~missing("`y_c0'") {;
		local truncation "above";
		if "`copula'" == "gaussian" {;
			scalar `theta1' = tanh(`atheta1');
			if "`margin'" == "normal" {;
				qui gen double `mill' = -`theta1'*normalden(invnormal(`Fs'))/(`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula'","`margin'","`truncation'", "`atheta1'", "`df1'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
			mata: _conditional("`Fs'","`mill'","`touse'","`copula'","`margin'","`truncation'", "`atheta1'", "`df1'") ;
		};
		generate `typlist' `varlist' = `xb1' + `sign1'*`sig1'*`mill' if `touse';
		exit;
	};
	if ~missing("`y_c1'") {;
		local truncation "below";
		if "`copula'" == "gaussian" {;
			scalar `theta1' = tanh(`atheta1');
			if "`margin'" == "normal" {;
				qui gen double `mill' = `theta1'*normalden(invnormal(`Fs'))/(1-`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula'","`margin'","`truncation'", "`atheta1'", "`df1'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
			mata: _conditional("`Fs'","`mill'","`touse'","`copula'","`margin'","`truncation'", "`atheta1'", "`df1'") ;
		};
		generate `typlist' `varlist' = `xb1' + `sign1'*`sig1'*`mill' if `touse';
		exit;
	};
	
	*** contribution to likelihood function *** ;
	
	if ~missing("`cll'") {;
	*** dependent variables *** ;
	tempvar S e1 z;
	
	capture sum `y_s'; 
	if _rc != 0 {; 
		qui gen `y_s' = (`y'!=.);
	};
	
	qui gen double `S' = `y_s';	//selection indicator
	*** Standarized continuous dependent variables *** ;
	qui gen double `e1' = `sign1'*(`y'-`xb1')/`sig1' ;
	
	*** Specify Margins ***;
	tempvar F1 f1 lnf1;
		
	* outcome *;
	if "`margin'"=="normal" {;
		qui gen double `F1' = normal(`e1');
		qui gen double `f1' = normalden(`e1');
		qui gen double `lnf1' = lnnormalden(`e1');
	};
	else if "`margin'" == "logistic" {;
		qui gen double `F1' = 1/(1+exp(-`e1'));
		qui gen double `f1' = `F1'*(1-`F1');
		qui gen double `lnf1' = -`e1'-2*ln(1+exp(-`e1'));
	};
	else if "`margin'" == "t" {;
		qui gen double `F1' = 1 - ttail(`df1',`e1');
		qui gen double `f1' = tden(`df1',`e1');
		qui	gen double `lnf1' = lngamma((`df1'+1)/2)-0.5*ln(`df1'*_pi)-lngamma(`df1'/2)
			-(`df1'+1)*ln(1+(`e1')^2/`df1')/2 ;
	};

	*** Copula-specific *** ;
	tempvar C1 ;
	
	*** Between Selection and Eq.1 *** ;
	if "`copula'" == "product" {;
		qui gen double `C1' = `Fs';
	};
	else if "`copula'" == "gaussian" {;
		tempvar vs v1 Z1;
		tempname sr_t1 ;
		scalar `theta1' = tanh(`atheta1');	//transform to dependent parameter
		scalar `sr_t1' = sqrt(1-(`theta1')^2);
		
		qui gen double `vs' = invnormal(`Fs');
		qui gen double `v1' = invnormal(`F1');
		qui gen double `Z1 ' = `vs'/`sr_t1' - `theta1'*`v1'/`sr_t1';
		qui gen double `C1' = normal(`Z1');
	};
	else if "`copula'" == "fgm" {;
		scalar `theta1' = tanh(`atheta1');
		qui gen double `C1' = `Fs'*(1+`theta1'*(1-`Fs')*(1-`F1')) - `theta1'*`F1'*`Fs'*(1-`Fs');
	};
	else if "`copula'" == "plackett" {;
		tempvar r root ;
		scalar `theta1' = exp(`atheta1');
		qui gen double `r' = 1+(`theta1'-1)*(`F1'+`Fs');
		qui gen double `root' = (`r')^2 - 4*`F1'*`Fs'*`theta1'*(`theta1'-1);
		qui gen double `C1' = 0.5 - 0.5*(`r'-2*`Fs'*`theta1')*(`root')^(-1/2);
		
	};
	else {;	//Acrhimedean copulas
		tempvar C phi1_F1 phi1_Fs phi1_C;
		if "`copula'" == "amh" {;
			scalar `theta1' = tanh(`atheta1');
			qui gen double `C' = `F1'*`Fs'/(1-`theta1'*(1-`Fs')*(1-`F1'));	
			qui gen double `phi1_F1' = `theta1'/(1-`theta1'*(1-`F1'))-1/`F1';
			qui gen double `phi1_Fs' = `theta1'/(1-`theta1'*(1-`Fs'))-1/`Fs';
			qui gen double `phi1_C' = `theta1'/(1-`theta1'*(1-`C'))-1/`C';
		};
		else if "`copula'"=="clayton" {;
			scalar `theta1' = exp(`atheta1');
			qui gen double `C' = (`Fs'^(-`theta1')+`F1'^(-`theta1')-1)^(-1/`theta1');
			qui gen double `phi1_F1' = -(`F1')^(-`theta1'-1);
			qui gen double `phi1_Fs' = -(`Fs')^(-`theta1'-1);
			qui gen double `phi1_C' =  -(`C')^(-`theta1'-1);
		};
		else if "`copula'"=="frank" {;
			scalar `theta1' = `atheta1';
			qui gen double `C' = -ln(1+(exp(-`theta1'*`Fs')-1)*(exp(-`theta1'*`F1')-1)/(exp(-`theta1')-1))/`theta1';
			qui gen double `phi1_F1' = `theta1'/(1-exp(`theta1'*`F1'));
			qui gen double `phi1_Fs' = `theta1'/(1-exp(`theta1'*`Fs'));
			qui gen double `phi1_C' = `theta1'/(1-exp(`theta1'*`C'));
		};
		else if "`copula'"=="gumbel" {;
			scalar `theta1' = 1 + exp(`atheta1');
			qui gen double `C' = exp(-((-ln(`F1'))^(`theta1')+(-ln(`Fs'))^(`theta1'))^(1/`theta1'));
			qui gen double `phi1_F1' = -`theta1'*(-ln(`F1'))^(`theta1'-1) /`F1';
			qui gen double `phi1_Fs' = -`theta1'*(-ln(`Fs'))^(`theta1'-1) /`Fs';
			qui gen double `phi1_C' = -`theta1'*(-ln(`C'))^(`theta1'-1) /`C';
		};
		else if "`copula'"=="joe" {;
			tempvar guard;
			scalar `theta1' = 1 + exp(`atheta1');
			qui gen double `C' = 1-((1-`Fs')^`theta1'+(1-`F1')^`theta1'-((1-`Fs')*(1-`F1'))^`theta1')^(1/`theta1');
			qui gen double `phi1_F1' = -`theta1'*(1-`F1')^(`theta1'-1)/(1-(1-`F1')^`theta1');
			qui gen double `phi1_Fs' = -`theta1'*(1-`Fs')^(`theta1'-1)/(1-(1-`Fs')^`theta1');
			qui gen double `phi1_C' = -`theta1'*(1-`C')^(`theta1'-1)/(1-(1-`C')^`theta1');
			qui gen double `guard' = 1-`Fs';
		};
		qui gen double `C1' = `phi1_F1'/`phi1_C';
	};
	qui replace `C1' = 0 if `Fs'==0 ;	// `F1'==0 | `F1'==1;	//?;
	qui replace `C1' = 1 if `Fs'==1 ;
	capture replace `C1' = 0 if `guard' == 1;
	
	*** Log-likelihood ***;
	* General * ;
	
	generate `typlist' `varlist' = cond(`S'==0, ln(`Fs'), ln(1-`C1') + `lnf1'-`lnsig1') if `touse' ;
	exit;
	}; //end of ~missing("`cll'");
	
end;
*include _conditional.mata;
*include num_inte.mata;
