*! version 1.2.0 August 7, 2012 Takuya Hasebe

* version 1.1.0 June 7, 2012
* version 1.0.1 January 2012 

*** predict command for switchcopula ***
#delimit ;

program define switchcopula_p;
	syntax anything [if] [in] [, psel xb0 xb1 xbsel y0_c0 y0_c1 y1_c0 y1_c1 cll ];
	
	marksample touse;
	
	syntax newvarname [if] [in] [, psel xb0 xb1 xbsel y0_c0 y0_c1 y1_c0 y1_c1 cll];
	
	if "`e(cmd)'" != "switchcopula" error 301;
	local check "`psel' `xb0' `xb1' `xbsel' `y0_c0' `y0_c1' `y1_c0' `y1_c1' `cll'";
	if wordcount("`check'")>1 {;
		dis as error "Only one statistic is allowed.";
	};
	if wordcount("`check'")==0 {;
		dis as text "Option psel is assumed";
		local psel "psel";
	};
	
	local y_s: word 1 of `e(depvar)';
	local y_0: word 2 of `e(depvar)';
	local y_1: word 3 of `e(depvar)';
	
	local margin1  "`e(margin1)'";
	local margin0  "`e(margin0)'";
	local margins "`e(margsel)'";
	local copula0  "`e(copula0)'";
	local copula1  "`e(copula1)'";
	
	tempname beta beta_s beta0 beta1 lnsig0 lnsig1 sig0 sig1 atheta0 atheta1 theta0 theta1 df0 df1;
	
	matrix `beta' = e(b);
	
	matrix `beta_s' = `beta'[1,"select:"];
	matrix `beta0' = `beta'[1,"regime0:"];
	matrix `beta1' = `beta'[1,"regime1:"];
	scalar `lnsig0' = _b[/lnsigma0];
	scalar `lnsig1' = _b[/lnsigma1];
	scalar `sig0' = exp(`lnsig0');
	scalar `sig1' = exp(`lnsig1');
	if "`copula0'"!="product" scalar `atheta0' = _b[/atheta0];
	if "`copula1'"!="product" scalar `atheta1' = _b[/atheta1];
	
	if "`margin0'"=="t" {; 
		if "`e(df0)'" == "" scalar `df0' = exp(_b[/lndf0]);
		else scalar `df0' = `e(df0)';
	};
	if "`margin1'"=="t" {; 
		if "`e(df1)'" == "" scalar `df1' = exp(_b[/lndf1]);
		else scalar `df1' = `e(df1)';
	};
	
	tempvar xb_s x_b0 x_b1 prob_sel;
	
	local sign0 = 2*`e(negative0)'-1 ;
	local sign1 = 2*`e(negative1)'-1 ;
	
	matrix score double `xb_s' = `beta_s' if `touse';
	
	if ~missing("`xbsel'") {;
		generate `typlist' `varlist' = `xb_s' if `touse';
		exit;
	};
	
	tempvar z Fs;
	qui gen double `z' = `xb_s';
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
	
	tempvar xb0_t xb1_t ;
	
	matrix score double `xb0_t' = `beta0' if `touse';
	matrix score double `xb1_t' = `beta1' if `touse';
	
	if ~missing("`xb0'") {;
		generate `typlist' `varlist' = `xb0_t' if `touse';
		exit ;
	};

	if ~missing("`xb1'") {;
		generate `typlist' `varlist' = `xb1_t' if `touse';
		exit ;
	};
	
	tempvar mill z_minus;
	qui gen double `z_minus' = - `xb_s' if `touse';
	

	if "`margin0'" != "t" {; scalar `df0' = 1;}; 
	if "`margin1'" != "t" {; scalar `df1' = 1;};
	
	if ~missing("`y0_c0'") {;
		local truncation "above";
		if "`copula0'" == "gaussian" {;
			scalar `theta0' = tanh(`atheta0');
			if "`margin0'" == "normal" {;
				qui gen double `mill' = -`theta0'*normalden(invnormal(`Fs'))/(`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula0'","`margin0'","`truncation'", "`atheta0'", "`df0'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
			mata: _conditional("`Fs'","`mill'","`touse'","`copula0'","`margin0'","`truncation'", "`atheta0'", "`df0'") ;
		};
		generate `typlist' `varlist' = `xb0_t' + `sign0'*`sig0'*`mill' if `touse';
		exit;
	};
	if ~missing("`y0_c1'") {;
		local truncation "below";
		if "`copula0'" == "gaussian" {;
			scalar `theta0' = tanh(`atheta0');
			if "`margin0'" == "normal" {;
				qui gen double `mill' = `theta0'*normalden(invnormal(`Fs'))/(1-`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula0'","`margin0'","`truncation'", "`atheta0'", "`df0'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula0'","`margin0'","`truncation'", "`atheta0'", "`df0'") ;
		};
		generate `typlist' `varlist' = `xb0_t' + `sign0'*`sig0'*`mill' if `touse';
		exit;
	};
	if ~missing("`y1_c0'") {;
		local truncation "above";
		if "`copula1'" == "gaussian" {;
			scalar `theta1' = tanh(`atheta1');
			if "`margin1'" == "normal" {;
				qui gen double `mill' = -`theta1'*normalden(invnormal(`Fs'))/(`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula1'","`margin1'","`truncation'", "`atheta1'", "`df1'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
			mata: _conditional("`Fs'","`mill'","`touse'","`copula1'","`margin1'","`truncation'", "`atheta1'", "`df1'") ;
		};
		generate `typlist' `varlist' = `xb1_t' + `sign1'*`sig1'*`mill' if `touse';
		exit;
	};
	if ~missing("`y1_c1'") {;
		local truncation "below";
		if "`copula1'" == "gaussian" {;
			scalar `theta1' = tanh(`atheta1');
			if "`margin1'" == "normal" {;
				qui gen double `mill' = `theta1'*normalden(invnormal(`Fs'))/(1-`Fs') if `touse';
			};
			else {;
				qui gen double `mill' = . if `touse';
				mata: _conditional("`Fs'","`mill'","`touse'","`copula1'","`margin1'","`truncation'", "`atheta1'", "`df1'") ;
			};
		};
		else {;
			qui gen double `mill' = . if `touse';
			mata: _conditional("`Fs'","`mill'","`touse'","`copula1'","`margin1'","`truncation'", "`atheta1'", "`df1'") ;
		};
		generate `typlist' `varlist' = `xb1_t' + `sign1'*`sig1'*`mill' if `touse';
		exit;
	};
	
	
	*** contribution to likelihood function *** ;
	
	if ~missing("`cll'") {;

	
	*** dependent variables *** ;
	tempvar S e0 e1 ;
	qui gen double `S' = `y_s';	//selection indicator
		
	*** Standarized continuous dependent variables *** ;
	qui gen double `e0' = `sign0'*(`y_0'-`xb0_t')/`sig0' ;
	qui gen double `e1' = `sign1'*(`y_1'-`xb1_t')/`sig1' ;
	
	*** Specify Margins ***;
	tempvar F0 f0 lnf0 F1 f1 lnf1;
	* outcomes * ;
	
	forvalue j = 0/1 {;
		if "`margin`j''" == "normal" {;
			qui gen double `F`j'' = normal(`e`j'');
			qui gen double `f`j'' = normalden(`e`j'');
			qui gen double `lnf`j'' = lnnormalden(`e`j'');
		};
		else if "`margin`j''" == "logistic" {;
			qui gen double `F`j'' = 1/(1+exp(-`e`j''));
			qui gen double `f`j'' = `F`j''*(1-`F`j'');
			qui gen double `lnf`j'' = -`e`j''-2*ln(1+exp(-`e`j''));
		};
		else if "`margin`j''" == "t" {;
			qui gen double `F`j'' = 1 - ttail(`df`j'',`e`j'');
			qui gen double `f`j'' = tden(`df`j'',`e`j'');
			qui	gen double `lnf`j'' = lngamma((`df`j''+1)/2)-0.5*ln(`df`j''*_pi)-lngamma(`df`j''/2)
				-(`df`j''+1)*ln(1+(`e`j'')^2/`df`j'')/2 ;
		};
	};
	
	*** Copula-specific *** ;
	tempvar CS0 CS1 dCS0 dCS1 p1_F0 p1_CS0 p1_F1 p1_CS1;
	
	
	*** Copula-specific *** ;
	tempvar C1 C0 ;
	
	forvalue j = 0/1 {;
		if "`copula`j''" == "product" {;
			qui gen double `C`j'' = `Fs';
		};
		else if "`copula`j''" == "gaussian" {;
			tempvar v`j' Z`j';
			tempname sr_t`j' ;
			
			if (`j'==0 | "`copula0'"!="gaussian") {;
				tempvar vs;
				qui gen double `vs' = invnormal(`Fs');
			};
			
			scalar `theta`j'' = tanh(`atheta`j'');	//transform to dependent parameter
			scalar `sr_t`j'' = sqrt(1-(`theta`j'')^2);
						
			qui gen double `v`j'' = invnormal(`F`j'');
			qui gen double `Z`j'' = `vs'/`sr_t`j'' - `theta`j''*`v`j''/`sr_t`j'';
			qui gen double `C`j'' = normal(`Z`j'');
		};
		else if "`copula`j''" == "fgm" {;
			scalar `theta`j'' = tanh(`atheta`j'');
			qui gen double `C`j'' = `Fs'*(1+`theta`j''*(1-`Fs')*(1-`F`j'')) - `theta`j''*`F`j''*`Fs'*(1-`Fs');
		};
		else if "`copula`j''" == "plackett" {;
			tempvar r`j' root`j';
			scalar `theta`j'' = exp(`atheta`j'');
			qui gen double `r`j'' = 1+(`theta`j''-1)*(`F`j''+`Fs');
			qui gen double `root`j'' = (`r`j'')^2 - 4*`F`j''*`Fs'*`theta`j''*(`theta`j''-1);
			qui gen double `C`j'' = 0.5 - 0.5*(`r`j''-2*`Fs'*`theta`j'')*(`root`j'')^(-1/2);
		};
		else {;	//Acrhimedean copulas
			tempvar C_`j' phi1_F`j' phi1_Fs_`j' phi1_C_`j';
			if "`copula`j''" == "amh" {;
				scalar `theta`j'' = tanh(`atheta`j'');
				qui gen double `C_`j'' = `F`j''*`Fs'/(1-`theta`j''*(1-`Fs')*(1-`F`j''));	
				qui gen double `phi1_F`j'' = `theta`j''/(1-`theta`j''*(1-`F`j''))-1/`F`j'';
				qui gen double `phi1_Fs_`j'' = `theta`j''/(1-`theta`j''*(1-`Fs'))-1/`Fs';
				qui gen double `phi1_C_`j'' = `theta`j''/(1-`theta`j''*(1-`C_`j''))-1/`C_`j'';
			};
			else if "`copula`j''"=="clayton" {;
				scalar `theta`j'' = exp(`atheta`j'');
				qui gen double `C_`j'' = (`Fs'^(-`theta`j'')+`F`j''^(-`theta`j'')-1)^(-1/`theta`j'');
				qui gen double `phi1_F`j'' = -(`F`j'')^(-`theta`j''-1);
				qui gen double `phi1_Fs_`j'' = -(`Fs')^(-`theta`j''-1);
				qui gen double `phi1_C_`j'' =  -(`C_`j'')^(-`theta`j''-1);
			};
			else if "`copula`j''"=="frank" {;
				scalar `theta`j'' = `atheta`j'';
				qui gen double `C_`j'' = -ln(1+(exp(-`theta`j''*`Fs')-1)*(exp(-`theta`j''*`F`j'')-1)/(exp(-`theta`j'')-1))/`theta`j'';
				qui gen double `phi1_F`j'' = `theta`j''/(1-exp(`theta`j''*`F`j''));
				qui gen double `phi1_Fs_`j'' = `theta`j''/(1-exp(`theta`j''*`Fs'));
				qui gen double `phi1_C_`j'' = `theta`j''/(1-exp(`theta`j''*`C_`j''));
			};
			else if "`copula`j''"=="gumbel" {;
				scalar `theta`j'' = 1 + exp(`atheta`j'');
				qui gen double `C_`j'' = exp(-((-ln(`F`j''))^(`theta`j'')+(-ln(`Fs'))^(`theta`j''))^(1/`theta`j''));
				qui gen double `phi1_F`j'' = -`theta`j''*(-ln(`F`j''))^(`theta`j''-1) /`F`j'';
				qui gen double `phi1_Fs_`j'' = -`theta`j''*(-ln(`Fs'))^(`theta`j''-1) /`Fs';
				qui gen double `phi1_C_`j'' = -`theta`j''*(-ln(`C_`j''))^(`theta`j''-1) /`C_`j'';
			};
			else if "`copula`j''"=="joe" {;
				tempvar guard`j';
				scalar `theta`j'' = 1 + exp(`atheta`j'');
				qui gen double `C_`j'' = 1-((1-`Fs')^`theta`j''+(1-`F`j'')^`theta`j''-((1-`Fs')*(1-`F`j''))^`theta`j'')^(1/`theta`j'');
				qui gen double `phi1_F`j'' = -`theta`j''*(1-`F`j'')^(`theta`j''-1)/(1-(1-`F`j'')^`theta`j'');
				qui gen double `phi1_Fs_`j'' = -`theta`j''*(1-`Fs')^(`theta`j''-1)/(1-(1-`Fs')^`theta`j'');
				qui gen double `phi1_C_`j'' = -`theta`j''*(1-`C_`j'')^(`theta`j''-1)/(1-(1-`C_`j'')^`theta`j'');
				qui gen double `guard`j'' = 1 - `Fs';
			};
			qui gen double `C`j'' = `phi1_F`j''/`phi1_C_`j'';
		};
	};	//end of forvalue;
	
	* General * ;
	qui replace `C0' = 0 if `Fs' == 0;
	qui replace `C0' = 1 if `Fs' == 1;
	qui replace `C1' = 0 if `Fs' == 0;
	qui replace `C1' = 1 if `Fs' == 1;
	capture replace `C1' = 0 if `guard1' == 1;
	capture replace `C0' = 0 if `guard0' == 1;
	
	*** Log-likelihood ***;
	generate `typlist' `varlist' = cond(`S'==0, cond(`C0'!=., ln(`C0') + `lnf0'-`lnsig0',ln(`Fs')),cond(`C1'!=.,ln(1-`C1') + `lnf1'-`lnsig1',ln(1-`Fs'))) if `touse' ;
	exit;
	}; //end of ~missing("`cll'");
	
end;
*include num_inte.mata;

*! version 1.1.0 June 11, 2012: revise numerical itegration
*! version 1.0.1 January 2012
