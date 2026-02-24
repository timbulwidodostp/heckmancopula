*! version 1.2.0 August 7, 2012 T.Hasebe

* version 1.2.0 August 7, 2012: allow negative dependence & minor changes;
* version 1.0.1 December 10, 2011:  


#delimit ;

program define switchcopula, sortpreserve;
	if replay() {;
		if ("`e(cmd)'" != "switchcopula") error 301 ;
		Replay `0' ;
	};
	else Estimate `0';
end;


program Estimate, eclass sortpreserve;
	version 11;
	syntax anything(id="equations" equalok) [aweight pweight iweight fweight] [if] [in], SELect(string) 
		[ copula0(string) copula1(string) margin0(string) margin1(string) margsel(string)
		df0(string) df1(string)
		CONsel
		negative0
		negative1
		noCONSTANT0
		noCONSTANT1
		vce(passthru) *];

	mlopts mlopts, `options';
	
	local margins `margsel';
	
	/** dealing with copulas and margins **/;
	if "`copula0'" == "" {; local copula0 "gaussian";}; if "`copula1'" == "" {; local copula1 "gaussian";};
	if "`margin0'" == "" {; local margin0 "normal";}; if "`margin1'" == "" {; local margin1 "normal";};
	if "`margins'" == "" | "`margins'" == "normal" {; local margins "probit";}; 
	if "`margins'" == "logistic" {; local margins "logit" ;};

	if "`margins'" != "" & "`margins'" != "probit" & "`margins'" != "logit" {;
		dis as error "`margins' is not an available margin."; exit 198;
	};
		
	foreach eq_num in 0 1 {;
		if  "`copula`eq_num''" != "gaussian" & "`copula`eq_num''" != "fgm" & "`copula`eq_num''" != "amh" & "`copula`eq_num''" != "product" &
			"`copula`eq_num''" != "clayton" & "`copula`eq_num''" != "frank" & "`copula`eq_num''" != "gumbel" & "`copula`eq_num''" != "joe" &
			"`copula`eq_num''" != "plackett"  {; //& "`copula`eq_num''" != "AP"//
			dis as error "`copula`eq_num'' is not an available copula."; exit 198;
		};	
	
		if "`margin`eq_num''" != "normal" & "`margin`eq_num''" != "logistic" & "`margin`eq_num''" != "t" {;
			dis as error "`margin`eq_num'' is not an available margin."; exit 198;
		};
	}; 
	
	* option no constant *;
	if "`constant0'"!="" {; local _constant0 , noconstant;};
	if "`constant1'"!="" {; local _constant1 , noconstant;};
	
	*** globals for ml estimation ***;
	local temp_global copula0 copula1 margin0 margin1 margins negative0 negative1 df0 df1 touseo;
	foreach v of local temp_global  {;
		if "$`v'" != "" {; local temp_`v' $`v'; };
		global `v' ``v'';
	};
	

	gettoken reg0 reg1_ : anything, match(parns) bind;
	if "`parns'" != "(" {;
		local reg0 "`anything'";
		local reg1_ "";
	};
	tokenize "`reg0'", parse("()");
	if "`1'" != "`reg0'" {;
		dis as error "If more than one, equations must be enclosed in brackets.";
		exit 198;
	};
	tokenize "`reg0'", parse("=");
	if "`2'" != "=" {;
		tokenize "`reg0'";
		local y0 `1';
		macro shift;
		local x0 `*';
	};
	else {;
		if "`4'" == "=" {;
			dis error "If more than one, equations must be enclosed in brackets.";
			exit 198;
		};
		local y0 `1';
		local x0 `3';
	};
	capture unab x0 : `x0';
	
	if "`reg1_'" != "" {;
		gettoken reg1 : reg1_, match(parns) bind;
		if "`parns'" != "(" {;
			local reg1 `reg1_';
		};
		tokenize "`reg1'", parse("=");
		if "`2'" != "=" {;
			tokenize "`reg1'";
			local y1 `1';
			macro shift;
			local x1 `*';
		};
		else {;
			if "`4'" == "=" {;
				dis as error "If more than one, equations must be enclosed in brackets.";
				exit 198;
			};
			local y1 `1';
			local x1 `3';
		};
	};
	else {;
		local y1 `y0'; local x1 `x0';
	};
	capture unab x1: `x1';
	
	tokenize `select', parse("=");
	if "`2'" != "=" {;
		tokenize `select'; local y_s `1'; macro shift; local x_s `*';
	};
	else {;	local y_s `1'; local x_s `3'; };
	capture unab x_s : `x_s';
	
	if "`consel'"!="" {;
		marksample touse;
		markout `touse' `y_s' `x_s';
		marksample touse0; markout `touse0' `y0' `x0' `y_s' `x_s';
		marksample touse1; markout `touse1' `y1' `x1' `y_s' `x_s'; 
		tempvar touseo touses;
		qui gen byte `touseo' = 1 if (`touse1'==1 & `y_s'==1) | (`touse0'==1 & `y_s'==0) ;
	};	
	else {;	
		marksample touse;
		markout `touse' `y0' `x0' `y1' `x1' `y_s' `x_s' ;
		tempvar touseo;
		qui gen byte `touseo' = `touse';
		tempvar touse1 touse0 touses;
		qui gen `touse1' = `touse';
		qui gen `touse0' = `touse';
		qui gen `touses' = `touse';
	};
	global touseo `touseo';

	qui tab `y_s' if `touse';
	if r(r) != 2 {; dis as error "`y_s' should be binary."; exit 198; };
	qui sum `y_s' if `touse';
	if r(max) != 1  & r(min) != 0 {;dis as error "`y_s' should be either 0 or 1."; exit 198; };
	
	*** checking collinearity ***;
	qui _rmdcoll `y_s' `x_s' if `touse';
	local result "`r(varlist)'";
	local coll_x: list x_s - result;
	if ~missing("`coll_x'") {;
		noisily display as text "note: `coll_x' omitted from selection equation because of colliearity";
		local x_s `result';
	};
		
	forvalue j = 0/1 {;
		qui _rmdcoll `y`j'' `x`j'' if `touse' `_constant`j'';
		local result "`r(varlist)'";
		local coll_x: list x`j' - result;
		if ~missing("`coll_x'") {;
			noisily display as text "note: `coll_x' omitted from outcome `j' because of colliearity";
			local x`j' `result';
		};
	};

	
	*** obtaining likelihood for the case where all errors are independent *** ;
	*** and use these estimates as starting values *** ;
	tempname start_s start_b0 start_b1 start_sig0 start_sig1 ;
	
	if "`margins'" == "probit" {;
		qui probit `y_s' `x_s' if `touse';	// , `mlopts'; 
	};
	else if "`margins'" == "logit" {;
		qui logit `y_s' `x_s' if `touse';	//, `mlopts';  
	};
	matrix `start_s' = e(b);
	local ll_s = `e(ll)';
	
	local eq = 5;
	
	foreach eq_num in 0 1 {;
		if "`margin`eq_num''" == "t" {;
			if "`df`eq_num''"=="" {;
				local ml_v`eq_num' /lndf`eq_num';
				tempname start_v`eq_num' ;
				matrix `start_v`eq_num'' = 10;
			};
			else {; global df = `df`eq_num''; };
		};
		
		if "`margin`eq_num''" == "normal" {;
			capture reg `y`eq_num'' `x`eq_num'' if `y`eq_num''!=. & `y_s' == `eq_num' & `touse' `_constant`eq_num''; 
			tempname beta; 
			matrix `start_b`eq_num'' = e(b);
			scalar `start_sig`eq_num'' = 0.5*ln(e(rss)/e(N));
		};
		else if "`margin`eq_num''" != "normal" {;
			capture reg `y`eq_num'' `x`eq_num'' if `y`eq_num''!=. & `y_s' == `eq_num' & `touse';
			tempname beta lnsig_i;
			matrix `beta' = e(b);
			matrix `lnsig_i' = 0.5*ln(e(rss)/e(N));
			
			qui ml model lf2 ml_reg_`margin`eq_num'' (eq`eq_num': `y`eq_num'' = `x`eq_num'' `_constant`eq_num'') /lnsig `ml_v`eq_num'' 
				if `y`eq_num''!=. & `y_s' == `eq_num' & `touse', 
				maximize init(`beta' `lnsig_i' `start_v`eq_num'', copy) search(off) `mlopts'
				; 
			tempname beta; 
			matrix `beta' = e(b);
			matrix `start_b`eq_num'' = `beta'[1,"eq`eq_num':"];
			capture scalar `start_sig`eq_num'' = _b[/lnsig];
		};
		local ll_`eq_num' = `e(ll)';
		
		if "`margin`eq_num''" == "t" {;
			if "`df`eq_num''" == "" {;
				matrix `start_v`eq_num'' = _b[`ml_v`eq_num'']; 
				if exp(_b[`ml_v`eq_num'']) > 200 matrix `start_v`eq_num'' == ln(120);
				local eq = `eq'+1;
			};	
		};
	};
	
	tempname start_value ;
	matrix `start_value' = (`start_s', `start_b0', `start_b1', `start_sig0', `start_sig1');
	
	foreach eq_num in 0 1 {;
		if "`copula`eq_num''" != "product" {;
			tempname start_t`eq_num' ;
			local ml_theta`eq_num' /atheta`eq_num';
			matrix `start_t`eq_num'' = 0.01 ;	//is this a good initial value?
			local eq = `eq'+1;
		};
		else local start_t`eq_num' "";
	};

	*** ML ESTIMATION *** ;
	*set trace on;
	ml model lf2 switchcopula_ml (select: `y_s' = `x_s') (regime0: `y0' = `x0' `_constant0') (regime1: `y1' = `x1' `_constant1') 
		/lnsigma0 /lnsigma1 `ml_v0' `ml_v1' `ml_theta0' `ml_theta1'  
		[`weight'`exp']
		if `touse'
		,
		title(Swithching Regression: Copulas `copula0'-`copula1', Margins `margins'-`margin0'-`margin1')
		missing
		init(`start_value' `start_v0' `start_v1' `start_t0' `start_t1' , copy) 
		maximize search(off)  `vce' `mlopts' 
		;
	*ml check;
	
	*** ereturn's ***;
	* scalar *;
	ereturn scalar ll0 = `ll_s' + `ll_0' + `ll_1';	//log likelihood of the case where all errors are independent;
	ereturn scalar AIC = -2*`e(ll)'+2*`e(k)';
	ereturn scalar BIC = -2*`e(ll)'+ln(`e(N)')*`e(k)';
	if "`margin0'"=="t" & "`df0'" != "" ereturn scalar df0 = `df0'; 
	if "`margin1'"=="t" & "`df1'" != "" ereturn scalar df1 = `df1'; 
	ereturn scalar negative0 = ("`negative0'"!="");
	ereturn scalar negative1 = ("`negative1'"!="");
	ereturn scalar consel = ("`consel'"!="");
	
	* macro *;
	ereturn local cmd switchcopula;
	ereturn local copula0 `copula0';
	ereturn local copula1 `copula1';
	ereturn local margin0 `margin0';
	ereturn local margin1 `margin1';
	ereturn local margsel `margins';	
	
	ereturn local predict "switchcopula_p";
	
	
	Replay, level(`level');

	*** replace global back ***; 
	foreach v of local temp_global  {;
		if "temp_`v'" != "" {; global `v' `temp_`v''; };
		else  macro drop `v';
	};
	
end;


************************************************************************;
program Replay;
	version 8.1;
	syntax [, Level(cilevel)];
	
	local n_equation = 5;	
	local num_d = 0;	//number of copulas that are not product;
	local num_b = 0;	//number of copulas that independence occurs at boundary;
	
	foreach j in 0 1 {;
		
		local copula`j' `e(copula`j')';
		local margin`j' `e(margin`j')';		
		
		if "`e(negative`j')'"=="" {; local sign`j' = 1; };
		else {; local sign`j' = -1;};	
				
		if "`copula`j''" == "gaussian" | "`copula`j''" == "amh" | "`copula`j''" == "fgm" {;
			local diparm_t`j' diparm(atheta`j', tanh label("theta`j'")) ;
			local test`j' (tanh(_b[/atheta`j'])=0);
				if "`copula`j''" == "gaussian" {;
					local diparm_tau`j' diparm(atheta`j', function(`sign`j''*2*asin(tanh(@))/_pi)
						derivative(`sign`j''*(2*sqrt(1-tanh(@)^2)/_pi)) label("tau`j'"));
				};
				else if "`copula`j''" == "amh" {;
					local diparm_tau`j' diparm(atheta`j', function(`sign`j''*((3*tanh(@)-2)/(3*tanh(@)) - (2/3)*(1-1/tanh(@))^2 *ln(1-tanh(@)))) 
						derivative(`sign`j''*(-(2/3)*(-2+tanh(@)+2*(1-1/tanh(@))*ln(1-tanh(@)))/(tanh(@)^2)*(1-tanh(@)^2))) label("tau`j'"));
				};
				else if "`copula`j''" == "fgm" {;
					local diparm_tau`j' diparm(atheta`j', function(`sign`j''*(2*tanh(@)/9))
						derivative(`sign`j''*(2*(1-tanh(@)^2)/9)) label("tau`j'"));
				};
		};
		else if "`copula`j''" == "clayton" | "`copula`j''" == "plackett" {;
			local diparm_t`j' diparm(atheta`j', exp label("theta`j'"));
				if "`copula`j''" == "clayton" {;
					local num_b = 1; 
					local diparm_tau`j' diparm(atheta`j', function(`sign`j''*(exp(@)/(exp(@)+2)))
						derivative(`sign`j''*(2*exp(@)/((exp(@)+2)^2))) label("tau`j'")) ;
					local test`j' (exp(_b[/atheta`j'])=0);
				};
				if "`copula`j''" == "plackett" {;
					local test`j' (exp(_b[/atheta`j'])=1);
				};
		};
		else if "`copula`j''" == "gumbel" | "`copula`j''" == "joe" {;
			local diparm_t`j' diparm(atheta`j', function(1+exp(@)) derivative(exp(@)) label("theta`j'"));
			local test`j' (1+exp(_b[/atheta`j'])=1);
			local num_b = 1;
				if "`copula`j''" == "gumbel" {;	
					local diparm_tau`j' diparm(atheta`j', function(`sign`j''*(1-1/(1+exp(@)))) derivative(`sign`j''*((1+exp(@))^(-2)*(exp(@)))) label("tau`j'")) ;
				};
				else if "`copula`j''" == "joe" {;
					local diparm_tau`j' ;
					qui _tau 1+exp(_b[/atheta`j']), copula(`copula`j'') dis;
					local tau`j' = `sign`j''*r(tau);
					local tau_display`j' dis in text _col(10) "tau`j'" _col(18) as result `tau`j'';
				};
		};
		else if "`copula`j''" == "frank" {;
			local diparm_t`j' diparm(atheta`j', function(@)
				derivative(1) label("theta`j'"));
			local diparm_tau`j' ;
			qui _tau _b[/atheta`j'], copula(`copula`j'') dis;
			local tau`j' = `sign`j''*r(tau);
			local tau_display`j' dis in text _col(10) "tau`j'" _col(18) as result `tau`j'';
			local test`j' (_b[/atheta`j']=0);
		};		

		local diparm_s`j' diparm(lnsigma`j', exp label("sigma`j'")) ;
		
		if "`margin`j''" == "t" {;
			if "`e(df`j')'"=="" {;
				local diparm_v`j' diparm(lndf`j', exp label("df`j'"));
				local n_equation = `n_equation' + 1;
			};
			else {;
				local fix_df`j' dis as res "df`j' " as text "is fixed at " as res "`e(df`j')'";
			};
		};
		if "`copula`j''" != "product" {;
			local n_equation = `n_equation' + 1;
			local num_d = `num_d' + 1;
		};
		
	};
	
	*** obtain test statistic ***;
	if "`e(vce)'" == "oim" | "`e(vce)'" == "opg" {; 
		local test_type LR;
		local test_stat = 2*(`e(ll)'-`e(ll0)');
	};
	else {; 
		local test_type Wald;
		qui testnl `test0' `test1'; 
		local test_stat = r(chi2);
	};
	
	if `num_d' == 0 {; /**no need to test**/ };
	if `num_d' == 1 & `num_b' == 0 {; local p_v = 1-chi2(1,`test_stat'); };
	if `num_d' == 1 & `num_b' == 1 {; local p_v = 1-(0.5*1+0.5*chi2(1,`test_stat')); };
	if `num_d' == 2 & `num_b' == 0 {; local p_v = 1-chi2(2,`test_stat'); };
	if `num_d' == 2 & `num_b' == 1 {; local p_v = 1-(0.5*chi2(1,`test_stat')+0.5*chi2(2,`test_stat')) ;};
	if `num_d' == 2 & `num_b' == 2 {; local p_v = 1-(0.25*1+ 0.5*chi2(1,`test_stat')+0.25*chi2(2,`test_stat')); };	

	ml display, level(`level') neq(`n_equation')
		`diparm_s0' `diparm_s1' `diparm_t0' `diparm_t1' `diparm_v0' `diparm_v1' `diparm_tau0' `diparm_tau1';
		`tau_display0'; `tau_display1';
	`fix_df0';
	`fix_df1';
	dis in text "`test_type' test of independence :" _col(34) "Test statistic " % 8.3f `test_stat' " with p-value" %8.4f `p_v'; 
	display	in smcl "{hline 78}";
end;	
