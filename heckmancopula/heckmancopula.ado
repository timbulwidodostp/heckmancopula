*! version 2.1.0  August 7, 2012 Takuya Hasebe 

* version 2.1.0 September 24, 2012: minor revision;
* version 1.2.0 August 7, 2012: allow negative dependence and minor changes;
* wrap file to excute copula-based maximum likelihood estimation of sample selection model;

#delimit ;

program define heckmancopula, sortpreserve;
	if replay() {;
		if ("`e(cmd)'" != "heckmancopula") error 301 ;
		Replay `0' ;
	};
	else Estimate `0';
end;


program Estimate, eclass sortpreserve;
	version 11;
	syntax anything(id="equations" equalok) [aweight pweight iweight fweight] [if] [in], SELect(string) 
		[ copula(string) margin1(string) margsel(string) negative NOCONStant 
		df(string)
		*];
	
	mlopts mlopts, `options';
	
	local margins `margsel';	//rename local in an orginal way;
	local margin `margin1';
	
	/** dealing with copulas and margins **/;
	if "`copula'" == "" {; local copula "gaussian";}; 
	if "`margin'" == "" {; local margin "normal";}; 
	if "`margins'" == "" | "`margins'" == "normal" {; local margins "probit";}; 
	if "`margins'" == "logistic" {; local margins "logit" ;};

	if "`margins'" != "" & "`margins'" != "probit" & "`margins'" != "logit" {;
		dis as error "`margins' is not an available margin."; exit 198;
	};
		
	if  "`copula'" != "gaussian" & "`copula'" != "fgm" & "`copula'" != "amh" & "`copula'" != "product" &
		"`copula'" != "clayton" & "`copula'" != "frank" & "`copula'" != "gumbel" & "`copula'" != "joe" &
		"`copula'" != "plackett"  {; 
		dis as error "`copula' is not an available copula."; exit 198;
	};	

	if "`margin'" != "normal" & "`margin'" != "logistic" & "`margin'" != "t" {;
		dis as error "`margin' is not an available margin."; exit 198;
	};
	
	if "`margin'"=="t" & "`df'"!="" {;
		if `df' <= 0 {;
			dis as error "Degree of freedom must be positive value.";
			exit 198;
		};
	};
	
	* option no constant *;
	if "`noconstant'"!="" {; local _constant , noconstant;};
	
	*** globals for ml estimation ***;
	local temp_global copula margin margins gradient negative df ;
	foreach v of local temp_global {;
		if "$`v'" != "" {; local temp_`v' $`v'; };
		global `v' ``v'';
	};
	
	tokenize `anything', parse("=");
	if "`2'" != "=" {;
		tokenize `anything';
		local y `1';
		macro shift;
		local x `*';
	};
	else {;
		local y `1';
		local x `3';
	};
	
	capture unab x : `x';
	
	tokenize `select', parse("=");
	if "`2'" != "=" {;
		local x_s `select';
		tempvar y_s;
		qui gen `y_s' = (`y'!=.);
	};
	else {;
		local y_s `1';
		local x_s `3';
	};
	
	capture unab x_s : `x_s';
	
	marksample touse;
	markout `touse' `y_s' `x_s';
	tempvar touse1;
	mark `touse1' if `y_s'==1;
	markout `touse1' `y' `x';
	qui replace `touse' = 0 if `touse1'==0 & `y_s'==1; //(`y_s'==1 & `y'==.);

	qui tab `y_s' if `touse';
	if r(r) != 2 {; dis as error "`y_s' should be binary."; exit 198; };
	qui sum `y_s' if `touse';
	if r(max) != 1  & r(min) != 0 {;dis as error "`y_s' should be either 0 or 1."; exit 198; };
	
	*** checking collinearity ***;
	qui _rmdcoll `y' `x' if `touse' `_constant';
	local result "`r(varlist)'";
	local coll_x: list x - result;
	if ~missing("`coll_x'") {;
		noisily display as text "note: `coll_x' omitted from outcome equation because of collinearity";
		local x `result';
	};
	qui _rmdcoll `y_s' `x_s' if `touse';
	local result "`r(varlist)'";
	local coll_x: list x_s - result;
	if ~missing("`coll_x'") {;
		noisily display as text "note: `coll_x' omitted from selection equation because of colliearity";
		local x_s `result';
	};
	
	*** obtaining likelihood for the case where all errors are independent *** ;
	*** and use these estimates as starting values *** ;
	tempname start_s start_b start_sig  ;
	
	*dis _newline as text "fitting initial values from `margins' model:";
	capture `margins' `y_s' `x_s' if `touse', `mlopts';  
	matrix `start_s' = e(b);
	local ll_s = `e(ll)';
		
	if "`margin'" == "t" {;
		if "`df'"=="" local ml_v /lndf;
	};
	if "`copula'" != "product" {; local ml_theta /atheta;};
	
	if "`margin'" == "normal" {;
		capture reg `y' `x' if `y'!=. & `y_s' == 1 & `touse' `_constant'; 
		tempname beta; 
		matrix `start_b' = e(b);
		scalar `start_sig' = 0.5*ln(e(rss)/e(N));
	};
	else if "`margin'" != "normal" {;
		capture reg `y' `x' if `y'!=. & `y_s' == 1 & `touse'; 
		capture ml model lf2 ml_reg_`margin' (eq: `y' = `x' `_constant') /lnsig `ml_v' 
			if `y'!=. & `y_s' == 1 & `touse', maximize 
			initial(e(b) 0 2, copy) search(off) `mlopts'; 
		tempname beta; 
		matrix `beta' = e(b);
		matrix `start_b' = `beta'[1,"eq:"];
		matrix `start_sig' = _b[/lnsig];
		
		if "`margin'" == "t" {;
			if "`df'" == "" {;
				tempname start_v;
				matrix `start_v' = _b[`ml_v']; 
			};
		};
	};

	local ll_r = `e(ll)';
		
	tempname start_value ;
	matrix `start_value' = (`start_s', `start_b', `start_sig');
	
	if "`copula'" != "product" {;
		tempname start_t;
		if "`init_theta'"=="" {;
			matrix `start_t' = 0.01;
		};
		else {; matrix `start_t' = `init_theta'; };
	};
			
	*** ML ESTIMATION *** ;
	ml model lf2 heckmancopula_ml (select: `y_s' = `x_s' `_constant') (`y': `y' = `x')  
		/lnsigma `ml_v' `ml_theta'  
		[`weight'`exp']
		if `touse', missing collinear
		title(Sample Selection Model: Copula `negative' `copula', Margins `margins'-`margin')
		init(`start_value' `start_v' `start_t'  , copy) 
		maximize `vce' `mlopts' search(off) 
		;
	*ml check; // 
		
	*** ereturn's ***;
	* scalar *;
	*ereturn scalar k_aux = `e(k_eq)'-`e(k_dv)';
	ereturn scalar ll0 = `ll_s' + `ll_r';	//log likelihood of the case where all errors are independent;
	ereturn scalar AIC = -2*`e(ll)'+2*`e(k)';
	ereturn scalar BIC = -2*`e(ll)'+ln(`e(N)')*`e(k)';
	if "`margin'"=="t" & "`df'" != "" ereturn scalar df = `df'; 
	ereturn scalar negative = ("`negative'"!="");
	
	* macro *;
	ereturn local cmd heckmancopula;
	ereturn local copula `copula';
	ereturn local margin1 `margin';
	ereturn local margsel `margins';	
	ereturn local predict "heckmancopula_p";
	
	Replay, level(`level');
	
	*** replace global back ***; 
	foreach v of local temp_global {;
		if "temp_`v'" != "" {; global `v' `temp_`v''; };
		else  macro drop `v';
	};
	
end;

program Replay;
	syntax [, Level(cilevel)];
	
	local n_equation = 3;	
	local num_d = 0;	//number of copulas that are not product;
	local num_b = 0;	//number of copulas that independence occurs at boundary;
		
	local copula `e(copula)';
	local margin `e(margin1)';		

	if "`e(negative)'"=="" {; local sign = 1; };
	else {; local sign = -1;};	
		
	if "`copula'" == "gaussian" | "`copula'" == "amh" | "`copula'" == "fgm" {;
		local diparm_t diparm(atheta, tanh label("theta")) ;
		local test (tanh(_b[/atheta])=0);
			if "`copula'" == "gaussian" {;
				local diparm_tau diparm(atheta, function(`sign'*2*asin(tanh(@))/_pi)
					derivative(`sign'*(2*sqrt(1-tanh(@)^2)/_pi)) label("tau"));
			};
			else if "`copula'" == "amh" {;
				local diparm_tau diparm(atheta, function(`sign'*((3*tanh(@)-2)/(3*tanh(@)) - (2/3)*(1-1/tanh(@))^2 *ln(1-tanh(@)))) 
					derivative(`sign'*(-(2/3)*(-2+tanh(@)+2*(1-1/tanh(@))*ln(1-tanh(@)))/(tanh(@)^2)*(1-tanh(@)^2))) label("tau"));
			};
			else if "`copula'" == "fgm" {;
				local diparm_tau diparm(atheta, function(`sign'*(2*tanh(@)/9))
					derivative(`sign'*(2*(1-tanh(@)^2)/9)) label("tau"));
			};
	};
	else if "`copula'" == "clayton" | "`copula'" == "plackett" {;
		local diparm_t diparm(atheta, exp label("theta"));
			if "`copula'" == "clayton" {;
				local num_b = 1; 
				local diparm_tau diparm(atheta, function(`sign'*(exp(@)/(exp(@)+2)))
					derivative(`sign'*(2*exp(@)/((exp(@)+2)^2))) label("tau")) ;
				local test (exp(_b[/atheta])=0);
			};
			if "`copula'" == "plackett" {;
				local test (exp(_b[/atheta])=1);
			};
	};
	else if "`copula'" == "gumbel" | "`copula'" == "joe" {;
		local diparm_t diparm(atheta, function(1+exp(@)) derivative(exp(@)) label("theta"));
		local test (1+exp(_b[/atheta])=1);
		local num_b = 1;
			if "`copula'" == "gumbel" {;	
				local diparm_tau diparm(atheta, function(`sign'*(1-1/(1+exp(@)))) derivative(`sign'*((1+exp(@))^(-2)*(exp(@)))) label("tau")) ;
			};
			else if "`copula'" == "joe" {;
				local diparm_tau ;
				qui _tau 1+exp(_b[/atheta]), copula(`copula') dis;
				local tau = `sign'*r(tau);
				local tau_display dis in text _col(10) "tau" _col(18) as result `tau';
			};
	};
	else if "`copula'" == "frank" {;
		local diparm_t diparm(atheta, function(@)
			derivative(1) label("theta"));
		local diparm_tau ;
		qui _tau _b[/atheta], copula(`copula') dis;
		local tau = `sign'*r(tau);
		local tau_display dis in text _col(10) "tau" _col(18) as result `tau';
		local test (_b[/atheta]=0);
	};
	
	if "`margin'" == "t" {;
		if "`e(df)'"=="" {;
			local diparm_v diparm(lndf, exp label("df"));
			local n_equation = `n_equation' + 1;
		};
		else {;
			local fix_df dis as res "df " as text "is fixed at " as res "`e(df)'";
		};
	};
		
	if "`copula'" != "product" {;
		local n_equation = `n_equation' + 1;
		local num_d = 1;
	};
		
	*** obtain test statistic ***;
	if "`copula'"! = "product" {;
		if "`e(vce)'" == "oim" | "`e(vce)'" == "opg" {; 
			local test_type LR;
			local test_stat = 2*(`e(ll)'-`e(ll0)');
		};
		else {; 
			local test_type Wald;
			qui testnl `test'; 
			local test_stat = r(chi2);
		};
	};
	
	if `num_d' == 0 {; /**no need to test**/ };
	if `num_d' == 1 & `num_b' == 0 {; local p_v = 1-chi2(1,`test_stat'); };
	if `num_d' == 1 & `num_b' == 1 {; local p_v = 1-(0.5*1+0.5*chi2(1,`test_stat')); };
	ml display, level(`level') neq(`n_equation')
		`diparm_s' `diparm_t' `diparm_v' `diparm_tau';
	`tau_display';
	`fix_df';	
	if "`copula'" != "product"  dis in text "`test_type' test of independence :" _col(34) "Test statistic " % 8.3f `test_stat' " with p-value" %8.4f `p_v'; 
	display	in smcl "{hline 78}";
	
	
end;	
