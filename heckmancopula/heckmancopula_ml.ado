*! version 1.2.0 August 7, 2012 Takuya Hasebe

* Version 1.2.0 August 7, 2012: allow negative dependence of Clayton, Gumbel, Joe (and others as well);

#delimit ;
*** Copula based Heckman Model with Flexible Margins and Copula ***;
program define heckmancopula_ml;

	local margin "$margin";
	local copula "$copula";
	local margins "$margins";
	local negative1 "$negative";
	
	local eq = 3;
	if "`margin'" == "t" & "$df"=="" {;
		local eq = `eq' + 1;
	};
	
	if "`copula'" != "product" {;
		local eq = `eq' + 1;
	};
	
	forvalue i = 1/`eq' {;
		local gradient `gradient' g`i';
	};
	
	args todo b lnf `gradient' Hessian;
	
		
	tempvar S z y xb1 e1 ;
	tempname lnsig1 sig1 atheta1 theta1 lndf1 df1 ;
	
	mleval `z' = `b', eq(1);	//selection equation;
	mleval `xb1' = `b', eq(2);	//outcome equation;
	mleval `lnsig1' = `b', eq(3) scalar; scalar `sig1' = exp(`lnsig1');	//scale parameter for outcome equation;
	
	local eq = 4;
	
	
	if "`margin'" == "t" {;
		if "$df" == "" {;
			mleval `lndf1' = `b', eq(`eq') scalar;	//ancillary degree of freedom;
			scalar `df1' = exp(`lndf1');
			local eq  = `eq' + 1;
		};
		else scalar `df1' = $df;
	};
	
	if "`copula'" != "product" {;
		mleval `atheta1' = `b', eq(`eq') scalar;	//ancillary depedence parameter;
		local eq = `eq' + 1;
	};
	

	*** dependent variables ***;
	qui gen double `S' = $ML_y1;	//selection indicator;
	qui gen double `y' = $ML_y2;		//outcome;
	
	*** standarized continuous dependent variable ***;
	if "`negative1'"=="" {; local sign1 = 1; };
	else {; local sign1 = -1; };
	
	qui gen double `e1' = `sign1'*(`y'-`xb1')/`sig1';
	
	*** Specify Margins ***;
	tempvar Fs fs F1 f1 lnf1;
	
	* selection *;
	if "`margins'"=="probit" {;
		qui gen double `Fs' = normal(-`z');
		qui gen double `fs' = normalden(`z');
	};
	else if "`margins'"=="logit" {;
		qui gen double `Fs' = 1/(1+exp(`z'));
		qui gen double `fs' = `Fs'*(1-`Fs');
	};
	
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
	qui replace `lnf' = cond(`S'==0, ln(`Fs'), ln(1-`C1') + `lnf1'-`lnsig1');
	
	if (`todo'==0) exit ;
	
	*********************************************************************************************************************;
	*** Gradient *** ;
	
	//margin-specific ;
	tempvar dlnf1db dlnf1ds dF1db dF1ds dFsdb dF1dv dlnf1dv;

	* outcome *;
	if "`margin'" == "normal" {;
		qui gen double `dF1db' = -`f1'/`sig1';
		qui gen double `dF1ds' = -`f1'*`e1';
		qui gen double `dlnf1db' = `e1'/`sig1';
		qui gen double `dlnf1ds' = (`e1')^2;
	};
	else if "`margin'" == "logistic" {;
		qui gen double `dF1db' = -`f1'/`sig1';
		qui gen double `dF1ds' = -`f1'*`e1';
		qui gen double `dlnf1db' = 1/`sig1' - 2*(1-`F1')/`sig1';
		qui gen double `dlnf1ds' = `e1' - 2*(1-`F1')*`e1';
	};
	else if "`margin'" == "t" {;
		qui gen double `dlnf1db' = ((`df1'+1)/(`df1'+(`e1')^2))*`e1'/`sig1';
		qui gen double `dlnf1ds' = ((`df1'+1)/(`df1'+(`e1')^2))*(`e1')^2;
		qui gen double `dlnf1dv' = 0.5*(digamma((`df1'+1)/2)-1/`df1'-digamma(`df1'/2)-ln(1+(`e1')^2/`df1')+(`df1'+1)*(`e1'/`df1')^2 /(1+(`e1')^2/`df1'));
		qui replace `dlnf1dv' = `dlnf1dv'*`df1';
		qui gen double `dF1db' = -`f1'/`sig1';
		qui gen double `dF1ds' = -`f1'*`e1';
		* numerical derivative *;
		tempname d ddf1;
		tempvar P_F1 M_F1;
		scalar `d' = (abs(`df1')+10^(-8))/(10^8);
		scalar `ddf1' = `df1' + `d'/2;
		qui gen double `P_F1' = 1-ttail(`ddf1',`e1');
		scalar `ddf1' = `df1' - `d'/2;
		qui gen double `M_F1' = 1-ttail(`ddf1',`e1');
		qui gen double `dF1dv' = (`P_F1'-`M_F1')/(`d');
		qui replace `dF1dv' = `dF1dv'*`df1';
	};
	
	* selection *;
	if "`margins'"=="probit" {;
		qui gen double `dFsdb' = -`fs';
	};
	else if "`margins'"=="logit" {;
		qui gen double `dFsdb' = -`fs';
	};
	
	//copula-specific ;
	tempvar C11 C1s dC1dt1;
	tempname dt1dat1;
	
	*** Between Selection and Eq.1 ***;	
	if "`copula'" == "product" {;
		scalar `dt1dat1'=0;
		qui gen double `C11' = 0;
		qui gen double `C1s' = 1;
		qui gen double `dC1dt1' = 0;
	};
	else if "`copula'"=="gaussian" {;
		tempvar dZ1dt1;
		
		scalar `dt1dat1' = (`sr_t1')^2;
		qui gen double `C11' = -`theta1'*normalden(`Z1')/(`sr_t1'*normalden(`v1'));
		qui gen double `C1s' = normalden(`Z1')/(`sr_t1'*normalden(`vs'));
		qui gen double `dZ1dt1' = `theta1'*`vs'/(`sr_t1'^3) - `v1'*(1-(`theta1')^2)^(-3/2);
		qui gen double `dC1dt1' = normalden(`Z1')*`dZ1dt1';
	};
	else if "`copula'" == "fgm" {;
		scalar `dt1dat1' = 1-(`theta1')^2;
		qui gen double `C11' = -2*`theta1'*(1-`Fs')*`Fs';
		qui gen double `C1s' = (4*`theta1'*`F1'-2*`theta1')*`Fs' - 2*`theta1'*`F1'+`theta1'+1;
		qui gen double `dC1dt1' = (1-`F1')*(1-`Fs')*`Fs' - `F1'*(1-`Fs')*`Fs';
	};	
	else if "`copula'" == "plackett" {;
		scalar `dt1dat1' = `theta1';
		qui gen double `C11' = -(0.5*(`theta1'-1))/sqrt(`root')+(0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1'))/((`root')^(3/2));
		qui gen double `C1s' = -(0.5*(-`theta1'-1))/sqrt(`root')+(0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1'))/((`root')^(3/2));
		qui gen double `dC1dt1' = -(0.5*(`F1'-`Fs'))/sqrt(`root')+(0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(3/2));
	};
	else {;
		tempvar phi2_F1 phi2_Fs phi2_C dphi1dt1_F1 dphi1dt1_Fs dphi1dt1_C dCdt1;
		
		if "`copula'" == "amh" {;
			scalar `dt1dat1' = (1-(`theta1')^2);
			qui gen double `phi2_Fs' = -(`theta1'/(1-`theta1'*(1-`Fs')))^2+1/(`Fs'^2);
			qui gen double `phi2_F1' = -(`theta1'/(1-`theta1'*(1-`F1')))^2+1/(`F1'^2);
			qui gen double `phi2_C' = -(`theta1'/(1-`theta1'*(1-`C')))^2+1/(`C'^2);
		
			qui gen double `dphi1dt1_F1' = (1-`theta1'*(1-`F1'))^(-2);
			qui gen double `dphi1dt1_Fs' = (1-`theta1'*(1-`Fs'))^(-2);
			qui gen double `dphi1dt1_C' = (1-`theta1'*(1-`C'))^(-2);
			qui gen double `dCdt1' = `C'*(1-`Fs')*(1-`F1')/(1-`theta1'*(1-`Fs')*(1-`F1'));
		};
		else if "`copula'"=="clayton" {;
			scalar `dt1dat1' = `theta1';
			qui gen double `phi2_Fs' = (`theta1'+1)*`Fs'^(-`theta1'-2);
			qui gen double `phi2_F1' = (`theta1'+1)*`F1'^(-`theta1'-2);
			qui gen double `phi2_C' = (`theta1'+1)*`C'^(-`theta1'-2);
			
			qui gen double `dphi1dt1_F1' = ln(`F1')*(`F1')^(-`theta1'-1);
			qui gen double `dphi1dt1_Fs' = ln(`Fs')*(`Fs')^(-`theta1'-1);
			qui gen double `dphi1dt1_C' =  ln(`C')*(`C')^(-`theta1'-1);
			qui gen double `dCdt1' = `C'*(ln(`F1'^(-`theta1')+`Fs'^(-`theta1')-1)*(`theta1')^(-2)+(ln(`F1')*`F1'^(-`theta1') + ln(`Fs')*`Fs'^(-`theta1'))/(`theta1'*(`F1'^(-`theta1')+`Fs'^(-`theta1')-1)));
		};
		else if "`copula'" == "frank" {;
			scalar `dt1dat1' = 1;
			qui gen double `phi2_F1' = (`theta1'/(1-exp(`theta1'*`F1')))^2 *exp(`theta1'*`F1');
			qui gen double `phi2_Fs' = (`theta1'/(1-exp(`theta1'*`Fs')))^2 *exp(`theta1'*`Fs');
			qui gen double `phi2_C' = (`theta1'/(1-exp(`theta1'*`C')))^2 *exp(`theta1'*`C');
			
			qui gen double `dphi1dt1_F1' = (1+(`theta1'*`F1'-1)*exp(`theta1'*`F1'))/((1-exp(`theta1'*`F1'))^2) ;
			qui gen double `dphi1dt1_Fs' = (1+(`theta1'*`Fs'-1)*exp(`theta1'*`Fs'))/((1-exp(`theta1'*`Fs'))^2) ;
			qui gen double `dphi1dt1_C' = (1+(`theta1'*`C'-1)*exp(`theta1'*`C'))/((1-exp(`theta1'*`C'))^2) ;
			qui gen double `dCdt1' = -`C'/`theta1' + (`Fs'*exp(-`theta1'*`Fs')*(exp(-`theta1'*`F1')-1)*(exp(-`theta1')-1)+`F1'*exp(-`theta1'*`F1')*
				(exp(-`theta1'*`Fs')-1)*(exp(-`theta1')-1)-exp(-`theta1')*(exp(-`theta1'*`Fs')-1)*(exp(-`theta1'*`F1')-1))
				/((1+(exp(-`theta1'*`Fs')-1)*(exp(-`theta1'*`F1')-1)/(exp(-`theta1')-1))*((exp(-`theta1')-1)^2)*(`theta1'));
		};
		else if "`copula'"=="gumbel" {;
			tempvar Fs1;
			scalar `dt1dat1' = exp(`atheta1');
			qui gen double `phi2_F1' = `theta1'*(-ln(`F1'))^(`theta1'-1)*((`theta1'-1)/(-ln(`F1'))+1)/(`F1'^2);
			qui gen double `phi2_Fs' = `theta1'*(-ln(`Fs'))^(`theta1'-1) *((`theta1'-1)/(-ln(`Fs'))+ 1)/(`Fs'^2);
			qui gen double `phi2_C' = `theta1'*(-ln(`C'))^(`theta1'-1) *((`theta1'-1)/(-ln(`C'))+ 1)/(`C'^2);
			
			qui gen double `dphi1dt1_F1' = -(-ln(`F1'))^(`theta1'-1)*(1+`theta1'*ln(-ln(`F1')))/`F1';
			qui gen double `dphi1dt1_Fs' = -(-ln(`Fs'))^(`theta1'-1)*(1+`theta1'*ln(-ln(`Fs')))/`Fs';
			qui gen double `dphi1dt1_C' = -(-ln(`C'))^(`theta1'-1)*(1+`theta1'*ln(-ln(`C')))/`C';
			qui gen double `Fs1' = ((-ln(`Fs'))^`theta1' + (-ln(`F1'))^`theta1');
			qui gen double `dCdt1' = `C'*((ln(`Fs1')/(`theta1'^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`theta1'+
				ln(-ln(`F1'))*(-ln(`F1'))^`theta1')/(`Fs1'*`theta1'))*`Fs1'^(1/`theta1'));
		};		
		else if "`copula'" == "joe" {;
			scalar `dt1dat1' = exp(`atheta1');
			qui gen double `phi2_F1' = (`phi1_F1')^2 - (`theta1'-1)*`phi1_F1'/(1-`F1');
			qui gen double `phi2_Fs' = (`phi1_Fs')^2 - (`theta1'-1)*`phi1_Fs'/(1-`Fs');
			qui gen double `phi2_C' = (`phi1_C')^2 - (`theta1'-1)*`phi1_C'/(1-`C');
			
			qui gen double `dphi1dt1_F1' = -(((1-`F1')^(`theta1'-1)*(1+`theta1'*ln(1-`F1')))/(1-(1-`F1')^`theta1') +
				(`theta1'*(1-`F1')^(2*`theta1'-1)*ln(1-`F1'))/((1-(1-`F1')^`theta1')^2));
			qui gen double `dphi1dt1_Fs' = -(((1-`Fs')^(`theta1'-1)*(1+`theta1'*ln(1-`Fs')))/(1-(1-`Fs')^`theta1') +
				(`theta1'*(1-`Fs')^(2*`theta1'-1)*ln(1-`Fs'))/((1-(1-`Fs')^`theta1')^2));
			qui gen double `dphi1dt1_C' = -(((1-`C')^(`theta1'-1)*(1+`theta1'*ln(1-`C')))/(1-(1-`C')^`theta1') +
				(`theta1'*(1-`C')^(2*`theta1'-1)*ln(1-`C'))/((1-(1-`C')^`theta1')^2));	
			qui gen double `dCdt1' = (`C'-1)*((ln(1-`Fs')*(1-`Fs')^`theta1' + ln(1-`F1')*(1-`F1')^`theta1' - ln((1-`Fs')*(1-`F1'))*((1-`Fs')*(1-`F1'))^`theta1')
				/((1-`Fs')^`theta1'+(1-`F1')^`theta1'-((1-`Fs')*(1-`F1'))^`theta1') /`theta1'-ln((1-`Fs')^`theta1'+(1-`F1')^`theta1'-((1-`Fs')*(1-`F1'))^`theta1')/(`theta1'^2));	
		};
		
		*** general to Archimedean copulas ***;
		qui gen double `C11'=`phi2_F1'/`phi1_C' - (`phi2_C'*(`phi1_F1')^2)/((`phi1_C')^3) ;
		qui gen double `C1s' = -`phi2_C'*`phi1_F1'*`phi1_Fs'/((`phi1_C')^3) ;
		qui gen double `dC1dt1' = `dphi1dt1_F1'/`phi1_C' - `phi1_F1'*(`dphi1dt1_C'+`phi2_C'*`dCdt1')/((`phi1_C')^2);
	};
	qui replace `C11' = 0 if `Fs'==0 | `Fs'==1;
	qui replace `C1s' = 0 if `Fs'==0 | `Fs'==1;
	qui replace `dC1dt1' = 0 if `Fs'== 0 | `Fs'==1;
	
	capture replace `C11' = 0 if `guard' == 1; capture replace `C1s' = 0 if `guard' == 1; capture replace `dC1dt1' = 0 if `guard' == 1;
	
	*** General *** ;
	qui replace `g1' = cond(`S'==0, `dFsdb'/`Fs', -`C1s'*`dFsdb'/(1-`C1') );
	qui replace `g2' = `sign1'*cond(`S'==0, 0, -`C11'*`dF1db'/(1-`C1') + `dlnf1db');
	qui replace `g3' = cond(`S'==0, 0, -`C11'*`dF1ds'/(1-`C1') + `dlnf1ds' - 1);
	
	local eq = 4;
	
	if "`margin'" == "t" {;
		if "$df" == "" {;
			qui replace `g`eq'' = cond(`S'==0, 0, -`C11'*`dF1dv'/(1-`C1') + `dlnf1dv');
			local eq = `eq' + 1;
		};
	};
	
	if "`copula'" != "product" {;
		qui replace `g`eq'' = cond(`S'==0, 0, -`dC1dt1'/(1-`C1'));
		qui replace `g`eq'' = `g`eq'' * `dt1dat1';		
	};
	
	
	
	if `todo'==1 exit;
	
	**************************************************************************************;
	*** HESSIAN ***;
	* selection *;
	tempvar d2Fsdb2 ;
	
	if "`margins'"=="probit" {;
		qui gen double `d2Fsdb2' = `z'*`fs';
	};
	else if "`margins'" == "logit" {;
		qui gen double `d2Fsdb2' = `dFsdb'*(2*`Fs'-1);
	};
	
	* outcome ;
	tempvar d2F1db2 d2F1dbds d2F1ds2 d2lnf1db2 d2lnf1dbds d2lnf1ds2 d2F1dbdv d2F1dsdv d2F1dv2 d2lnf1dbdv d2lnf1dsdv d2lnf1dv2;
	if "`margin'" == "normal" {;
		qui gen double `d2F1db2' = -`e1'*normalden(`e1')/(`sig1'^2);
		qui gen double `d2F1dbds' = (1-(`e1')^2)*normalden(`e1')/`sig1';
		qui gen double `d2F1ds2' = (1-(`e1')^2)*normalden(`e1')*`e1';
		qui gen double `d2lnf1db2' = -`sig1'^(-2);
		qui gen double `d2lnf1dbds' = -2*`e1'/`sig1';
		qui gen double `d2lnf1ds2' = -2*(`e1')^2;
		
	};
	if "`margin'" == "logistic" {;
		qui gen double `d2F1db2' = `dF1db'*(2*`F1'-1)/`sig1';
		qui gen double `d2F1dbds' = `dF1ds'*(2*`F1'-1)/`sig1' + `f1'/`sig1';
		qui gen double `d2F1ds2' = `dF1ds'*(2*`F1'-1)*`e1' + `f1'*`e1';
		
		qui gen double `d2lnf1db2' = 2*`dF1db'/`sig1';
		qui gen double `d2lnf1dbds' = -1/`sig1'+2*`dF1ds'/`sig1'+2*(1-`F1')/`sig1';
		qui gen double `d2lnf1ds2' = -`e1'+2*`dF1ds'*`e1'+2*(1-`F1')*`e1';
	};
	if "`margin'" == "t" {;
		qui gen double `d2lnf1db2' = (`df1'+1)*((`e1')^2-`df1')/((`sig1'*(`df1'+(`e1')^2))^2);
		qui gen double `d2lnf1dbds' = -2*`df1'*`e1'*(`df1'+1)/(`sig1'*(`df1'+(`e1')^2)^2);
		qui gen double `d2lnf1dbdv' = ((`e1')^2-1)*`e1'/(`sig1'*(`df1'+(`e1')^2)^2) *`df1';
		qui gen double `d2lnf1ds2' = -2*`df1'*(`df1'+1)*(`e1')^2/((`df1'+(`e1')^2)^2);
		qui gen double `d2lnf1dsdv' =  ((`e1')^2-1)*(`e1')^2/((`df1'+(`e1')^2)^2) *`df1';
		qui gen double `d2lnf1dv2' = 0.5*digamma((`df1'+1)/2)+0.25*`df1'*trigamma((`df1'+1)/2)
			-0.5*digamma(`df1'/2)-0.25*`df1'*trigamma(`df1'/2)-0.5*ln(1+(`e1')^2/`df1')
			+0.5*(`e1')^2 *(2-(`df1'+1)/(`df1'+(`e1')^2))/(`df1'+(`e1')^2);
		qui replace `d2lnf1dv2' = `d2lnf1dv2' * `df1';
		
		qui gen double `d2F1db2' = -`dlnf1db'*`f1'/`sig1';
		qui gen double `d2F1dbds' = -`dlnf1ds'*`f1'/`sig1' + `f1'/`sig1';
		qui gen double `d2F1dbdv' = -`dlnf1dv'*`f1'/`sig1';
		qui gen double `d2F1ds2' =  -`dlnf1ds'*`f1'*`e1' + `f1'*`e1';
		qui gen double `d2F1dsdv' = -`dlnf1dv'*`f1'*`e1';
		
		*** numerica derivative ***;
		* numerical derivative *;
		scalar `d' = (abs(`df1')+10^(-5))/(10^5);
		scalar `ddf1' = `df1' + `d';
		qui replace `P_F1' = 1-ttail(`ddf1',`e1');
		scalar `ddf1' = `df1' - `d';
		qui replace `M_F1' = 1-ttail(`ddf1',`e1');
		qui gen double `d2F1dv2' = (`P_F1'+`M_F1'-2*`F1')/(`d'^2);
		qui replace `d2F1dv2' = `d2F1dv2'*(`df1')^2 + `dF1dv';
	};

	*** Copula-specific ***;
	tempvar C111 C11s C1ss d2C1dt12 dC11dt1 dC1sdt1;
	tempname d2t1dat12;
	if "`copula'" == "product" {;
		scalar `d2t1dat12' = 0;
		qui gen double `C111' = 0; qui gen double `C11s' = 0;
		qui gen double `C1ss' = 0; qui gen double `d2C1dt12' = 0;
		qui gen double `dC11dt1' =0; qui gen double `dC1sdt1' = 0;
	};
	else if "`copula'" == "gaussian" {;
		tempvar dZ1dt1 d2Z1dt12 ;
		
		scalar `d2t1dat12' = -2*`theta1'*(`dt1dat1');
		qui gen double `C111' = -(`theta1')^2*`Z1'*normalden(`Z1')/((`sr_t1'*normalden(`v1'))^2)
			-`theta1'*normalden(`Z1')*`v1'/(`sr_t1'*normalden(`v1')^2);
		qui gen double `C11s' = `theta1'*`Z1'*normalden(`Z1')/(`sr_t1'^2 * normalden(`v1')*normalden(`vs'));	
		qui gen double `C1ss' = -`Z1'*normalden(`Z1')/((`sr_t1'*normalden(`vs'))^2) + normalden(`Z1')*`vs'/(`sr_t1'*normalden(`vs')^2 );
		qui gen double `dZ1dt1' = (`theta1'*`vs'-`v1')/(`sr_t1'^3);
		qui gen double `d2Z1dt12' = 3*`theta1'*(`theta1'*`vs'-`v1')/(`sr_t1'^5) + `vs'/(`sr_t1'^3);

		qui gen double `dC11dt1' = (normalden(`Z1')/normalden(`v1'))*(-(`sr_t1')^(-3) + `theta1'*`Z1'*`dZ1dt1'/`sr_t1');
		qui gen double `dC1sdt1' = (normalden(`Z1')/normalden(`vs'))*(`theta1'/(`sr_t1'^3) -`Z1'* `dZ1dt1'/`sr_t1');
		qui gen double `d2C1dt12' = -`Z1'*normalden(`Z1')*(`dZ1dt1')^2 + normalden(`Z1')*`d2Z1dt12';
	};
	else if "`copula'" == "fgm" {;
		scalar `d2t1dat12' = -2*`theta1'*(`dt1dat1');
		qui gen double `C111' = 0; 
		qui gen double `C11s' = 4*`theta1'*`Fs'-2*`theta1';
		qui gen double `C1ss' = 4*`theta1'*`F1'-2*`theta1';
		qui gen double `dC11dt1' = -2*`Fs'*(1-`Fs');
		qui gen double `dC1sdt1' = (4*`F1'-2)*`Fs' - 2*`F1' + 1;
		qui gen double `d2C1dt12' = 0;
	};
	else if "`copula'" == "plackett" {;
		scalar `d2t1dat12' = `theta1';
		qui gen double `C111' = (0.5*(`theta1'-1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1'))/((`root')^(3/2))
			+(0.5*(`theta1'-1)^2*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1')^2)/((`root')^(5/2));

		qui gen double `C11s' = (0.25*(`theta1'-1)*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1'))/((`root')^(3/2))
			+(0.25*(-`theta1'-1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1'))/((`root')^(3/2))
			+(0.25*(2*(`theta1'-1)^2-4*(`theta1'-1)*`theta1')*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1')*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1'))/((`root')^(5/2));
		qui gen double `C1ss' =  (0.5*(-`theta1'-1)*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1'))/((`root')^(3/2))
			+(0.5*(`theta1'-1)^2*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1')^2)/((`root')^(5/2));
		
		qui gen double `dC11dt1' = -0.5/sqrt(`root')
			+(0.25*(`theta1'-1)*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(3/2))
			+(0.25*(`F1'-`Fs')*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1'))/((`root')^(3/2))
			+(0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`r')+2*(`theta1'-1)*(`F1'+`Fs')-4*`Fs'*`theta1'-4*`Fs'*(`theta1'-1)))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*`Fs'*(`theta1'-1)*`theta1')*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(5/2));
		qui gen double `dC1sdt1' = 0.5/sqrt(`root')
			+(0.25*(-`theta1'-1)*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(3/2))
			+(0.25*(`F1'-`Fs')*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1'))/((`root')^(3/2))
			+(0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`r')+2*(`theta1'-1)*(`F1'+`Fs')-4*`theta1'*`F1'-4*(`theta1'-1)*`F1'))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`theta1'-1)*(`r')-4*(`theta1'-1)*`theta1'*`F1')*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(5/2));
		qui gen double `d2C1dt12' = (0.25*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`F1'+`Fs')^2-8*`Fs'*`F1'))/((`root')^(3/2))
			+(0.5*(`F1'-`Fs')*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1'))/((`root')^(3/2))
			-(0.375*((`theta1'-1)*(`F1'+`Fs')-2*`Fs'*`theta1'+1)*(2*(`F1'+`Fs')*(`r')-4*`Fs'*`theta1'*`F1'-4*`Fs'*(`theta1'-1)*`F1')^2)/((`root')^(5/2));

	};
	else {;	//Archimedean copulas;
		tempvar phi3_F1 phi3_Fs phi3_C dphi2dt1_F1 dphi2dt1_Fs dphi2dt1_C d2phi1dt12_F1 d2phi1dt12_C d2Cdt12;
		
		if "`copula'" == "amh" {;
			scalar `d2t1dat12' = -2*`theta1'*(`dt1dat1');
			qui gen double `phi3_F1' = 2*(`theta1'/(1-`theta1'*(1-`F1')))^3 - 2*`F1'^(-3);
			qui gen double `phi3_Fs' = 2*(`theta1'/(1-`theta1'*(1-`Fs')))^3 - 2*`Fs'^(-3);
			qui gen double `phi3_C' = 2*(`theta1'/(1-`theta1'*(1-`C')))^3 - 2*`C'^(-3);
			
			qui gen double `dphi2dt1_F1' = (2*(`F1'-1)*(`theta1')^2)/((1-(1-`F1')*`theta1')^3)-(2*`theta1')/((1-(1-`F1')*`theta1')^2);
			qui gen double `dphi2dt1_Fs' = (2*(`Fs'-1)*(`theta1')^2)/((1-(1-`Fs')*`theta1')^3)-(2*`theta1')/((1-(1-`Fs')*`theta1')^2);
			qui gen double `dphi2dt1_C' = (2*(`C'-1)*(`theta1')^2)/((1-(1-`C')*`theta1')^3)-(2*`theta1')/((1-(1-`C')*`theta1')^2);
			
			qui gen double `d2phi1dt12_F1' =(2*(`F1'-1)^2*`theta1')/((1-(1-`F1')*`theta1')^3)-(2*(`F1'-1))/((1-(1-`F1')*`theta1')^2);
			qui gen double `d2phi1dt12_C' = (2*(`C'-1)^2*`theta1')/((1-(1-`C')*`theta1')^3)-(2*(`C'-1))/((1-(1-`C')*`theta1')^2);
			
			qui gen double `d2Cdt12' = (2*(1-`F1')^2*`F1'*(1-`Fs')^2*`Fs')/((1-`theta1'*(1-`F1')*(1-`Fs'))^3) ;
		};
		else if "`copula'" == "clayton" {;
			scalar `d2t1dat12' = `theta1';
			qui gen double `phi3_F1' = -(`theta1'+1)*(`theta1'+2)*(`F1')^(-`theta1'-3);
			qui gen double `phi3_Fs' = -(`theta1'+1)*(`theta1'+2)*(`Fs')^(-`theta1'-3);
			qui gen double `phi3_C' = -(`theta1'+1)*(`theta1'+2)*(`C')^(-`theta1'-3);
			
			qui gen double `dphi2dt1_F1' = (`F1')^(-`theta1'-2) - (`theta1'+1)*ln(`F1')*(`F1')^(-`theta1'-2);
			qui gen double `dphi2dt1_Fs' = (`Fs')^(-`theta1'-2) - (`theta1'+1)*ln(`Fs')*(`Fs')^(-`theta1'-2);
			qui gen double `dphi2dt1_C' = (`C')^(-`theta1'-2) - (`theta1'+1)*ln(`C')*(`C')^(-`theta1'-2);
			
			qui gen double `d2phi1dt12_F1' = -ln(`F1')^2*(`F1')^(-`theta1'-1);
			qui gen double `d2phi1dt12_C' = -ln(`C')^2*(`C')^(-`theta1'-1);
			
			*** numerical derivative ***;
			tempname d dtheta1 ;
			tempvar PdCdt1 MdCdt1 PC MC ;
			scalar `d' = 0.5*(abs(`theta1')+ 10^(-8))/(10^8);
			scalar `dtheta1' = `theta1' + `d';
			qui gen double `PC' = (`Fs'^(-`dtheta1')+`F1'^(-`dtheta1')-1)^(-1/`dtheta1');
			qui gen double `PdCdt1' = `PC'*(ln(`F1'^(-`dtheta1')+`Fs'^(-`dtheta1')-1)*(`dtheta1')^(-2)+(ln(`F1')*`F1'^(-`dtheta1') + ln(`Fs')*`Fs'^(-`dtheta1'))/(`dtheta1'*(`F1'^(-`dtheta1')+`Fs'^(-`dtheta1')-1)));
			scalar `dtheta1' = `theta1' - `d';
			qui gen double  `MC' = (`Fs'^(-`dtheta1')+`F1'^(-`dtheta1')-1)^(-1/`dtheta1');
			qui gen double `MdCdt1' = `MC'*(ln(`F1'^(-`dtheta1')+`Fs'^(-`dtheta1')-1)*(`dtheta1')^(-2)+(ln(`F1')*`F1'^(-`dtheta1') + ln(`Fs')*`Fs'^(-`dtheta1'))/(`dtheta1'*(`F1'^(-`dtheta1')+`Fs'^(-`dtheta1')-1)));
			qui gen double `d2Cdt12' = (`PdCdt1'-`MdCdt1')/(2*`d');
		};	
		else if "`copula'" == "frank" {;
			scalar `d2t1dat12' = 0;
			qui gen double `phi3_F1' = (`theta1'/(1-exp(`theta1'*`F1')))^3 *exp(`theta1'*`F1')*(1+exp(`theta1'*`F1'));
			qui gen double `phi3_Fs' = (`theta1'/(1-exp(`theta1'*`Fs')))^3 *exp(`theta1'*`Fs')*(1+exp(`theta1'*`Fs'));
			qui gen double `phi3_C' = (`theta1'/(1-exp(`theta1'*`C')))^3 *exp(`theta1'*`C')*(1+exp(`theta1'*`C'));
			
			qui gen double `dphi2dt1_F1' = `theta1'*exp(`theta1'*`F1')*(2+`F1'*`theta1'+(2*`theta1'*`F1'-`theta1'*`F1'-2)*exp(`theta1'*`F1'))/((1-exp(`theta1'*`F1'))^3);
			qui gen double `dphi2dt1_Fs' = `theta1'*exp(`theta1'*`Fs')*(2+`Fs'*`theta1'+(2*`theta1'*`Fs'-`theta1'*`Fs'-2)*exp(`theta1'*`Fs'))/((1-exp(`theta1'*`Fs'))^3);
			qui gen double `dphi2dt1_C' = `theta1'*exp(`theta1'*`C')*(2+`C'*`theta1'+(2*`theta1'*`C'-`theta1'*`C'-2)*exp(`theta1'*`C'))/((1-exp(`theta1'*`C'))^3);
			
			qui gen double `d2phi1dt12_F1' = 2*`F1'^2*`theta1'*exp(2*`F1'*`theta1')/((1-exp(`theta1'*`F1'))^3) + `F1'^2 *`theta1'*exp(`theta1'*`F1')/((1-exp(`theta1'*`F1'))^2)
				+ 2*`F1'*exp(`theta1'*`F1')/((1-exp(`theta1'*`F1'))^2);
			qui gen double `d2phi1dt12_C' = 2*`C'^2*`theta1'*exp(2*`C'*`theta1')/((1-exp(`theta1'*`C'))^3) + `C'^2 *`theta1'*exp(`theta1'*`C')/((1-exp(`theta1'*`C'))^2)
				+ 2*`C'*exp(`theta1'*`C')/((1-exp(`theta1'*`C'))^2);
			
			*** numerical derivatives ***;
			tempname d dtheta1;
			tempvar PdCdt1 MdCdt1 dC ;
			scalar `d' = 0.5*(abs(`theta1')+ 10^(-8))/(10^8);
			scalar `dtheta1' = `theta1' + `d';
			qui gen double `dC' = -ln(1+(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1)/(exp(-`dtheta1')-1))/`dtheta1';
			qui gen double `PdCdt1' = -`dC'/`dtheta1' + (`Fs'*exp(-`dtheta1'*`Fs')*(exp(-`dtheta1'*`F1')-1)*(exp(-`dtheta1')-1)+`F1'*exp(-`dtheta1'*`F1')*
				(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1')-1)-exp(-`dtheta1')*(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1))
				/((1+(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1)/(exp(-`dtheta1')-1))*((exp(-`dtheta1')-1)^2)*(`dtheta1'));
			scalar `dtheta1' = `theta1' - `d';
			qui replace `dC' = -ln(1+(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1)/(exp(-`dtheta1')-1))/`dtheta1';
			qui gen double `MdCdt1' = -`dC'/`dtheta1' + (`Fs'*exp(-`dtheta1'*`Fs')*(exp(-`dtheta1'*`F1')-1)*(exp(-`dtheta1')-1)+`F1'*exp(-`dtheta1'*`F1')*
				(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1')-1)-exp(-`dtheta1')*(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1))
				/((1+(exp(-`dtheta1'*`Fs')-1)*(exp(-`dtheta1'*`F1')-1)/(exp(-`dtheta1')-1))*((exp(-`dtheta1')-1)^2)*(`dtheta1'));
			qui gen double `d2Cdt12' = (`PdCdt1'-`MdCdt1')/(2*`d');	
		};
		else if "`copula'" == "gumbel" {;
			scalar `d2t1dat12' = exp(`atheta1');
			qui gen double `phi3_F1' = ((-ln(`F1'))^`theta1'*`theta1'^3+(-3*ln(`F1')-3)*(-ln(`F1'))^`theta1'*`theta1'^2+(-log(`F1'))^`theta1'*(2*log(`F1')^2+3*log(`F1')+2)*`theta1')/(`F1'^3*ln(`F1')^3);
			qui gen double `phi3_Fs' = ((-ln(`Fs'))^`theta1'*`theta1'^3+(-3*ln(`Fs')-3)*(-ln(`Fs'))^`theta1'*`theta1'^2+(-log(`Fs'))^`theta1'*(2*log(`Fs')^2+3*log(`Fs')+2)*`theta1')/(`Fs'^3*ln(`Fs')^3);
			qui gen double `phi3_C' = ((-ln(`C'))^`theta1'*`theta1'^3+(-3*ln(`C')-3)*(-ln(`C'))^`theta1'*`theta1'^2+(-log(`C'))^`theta1'*(2*log(`C')^2+3*log(`C')+2)*`theta1')/(`C'^3*ln(`C')^3);
			
			qui gen double `dphi2dt1_F1' = ((-ln(`F1'))^`theta1'*ln(-ln(`F1'))*`theta1'^2+(-ln(`F1'))^`theta1'*((-ln(`F1')-1)*ln(-ln(`F1'))+2)*`theta1'+(-ln(`F1')-1)*(-ln(`F1'))^`theta1')/(`F1'^2*ln(`F1')^2);
			qui gen double `dphi2dt1_Fs' = ((-ln(`Fs'))^`theta1'*ln(-ln(`Fs'))*`theta1'^2+(-ln(`Fs'))^`theta1'*((-ln(`Fs')-1)*ln(-ln(`Fs'))+2)*`theta1'+(-ln(`Fs')-1)*(-ln(`Fs'))^`theta1')/(`Fs'^2*ln(`Fs')^2);
			qui gen double `dphi2dt1_C' = ((-ln(`C'))^`theta1'*ln(-ln(`C'))*`theta1'^2+(-ln(`C'))^`theta1'*((-ln(`C')-1)*ln(-ln(`C'))+2)*`theta1'+(-ln(`C')-1)*(-ln(`C'))^`theta1')/(`C'^2*ln(`C')^2);
			
			qui gen double `d2phi1dt12_F1' = -((-ln(`F1'))^(`theta1'-1)*ln(-ln(`F1'))^2*`theta1')/`F1'-(2*(-ln(`F1'))^(`theta1'-1)*ln(-ln(`F1')))/`F1';
			qui gen double `d2phi1dt12_C' = -((-ln(`C'))^(`theta1'-1)*ln(-ln(`C'))^2*`theta1')/`C'-(2*(-ln(`C'))^(`theta1'-1)*ln(-ln(`C')))/`C';
			
			*** numerica derivative ***;
			tempname d dtheta1 ;
			tempvar PdCdt1 MdCdt1 dFs1 dC;
			scalar `d' = 0.5*(abs(`theta1')+ 10^(-8))/(10^8);
			scalar `dtheta1' = `theta1' + `d';
			qui gen double `dC' = exp(-((-ln(`F1'))^(`dtheta1')+(-ln(`Fs'))^(`dtheta1'))^(1/`dtheta1'));
			qui gen double `dFs1' = ((-ln(`Fs'))^`dtheta1' + (-ln(`F1'))^`dtheta1');
			qui gen double `PdCdt1' = `dC'*((ln(`dFs1')/(`dtheta1'^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`dtheta1'+
				ln(-ln(`F1'))*(-ln(`F1'))^`dtheta1')/(`dFs1'*`dtheta1'))*`dFs1'^(1/`dtheta1'));
			scalar `dtheta1' = `theta1' - `d';
			qui replace `dC' = exp(-((-ln(`F1'))^(`dtheta1')+(-ln(`Fs'))^(`dtheta1'))^(1/`dtheta1'));
			qui replace `dFs1' = ((-ln(`Fs'))^`dtheta1' + (-ln(`F1'))^`dtheta1');
			qui gen double `MdCdt1' = `dC'*((ln(`dFs1')/(`dtheta1'^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`dtheta1'+
				ln(-ln(`F1'))*(-ln(`F1'))^`dtheta1')/(`dFs1'*`dtheta1'))*`dFs1'^(1/`dtheta1'));
			qui gen double `d2Cdt12' = (`PdCdt1'-`MdCdt1')/(2*`d');
		};
		else if "`copula'"=="joe" {;
			scalar `d2t1dat12' = exp(`atheta1');
			qui gen double `phi3_F1' = 2*`phi1_F1'*`phi2_F1' - (`theta1'-1)*`phi2_F1'/(1-`F1') - (`theta1'-1)*`phi1_F1'/((1-`F1')^2);
			qui gen double `phi3_Fs' = 2*`phi1_Fs'*`phi2_Fs' - (`theta1'-1)*`phi2_Fs'/(1-`Fs') - (`theta1'-1)*`phi1_Fs'/((1-`Fs')^2);
			qui gen double `phi3_C' = 2*`phi1_C'*`phi2_C' - (`theta1'-1)*`phi2_C'/(1-`C') - (`theta1'-1)*`phi1_C'/((1-`C')^2);
			
			tempvar phi3;
			qui gen double `phi3' = -(2*(1-`F1')^(3*`theta1'-3)*`theta1'^3)/((1-(1-`F1')^`theta1')^3)-((1-`F1')^(2*`theta1'-3)*(`theta1'-1)*`theta1'^2)/((1-(1-`F1')^`theta1')^2)
				-((1-`F1')^(2*`theta1'-3)*`theta1'^2*(2*`theta1'-2))/((1-(1-`F1')^`theta1')^2)-((1-`F1')^(`theta1'-3)*(`theta1'-2)*(`theta1'-1)*`theta1')/(1-(1-`F1')^`theta1');
			
			qui gen double `dphi2dt1_F1' = 2*`phi1_F1'*`dphi1dt1_F1' - `phi1_F1'/(1-`F1') - (`theta1'-1)*`dphi1dt1_F1'/(1-`F1');
			qui gen double `dphi2dt1_Fs' = 2*`phi1_Fs'*`dphi1dt1_Fs' - `phi1_Fs'/(1-`Fs') - (`theta1'-1)*`dphi1dt1_Fs'/(1-`Fs');
			qui gen double `dphi2dt1_C' = 2*`phi1_C'*`dphi1dt1_C' - `phi1_C'/(1-`C') - (`theta1'-1)*`dphi1dt1_C'/(1-`C');
			
			tempvar dphi2;
			qui gen double `dphi2' = -((log(1-`F1')*(1-`F1')^(2*`theta1')+log(1-`F1')*(1-`F1')^`theta1')*`theta1'^2+((log(1-`F1')-2)*(1-`F1')^(2*`theta1')+(2-log(1-`F1'))*(1-`F1')^`theta1')*`theta1'-(1-`F1')^(3*`theta1')+2*(1-`F1')^(2*`theta1')-(1-`F1')^`theta1')
				/(((1-`F1')^(3*`theta1')-3*(1-`F1')^(2*`theta1')+3*(1-`F1')^`theta1'-1)*`F1'^2+(-2*(1-`F1')^(3*`theta1')+6*(1-`F1')^(2*`theta1')-6*(1-`F1')^`theta1'+2)*`F1'+(1-`F1')^(3*`theta1')-3*(1-`F1')^(2*`theta1')+3*(1-`F1')^`theta1'-1);

			qui gen double `d2phi1dt12_F1' = -(2*ln(1-`F1')^2*(1-`F1')^(3*`theta1'-1)*`theta1')/((1-(1-`F1')^`theta1')^3)-(3*ln(1-`F1')^2*(1-`F1')^(2*`theta1'-1)*`theta1')/((1-(1-`F1')^`theta1')^2)-(ln(1-`F1')^2*(1-`F1')^(`theta1'-1)*`theta1')/(1-(1-`F1')^`theta1')
				-(2*ln(1-`F1')*(1-`F1')^(2*`theta1'-1))/((1-(1-`F1')^`theta1'))^2-(2*ln(1-`F1')*(1-`F1')^(`theta1'-1))/(1-(1-`F1')^`theta1');
			qui gen double `d2phi1dt12_C' =  -(2*ln(1-`C')^2*(1-`C')^(3*`theta1'-1)*`theta1')/((1-(1-`C')^`theta1')^3)-(3*ln(1-`C')^2*(1-`C')^(2*`theta1'-1)*`theta1')/((1-(1-`C')^`theta1')^2)-(ln(1-`C')^2*(1-`C')^(`theta1'-1)*`theta1')/(1-(1-`C')^`theta1')
				-(2*ln(1-`C')*(1-`C')^(2*`theta1'-1))/((1-(1-`C')^`theta1'))^2-(2*ln(1-`C')*(1-`C')^(`theta1'-1))/(1-(1-`C')^`theta1');

			*** numerica derivative ***;
			tempname d dtheta1 ;
			tempvar PdCdt1 MdCdt1 dFs1 dC;
			scalar `d' = 0.5*(abs(`theta1')+ 10^(-8))/(10^8);
			scalar `dtheta1' = `theta1' + `d';
			qui gen double `dC' = 1-((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1')^(1/`dtheta1');
			qui gen double `PdCdt1' = (`dC'-1)*((ln(1-`Fs')*(1-`Fs')^`dtheta1' + ln(1-`F1')*(1-`F1')^`dtheta1' - ln((1-`Fs')*(1-`F1'))*((1-`Fs')*(1-`F1'))^`dtheta1')
				/((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1') /`dtheta1'-ln((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1')/(`dtheta1'^2));	
			scalar `dtheta1' = `theta1' - `d';
			qui replace `dC' = 1-((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1')^(1/`dtheta1');
			qui gen double `MdCdt1' = (`dC'-1)*((ln(1-`Fs')*(1-`Fs')^`dtheta1' + ln(1-`F1')*(1-`F1')^`dtheta1' - ln((1-`Fs')*(1-`F1'))*((1-`Fs')*(1-`F1'))^`dtheta1')
				/((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1') /`dtheta1'-ln((1-`Fs')^`dtheta1'+(1-`F1')^`dtheta1'-((1-`Fs')*(1-`F1'))^`dtheta1')/(`dtheta1'^2));	
			qui gen double `d2Cdt12' = (`PdCdt1'-`MdCdt1')/(2*`d');
		};
		
		*** general to Archimedean copula ***;
		qui gen double `C111' = `phi3_F1'/`phi1_C' - `phi1_F1'*`phi2_F1'*`phi2_C'/((`phi1_C')^3) - 2*`phi1_F1'*`phi2_F1'*`phi2_C'/((`phi1_C')^3)
			+ 3*(`phi1_F1')^3*(`phi2_C')^2 /((`phi1_C')^5) - (`phi1_F1')^3*`phi3_C'/((`phi1_C')^4);
		qui gen double `C11s' = -`phi1_Fs'*`phi2_F1'*`phi2_C'/((`phi1_C')^3) - (`phi1_F1')^2*`phi1_Fs'*`phi3_C'/((`phi1_C')^4)
			+ 3*(`phi1_F1'*`phi2_C')^2*`phi1_Fs'/((`phi1_C')^5);
		qui gen double `C1ss' = -`phi1_F1'*`phi2_Fs'*`phi2_C'/((`phi1_C')^3) - `phi1_F1'*`phi3_C'*(`phi1_Fs')^2 /((`phi1_C')^4)
			+3*`phi1_F1'*(`phi1_Fs'*`phi2_C')^2 /((`phi1_C')^5);
			
		qui gen double `dC11dt1' = `dphi2dt1_F1'/`phi1_C'
			-`phi2_F1'*(`dphi1dt1_C'+`phi2_C'*`dCdt1')/((`phi1_C')^2)
			-2*`phi1_F1'*`dphi1dt1_F1'*`phi2_C'/((`phi1_C')^3)
			-(`dphi2dt1_C' + `phi3_C'*`dCdt1')*((`phi1_F1')^2)/((`phi1_C')^3)
			+3*(`dphi1dt1_C'+`phi2_C'*`dCdt1')*`phi2_C'*((`phi1_F1')^2)/((`phi1_C')^4);
			
		qui gen double `dC1sdt1' = -(`dphi1dt1_F1'*`phi1_Fs'*`phi2_C'+`dphi1dt1_Fs'*`phi1_F1'*`phi2_C'
			+(`dphi2dt1_C' + `phi3_C'*`dCdt1')*`phi1_F1'*`phi1_Fs')/((`phi1_C')^3) 
			+3*(`phi1_F1'*`phi1_Fs'*`phi2_C')*(`dphi1dt1_C'+`phi2_C'*`dCdt1')/((`phi1_C')^4);
			
		qui gen double `d2C1dt12' = `d2phi1dt12_F1'/`phi1_C'-2*`dphi1dt1_F1'*(`dphi1dt1_C'+`phi2_C'*`dCdt1')/((`phi1_C')^2)
			+2*`phi1_F1'*(`dphi1dt1_C'+`phi2_C'*`dCdt1')^2/((`phi1_C')^3)
			-`phi1_F1'*(`d2phi1dt12_C'+`dphi2dt1_C'*`dCdt1'+(`dphi2dt1_C'+`phi3_C'*`dCdt1')*`dCdt1'+`phi2_C'*`d2Cdt12')/((`phi1_C')^2);
	};
	
	foreach x in C111 C11s C1ss dC11dt1 dC1sdt1 d2C1dt12 {;
		qui replace ``x'' = 0 if `Fs'==0 | `Fs'==1; 
	};
	
	capture {;
		foreach x in C111 C11s C1ss dC11dt1 dC1sdt1 d2C1dt12 {;
			qui replace ``x'' = 0 if `guard' == 1; 
		};
	};
	*** GENERAL ***;
	tempvar d2Ldgdg d2Ldgdb d2Ldgds d2Ldgdt d2Ldbdb d2Ldbds d2Ldbdt d2Ldsds d2Ldsdt d2Ldtdt d2Ldgdv d2Ldbdv d2Ldsdv d2Ldvdt d2Ldvdv;
	
	qui gen double `d2Ldgdg' = cond(`S'==0, -(`dFsdb'/`Fs')^2+`d2Fsdb2'/`Fs', 
		-(`dFsdb'*`C1s'/(1-`C1'))^2 - (`C1ss'/(1-`C1'))*(`dFsdb')^2- `C1s'*`d2Fsdb2'/(1-`C1') );
	qui gen double `d2Ldgdb' = `sign1'*cond(`S'==0, 0, -`C1s'*`C11'*`dFsdb'*`dF1db'/((1-`C1')^2)-`C11s'*`dFsdb'*`dF1db'/(1-`C1') );
	qui gen double `d2Ldgds' = cond(`S'==0, 0, -`C1s'*`C11'*`dFsdb'*`dF1ds'/((1-`C1')^2)-`C11s'*`dFsdb'*`dF1ds'/(1-`C1') ) ;
	qui gen double `d2Ldgdt' = cond(`S'==0, 0, -`C1s'*`dFsdb'*`dC1dt1'/((1-`C1')^2) - `dFsdb'*`dC1sdt1'/(1-`C1') );
	qui replace `d2Ldgdt' = `d2Ldgdt'*`dt1dat1'; 
	capture qui gen double `d2Ldgdv' = cond(`S'==0, 0, -`C1s'*`C11'*`dFsdb'*`dF1dv'/((1-`C1')^2)-`C11s'*`dFsdb'*`dF1dv'/(1-`C1'));
		
	qui gen double `d2Ldbdb' = cond(`S'==0, 0, -(`C11'*`dF1db'/(1-`C1'))^2 - `C111'/(1-`C1')*(`dF1db')^2 -`C11'*`d2F1db2'/(1-`C1') + `d2lnf1db2' );
	qui gen double `d2Ldbds' = `sign1'*cond(`S'==0, 0, -((`C11'/(1-`C1'))^2+`C111'/(1-`C1'))*`dF1db'*`dF1ds'-`C11'*`d2F1dbds'/(1-`C1')+ `d2lnf1dbds' );
	qui gen double `d2Ldbdt' = `sign1'*cond(`S'==0, 0, -`C11'*`dF1db'*`dC1dt1'/((1-`C1')^2) - `dF1db'*`dC11dt1'/(1-`C1') );
	qui replace `d2Ldbdt' = `d2Ldbdt'*`dt1dat1';
	capture qui gen double `d2Ldbdv' = `sign1'*cond(`S'==0, 0, -((`C11'/(1-`C1'))^2+`C111'/(1-`C1'))*`dF1db'*`dF1dv'-`C11'*`d2F1dbdv'/(1-`C1')+ `d2lnf1dbdv');
	
	qui gen double `d2Ldsds' = cond(`S'==0, 0, -(`C11'*`dF1ds'/(1-`C1'))^2 - `C111'/(1-`C1')*(`dF1ds')^2 -`C11'*`d2F1ds2'/(1-`C1') + `d2lnf1ds2' );
	qui gen double `d2Ldsdt' = cond(`S'==0, 0, -`C11'*`dF1ds'*`dC1dt1'/((1-`C1')^2) - `dF1ds'*`dC11dt1'/(1-`C1') );
	qui replace `d2Ldsdt' = `d2Ldsdt' * `dt1dat1';
	capture qui gen double `d2Ldsdv' = cond(`S'==0, 0, -((`C11'/(1-`C1'))^2+`C111'/(1-`C1'))*`dF1ds'*`dF1dv'-`C11'*`d2F1dsdv'/(1-`C1')+ `d2lnf1dsdv');
	
	qui gen double `d2Ldtdt' = cond(`S'==0, 0, -(`dC1dt1'*`dt1dat1'/(1-`C1'))^2 - `d2C1dt12'/(1-`C1')*(`dt1dat1')^2 -`dC1dt1'/(1-`C1')*`d2t1dat12');

	capture qui gen double `d2Ldvdt' = cond(`S'==0, 0, -`C11'*`dF1dv'*`dC1dt1'/((1-`C1')^2) - `dF1dv'*`dC11dt1'/(1-`C1'));
	capture qui replace `d2Ldvdt' = `d2Ldvdt'* `dt1dat1'; // ;
	
	capture qui gen double `d2Ldvdv' = cond(`S'==0, 0, -(`C11'*`dF1dv'/(1-`C1'))^2-`C111'/(1-`C1')*(`dF1dv')^2-`C11'*`d2F1dv2'/(1-`C1')+`d2lnf1dv2'); 
	
	tempname H h d11 d12 d13 d14 d15 d22 d23 d24 d25 d33 d34 d35 d44 d45 d55 ;
	
	
	tempvar touse ; mark `touse'; markout `touse' `C1';
	
	mlmatsum `lnf' `d11' = `d2Ldgdg' , eq(1);
	mlmatsum `lnf' `d12' = `d2Ldgdb' if `touse', eq(1,2);
	mlmatsum `lnf' `d13' = `d2Ldgds' if `touse', eq(1,3);
	
	mlmatsum `lnf' `d22' = `d2Ldbdb' if `touse', eq(2);
	mlmatsum `lnf' `d23' = `d2Ldbds' if `touse', eq(2,3);
	
	mlmatsum `lnf' `d33' = `d2Ldsds', eq(3);
	
	matrix `H' = (`d11', `d12', `d13' \ `d12'', `d22', `d23' \ `d13'', `d23'', `d33');
	
	local deriv_2nd dg db ds;
	local k = 4;

	if "`margin'" == "t" {;
		if "$df" == "" {;
			local deriv_2nd `deriv_2nd' dv;
			local i = 1;
			foreach j of local deriv_2nd {;
				if "`j'" != "dv" {;
					mlmatsum `lnf' `d`i'`k'' = `d2L`j'dv' if `touse', eq(`i',`k');
					if `i' == 1 {;
						matrix `h' = `d`i'`k'';
					};
					else {; matrix `h' = (`h' \ `d`i'`k''); };
					local i = `i' + 1;
				};
				else {;
					mlmatsum `lnf' `d`k'`k'' = `d2Ldvdv' if `touse', eq(`k');
				};
			};
			matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
			local k = `k' + 1;
		};
	};
	
	if "`copula'" != "product" {;
		local deriv_2nd `deriv_2nd' dt;
		local i = 1;
		foreach j of local deriv_2nd {;
			if "`j'" != "dt" {;
				mlmatsum `lnf' `d`i'`k'' = `d2L`j'dt' if `touse', eq(`i',`k');
				if `i' == 1 {;
					matrix `h' = `d`i'`k'';
				};
				else {; matrix `h' = (`h' \ `d`i'`k''); };
				local i = `i' + 1;
			};
			else {;
				mlmatsum `lnf' `d`k'`k'' = `d2Ldtdt' if `touse', eq(`k');
			};
		};
		matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
		local k = `k' + 1;
	};
	
	matrix `Hessian' = `H';
	
end;
	
*** Analytical Hessian version 1.1.0;
*** minor change in version 1.1.1;
*** allow negative dependence varsion 1.2.1;
