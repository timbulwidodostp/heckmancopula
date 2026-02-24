*! version 1.2.0 August 7, 2012 T.Hasebe

* version 1.2.0 August 7, 2012: allow negative dependence of Clayton, Gumbel, Joe (and others as well);

*** Copula Switching Model with normal, logistic and t-distribution as margins ***;
program define switchcopula_ml
	#delimit ;
	*** option ***;
	local margin0	"$margin0";
	local margin1	"$margin1";
	local copula0	"$copula0";
	local copula1	"$copula1";
	local margins	"$margins";
	local negative1 "$negative1";
	local negative0 "$negative0";
	
	tempvar touseo;
	qui gen byte `touseo' = $touseo;
	

	local k = 5 ;
	
	forvalue i = 0/1 {;
		if "`margin`i''"=="t" & "${df`i'}"=="" {;
			local k = `k' + 1;
		};
		
		if "`copula`i''"!="product" {;
			local k = `k' + 1;
		};
	};
	
	forvalue i = 1/`k' {;
		local gradient `gradient' g`i';
	};
	
	
	args todo b lnf `gradient' Hessian ;
	
	tempvar S z y0 y1 xb0 xb1 e0 e1;
	tempname lnsig0 lnsig1 sig0 sig1 atheta0 atheta1 theta0 theta1 lndf0 df0 lndf1 df1;
	
	mleval `z' = `b', eq(1);	//selection equation
	mleval `xb0' = `b', eq(2);	//equation 0 if sel. dummy is 0
	mleval `xb1' = `b', eq(3);	//eq. 1 if sel. dummy is 1
	mleval `lnsig0' = `b', eq(4) scalar; scalar `sig0' = exp(`lnsig0');	//scale parameter for eq. 0
	mleval `lnsig1' = `b', eq(5) scalar; scalar `sig1' = exp(`lnsig1'); //s.p. for eq.1
	
	local eq = 6;
	
	if "`margin0'" == "t" {;
		if "$df0" == "" {;
			mleval `lndf0' = `b', eq(`eq') scalar; scalar `df0' = exp(`lndf0');
			local eq = `eq'+ 1;		
		};
		else  scalar `df0' = $df0 ; 
	};
	
	if "`margin1'" == "t" {;
		if "$df1" == "" {;
			mleval `lndf1' = `b', eq(`eq') scalar; scalar `df1' = exp(`lndf1');
			local eq = `eq'+1;
		};
		else scalar `df1' = $df1 ; 
	};
	
	
	if "`copula0'" != "product" {;
		mleval `atheta0' = `b', eq(`eq') scalar;
		local eq = `eq' + 1;
	};
	
	if "`copula1'" != "product" {;
		mleval `atheta1' = `b', eq(`eq') scalar;
	};
	
	
	*** dependent variables *** ;
	qui gen double `S' = $ML_y1;	//selection indicator
	qui gen double `y0' = $ML_y2 ;
	qui gen double `y1' = $ML_y3 ;
	
	*** Standarized continuous dependent variables *** ;
	if "`negative0'"=="" {; local sign0 = 1 ; };
	else {; local sign0 = -1; };
	qui gen double `e0' = `sign0'*(`y0'-`xb0')/`sig0' ;
	
	if "`negative1'"=="" {; local sign1 = 1 ; };
	else {; local sign1 = -1; };
	qui gen double `e1' = `sign1'*(`y1'-`xb1')/`sig1' ;
	
	*** Specify Margins ***;
	tempvar Fs fs F0 f0 lnf0 F1 f1 lnf1;
	
	* selection * ;
	if "`margins'" == "probit" {;
		qui gen double `Fs' = normal(-`z');
		qui gen double `fs' = normalden(`z');
	};
	else if "`margins'" == "logit" {;
		qui gen double `Fs' = 1/(1+exp(`z'));
		qui gen double `fs' = (exp(`z')/(1+exp(`z'))^2);
	};
	
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
			/**dis "`j''" ; sum `e`j''; dis `df`j'';**/;
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
		
		qui replace `C`j'' = 0 if `Fs'==0 | `F`j''==0;
		qui replace `C`j'' = 1 if `Fs'==1 ;
		cap replace `C`j'' = 0 if `guard`j'' == 1;
		
	};	//end of forvalue;
	
	*** Log-likelihood ***;

	qui replace `lnf' = cond(`S'==0, cond(`touseo'==1, ln(`C0') + `lnf0'-`lnsig0',ln(`Fs')),
		cond(`touseo'==1,ln(1-`C1') + `lnf1'-`lnsig1',ln(1-`Fs')));
	
	
	*sum `Fs' `F1' `F0' `e1' `e0' `S' if `lnf' == .;
	
	if (`todo'==0) exit ;
	
	*********************************************************************************************************************;
	*** Gradient *** ;
	
	//margin-specific ;
	tempvar dFsdb ;
	* selection *;
	if "`margins'"=="probit" {;
		qui gen double `dFsdb' = -normalden(-`z');
	};
	else if "`margins'"=="logit" {;
		qui gen double `dFsdb' = -`fs';
	};
	
	* outcomes *;
	forvalue j = 0/1 {;
		tempvar dF`j'db dF`j'ds dF`j'dv dlnf`j'db dlnf`j'ds dlnf`j'dv;
		if "`margin`j''" == "normal" {;
			qui gen double `dF`j'db' = -`f`j''/`sig`j'';
			qui gen double `dF`j'ds' = -`f`j''*`e`j'';
			qui gen double `dlnf`j'db' = `e`j''/`sig`j'';
			qui gen double `dlnf`j'ds' = (`e`j'')^2;
		};
		else if "`margin`j''" == "logistic" {;
			qui gen double `dF`j'db' = -`f`j''/`sig`j'';
			qui gen double `dF`j'ds' = -`f`j''*`e`j'';
			qui gen double `dlnf`j'db' = 1/`sig`j'' - 2*(1-`F`j'')/`sig`j'';
			qui gen double `dlnf`j'ds' = `e`j'' - 2*(1-`F`j'')*`e`j'';
		};
		else if "`margin`j''" == "t" {;
			qui gen double `dlnf`j'db' = ((`df`j''+1)/(`df`j''+(`e`j'')^2))*`e`j''/`sig`j'';
			qui gen double `dlnf`j'ds' = ((`df`j''+1)/(`df`j''+(`e`j'')^2))*(`e`j'')^2;
			qui gen double `dlnf`j'dv' = 0.5*(digamma((`df`j''+1)/2)-1/`df`j''-digamma(`df`j''/2)-ln(1+(`e`j'')^2/`df`j'')+(`df`j''+1)*(`e`j''/`df`j'')^2 /(1+(`e`j'')^2/`df`j''));
			qui replace `dlnf`j'dv' = `dlnf`j'dv'*`df`j'';
			qui gen double `dF`j'db' = -`f`j''/`sig`j'';
			qui gen double `dF`j'ds' = -`f`j''*`e`j'';
			* numerical derivative *;
			tempname d`j' ddf`j';
			tempvar P_F`j' M_F`j';
			scalar `d`j'' = (abs(`df`j'')+10^(-8))/(10^8);
			scalar `ddf`j'' = `df`j'' + `d`j''/2;
			qui gen double `P_F`j'' = 1-ttail(`ddf`j'',`e`j'');
			scalar `ddf`j'' = `df`j'' - `d`j''/2;
			qui gen double `M_F`j'' = 1-ttail(`ddf`j'',`e`j'');
			qui gen double `dF`j'dv' = (`P_F`j''-`M_F`j'')/(`d`j'');
			qui replace `dF`j'dv' = `dF`j'dv'*`df`j'';
		};
		
	}; //end of forvalue ;
	
	
	//copula-specific ;
	tempvar C11 C1s dC1dt1 C00 C0s dC0dt0;
	tempname dt1dat1 dt0dat0;
	
	forvalue j = 0/1 {;
	if "`copula`j''" == "product" {;
		scalar `dt`j'dat`j'' = 0;
		qui gen double `C`j'`j'' = 0;
		qui gen double `C`j's' = 1;
		qui gen double `dC`j'dt`j'' = 0;
	};
	else if "`copula`j''"=="gaussian" {;
		tempvar dZ`j'dt`j';
		scalar `dt`j'dat`j'' = (`sr_t`j'')^2;
		qui gen double `C`j'`j'' = -`theta`j''*normalden(`Z`j'')/(`sr_t`j''*normalden(`v`j''));
		qui gen double `C`j's' = normalden(`Z`j'')/(`sr_t`j''*normalden(`vs'));
		qui gen double `dZ`j'dt`j'' = `theta`j''*`vs'/(`sr_t`j''^3) - `v`j''*(1-(`theta`j'')^2)^(-3/2);
		qui gen double `dC`j'dt`j'' = normalden(`Z`j'')*`dZ`j'dt`j'';
	};
	else if "`copula`j''" == "fgm" {;
		scalar `dt`j'dat`j'' = 1-(`theta`j'')^2;
		qui gen double `C`j'`j'' = -2*`theta`j''*(1-`Fs')*`Fs';
		qui gen double `C`j's' = (4*`theta`j''*`F`j''-2*`theta`j'')*`Fs' - 2*`theta`j''*`F`j''+`theta`j''+1;
		qui gen double `dC`j'dt`j'' = (1-`F`j'')*(1-`Fs')*`Fs' - `F`j''*(1-`Fs')*`Fs';
	};	
	else if "`copula`j''" == "plackett" {;
		scalar `dt`j'dat`j'' = `theta`j'';
		qui gen double `C`j'`j'' = -(0.5*(`theta`j''-1))/sqrt(`root`j'')+(0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j''))/((`root`j'')^(3/2));
		qui gen double `C`j's' = -(0.5*(-`theta`j''-1))/sqrt(`root`j'')+(0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j''))/((`root`j'')^(3/2));
		qui gen double `dC`j'dt`j'' = -(0.5*(`F`j''-`Fs'))/sqrt(`root`j'')+(0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(3/2));
	};
	else {;
		tempvar phi2_F`j' phi2_Fs_`j' phi2_C_`j' dphi1dt`j'_F`j' dphi1dt`j'_Fs dphi1dt`j'_C dCdt`j';
		
		if "`copula`j''" == "amh" {;
			scalar `dt`j'dat`j'' = (1-(`theta`j'')^2);
			qui gen double `phi2_Fs_`j'' = -(`theta`j''/(1-`theta`j''*(1-`Fs')))^2+1/(`Fs'^2);
			qui gen double `phi2_F`j'' = -(`theta`j''/(1-`theta`j''*(1-`F`j'')))^2+1/(`F`j''^2);
			qui gen double `phi2_C_`j'' = -(`theta`j''/(1-`theta`j''*(1-`C_`j'')))^2+1/(`C_`j''^2);
		
			qui gen double `dphi1dt`j'_F`j'' = (1-`theta`j''*(1-`F`j''))^(-2);
			qui gen double `dphi1dt`j'_Fs' = (1-`theta`j''*(1-`Fs'))^(-2);
			qui gen double `dphi1dt`j'_C' = (1-`theta`j''*(1-`C_`j''))^(-2);
			qui gen double `dCdt`j'' = `C_`j''*(1-`Fs')*(1-`F`j'')/(1-`theta`j''*(1-`Fs')*(1-`F`j''));
		};
		else if "`copula`j''"=="clayton" {;
			scalar `dt`j'dat`j'' = `theta`j'';
			qui gen double `phi2_Fs_`j'' = (`theta`j''+1)*`Fs'^(-`theta`j''-2);
			qui gen double `phi2_F`j'' = (`theta`j''+1)*`F`j''^(-`theta`j''-2);
			qui gen double `phi2_C_`j'' = (`theta`j''+1)*`C_`j''^(-`theta`j''-2);
			
			qui gen double `dphi1dt`j'_F`j'' = ln(`F`j'')*(`F`j'')^(-`theta`j''-1);
			qui gen double `dphi1dt`j'_Fs' = ln(`Fs')*(`Fs')^(-`theta`j''-1);
			qui gen double `dphi1dt`j'_C' =  ln(`C_`j'')*(`C_`j'')^(-`theta`j''-1);
			qui gen double `dCdt`j'' = `C_`j''*(ln(`F`j''^(-`theta`j'')+`Fs'^(-`theta`j'')-1)*(`theta`j'')^(-2)+(ln(`F`j'')*`F`j''^(-`theta`j'') + ln(`Fs')*`Fs'^(-`theta`j''))/(`theta`j''*(`F`j''^(-`theta`j'')+`Fs'^(-`theta`j'')-1)));
		};
		else if "`copula`j''" == "frank" {;
			scalar `dt`j'dat`j'' = 1;
			qui gen double `phi2_F`j'' = (`theta`j''/(1-exp(`theta`j''*`F`j'')))^2 *exp(`theta`j''*`F`j'');
			qui gen double `phi2_Fs_`j'' = (`theta`j''/(1-exp(`theta`j''*`Fs')))^2 *exp(`theta`j''*`Fs');
			qui gen double `phi2_C_`j'' = (`theta`j''/(1-exp(`theta`j''*`C_`j'')))^2 *exp(`theta`j''*`C_`j'');
			
			qui gen double `dphi1dt`j'_F`j'' = (1+(`theta`j''*`F`j''-1)*exp(`theta`j''*`F`j''))/((1-exp(`theta`j''*`F`j''))^2) ;
			qui gen double `dphi1dt`j'_Fs' = (1+(`theta`j''*`Fs'-1)*exp(`theta`j''*`Fs'))/((1-exp(`theta`j''*`Fs'))^2) ;
			qui gen double `dphi1dt`j'_C' = (1+(`theta`j''*`C_`j''-1)*exp(`theta`j''*`C_`j''))/((1-exp(`theta`j''*`C_`j''))^2) ;
			qui gen double `dCdt`j'' = -`C_`j''/`theta`j'' + (`Fs'*exp(-`theta`j''*`Fs')*(exp(-`theta`j''*`F`j'')-1)*(exp(-`theta`j'')-1)+`F`j''*exp(-`theta`j''*`F`j'')*
				(exp(-`theta`j''*`Fs')-1)*(exp(-`theta`j'')-1)-exp(-`theta`j'')*(exp(-`theta`j''*`Fs')-1)*(exp(-`theta`j''*`F`j'')-1))
				/((1+(exp(-`theta`j''*`Fs')-1)*(exp(-`theta`j''*`F`j'')-1)/(exp(-`theta`j'')-1))*((exp(-`theta`j'')-1)^2)*(`theta`j''));
		};
		else if "`copula`j''"=="gumbel" {;
			tempvar Fs`j';
			scalar `dt`j'dat`j'' = exp(`atheta`j'');
			qui gen double `phi2_F`j'' = `theta`j''*(-ln(`F`j''))^(`theta`j''-1)*((`theta`j''-1)/(-ln(`F`j''))+1)/(`F`j''^2);
			qui gen double `phi2_Fs_`j'' = `theta`j''*(-ln(`Fs'))^(`theta`j''-1) *((`theta`j''-1)/(-ln(`Fs'))+ 1)/(`Fs'^2);
			qui gen double `phi2_C_`j'' = `theta`j''*(-ln(`C_`j''))^(`theta`j''-1) *((`theta`j''-1)/(-ln(`C_`j''))+ 1)/(`C_`j''^2);
			
			qui gen double `dphi1dt`j'_F`j'' = -(-ln(`F`j''))^(`theta`j''-1)*(1+`theta`j''*ln(-ln(`F`j'')))/`F`j'';
			qui gen double `dphi1dt`j'_Fs' = -(-ln(`Fs'))^(`theta`j''-1)*(1+`theta`j''*ln(-ln(`Fs')))/`Fs';
			qui gen double `dphi1dt`j'_C' = -(-ln(`C_`j''))^(`theta`j''-1)*(1+`theta`j''*ln(-ln(`C_`j'')))/`C_`j'';
			qui gen double `Fs`j'' = ((-ln(`Fs'))^`theta`j'' + (-ln(`F`j''))^`theta`j'');
			qui gen double `dCdt`j'' = `C_`j''*((ln(`Fs`j'')/(`theta`j''^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`theta`j''+
				ln(-ln(`F`j''))*(-ln(`F`j''))^`theta`j'')/(`Fs`j''*`theta`j''))*(`Fs`j'')^(1/`theta`j''));	
		};		
		else if "`copula`j''" == "joe" {;
			scalar `dt`j'dat`j'' = exp(`atheta`j'');
			
			qui gen double `phi2_F`j'' = (`phi1_F`j'')^2 - (`theta`j''-1)*`phi1_F`j''/(1-`F`j'');
			qui gen double `phi2_Fs_`j'' = (`phi1_Fs_`j'')^2 - (`theta`j''-1)*`phi1_Fs_`j''/(1-`Fs');
			qui gen double `phi2_C_`j'' = (`phi1_C_`j'')^2 - (`theta`j''-1)*`phi1_C_`j''/(1-`C_`j'');
			
			qui gen double `dphi1dt`j'_F`j'' = -(((1-`F`j'')^(`theta`j''-1)*(1+`theta`j''*ln(1-`F`j'')))/(1-(1-`F`j'')^`theta`j'')+
				(`theta`j''*(1-`F`j'')^(2*`theta`j''-1)*ln(1-`F`j''))/((1-(1-`F`j'')^`theta`j'')^2));
			qui gen double `dphi1dt`j'_Fs' = -(((1-`Fs')^(`theta`j''-1)*(1+`theta`j''*ln(1-`Fs')))/(1-(1-`Fs')^`theta`j'') +
				(`theta`j''*(1-`Fs')^(2*`theta`j''-1)*ln(1-`Fs'))/((1-(1-`Fs')^`theta`j'')^2));
			qui gen double `dphi1dt`j'_C' = -(((1-`C_`j'')^(`theta`j''-1)*(1+`theta`j''*ln(1-`C_`j'')))/(1-(1-`C_`j'')^`theta`j'') +
				(`theta`j''*(1-`C_`j'')^(2*`theta`j''-1)*ln(1-`C_`j''))/((1-(1-`C_`j'')^`theta`j'')^2));	
			qui gen double `dCdt`j'' = (`C_`j''-1)*((ln(1-`Fs')*(1-`Fs')^`theta`j'' + ln(1-`F`j'')*(1-`F`j'')^`theta`j'' - ln((1-`Fs')*(1-`F`j''))*((1-`Fs')*(1-`F`j''))^`theta`j'')
				/((1-`Fs')^`theta`j''+(1-`F`j'')^`theta`j''-((1-`Fs')*(1-`F`j''))^`theta`j'')/`theta`j''-ln((1-`Fs')^`theta`j''+(1-`F`j'')^`theta`j''-((1-`Fs')*(1-`F`j''))^`theta`j'')/(`theta`j''^2));	
		};
		*** general to Archimedean copulas ***;
		qui gen double `C`j'`j''=`phi2_F`j''/`phi1_C_`j'' - (`phi2_C_`j''*(`phi1_F`j'')^2)/((`phi1_C_`j'')^3) ;
		qui gen double `C`j's' = -`phi2_C_`j''*`phi1_F`j''*`phi1_Fs_`j''/((`phi1_C_`j'')^3) ;
		qui gen double `dC`j'dt`j'' = `dphi1dt`j'_F`j''/`phi1_C_`j'' - `phi1_F`j''*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')/((`phi1_C_`j'')^2);
	};
	
	qui replace `C`j'`j'' = 0 if `Fs'==0 | `Fs'==1;
	qui replace `C`j's' = 0 if `Fs'==0 | `Fs'==1;
	qui replace `dC`j'dt`j'' = 0 if `Fs' ==0 | `Fs'==1;
	
	capture replace `C`j'`j'' = 0 if `guard`j''==1;
	capture replace `C`j's' = 0 if `guard`j''==1;
	capture replace `dC`j'dt`j'' = 0 if `guard`j''==1;
	};	//end of forvalue 
	
	
	*** General *** ;
	qui replace `g1' = cond(`S'==0, cond(`touseo'==1, `C0s'*`dFsdb'/`C0',`dFsdb'/`Fs'), cond(`touseo'==1,-`C1s'*`dFsdb'/(1-`C1'),-`dFsdb'/(1-`Fs')));
	qui replace `g2' = `sign0'*cond(`S'==0, cond(`touseo'==1, `C00'*`dF0db'/`C0'+`dlnf0db', 0),0);
	qui replace `g3' = `sign1'*cond(`S'==0, 0, cond(`touseo'==1, -`C11'*`dF1db'/(1-`C1')+`dlnf1db', 0));
	qui replace `g4' = cond(`S'==0, cond(`touseo'==1, `C00'*`dF0ds'/`C0'+`dlnf0ds'-1, 0),0);
	qui replace `g5' = cond(`S'==0, 0, cond(`touseo'==1, -`C11'*`dF1ds'/(1-`C1')+`dlnf1ds'-1, 0));
	
	local eq = 6;
	
	
	if "`margin0'" == "t" & "$df0" == "" {;
			qui replace `g`eq'' =  cond(`S'==0, cond(`touseo'==1, `C00'*`dF0dv'/`C0'+`dlnf0dv', 0),0);
			local eq = `eq' + 1;
		};
	
	if "`margin1'" == "t" & "$df1"== "" {;
		qui replace `g`eq'' =  cond(`S'==0, 0, cond(`touseo'==1, -`C11'*`dF1dv'/(1-`C1')+`dlnf1dv', 0));
		local eq = `eq' + 1;
	};
	
	
	
	if "`copula0'" != "product" {;
		qui replace `g`eq'' = cond(`S'==0, cond(`touseo'==1,`dC0dt0'/`C0',0), 0);
		qui replace `g`eq'' = `g`eq''*`dt0dat0';
		local eq = `eq' + 1;
	};
	
	if "`copula1'" != "product" {;
		qui replace `g`eq'' = cond(`S'==0, 0, cond(`touseo'==1,-`dC1dt1'/(1-`C1'),0));
		qui replace `g`eq'' = `g`eq''*`dt1dat1';
	};
	
	if `todo'==1 exit;
	
	**************************************************************************************;
	*** Hessian ***;
	* selection *;
	tempvar d2Fsdb2;
	
	if "`margins'" == "probit" {;
		qui gen double `d2Fsdb2' = `z'*`fs';
	};
	else if "`margins'" == "logit" {;
		qui gen double `d2Fsdb2' = `dFsdb'*(2*`Fs'-1);
	};
	* outcome *;
	forvalue j = 0/1 {;
		tempvar d2F`j'db2 d2F`j'dbds d2F`j'dbdv d2F`j'ds2 d2F`j'dsdv d2F`j'dv2 d2lnf`j'db2 d2lnf`j'dbds d2lnf`j'dbdv d2lnf`j'ds2 d2lnf`j'dsdv d2lnf`j'dv2 ;
		if "`margin`j''" == "normal" {;
			qui gen double `d2F`j'db2' = -`e`j''*normalden(`e`j'')/(`sig`j''^2);
			qui gen double `d2F`j'dbds' = (1-(`e`j'')^2)*normalden(`e`j'')/`sig`j'';
			qui gen double `d2F`j'ds2' = (1-(`e`j'')^2)*normalden(`e`j'')*`e`j'';
			qui gen double `d2lnf`j'db2' = -`sig`j''^(-2);
			qui gen double `d2lnf`j'dbds' = -2*`e`j''/`sig`j'';
			qui gen double `d2lnf`j'ds2' = -2*(`e`j'')^2;
		};
		if "`margin`j''" == "logistic" {;
			qui gen double `d2F`j'db2' = `dF`j'db'*(2*`F`j''-1)/`sig`j'';
			qui gen double `d2F`j'dbds' = `dF`j'ds'*(2*`F`j''-1)/`sig`j'' + `f`j''/`sig`j'';
			qui gen double `d2F`j'ds2' = `dF`j'ds'*(2*`F`j''-1)*`e`j'' + `f`j''*`e`j'';
			qui gen double `d2lnf`j'db2' = 2*`dF`j'db'/`sig`j'';
			qui gen double `d2lnf`j'dbds' = -1/`sig`j''+2*`dF`j'ds'/`sig`j''+2*(1-`F`j'')/`sig`j'';
			qui gen double `d2lnf`j'ds2' = -`e`j''+2*`dF`j'ds'*`e`j''+2*(1-`F`j'')*`e`j'';
		};
		if "`margin`j''" == "t" {;
			qui gen double `d2lnf`j'db2' = (`df`j''+1)*((`e`j'')^2-`df`j'')/((`sig`j''*(`df`j''+(`e`j'')^2))^2);
			qui gen double `d2lnf`j'dbds' = -2*`df`j''*`e`j''*(`df`j''+1)/(`sig`j''*(`df`j''+(`e`j'')^2)^2);
			qui gen double `d2lnf`j'dbdv' = ((`e`j'')^2-1)*`e`j''/(`sig`j''*(`df`j''+(`e`j'')^2)^2)*`df`j'';
			qui gen double `d2lnf`j'ds2' = -2*`df`j''*(`df`j''+1)*(`e`j'')^2/((`df`j''+(`e`j'')^2)^2);
			qui gen double `d2lnf`j'dsdv' =  ((`e`j'')^2-1)*(`e`j'')^2/((`df`j''+(`e`j'')^2)^2) *`df`j'';
			qui gen double `d2lnf`j'dv2' = 0.5*digamma((`df`j''+1)/2)+0.25*`df`j''*trigamma((`df`j''+1)/2)
				-0.5*digamma(`df`j''/2)-0.25*`df`j''*trigamma(`df`j''/2)-0.5*ln(1+(`e`j'')^2/`df`j'')
				+0.5*(`e`j'')^2 *(2-(`df`j''+1)/(`df`j''+(`e`j'')^2))/(`df`j''+(`e`j'')^2);
			qui replace `d2lnf`j'dv2' = `d2lnf`j'dv2' * `df`j'';
			
			qui gen double `d2F`j'db2' = -`dlnf`j'db'*`f`j''/`sig`j'';
			qui gen double `d2F`j'dbds' = -`dlnf`j'ds'*`f`j''/`sig`j'' + `f`j''/`sig`j'';
			qui gen double `d2F`j'dbdv' = -`dlnf`j'dv'*`f`j''/`sig`j'';
			qui gen double `d2F`j'ds2' =  -`dlnf`j'ds'*`f`j''*`e`j'' + `f`j''*`e`j'';
			qui gen double `d2F`j'dsdv' = -`dlnf`j'dv'*`f`j''*`e`j'';
			
			*** numerica derivative ***;
			* numerical derivative *;
			scalar `d`j'' = (abs(`df`j'')+10^(-5))/(10^5);
			scalar `ddf`j'' = `df`j'' + `d`j'';
			qui replace `P_F`j'' = 1-ttail(`ddf`j'',`e`j'');
			scalar `ddf`j'' = `df`j'' - `d`j'';
			qui replace `M_F`j'' = 1-ttail(`ddf`j'',`e`j'');
			qui gen double `d2F`j'dv2' = (`P_F`j''+`M_F`j''-2*`F`j'')/(`d`j''^2);
			qui replace `d2F`j'dv2' = `d2F`j'dv2'*(`df`j'')^2 + `dF`j'dv';
			
		};

	};	//end of forvalue ;

	* copula specific *;
	forvalue j = 0/1 {;
		tempvar C`j'`j'`j' C`j'`j's C`j'ss d2C`j'dt`j'2 dC`j'`j'dt`j' dC`j'sdt`j';
		tempname d2t`j'dat`j'2;
		if "`copula`j''" == "product" {;
			scalar `d2t`j'dat`j'2' = 0;
			qui gen double `C`j'`j'`j'' = 0; qui gen double `C`j'`j's' = 0;
			qui gen double `C`j'ss' = 0; qui gen double `d2C`j'dt`j'2' = 0;
			qui gen double `dC`j'`j'dt`j'' =0; qui gen double `dC`j'sdt`j'' = 0;
		};
		else if "`copula`j''" == "gaussian" {;
			tempvar dZ`j'dt`j' d2Z`j'dt`j'2 ;
			
			scalar `d2t`j'dat`j'2' = -2*`theta`j''*(`dt`j'dat`j'');
			qui gen double `C`j'`j'`j'' = -(`theta`j'')^2*`Z`j''*normalden(`Z`j'')/((`sr_t`j''*normalden(`v`j''))^2)
				-`theta`j''*normalden(`Z`j'')*`v`j''/(`sr_t`j''*normalden(`v`j'')^2);
			qui gen double `C`j'`j's' = `theta`j''*`Z`j''*normalden(`Z`j'')/(`sr_t`j''^2 * normalden(`v`j'')*normalden(`vs'));	
			qui gen double `C`j'ss' = -`Z`j''*normalden(`Z`j'')/((`sr_t`j''*normalden(`vs'))^2) + normalden(`Z`j'')*`vs'/(`sr_t`j''*normalden(`vs')^2 );
			qui gen double `dZ`j'dt`j'' = (`theta`j''*`vs'-`v`j'')/(`sr_t`j''^3);
			qui gen double `d2Z`j'dt`j'2' = 3*`theta`j''*(`theta`j''*`vs'-`v`j'')/(`sr_t`j''^5) + `vs'/(`sr_t`j''^3);

			qui gen double `dC`j'`j'dt`j'' = (normalden(`Z`j'')/normalden(`v`j''))*(-(`sr_t`j'')^(-3) + `theta`j''*`Z`j''*`dZ`j'dt`j''/`sr_t`j'');
			qui gen double `dC`j'sdt`j'' = (normalden(`Z`j'')/normalden(`vs'))*(`theta`j''/(`sr_t`j''^3) -`Z`j''* `dZ`j'dt`j''/`sr_t`j'');
			qui gen double `d2C`j'dt`j'2' = -`Z`j''*normalden(`Z`j'')*(`dZ`j'dt`j'')^2 + normalden(`Z`j'')*`d2Z`j'dt`j'2';
		};
		else if "`copula`j''" == "fgm" {;
			scalar `d2t`j'dat`j'2' = -2*`theta`j''*(`dt`j'dat`j'');
			qui gen double `C`j'`j'`j'' = 0; 
			qui gen double `C`j'`j's' = 4*`theta`j''*`Fs'-2*`theta`j'';
			qui gen double `C`j'ss' = 4*`theta`j''*`F`j''-2*`theta`j'';
			qui gen double `dC`j'`j'dt`j'' = -2*`Fs'*(1-`Fs');
			qui gen double `dC`j'sdt`j'' = (4*`F`j''-2)*`Fs' - 2*`F`j'' + 1;
			qui gen double `d2C`j'dt`j'2' = 0;
		};
		else if "`copula`j''" == "plackett" {;
			scalar `d2t`j'dat`j'2' = `theta`j'';
			qui gen double `C`j'`j'`j'' = (0.5*(`theta`j''-1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j''))/((`root`j'')^(3/2))
				+(0.5*(`theta`j''-1)^2*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j'')^2)/((`root`j'')^(5/2));

			qui gen double `C`j'`j's' = (0.25*(`theta`j''-1)*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j''))/((`root`j'')^(3/2))
				+(0.25*(-`theta`j''-1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j''))/((`root`j'')^(3/2))
				+(0.25*(2*(`theta`j''-1)^2-4*(`theta`j''-1)*`theta`j'')*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j'')*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j''))/((`root`j'')^(5/2));
			qui gen double `C`j'ss' =  (0.5*(-`theta`j''-1)*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j''))/((`root`j'')^(3/2))
				+(0.5*(`theta`j''-1)^2*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j'')^2)/((`root`j'')^(5/2));
			
			qui gen double `dC`j'`j'dt`j'' = -0.5/sqrt(`root`j'')
				+(0.25*(`theta`j''-1)*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(3/2))
				+(0.25*(`F`j''-`Fs')*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j''))/((`root`j'')^(3/2))
				+(0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`r`j'')+2*(`theta`j''-1)*(`F`j''+`Fs')-4*`Fs'*`theta`j''-4*`Fs'*(`theta`j''-1)))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*`Fs'*(`theta`j''-1)*`theta`j'')*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(5/2));
			qui gen double `dC`j'sdt`j'' = 0.5/sqrt(`root`j'')
				+(0.25*(-`theta`j''-1)*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(3/2))
				+(0.25*(`F`j''-`Fs')*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j''))/((`root`j'')^(3/2))
				+(0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`r`j'')+2*(`theta`j''-1)*(`F`j''+`Fs')-4*`theta`j''*`F`j''-4*(`theta`j''-1)*`F`j''))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`theta`j''-1)*(`r`j'')-4*(`theta`j''-1)*`theta`j''*`F`j'')*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(5/2));
			qui gen double `d2C`j'dt`j'2' = (0.25*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`F`j''+`Fs')^2-8*`Fs'*`F`j''))/((`root`j'')^(3/2))
				+(0.5*(`F`j''-`Fs')*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j''))/((`root`j'')^(3/2))
				-(0.375*((`theta`j''-1)*(`F`j''+`Fs')-2*`Fs'*`theta`j''+1)*(2*(`F`j''+`Fs')*(`r`j'')-4*`Fs'*`theta`j''*`F`j''-4*`Fs'*(`theta`j''-1)*`F`j'')^2)/((`root`j'')^(5/2));

		};
		else {;	//Archimedean copulas;
			tempvar phi3_F`j' phi3_Fs_`j' phi3_C_`j' dphi2dt`j'_F`j' dphi2dt`j'_Fs dphi2dt`j'_C d2phi1dt`j'2_F`j' d2phi1dt`j'2_C d2Cdt`j'2;
			
			if "`copula`j''" == "amh" {;
				scalar `d2t`j'dat`j'2' = -2*`theta`j''*(`dt`j'dat`j'');
				qui gen double `phi3_F`j'' = 2*(`theta`j''/(1-`theta`j''*(1-`F`j'')))^3 - 2*`F`j''^(-3);
				qui gen double `phi3_Fs_`j'' = 2*(`theta`j''/(1-`theta`j''*(1-`Fs')))^3 - 2*`Fs'^(-3);
				qui gen double `phi3_C_`j'' = 2*(`theta`j''/(1-`theta`j''*(1-`C_`j'')))^3 - 2*(`C_`j'')^(-3);
				
				qui gen double `dphi2dt`j'_F`j'' = (2*(`F`j''-1)*(`theta`j'')^2)/((1-(1-`F`j'')*`theta`j'')^3)-(2*`theta`j'')/((1-(1-`F`j'')*`theta`j'')^2);
				qui gen double `dphi2dt`j'_Fs' = (2*(`Fs'-1)*(`theta`j'')^2)/((1-(1-`Fs')*`theta`j'')^3)-(2*`theta`j'')/((1-(1-`Fs')*`theta`j'')^2);
				qui gen double `dphi2dt`j'_C' = (2*(`C_`j''-1)*(`theta`j'')^2)/((1-(1-`C_`j'')*`theta`j'')^3)-(2*`theta`j'')/((1-(1-`C_`j'')*`theta`j'')^2);
				
				qui gen double `d2phi1dt`j'2_F`j'' =(2*(`F`j''-1)^2*`theta`j'')/((1-(1-`F`j'')*`theta`j'')^3)-(2*(`F`j''-1))/((1-(1-`F`j'')*`theta`j'')^2);
				qui gen double `d2phi1dt`j'2_C' = (2*(`C_`j''-1)^2*`theta`j'')/((1-(1-`C_`j'')*`theta`j'')^3)-(2*(`C_`j''-1))/((1-(1-`C_`j'')*`theta`j'')^2);
				
				qui gen double `d2Cdt`j'2' = (2*(1-`F`j'')^2*`F`j''*(1-`Fs')^2*`Fs')/((1-`theta`j''*(1-`F`j'')*(1-`Fs'))^3) ;
			};
			else if "`copula`j''" == "clayton" {;
				scalar `d2t`j'dat`j'2' = `theta`j'';
				qui gen double `phi3_F`j'' = -(`theta`j''+1)*(`theta`j''+2)*(`F`j'')^(-`theta`j''-3);
				qui gen double `phi3_Fs_`j'' = -(`theta`j''+1)*(`theta`j''+2)*(`Fs')^(-`theta`j''-3);
				qui gen double `phi3_C_`j'' = -(`theta`j''+1)*(`theta`j''+2)*(`C_`j'')^(-`theta`j''-3);
				
				qui gen double `dphi2dt`j'_F`j'' = (`F`j'')^(-`theta`j''-2) - (`theta`j''+1)*ln(`F`j'')*(`F`j'')^(-`theta`j''-2);
				qui gen double `dphi2dt`j'_Fs' = (`Fs')^(-`theta`j''-2) - (`theta`j''+1)*ln(`Fs')*(`Fs')^(-`theta`j''-2);
				qui gen double `dphi2dt`j'_C' = (`C_`j'')^(-`theta`j''-2) - (`theta`j''+1)*ln(`C_`j'')*(`C_`j'')^(-`theta`j''-2);
				
				qui gen double `d2phi1dt`j'2_F`j'' = -ln(`F`j'')^2*(`F`j'')^(-`theta`j''-1);
				qui gen double `d2phi1dt`j'2_C' = -ln(`C_`j'')^2*(`C_`j'')^(-`theta`j''-1);
				
				*** numerical derivative ***;
				tempname d`j' dtheta`j' ;
				tempvar PdCdt`j' MdCdt`j' PC`j' MC`j' ;
				scalar `d`j'' = 0.5*(abs(`theta`j'')+ 10^(-8))/(10^8);
				scalar `dtheta`j'' = `theta`j'' + `d`j'';
				qui gen double `PC`j'' = (`Fs'^(-`dtheta`j'')+`F`j''^(-`dtheta`j'')-1)^(-1/`dtheta`j'');
				qui gen double `PdCdt`j'' = `PC`j''*(ln(`F`j''^(-`dtheta`j'')+`Fs'^(-`dtheta`j'')-1)*(`dtheta`j'')^(-2)+(ln(`F`j'')*`F`j''^(-`dtheta`j'') + ln(`Fs')*`Fs'^(-`dtheta`j''))/(`dtheta`j''*(`F`j''^(-`dtheta`j'')+`Fs'^(-`dtheta`j'')-1)));
				scalar `dtheta`j'' = `theta`j'' - `d`j'';
				qui gen double  `MC`j'' = (`Fs'^(-`dtheta`j'')+`F`j''^(-`dtheta`j'')-1)^(-1/`dtheta`j'');
				qui gen double `MdCdt`j'' = `MC`j''*(ln(`F`j''^(-`dtheta`j'')+`Fs'^(-`dtheta`j'')-1)*(`dtheta`j'')^(-2)+(ln(`F`j'')*`F`j''^(-`dtheta`j'') + ln(`Fs')*`Fs'^(-`dtheta`j''))/(`dtheta`j''*(`F`j''^(-`dtheta`j'')+`Fs'^(-`dtheta`j'')-1)));
				qui gen double `d2Cdt`j'2' = (`PdCdt`j''-`MdCdt`j'')/(2*`d`j'');
			};	
			else if "`copula`j''" == "frank" {;
				scalar `d2t`j'dat`j'2' = 0;
				qui gen double `phi3_F`j'' = (`theta`j''/(1-exp(`theta`j''*`F`j'')))^3 *exp(`theta`j''*`F`j'')*(1+exp(`theta`j''*`F`j''));
				qui gen double `phi3_Fs_`j'' = (`theta`j''/(1-exp(`theta`j''*`Fs')))^3 *exp(`theta`j''*`Fs')*(1+exp(`theta`j''*`Fs'));
				qui gen double `phi3_C_`j'' = (`theta`j''/(1-exp(`theta`j''*`C_`j'')))^3 *exp(`theta`j''*`C_`j'')*(1+exp(`theta`j''*`C_`j''));
				
				qui gen double `dphi2dt`j'_F`j'' = `theta`j''*exp(`theta`j''*`F`j'')*(2+`F`j''*`theta`j''+(2*`theta`j''*`F`j''-`theta`j''*`F`j''-2)*exp(`theta`j''*`F`j''))/((1-exp(`theta`j''*`F`j''))^3);
				qui gen double `dphi2dt`j'_Fs' = `theta`j''*exp(`theta`j''*`Fs')*(2+`Fs'*`theta`j''+(2*`theta`j''*`Fs'-`theta`j''*`Fs'-2)*exp(`theta`j''*`Fs'))/((1-exp(`theta`j''*`Fs'))^3);
				qui gen double `dphi2dt`j'_C' = `theta`j''*exp(`theta`j''*`C_`j'')*(2+`C_`j''*`theta`j''+(2*`theta`j''*`C_`j''-`theta`j''*`C_`j''-2)*exp(`theta`j''*`C_`j''))/((1-exp(`theta`j''*`C_`j''))^3);
				
				qui gen double `d2phi1dt`j'2_F`j'' = 2*`F`j''^2*`theta`j''*exp(2*`F`j''*`theta`j'')/((1-exp(`theta`j''*`F`j''))^3) + `F`j''^2 *`theta`j''*exp(`theta`j''*`F`j'')/((1-exp(`theta`j''*`F`j''))^2)
					+ 2*`F`j''*exp(`theta`j''*`F`j'')/((1-exp(`theta`j''*`F`j''))^2);
				qui gen double `d2phi1dt`j'2_C' = 2*`C_`j''^2*`theta`j''*exp(2*`C_`j''*`theta`j'')/((1-exp(`theta`j''*`C_`j''))^3) + `C_`j''^2 *`theta`j''*exp(`theta`j''*`C_`j'')/((1-exp(`theta`j''*`C_`j''))^2)
					+ 2*`C_`j''*exp(`theta`j''*`C_`j'')/((1-exp(`theta`j''*`C_`j''))^2);
				
				*** numerical derivatives ***;
				tempname d`j' dtheta`j';
				tempvar PdCdt`j' MdCdt`j' dC`j' ;
				scalar `d`j'' = 0.5*(abs(`theta`j'')+ 10^(-8))/(10^8);
				scalar `dtheta`j'' = `theta`j'' + `d`j'';
				qui gen double `dC`j'' = -ln(1+(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1)/(exp(-`dtheta`j'')-1))/`dtheta`j'';
				qui gen double `PdCdt`j'' = -`dC`j''/`dtheta`j'' + (`Fs'*exp(-`dtheta`j''*`Fs')*(exp(-`dtheta`j''*`F`j'')-1)*(exp(-`dtheta`j'')-1)+`F`j''*exp(-`dtheta`j''*`F`j'')*
					(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j'')-1)-exp(-`dtheta`j'')*(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1))
					/((1+(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1)/(exp(-`dtheta`j'')-1))*((exp(-`dtheta`j'')-1)^2)*(`dtheta`j''));
				scalar `dtheta`j'' = `theta`j'' - `d`j'';
				qui replace `dC`j'' = -ln(1+(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1)/(exp(-`dtheta`j'')-1))/`dtheta`j'';
				qui gen double `MdCdt`j'' = -`dC`j''/`dtheta`j'' + (`Fs'*exp(-`dtheta`j''*`Fs')*(exp(-`dtheta`j''*`F`j'')-1)*(exp(-`dtheta`j'')-1)+`F`j''*exp(-`dtheta`j''*`F`j'')*
					(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j'')-1)-exp(-`dtheta`j'')*(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1))
					/((1+(exp(-`dtheta`j''*`Fs')-1)*(exp(-`dtheta`j''*`F`j'')-1)/(exp(-`dtheta`j'')-1))*((exp(-`dtheta`j'')-1)^2)*(`dtheta`j''));
				qui gen double `d2Cdt`j'2' = (`PdCdt`j''-`MdCdt`j'')/(2*`d`j'');	
				*sum `d2Cdt12';
			};
			else if "`copula`j''" == "gumbel" {;
				scalar `d2t`j'dat`j'2' = exp(`atheta`j'');
				
				qui gen double `phi3_F`j'' = ((-ln(`F`j''))^`theta`j''*`theta`j''^3+(-3*ln(`F`j'')-3)*(-ln(`F`j''))^`theta`j''*`theta`j''^2+(-log(`F`j''))^`theta`j''*(2*log(`F`j'')^2+3*log(`F`j'')+2)*`theta`j'')/(`F`j''^3*ln(`F`j'')^3);
				qui gen double `phi3_Fs_`j'' = ((-ln(`Fs'))^`theta`j''*`theta`j''^3+(-3*ln(`Fs')-3)*(-ln(`Fs'))^`theta`j''*`theta`j''^2+(-log(`Fs'))^`theta`j''*(2*log(`Fs')^2+3*log(`Fs')+2)*`theta`j'')/(`Fs'^3*ln(`Fs')^3);
				qui gen double `phi3_C_`j'' = ((-ln(`C_`j''))^`theta`j''*`theta`j''^3+(-3*ln(`C_`j'')-3)*(-ln(`C_`j''))^`theta`j''*`theta`j''^2+(-log(`C_`j''))^`theta`j''*(2*log(`C_`j'')^2+3*log(`C_`j'')+2)*`theta`j'')/(`C_`j''^3*ln(`C_`j'')^3);
				
				qui gen double `dphi2dt`j'_F`j'' = ((-ln(`F`j''))^`theta`j''*ln(-ln(`F`j''))*`theta`j''^2+(-ln(`F`j''))^`theta`j''*((-ln(`F`j'')-1)*ln(-ln(`F`j''))+2)*`theta`j''+(-ln(`F`j'')-1)*(-ln(`F`j''))^`theta`j'')/(`F`j''^2*ln(`F`j'')^2);
				qui gen double `dphi2dt`j'_Fs' = ((-ln(`Fs'))^`theta`j''*ln(-ln(`Fs'))*`theta`j''^2+(-ln(`Fs'))^`theta`j''*((-ln(`Fs')-1)*ln(-ln(`Fs'))+2)*`theta`j''+(-ln(`Fs')-1)*(-ln(`Fs'))^`theta`j'')/(`Fs'^2*ln(`Fs')^2);
				qui gen double `dphi2dt`j'_C' = ((-ln(`C_`j''))^`theta`j''*ln(-ln(`C_`j''))*`theta`j''^2+(-ln(`C_`j''))^`theta`j''*((-ln(`C_`j'')-1)*ln(-ln(`C_`j''))+2)*`theta`j''+(-ln(`C_`j'')-1)*(-ln(`C_`j''))^`theta`j'')/(`C_`j''^2*ln(`C_`j'')^2);
				
				qui gen double `d2phi1dt`j'2_F`j'' = -((-ln(`F`j''))^(`theta`j''-1)*ln(-ln(`F`j''))^2*`theta`j'')/`F`j''-(2*(-ln(`F`j''))^(`theta`j''-1)*ln(-ln(`F`j'')))/`F`j'';
				qui gen double `d2phi1dt`j'2_C' = -((-ln(`C_`j''))^(`theta`j''-1)*ln(-ln(`C_`j''))^2*`theta`j'')/`C_`j''-(2*(-ln(`C_`j''))^(`theta`j''-1)*ln(-ln(`C_`j'')))/`C_`j'';
				
				*** numerica derivative ***;
				tempname d`j' dtheta`j' ;
				tempvar PdCdt`j' MdCdt`j' dFs`j' dC`j';
				scalar `d`j'' = 0.5*(abs(`theta`j'')+ 10^(-8))/(10^8);
				scalar `dtheta`j'' = `theta`j'' + `d`j'';
				qui gen double `dC`j'' = exp(-((-ln(`F`j''))^(`dtheta`j'')+(-ln(`Fs'))^(`dtheta`j''))^(1/`dtheta`j''));
				qui gen double `dFs`j'' = ((-ln(`Fs'))^`dtheta`j'' + (-ln(`F`j''))^`dtheta`j'');
				qui gen double `PdCdt`j'' = `dC`j''*((ln(`dFs`j'')/(`dtheta`j''^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`dtheta`j''+
					ln(-ln(`F`j''))*(-ln(`F`j''))^`dtheta`j'')/(`dFs`j''*`dtheta`j''))*`dFs`j''^(1/`dtheta`j''));
				scalar `dtheta`j'' = `theta`j'' - `d`j'';
				qui replace `dC`j'' = exp(-((-ln(`F`j''))^(`dtheta`j'')+(-ln(`Fs'))^(`dtheta`j''))^(1/`dtheta`j''));
				qui replace `dFs`j'' = ((-ln(`Fs'))^`dtheta`j'' + (-ln(`F`j''))^`dtheta`j'');
				qui gen double `MdCdt`j'' = `dC`j''*((ln(`dFs`j'')/(`dtheta`j''^2)-(ln(-ln(`Fs'))*(-ln(`Fs'))^`dtheta`j''+
					ln(-ln(`F`j''))*(-ln(`F`j''))^`dtheta`j'')/(`dFs`j''*`dtheta`j''))*`dFs`j''^(1/`dtheta`j''));
				qui gen double `d2Cdt`j'2' = (`PdCdt`j''-`MdCdt`j'')/(2*`d`j'');
			};
			else if "`copula`j''"=="joe" {;
				scalar `d2t`j'dat`j'2' = exp(`atheta`j'');
				qui gen double `phi3_F`j'' = 2*`phi1_F`j''*`phi2_F`j'' - (`theta`j''-1)*`phi2_F`j''/(1-`F`j'') - (`theta`j''-1)*`phi1_F`j''/((1-`F`j'')^2);
				qui gen double `phi3_Fs_`j'' = 2*`phi1_Fs_`j''*`phi2_Fs_`j'' - (`theta`j''-1)*`phi2_Fs_`j''/(1-`Fs') - (`theta`j''-1)*`phi1_Fs_`j''/((1-`Fs')^2);
				qui gen double `phi3_C_`j'' = 2*`phi1_C_`j''*`phi2_C_`j'' - (`theta`j''-1)*`phi2_C_`j''/(1-`C_`j'') - (`theta`j''-1)*`phi1_C_`j''/((1-`C_`j'')^2);
				
				qui gen double `dphi2dt`j'_F`j'' = 2*`phi1_F`j''*`dphi1dt`j'_F`j'' - `phi1_F`j''/(1-`F`j'') - (`theta`j''-1)*`dphi1dt`j'_F`j''/(1-`F`j'');
				qui gen double `dphi2dt`j'_Fs' = 2*`phi1_Fs_`j''*`dphi1dt`j'_Fs' - `phi1_Fs_`j''/(1-`Fs') - (`theta`j''-1)*`dphi1dt`j'_Fs'/(1-`Fs');
				qui gen double `dphi2dt`j'_C' = 2*`phi1_C_`j''*`dphi1dt`j'_C' - `phi1_C_`j''/(1-`C_`j'') - (`theta`j''-1)*`dphi1dt`j'_C'/(1-`C_`j'');
				
				qui gen double `d2phi1dt`j'2_F`j'' = -(2*ln(1-`F`j'')^2*(1-`F`j'')^(3*`theta`j''-1)*`theta`j'')/((1-(1-`F`j'')^`theta`j'')^3)-(3*ln(1-`F`j'')^2*(1-`F`j'')^(2*`theta`j''-1)*`theta`j'')/((1-(1-`F`j'')^`theta`j'')^2)-(ln(1-`F`j'')^2*(1-`F`j'')^(`theta`j''-1)*`theta`j'')/(1-(1-`F`j'')^`theta`j'')
					-(2*ln(1-`F`j'')*(1-`F`j'')^(2*`theta`j''-1))/((1-(1-`F`j'')^`theta`j''))^2-(2*ln(1-`F`j'')*(1-`F`j'')^(`theta`j''-1))/(1-(1-`F`j'')^`theta`j'');
				qui gen double `d2phi1dt`j'2_C' =  -(2*ln(1-`C_`j'')^2*(1-`C_`j'')^(3*`theta`j''-1)*`theta`j'')/((1-(1-`C_`j'')^`theta`j'')^3)-(3*ln(1-`C_`j'')^2*(1-`C_`j'')^(2*`theta`j''-1)*`theta`j'')/((1-(1-`C_`j'')^`theta`j'')^2)-(ln(1-`C_`j'')^2*(1-`C_`j'')^(`theta`j''-1)*`theta`j'')/(1-(1-`C_`j'')^`theta`j'')
					-(2*ln(1-`C_`j'')*(1-`C_`j'')^(2*`theta`j''-1))/((1-(1-`C_`j'')^`theta`j''))^2-(2*ln(1-`C_`j'')*(1-`C_`j'')^(`theta`j''-1))/(1-(1-`C_`j'')^`theta`j'');

				*** numerica derivative ***;
				tempname d`j' dtheta`j' ;
				tempvar PdCdt`j' MdCdt`j' dFs`j' dC`j';
				scalar `d`j'' = 0.5*(abs(`theta`j'')+ 10^(-8))/(10^8);
				scalar `dtheta`j'' = `theta`j'' + `d`j'';
				qui gen double `dC`j'' = 1-((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')^(1/`dtheta`j'');
				qui gen double `PdCdt`j'' = (`dC`j''-1)*((ln(1-`Fs')*(1-`Fs')^`dtheta`j'' + ln(1-`F`j'')*(1-`F`j'')^`dtheta`j'' - ln((1-`Fs')*(1-`F`j''))*((1-`Fs')*(1-`F`j''))^`dtheta`j'')
					/((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')/`dtheta`j''-ln((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')/(`dtheta`j''^2));	
				scalar `dtheta`j'' = `theta`j'' - `d`j'';
				qui replace `dC`j'' = 1-((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')^(1/`dtheta`j'');
				qui gen double `MdCdt`j'' = (`dC`j''-1)*((ln(1-`Fs')*(1-`Fs')^`dtheta`j'' + ln(1-`F`j'')*(1-`F`j'')^`dtheta`j'' - ln((1-`Fs')*(1-`F`j''))*((1-`Fs')*(1-`F`j''))^`dtheta`j'')
					/((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')/`dtheta`j''-ln((1-`Fs')^`dtheta`j''+(1-`F`j'')^`dtheta`j''-((1-`Fs')*(1-`F`j''))^`dtheta`j'')/(`dtheta`j''^2));	
				qui gen double `d2Cdt`j'2' = (`PdCdt`j''-`MdCdt`j'')/(2*`d`j'');
			};
			
			*** general to Archimedean copula ***;
			qui gen double `C`j'`j'`j'' = `phi3_F`j''/`phi1_C_`j'' - `phi1_F`j''*`phi2_F`j''*`phi2_C_`j''/((`phi1_C_`j'')^3) - 2*`phi1_F`j''*`phi2_F`j''*`phi2_C_`j''/((`phi1_C_`j'')^3)
				+ 3*(`phi1_F`j'')^3*(`phi2_C_`j'')^2 /((`phi1_C_`j'')^5) - (`phi1_F`j'')^3*`phi3_C_`j''/((`phi1_C_`j'')^4);
			qui gen double `C`j'`j's' = -`phi1_Fs_`j''*`phi2_F`j''*`phi2_C_`j''/((`phi1_C_`j'')^3) - (`phi1_F`j'')^2*`phi1_Fs_`j''*`phi3_C_`j''/((`phi1_C_`j'')^4)
				+ 3*(`phi1_F`j''*`phi2_C_`j'')^2*`phi1_Fs_`j''/((`phi1_C_`j'')^5);
			qui gen double `C`j'ss' = -`phi1_F`j''*`phi2_Fs_`j''*`phi2_C_`j''/((`phi1_C_`j'')^3) - `phi1_F`j''*`phi3_C_`j''*(`phi1_Fs_`j'')^2 /((`phi1_C_`j'')^4)
				+3*`phi1_F`j''*(`phi1_Fs_`j''*`phi2_C_`j'')^2 /((`phi1_C_`j'')^5);
			
			qui gen double `dC`j'`j'dt`j'' = `dphi2dt`j'_F`j''/`phi1_C_`j''
				-`phi2_F`j''*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')/((`phi1_C_`j'')^2)
				-2*`phi1_F`j''*`dphi1dt`j'_F`j''*`phi2_C_`j''/((`phi1_C_`j'')^3)
				-(`dphi2dt`j'_C' + `phi3_C_`j''*`dCdt`j'')*((`phi1_F`j'')^2)/((`phi1_C_`j'')^3)
				+3*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')*`phi2_C_`j''*((`phi1_F`j'')^2)/((`phi1_C_`j'')^4);
				
			qui gen double `dC`j'sdt`j'' = -(`dphi1dt`j'_F`j''*`phi1_Fs_`j''*`phi2_C_`j''+`dphi1dt`j'_Fs'*`phi1_F`j''*`phi2_C_`j''
				+(`dphi2dt`j'_C' + `phi3_C_`j''*`dCdt`j'')*`phi1_F`j''*`phi1_Fs_`j'')/((`phi1_C_`j'')^3) 
				+3*(`phi1_F`j''*`phi1_Fs_`j''*`phi2_C_`j'')*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')/((`phi1_C_`j'')^4);
				
			qui gen double `d2C`j'dt`j'2' = `d2phi1dt`j'2_F`j''/`phi1_C_`j''-2*`dphi1dt`j'_F`j''*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')/((`phi1_C_`j'')^2)
				+2*`phi1_F`j''*(`dphi1dt`j'_C'+`phi2_C_`j''*`dCdt`j'')^2/((`phi1_C_`j'')^3)
				-`phi1_F`j''*(`d2phi1dt`j'2_C'+`dphi2dt`j'_C'*`dCdt`j''+(`dphi2dt`j'_C'+`phi3_C_`j''*`dCdt`j'')*`dCdt`j''+`phi2_C_`j''*`d2Cdt`j'2')/((`phi1_C_`j'')^2);
		
		};	//end of Archimedean copual
		
		
		foreach x in C`j'`j'`j' C`j'`j's C`j'ss dC`j'`j'dt`j' dC`j'sdt`j' d2C`j'dt`j'2 {;
			qui replace ``x'' = 0 if `Fs'==0 | `Fs'==1; 
		};
	
		capture {;
			foreach x in C`j'`j'`j' C`j'`j's C`j'ss dC`j'`j'dt`j' dC`j'sdt`j' d2C`j'dt`j'2 {;
				qui replace ``x'' = 0 if `guard' == 1; 
			}; 
		}; //end of capture ;
	};	//end of forvalue
	
	*** GENERAL ***;
	local deriv dg db0 db1 ds0 ds1 dv0 dv1 dt0 dt1 ;

	foreach j of local deriv {;
		local rest `ferest()';
		tempvar d2L`j'`j';
		foreach k of local rest {;
			tempvar d2L`j'`k';
		};
	};

	qui gen double `d2Ldgdg' = cond(`S'==0,cond(`touseo'==1,-((`C0s'/`C0')^2-`C0ss'/`C0')*(`dFsdb')^2+(`C0s'/`C0')*`d2Fsdb2',-(`dFsdb'/`Fs')^2+`d2Fsdb2'/`Fs')
		, cond(`touseo'==1,-((`C1s'/(1-`C1'))^2+`C1ss'/(1-`C1'))*(`dFsdb')^2-(`C1s'/(1-`C1'))*`d2Fsdb2', -(`dFsdb'/(1-`Fs'))^2-`d2Fsdb2'/(1-`Fs')));
	qui gen double `d2Ldgdb0' = `sign0'*cond(`S'==0, cond(`touseo'==1,-(`C0s'*`C00'/((`C0')^2)-`C00s'/`C0')*`dFsdb'*`dF0db',0),0);	
	qui gen double `d2Ldgdb1' = `sign1'*cond(`S'==0, 0, cond(`touseo'==1,-(`C1s'*`C11'/((1-`C1')^2)+`C11s'/(1-`C1'))*`dFsdb'*`dF1db',0));
	qui gen double `d2Ldgds0' = cond(`S'==0, cond(`touseo'==1,-(`C0s'*`C00'/((`C0')^2)-`C00s'/`C0')*`dFsdb'*`dF0ds',0),0);	
	qui gen double `d2Ldgds1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C1s'*`C11'/((1-`C1')^2)+`C11s'/(1-`C1'))*`dFsdb'*`dF1ds',0));
	qui gen double `d2Ldgdt0' = cond(`S'==0, cond(`touseo'==1,-(`C0s'*`dC0dt0'/(`C0'^2)-`dC0sdt0'/`C0')*`dFsdb',0),0);
	qui replace `d2Ldgdt0' = `d2Ldgdt0'*`dt0dat0';
	qui gen double `d2Ldgdt1' = cond(`S'==0, 0, cond(`touseo'==1, -(`C1s'*`dC1dt1'/((1-`C1')^2)+`dC1sdt1'/(1-`C1'))*`dFsdb',0));
	qui replace `d2Ldgdt1' = `d2Ldgdt1'*`dt1dat1';
	capture gen double `d2Ldgdv0' = cond(`S'==0, cond(`touseo'==1,-(`C0s'*`C00'/((`C0')^2)-`C00s'/`C0')*`dFsdb'*`dF0dv',0),0);	 
	capture gen double `d2Ldgdv1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C1s'*`C11'/((1-`C1')^2)+`C11s'/(1-`C1'))*`dFsdb'*`dF1dv',0));
	
	qui gen double `d2Ldb0db0' = cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0db')^2+`C00'*`d2F0db2'/`C0'+`d2lnf0db2',0),0);
	qui gen double `d2Ldb0db1' = 0;
	qui gen double `d2Ldb0ds0' = `sign0'*cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0db')*(`dF0ds')+`C00'*`d2F0dbds'/`C0'+`d2lnf0dbds',0),0);
	qui gen double `d2Ldb0ds1' = 0;
	qui gen double `d2Ldb0dt0' = `sign0'*cond(`S'==0, cond(`touseo'==1,-(`C00'*`dC0dt0'/((`C0')^2)-`dC00dt0'/`C0')*`dF0db' ,0), 0);
	qui replace `d2Ldb0dt0' = `d2Ldb0dt0'*`dt0dat0';
	qui gen double `d2Ldb0dt1' = 0;
	capture gen double `d2Ldb0dv0' = `sign0'*cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0db')*(`dF0dv')+`C00'*`d2F0dbdv'/`C0'+`d2lnf0dbdv',0),0);
	capture gen double `d2Ldb0dv1' = 0;
	
	qui gen double `d2Ldb1db1' = cond(`S'==0,0,cond(`touseo'==1,-(`C11'*`dF1db'/(1-`C1'))^2 - `C111'/(1-`C1')*(`dF1db')^2-`C11'*`d2F1db2'/(1-`C1')+`d2lnf1db2',0));
	qui gen double `d2Ldb1ds0' = 0;
	qui gen double `d2Ldb1ds1' = `sign1'*cond(`S'==0,0,cond(`touseo'==1,-((`C11'/(1-`C1'))^2 + `C111'/(1-`C1'))*(`dF1db')*(`dF1ds')-`C11'*`d2F1dbds'/(1-`C1')+`d2lnf1dbds',0));
	qui gen double `d2Ldb1dt0' = 0;
	qui gen double `d2Ldb1dt1' = `sign1'*cond(`S'==0, 0, cond(`touseo'==1,-(`C11'*`dC1dt1'/((1-`C1')^2)+`dC11dt1'/(1-`C1'))*`dF1db',0));
	qui replace `d2Ldb1dt1' = `d2Ldb1dt1'*`dt1dat1';
	capture gen double `d2Ldb1dv0' = 0;
	capture gen double `d2Ldb1dv1' = `sign1'*cond(`S'==0,0,cond(`touseo'==1,-((`C11'/(1-`C1'))^2 + `C111'/(1-`C1'))*(`dF1db')*(`dF1dv')-`C11'*`d2F1dbdv'/(1-`C1')+`d2lnf1dbdv',0));
	
	qui gen double `d2Lds0ds0' = cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0ds')^2+`C00'*`d2F0ds2'/`C0'+`d2lnf0ds2',0),0);
	qui gen double `d2Lds0ds1' = 0;
	qui gen double `d2Lds0dt0' = cond(`S'==0, cond(`touseo'==1,-(`C00'*`dC0dt0'/((`C0')^2)-`dC00dt0'/`C0')*`dF0ds',0), 0);
	qui replace `d2Lds0dt0' = `d2Lds0dt0'*`dt0dat0';
	qui gen double `d2Lds0dt1' = 0;
	capture gen double `d2Lds0dv0' = cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0ds')*(`dF0dv')+`C00'*`d2F0dsdv'/`C0'+`d2lnf0dsdv',0),0);
	capture gen double `d2Lds0dv1' = 0;
	
	qui gen double `d2Lds1ds1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C11'*`dF1ds'/(1-`C1'))^2 - `C111'/(1-`C1')*(`dF1ds')^2-`C11'*`d2F1ds2'/(1-`C1')+`d2lnf1ds2',0));
	qui gen double `d2Lds1dt0' = 0;
	qui gen double `d2Lds1dt1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C11'*`dC1dt1'/((1-`C1')^2)+`dC11dt1'/(1-`C1'))*`dF1ds',0));
	qui replace `d2Lds1dt1' = `d2Lds1dt1'*`dt1dat1';
	capture gen double `d2Lds1dv0' = 0;
	capture gen double `d2Lds1dv1' = cond(`S'==0,0,cond(`touseo'==1,-((`C11'/(1-`C1'))^2 + `C111'/(1-`C1'))*(`dF1ds')*(`dF1dv')-`C11'*`d2F1dsdv'/(1-`C1')+`d2lnf1dsdv',0));
	
	capture gen double `d2Ldv0dv0' = cond(`S'==0, cond(`touseo'==1,-((`C00'/`C0')^2-`C000'/`C0')*(`dF0dv')^2+`C00'*`d2F0dv2'/`C0'+`d2lnf0dv2',0),0);
	capture gen double `d2Ldv0dv1' = 0;
	capture gen double `d2Ldv0dt0' =  cond(`S'==0, cond(`touseo'==1,-(`C00'*`dC0dt0'/((`C0')^2)-`dC00dt0'/`C0')*`dF0dv',0), 0);
	capture replace `d2Ldv0dt0' = `d2Ldv0dt0'*`dt0dat0';
	capture gen double `d2Ldv0dt1' = 0;
	
	
	capture gen double `d2Ldv1dv1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C11'*`dF1dv'/(1-`C1'))^2 - `C111'/(1-`C1')*(`dF1dv')^2-`C11'*`d2F1dv2'/(1-`C1')+`d2lnf1dv2',0));
	capture gen double `d2Ldv1dt0' = 0;
	capture gen double `d2Ldv1dt1' = cond(`S'==0, 0, cond(`touseo'==1,-(`C11'*`dC1dt1'/((1-`C1')^2)+`dC11dt1'/(1-`C1'))*`dF1dv',0));
	capture replace `d2Ldv1dt1' = `d2Ldv1dt1' * `dt1dat1';
	
	
	qui gen double `d2Ldt0dt0' = cond(`S'==0,cond(`touseo'==1,-(`dC0dt0'*`dt0dat0'/`C0')^2+(`d2C0dt02'/`C0')*(`dt0dat0')^2+(`dC0dt0'/`C0')*`d2t0dat02' ,0),0);
	qui gen double `d2Ldt0dt1' = 0;
	
	qui gen double `d2Ldt1dt1' = cond(`S'==0, 0, cond(`touseo'==1,-(`dC1dt1'*`dt1dat1'/(1-`C1'))^2-`d2C1dt12'/(1-`C1')*(`dt1dat1')^2-`dC1dt1'/(1-`C1')*`d2t1dat12',0));
	
	forvalue j = 1/`eq' {;
		forvalue k = 1/`eq' {;
			tempname d`j'`k';
		};
	};
	tempname h H;
	

	tempvar touse1 touse0 ; mark `touse1'; mark `touse0'; markout `touse0' `C0'; markout `touse1' `C1';
	
	mlmatsum `lnf' `d11' = `d2Ldgdg', eq(1);
	mlmatsum `lnf' `d12' = `d2Ldgdb0' if `touseo'==1 & `S'==0, eq(1,2);
	mlmatsum `lnf' `d13' = `d2Ldgdb1' if `touseo'==1 & `S'==1, eq(1,3);
	mlmatsum `lnf' `d14' = `d2Ldgds0' if `touseo'==1 & `S'==0, eq(1,4);
	mlmatsum `lnf' `d15' = `d2Ldgds1' if `touseo'==1 & `S'==1, eq(1,5);
	mlmatsum `lnf' `d22' = `d2Ldb0db0' if `touseo'==1 & `S'==0, eq(2);
	mlmatsum `lnf' `d23' = `d2Ldb0db1' if `touseo'==1 & `C0'!=. & `C1'!=., eq(2,3);
	mlmatsum `lnf' `d24' = `d2Ldb0ds0' if `touseo'==1 & `S'==0, eq(2,4);
	mlmatsum `lnf' `d25' = `d2Ldb0ds1' if `touseo'==1 & `C0'!=. & `C1'!=., eq(2,5);
	mlmatsum `lnf' `d33' = `d2Ldb1db1' if `touseo'==1 & `S'==1, eq(3);
	mlmatsum `lnf' `d34' = `d2Ldb1ds0' if `touseo'==1 & `C0'!=. & `C1'!=., eq(3,4);
	mlmatsum `lnf' `d35' = `d2Ldb1ds1' if `touseo'==1 & `S'==1, eq(3,5);
	mlmatsum `lnf' `d44' = `d2Lds0ds0' if `touseo'==1 & `S'==0, eq(4);
	mlmatsum `lnf' `d45' = `d2Lds0ds1' if `touseo'==1 & `C0'!=. & `C1'!=., eq(4,5);
	mlmatsum `lnf' `d55' = `d2Lds1ds1' if `touseo'==1 & `S'==1, eq(5);
	
	matrix `H' = (`d11', `d12', `d13', `d14', `d15' \ `d12'', `d22', `d23', `d24', `d25' \ `d13'', `d23'', `d33', `d34', `d35' 
				\ `d14'', `d24'', `d34'', `d44', `d45' \ `d15'', `d25'', `d35'', `d45'', `d55');

	local deriv_2nd dg db0 db1 ds0 ds1;
	local k = 6;
		
	if ("`margin0'" == "t") & "$df0"=="" {;
		local deriv_2nd `deriv_2nd' dv0;
		local i = 1;
		foreach j of local deriv_2nd {;
			if "`j'" != "dv0" {;
				if "`j'" == "db1" local condition `C0'!=. & `C1' !=. ;
				else local condition `S'==0 ;
				mlmatsum `lnf' `d`i'`k'' = `d2L`j'dv0' if `touseo'==1 & `condition', eq(`i',`k');
				if `i' == 1 {;
					matrix `h' = `d`i'`k'';
				};
				else {;	matrix `h' = (`h' \ `d`i'`k''); };
				local i = `i' + 1;
			};
			else {;
				mlmatsum `lnf' `d`k'`k'' = `d2Ldv0dv0' if `touseo'==1 & `S'==0, eq(`k');
			};
		};
		matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
		local k = `k'+1;
	};
	
	if ("`margin1'" == "t") & "$df1"=="" {;
		local deriv_2nd `deriv_2nd' dv1;
		local i = 1;
		foreach j of local deriv_2nd {;
			if "`j'" != "dv1" {;
				if "`j'" == "db0" local condition `C0'!=. & `C1' !=. ;
				else local condition `S'==1 ;
				mlmatsum `lnf' `d`i'`k'' = `d2L`j'dv1' if `touseo'==1 & `condition', eq(`i',`k');
				if `i' == 1 {;
					matrix `h' = `d`i'`k'';
				};
				else {;	matrix `h' = (`h' \ `d`i'`k''); };
				local i = `i' + 1;
			};
			else {;
				mlmatsum `lnf' `d`k'`k'' = `d2Ldv1dv1' if `touseo'==1 & `S'==1, eq(`k');
			};
		};
		matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
		local k = `k'+1;
	};
		
	if ("`copula0'"!="product") {;
		local deriv_2nd `deriv_2nd' dt0;
		local i = 1;
		foreach j of local deriv_2nd {;
			if "`j'" != "dt0" {;
				if "`j'" == "db1" local condition `C0'!=. & `C1' !=. ;
				else local condition `S'==0 ;
				mlmatsum `lnf' `d`i'`k'' = `d2L`j'dt0' if `touseo'==1 & `condition', eq(`i',`k');
				if `i' == 1 {;
					matrix `h' = `d`i'`k'';
				};
				else {;	matrix `h' = (`h' \ `d`i'`k''); };
				local i = `i' + 1;
			};
			else {;
				mlmatsum `lnf' `d`k'`k'' = `d2Ldt0dt0' if `touseo'==1 & `S'==0, eq(`k');
			};
		};
		matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
		local k = `k'+1;
	};
	
	if ("`copula1'"!="product") {;
		local deriv_2nd `deriv_2nd' dt1;
		local i = 1;
		foreach j of local deriv_2nd {;
			if "`j'" != "dt1" {;
				if "`j'" == "db0" local condition `C0'!=. & `C1' !=. ;
				else local condition `S'==1 ;
				mlmatsum `lnf' `d`i'`k'' = `d2L`j'dt1' if `touseo'==1 & `S'==1, eq(`i',`k');
				if `i' == 1 {;
					matrix `h' = `d`i'`k'';
				};
				else {;	matrix `h' = (`h' \ `d`i'`k''); };
				local i = `i' + 1;
			};
			else {;
				mlmatsum `lnf' `d`k'`k'' = `d2Ldt1dt1' if `touseo'==1 & `S'==1, eq(`k');
			};
		};
		matrix `H' = (`H', `h' \ `h'', `d`k'`k'');
		local k = `k'+1;
	};
	

	
	matrix `Hessian' = `H';
	
end;
	
