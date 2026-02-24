*! version1.1.0 March 2012 T.Hasebe

*** linear regression model with t distribution ***

#delimit ;
program define ml_reg_t ;

	local gradient g1 g2 ;
	if "$df" == "" local gradient `gradient' g3;
	args todo b lnf `gradient' H ;
	
	#delimit ;
	tempvar fi xb ;
	tempvar lnsig sig p_v v;
	
	mleval `xb' = `b', eq(1);
	mleval `lnsig' = `b', eq(2) scalar; scalar `sig' = exp(`lnsig');
	
	if "$df" == "" {;	
		mleval `p_v' = `b', eq(3) scalar; scalar `v' = exp(`p_v');
	};
	else scalar `v' = $df;
	
	tempvar y e;
	qui gen double `y' = $ML_y1;
	qui gen double `e' = (`y'-`xb')/`sig';
	
	
	qui replace `lnf' = lngamma((`v'+1)/2)-0.5*ln(`v'*_pi)-lngamma(`v'/2)-`lnsig'
		-((`v'+1)/2)*ln(1+(`e')^2 /`v') ;
	*qui replace `lnf' = ln(tden(`v',`e'))-`lnsig';
	
	if (`todo' ==0) exit;
	
	*** Gradient *** ;
	tempvar dlnf1db dlnf1ds dlnf1dv;
	qui gen double `dlnf1db' = (`v'+1)*`e'/(`sig'*(`v'+(`e')^2));
	qui gen double `dlnf1ds' = (`v'+1)*(`e')^2/(`v'+(`e')^2);
	qui gen double `dlnf1dv' = 0.5*(digamma((`v'+1)/2)-1/`v'-digamma(`v'/2)-ln(1+(`e')^2/`v')+(`v'+1)*(`e'/`v')^2 /(1+(`e')^2/`v'));
	
	qui replace `g1' = `dlnf1db';
	qui replace `g2' = `dlnf1ds'-1;
	
	if "$df" == "" qui replace `g3' = `dlnf1dv'* `v';
	
	if (`todo'==1) exit ;
	
	*** Hessian *** ;
	tempvar d2lnf1db2 d2lnf1dbds d2lnf1dbdv d2lnf1ds2 d2lnf1dsdv d2lnf1dv2;
	
	qui gen double `d2lnf1db2' = (`v'+1)*((`e')^2-`v')/((`sig'*(`v'+(`e')^2))^2);
	qui gen double `d2lnf1dbds' = -2*`v'*`e'*(`v'+1)/(`sig'*(`v'+(`e')^2)^2);
	qui gen double `d2lnf1dbdv' = ((`e')^2-1)*`e'/(`sig'*(`v'+(`e')^2)^2) *`v';
	
	qui gen double `d2lnf1ds2' = -2*`v'*(`v'+1)*(`e')^2/((`v'+(`e')^2)^2);
	qui gen double `d2lnf1dsdv' = ((`e')^2-1)*(`e')^2/((`v'+(`e')^2)^2) *`v';
	
	qui gen double `d2lnf1dv2' = 0.5*digamma((`v'+1)/2)+0.25*`v'*trigamma((`v'+1)/2)
		-0.5*digamma(`v'/2)-0.25*`v'*trigamma(`v'/2)-0.5*ln(1+(`e')^2/`v')
		+0.5*(`e')^2 *(2-(`v'+1)/(`v'+(`e')^2))/(`v'+(`e')^2) ;
	qui replace `d2lnf1dv2' = `d2lnf1dv2'*`v';	
		
	tempname d11 d12 d13 d22 d23 d33;
	
	mlmatsum `lnf' `d11' = `d2lnf1db2', eq(1);
	mlmatsum `lnf' `d12' = `d2lnf1dbds', eq(1,2);
	mlmatsum `lnf' `d22' = `d2lnf1ds2', eq(2);
	
	
	matrix `H' = (`d11', `d12' \
					`d12'', `d22');
	
	if "$df" == "" {;
		mlmatsum `lnf' `d13' = `d2lnf1dbdv', eq(1,3);
		mlmatsum `lnf' `d23' = `d2lnf1dsdv', eq(2,3);
		mlmatsum `lnf' `d33' = `d2lnf1dv2', eq(3);
		
		matrix `H' = (`d11', `d12', `d13' \
						`d12'', `d22', `d23' \
						`d13'', `d23'', `d33');
	};
	
end;
