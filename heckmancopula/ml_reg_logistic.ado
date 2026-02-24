*! version1.0.0 Dec2011 T.Hasebe

*** linear regression model with logistic distribution ***

#delimit;
program define ml_reg_logistic ;
	args todo b lnf g1 g2 H; 
	
	tempvar fi xb ;
	tempvar lnsig sig;
	
	mleval `xb' = `b', eq(1);
	mleval `lnsig' = `b', eq(2) scalar; scalar `sig' = exp(`lnsig');
	
	tempvar y e;
	qui gen double `y' = $ML_y1;
	qui gen double `e' = (`y'-`xb')/`sig';
	
	qui replace `lnf' = -`e'-2*ln(1+exp(-`e'))-`lnsig' ;
	
	if (`todo' ==0) exit;
	
	*** Gradient *** ;
	tempvar de ;
	qui gen double `de' = (exp(-`e')-1)/(exp(-`e')+1);
	
	qui replace `g1' = `de'*(-1/`sig');
	qui replace `g2' = `de'*(-`e')-1;
	
	if (`todo' == 1) exit ;
	
	*** Hessian *** ;
	tempvar de2 g11 g12 g22;
	qui gen double `de2' = -2*exp(-`e')/(exp(-`e')+1)^2 ;
	
	tempname d11 d12 d22 ;
	mlmatsum `lnf' `d11' = `de2'/`sig'^2, eq(1);
	mlmatsum `lnf' `d12' = `de2'*`e'/`sig'+`de'/`sig' , eq(1,2);
	mlmatsum `lnf' `d22' = `de2'*(`e')^2   + `de'*`e', eq(2) ;
	
	matrix `H' = (`d11', `d12' \
					`d12'', `d22');	
	
end;
