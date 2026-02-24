mata:
version 10.1
mata set matastrict on 
mata set matalnum on 
void _conditional(string scalar probability,
				string scalar newvarlist,
				string scalar touse,
				string scalar copula,
				string scalar margin,
				string scalar truncation,
				string scalar atheta_t,
				string scalar df_t)
{
	real matrix P, mill, u1, us, vs, v1, C1s, e1, u, r, root, C, phi1_u1, phi1_us, phi1_C, phi2_C
	real scalar N, S, atheta, theta, df
	
	st_view(P=.,.,tokens(probability),touse)
	st_view(mill=.,.,tokens(newvarlist),touse)
	
	atheta = st_numscalar(atheta_t)
	df = st_numscalar(df_t)
	
	N = rows(P)
	S = 500
	
	u = halton(N*S,2,100)
	//u = uniform(N*S,2)

	u1 = J(N,4*S,.); us = J(N,4*S,.)
	u1[.,1::S] = colshape(u[.,2],S)
	u1[.,(S+1)::(2*S)] = 1:-u1[.,1::S]
	u1[.,(2*S+1)::(3*S)]=u1[.,1::S]
	u1[.,(3*S+1)::(4*S)]=1:-u1[.,1::S]
	us[.,1::S] = colshape(u[.,1],S)
	us[.,(S+1)::(2*S)] = 1:-us[.,1::S]
	us[.,(2*S+1)::(3*S)]=1:-us[.,1::S]
	us[.,(3*S+1)::(4*S)]=us[.,1::S]
	
	if (truncation == "below") {
		us = us:*(1:-P):+P
	}
	else {
		us = us:*P
	}

	if (copula == "gaussian") {
		theta = tanh(atheta)
		vs = invnormal(us)
		v1 = invnormal(u1)
		C1s = exp(-0.5*((theta*v1):^2+(theta*vs):^2-2*theta*v1:*vs)/(1-(theta)^2))/sqrt(1-(theta)^2);
	}
	else if (copula == "fgm" ) {
		theta = tanh(atheta)
		C1s = 1:+theta*(1:-2*u1):*(1:-2*us);
	}
	else if (copula == "plackett") {
		theta = exp(atheta)
		r = 1:+(theta-1)*(u1+us)
		root = r:^2 - 4:*u1:*us*theta*(theta-1)
		C1s = -(0.5*(-theta-1)):/sqrt(root)+(0.25*((theta-1)*(u1+us)-2*us*theta:+1):*(2*(theta-1)*(r)-4*(theta-1)*theta*u1)):/((root):^(3/2))
	}
	else {
		if (copula == "amh") {
			theta = tanh(atheta)
			C = u1:*us:/(1:-theta*(1:-us):*(1:-u1));	
			phi1_u1 = theta:/(1:-theta*(1:-u1))-1:/u1;
			phi1_us = theta:/(1:-theta*(1:-us))-1:/us;
			phi1_C = theta:/(1:-theta*(1:-C))-1:/C;
			phi2_C = -(theta:/(1:-theta*(1:-C))):^2+1:/(C:^2);
		}
		else if (copula == "clayton") {
			theta = exp(atheta)
			C = (us:^(-theta)+u1:^(-theta):-1):^(-1/theta)
			phi1_u1 = (-1):*(u1):^(-theta-1)
			phi1_us = (-1):*(us):^(-theta-1)
			phi1_C = (-1):*(C):^(-theta-1)
			phi2_C = (theta+1)*C:^(-theta-2)
		}
		else if (copula=="frank") {
			theta = atheta
			C = -ln(1:+(exp(-theta*us):-1):*(exp(-theta*u1):-1):/(exp(-theta)-1))/theta
			phi1_u1 = theta:/(1:-exp(theta*u1))
			phi1_us = theta:/(1:-exp(theta*us))
			phi1_C = theta:/(1:-exp(theta*C))
			phi2_C = (theta:/(1:-exp(theta*C))):^2:*exp(theta*C)
		}
		else if (copula=="gumbel") {
			theta = 1 + exp(atheta)
			C = exp(-((-ln(u1)):^(theta)+(-ln(us)):^(theta)):^(1/theta))
			phi1_u1 = -theta*(-ln(u1)):^(theta-1):/u1
			phi1_us = -theta*(-ln(us)):^(theta-1):/us
			phi1_C = -theta*(-ln(C)):^(theta-1):/C
			phi2_C = theta*(-ln(C)):^(theta-1):*((theta-1):/(-ln(C)):+ 1):/(C:^2)
		}
		else if (copula=="joe") {
			theta = 1 + exp(atheta)
			C = 1:-((1:-us):^theta+(1:-u1):^(theta) -((1:-us):*(1:-u1)):^theta):^(1/theta)
			phi1_u1 = -theta*(1:-u1):^(theta-1):/(1:-(1:-u1):^theta)
			phi1_us = -theta*(1:-us):^(theta-1):/(1:-(1:-us):^theta)
			phi1_C = -theta*(1:-C):^(theta-1):/(1:-(1:-C):^theta)
			phi2_C = (phi1_C):^2 - (theta-1):*phi1_C:/(1:-C)
		}
		
		C1s = -phi2_C:*phi1_u1:*phi1_us:/(phi1_C:^3)
	}
	if (margin=="normal") {
		e1 = invnormal(u1)
	}
	else if (margin=="logistic") {
		e1 = ln(u1:/(1:-u1))
	}
	else if (margin == "t") {
		e1 = (invttail(df,1:-u1)) 
	}	
	e1 = C1s:*e1
	mill[.,.] = rowsum(e1)/cols(e1)
}
mata set matalnum off
end
mata: mata mosave _conditional(), dir(PERSONAL) replace
