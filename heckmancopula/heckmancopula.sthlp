{smcl}
{* *! version 1.0.0 January 2012}{...}
{cmd:help heckmancopula}{right: ({browse "http://www.stata-journal.com/article.html?article=st0308":SJ13-3: st0308})}
{hline}

{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{hi:heckmancopula} {hline 2}}Maximum likelihood estimation of the copula-based Heckman sample-selection model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:heckmancopula}
{depvar} [{cmd:=}] {varlist}
{ifin} 
{weight}{cmd:,}
{cmdab:sel:ect(}[{it:{help depvar:depvar_s}} {cmd:=}]
   {it:{help varlist:varlist_s}}{cmd:)} 
[{cmd:copula(}{it:{help heckmancopula##copula:copula}}{cmd:)} 
{cmd:margsel(}{it:{help heckmancopula##margsel:margin}}{cmd:)} 
{cmd:margin1(}{it:{help heckmancopula##margin1:margin}}{cmd:)} 
{cmd:df(}{it:#}{cmd:)}
{cmd:negative}
{cmdab:nocons:tant}
{cmd:vce(}{it:{help heckmancopula##vcetype:vcetype}}{cmd:)}
{it:{help maximize:maximize_options}}]


{title:Description}

{pstd}{cmd:heckmancopula} estimates a copula-based Heckman-type
sample-selection model by the maximum likelihood method. 


{title:Options}

{phang}
{cmd:select(}[{it:{help depvar:depvar_s}} {cmd:=}]
{it:{help varlist:varlist_s}}{cmd:)}
specifies a selection equation.  If {it:depvar_s} is specified, it
should be coded as 0 and 1, with 0 indicating an outcome not observed
for an observation and 1 indicating an outcome observed for an
observation.  {cmd:select()} is required.

{marker copula}{...}
{phang}
{opt copula(copula)} specifies a copula function governing the
dependence between the errors in the outcome equation and selection
equation.  {it:copula} may be one of the following: 

{p 8 8 2}{cmd:product}, {cmd:gaussian}, {cmd:fgm}, {cmd:plackett},
{cmd:amh}, {cmd:clayton}, {cmd:frank}, {cmd:gumbel}, {cmd:joe}

{pmore}
The default is {cmd:copula(gaussian)}.  Except for {cmd:product}, the result
table displays estimates of an ancillary dependence parameter ({bf:atheta})
and a dependence parameter ({bf:theta}).

{marker margsel}{...}
{phang}
{opt margsel(margin)} specifies the marginal distribution of the
error term in the selection equation.  {it:margin} may be
{cmd:normal} (or {cmd:probit})  or {cmd:logistic} (or {cmd:logit}).
The default is {cmd:margsel(normal)}.

{marker margin1}{...}
{phang}
{opt margin1(margin)} specifies the marginal distribution of the error term in
the outcome equation.  {it:margin} may be {cmd:normal}, {cmd:logistic}, or
{cmd:t}.  The default is {cmd:margin1(normal)}.

{phang}
{opt df(#)} fixes the degrees of freedom if {cmd:margin1()} is
{cmd:t}.  The specified value must be greater than 0.  When
{cmd:margin1()} is {cmd:t} and {cmd:df()} is not specified, the degrees
of freedom will be a parameter to estimate.  The result table reports an
ancillary parameter ({cmd:lndf}, log of degrees of freedom) and an
estimated degree of freedom, {cmd:df()}.  If {cmd:margin1()} is not
{cmd:t}, this option will be ignored.

{phang}
{opt negative} makes the error term of the outcome equation
negative.  This option allows a negative dependence between the
selection and outcome equations.

{phang}
{opt nonconstant} suppresses a constant term of the outcome
equation.

{marker vcetype}{...}
{phang}
{opt vce(vcetype)} specifies the type of standard errors reported; see
{it:{help vce_option}}.

{phang}
{it:maximization_options} control the maximization process; see
{helpb maximize}.


{title:Syntax for predict}

{p 8 17 2}
{cmd:predict}
{dtype} 
{newvar} {ifin}
[{cmd:,} {it:options}]


{title:Options for predict}

{phang}
{cmd:psel} computes the probability of the outcome being
observed.  This is a default.

{phang}
{cmd:xbsel} computes the linear prediction for the selection
equation.

{phang}
{cmd:xb} computes the linear prediction of the dependent variable
in the outcome equation.

{phang}
{cmd:cll} computes the contribution to the log-likelihood
function of each observation.  This will be useful to conduct Vuong's
test.

{phang}
{cmd:y_c0} computes the prediction of the dependent variable in
the outcome equation, conditional on not being observed: E(y|S=0).  If
{cmd:copula()} is {cmd:gaussian} and {cmd:margin1()} is {cmd:normal}, it
is computed analytically; otherwise, it is computed numerically.

{phang}
{cmd:y_c1} computes the prediction of the dependent variable in
the outcome equation, conditional on being observed: E(y|S=1).  If
{cmd:copula()} is {cmd:gaussian} and {cmd:margin1()} is {cmd:normal}, it
is computed analytically; otherwise, it is computed numerically.


{title:Stored results}

{pstd}
{cmd:heckmancopula} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_aux)}}number of auxiliary parameters{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(p)}}significance {p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(ll0)}}log likelihood, independent model{p_end}
{synopt:{cmd:e(AIC)}}Akaike information criterion{p_end}
{synopt:{cmd:e(BIC)}}Bayesian information criterion{p_end}
{synopt:{cmd:e(df)}}fixed value of {cmd:df}; only when option {cmd:df} is specified{p_end}
{synopt:{cmd:e(negative)}}{cmd:1} if option {cmd:negative} is specified, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:heckmancopula}{p_end}
{synopt:{cmd:e(depvar)}}names of dependent variables{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared
	test{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(crittype)}}optimization criterion{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(copula)}}specified {cmd:copula()}{p_end}
{synopt:{cmd:e(margsel)}}specified {cmd:margsel()}{p_end}
{synopt:{cmd:e(margin1)}}specified {cmd:margin1()}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:Author}

{pstd}Takuya Hasebe{p_end}
{pstd}Graduate Center, City University of New York{p_end}
{pstd}New York, NY{p_end}
{pstd}thasebe@gc.cuny.edu{p_end}


{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 13, number 3: {browse "http://www.stata-journal.com/article.html?article=st0308":st0308}

{p 7 14 2}
Help:  {manhelp heckman R}, {helpb switchcopula}
{p_end}
