{smcl}
{* *! version 1.0.0 Dec2011}{...}
{cmd: help switchcopula}{right: ({browse "http://www.stata-journal.com/article.html?article=st0308":SJ13-3: st0308})}
{hline}

{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{hi:switchcopula} {hline 2}}Maximum likelihood estimation of the copula-based endogenous switching regression model{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 20 2}
{cmd:switchcopula}
({it:{help depvar:depvar0}} [{cmd:=}] {it:{help varlist:varlist0}})
[({it:depvar1} [{cmd:=}] {it:varlist1})] 
{ifin} 
{weight}{cmd:,}
{cmdab:sel:ect(}{it:{help depvar:depvar_s}} [{cmd:=}]
     {it:{help varlist:varlist_s}}{cmd:)}
[{cmd:copula0(}{it:{help switchcopula##copula:copula}}{cmd:)} 
{cmd:copula1(}{it:{help switchcopula##copula:copula}}{cmd:)} 
{cmd:margsel(}{it:{help switchcopula##margsel:margin}}{cmd:)} 
{cmd:margin0(}{it:{help switchcopula##margin0:margin}}{cmd:)} 
{cmd:margin1(}{it:{help switchcopula##margin1:margin}}{cmd:)} 
{cmd:df0(}{it:#}{cmd:)}
{cmd:df1(}{it:#}{cmd:)}
{cmd:negative0}
{cmd:negative1}
{cmdab:con:sel}
{cmd:vce(}{it:{help switchcopula##vcetype:vcetype}}{cmd:)}
{it:{help maximize:maximize_options}}]


{title:Description}

{pstd}
{cmd:switchcopula} fits the copula-based endogenous switching
regression model by the maximum likelihood method.


{title:Options}

{phang}
{cmd:select(}{it:depvar_s} [{cmd:=}] {it:varlist_s}{cmd:)}
specifies the selection equation.  {it:depvar_s} should be coded as
0 and 1, with 0 indicating an observation being in regime 0 and 1
indicating an observation being in regime 1. {cmd:select()} is required.

{marker copula}{...}
{phang}
{opt copula0(copula)} specifies a copula function for the
dependence between the errors in the regime 0 equation and selection
equation.  {it:copula} may be one of the following:

{p 8 8 2}
{cmd:product}, {cmd:gaussian}, {cmd:fgm}, {cmd:plackett},
{cmd:amh}, {cmd:clayton}, {cmd:frank}, {cmd:gumbel}, {cmd:joe}

{pmore}
The default is {cmd:copula0(gaussian)}.  Except for {cmd:product}, the result
table displays estimates of an ancillary dependence parameter ({cmd:atheta0})
and a dependence parameter ({cmd:theta0}).

{phang}
{opt copula1(copula)} specifies a copula function for the
dependence between the errors in the regime 1 equation and selection
equation.  See {help switchcopula##copula:above} for the list of available
copulas.  The default is {cmd:copula1(gaussian)}.  Except for {cmd:product},
the result table displays estimates of an ancillary dependence parameter
({cmd:atheta1}) and a dependence parameter ({cmd:theta1}).

{marker margsel}{...}
{phang}
{opt margsel(margin)} specifies the marginal distribution of
the error term in the selection equation.  {it:margin} may be
{cmd:normal} (or {cmd:probit}) or {cmd:logistic} (or {cmd:logit}).
The default is {cmd:margsel(normal)}.

{marker margin0}{...}
{phang}
{opt margin0(margin)} specifies the marginal distribution of the 
error term in regime 0.  {it:margin} may be {cmd:normal}, {cmd:logistic},
or {cmd:t}.  The default is {cmd:margin0(normal)}.

{marker margin1}{...}
{phang}
{opt margin1(margin)} specifies the marginal distribution of the
error term in regime 1.  {it:margin} may be {cmd:normal}, {cmd:logistic},
or {cmd:t}.  The default is {cmd:margin1(normal)}.

{phang}
{opt df0(#)} fixes the degrees of freedom if {cmd:margin0()} is
{cmd:t}.  The specified value must be greater than 0.  When
{cmd:margin0()} is {cmd:t} and {cmd:df0()} is not specified, the degrees
of freedom will be a parameter to estimate.  The result table reports an
ancillary parameter ({cmd:lndf0}, log of degrees of freedom) and an
estimated degree of freedom, {cmd:df0}.  If {cmd:margin0} is not
{cmd:t}, this option will be ignored.

{phang}
{opt df1(#)} fixes the degrees of freedom if {cmd:margin1()} is
{cmd:t}; see {cmd:df0()}.

{phang}
{opt negative0} makes the error term of the regime 0 equation
negative.  This option allows a negative dependence between the regime 0
and selection equations.

{phang}
{opt negative1} makes the error term of the regime 1 equation
negative.  This option allows a negative dependence between the regime 1
and selection equations.

{phang}
{opt consel} allows contributions to the likelihood of the
selection equation by observations in which the selection decision is
observed but in which the outcome variables or some of the covariates in
the outcome equations are not observed.

{marker vcetype}{...}
{phang}
{opt vce(vcetype)} specifies the type of standard errors
reported; see {it:{help vce_option}}.

{phang}
{it:maximization_options} control the maximization process; see
{helpb maximize}.


{title:Syntax for predict}

{p 8 17 2}
{cmd:predict}
{dtype} 
{newvar} 
{ifin} 
[{cmd:,} {it:options}]


{title:Options for predict}

{phang}
{cmd:psel} computes the probability of being in regime 1. This is
a default.

{phang}
{cmd:xbsel} computes the linear prediction for the selection
equation.

{phang}
{cmd:xb0} computes the linear prediction of the dependent
variable in regime 0.

{phang}
{cmd:xb1} computes the linear prediction of the dependent
variable in regime 1.

{phang}
{cmd:cll} computes the contribution to the log likelihood by each
observation.

{phang}
{cmd:y0_c0} computes the expected value of the dependent variable
in regime 0 conditional on being in regime 0: E(y0|S=0).  If
{cmd:copula0()} is {cmd:gaussian} and {cmd:margin0()} is {cmd:normal},
it is computed analytically; otherwise, it is computed numerically.

{phang}
{cmd:y0_c1} computes the expected value of the dependent variable
in regime 0 conditional on being in regime 1: E(y0|S=1).  If
{cmd:copula0()} is {cmd:gaussian} and {cmd:margin0()} is {cmd:normal},
it is computed analytically; otherwise, it is computed numerically.

{phang}
{cmd:y1_c0} computes the expected value of the dependent variable
in regime 1 conditional on being in regime 0: E(y1|S=0).  If
{cmd:copula1()} is {cmd:gaussian} and {cmd:margin1()} is {cmd:normal},
it is computed analytically; otherwise, it is computed numerically.

{phang}
{cmd:y1_c1} computes the expected value of the dependent variable
in regime 1 conditional on being in regime 1: E(y1|S=1).  If
{cmd:copula1()} is {cmd:gaussian} and {cmd:margin1()} is {cmd:normal},
it is computed analytically; otherwise, it is computed numerically.


{title:Stored results}

{pstd}
{cmd:switchcopula} stores the following in {cmd:e()}:

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
{synopt:{cmd:e(df0)}}fixed value of {cmd:df0()}; only when option {cmd:df0()} is specified{p_end}
{synopt:{cmd:e(df1)}}fixed value of {cmd:df1()}; only when option {cmd:df1()} is specified{p_end}
{synopt:{cmd:e(negative0)}}{cmd:1} if option {cmd:negative0} is specified, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(negative1)}}{cmd:1} if option {cmd:negative1} is specified, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:switchcopula}{p_end}
{synopt:{cmd:e(depvar)}}names of dependent variable{p_end}
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
{synopt:{cmd:e(copula0)}}specified {cmd:copula0()}{p_end}
{synopt:{cmd:e(copula1)}}specified {cmd:copula1()}{p_end}
{synopt:{cmd:e(margsel)}}specified {cmd:margsel()}{p_end}
{synopt:{cmd:e(margin0)}}specified {cmd:margin0()}{p_end}
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
Help:  {helpb movestay}, {helpb heckmancopula} (if installed){p_end}
