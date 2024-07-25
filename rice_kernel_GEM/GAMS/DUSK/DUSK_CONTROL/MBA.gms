*************************************************************
*****************MBA for rice seed model*********************
*************************************************************
*****************Niaz Bahar Chowdhury************************
*************************************************************
$INLINECOM /*  */

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"
;

alias(j,j1);
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"

	proxy_v_max(j)				proxy of v_max set to reset v_max in the iteration counter
	
	proxy_v_min(j)				proxy of v_min set to reset v_min in the iteration counter
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance(i)				steady state mass balance
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	Z					objective value
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= v('Seed_Biomass[K]');

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

Model rice_kernel /all/;

****************Solving the model iteratively*******************

proxy_v_max(j) = v_max(j);

proxy_v_min(j) = v_min(j);

FILE RESULTS /rice_kernel_bottlenecks.txt/;

PUT RESULTS;

PUT "Nominal Biomass:"

solve rice_kernel using lp maximizing Z;

PUT Z.l:20:5/;

PUT "reaction      BIOMASS"/;

LOOP(j1,

	v_min(j1) = -1000;
	
	v_max(j1) =  1000;
	
	PUT j1.tl:0:100;
	
	solve rice_kernel using lp maximizing Z;

	PUT Z.l:20:5/;

	v_min(j1) = proxy_v_min(j1);
	
	v_max(j1) = proxy_v_max(j1);	
);

PUTCLOSE;
**********************************************