*************************************************************
*********Genome-scale metabolic model of rice panicle********
*************************************************************
*****************Niaz Bahar Chowdhury************************
*************************************************************
$INLINECOM /*  */

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"

;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	mass_balance(i)				steady state mass balance
	lower_bound(j)				lower bounds on reactions
	upper_bound(j)				upper bounds on reactions
	biomass_fix				fixing biomass value
	reform1(j)				reformulation 1
	reform2(j)				reformulation 2
	oxygen_uptake				oxygen uptake
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)					reaction flux
	Z					objective value
	D(j)					dummy variable
;

****************************************************************

***************Defining Model***********************************
objective..			Z =e= sum(j, D(j));

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

biomass_fix..			v('Seed_Biomass[K]') =e= 2.031;

reform1(j)..			v(j) =l= D(j);

reform2(j)..			- v(j) =l= D(j);

oxygen_uptake..			v('Exe2[K]') =e= -36;

Model rice_kernel /all/;
******************************************************************

**********Solving Model*********************
rice_kernel.optfile = 1;

rice_kernel.holdfixed = 1;

solve rice_kernel using lp minimizing Z;
********************************************
****************Output File*****************
FILE RESULTS /FLUX_DISTRIBUTION_KERNEL.txt/;

PUT RESULTS;

PUT "reaction      FLUX"/;

LOOP(j,
	
	PUT j.tl:0:100,"    ", v.l(j):20:5/;
		
);

PUTCLOSE;
**********************************************