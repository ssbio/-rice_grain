*************************************************************
******************FLUX SUM ANALYSIS**************************
*************************************************************
*****************Niaz Bahar Chowdhury************************
*************************************************************
$INLINECOM /*  */

Options

       	optCR = 0
        optCA = 0
;

*******************Defining Sets*****************************
SETS

	i				set of metabolites

$include "metabolites.txt"

	SM(i)				selected metabolites from metabolomics data

$include "selected_metabolites.txt"

	j				set of reactions

$include "reactions.txt"
;

alias(i,i1);

**************************************************************

***********Defining Parameters********************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "v_max.txt"

	v_min(j)				minimum flux of v(j)

$include "v_min.txt"

	c(i)					iteration counter

;
**************************************************************

**********Defining Scalars************************************
SCALARS
	
	v_wild_type				wild type biomass
	
/2.031/

	coeff					fraction of biomass

/0.99/

	M					a large number

/100000/

;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective			objective function
	mass_balance(i)			steady state mass balance
	rearrangement_1(i,j)		rearrangement of absolute objective
	rearrangement_2(i,j)		rearrangement of absolute objective
	rearrangement_3(i,j)		rearrangement of absolute objective
	rearrangement_4(i,j)		rearrangement of absolute objective
	lower_bound(j)			reaction flux lower bound
	upper_bound(j)			reaction flux upper bound
	biomass_enforcement(j)		biomass enforcement
;
**************************************************************

*********Defining Variables***********************************
FREE VARIABLES

	v(j)				reaction flux
	Z				objective value

POSITIVE VARIABLES

	f_plus(i,j)			absolute value rearrangement
	f_minus(i,j)			absolute value rearrangement

BINARY VARIABLES

	i_plus(i,j)			binary varaible for absolute value rearrangement
	i_minus(i,j)			binary varaible for absolute value rearrangement
;				
****************************************************************

***************Defining Model***********************************

objective..			Z =e= 0.5 * sum((i,j)$SM(i), (c(i) * (f_plus(i,j) + f_minus(i,j))));

mass_balance(i)..		sum(j, S(i,j) * v(j)) =e= 0;

rearrangement_1(i,j)$SM(i)..	S(i,j) * v(j) =e= f_plus(i,j) - f_minus(i,j);

rearrangement_2(i,j)$SM(i)..	f_plus(i,j) =l= i_plus(i,j) * M;
	
rearrangement_3(i,j)$SM(i)..	f_minus(i,j) =l= i_minus(i,j) * M;

rearrangement_4(i,j)$SM(i)..	i_plus(i,j) + i_minus(i,j) =e= 1;	

lower_bound(j)..		v_min(j) =l= v(j);

upper_bound(j)..		v(j) =l= v_max(j);

biomass_enforcement(j)..	v('Seed_Biomass[K]') =e= coeff * v_wild_type;

Model flux_sum_analysis /all/;

*flux_sum_analysis.optfile = 1;

*flux_sum_analysis.holdfixed = 1;

******************************************************************

***********************Iteration**********************************

File results /FSA.txt/;

Put results;

Put " Metabolite	Max N-		Min N-"/;

LOOP(i1$SM(i1),

	c(i) = 0;

	c(i1) = 1;

	Put i1.tl:0:100;
		
	Solve flux_sum_analysis using mip maximizing Z;	

	Put Z.l:10:5
	
	Solve flux_sum_analysis using mip minimizing Z

	Put Z.l:10:5/;
);

Putclose;
	
