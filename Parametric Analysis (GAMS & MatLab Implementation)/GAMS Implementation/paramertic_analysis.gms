$OFFLISTING
$EOLCOM //
$INLINECOM { }
$TITLE Parametric Analysis for determining optimal solutions when pertubations are done on the cost vector or on the right-hand-side vector.

* Options:
OPTION LIMROW = 0; OPTION LIMCOL = 0;
OPTION SOLPRINT = OFF; OPTION SYSOUT = OFF;

SETS
J indices of variables / 1 * 2 /,
L indices of pertubation rate / L01 * L29 /;

// Pertubation on the cost vector
PARAMETER
C(J) cost vector / 1 -1, 2 -3 /, // c = (-1, -3)
Cp(J) cost direction of pertubation / 1 2, 2 1 /, // c' = (2, 1)
Cof2(J) second constraint coefficients /1 -1, 2 2 /,
lambdaL(L) pertubation rate scenarios
/ L01 0, L02 0.1, L03 0.5, L04 0.9,
  L05 1, L06 1.1, L07 1.5, L08 1.9,
  L09 2, L10 2.1, L11 2.5, L12 2.9,
  L13 3, L14 3.1, L15 3.5, L16 3.9,
  L17 4, L18 4.1, L19 4.5, L20 4.9,
  L21 5, L22 5.1, L23 5.5, L24 5.9,
  L25 6, L26 6.1, L27 6.5, L28 6.9,
  L29 7 /
Cf(J) final pertubed cost vector / 1 -1, 2 -3 /;

SCALAR
lambda pertubation rate;

VARIABLE
Z objective function;

POSITIVE VARIABLE X(J) optimal solution variable;

EQUATIONS objfuncCost, constr1Cost, constr2Cost;
//objfuncCost.. Z =e= SUM(J, (C(J) + (lambda * Cp(J))) * X(J));
objfuncCost.. Z =e= SUM(J, Cf(J) * X(J));
constr1Cost.. SUM(J, X(J)) =l= 6;
constr2Cost.. SUM(J, Cof2(J) * X(J)) =l= 6;

MODEL ParmAnalysisCost / objfuncCost, constr1Cost, constr2Cost /;

// Output File
FILE OUT / "OUT-PMTR-ANLYS.TXT" /;
PUT OUT;
// File Heading
PUT @1 "Parametric Analysis Results" /;
PUT @1 "==============================" / /;
PUT @1 "Pertubation of the Cost Vector:" /;
PUT @1 "----------------------------------" /;
// Column Headers
PUT @1 "Scen.:",@10"Lambda:",@20"Obj.f:",@30"Type",@35"Sol?",@43;
LOOP(J, PUT" X(",J.TL:1:0,'):    '; );
LOOP(J, PUT" C(",J.TL:1:0,'):   '; );
PUT /;
// Solves model and enter values
LOOP(L, // for different values of lambda
lambda = lambdaL(L);
Cf(J) = C(J) + (lambda * Cp(J));
SOLVE ParmAnalysisCost MINIMIZING Z using LP;
PUT @1 L.TL:8,@10 lambda:8:6, @20 Z.L:9:4, @30 ParmAnalysisCost.MODELSTAT:3:0,
@35 ParmAnalysisCost.SOLVESTAT:3:0, @40;
LOOP(J, PUT X.L(J):9:4, ' ';);
LOOP(J, PUT Cf(J):9:4, ' ';); PUT /;
);

// Additional Codes for Pertubation on the Right-hand-side vector b
PARAMETER
b(J) initial right-hand-side vector / 1 6, 2 6/,
bp(J) right-hand-side direction of pertubation / 1 -1, 2 1 /,
bf(J) final right-hand-side vector / 1 6, 2 6/;

EQUATIONS objfuncRHS, constr1RHS, constr2RHS;
objfuncRHS.. Z =e= SUM(J, C(J) * X(J));
constr1RHS.. SUM(J, X(J)) =l= bf("1");
constr2RHS.. SUM(J, Cof2(J) * X(J)) =l= bf("2");

MODEL ParmAnalysisRHS / objfuncRHS, constr1RHS, constr2RHS /;

// Output File
PUT / @1 "Pertubation of the Right-hand-side Vector:" /;
PUT @1 "----------------------------------" /;
// Column Headers
PUT @1 "Scen.:",@10"Lambda:",@20"Obj.f:",@30"Type",@35"Sol?",@43;
LOOP(J, PUT" X(",J.TL:1:0,'):    '; );
LOOP(J, PUT" b(",J.TL:1:0,'):   '; );
PUT /;
// Solves model and enter values
LOOP(L, // for different values of lambda
lambda = lambdaL(L);
bf(J) = b(J) + (lambda * bp(J));
SOLVE ParmAnalysisRHS MINIMIZING Z using LP;
PUT @1 L.TL:8,@10 lambda:8:6, @20 Z.L:9:4, @30 ParmAnalysisRHS.MODELSTAT:3:0,
@35 ParmAnalysisRHS.SOLVESTAT:3:0, @40;
LOOP(J, PUT X.L(J):9:4, ' ';);
LOOP(J, PUT bf(J):9:4, ' ';); PUT /;
);