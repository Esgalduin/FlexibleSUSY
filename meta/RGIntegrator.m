BeginPackage["RGIntegrator`"];
EndPackage[];

RGIntegrate::usage = "Integrates a given list of renormalization group
 equations from a starting scale to a destination scale.

 The first parameter of RGIntegrate[] is the list of beta functions.
 Each list element is again a list, where the first entry is the name
 of the parameter and the rest are the 1-loop, 2-loop, ... beta
 functions.

 Example:

 betas = { {g, a g, b g^2 l^2},
           {l, c l + d g, e l^2 + f g^2 + k l g} }

 RGIntegrate[betas, Q1, Q2]
";

Begin["RGIntegrator`Private`"];

(* loop factor for expansion *)
{h};

AddScale[par_, Q_] := Rule[par, par[Q]];

MultiplyLoopFactor[{par_, betas___}, h_] :=
    Join[{par}, MapIndexed[#1 h^(First[#2])&, {betas}]];

PerformIntegrals[expr_] := expr /. Integral -> Integrate;

IntegrateSingleRHS[{par_, betas___}, Q1_, Q2_, Qp_, addScales_, sol_] :=
    par[Q1] -> par[Q2] - Integral[
        Total[{betas} /. addScales /. (sol /. Q1 -> Qp)]/Qp,
        {Qp, Q1, Q2}, Assumptions :> Q1 > 0 && Q2 > 0 && Qp > 0 && Q2 > Q1];

IntegrateRHS[betas_List, Q1_, Q2_, sol_] :=
    Module[{Qp, addScales, ints},
           addScales = AddScale[First[#], Qp]& /@ betas;
           ints = IntegrateSingleRHS[#, Q1, Q2, Qp, addScales, sol]& /@ betas;
           PerformIntegrals /@ ints
          ];

(* integrate beta functions of order h^N,
   given a solution of the (N-1)th order beta functions *)
RGIntegrateRecursively[betas_List, Q1_, Q2_, sol_] :=
    IntegrateRHS[betas, Q1, Q2, sol];

(* solution of empty beta functions *)
RGIntegrateRecursively[{}, Q1_, Q2_] := {};

(* solution of trivial beta functions *)
RGIntegrateRecursively[betas : {{par_}, ___}, Q1_, Q2_] :=
    Rule[First[#][Q1], First[#][Q2]]& /@ betas;

RGIntegrateRecursively[betas_List, Q1_, Q2_] :=
    Module[{sol, betasLower},
           betasLower = Drop[#, -1]& /@ betas;
           sol = RGIntegrateRecursively[betasLower, Q1, Q2];
           RGIntegrateRecursively[betas, Q1, Q2, sol]
      ];

RGIntegrate[beta_List, Q1_, Q2_] :=
    Module[{lbeta},
           lbeta = MultiplyLoopFactor[#, h]& /@ beta;
           RGIntegrateRecursively[lbeta, Q1, Q2] /. h -> 1
          ];

End[];
