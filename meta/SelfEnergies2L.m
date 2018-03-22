BeginPackage["SelfEnergies2L`", {"SARAH`","TreeMasses`","CConversion`"}];

Unprotect["SelfEnergies2L`*"];
ClearAll["SelfEnergies2L`*"];
ClearAll["SelfEnergies2L`Private`*"];


ConvertSarah2LDiagramList::usage = "Converts SARAH's list with 2-loop
 self-energy and tadpole diagrams."

(* GetTadpoleField::usage = "Returns the field to which the tadpole diagram expression belongs to." *)

Make1L2LShifts::usage = "takes in list of SARAH tadpoles and eigenstates, and gives all 1L2L Tadpole shifts"

CreateEnterGauglessLimitFunction::usage = "writes the entergaugelesslimit function in mass_eigenstates.cpp.in"

Begin["`Private`"];

AppendIndex[p_, idx_, 1] := p;
AppendIndex[p_, idx_, range_] := AppendIndex[p, idx];
AppendIndex[SARAH`bar[p_], idx_] := SARAH`bar[AppendIndex[p, idx]];
AppendIndex[Susyno`LieGroups`conj[p_], idx_] := Susyno`LieGroups`conj[AppendIndex[p, idx]];
AppendIndex[p_, idx_] := p[{idx}];

DistChir[{}, {}] := 1;
DistChir[{0, rest1___}, {c_, rest2___}] := c     DistChir[{rest1}, {rest2}];
DistChir[{L, rest1___}, {c_, rest2___}] := c[PL] DistChir[{rest1}, {rest2}];
DistChir[{R, rest1___}, {c_, rest2___}] := c[PR] DistChir[{rest1}, {rest2}];

DistributeChiralities[chiralities_List, couplings_List] :=
    If[chiralities === (chiralities /. {L -> R, R -> L}),
       DistChir[chiralities, couplings],
       DistChir[chiralities, couplings] + (DistChir[chiralities /. {L -> R, R -> L}, couplings])
      ];

MultLF[{factor_, func_, chiralities_List}, particles_List, expr_] :=
    factor func[Sequence @@ (
        If[#3 === 1, SARAH`Mass2[#1],
           SARAH`Mass2[#1, #2]]& @@@ particles
    )] DistributeChiralities[chiralities, expr];

MultiplyLoopFunction[loopfuncs_List, particles_List, expr_] :=
    Total[MultLF[#, particles, expr]& /@ loopfuncs];

(* check for ambiguous contraction of indices *)
IsAmbiguousContraction[indices_List] :=
    Or @@ ((#[[2]] > 2)& /@
           Tally[First /@ (DeleteCases[indices, {_, SARAH`gE1 | SARAH`gE2}] /.
                           SARAH`bar -> Identity /. Susyno`LieGroups`conj -> Identity)]);

(* one or more fields with indices appear more once *)
IsAmbiguousIndex[indices_List] :=
    Or @@ ((#[[2]] > 1)& /@ Tally[First /@ indices]);

DistributeIndices[{}, coupling_] := coupling;
DistributeIndices[indices_List, c : C[particles__]] :=
    If[IsAmbiguousIndex[indices],
       (* replace indices in order from indices list *)
       Assert[First /@ indices === {particles}];
       C[Sequence @@
         MapThread[(#1 /. #2)&, {{particles}, (Rule[#1, AppendIndex[#1, #2]]& @@@ indices)}]]
       ,
       (* use one common replacement rule for all fields *)
       c /. (Rule[#1, AppendIndex[#1, #2]]& @@@ indices)
    ];

(* multiplies all given couplings *)
MultiplyCouplings[lst_List] := DistributeIndices @@@ lst;

SumOverIndices[{}, expr_] := expr;
SumOverIndices[{{_, _, 1}, rest___}, expr_] :=
    SumOverIndices[{rest}, expr];
SumOverIndices[{{_, idx_, range_}, rest___}, expr_] :=
    SumOverIndices[{rest}, SARAH`sum[idx, 1, range, expr]];

(* combines list of particles/loop functions/couplings to a sum *)
CreateTadpoleDiag[{particles_List, loopFuncs_List, couplings_List}] :=
    SumOverIndices[particles,
                   MultiplyLoopFunction[loopFuncs, particles,
                                        MultiplyCouplings[couplings]]];

(* sums all diagrams that for a given type *)
SumTadpoleType[{type_, diags_List}] := Total[CreateTadpoleDiag /@ diags];

(*returns the field the tadpole diagram corresponds to*)
(* GetTadpoleField[expr_] := Module[{fields, prep1, prep2},
   If[Head[expr] === Plus, prep1 = (List @@ expr)[[1]],prep1 = expr];
	 prep2 = prep1 //. {Conj[a_] -> a, conj[a_] -> a, bar[b_] -> b, c_[{b__}] /; c =!= Cp -> c, Cp[x___][a_] -> Cp[x]};
   fields = Tally[Cases[List[prep2], Cp[a___] -> a, Infinity]];
   DeleteCases[fields, {a_, _?EvenQ}][[1, 1]]
   ]; *)

(* if there is more than one generation of HiggsBosons or Pseudoscalars, respectively attach a "U" to their name,*)
(* plus exchange the index of form gE1 with gO1 *)
UnrotateRules[]:=Module[{rules={}},
          If[TreeMasses`GetDimension[SARAH`HiggsBoson]>1,
            rules=Append[rules,(SARAH`HiggsBoson[{i_}] /; StringMatchQ[ToString[i], RegularExpression["gE[0-9]*"]]) :> ToExpression["U" <> ToString[SARAH`HiggsBoson[{Symbol[StringReplace[ToString[i], "gE" -> "gO"]]}]]]];
          ];
          If[TreeMasses`GetDimension[SARAH`PseudoScalar]>1,
            rules=Append[rules,(SARAH`PseudoScalar[{i_}] /; StringMatchQ[ToString[i], RegularExpression["gE[0-9]*"]]) :> ToExpression["U" <> ToString[SARAH`PseudoScalar[{Symbol[StringReplace[ToString[i], "gE" -> "gO"]]}]]]];
          ];
          rules
        ];

simpleUnrotateRules[]:={Symbol["U"<>ToString[SARAH`HiggsBoson]][x_]->SARAH`HiggsBoson[x],Symbol["U"<>ToString[SARAH`PseudoScalar]][x_]->SARAH`PseudoScalar[x]};


AddSEMomDep[] := {Symbol["WfSSSS"][masses___] -> Symbol["WfSSSS"][p^2,masses],  Symbol["XfSSS"][masses___] -> Symbol["XfSSS"][p^2,masses],  Symbol["YfSSSS"][masses___] -> Symbol["YfSSSS"][p^2,masses],
                  Symbol["SfSSS"][masses___] -> Symbol["SfSSS"][p^2,masses],  Symbol["UfSSSS"][masses___] -> Symbol["UfSSSS"][p^2,masses],  Symbol["VfSSSSS"][masses___] -> Symbol["VfSSSSS"][p^2,masses],
                  Symbol["WfSSSV"][masses___] -> Symbol["WfSSSV"][p^2,masses],  Symbol["MfSSSSV"][masses___] -> Symbol["MfSSSSV"][p^2,masses], Symbol["WfSSFF"][masses___] -> Symbol["WfSSFF"][p^2,masses],
                  Symbol["WfSSFbFb"][masses___] -> Symbol["WfSSFbFb"][p^2,masses],  Symbol["MfFbFbFbFbS"][masses___] -> Symbol["MfFbFbFbFbS"][p^2,masses],  Symbol["MfFFbFbFS"][masses___] -> Symbol["MfFFbFbFS"][p^2,masses],
                  Symbol["MfFFbFFbS"][masses___] -> Symbol["MfFFbFFbS"][p^2,masses],  Symbol["MfFFFbFbS"][masses___] -> Symbol["MfFFFbFbS"][p^2,masses],  Symbol["MfFFFFS"][masses___] -> Symbol["MfFFFFS"][p^2,masses],
                  Symbol["MfSFbSFbFb"][masses___] -> Symbol["MfSFbSFbFb"][p^2,masses],  Symbol["MfSFSFbF"][masses___] -> Symbol["MfSFSFbF"][p^2,masses], Symbol["MfSFSFFb"][masses___] -> Symbol["MfSFSFFb"][p^2,masses],
                  Symbol["VfSSSFbFb"][masses___] -> Symbol["VfSSSFbFb"][p^2,masses],  Symbol["VfSSSFF"][masses___] -> Symbol["VfSSSFF"][p^2,masses],  Symbol["VfFbFbFbFbS"][masses___] -> Symbol["VfFbFbFbFbS"][p^2,masses],
                  Symbol["VfFbFFbFS"][masses___] -> Symbol["VfFbFFbFS"][p^2,masses],  Symbol["VfFbFFFbS"][masses___] -> Symbol["VfFbFFFbS"][p^2,masses],  Symbol["VfFFbFbFS"][masses___] -> Symbol["VfFFbFbFS"][p^2,masses],
                  Symbol["VfFFFbFbS"][masses___] -> Symbol["VfFFFbFbS"][p^2,masses],  Symbol["VfFFFFS"][masses___] -> Symbol["VfFFFFS"][p^2,masses],  Symbol["GfFFV"][masses___] -> Symbol["GfFFV"][p^2,masses],
                  Symbol["GfFbFbV"][masses___] -> Symbol["GfFbFbV"][p^2,masses],Symbol["ZfSSSS"][masses___] -> Symbol["ZfSSSS"][p^2,masses],Symbol["MfSSSSS"][masses___] -> Symbol["MfSSSSS"][p^2,masses]};

(* result as from CalculatePi2S in SPheno except for a global (-1),
   to be consistent with FS 1L expressions. See comment at function Make1L2LShifts
   at the bottom. *)
ConvertSarah2LDiagramList[tad_List, head_:Total] :=
    head[SumTadpoleType /@ tad]*(-1) //. {
        (m : (SARAH`Mass | SARAH`Mass2))[(SARAH`bar | Susyno`LieGroups`conj)[p_], idx___] :> m[p, idx],
        (m : (SARAH`Mass | SARAH`Mass2))[p_, idx__] :> m[p[{idx}]],
        C[p__] :> Cp[p]} //. {Symbol["i1"]->SARAH`gI1, Symbol["i2"]->SARAH`gI2, Symbol["i3"]->SARAH`gI3,
        Symbol["i4"]->SARAH`gI4,Symbol["i5"]->SARAH`gI5,Symbol["i6"]->SARAH`gI6} //. UnrotateRules[] /. AddSEMomDep[];




(* given a 1-loop nPoint list of diagrams, determines the field that these diagrams correspond to *)
GetnPointField[tempdiags_List]:=Module[{testdiag=tempdiags[[1]],coupfields,loopfields,tempnPointField},
  loopfields={testdiag[[1]],testdiag[[2]]}/. {Conj[a_] -> a, conj[a_] -> a, bar[b_] -> b}/.{x_[b_]->x};
  coupfields=testdiag[[3]]/.{Cp[tempfields__]->List[tempfields]}/. {Conj[a_] -> a, conj[a_] -> a, bar[b_] -> b}/.{x_[b_]->x};
  For[i = 1, i <= Length[loopfields], i++, coupfields = Delete[coupfields, FirstPosition[coupfields, loopfields[[i]]]];];
  coupfields=DeleteDuplicates[coupfields];
  If[Length[coupfields] != 1,Print["Problem while determining to which field the diagram "];
                           Print[testdiag];
                           Print["corresponds to."]
    ];
  coupfields[[1]]
  ];

GetUsedParameters[expr_] := Select[Transpose[parameters][[1]], ! FreeQ[expr, #] &];

GetShiftedExpressions[SarahList_List,ewsbEqparameters_,eigenstates_]:=Select[#,
   With[{unused = #},ContainsAny[Join[GetUsedParameters[TreeMass[#[[1]], eigenstates]],
      GetUsedParameters[TreeMass[#[[2]], eigenstates]]],ewsbEqparameters]] &] & /@ SarahList;
    (* the 'unused = #'' is set, so that we can use a pure function inside another one *)


VertexZeroQ[fields_List,subs_List,mod_:-1]:= If[mod===-1,
      SameQ[0, Simplify[(Drop[Vertex[fields//.simpleUnrotateRules[]],1][[1,1]])//.subs]],
      SameQ[0, Simplify[(Select[Drop[Vertex[fields//.simpleUnrotateRules[]],1],#[[2]]===mod&][[1,1]])//.subs]]];

(*TadpoleReplacements1L2L:={Cp[a___, Symbol["Uhh"][{SARAH`gI3}] , b___]*Symbol["tadpole"][i_Integer] :> Cp[a, Symbol["Uhh"][{i}], b]*Symbol["tadpole"][i],
Cp[a___, Symbol["Uhh"][{SARAH`gI3}] , b___][c___]*Symbol["tadpole"][i_Integer] :> Cp[a, Symbol["Uhh"][{i}], b][c]*Symbol["tadpole"][i]};*)


(* we use the convention, that the mass shifts have the indices gI1,gI2, and the extra coupling(s) have gI4,gI5, if there are multiple generations *)
(* The contraction of tadpole[i] with the correct tadpole expression happens via a Kronckerdelta(i,gI3), where*)
(* gI3 is always going to be the index of the Higgs tadpole we are attaching *)


(* Note: the A0=-x(log(x)-1) function used by FlexibleSUSY is the same as in SPheno (though they both differ from the paper, where A=x(log(x)-1))*)
Calc1L2LTadShiftExpr[diags_,massshifts_List] := Module[{result={}},
  For[n=1,n <= Length[diags],n++,
    Switch[diags[[n,4]],
      Symbol["SSS"],
        result=Append[result,Calc1L2LTadShiftSSS[diags[[n]],massshifts[[n]]]];,
      Symbol["FFS"],
        result=Append[result,Calc1L2LTadShiftFFS[diags[[n]],massshifts[[n]]]];,
      _,Null;];
  ];
  result
];

(* Tadpole shift functions for the different topologies *)

Calc1L2LTadShiftSSS[diag_List,massshift_]:=Module[{tempexpr,loopfield,loopfunction,nField,couplings,prefactor},
    loopfield=diag[[2]]/.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nField = TreeMasses`GetDimension[loopfield];

    couplings = diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};
    prefactor=4*0.5*diag[[5]]*diag[[6]]*couplings*massshift//.{Symbol["generation"]->SARAH`gI4}; (* SARAH diag basefactor: 4 *)
    loopfunction=-Symbol["BB"][Mass2[loopfield[{SARAH`gI4}]],Mass2[loopfield[{SARAH`gI4}]]]; (* -BB is P_{SS} *)

    If[nField==1, loopfunction=loopfunction//.{x_[{SARAH`gI4}]->x}];
    tempexpr=prefactor*loopfunction;
    If[nField > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nField,tempexpr];];
    tempexpr
];

Calc1L2LTadShiftFFS[diag_List,massshift_]:=Module[{tempexpr,loopfield,loopfunction,nField,couplings,prefactor},
  If[FlexibleSUSY`Exclude1L2LFermionShifts === False,
    loopfield = diag[[2]]/.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nField = TreeMasses`GetDimension[loopfield];

    couplings = diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};
    prefactor = -4*diag[[5]]*diag[[6]]*couplings*{Sqrt[massshift],massshift*SARAH`Mass[loopfield[{SARAH`gI4}]]}//.{Symbol["generation"]->SARAH`gI4}; (* SARAH diag basefactor: 4*)
    loopfunction = {-Symbol[A0][Mass2[loopfield[{SARAH`gI4}]]],-Symbol["BB"][Mass2[loopfield[{SARAH`gI4}]],Mass2[loopfield[{SARAH`gI4}]]]}; (* -BB is P_{SS}, the A0 implemented has opposite sign from the paper*)

    If[nField == 1, loopfunction = loopfunction //. {x_[{SARAH`gI4}]->x}];
    tempexpr = prefactor * loopfunction;
    If[nField > 1,tempexpr = SARAH`sum[SARAH`gI4,1,nField,tempexpr];];
    tempexpr
  ,Nothing]
];


(* calculates the 1L2L shifts to the given selfenergies *)
Calc1L2LSEShiftExpr[diags_,massshifts_List]:=Module[{result={}},
  For[n=1,n<=Length[diags],n++,
    Switch[diags[[n,4]],
      Symbol["SSS"],
        result=Append[result,Calc1L2LSEShiftSSS[diags[[n]],massshifts[[n]]]];,
      Symbol["SSSS"],
        result=Append[result,Calc1L2LSEShiftSSSS[diags[[n]],massshifts[[n]]]];,
      Symbol["FFS"],
        result=Append[result,Calc1L2LSEShiftFFS[diags[[n]],massshifts[[n]]]];,
      _,Null;];
  ];
  result
];

(* SelfEnergy shift functions for the different topologies *)

Calc1L2LSEShiftSSS[diag_List,massshifts_]:=Module[{tempexpr,loopfields,loopfunction,nFields,couplings,prefactors,loopfunctions},
  loopfields={diag[[1]],diag[[2]]}//.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
  nFields = TreeMasses`GetDimension[#]& /@ loopfields;

  couplings= diag[[3]] * (diag[[3]]/.{SARAH`gO1->SARAH`gO2}/. {Cp[tempfields__]:>Cp[Sequence @@ (AntiField /@ List[tempfields])]})/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI5};
  prefactors= -2*diag[[5]]*diag[[6]]*couplings*{massshifts[[1]]//.{Symbol["generation"]->SARAH`gI4},massshifts[[2]]//.{Symbol["generation"]->SARAH`gI5}};
  loopfunctions={Symbol["CCtilde"][Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[1]][{SARAH`gI4}]]],
                 Symbol["CCtilde"][Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[2]][{SARAH`gI5}]]]}; (* C in Braathen paper is -CCtilde, sort of*)

  If[nFields[[1]]==1, loopfunctions=loopfunctions//.{x_[{SARAH`gI4}]->x}];
  If[nFields[[2]]==1, loopfunctions=loopfunctions//.{x_[{SARAH`gI5}]->x}];
  tempexpr=Plus @@ (prefactors*loopfunctions);
  If[nFields[[1]] > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
  If[nFields[[2]] > 1,tempexpr=SARAH`sum[SARAH`gI5,1,nFields[[2]],tempexpr];];
  tempexpr
];

Calc1L2LSEShiftSSSS[diag_List,massshifts_]:=Module[{tempexpr,loopfields,loopfunction,nFields,couplings,prefactors,loopfunctions},
  If[!(MatchQ[GetnPointField[{diag}],SARAH`PseudoScalar] && FlexibleSUSY`Exclude1L2LAhShiftSSSS===True),
    loopfields={diag[[2]]}//.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nFields = TreeMasses`GetDimension[#]& /@ loopfields;

    couplings= diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};
    If[MatchQ[couplings, x_[{SARAH`gO1}]],couplings=ReplacePart[couplings, FirstPosition[couplings, gO1] -> SARAH`gO2];];
    prefactors=2*0.5*diag[[5]]*diag[[6]]*couplings*{massshifts[[1]]//.{Symbol["generation"]->SARAH`gI4}};
    loopfunctions={-Symbol["BB"][Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[1]][{SARAH`gI4}]]]};

    If[nFields[[1]]==1, loopfunctions=loopfunctions//.{x_[{SARAH`gI4}]->x}];
    tempexpr=Plus @@ (prefactors*loopfunctions);
    If[nFields[[1]] > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
    tempexpr
  ,Nothing]
];

Calc1L2LSEShiftFFS[diag_List,massshifts_]:=Module[{tempexpr,loopfields,loopfunction,nFields,couplings,prefactors,loopfunctions},
   If[FlexibleSUSY`Exclude1L2LFermionShifts === False,
     loopfields = {diag[[1]],diag[[2]]}//.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
     nFields = TreeMasses`GetDimension[#]& /@ loopfields;

     couplings = diag[[3]] * (diag[[3]]/.{SARAH`gO1->SARAH`gO2}/. {Cp[tempfields__]:>Cp[Sequence @@ (AntiField /@ List[tempfields])]})/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI5};
     prefactors = diag[[5]]*diag[[6]]*couplings*{-2*Sqrt[massshifts[[1]]]/.{Symbol["generation"]->SARAH`gI4}*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]], (* the first 4 parts correspond to the -2*m1*m2*B(p2,m12,m22) part of the SE *)
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*Sqrt[massshifts[[2]]]/.{Symbol["generation"]->SARAH`gI5},
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]*massshift[[1]]/.{Symbol["generation"]->SARAH`gI4},
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]*massshift[[2]]/.{Symbol["generation"]->SARAH`gI5},
                                                +massshift[[1]]/.{Symbol["generation"]->SARAH`gI4},  (* these last parts correspond to the +G0(p2,m12,m22) part of the SE *)
                                                +massshift[[2]]/.{Symbol["generation"]->SARAH`gI5},
                                                -((massshift[[1]]/.{Symbol["generation"]->SARAH`gI4})+(massshift[[2]]/.{Symbol["generation"]->SARAH`gI5})),
                                                -(Mass2[loopfields[[1]][{SARAH`gI4}]] + Mass2[loopfields[[2]][{SARAH`gI5}]]-Symbol["p2"])*massshift[[1]]/.{Symbol["generation"]->SARAH`gI4},
                                                -(Mass2[loopfields[[1]][{SARAH`gI4}]] + Mass2[loopfields[[2]][{SARAH`gI5}]]-Symbol["p2"])*massshift[[2]]/.{Symbol["generation"]->SARAH`gI5}};

     loopfunctions = {+Symbol["BBs"][Symbol["p2"], Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["BBs"][Symbol["p2"], Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["CCtilde"][Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      +Symbol["CCtilde"][Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      -Symbol["BB"][Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      -Symbol["BB"][Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["BBs"][Symbol["p2"], Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["CCtilde"][Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      +Symbol["CCtilde"][Mass2[loopfields[[1]][{SARAH`gI4}]],Mass2[loopfields[[2]][{SARAH`gI5}]],Mass2[loopfields[[2]][{SARAH`gI5}]]]};

     If[nFields[[1]] == 1, loopfunctions=loopfunctions//.{x_[{SARAH`gI4}]->x}];
     If[nFields[[2]] == 1, loopfunctions=loopfunctions//.{x_[{SARAH`gI5}]->x}];
     tempexpr = Plus @@ (prefactors*loopfunctions);
     If[nFields[[1]] > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
     If[nFields[[2]] > 1,tempexpr=SARAH`sum[SARAH`gI5,1,nFields[[2]],tempexpr];];
     tempexpr
   ,Nothing]
];

GetTadpolesfromNPointFunctions[nPointFunctions_List]:=If[Length[nPointFunctions]==0,Print["No nPointFunctions available."];{},Cases[nPointFunctions,_SelfEnergies`Tadpole]];

(* replaces tadpole[i] expressions with the explicit 1-loop tadpole expression *)
tadpoleReplacementRules[n_, assoc_, tadexpr_] :=
    Module[{},
     If[n == 1,
     {Rule[Symbol["tadpole"][1], Cases[tadexpr,SelfEnergies`Tadpole[assoc[[1, 1]],expr1L_] -> (expr1L)][[1]]]},
      Table[Rule[Symbol["tadpole"][i],
      ReleaseHold[
         Cases[tadexpr,SelfEnergies`Tadpole[assoc[[i, 1]][ind_], expr1L_] ->
            Hold[SARAH`sum[SARAH`gI3,1,n,(expr1L //. {ind :> SARAH`gI3})*Symbol["KroneckerDelta"][SARAH`gI3,assoc[[i, 2]]-1]]]][[1]]
            ]
          ], {i, 1, n}]
        ]
      ];

(* replaces vertices with zero, if they are in fact zero given the specified substiutions *)
gaugelessVertexRules[subs_List] := {Symbol["Cp"][fields__] /; VertexZeroQ[List[fields], subs] -> 0,
                             Symbol["Cp"][fields__][mod_] /; VertexZeroQ[List[fields], subs, mod] -> 0};


(* shamelessly stolen from SelfEnergies` *)
AppendShiftFieldIndices[lst_List, idx__] :=
   AppendShiftFieldIndicesTo[#,idx]& /@ lst;

AppendShiftFieldIndicesTo[f_[field_, ex___], idx__] :=
   If[SARAH`getGen[field] > 1,
      f[field[idx], ex ],
      f[field, ex]
     ];

(* these are the fields, whose selfenergies will get a shift due to the 1L tadpoles *)
RelevantFieldsSelfEnergies:={ToExpression["U"<>ToString[SARAH`HiggsBoson]],
                        ToExpression["U"<>ToString[SARAH`PseudoScalar]],
                        SARAH`HiggsBoson,SARAH`PseudoScalar};

GetHiggsSelfEnergy[energies_List]:=Select[energies,MemberQ[RelevantFieldsSelfEnergies,GetnPointField[#]]&];


(* creates expressions for shifts of 1L selfenergies and tadpoles through
   1L corrections to the EWSB parameters. The sign is opposite to that of
   the values given by CalculatePi2S in SPheno, since there they multiply
   their result with a -1 later on. The signs used here should be consistent
   with both the rest of FlexibleSUSY *)
Make1L2LShifts[SarahTads_List,SarahSelfenergies_List,nPointFuncs_List,tadHiggsassoc_,EWSBpars_List,treeEWSBsol_List,sub_List,eigenstates_] :=
Module[{gaugelesssub,relevantTadpoles,tad1Lexpr,SE1Lexpr,treelevelsolution,higgstoewsb,
     nHiggs,tadpoleSeriesParameters,tadpolefields,massshiftsintadpole,shifts,
     tadsnPointform={},SEnPointform={},nPointform={}},

     gaugelesssub = sub /. {Rule[gaugelesscoup_,probzero_]:>Rule[Symbol[SymbolName[gaugelesscoup]],probzero]}; (* for some reason, the couplings have context FlexibleSUSY`Private` for no obvious reason. This fixes that. *)
     tad1Lexpr = GetTadpolesfromNPointFunctions[nPointFuncs] //. {SelfEnergies`Tadpole[f_,L1_,___]->SelfEnergies`Tadpole[f,L1]};
     tad1Lexpr = tad1Lexpr //. gaugelessVertexRules[gaugelesssub];
     If[Length[tad1Lexpr] == 1,

       relevantTadpoles = GetShiftedExpressions[SarahTads,EWSBpars,eigenstates];
       relevantSelfEnergies = GetShiftedExpressions[GetHiggsSelfEnergy[SarahSelfenergies],EWSBpars,eigenstates];
       (*We don't have to worry about ignoring Goldstones, since those contributions will
         be set to zero in the loop functions (there, e.g. BB(small,small,scale)=0)
         given that their masses have been properly set to zero *)

       treelevelsolution = treeEWSBsol//.gaugelesssub;
       nHiggs = Length[tadHiggsassoc];

       tadpoleSeriesParametersFirstOrder = Sequence @@ Table[{Symbol["tadpole"][i], 0, 1}, {i, 1, nHiggs}];
       tadpoleSeriesParametersZeroOrder = Sequence @@ Table[{Symbol["tadpole"][i], 0, 0}, {i, 1, nHiggs}];
       tadpolefields = (GetnPointField[#] & /@ relevantTadpoles)/.{Symbol["U"<>ToString[SARAH`HiggsBoson]]->SARAH`HiggsBoson,Symbol["U"<>ToString[SARAH`PseudoScalar]]->SARAH`PseudoScalar};
       selfenergyfields = (GetnPointField[#] & /@ relevantSelfEnergies)/.{Symbol["U"<>ToString[SARAH`HiggsBoson]]->SARAH`HiggsBoson,Symbol["U"<>ToString[SARAH`PseudoScalar]]->SARAH`PseudoScalar};

       massshiftsintadpole = Map[TreeMass[#[[2]], eigenstates] &, relevantTadpoles, {2}] //.gaugelesssub /.treelevelsolution;
       massshiftsintadpole = TreeMasses`StripGenerators[massshiftsintadpole,{SARAH`ct1,SARAH`ct2,SARAH`ct3,SARAH`ct4}]; (* get rid of all colour indices and any generators, that might be present *)
       massshiftsintadpole = Map[Normal[Series[#,tadpoleSeriesParametersFirstOrder]]-Normal[Series[#,tadpoleSeriesParametersZeroOrder]]& , massshiftsintadpole, {2}]; (* subbing the treelevel solution into the masses and throwing out the leading order part *)

       massshiftsinSE = Map[{TreeMass[#[[1]], eigenstates],TreeMass[#[[2]], eigenstates]} &, relevantSelfEnergies, {2}] //.gaugelesssub /.treelevelsolution;
       massshiftsinSE = TreeMasses`StripGenerators[massshiftsinSE,{SARAH`ct1,SARAH`ct2,SARAH`ct3,SARAH`ct4}]; (* get rid of all colour indices and any generators, that might be present *)
       massshiftsinSE = Map[Normal[Series[#,tadpoleSeriesParametersFirstOrder]]-Normal[Series[#,tadpoleSeriesParametersZeroOrder]]& , massshiftsinSE, {3}];

       tadpoleshifts = Plus @@@ (Thread[noEvalfunc[relevantTadpoles,massshiftsintadpole]]//.{noEvalfunc[pars___]->Calc1L2LTadShiftExpr[pars]}); (* Function evaluation with a list as parameter has higher priority than the distribution of lists via Thread, therefore I am using this workaround. Not pretty, but gets the job done. *)
       SEshifts = Plus @@@ (Thread[noEvalfunc[relevantSelfEnergies,massshiftsinSE]]//.{noEvalfunc[pars___]->Calc1L2LSEShiftExpr[pars]});

       tadsnPointform = Thread[SelfEnergies`TadpoleShift1L[tadpolefields,0,tadpoleshifts]]//.tadpoleReplacementRules[nHiggs,tadHiggsassoc,tad1Lexpr]/. SARAH`Mass -> FlexibleSUSY`M;
       SEnPointform = Thread[SelfEnergies`FSSelfEnergyShift1L[selfenergyfields,0,SEshifts]]//.tadpoleReplacementRules[nHiggs,tadHiggsassoc,tad1Lexpr]/. SARAH`Mass -> FlexibleSUSY`M;

       nPointform=Join[nPointform,AppendShiftFieldIndices[tadsnPointform,SARAH`gO1],AppendShiftFieldIndices[SEnPointform,SARAH`gO1,SARAH`gO2]];,
       Print["More than one tadpole present, throwing the towel, no tadpole shifts were calculated. Fix SelfEnergies2L to handle this case."];
     ];
     nPointform
];

CreateEnterGauglessLimitFunction[brokencouplings_]:=Module[{output="",goldstones=SARAH`GoldstoneGhost,goldstonemasses,
   couplingnames=Symbol[SymbolName[#]]& /@ brokencouplings},
   goldstonemasses=FlexibleSUSY`M[#]& /@ (Transpose[goldstones][[2]] /. {field_[{Ind_}] -> field[Ind-1]}) ;

   For[nm=1,nm <= Length[couplingnames], nm++,
      output = output <> CConversion`RValueToCFormString[couplingnames[[nm]]] <> " = 0;\n";
   ];

   output = output <> "solve_ewsb_tree_level();\ncalculate_DRbar_masses();\n";

   For[nm=1,nm <= Length[goldstonemasses], nm++,
      output = output <> CConversion`RValueToCFormString[goldstonemasses[[nm]]] <> " = 0;\n";
   ];

   output
];

End[];
Protect["SelfEnergies2L`*"];
EndPackage[];
