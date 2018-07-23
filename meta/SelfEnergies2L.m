BeginPackage["SelfEnergies2L`", {"SARAH`","TreeMasses`","CConversion`","EWSB`"}];

Unprotect["SelfEnergies2L`*"];
ClearAll["SelfEnergies2L`*"];
ClearAll["SelfEnergies2L`Private`*"];


ConvertSarah2LDiagramList::usage = "Converts SARAH's list with 2-loop
 self-energy and tadpole diagrams."

CreateCHKZEROMULTWrapper::usage = "Adds CHKZEROMULT Macro to expressions with 2-loop functions for improved efficiency"

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


DistributeIndices[{}, coupling_] := coupling;
DistributeIndices[indices_List, c : C[particles__]] := ReplaceFirst[c, Rule[#1, AppendIndex[#1, #2]]& @@@ indices];


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

simpleUnrotateRules:={Symbol["U"<>ToString[SARAH`HiggsBoson]][x_]->SARAH`HiggsBoson[x],
                        Symbol["U"<>ToString[SARAH`PseudoScalar]][x_]->SARAH`PseudoScalar[x]};


AddSEMomDep[] := {Symbol["WfSSSS"][masses__] -> Symbol["WfSSSS"][p^2,masses],  Symbol["XfSSS"][masses__] -> Symbol["XfSSS"][p^2,masses],  Symbol["YfSSSS"][masses__] -> Symbol["YfSSSS"][p^2,masses],
                  Symbol["SfSSS"][masses__] -> Symbol["SfSSS"][p^2,masses],  Symbol["UfSSSS"][masses__] -> Symbol["UfSSSS"][p^2,masses],  Symbol["VfSSSSS"][masses__] -> Symbol["VfSSSSS"][p^2,masses],
                  Symbol["WfSSSV"][masses__] -> Symbol["WfSSSV"][p^2,masses],  Symbol["MfSSSSV"][masses__] -> Symbol["MfSSSSV"][p^2,masses], Symbol["WfSSFF"][masses__] -> Symbol["WfSSFF"][p^2,masses],
                  Symbol["WfSSFbFb"][masses__] -> Symbol["WfSSFbFb"][p^2,masses],  Symbol["MfFbFbFbFbS"][masses__] -> Symbol["MfFbFbFbFbS"][p^2,masses],  Symbol["MfFFbFbFS"][masses__] -> Symbol["MfFFbFbFS"][p^2,masses],
                  Symbol["MfFFbFFbS"][masses__] -> Symbol["MfFFbFFbS"][p^2,masses],  Symbol["MfFFFbFbS"][masses__] -> Symbol["MfFFFbFbS"][p^2,masses],  Symbol["MfFFFFS"][masses__] -> Symbol["MfFFFFS"][p^2,masses],
                  Symbol["MfSFbSFbFb"][masses__] -> Symbol["MfSFbSFbFb"][p^2,masses],  Symbol["MfSFSFbF"][masses__] -> Symbol["MfSFSFbF"][p^2,masses], Symbol["MfSFSFFb"][masses__] -> Symbol["MfSFSFFb"][p^2,masses],
                  Symbol["VfSSSFbFb"][masses__] -> Symbol["VfSSSFbFb"][p^2,masses],  Symbol["VfSSSFF"][masses__] -> Symbol["VfSSSFF"][p^2,masses],  Symbol["VfFbFbFbFbS"][masses__] -> Symbol["VfFbFbFbFbS"][p^2,masses],
                  Symbol["VfFbFFbFS"][masses__] -> Symbol["VfFbFFbFS"][p^2,masses],  Symbol["VfFbFFFbS"][masses__] -> Symbol["VfFbFFFbS"][p^2,masses],  Symbol["VfFFbFbFS"][masses__] -> Symbol["VfFFbFbFS"][p^2,masses],
                  Symbol["VfFFFbFbS"][masses__] -> Symbol["VfFFFbFbS"][p^2,masses],  Symbol["VfFFFFS"][masses__] -> Symbol["VfFFFFS"][p^2,masses],  Symbol["GfFFV"][masses__] -> Symbol["GfFFV"][p^2,masses],
                  Symbol["GfFbFbV"][masses__] -> Symbol["GfFbFbV"][p^2,masses],Symbol["ZfSSSS"][masses__] -> Symbol["ZfSSSS"][p^2,masses],Symbol["MfSSSSS"][masses__] -> Symbol["MfSSSSS"][p^2,masses]};

GetFieldType[x_] := SARAH`getType[x, False, FlexibleSUSY`FSEigenstates];

ReplaceFirst[expr_, rule_Rule] := ReplacePart[expr, (FirstPosition[expr, #1] -> #2) &[Sequence @@ rule]];

ReplaceFirst[expr_, rules_List] := Fold[ReplaceFirst, expr, rules];

ReFields[part_] := part /. {SARAH`bar[x_] -> x, Susyno`LieGroups`conj[x_] -> x, SARAH`Conj[x_] -> x};

AllInternalFieldsQ[fieldslist_List] := ContainsNone[fieldslist /. {fd_[{indx__}] -> indx},{SARAH`gE1,SARAH`gE2}];

FourScalarFieldsQ[fieldslist_List] := (Length[fieldslist] === 4) && (GetFieldType /@ fieldslist === {S,S,S,S});

AllUnbrokenIndicesQ[fieldslist_List,unbrokesymmetries_List] :=
   Module[{modelParticles = SARAH`Particles[FlexibleSUSY`FSEigenstates],fieldIndices,indextypelist},
   fieldIndices = Map[Function[par,Select[modelParticles,#[[1]] === par &]], ReFields[fieldslist] /. {y_[{__}] -> y}][[All, 1, 5]];
   indextypelist = (Transpose /@ (fieldIndices/. {{} -> {{noIndex, 1}}}))[[All, 1]];
   Or @@ Map[Function[indtype, And @@ (MemberQ[#, indtype] & /@ indextypelist)], unbrokesymmetries]
];

(* code for finding unbroken symmetries inspired by SPhenoCoupling.m *)
MarkColorSummableScalarVertices[] :=
   Module[{unbrokesymmetries = SARAH`Gauge[[#, 3]]& /@ (Position[SARAH`SGauge /. A_[{b__}] -> A, #][[1, 1]] & /@
   Select[SARAH`SGauge /. A_[{b__}] -> A, FreeQ[Particles[FlexibleSUSY`FSEigenstates], #] == False &])},
   {Cp[fieldlist__] /; And[FourScalarFieldsQ[List[fieldlist]],AllInternalFieldsQ[List[fieldlist]],AllUnbrokenIndicesQ[List[fieldlist],unbrokesymmetries]] -> Cp[fieldlist][Symbol["L2"]]}
   ];


(* result as from CalculatePi2S in SPheno except for a global (-1),
   to be consistent with FS 1L expressions. See comment at function Make1L2LShifts
   at the bottom. *)
ConvertSarah2LDiagramList[tad_List, head_:Total] :=
    head[SumTadpoleType /@ tad]*(-1) //. {
        (m : (SARAH`Mass | SARAH`Mass2))[(SARAH`bar | Susyno`LieGroups`conj)[p_], idx___] :> m[p, idx],
        (m : (SARAH`Mass | SARAH`Mass2))[p_, idx__] :> m[p[{idx}]],
        C[p__] :> Cp[p]} /. MarkColorSummableScalarVertices[] /. {Symbol["i1"]->SARAH`gI1, Symbol["i2"]->SARAH`gI2, Symbol["i3"]->SARAH`gI3,
        Symbol["i4"]->SARAH`gI4,Symbol["i5"]->SARAH`gI5,Symbol["i6"]->SARAH`gI6} /. UnrotateRules[] /. AddSEMomDep[];




(* given a 1-loop nPoint list of diagrams, determines the field that these diagrams correspond to *)
GetnPointField[tempdiags_List]:=Module[{testdiag=tempdiags[[1]],coupfields,loopfields,tempnPointField},
  loopfields={testdiag[[1]],testdiag[[2]]}/. {Conj[a_] -> a, conj[a_] -> a, bar[b_] -> b}/.{x_[b_]->x};
  coupfields=testdiag[[3]]/.{Cp[tempfields__]->List[tempfields]}/. {Conj[a_] -> a, conj[a_] -> a, bar[b_] -> b}/.{x_[b_]->x};
  For[i = 1, i <= Length[loopfields], i++, coupfields = Delete[coupfields, FirstPosition[coupfields, loopfields[[i]]]];];
  coupfields=DeleteDuplicates[coupfields];
  If[Length[coupfields] != 1,Print["Problem while determining, which field the diagram "];
                           Print[testdiag];
                           Print["corresponds to."]
    ];
  coupfields[[1]]
  ];

GetUsedParameters[expr_,extraPars_:{}] := Select[Join[Transpose[parameters][[1]],extraPars], ! FreeQ[expr, #] &];

GetMassShiftedExpressions[SarahList_List,ewsbEqparameters_,EWSBSubst_,eigenstates_]:=
Select[#, With[{unused = #},ContainsAny[Join[GetUsedParameters[TreeMass[#[[1]], eigenstates] /. EWSBSubst,ewsbEqparameters],
      GetUsedParameters[TreeMass[#[[2]], eigenstates] /. EWSBSubst],ewsbEqparameters],ewsbEqparameters]] &] & /@ SarahList;
    (* the 'unused = #'' is set, so that we can use a pure function inside another one *)


VertexZeroQ[fields_List,subs_List,mod_:-1]:= If[mod===-1,
      SameQ[0, Simplify[(Drop[Vertex[fields//.simpleUnrotateRules],1][[1,1]])//.subs]],
      SameQ[0, Simplify[(Select[Drop[Vertex[fields//.simpleUnrotateRules],1],#[[2]]===mod&][[1,1]])//.subs]]];

(*TadpoleReplacements1L2L:={Cp[a___, Symbol["Uhh"][{SARAH`gI3}] , b___]*Symbol["tadpole"][i_Integer] :> Cp[a, Symbol["Uhh"][{i}], b]*Symbol["tadpole"][i],
Cp[a___, Symbol["Uhh"][{SARAH`gI3}] , b___][c___]*Symbol["tadpole"][i_Integer] :> Cp[a, Symbol["Uhh"][{i}], b][c]*Symbol["tadpole"][i]};*)


(* we use the convention, that the mass shifts have the indices gI1,gI2, and the extra coupling(s) have gI4,gI5, if there are multiple generations *)
(* The contraction of tadpole[i] with the correct tadpole expression happens via a Kronckerdelta(i,gI3), where*)
(* gI3 is always going to be the index of the Higgs tadpole we are attaching *)


(* Note: the A0=-x(log(x)-1) function used by FlexibleSUSY is the same as in SPheno (though they both differ from the paper, where A=x(log(x)-1))*)
Calc1L2LTadShiftExpr[diags_,massshifts_List] := Module[{result={}},
  For[n = 1,n <= Length[diags],n++,
    Switch[diags[[n,4]],
      Symbol["SSS"],
        result = Append[result,Calc1L2LTadShiftSSS[diags[[n]],massshifts[[n]]]];,
      Symbol["FFS"],
        result = Append[result,Calc1L2LTadShiftFFS[diags[[n]],massshifts[[n]]]];,
      _,Null;];
  ];
  result
];

(* Tadpole shift functions for the different topologies *)

Calc1L2LTadShiftSSS[diag_List,massshift_]:=Module[{tempexpr,loopfield,loopfunction,nField,couplings,prefactor},
    loopfield = diag[[2]]/.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nField = TreeMasses`GetDimension[loopfield];

    couplings = diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};

    prefactor = 4*0.5*diag[[5]]*diag[[6]]*couplings*massshift//.{Symbol["generation"]->SARAH`gI4}; (* SARAH diag basefactor: 4 *)

    loopfunction =-Symbol["BB"][SARAH`Mass2[loopfield[{SARAH`gI4}]],SARAH`Mass2[loopfield[{SARAH`gI4}]]]; (* -BB is P_{SS} *)

    If[nField == 1, loopfunction=loopfunction//.{x_[{SARAH`gI4}]->x}];
    tempexpr = prefactor*loopfunction;
    If[nField > 1, tempexpr = SARAH`sum[SARAH`gI4,1,nField,tempexpr];];
    tempexpr
];

Calc1L2LTadShiftFFS[diag_List,massshift_]:=Module[{tempexpr,loopfield,loopfunctions,nField,couplings,prefactors},
  If[FlexibleSUSY`Exclude1L2LFermionShifts === False,
    loopfield = diag[[2]]/.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nField = TreeMasses`GetDimension[loopfield];

    couplings = diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};

    prefactors = -4*diag[[5]]*diag[[6]]*couplings*{AbsSqrt[Re[massshift]],
                                                  massshift*SARAH`Mass[loopfield[{SARAH`gI4}]]}/.{Symbol["generation"]->SARAH`gI4};

    loopfunctions = {-Symbol["A0"][SARAH`Mass2[loopfield[{SARAH`gI4}]]],
                    -Symbol["BB"][SARAH`Mass2[loopfield[{SARAH`gI4}]],SARAH`Mass2[loopfield[{SARAH`gI4}]]]};

    If[nField == 1, loopfunctions = loopfunctions /. {x_[{SARAH`gI4}]->x}];
    tempexpr = Plus @@ (prefactors * loopfunctions);
    If[nField > 1,tempexpr = SARAH`sum[SARAH`gI4,1,nField,tempexpr];];
    tempexpr
  ,Nothing]
];


(* calculates the 1L2L shifts to the given selfenergies *)
Calc1L2LSEShiftExpr[diags_,massshifts_List]:=Module[{result={}},
  For[n=1,n<=Length[diags],n++,
    Switch[diags[[n,4]],
      Symbol["SSS"],
        result = Append[result,Calc1L2LselfEnergyShiftsSS[diags[[n]],massshifts[[n]]]];,
      Symbol["SSSS"],
        result = Append[result,Calc1L2LselfEnergyShiftsSSS[diags[[n]],massshifts[[n]]]];,
      Symbol["FFS"],
        result = Append[result,Calc1L2LSEShiftFFS[diags[[n]],massshifts[[n]]]];,
      _,Null;];
  ];
  result
];

(* SelfEnergy shift functions for the different topologies *)

Calc1L2LselfEnergyShiftsSS[diag_List,massshifts_]:=Module[{tempexpr,loopfields,loopfunction,nFields,couplings,prefactors,loopfunctions},
  loopfields = {diag[[1]],diag[[2]]}//.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
  nFields = TreeMasses`GetDimension[#]& /@ loopfields;

  couplings = diag[[3]] * (diag[[3]]/.{SARAH`gO1->SARAH`gO2}/. {Cp[tempfields__]:>Cp[Sequence @@ (AntiField /@ List[tempfields])]})/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI5};

  prefactors = 2*diag[[5]]*diag[[6]]*couplings*{massshifts[[1]]/.{Symbol["generation"]->SARAH`gI4},
                                                massshifts[[2]]/.{Symbol["generation"]->SARAH`gI5}};

  loopfunctions = {-Symbol["CCtilde"][SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]]],
                 -Symbol["CCtilde"][SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]]};

  If[nFields[[1]] == 1, loopfunctions = loopfunctions/.{x_[{SARAH`gI4}]->x}];
  If[nFields[[2]] == 1, loopfunctions = loopfunctions/.{x_[{SARAH`gI5}]->x}];
  tempexpr = Plus @@ (prefactors * loopfunctions);
  If[nFields[[1]] > 1,tempexpr = SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
  If[nFields[[2]] > 1,tempexpr = SARAH`sum[SARAH`gI5,1,nFields[[2]],tempexpr];];
  tempexpr
];

Calc1L2LselfEnergyShiftsSSS[diag_List,massshifts_]:=Module[{tempexpr,loopfields,loopfunction,nFields,couplings,prefactors,loopfunctions},
  If[!(MatchQ[GetnPointField[{diag}],SARAH`PseudoScalar] && FlexibleSUSY`Exclude1L2LAhShiftSSSS === True),
    loopfields = {diag[[2]]}/.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
    nFields = TreeMasses`GetDimension[#]& /@ loopfields;

    couplings= diag[[3]]/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI4};
    If[!FreeQ[couplings, x_[{SARAH`gO1}]],couplings = ReplaceFirst[couplings,SARAH`gO1->SARAH`gO2];];

    prefactors = 2*0.5*diag[[5]]*diag[[6]]*couplings*{massshifts[[1]]/.{Symbol["generation"]->SARAH`gI4}};

    loopfunctions = {-Symbol["BB"][SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]]]};

    If[nFields[[1]] == 1, loopfunctions=loopfunctions//.{x_[{SARAH`gI4}]->x}];
    tempexpr=Plus @@ (prefactors*loopfunctions);
    If[nFields[[1]] > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
    tempexpr
  ,Nothing]
];

Calc1L2LSEShiftFFS[diag_List,massshiftsgen_]:=Module[{tempexpr,loopfields,loopfunction,nFields,
   couplings,prefactors,loopfunctions,massshifts={massshiftsgen[[1]]/.{Symbol["generation"]->SARAH`gI4},massshiftsgen[[2]]/.{Symbol["generation"]->SARAH`gI5}}},
   If[FlexibleSUSY`Exclude1L2LFermionShifts === False,
     loopfields = {diag[[1]],diag[[2]]}//.{bar[x_]->x,conj[x_]->x,Conj[x_]->x};
     nFields = TreeMasses`GetDimension[#]& /@ loopfields;

     couplings = diag[[3]] * (diag[[3]]/.{SARAH`gO1->SARAH`gO2}/. {Cp[tempfields__]:>Cp[Sequence @@ (AntiField /@ List[tempfields])]})/.{SARAH`gI1->SARAH`gI4,SARAH`gI2->SARAH`gI5};
     prefactors = diag[[5]]*diag[[6]]*couplings*{-2*massshifts[[1]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]], (* the first 4 parts correspond to the -2*m1*m2*B(p2,m12,m22) part of the SE *)
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*massshifts[[2]],
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]*2*massshifts[[1]]*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]],
                                                -2*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]*2*massshifts[[2]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]],
                                                +2*massshifts[[1]]*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]],  (* these last parts correspond to the +G0(p2,m12,m22) part of the SE *)
                                                +2*massshifts[[2]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]],
                                                -(2*massshifts[[1]]*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]]+2*massshifts[[2]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]),
                                                -(SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]] + SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]-SARAH`p^2)*2*massshifts[[1]]*SARAH`Mass[loopfields[[1]][{SARAH`gI4}]],
                                                -(SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]] + SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]-SARAH`p^2)*2*massshifts[[2]]*SARAH`Mass[loopfields[[2]][{SARAH`gI5}]]};

     loopfunctions = {+Symbol["BBs"][SARAH`p^2, SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["BBs"][SARAH`p^2, SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      -Symbol["CCtilde"][SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      -Symbol["CCtilde"][SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      -Symbol["BB"][SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      -Symbol["BB"][SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      +Symbol["BBs"][SARAH`p^2, SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]],
                      -Symbol["CCtilde"][SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]]],
                      -Symbol["CCtilde"][SARAH`Mass2[loopfields[[1]][{SARAH`gI4}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]],SARAH`Mass2[loopfields[[2]][{SARAH`gI5}]]]};

     If[nFields[[1]] == 1, loopfunctions=loopfunctions//.{x_[{SARAH`gI4}]->x}];
     If[nFields[[2]] == 1, loopfunctions=loopfunctions//.{x_[{SARAH`gI5}]->x}];
     tempexpr = Plus @@ (prefactors*loopfunctions);
     If[nFields[[1]] > 1,tempexpr=SARAH`sum[SARAH`gI4,1,nFields[[1]],tempexpr];];
     If[nFields[[2]] > 1,tempexpr=SARAH`sum[SARAH`gI5,1,nFields[[2]],tempexpr];];
     tempexpr
   ,Nothing]
];

GetTadpolesfromNPointFunctions[nPointFunctions_List]:=If[Length[nPointFunctions]==0,Print["No nPointFunctions available."];{},Cases[nPointFunctions,_SelfEnergies`Tadpole]];
GetHiggsSEfromNPointFunctions[nPointFunctions_List]:=If[Length[nPointFunctions]==0,Print["No nPointFunctions available."];{},Cases[nPointFunctions,SelfEnergies`FSSelfEnergy[field_, __]/;IsRelevantSelfEnergyFieldQ[field]]];


(* replaces tadpole[i] expressions with the explicit 1-loop tadpole expression *)
tadpoleReplacementRules[assoc_, tadexpr_] :=
    Module[{n=Length[assoc]},
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

(* replaces vertices with zero, if they are in fact zero given the specified substitutions *)
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
RelevantselfenergyFields[] := {ToExpression["U"<>ToString[SARAH`HiggsBoson]],
                        ToExpression["U"<>ToString[SARAH`PseudoScalar]],
                        SARAH`HiggsBoson,SARAH`PseudoScalar};

IsRelevantSelfEnergyFieldQ[field_[___]] := MemberQ[RelevantselfenergyFields[],field];
IsRelevantSelfEnergyFieldQ[field_] := MemberQ[RelevantselfenergyFields[],field];


GetHiggsSelfEnergy[energies_List]:=Select[energies,MemberQ[RelevantselfenergyFields[],GetnPointField[#]]&];

ReduceExplicitGenIndices := Module[{generationMatrices = First /@ Select[SARAH`parameters, (SequenceCount[Flatten[#], {generation, generation}] === 1) &]},
                                    {matx_[a_Symbol, b_Integer] /; MemberQ[generationMatrices, matx] -> matx[a, b - 1],
                                     matx_[a_Integer, b_Symbol] /; MemberQ[generationMatrices, matx] -> matx[a - 1, b],
                                     matx_[a_Integer, b_Integer] /; MemberQ[generationMatrices, matx] -> matx[a - 1, b - 1]}
   ];

ReplaceSARAHMassHeads :=
    {
        SARAH`Mass2[p_] /; TreeMasses`IsFermion[p] :> FlexibleSUSY`M[p]^2,
        SARAH`Mass -> FlexibleSUSY`M,
        SARAH`Mass2 -> FlexibleSUSY`M2,
        SARAH`p^2 -> FlexibleSUSY`p2
    };

ExtractFieldName[field_[idx1_,idx2_]] := ExtractFieldName[field];
ExtractFieldName[field_[PL]]          := ExtractFieldName[field];
ExtractFieldName[field_[PR]]          := ExtractFieldName[field];
ExtractFieldName[field_[1]]           := ExtractFieldName[field];
ExtractFieldName[field_[idx_]]        := ExtractFieldName[field];
ExtractFieldName[field_]              := ToValidCSymbolString[field];

StripFieldRotation := {Symbol["U"<>ToString[SARAH`HiggsBoson]]->SARAH`HiggsBoson,Symbol["U"<>ToString[SARAH`PseudoScalar]]->SARAH`PseudoScalar};

(* creates expressions for shifts of 1L selfenergies and tadpoles through
   1L corrections to the EWSB parameters. The sign is opposite to that of
   the values given by CalculatePi2S in SPheno. The signs used here should be
   consistent with the rest of FlexibleSUSY *)
Make1L2LShifts[SarahTads_List,SarahSelfenergies_List,nPointFuncs_List,tadHiggsassoc_,EWSBpars_List,treeEWSBsol_List,sub_List,EWSBSubst_List,eigenstates_] :=
Module[{gaugelessSub,relevantMassTadpoles,relevantMassSelfEnergies,tadpole1LoopExpr,selfEnergy1LoopExpr,treeSol,higgstoewsb,
     nHiggs,tadpoleSeriesParameters,tadpoleFields,tadpoleShifts,selfEnergyShifts,couplingTadpoleShift = 0,couplingSelfEnergyShift = 0,
     massTadpoleShift,massSelfEnergyShift,tadpoleNPointForm={},selfEnergyNPointForm={},nPointform={}},

     gaugelessSub = sub /. {Rule[gaugelesscoup_,probzero_]:>Rule[Symbol[SymbolName[gaugelesscoup]],probzero]}; (* for some reason, the couplings have context FlexibleSUSY`Private` for no obvious reason. This fixes that. *)
     tadpole1LoopExpr = GetTadpolesfromNPointFunctions[nPointFuncs] /. {SelfEnergies`Tadpole[f_,L1_,___]->SelfEnergies`Tadpole[f,L1]};
     tadpole1LoopExpr = tadpole1LoopExpr /. gaugelessVertexRules[gaugelessSub];
     selfEnergy1LoopExpr = GetHiggsSEfromNPointFunctions[nPointFuncs] /. {SelfEnergies`FSSelfEnergy[f_,L1_,___]->SelfEnergies`SelfEnergies`FSSelfEnergy[f,L1]};
     selfEnergy1LoopExpr = selfEnergy1LoopExpr /. gaugelessVertexRules[gaugelessSub];

     If[Length[tadpole1LoopExpr] == 1,

       tadpoleFields = (GetnPointField[#] & /@ SarahTads) /. StripFieldRotation;
       selfenergyFields = (GetnPointField[#] & /@ GetHiggsSelfEnergy[SarahSelfenergies]) /. StripFieldRotation;

       relevantMassTadpoles = GetMassShiftedExpressions[SarahTads,EWSBpars,EWSBSubst,eigenstates];
       relevantMassSelfEnergies = GetMassShiftedExpressions[GetHiggsSelfEnergy[SarahSelfenergies],EWSBpars,EWSBSubst,eigenstates];
       (*We don't have to worry about ignoring Goldstones, since those contributions will
         be set to zero in the loop functions (there, e.g. BB(small,small,scale)=0)
         given that their masses have been properly set to zero *)

       treeSol = EWSB`ReplaceFixedParametersBySymbolsInTarget[treeEWSBsol /. gaugelessSub /. EWSBSubst];
       nHiggs = Length[tadHiggsassoc];



       massTadpoleShift = GenerateTadpoleMassShifts[relevantMassTadpoles,gaugelessSub,treeSol,nHiggs,EWSBSubst,eigenstates];
       massSelfEnergyShift = GenerateSelfEnergyMassShifts[relevantMassSelfEnergies,gaugelessSub,treeSol,nHiggs,EWSBSubst,eigenstates];

       (* couplingTadpoleShift = GenerateCouplingShifts[tadpole1LoopExpr,gaugelessSub,treeSol,nHiggs,eigenstates];
       couplingSelfEnergyShift = GenerateCouplingShifts[selfEnergy1LoopExpr,gaugelessSub,treeSol,nHiggs,eigenstates]; *)

       (* couplingTadpoleShift = OrderingToTarget[Replace[couplingTadpoleShift,{_,expr_}->expr,{1}],Replace[couplingTadpoleShift,{field_,_}->ExtractFieldName[field],{1}],tadpoleFields];
       couplingSelfEnergyShift = OrderingToTarget[Replace[couplingSelfEnergyShift,{_,expr_}->expr,{1}],Replace[couplingSelfEnergyShift,{field_,_}->ExtractFieldName[field],{1}],selfenergyFields]; *)

       tadpoleShifts = massTadpoleShift + couplingTadpoleShift;
       selfEnergyShifts = massSelfEnergyShift + couplingSelfEnergyShift;

       tadpoleNPointForm = Thread[SelfEnergies`TadpoleShift1L[tadpoleFields,0,tadpoleShifts]]/.tadpoleReplacementRules[tadHiggsassoc,tadpole1LoopExpr]/. SARAH`Mass -> FlexibleSUSY`M  /. {xy_^(-1/2) -> 1/AbsSqrt[xy], Sqrt[xy_] -> AbsSqrt[xy]};
       selfEnergyNPointForm = Thread[SelfEnergies`FSSelfEnergyShift1L[selfenergyFields,0,selfEnergyShifts]]/.tadpoleReplacementRules[tadHiggsassoc,tadpole1LoopExpr]/. SARAH`Mass -> FlexibleSUSY`M  /. {xy_^(-1/2) -> 1/AbsSqrt[xy], Sqrt[xy_] -> AbsSqrt[xy]};

       tadpoleNPointForm = tadpoleNPointForm /. ReduceExplicitGenIndices /. ReplaceSARAHMassHeads;
       selfEnergyNPointForm = selfEnergyNPointForm /. ReduceExplicitGenIndices /. ReplaceSARAHMassHeads;

       nPointform=Join[nPointform,AppendShiftFieldIndices[tadpoleNPointForm,SARAH`gO1],AppendShiftFieldIndices[selfEnergyNPointForm,SARAH`gO1,SARAH`gO2]];,
       Print["Whacky model with more than one tadpole. Tell someone to fix SelfEnergies2L to handle this case."];
       Quit[];
     ];
     nPointform
];

ReplaceMassInternalIndices := {SARAH`sum[idx_, bndLow_, bndHigh_, Expr_] /;
   StringMatchQ[ToString[idx],
    RegularExpression["j[0-9]*"]] :> (SARAH`sum[
    ToExpression[
     StringReplace[
      ToString[idx], {"j" -> "SARAH`gI",
       idxnum : DigitCharacter .. :>
        ToString[ToExpression[idxnum] + 5]}]], bndLow, bndHigh,
    Expr /. {idx ->
       ToExpression[
        StringReplace[
         ToString[idx], {"j" -> "SARAH`gI",
          idxnum : DigitCharacter .. :>
           ToString[ToExpression[idxnum] + 5]}]]}])};

GenerateCouplingShifts[selfEnergiesFormat_,gaugelessSub_,treeSol_,nHiggs_,eigenstates_] := Module[
   {},

   tadpoleSeriesParametersFirstOrder = Sequence @@ Table[{Symbol["tadpole"][i], 0, 1}, {i, 1, nHiggs}];
   tadpoleSeriesParametersZeroOrder = Sequence @@ Table[{Symbol["tadpole"][i], 0, 0}, {i, 1, nHiggs}];

   (* Necessary output format: {{field,expr},{field,expr},...} *)
];


(* the following command plus the 'SetOptions[D,...]' are necessary to get SARAH`sum to work correctly
   with differentiation *)
SARAH`sum /: D[SARAH`sum[idx_, a_, b_, exprs_], y_, c___] :=
   SARAH`sum[idx, a, b, D[exprs, y, c]];


FirstOrderSeries[seriesExpr_, seriesPars_List] :=
    Plus @@ (((D[seriesExpr, #] & /@ seriesPars[[All, 1]]) /. (Rule[#[[1]], #[[2]]] & /@ seriesPars))*seriesPars[[All, 1]]);

GenerateTadpoleMassShifts[relevantMassTadpoles_,gaugelessSub_,treeSol_,nHiggs_,EWSBSubst_,eigenstates_] := Module[{massshiftsTadpole,
   tadpoleSeriesParameters,makeParametersUnique},

   makeParametersUnique = EWSB`MakeParametersUnique[FlexibleSUSY`EWSBOutputParameters];

   tadpoleSeriesParameters =  Table[{Symbol["tadpole"][i], 0}, {i, 1, nHiggs}];

   massshiftsTadpole = Map[(TreeMass[#[[2]], eigenstates] /. EWSBSubst) &, relevantMassTadpoles, {2}];
   massshiftsTadpole = massshiftsTadpole /. makeParametersUnique;
   massshiftsTadpole = massshiftsTadpole //. ReplaceMassInternalIndices /. gaugelessSub /. treeSol;

   massshiftsTadpole = TreeMasses`StripGenerators[massshiftsTadpole,{SARAH`ct1,SARAH`ct2,SARAH`ct3,SARAH`ct4}]; (* get rid of all colour indices and any generators, that might be present *)
   SetOptions[D, NonConstants -> {SARAH`sum}];
   massshiftsTadpole = Map[FirstOrderSeries[#,tadpoleSeriesParameters]& , massshiftsTadpole, {2}]; (* subbing the treelevel solution into the masses and throwing out the leading order part *)
   SetOptions[D, NonConstants -> {}];

   makeParametersUnique = Reverse /@ makeParametersUnique;
   massshiftsTadpole = massshiftsTadpole /. makeParametersUnique;

   Plus @@@ (Thread[noEvalfunc[relevantMassTadpoles,massshiftsTadpole]]/.{noEvalfunc[pars___]->Calc1L2LTadShiftExpr[pars]}) (* Function evaluation with a list as parameter has higher priority than the distribution of lists via Thread, therefore I am using this workaround. Not pretty, but gets the job done. *)
];

GenerateSelfEnergyMassShifts[relevantMassSelfEnergies_,gaugelessSub_,treeSol_,nHiggs_,EWSBSubst_,eigenstates_] := Module[{massshiftsSelfEnergy,
   tadpoleSeriesParameters,makeParametersUnique},

   makeParametersUnique = EWSB`MakeParametersUnique[FlexibleSUSY`EWSBOutputParameters];

   tadpoleSeriesParameters =  Table[{Symbol["tadpole"][i], 0}, {i, 1, nHiggs}];

   massshiftsSelfEnergy = Map[{TreeMass[#[[1]], eigenstates] /. EWSBSubst,TreeMass[#[[2]], eigenstates] /. EWSBSubst} &, relevantMassSelfEnergies, {2}];
   massshiftsSelfEnergy = massshiftsSelfEnergy /. makeParametersUnique;
   massshiftsSelfEnergy = massshiftsSelfEnergy //. ReplaceMassInternalIndices /. gaugelessSub /. treeSol;

   massshiftsSelfEnergy = TreeMasses`StripGenerators[#,{SARAH`ct1,SARAH`ct2,SARAH`ct3,SARAH`ct4}] & /@ massshiftsSelfEnergy; (* get rid of all colour indices and any generators, that might be present *)
   SetOptions[D, NonConstants -> {SARAH`sum}];
   massshiftsSelfEnergy = Map[FirstOrderSeries[#,tadpoleSeriesParameters]& , massshiftsSelfEnergy, {3}];
   SetOptions[D, NonConstants -> {}];

   makeParametersUnique = Reverse /@ makeParametersUnique;
   massshiftsSelfEnergy = massshiftsSelfEnergy /. makeParametersUnique;

   Plus @@@ (Thread[noEvalfunc[relevantMassSelfEnergies,massshiftsSelfEnergy]]/.{noEvalfunc[pars___]->Calc1L2LSEShiftExpr[pars]})
];

OrderingToTarget[list_, sourceIds_, targetIds_] :=
  list[[Ordering @ sourceIds]][[Ordering @ Ordering@targetIds]];

CreateEnterGauglessLimitFunction[brokencouplings_]:=Module[{output="",
   couplingnames = Symbol[SymbolName[#]]& /@ brokencouplings},

   For[nm = 1,nm <= Length[couplingnames], nm++,
      output = output <> CConversion`RValueToCFormString[couplingnames[[nm]]] <> " = 0;\n";
   ];

   output = output <> "if(ewsb_solve_consistently){solve_ewsb_tree_level();}\n";
   output = output <> "problems.set_suppress_running_tachyon_warning(true);\n";
   output = output <> "calculate_DRbar_masses();\n";
   output = output <> "calculate_vertices();";
   output
];

NotLFFreeQ[x__] := And @@ (! FreeQ[#, Alternatives @@ loopFunctions] & /@ (List[x]));
LFFreeQ[x__] :=  And @@ (FreeQ[#, Alternatives @@ loopFunctions] & /@ (List[x]));
AllNumericQ[x__] := And @@ (NumericQ[#] & /@ List[x]);

addCHKZEROMULTWrapper := {Times[b__,aft__] /; (NotLFFreeQ[b] &&
   LFFreeQ[aft] && !AllNumericQ[aft]) :>
   FlexibleSUSY`CHKZEROMULT[Times[aft],Times[b] /. addCHKZEROMULTWrapper]};

loopFunctions := {Symbol["TfSS"], Symbol["TfSSS"], Symbol["TfSSSS"],
   Symbol["TfSSFF"], Symbol["TfSSFbFb"], Symbol["TfFFFbS"],
   Symbol["TfFFbFS"], Symbol["TfFbFbFbS"], Symbol["TfSV"],
   Symbol["TfFV"], Symbol["WfSSSS"], Symbol["XfSSS"],
   Symbol["YfSSSS"], Symbol["SfSSS"], Symbol["UfSSSS"],
   Symbol["VfSSSSS"], Symbol["WfSSSV"], Symbol["MfSSSSV"],
   Symbol["WfSSFF"], Symbol["WfSSFbFb"], Symbol["MfFbFbFbFbS"],
   Symbol["MfFFbFbFS"], Symbol["MfFFbFFbS"], Symbol["MfFFFbFbS"],
   Symbol["MfFFFFS"], Symbol["MfSFbSFbFb"], Symbol["MfSFSFbF"],
   Symbol["MfSFSFFb"], Symbol["VfSSSFbFb"], Symbol["VfSSSFF"],
   Symbol["VfFbFbFbFbS"], Symbol["VfFbFFbFS"], Symbol["VfFbFFFbS"],
   Symbol["VfFFbFbFS"], Symbol["VfFFFbFbS"], Symbol["VfFFFFS"],
   Symbol["GfFFV"], Symbol["GfFbFbV"], Symbol["ZfSSSS"],
   Symbol["MfSSSSS"], Symbol["CCtilde"], Symbol["BB"], Symbol["BBs"]};

splitLoopFunctionSum = {SARAH`sum[idx_, start_, stop_,
   Plus[c___, a_, b_,d___]] /; (!FreeQ[a, Alternatives[Sequence @@ loopFunctions], Heads -> True] &&
   !FreeQ[b, Alternatives[Sequence @@ loopFunctions], Heads -> True]) :> SARAH`sum[idx, start, stop, Plus[c, a]] + SARAH`sum[idx, start, stop, Plus[b, d]]};

distributeNumericFactors = {a_?NumericQ*Plus[b_, c__] :> Plus[a*b, a*Plus[c]]};

CreateCHKZEROMULTWrapper[wrapExpr_] := wrapExpr //. splitLoopFunctionSum //. distributeNumericFactors /. addCHKZEROMULTWrapper;


End[];
Protect["SelfEnergies2L`*"];
EndPackage[];
