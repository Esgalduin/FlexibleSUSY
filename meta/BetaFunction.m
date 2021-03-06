(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["BetaFunction`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "Traces`", "Utils`"}];

BetaFunction[];

ConvertSarahRGEs::usage="converts SARAH's beta functions";
CreateBetaFunction::usage="";
GetAllBetaFunctions::usage="";
CountNumberOfParameters::usage="";
CreateDisplayFunction::usage="";
CreateSetFunction::usage="";
CreateSetters::usage="";
CreateGetters::usage="";
CreateParameterDefinitions::usage="";
CreateCCtorInitialization::usage="";
CreateCCtorParameterList::usage="";
CreateParameterList::usage="";
ClearParameters::usage="";

CreateParameterEnum::usage="";
CreateParameterNames::usage="";

GetName::usage="returns parameter name from beta function";
GetBeta::usage="returns beta function expression";

CreateSingleBetaFunctionDecl::usage="";
CreateSingleBetaFunctionDefs::usage="";

Begin["`Private`"];

GetName[BetaFunction[name_, type_, beta_List]] :=
    name /. Parameters`StripSARAHIndicesRules[1] /.
            Parameters`StripSARAHIndicesRules[2] /.
            Parameters`StripSARAHIndicesRules[3] /.
            Parameters`StripSARAHIndicesRules[4];

GetType[BetaFunction[name_, type_, beta_List]] := type;

GetBeta[BetaFunction[name_, type_, beta_List]] := beta;

GetBeta[BetaFunction[name_, type_, beta_List], loopOrder_Integer] :=
    If[Length[beta] < loopOrder, 0, beta[[loopOrder]]];

GetAllBetaFunctions[BetaFunction[name_, type_, beta_List]] := beta;

betaIndices = {
    Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4
};

IsBetaIdx[i_] := MemberQ[betaIndices, i];

GuessType[sym_[_?IsBetaIdx, _?IsBetaIdx]] :=
    Parameters`GetType[sym];

GuessType[sym_[_?IsBetaIdx]] :=
    Parameters`GetType[sym];

GuessType[sym_] :=
    Parameters`GetType[sym];

CreateSingleBetaFunctionDecl[par_, loops_] :=
    CConversion`CreateCType[GetType[par]] <>
    " calc_beta_" <> CConversion`ToValidCSymbolString[GetName[par]] <>
    "_" <> ToString[loops] <> "_loop(const TRACE_STRUCT_TYPE&) const;\n";

CreateSingleBetaFunctionDecl[betaFun_List] :=
    StringJoin[
        (
            CreateSingleBetaFunctionDecl[#, 1] <> 
            CreateSingleBetaFunctionDecl[#, 2] <> 
            CreateSingleBetaFunctionDecl[#, 3] <> 
            CreateSingleBetaFunctionDecl[#, 4]
        )& /@ betaFun
    ];

CreateSingleBetaFunctionDefs[betaFun_List, templateFile_String, sarahTraces_List] :=
    Module[{b, para, type, paraStr, typeStr, files = {},
            inputFile, outputFile,
            localDeclOneLoop, localDeclTwoLoop, localDeclThreeLoop, localDeclFourLoop,
            betaOneLoop, betaTwoLoop, betaThreeLoop, betaFourLoop},
           For[b = 1, b <= Length[betaFun], b++,
               para = GetName[betaFun[[b]]];
               type = GetType[betaFun[[b]]];
               paraStr = CConversion`ToValidCSymbolString[para];
               typeStr = CConversion`CreateCType[type];
               inputFile  = FileNameJoin[{FlexibleSUSY`$flexiblesusyTemplateDir, templateFile}];
               outputFile = FileNameJoin[{FlexibleSUSY`FSOutputDir,
                                          FlexibleSUSY`FSModelName <> "_" <>
                                          StringReplace[templateFile,
                                                        {".cpp.in" -> paraStr <> ".cpp"}]}];
               {localDeclOneLoop, betaOneLoop} = CreateBetaFunction[betaFun[[b]], 1, sarahTraces];
               {localDeclTwoLoop, betaTwoLoop} = CreateBetaFunction[betaFun[[b]], 2, sarahTraces];
               {localDeclThreeLoop, betaThreeLoop} = CreateBetaFunction[betaFun[[b]], 3, sarahTraces];
               {localDeclFourLoop, betaFourLoop}   = CreateBetaFunction[betaFun[[b]], 4, sarahTraces];
               WriteOut`ReplaceInFiles[{{inputFile, outputFile}},
                     { "@ModelName@"     -> FlexibleSUSY`FSModelName,
                       "@parameterType@" -> typeStr,
                       "@parameterName@" -> paraStr,
                       "@localDeclOneLoop@" -> WrapLines[IndentText[localDeclOneLoop]],
                       "@localDeclTwoLoop@" -> WrapLines[IndentText[localDeclTwoLoop]],
                       "@localDeclThreeLoop@" -> WrapLines[IndentText[localDeclThreeLoop]],
                       "@localDeclFourLoop@"  -> WrapLines[IndentText[localDeclFourLoop]],
                       "@betaOneLoop@"   -> WrapLines[IndentText[betaOneLoop]],
                       "@betaTwoLoop@"   -> WrapLines[IndentText[betaTwoLoop]],
                       "@betaThreeLoop@" -> WrapLines[IndentText[betaThreeLoop]],
                       "@betaFourLoop@"  -> WrapLines[IndentText[betaFourLoop]],
                       "@DateAndTime@"   -> DateString[]
                     } ];
               AppendTo[files, outputFile];
              ];
           Return[files];
          ];

(* expand expression and replace given head (usually Plus) by List *)
ToList[expr_, head_] :=
    Module[{exp = Expand[expr]},
           If[Head[exp] === head,
              List @@ exp,
              {expr}
             ]
          ];

FactorOutLoopFactor[expr_] :=
    Module[{i, prefactors = {CConversion`oneOver16PiSqr, CConversion`twoLoop,
                             CConversion`threeLoop, CConversion`oneOver16PiSqr^4}},
           For[i = 1, i <= Length[prefactors], i++,
               If[Coefficient[expr, prefactors[[i]]] =!= 0,
                  Return[TimeConstrainedSimplify[prefactors[[i]] Expand[expr / prefactors[[i]]]]];
                 ];
              ];
           expr
          ];

TimeConstrainedSimplify[expr_] :=
    TimeConstrained[Simplify[expr], FlexibleSUSY`FSSimplifyBetaFunctionsTimeConstraint, expr];

CollectMatMul[expr_] :=
    TimeConstrained[Collect[expr, SARAH`MatMul[___], TimeConstrainedSimplify],
                    FlexibleSUSY`FSSimplifyBetaFunctionsTimeConstraint,
                    expr];

(* split expression into sub-expressions of given maximum size *)
SplitExpression[expr_, size_Integer] :=
    CollectMatMul[FactorOutLoopFactor /@ (Plus @@@ Utils`SplitList[ToList[expr, Plus], size])];

NeedToSplitExpression[expr_, threshold_Integer] :=
    Length[ToList[expr, Plus]] > threshold;

ConvertSingleExprToC[expr_, type_, target_String] :=
    "const " <> CConversion`CreateCType[type] <> " " <> target <>
    " = " <> CastTo[RValueToCFormString[expr], type] <> ";\n"

TryCreateUnitMatrix[CConversion`MatrixType[_,m_,n_] /; m =!= n] := 1;
TryCreateUnitMatrix[type_] := CConversion`CreateUnitMatrix[type];

ConvertExprToC[expr_, type_, target_String] :=
    Module[{result, splitExpr},
           If[NeedToSplitExpression[expr, FlexibleSUSY`FSMaximumExpressionSize],
              splitExpr = SplitExpression[expr, FlexibleSUSY`FSMaximumExpressionSize];
              result = MapIndexed[
                  ConvertSingleExprToC[
                      #1 * TryCreateUnitMatrix[type] /. {
                          CConversion`UNITMATRIX[r_]^_        :> CConversion`UNITMATRIX[r],
                          CConversion`UNITMATRIXCOMPLEX[r_]^_ :> CConversion`UNITMATRIXCOMPLEX[r]
                      },
                      type, target <> "_" <> ToString[#2[[1]]]
                  ]&,
                  splitExpr
              ];
              result = StringJoin[result] <> "\n" <>
                       target <> " = " <>
                       StringJoin[Riffle[MapIndexed[(target <> "_" <> ToString[#2[[1]]])&, splitExpr], " + "]] <>
                       ";\n";
              ,
              result = target <> " = " <> CastTo[RValueToCFormString[expr], type] <> ";\n";
             ];
           result
          ];

(*
 * Create one-loop and two-loop beta function assignments and local definitions.
 *)
CreateBetaFunction[betaFunction_BetaFunction, loopOrder_Integer, sarahTraces_List] :=
     Module[{beta, betaName, name, betaStr,
             type = ErrorType, localDecl, traceRules, expr, loopFactor},
            Switch[loopOrder,
                   1, loopFactor = CConversion`oneOver16PiSqr;,
                   2, loopFactor = CConversion`twoLoop;,
                   3, loopFactor = CConversion`threeLoop;,
                   _, loopFactor = CConversion`oneOver16PiSqr^loopOrder];
            name      = ToValidCSymbolString[GetName[betaFunction]];
            betaName  = "beta_" <> name;
            type = GetType[betaFunction];
            expr = GetAllBetaFunctions[betaFunction];
            If[loopOrder > Length[expr],
               Return[{"",
                       betaName <> " = " <>
                       RValueToCFormString[CConversion`CreateZero[type]] <> ";\n"}];
              ];
            expr       = expr[[loopOrder]];
            (* convert beta function expressions to C form *)
            beta       = (loopFactor * expr) /.
                            { Kronecker[Susyno`LieGroups`i1,SARAH`i2] :> CreateUnitMatrix[type],
                              a_[Susyno`LieGroups`i1] :> a,
                              a_[Susyno`LieGroups`i1,SARAH`i2] :> a,
                              a_[SARAH`i2,Susyno`LieGroups`i1] :> SARAH`Tp[a] };
            traceRules = Traces`CreateTraceRules[{expr}];
            localDecl  = Traces`CreateLocalCopiesOfTraces[{expr}, "TRACE_STRUCT"];
            beta = beta /. traceRules;
            (* replace SARAH traces in expr *)
            traceRules = Rule[#,ToValidCSymbol[#]]& /@ (Traces`FindSARAHTraces[expr, sarahTraces]);
            beta = beta /. traceRules;
            (* collecting complicated matrix multiplications *)
            beta = loopFactor CollectMatMul[beta / loopFactor];
            (* declare SARAH traces locally *)
            localDecl  = localDecl <> Traces`CreateLocalCopiesOfSARAHTraces[expr, sarahTraces, "TRACE_STRUCT"];
            If[beta === 0 || PossibleZeroQ[beta],
               beta = CConversion`CreateZero[type];
              ];
            beta       = Parameters`DecreaseIndexLiterals[beta];
            betaStr    = ConvertExprToC[beta, type, betaName];
            localDecl  = Parameters`CreateLocalConstRefsForInputParameters[expr] <>
                         localDecl;
            Return[{localDecl, betaStr}];
           ];

DeclareBetaFunction[betaFunction_BetaFunction] :=
    Module[{name, betaName, type, ctype},
           name = ToValidCSymbolString[GetName[betaFunction]];
           betaName = "beta_" <> name;
           type = GetType[betaFunction];
           ctype = CConversion`CreateCType[type];
           ctype <> " " <> CConversion`SetToDefault[betaName, type]
          ];

CreateBetaFunctionCallSequential[betaFunction_BetaFunction, loopOrder_Integer] :=
    Module[{name, betaName, result = ""},
           If[Length[GetAllBetaFunctions[betaFunction]] >= loopOrder,
              name = ToValidCSymbolString[GetName[betaFunction]];
              betaName = "beta_" <> name;
              result = betaName <> " += " <>
                       "calc_" <> betaName <> "_" <> ToString[loopOrder] <> "_loop(TRACE_STRUCT);\n";
             ];
           result
          ];

CreateBetaFunctionCallsSequential[betaFunctions_List, loopOrder_Integer] :=
    StringJoin[CreateBetaFunctionCallSequential[#,loopOrder]& /@ betaFunctions];

CreateBetaFunctionCallParallel[betaFunction_BetaFunction, loopOrder_Integer] :=
    Module[{name, betaName, result = ""},
           If[Length[GetAllBetaFunctions[betaFunction]] >= loopOrder,
              name = ToValidCSymbolString[GetName[betaFunction]];
              betaName = "beta_" <> name;
              result = "auto fut_" <> name <>
                       " = global_thread_pool().run_packaged_task([this, &TRACE_STRUCT](){ return calc_" <>
                       betaName <> "_" <> ToString[loopOrder] <> "_loop(TRACE_STRUCT); });\n";
             ];
           result
          ];

CollectBetaFunctionFutures[betaFunction_BetaFunction, loopOrder_Integer] :=
    Module[{name, betaName, result = ""},
           If[Length[GetAllBetaFunctions[betaFunction]] >= loopOrder,
              name = ToValidCSymbolString[GetName[betaFunction]];
              betaName = "beta_" <> name;
              result = betaName <> " += fut_" <> name <> ".get();\n";
             ];
           result
          ];

CreateBetaFunctionCallsParallel[betaFunctions_List, loopOrder_Integer] :=
    "{\n" <>
    TextFormatting`IndentText[
        StringJoin[CreateBetaFunctionCallParallel[#,loopOrder]& /@ betaFunctions] <> "\n" <>
        StringJoin[CollectBetaFunctionFutures[#,loopOrder]& /@ betaFunctions]
    ] <>
    "\n}\n";

CreateBetaFunctionCalls[betaFunctions_List, loopOrder_Integer, parallel_:False] :=
    If[parallel,
       CreateBetaFunctionCallsParallel[betaFunctions, loopOrder],
       CreateBetaFunctionCallsSequential[betaFunctions, loopOrder]
      ];

CreateBetaFunction[betaFunctions_List] :=
    Module[{allBeta1L, allBeta2L, allBeta3L, allBeta3LParallel, allBeta4L},
           allBeta1L = CreateBetaFunctionCalls[betaFunctions,1];
           allBeta2L = TextFormatting`IndentText @ CreateBetaFunctionCalls[betaFunctions,2];
           allBeta3L = TextFormatting`IndentText @ CreateBetaFunctionCalls[betaFunctions,3];
           allBeta4L = TextFormatting`IndentText @ CreateBetaFunctionCalls[betaFunctions,4];
           (* only parallelization of 3L betas leads to a speed-up *)
           allBeta3LParallel = TextFormatting`IndentText @ CreateBetaFunctionCalls[betaFunctions,3,True];
           StringJoin[DeclareBetaFunction /@ betaFunctions] <> "\n" <>
           "if (loops > 0) {\n" <>
           TextFormatting`IndentText[
               "const auto TRACE_STRUCT = CALCULATE_TRACES(loops);\n\n" <>
               allBeta1L <> "\n" <>
               "if (loops > 1) {\n" <>
               allBeta2L <> "\n" <>
               TextFormatting`IndentText[
                   "if (loops > 2) {\n" <>
                   "#ifdef ENABLE_THREADS\n" <>
                   allBeta3LParallel <>
                   "#else\n" <>
                   allBeta3L <>
                   "#endif\n" <>
                   "\n" <>
                   TextFormatting`IndentText[
                       "if (loops > 3) {\n" <>
                       allBeta4L <>
                       "\n}"
                   ] <>
                   "\n}"
               ] <>
               "\n}"
           ] <>
           "\n}\n"
          ];

(* Converts SARAH beta functions to our own format.
 *
 * SARAH format:
 *   { name, one-loop beta, two-loop beta }
 *
 * Our format:
 *   BetaFunction[name, type, {one-loop beta, two-loop beta}]
 *
 * @param betaFunctions list of SARAH-like formated beta functions
 *)
ConvertSarahRGEs[beta_List] :=
    Module[{lst = {}, k, name, type, expr},
           (* extract all beta functions and guess type *)
           For[k = 1, k <= Length[beta], k++,
               If[Length[beta[[k]] < 2], Continue[];];
               (* beta[[k,1]] == name, beta[[k,2]] == 1 loop beta function *)
               name = beta[[k,1]];
               type = GuessType[name];
               expr = Drop[beta[[k]], 1];
               (* protect tensor products *)
               expr = CConversion`ProtectTensorProducts[#, name]& /@ expr;
               (* simplify expressions *)
               expr = TimeConstrainedSimplify /@ expr;
               AppendTo[lst, BetaFunction[name, type, expr]];
              ];
           lst
          ];

(* count number of parameters in beta functions list *)
CountNumberOfParameters[CConversion`ScalarType[CConversion`realScalarCType]] := 1;
CountNumberOfParameters[CConversion`ScalarType[CConversion`complexScalarCType]] := 2;
CountNumberOfParameters[CConversion`ArrayType[CConversion`realScalarCType, entries_]] := entries;
CountNumberOfParameters[CConversion`ArrayType[CConversion`complexScalarCType, entries_]] := 2 * entries;
CountNumberOfParameters[CConversion`VectorType[CConversion`realScalarCType, entries_]] := entries;
CountNumberOfParameters[CConversion`VectorType[CConversion`complexScalarCType, entries_]] := 2 * entries;
CountNumberOfParameters[CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] := rows * cols;
CountNumberOfParameters[CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] := 2 * rows * cols;
CountNumberOfParameters[CConversion`TensorType[CConversion`realScalarCType, dims__]] := Times[dims];
CountNumberOfParameters[CConversion`TensorType[CConversion`complexScalarCType, dims__]] := 2 * Times[dims];

CountNumberOfParameters[betaFunctions_List] :=
    Module[{num = 0},
           (num += CountNumberOfParameters[GetType[#]])& /@ betaFunctions;
           Return[num];
          ];

(* creating set function *)
CreateSetFunction[betaFunctions_List, parameterNumberOffset_:0] :=
    Module[{set = "", paramCount = parameterNumberOffset, name = "", beta = {},
            type = ErrorType, i, numberOfParameters, assignment = "",
            nAssignments = 0},
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = GetAllBetaFunctions[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = ToValidCSymbolString[GetName[betaFunctions[[i]]]];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           numberOfParameters = CountNumberOfParameters[betaFunctions] + parameterNumberOffset;
           If[paramCount != numberOfParameters,
              Print["Error: CreateSetFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[1];];
           Return[set];
          ];

(* creating display function *)
CreateDisplayFunction[betaFunctions_List, parameterNumberOffset_:0] :=
    Module[{display = "", paramCount = parameterNumberOffset, name = "",
            beta = {}, type = ErrorType, i, numberOfParameters, assignment = "",
            nAssignments = 0},
           numberOfParameters = CountNumberOfParameters[betaFunctions] + parameterNumberOffset;
           For[i = 1, i <= Length[betaFunctions], i++,
               beta = GetAllBetaFunctions[betaFunctions[[i]]];
               type = GetType[betaFunctions[[i]]];
               name = GetName[betaFunctions[[i]]];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               display = display <> assignment;
               paramCount += nAssignments;
              ];
           (* sanity check *)
           If[paramCount != numberOfParameters,
              Print["Error: CreateDisplayFunction: number of parameters does not match: ", paramCount,
                    " != ", numberOfParameters]; Quit[1];];
           Return[display];
          ];

CreateParameterNames[betaFunctions_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[
               Parameters`CreateParameterSARAHNames[GetName[#],GetType[#]]& /@ betaFunctions, ", "];
           "const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {" <>
           result <> "};\n"
          ];

CreateParameterEnum[betaFunctions_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[
               Parameters`CreateParameterEnums[GetName[#],GetType[#]]& /@ betaFunctions, ", "];
           If[Length[betaFunctions] > 0, result = result <> ", ";];
           "enum Parameters : int { " <> result <> "NUMBER_OF_PARAMETERS };\n"
          ];

(* create setters *)
CreateElementSetter[name_String, CConversion`ScalarType[_]] := "";

CreateElementSetter[name_String, type_] :=
    CConversion`CreateInlineElementSetter[name, type];

CreateSetters[betaFunction_BetaFunction] :=
    Module[{setter = "", name = "", type},
           name = ToValidCSymbolString[GetName[betaFunction]];
           type = GetType[betaFunction];
           setter = setter <> CConversion`CreateInlineSetter[name, type];
           setter = setter <> CreateElementSetter[name, type];
           Return[setter];
          ];

CreateSetters[betaFunctions_List] :=
    Module[{setter = ""},
           (setter = setter <> CreateSetters[#])& /@ betaFunctions;
           Return[setter];
          ];

(* create getters *)
CreateGetters[betaFunction_BetaFunction] :=
    Module[{getter = "", name = "", type},
           name = ToValidCSymbolString[GetName[betaFunction]];
           type = GetType[betaFunction];
           getter = getter <> CConversion`CreateInlineGetters[name, name, type];
           Return[getter];
          ];

CreateGetters[betaFunctions_List] :=
    Module[{getter = ""},
           (getter = getter <> CreateGetters[#])& /@ betaFunctions;
           Return[getter];
          ];

(* create parameter definition in C++ class *)
CreateParameterDefinitions[betaFunction_BetaFunction] :=
    Parameters`CreateParameterDefinitionAndDefaultInitialize[{GetName[betaFunction], GetType[betaFunction]}];

CreateParameterDefinitions[betaFunctions_List] :=
    StringJoin[CreateParameterDefinitions /@ betaFunctions];

(* create copy constructor initialization list *)
CreateCCtorInitialization[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def  = def <> ", " <> name <> "(" <> name <> "_)";
           Return[def];
          ];

CreateCCtorInitialization[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> CreateCCtorInitialization[#])& /@ betaFunctions;
           Return[def];
          ];

(* create copy constructor initialization list *)
CreateCCtorParameterList[betaFunction_BetaFunction] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CreateGetterReturnType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           def = def <> ", " <> dataType <> " " <> name <> "_";
           Return[def];
          ];

CreateCCtorParameterList[betaFunctions_List] :=
    Module[{def = "", i},
           For[i = 1, i <= Length[betaFunctions], i++,
               def = def <> CreateCCtorParameterList[betaFunctions[[i]]];
              ];
           Return[def];
          ];

(* create parameter list *)
CreateParameterList[betaFunction_BetaFunction, prefix_] :=
    Module[{def = "", name = "", dataType = ""},
           dataType = CConversion`CreateCType[GetType[betaFunction]];
           name = ToValidCSymbolString[GetName[betaFunction]];
           If[def != "", def = def <> ", "];
           def = def <> prefix <> name;
           Return[def];
          ];

CreateParameterList[betaFunctions_List, prefix_:""] :=
    Module[{def = "", i},
           For[i = 1, i <= Length[betaFunctions], i++,
               If[def != "", def = def <> ", "];
               def = def <> CreateParameterList[betaFunctions[[i]], prefix];
              ];
           Return[def];
          ];

ClearParameter[betaFunction_BetaFunction] :=
    Module[{def, type},
           def = CConversion`SetToDefault[ToValidCSymbolString[GetName[betaFunction]],
                                          GetType[betaFunction]];
           Return[def];
          ];

ClearParameters[betaFunctions_List] :=
    Module[{def = ""},
           (def = def <> ClearParameter[#])& /@ betaFunctions;
           Return[def];
          ];

End[];

EndPackage[];
