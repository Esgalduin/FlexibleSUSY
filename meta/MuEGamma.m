(* ::Package:: *)

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

BeginPackage["MuEGamma`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}];

MuEGammaCreateInterfaceFunctionForField::usage="";
MuEGammaContributingDiagramsForFieldAndGraph::usage="";
MuEGammaContributingGraphs::usage="";

(*Begin["Private`"];*)

(* The graphs that contribute to the EDM are precisely those with three
   external lines given by the field in question, its Lorentz conjugate
   and a photon.
   They are given as a List of undirected adjacency matrices where
    1 is the field itself
    2 is its Lorentz conjugate
    3 is the photon
   and all other indices unspecified. *)
vertexCorrectionGraph = {{0,0,0,1,0,0},
                         {0,0,0,0,1,0},
                         {0,0,0,0,0,1},
                         {1,0,0,0,1,1},
                         {0,1,0,1,0,1},
                         {0,0,1,1,1,0}};
contributingGraphs = {vertexCorrectionGraph};

MuEGammaContributingGraphs[] := contributingGraphs

GetPhoton[] := SARAH`Photon

MuEGammaContributingDiagramsForLeptonPairAndGraph[{inLepton_, outLepton_}, graph_] :=
  Module[{diagrams},
    diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph,
         {1 ->CXXDiagrams`LorentzConjugate[inLepton], 2 -> outLepton,
          3 -> CXXDiagrams`LorentzConjugate[GetPhoton[]]}];

    Select[diagrams,IsDiagramSupported[inLepton,outLepton,graph,#] &]
 ]

IsDiagramSupported[inLepton_,outLepton_,vertexCorrectionGraph,diagram_] :=
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = CXXDiagrams`LorentzConjugate[diagram[[4,3]]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {GetPhoton[],photonEmitter,CXXDiagrams`LorentzConjugate[photonEmitter]},
       Return["(unknown diagram)"]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[True]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[True]];
    
    Return[False];
  ]

MuEGammaCreateInterfaceFunctionForLeptonPair[{inLepton_, outLepton_},gTaggedDiagrams_List] :=
  Module[{prototype,definition,
          numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inLepton],
numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outLepton]},
  
    prototype = "double calculate_" <> CXXNameOfField[inLepton] <> "_to_" 
                  <> CXXNameOfField[outLepton] <> "_gamma" <>
                  "(" <> If[TreeMasses`GetDimension[inLepton] =!= 1,
                           " int generationIndex1, ",
                           " "] <>
                        If[TreeMasses`GetDimension[outLepton] =!= 1,
                           " int generationIndex2, ",
                           " "] <>
                 "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model );";
                 
    definition = "double calculate_" <> CXXNameOfField[inLepton] <> "_to_" 
                  <> CXXNameOfField[outLepton] <> "_gamma" <>
                  "(" <> If[TreeMasses`GetDimension[inLepton] =!= 1,
                           " int generationIndex1, ",
                           " "] <>
                        If[TreeMasses`GetDimension[outLepton] =!= 1,
                           " int generationIndex2, ",
                           " "] <>
                 "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )\n" <>
                 "{\n" <>
                 IndentText[
                   FlexibleSUSY`FSModelName <> "_mass_eigenstates model_ = model;\n" <>
                   "EvaluationContext context{ model_ };\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[inLepton] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[outLepton] =!= 1,
                          " generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
                             ""]
                         ] <> "};\n\n" <>
                                 
                   "std::valarray<std::complex<double>> val {0.0, 0.0};\n\n" <>
                   
                   StringJoin @ Riffle[("val += " <> ToString @ # <> "::value(indices1, indices2, context);") & /@ 
                     Flatten[CXXEvaluatorsForLeptonPairAndDiagramsFromGraph[{inLepton, outLepton},#[[2]],#[[1]]] & /@ gTaggedDiagrams],
                                       "\n"] <> "\n\n" <>
                  "const auto leptonInMass = context.mass<Fe>(indices1);\n" <> 
                  "const double width = pow(leptonInMass,5)/(16.0*Pi) * (abs(val)*abs(val)).sum().real();\n" <>
                  "return width/(width + sm_width(generationIndex1, generationIndex2, model));"
                 ] <> "\n}";
    
    {prototype, definition}
  ];

CXXEvaluatorsForLeptonPairAndDiagramsFromGraph[{inLepton_, outLepton_},diagrams_,graph_] :=
  CXXEvaluatorsForLeptonPairAndDiagramFromGraph[inLepton,outLepton,#,graph] & /@ diagrams;
CXXEvaluatorsForLeptonPairAndDiagramFromGraph[inLepton_,outLepton_,diagram_,vertexCorrectionGraph] := 
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = CXXDiagrams`LorentzConjugate[diagram[[4,3]]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[CXXEvaluatorFS[inLepton,outLepton,photonEmitter,exchangeParticle]]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[CXXEvaluatorSF[inLepton,outLepton,photonEmitter,exchangeParticle]]];
    
    Return["(unknown diagram)"];
  ]

CXXEvaluatorFS[inLepton_,outLepton_,photonEmitter_,exchangeParticle_] :=
  "EDMVertexCorrectionFS<" <> CXXDiagrams`CXXNameOfField[inLepton] <> ", " <>
  CXXDiagrams`CXXNameOfField[outLepton] <> ", " <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"

CXXEvaluatorSF[inLepton_,outLepton_,photonEmitter_,exchangeParticle_] :=
  "EDMVertexCorrectionSF<" <> CXXDiagrams`CXXNameOfField[inLepton] <> ", " <>
  CXXDiagrams`CXXNameOfField[outLepton] <> ", " <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"

(*End[];*)
EndPackage[];
