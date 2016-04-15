
BeginPackage["Utils`"];

ApplyAndConcatenate::usage = "Applies a function to a list and
concatenates the resulting list.";

InputFormOfNonStrings::usage = "apply InputForm to non-strings";

MaxRelDiff::usage="Returns maximum relative difference between numbers
 in a given list.  The second argument defines the numeric underflow.

In[]:= MaxRelDiff[{1, 1.1, 0.9}]
Out[]= 0.181818

In[]:= MaxRelDiff[{0, 10^(-100)}]
Out[]= 0

In[]:= MaxRelDiff[{0, 10^(-100)}, 10^(-101)]
Out[]= 1

In[]:= MaxRelDiff[{1, -1}]
Out[]= 2
";

StringJoinWithSeparator::usage = "Joins a list of strings with a given separator string";

Zip::usage = "Combines two lists to a list of touples.
Example:

In[]:= Zip[{a,b,c},{d,e,f}]
Out[]= {{a, d}, {b, e}, {c, f}}
";

StringZip::usage = "Combines lists of strings to a list of concated
strings.  Example:

In[]:= StringZip[{\"a\",\"b\"},{\"d\",\"e\"}]
Out[]= {\"ad\", \"be\"}}
";

StringZipWithSeparator::usage = "Combines lists of strings to a list
of concated strings using a separator.  Example:

In[]:= StringZipWithSeparator[{\"a\",\"b\"},{\"d\",\"e\"}, \"_\"]
Out[]= {\"a_d\", \"b_e\"}}
";

SplitList::usage = "split list into list of sub-lists of given maximum
 size";

ForceJoin::usage = "Joins the given arguments if they are lists.";

FSGetOption::usage = "Returns the value of an option from a list of
options.

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSGetOption[opts, o1]
Out[]= x

In[]:= FSGetOption[opts, o2]
Out[]= y

In[]:= FSGetOption[opts, o3]
Error: option o3 not found
";

FSSetOption::usage = "Returns given list of options, where the
occurrence of the given rule is replaced (if it exists) or added (if
it does not exist).

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSSetOption[opts, o1 -> z]
Out[]= {o2 -> y, o1 -> z}

In[]:= FSSetOption[opts, o3 -> z]
Out[]= {o1 -> x, o2 -> y, o3 -> z}
";

FSImportString::usage = "Returns the content of a file in form of a string.  If the file does not exist, \"unknown\" is returned.";

Begin["`Private`"];

ApplyAndConcatenate[Func_, l_List] :=
    Module[{result = ""},
           (result = result <> Evaluate[Func[#]])& /@ l;
           result
          ];

ApplyAndConcatenate[Func_, l_] := Evaluate[Func[l]];

SetAttributes[ApplyAndConcatenate, HoldFirst];

StringJoinWithSeparator[list_List, separator_String, transformer_:ToString] :=
    StringJoin[Riffle[transformer /@ list, separator]];

Zip[list1_List, list2_List] :=
    MapThread[List, {list1, list2}];

StringZip[lists___List] :=
    MapThread[StringJoin, {lists}];

StringZipWithSeparator[lists___List, separator_String] :=
    MapThread[StringJoinWithSeparator[{##},separator]&, {lists}];

SplitList[lst_List, 0] := {lst};

SplitList[lst_List, size_Integer] :=
    Module[{result = {}, list = lst, drops},
           While[list =!= {},
                 drops = Min[size,Length[list]];
                 AppendTo[result, Take[list, drops]];
                 list = Drop[list, drops];
                ];
           result
          ];

FSGetOption[opts_List, opt_] :=
    Module[{values},
           values = Cases[opts, (Rule[opt, value_] | RuleDelayed[opt, value_]) :> value];
           Switch[Length[values],
                  0, Print["Error: option ", opt, " not found"];
                     Null,
                  1, values[[1]],
                  _, Print["Warning: option ", opt, " not unique"];
                     values[[1]]
                 ]
          ];

FSSetOption[opts_List, rule:Rule[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSSetOption[opts_List, rule:RuleDelayed[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSImportString[fileName_String] :=
    Module[{str = Import[fileName, "String"]},
           If[str =!= $Failed,
              str,
              "unknown"
             ]
          ];

ForceJoin[elem___] :=
    Join[Sequence @@ Select[{elem}, (Head[#] === List)&]];

InputFormOfNonStrings[a_String] := a;
InputFormOfNonStrings[a_] := InputForm[a];

MaxRelDiff[{}, _] := 0;

MaxRelDiff[{a_, b_}, underflow_:10^(-16)] :=
    If[Abs[a] < underflow && Abs[b] < underflow,
       0,
       Abs[(a-b)/Max[Abs[{a,b}]]]
      ];

MaxRelDiff[numbers_List, underflow_:10^(-16)] :=
    Max[MaxRelDiff[#,underflow]& /@ Tuples[numbers, 2]];

End[];

EndPackage[];
