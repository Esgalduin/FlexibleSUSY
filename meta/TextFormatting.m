$flexiblesusyCSrcChunkSize = 2^21;

BeginPackage["TextFormatting`"];

WrapLines::usage="breaks text lines";
WrapText::usage="WrapText[text, maxWidth, indentation] breaks text lines.  It tries to wrap a line at a blank or a special character that maximizes the line length within maxWidth characters.";
IndentText::usage="indents text by a given number of spaces";

Begin["`Private`"];

GetBestSplitPoint[line_String, maxWidth_:79] :=
    Module[{chars, lastSplitPoint, i, char, nextChar, numberOfQuotes = 0,
            splitChars = StringSplit["(){}[]*+,;<> ", ""]},
           If[StringLength[line] <= maxWidth, Return[StringLength[line]]];
           chars = StringSplit[line, ""];
           lastSplitPoint = StringLength[line];
           For[i = 1, i < Length[chars] && i <= maxWidth, i++,
               char = chars[[i]];
               If[char == "\"", numberOfQuotes++;];
               nextChar = chars[[i+1]];
               If[(MemberQ[splitChars, char] || MemberQ[splitChars, nextChar])
                  && EvenQ[numberOfQuotes],
                  lastSplitPoint = i;
                 ];
              ];
           Return[lastSplitPoint];
          ];

SplitLine[line_String, maxWidth_:79] :=
    Module[{bestSplitPoint, first, rest, result = {}},
           bestSplitPoint = GetBestSplitPoint[line, maxWidth];
           first = StringTake[line, bestSplitPoint];
           AppendTo[result, first];
           rest = StringDrop[line, bestSplitPoint];
           While[rest != "",
                 bestSplitPoint = GetBestSplitPoint[rest, maxWidth];
                 first = StringTake[rest, bestSplitPoint];
                 AppendTo[result, first];
                 rest = StringDrop[rest, bestSplitPoint];
                ];
           Return[result];
          ];

WrapLines[text_String, maxWidth_:79, offset_:"   "] :=
    Module[{result = "", lines, line, i, k, indent, splitLines},
           lines = StringSplit[text, "\n"];
           For[i = 1, i <= Length[lines], i++,
               line = lines[[i]];
               (* get intentation of line *)
               indent = StringReplace[line, x:(StartOfString ~~ Whitespace...) ~~ ___ -> x];
               splitLines = SplitLine[line, maxWidth - StringLength[indent]
                                      - StringLength[offset]];
               (* remove strings that consist of whitespace only *)
               splitLines = Cases[StringReplace[#,StartOfString ~~ Whitespace.. ~~ EndOfString -> EmptyString]& /@ splitLines, _String];
               For[k = 1, k <= Length[splitLines], k++,
                   (* strip whitespace *)
                   splitLines[[k]] = StringTrim[splitLines[[k]]];
                   result = result <> indent;
                   If[k > 1, result = result <> offset];
                   result = result <> splitLines[[k]] <> "\n";
                  ];
              ];
           Return[result];
          ];

IndentText[text_String, spaces_Integer:3] :=
    Module[{i, whiteSpace = ""},
           If[text == "", Return[text];];
           For[i = 0, i < spaces, i++, whiteSpace = whiteSpace <> " "];
           Return[StringReplace[whiteSpace <> text,
                                { x:("\n" ~~ ("\n"..)) :> x <> whiteSpace,
                                  "\n" ~~ EndOfString -> "\n",
                                  "\n" -> "\n" <> whiteSpace } ]];
          ];

WrapText[text_String, maxWidth_Integer:79, indentation_Integer:2] := Block[{
    maxLength = maxWidth,
    relSkip = indentation
  },
  StringJoin[WrapLine /@ StringSplit[text, "\n"]]
];

WrapLine[blank_String] := "\n" /;
  StringMatchQ[blank, RegularExpression["^[[:space:]]*$"]];

WrapLine[line_String] := Block[{
    nLeadingSpaces = NLeadingSpaces[line],
    absSkip,
    lst,
    firstBlankStr,
    otherBlankStr
  },
  absSkip = nLeadingSpaces + relSkip;
  lst = SplitLine[nLeadingSpaces, ProtectCTokens@StringTrim[line]];
  firstBlankStr = StringJoin@Table[" ", {nLeadingSpaces}];
  otherBlankStr = StringJoin@Table[" ", {absSkip}];
  If[lst === {}, {},
     {firstBlankStr, First[lst], "\n", {otherBlankStr, #, "\n"}& /@ Rest[lst]}]
];

NLeadingSpaces[line_String] := Module[{
    sp = StringPosition[line, RegularExpression["[^ ]"], 1]
  },
  If[sp === {}, 0, First@First[sp] - 1]
];

ProtectCTokens[line_String] :=
  DeleteCases[StringSplit[
      line,
      s:RegularExpression["\".*?(?<!\\\\)\"|##|<:|:>|<%|%>|(%:){1,2}|\\.\\.\\.|::|\\.\\*|[-+*/%&|^]=|(<<|>>)=?|[=!<>]=|&&|\\|\\||\\+\\+|--|->\\*?"]
      :> Hold[s]], ""];

SplitLine[fstLen_Integer, strs_List] :=
  StringJoin /@ Last@Reap[Fold[SplitString, {1, fstLen}, strs]];

SplitString[{lineN_, curLen_}, Hold[str_String]] := Module[{
    strLen = StringLength[str],
    nextLen, nextLineN = lineN
  },
  nextLen = curLen + strLen;
  If[nextLen > maxLength, nextLineN++; nextLen = absSkip + strLen];
  Sow[str, nextLineN];
  If[nextLen < maxLength, {nextLineN, nextLen}, {nextLineN + 1, absSkip}]
];

SplitString[{lineN_, curLen_}, str_String] := Block[{
    split,
    strLen, nextLen, nextLineN = lineN,
    remainder
  },
  split = Flatten@Last@Reap@SowStrings[curLen, str];
  If[split === {}, {lineN, curLen},
    nextLen = curLen + StringLength@First[split];
    If[nextLen > maxLength, nextLineN++];
    Sow[#, nextLineN++]& /@ split;
    If[Length[split] > 1, nextLen = absSkip + StringLength@Last[split]];
    If[nextLen < maxLength, {nextLineN - 1, nextLen}, {nextLineN, absSkip}]]
];

splitRegex = "[(){}\\[\\]<>,;:?~!=%^&|*+/-]";

SowStrings[_Integer, ""] := True;

SowStrings[curLen_Integer, str_String] :=
  Which[
    curLen + StringLength[str] <= maxLength, Sow[str],
    Head[StringReplace[str,
		       RegularExpression[
			   "^((?:.{0," <> ToString[maxLength - curLen - 1] <>
			   "}|.*?)(?:" <> splitRegex <>
			   "|[^[:space:]](?=[[:space:]])))[[:space:]]*(.*)"] :>
		       (Sow["$1"]; remainder = "$2"; True)]] =!= String,
      SowStrings[absSkip, remainder],
    True, Sow[str]];

End[];

EndPackage[];
