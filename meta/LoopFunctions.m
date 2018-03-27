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

BeginPackage["LoopFunctions`"];
EndPackage[];

(* loop functions *)
{A0, B0, B1, B00, B11, B22, B22tilde, C0, C1, C2, D0, D27, F, G, H};

LFFull::usage = "Returns explicit form of loop functions with full
 momentum dependence.  The convention depends on the value of
 $BPMZSign .

 The following functions are implemented:

   A0[m,Q]                     [arxiv:hep-ph/9606211, Eq. (B.5)]
   B0[p,m1,m2,Q]               [arxiv:hep-ph/9606211, Eq. (B.7)]
   B1[p,m1,m2,Q]               [arxiv:hep-ph/9606211, Eq. (B.5)]
   B00[p,m1,m2,Q]              [arxiv:hep-ph/9606211, Eq. (B.10)]
   B11[p,m1,m2,Q]              [arxiv:hep-ph/0709.1075, Eq. (4.6)]
   B22[p,m1,m2,Q]              [arxiv:hep-ph/9606211, Eq. (B.10)]
   B22tilde[p,m1,m2,Q]         [arxiv:hep-ph/9606211, Eq. (B.14)]
   F[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.11)]
   G[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.12)]
   H[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.13)]
   C0[0,0,m1,m2,m3,Q]          [arxiv:hep-ph/9606211, Eq. (C.19)]
   C1[0,0,m1,m2,m3,Q]
   C2[0,0,m1,m2,m3,Q]
   D0[0,0,0,m1,m2,m3,m4,Q]     [arxiv:hep-ph/9606211, Eq. (C.21)]
   D27[0,0,0,m1,m2,m3,m4,Q]    [arxiv:hep-ph/9606211, Eq. (C.22)]

B00[p,m1,m2,Q] is an alias for B22[p,m1,m2,Q].
";

LFZeroMomentum::usage = "Returns loop functions at zero momentum.  The
 convention depends on the value of $BPMZSign .";

LFDivergence::usage = "Returns divergent part of loop functions.  The
 convention depends on the value of $BPMZSign .";

LFScaleDependence::usage = "Returns the renormalization scale
 dependent part of loop functions.  The convention depends on the
 value of $BPMZSign .";

Delta::usage = "1/eps - \[Gamma]_E + Log[4 Pi]";

$BPMZSign::usage = "If set to 1, loop functions are returned in BPMZ
 convention (arxiv:hep-ph/9606211).  If set to -1, loop functions are
 returned in Denner convention (arxiv:hep-ph/0709.1075).  The default
 is -1.";

$BPMZSign = -1;

Begin["LoopFunctions`Private`"];

LFZeroMomentum[] := {
    A0[m_,mu_]                                               :> A0impl[m,mu],
    B0[0,m1_,m2_,mu_]             | B0[m1_,m2_,mu_]          :> B0zero[m1,m2,mu],
    B1[0,m1_,m2_,mu_]             | B1[m1_,m2_,mu_]          :> B1zero[m1,m2,mu],
    B00[0,m1_,m2_,mu_]            | B00[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B11[0,m1_,m2_,mu_]            | B11[m1_,m2_,mu_]         :> B11zero[m1,m2,mu],
    B22[0,m1_,m2_,mu_]            | B22[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B22tilde[0,m1_,m2_,mu_]       | B22tilde[m1_,m2_,mu_]    :> B22tildezero[m1,m2,mu],
    C0[0,0,m1_,m2_,m3_,mu_]       | C0[m1_,m2_,m3_,mu_]      :> C0zero[m1,m2,m3],
    C1[0,0,m1_,m2_,m3_,mu_]       | C1[m1_,m2_,m3_,mu_]      :> C1zero[m1,m2,m3],
    C2[0,0,m1_,m2_,m3_,mu_]       | C2[m1_,m2_,m3_,mu_]      :> C2zero[m1,m2,m3],
    D0[0,0,0,m1_,m2_,m3_,m4_,mu_] | D0[m1_,m2_,m3_,m4_,mu_]  :> D0zero[m1,m2,m3,m4],
    D27[0,0,0,m1_,m2_,m3_,m4_,mu_]| D27[m1_,m2_,m3_,m4_,mu_] :> D27zero[m1,m2,m3,m4],
    F[0,m1_,m2_,mu_]              | F[m1_,m2_,mu_]           :> Fzero[m1,m2,mu],
    G[0,m1_,m2_,mu_]              | G[m1_,m2_,mu_]           :> Gzero[m1,m2,mu],
    H[0,m1_,m2_,mu_]              | H[m1_,m2_,mu_]           :> Hzero[m1,m2,mu]
};

LFFull[] := {
    A0[m_,mu_]                          :> A0impl[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> B0impl[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> B1impl[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> B11impl[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> B22tildeimpl[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> C0impl[p1,p2,m1,m2,m3,mu],
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> C1impl[p1,p2,m1,m2,m3],
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> C2impl[p1,p2,m1,m2,m3],
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> D0impl[p1,p2,p3,m1,m2,m3,m4,mu],
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> D27impl[p1,p2,p3,m1,m2,m3,m4,mu],
    F[p_,m1_,m2_,mu_]                   :> Fimpl[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> Gimpl[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> Himpl[p,m1,m2,mu]
};

LFDivergence[] := {
    A0[m_,mu_]                          :> DivA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> DivB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> DivB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> DivB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> DivB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> DivF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> DivG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> DivH[p,m1,m2,mu]
};

LFScaleDependence[] := {
    A0[m_,mu_]                          :> LogA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> LogB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> LogB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> LogB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> LogB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> LogF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> LogG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> LogH[p,m1,m2,mu]
};

(********************* A0 *********************)

(* A0 [arxiv:hep-ph/9606211 Eq. (B.5)] *)
A0impl[m2_, mu2_] :=
    If[PossibleZeroQ[m2],
       0,
       m2 (Delta + 1 + Log[mu2/m2])
      ];

DivA0[m2_, _] := m2 Delta;

LogA0[m2_, mu2_] :=
    If[PossibleZeroQ[m2],
       0,
       m2 Log[mu2/m2]
      ];

(********************* B0 *********************)

(* B0 [arxiv:hep-ph/9606211 Eq. (B.6)] *)
B0impl[p2_, m12_, m22_, mu2_] :=
    Which[PossibleZeroQ[p2],
          B0zero[m12,m22,mu2],
          PossibleZeroQ[p2 - m22] && PossibleZeroQ[m12],
          Delta + Log[mu2/m22] + 2,
          PossibleZeroQ[p2 - m12] && PossibleZeroQ[m22],
          Delta + Log[mu2/m12] + 2,
          PossibleZeroQ[m12] && PossibleZeroQ[m22],
          Delta - Log[(-p2 - I eps)/mu] + 2,
          True,
          B0analytic[p2,m12,m22,mu2]
      ];

(* B0 for p^2 = 0 *)
B0zero[m12_, m22_, mu2_] :=
    Which[PossibleZeroQ[m12 - m22],
          Delta + Log[mu2/m22],
          PossibleZeroQ[m12],
          Delta + 1 + Log[mu2/m22],
          PossibleZeroQ[m22],
          Delta + 1 + Log[mu2/m12],
          True,
          Delta + 1 + (m12 Log[mu2/m12] - m22 Log[mu2/m22])/(m12 - m22)
   ];

(* B0 general (symmetric) form [arxiv:hep-ph/9606211 Eq. (B.7)] *)
B0analytic[p2_, m12_, m22_, mu2_] :=
    Module[{eps, fB, s, xp, xm},
           s = p2 - m22 + m12;
           xp = (s + Sqrt[s^2 - 4 p2 (m12 - I eps)]) / (2 p2);
           xm = (s - Sqrt[s^2 - 4 p2 (m12 - I eps)]) / (2 p2);
           fB[x_] := Log[1-x] - x Log[1 - 1/x] - 1;
           Limit[Delta - Log[p2/mu2] - fB[xp] - fB[xm],
                 eps -> 0, Direction -> -1,
                 Assumptions :> p2 > 0 && m12 >= 0 && m22 >= 0 && mu2 > 0]
          ];

(* B0 with explicit integration [arxiv:hep-ph/9606211 Eq. (B.6)] *)
B0integral[p2_, m12_, m22_, mu2_] :=
    Module[{eps},
           Limit[
               Delta - Integrate[Log[((1-x) m12 + x m22 - x (1-x) p2 - I eps)/mu2], {x,0,1}],
               eps -> 0]
          ];

DivB0[_, _, _, _] := Delta;

LogB0[p2_, _, _, mu2_] := Log[mu2/p2];

(********************* B1 *********************)

(* B1 [arxiv:hep-ph/9606211 Eq. (B.9)] *)
B1impl[p2_, m12_, m22_, mu2_] :=
    If[PossibleZeroQ[p2],
       B1zero[m12,m22,mu2],
       $BPMZSign 1/(2 p2) (A0impl[m22,mu2] - A0impl[m12,mu2]
                           + (p2 + m12 - m22) B0impl[p2,m12,m22,mu2])
      ];

(* B1 for p = 0 *)
B1zero[m12_, m22_, mu2_] :=
    Which[PossibleZeroQ[m12],
          $BPMZSign (1/4 + Delta/2 + Log[mu2/m22]/2),
          PossibleZeroQ[m22],
          $BPMZSign (3/4 + Delta/2 + Log[mu2/m12]/2),
          PossibleZeroQ[m12 - m22],
          $BPMZSign (Delta + Log[mu2/m22])/2,
          True,
          $BPMZSign 1/2 (Delta + 1 + Log[mu2/m22]
                         + (m12/(m12 - m22))^2 Log[m22/m12]
                         + 1/2 (m12 + m22)/(m12 - m22))
         ];

DivB1[_, _, _, _] := $BPMZSign Delta / 2;

LogB1[p2_, m12_, m22_, mu2_] :=
    If[PossibleZeroQ[p2],
       $BPMZSign Log[mu2/m22]/2,
       $BPMZSign 1/(2 p2) (LogA0[m22,mu2] - LogA0[m12,mu2]
                           + (p2 + m12 - m22) LogB0[p2,m12,m22,mu2])
      ];

(********************* B11 *********************)

(* B11 *)
B11impl[p2_, m12_, m22_, mu2_] :=
    If[PossibleZeroQ[p2],
       B11zero[m12,m22,mu2],
       1/(6 p2) (
           2 A0impl[m22,mu2]
           - 2 m12 B0impl[p2,m12,m22,mu2]
           + $BPMZSign 4 (p2 - m22 + m12) B1impl[p2,m12,m22,mu2]
           - m12 - m22 + p2/3)
      ];

(* B11 for p = 0 *)
(* DirectedInfinity[(m1^2 - m2^2)*(-2*m1^2 + 2*m2^2)] *)
B11zero[__] := Undefined;

DivB11[__] := Delta / 3;

LogB11[p2_, m12_, m22_, mu2_] :=
    If[PossibleZeroQ[p2],
       Undefined,
       1/(6 p2) (
           2 LogA0[m22,mu2]
           - 2 m12 LogB0[p2,m12,m22,mu2]
           + $BPMZSign 4 (p2 - m22 + m12) LogB1[p2,m12,m22,mu2])
      ];

(********************* B22 (= B00) *********************)

(* B22 [arxiv:hep-ph/9606211 Eq. (B.10)],
   identical to B00[p,m1,m2,mu] *)
B22impl[p2_, m12_, m22_, mu2_] :=
    If[PossibleZeroQ[p2],
       B22zero[m12,m22,mu2],
       1/6 (1/2 (A0impl[m12,mu2] + A0impl[m22,mu2])
            + (m12 + m22 - p2/2) B0impl[p2,m12,m22,mu2]
            + (m22 - m12)/(2 p2) (A0impl[m22,mu2] - A0impl[m12,mu2]
                                     - (m22 - m12) B0impl[p2,m12,m22,mu2])
            + m12 + m22 - p2/3)
      ];

(* B22 for p = 0 *)
B22zero[m12_, m22_, mu2_] :=
    Which[PossibleZeroQ[m12] && PossibleZeroQ[m22],
          0,
          PossibleZeroQ[m12],
          (m22*(5 + 3*Delta + 3*Log[mu2/m22]))/12,
          PossibleZeroQ[m22],
          (m12*(5 + 3*Delta + 3*Log[mu2/m12]))/12,
          PossibleZeroQ[m12 - m22],
          (m22*(1 + Delta + Log[mu2/m22]))/2,
          True,
          ((5 + 3*Delta)*(m12^2 - m22^2) + m12*(3*m12 + m22)*Log[mu2/m12] -
           m22*(m12 + 3*m22)*Log[mu2/m22])/(12*(m12 - m22))
      ];

DivB22[p2_, m12_, m22_, _] := Delta (3*m12 + 3*m22 - p2)/12;

LogB22[p2_, m12_, m22_, mu2_] :=
    Which[PossibleZeroQ[p2] && PossibleZeroQ[m12 - m22],
          m22 Log[mu2/m22]/2,
          PossibleZeroQ[p2],
          (m12*(3*m12 + m22)*Log[mu2/m12] -
           m22*(m12 + 3*m22)*Log[mu2/m22])/(12*(m12 - m22)),
          True,
          1/6 (1/2 (LogA0[m12,mu2] + LogA0[m22,mu2])
               + (m12 + m22 - p2/2) LogB0[p2,m12,m22,mu2]
               + (m22 - m12)/(2 p2) (LogA0[m22,mu2] - LogA0[m12,mu2]
                                        - (m22 - m12) LogB0[p2,m12,m22,mu2]))
         ];

(********************* C0 *********************)

C0impl[p12_, p22_, m12_, m22_, m32_, mu2_] :=
    If[PossibleZeroQ[p12] && PossibleZeroQ[p22],
       C0zero[m12,m22,m32],
       C0analytic[Sqrt[p12], Sqrt[p22], Sqrt[m12], Sqrt[m22], Sqrt[m32], Sqrt[mu2]]
      ];

(* C0 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.19)] *)
C0zero[m12_, m22_, m32_] := Which[
    PossibleZeroQ[m12 - m22] && PossibleZeroQ[m12 - m32],
    -1/(2*m32),
    PossibleZeroQ[m12 - m22],
    (-m22 + m32 - m32*Log[m32/m22])/(m22 - m32)^2,
    PossibleZeroQ[m12 - m32],
    (m22 - m32 - m22*Log[m22/m32])/(m22 - m32)^2,
    PossibleZeroQ[m22 - m32],
    (m12 - m32 + m12*Log[m32/m12])/(m12 - m32)^2,
    (* general case *)
    True,
    (+ m22 / (m12 - m22) Log[m22/m12]
     - m32 / (m12 - m32) Log[m32/m12]) / (m22 - m32)
];

(* C0 for complex momenta [arxiv:hep-ph/0709.1075, Eq. (4.26)] *)
(* Note: This implementation is divergent for real momenta.    *)
(* For real momenta a different implementation should be used. *)
(* Note: takes non-squared mass arguments! *)
C0analytic[p1_, p2_, m1_, m2_, m3_, mu_] :=
    Module[{p21 = (p2 - p1), p12 = (p1 - p2), pjk, pki, pij, mi, mj, mk,
            Dilogs, y0, xi, yi, alpha, alphai, eps, kappa, eta, result},
           (* Källén function, [arxiv:hep-ph/0709.1075, Eq. (4.28)] *)
           kappa[x_, y_, z_] := Sqrt[x^2 + y^2 + z^2 - 2 (x y + y z + z x)];
           (* [arxiv:hep-ph/0709.1075, Eq. (4.30)] *)
           eta[a_, b_] := Log[a b] - Log[a] - Log[b];

           pjk[i_?IntegerQ] :=
               Which[i == 0, p12, (* j = 1, k = 2 *)
                     i == 1, p2,  (* j = 2, k = 0 *)
                     i == 2, p1   (* j = 0, k = 1 *)];
           pki[i_?IntegerQ] :=
               Which[i == 0, p2,  (* j = 1, k = 2 *)
                     i == 1, p1,  (* j = 2, k = 0 *)
                     i == 2, p12  (* j = 0, k = 1 *)];
           pij[i_?IntegerQ] :=
               Which[i == 0, p1,  (* j = 1, k = 2 *)
                     i == 1, p12, (* j = 2, k = 0 *)
                     i == 2, p2   (* j = 0, k = 1 *)];
           mi[i_?IntegerQ] := Which[i == 0, m1,
                                    i == 1, m2,
                                    i == 2, m3];
           mj[i_?IntegerQ] := Which[i == 0, m2, (* j = 1, k = 2 *)
                                    i == 1, m3, (* j = 2, k = 0 *)
                                    i == 2, m1  (* j = 0, k = 1 *)];
           mk[i_?IntegerQ] := Which[i == 0, m3, (* j = 1, k = 2 *)
                                    i == 1, m1, (* j = 2, k = 0 *)
                                    i == 2, m2  (* j = 0, k = 1 *)];
           y0[i_] := (pjk[i]^2 (pjk[i]^2 - pki[i]^2 - pij[i]^2 + 2 mi[i]^2 - mj[i]^2 - mk[i]^2)
                      - (pki[i]^2 - pij[i]^2) (mj[i]^2 - mk[i]^2)
                      + alpha (pjk[i]^2 - mj[i]^2 + mk[i]^2)) / (2 alpha pjk[i]^2);
           xi[i_, s_] := (pjk[i]^2 - mj[i]^2 - mk[i]^2 + s alphai[i]) / (2 pjk[i]^2);
           yi[i_, s_] := y0[i] - xi[i,s];
           alpha = kappa[p1^2, p21^2, p2^2];
           alphai[i_] := kappa[pjk[i]^2, mj[i]^2, mk[i]^2] (1 + I eps pjk[i]^2);
           Dilogs[i_, s_] := (
               PolyLog[2,(y0[i] - 1)/yi[i,s]] - PolyLog[2,y0[i]/yi[i,s]]
               + eta[1 - xi[i,s], 1/yi[i,s]] Log[(y0[i] - 1)/yi[i,s]]
               - eta[-xi[i,s], 1/yi[i,s]] Log[y0[i]/yi[i,s]]);

           result = Sum[Dilogs[i,+1] + Dilogs[i,-1]
                        - (eta[-xi[i,+1],-xi[i,-1]]
                           - eta[yi[i,+1],yi[i,-1]]
                           - 2 Pi I UnitStep[-Re[pjk[i]^2]] UnitStep[-Im[yi[i,+1] yi[i,-1]]]
                          ) Log[(1 - y0[i])/(-y0[i])],
                        {i,0,2}] / alpha;

           Limit[result, eps -> 0, Direction -> -1,
                 Assumptions :> Element[p1, Complexes] || \
                                Element[p2, Complexes]]
          ];

(********************* C1 *********************)

C1impl[p12_, p22_, m12_, m22_, m32_] :=
    If[PossibleZeroQ[p12] && PossibleZeroQ[p22],
       C1zero[m12,m22,m32],
       NotImplemented
      ];

(* C1 for p = 0 *)
C1zero[m12_, m22_, m32_] :=
    Module[{t1 = m22/m12, t2 = m32/m12},
           Which[
               PossibleZeroQ[m12 - m22] && PossibleZeroQ[m12 - m32],
               -1/(6 m12),
               PossibleZeroQ[m12 - m22],
               -(m22^2 - 4*m22*m32 + 3*m32^2 - 2*m32^2*Log[m32/m22])/(4*(m22 - m32)^3),
               PossibleZeroQ[m22 - m32],
               (3*m12^2 - 4*m12*m32 + m32^2 + 2*m12^2*Log[m32/m12])/(4*(m12 - m32)^3),
               PossibleZeroQ[m12 - m32],
               (-m22^2 + m32^2 + 2*m22*m32*Log[m22/m32])/(2*(m22 - m32)^3),
               (* general case *)
               True,
               -(t1 / (2 (t1 - 1) (t1 - t2))
                 - t1 (t1 - 2 t2 + t1 t2) Log[t1] / (2 (t1 - 1)^2 (t1 - t2)^2)
                 + (t2^2 - 2 t1 t2^2 + t1^2 t2^2) Log[t2] / (2 (t1 - 1)^2 (t1 - t2)^2 (t2 - 1))
                )/m12
                ]
          ];

(********************* C2 *********************)

C2impl[p12_, p22_, m12_, m22_, m32_] :=
    If[PossibleZeroQ[p12] && PossibleZeroQ[p22],
       C2zero[m12,m22,m32],
       NotImplemented
      ];

(* C2 for p = 0 *)
C2zero[m12_, m22_, m32_] :=
    Module[{t1 = m22/m12, t2 = m32/m12},
           Which[
               PossibleZeroQ[m12 - m22] && PossibleZeroQ[m12 - m32],
               -1/(6 m12),
               PossibleZeroQ[m12 - m22],
               (-m22^2 + m32^2 + 2*m22*m32*Log[m22/m32])/(2*(m22 - m32)^3),
               PossibleZeroQ[m22 - m32],
               (3*m12^2 - 4*m12*m32 + m32^2 + 2*m12^2*Log[m32/m12])/(4*(m12 - m32)^3),
               PossibleZeroQ[m12 - m32],
               (3*m22^2 - 4*m22*m32 + m32^2 - 2*m22^2*Log[m22/m32])/(4*(m22 - m32)^3),
               (* general case *)
               True,
               -(- t2 / (2 (t1 - t2) (t2 - 1))
                 + Log[t1] / (2 (t1 - 1) (t2 - 1)^2)
                 + (2 t1 t2 - 2 t1^2 t2 - t2^2 + t1^2 t2^2) Log[t1/t2] / (2 (t1 - 1) (t1 - t2)^2 (t2 - 1)^2)
                )/m12
                ]
          ];

(********************* D0 *********************)

(* D0 *)
D0impl[p12_, p22_, p32_, m12_, m22_, m32_, m42_, mu2_] :=
    If[PossibleZeroQ[p12] && PossibleZeroQ[p22] && PossibleZeroQ[p32],
       D0zero[m12, m22, m32, m42],
       NotImplemented
      ];

(* D0 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.21)] *)
D0zero[m12_, m22_, m32_, m42_] := Which[
   PossibleZeroQ[m12 - m22] && PossibleZeroQ[m12 - m32] && 
    PossibleZeroQ[m12 - m42],
   1/(6 m12^2),
   PossibleZeroQ[m12 - m22] && PossibleZeroQ[m22 - m32],
   (8*(m22 - m42)2*(m22 + m42) + (-23*m22^2*m42 - 10*m22*m42^2 + 
         m42^3)*Log[m22/m42] + (-7*m22^2*m42 - 26*m22*m42^2 + m42^3)*
       Log[m42/m22])/(16*m22*(m22 - m42)2^2),
   PossibleZeroQ[m12 - m22] && PossibleZeroQ[m22 - m42],
   (8*(m22 - m32)2*(m22 + m32) + (7*m22^2*m32 + 26*m22*m32^2 - 
         m32^3)*Log[m22/m32] + (23*m22^2*m32 + 10*m22*m32^2 - m32^3)*
       Log[m32/m22])/(16*m22*(m22 - m32)2^2),
   PossibleZeroQ[m12 - m22] && PossibleZeroQ[m32 - m42],
   (-8*(m22 - m42)2 + (m22^2 - 8*m22*m42 - 5*m42^2)*
       Log[m22/m42] - (3*m22^2 + 8*m22*m42 + m42^2)*
       Log[m42/m22])/(4*(m22 - m42)2^2),
   PossibleZeroQ[m12 - m22],
   (-((m2 - m3)*(m2 + m3)*(m2 - m4)*(m3 - m4)*(m2 + m4)*(m3 + 
           m4)) - (m22^2 - m32*m42)*(m32*Log[m32/m22] + 
         m42*Log[m22/m42]) + 
      m32*m42*(-2*m22 + m32 + m42)*
       Log[m42/m32])/((m22 - m32)2*(m22 - m42)2*(m32 - m42)),
   True,
   (C0zero[m12, m32, m42] - C0zero[m22, m32, m42])/(m12 - m22)];

(********************* D27 *********************)

(* D27 *)
D27impl[p12_, p22_, p32_, m12_, m22_, m32_, m42_, mu2_] :=
    If[PossibleZeroQ[p12] && PossibleZeroQ[p22] && PossibleZeroQ[p32],
       D27zero[m12, m22, m32, m42],
       NotImplemented
      ];

(* D27 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.22)] *)
D27zero[m12_, m22_, m32_, m42_] :=
   (m12 C0zero[m12, m32, m42] - m22 C0zero[m22, m32, m42])/(4 (m12 - m22));

(********************* F *********************)

(* F [arxiv:hep-ph/9606211 Eq. (B.11)] *)
Fimpl[p2_, m12_, m22_, mu2_] :=
    A0impl[m12,mu2] - 2 A0impl[m22,mu2] - (2 p2 + 2 m12 - m22) B0impl[p2,m12,m22,mu2];

Fzero[m12_, m22_, mu2_] := Fimpl[0,m12,m22,mu2];

DivF[p2_, m12_, m22_, mu2_] := Delta (-m12 - m22 - 2*p2);

LogF[p2_, m12_, m22_, mu2_] :=
    LogA0[m12,mu2] - 2 LogA0[m22,mu2] - (2 p2 + 2 m12 - m22) LogB0[p2,m12,m22,mu2];

(********************* G *********************)

(* G [arxiv:hep-ph/9606211 Eq. (B.12)] *)
Gimpl[p2_, m12_, m22_, mu2_] :=
    (p2 - m12 - m22) B0impl[p2,m12,m22,mu2] - A0impl[m12,mu2] - A0impl[m22,mu2];

Gzero[m12_, m22_, mu2_] := Gimpl[0,m12,m22,mu2];

DivG[p2_, m12_, m22_, mu2_] := Delta (-2*m12 - 2*m22 + p2);

LogG[p2_, m12_, m22_, mu2_] :=
    (p2 - m12 - m22) LogB0[p2,m12,m22,mu2] - LogA0[m12,mu2] - LogA0[m22,mu2];

(********************* H *********************)

(* H [arxiv:hep-ph/9606211 Eq. (B.13)] *)
Himpl[p2_, m12_, m22_, mu2_] := 4 B22impl[p2,m12,m22,mu2] + Gimpl[p2,m12,m22,mu2];

Hzero[m12_, m22_, mu2_] := Himpl[0,m12,m22,mu2];

DivH[p2_, m12_, m22_, mu2_] := Delta (-m12 - m22 + (2*p2)/3);

LogH[p2_, m12_, m22_, mu2_] := 4 LogB22[p2,m12,m22,mu2] + LogG[p2,m12,m22,mu2];

(********************* ~B22 *********************)

(* ~B22 [arxiv:hep-ph/9606211 Eq. (B.14)] *)
B22tildeimpl[p2_, m12_, m22_, mu2_] :=
    B22impl[p2,m12,m22,mu2] - A0impl[m12,mu2]/4 - A0impl[m22,mu2]/4;

B22tildezero[m12_, m22_, mu2_] := B22tildeimpl[0,m12,m22,mu2];

DivB22tilde[p2_, _, _, _] := Delta (-p2/12);

LogB22tilde[p2_, m12_, m22_, mu2_] :=
    LogB22[p2,m12,m22,mu2] - LogA0[m12,mu2]/4 - LogA0[m22,mu2]/4;

End[];
