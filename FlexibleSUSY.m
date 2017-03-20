FSModelName = "NMSSMtower";
FSEigenstates = SARAH`EWSB;
FSDefaultSARAHModel = NMSSM;
OnlyLowEnergyFlexibleSUSY = True;
FlexibleEFTHiggs = True;

EXTPAR = {
    {0, MSUSY},
    {1, M1Input},
    {2, M2Input},
    {3, M3Input},
    {4, MuInput},
    {25, TanBeta},
    {61, LambdaInput},
    {62, KappaInput},
    {63, ALambdaInput},
    {64, AKappaInput}
};

FSExtraInputParameters = {
    {mq2Input, MSQ2IN, {3,3}}, (* 3x3 matrix *)
    {mu2Input, MSU2IN, {3,3}}, (* 3x3 matrix *)
    {md2Input, MSD2IN, {3,3}}, (* 3x3 matrix *)
    {ml2Input, MSL2IN, {3,3}}, (* 3x3 matrix *)
    {me2Input, MSE2IN, {3,3}}, (* 3x3 matrix *)
    {AuInput, AUIN, {3,3}}, (* 3x3 matrix *)
    {AdInput, ADIN, {3,3}}, (* 3x3 matrix *)
    {AeInput, AEIN, {3,3}}  (* 3x3 matrix *)
};

EWSBOutputParameters = { mHd2, mHu2, ms2 };

SUSYScale = MSUSY;

SUSYScaleFirstGuess = MSUSY;

SUSYScaleInput = {
    {vu, Sqrt[vu^2 + vd^2] Sin[ArcTan[TanBeta]]},
    {vd, Sqrt[vu^2 + vd^2] Cos[ArcTan[TanBeta]]},
    {MassB, M1Input},
    {MassWB, M2Input},
    {MassG, M3Input},
    {mq2, mq2Input},
    {mu2, mu2Input},
    {md2, md2Input},
    {ml2, ml2Input},
    {me2, me2Input},
    {T[Yu], AuInput Yu},
    {T[Yd], AdInput Yd},
    {T[Ye], AeInput Ye},
    {vS, Sqrt[2] MuInput/LambdaInput},
    {\[Kappa], KappaInput},
    {\[Lambda], LambdaInput},
    {T[\[Kappa]], AKappaInput KappaInput},
    {T[\[Lambda]], ALambdaInput LambdaInput}
};

InitialGuessAtSUSYScale = {
    {vu, LowEnergyConstant[vev] Sin[ArcTan[TanBeta]]},
    {vd, LowEnergyConstant[vev] Cos[ArcTan[TanBeta]]},
    {MassB, M1Input},
    {MassWB, M2Input},
    {MassG, M3Input},
    {mq2, mq2Input},
    {mu2, mu2Input},
    {md2, md2Input},
    {ml2, ml2Input},
    {me2, me2Input},
    {T[Yu], AuInput Yu},
    {T[Yd], AdInput Yd},
    {T[Ye], AeInput Ye},
    {vS, Sqrt[2] MuInput/LambdaInput},
    {\[Kappa], KappaInput},
    {\[Lambda], LambdaInput},
    {T[\[Kappa]], AKappaInput KappaInput},
    {T[\[Lambda]], ALambdaInput LambdaInput}
};

(* VEV is the SM-like VEV in the MSSM *)
SUSYScaleMatching = {
    {vu, VEV Sin[ArcTan[vu/vd]]},
    {vd, VEV Cos[ArcTan[vu/vd]]}
};

LowScale = LowEnergyConstant[MZ];

LowScaleFirstGuess = LowEnergyConstant[MZ];


UseHiggs2LoopNMSSM = True;
EffectiveMu = \[Lambda] vS / Sqrt[2];
EffectiveMASqr = (T[\[Lambda]] vS / Sqrt[2] + 0.5 \[Lambda] \[Kappa] vS^2) (vu^2 + vd^2) / (vu vd);

PotentialLSPParticles = { Chi, Sv, Su, Sd, Se, Cha, Glu };

DefaultPoleMassPrecision = HighPrecision;
HighPoleMassPrecision    = {hh, Ah, Hpm};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

ExtraSLHAOutputBlocks = {
   {FlexibleSUSYOutput,
           {{0, Hold[HighScale]},
            {1, Hold[SUSYScale]},
            {2, Hold[LowScale]} } },
   {ALPHA, {{ArcSin[Pole[ZH[2,2]]]}}},
   {HMIX , {{2, vu / vd},
            {3, Sqrt[vu^2 + vd^2]},
            {4, M[Ah[2]]^2},
            {102, vd},
            {103, vu} } },
   {Au,    {{1, 1, T[Yu][1,1] / Yu[1,1]},
            {2, 2, T[Yu][2,2] / Yu[2,2]},
            {3, 3, T[Yu][3,3] / Yu[3,3]} } },
   {Ad,    {{1, 1, T[Yd][1,1] / Yd[1,1]},
            {2, 2, T[Yd][2,2] / Yd[2,2]},
            {3, 3, T[Yd][3,3] / Yd[3,3]} } },
   {Ae,    {{1, 1, T[Ye][1,1] / Ye[1,1]},
            {2, 2, T[Ye][2,2] / Ye[2,2]},
            {3, 3, T[Ye][3,3] / Ye[3,3]} } },
   {MSOFT, {{1, MassB},
            {2, MassWB},
            {3, MassG},
            {21, mHd2},
            {22, mHu2},
            {31, SignedAbsSqrt[ml2[1,1]]},
            {32, SignedAbsSqrt[ml2[2,2]]},
            {33, SignedAbsSqrt[ml2[3,3]]},
            {34, SignedAbsSqrt[me2[1,1]]},
            {35, SignedAbsSqrt[me2[2,2]]},
            {36, SignedAbsSqrt[me2[3,3]]},
            {41, SignedAbsSqrt[mq2[1,1]]},
            {42, SignedAbsSqrt[mq2[2,2]]},
            {43, SignedAbsSqrt[mq2[3,3]]},
            {44, SignedAbsSqrt[mu2[1,1]]},
            {45, SignedAbsSqrt[mu2[2,2]]},
            {46, SignedAbsSqrt[mu2[3,3]]},
            {47, SignedAbsSqrt[md2[1,1]]},
            {48, SignedAbsSqrt[md2[2,2]]},
            {49, SignedAbsSqrt[md2[3,3]]} } },
   {NMSSMRUN,
           {{1, \[Lambda]},
            {2, \[Kappa]},
            {3, T[\[Lambda]] / \[Lambda]},
            {4, T[\[Kappa]] / \[Kappa]},
            {5, \[Lambda] vS / Sqrt[2]},
            {10, ms2} } }
};
