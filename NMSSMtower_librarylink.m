libNMSSMtower = FileNameJoin[{Directory[], "models", "NMSSMtower", "NMSSMtower_librarylink.so"}];

FSNMSSMtowerGetSettings = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerGetSettings", LinkObject, LinkObject];
FSNMSSMtowerGetSMInputParameters = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerGetSMInputParameters", LinkObject, LinkObject];
FSNMSSMtowerGetInputParameters = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerGetInputParameters", LinkObject, LinkObject];

FSNMSSMtowerOpenHandleLib = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerOpenHandle", {{Real,1}}, Integer];
FSNMSSMtowerCloseHandle = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerCloseHandle", {Integer}, Void];

FSNMSSMtowerSetLib = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerSet", {Integer, {Real,1}}, Void];

FSNMSSMtowerCalculateSpectrum = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerCalculateSpectrum", LinkObject, LinkObject];
FSNMSSMtowerCalculateObservables = LibraryFunctionLoad[libNMSSMtower, "FSNMSSMtowerCalculateObservables", LinkObject, LinkObject];

FSNMSSMtowerCalculateSpectrum::error = "`1`";
FSNMSSMtowerCalculateSpectrum::warning = "`1`";

FSNMSSMtowerCalculateObservables::error = "`1`";
FSNMSSMtowerCalculateObservables::warning = "`1`";

FSNMSSMtower::info = "`1`";
FSNMSSMtower::nonum = "Error: `1` is not a numeric input value!";
FSNMSSMtowerMessage[s_] := Message[FSNMSSMtower::info, s];

FSNMSSMtowerCheckIsNumeric[a_?NumericQ] := a;
FSNMSSMtowerCheckIsNumeric[a_] := (Message[FSNMSSMtower::nonum, a]; Abort[]);

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
      poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
      ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
      betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
      thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
      higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
      higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
      higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
      higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
      forceOutput -> 0,                  (* FlexibleSUSY[12] *)
      topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
      betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
      forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
      poleMassScale -> 0,                (* FlexibleSUSY[17] *)
      eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
      eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
      eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
      eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
      eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
      calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
      parameterOutputScale -> 0          (* MODSEL[12] *)
};

fsDefaultSMParameters = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.16637*^-5,      (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.777,          (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

fsNMSSMtowerDefaultInputParameters = {
   MSUSY -> 0,
   M1Input -> 0,
   M2Input -> 0,
   M3Input -> 0,
   MuInput -> 0,
   TanBeta -> 0,
   LambdaInput -> 0,
   KappaInput -> 0,
   ALambdaInput -> 0,
   AKappaInput -> 0,
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AuInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AdInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   AeInput -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
};

Options[FSNMSSMtowerOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsNMSSMtowerDefaultInputParameters
};

FSNMSSMtowerOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSNMSSMtowerOpenHandle[a, Sequence @@ s, r];

FSNMSSMtowerOpenHandle[OptionsPattern[]] :=
    FSNMSSMtowerOpenHandleLib[
        FSNMSSMtowerCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* NMSSMtower input parameters *)
            ,
            OptionValue[MSUSY],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[M3Input],
            OptionValue[MuInput],
            OptionValue[TanBeta],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[ALambdaInput],
            OptionValue[AKappaInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]]
        }
];

Options[FSNMSSMtowerSet] = Options[FSNMSSMtowerOpenHandle];

FSNMSSMtowerSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSNMSSMtowerSet[handle, a, Sequence @@ s, r];

FSNMSSMtowerSet[handle_Integer, p:OptionsPattern[]] :=
    FSNMSSMtowerSetLib[
        handle,
        ReleaseHold[Hold[{
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* NMSSMtower input parameters *)
            ,
            OptionValue[MSUSY],
            OptionValue[M1Input],
            OptionValue[M2Input],
            OptionValue[M3Input],
            OptionValue[MuInput],
            OptionValue[TanBeta],
            OptionValue[LambdaInput],
            OptionValue[KappaInput],
            OptionValue[ALambdaInput],
            OptionValue[AKappaInput],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[AuInput][[1,1]],
            OptionValue[AuInput][[1,2]],
            OptionValue[AuInput][[1,3]],
            OptionValue[AuInput][[2,1]],
            OptionValue[AuInput][[2,2]],
            OptionValue[AuInput][[2,3]],
            OptionValue[AuInput][[3,1]],
            OptionValue[AuInput][[3,2]],
            OptionValue[AuInput][[3,3]],
            OptionValue[AdInput][[1,1]],
            OptionValue[AdInput][[1,2]],
            OptionValue[AdInput][[1,3]],
            OptionValue[AdInput][[2,1]],
            OptionValue[AdInput][[2,2]],
            OptionValue[AdInput][[2,3]],
            OptionValue[AdInput][[3,1]],
            OptionValue[AdInput][[3,2]],
            OptionValue[AdInput][[3,3]],
            OptionValue[AeInput][[1,1]],
            OptionValue[AeInput][[1,2]],
            OptionValue[AeInput][[1,3]],
            OptionValue[AeInput][[2,1]],
            OptionValue[AeInput][[2,2]],
            OptionValue[AeInput][[2,3]],
            OptionValue[AeInput][[3,1]],
            OptionValue[AeInput][[3,2]],
            OptionValue[AeInput][[3,3]]
        }] /. HoldPattern[OptionValue[param_]] :> FSNMSSMtowerCheckIsNumeric[param] /.
        { p } /.
        FSNMSSMtowerGetSettings[handle] /.
        FSNMSSMtowerGetSMInputParameters[handle] /.
        FSNMSSMtowerGetInputParameters[handle]]];
