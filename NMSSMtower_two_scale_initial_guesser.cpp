// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Mon 20 Mar 2017 08:24:43

#include "NMSSMtower_two_scale_initial_guesser.hpp"
#include "NMSSMtower_two_scale_model.hpp"
#include "NMSSMtower_standard_model_two_scale_matching.hpp"
#include "standard_model_two_scale_model.hpp"
#include "lowe.h"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <cassert>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define SMPARAMETER(p) eft->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

NMSSMtower_standard_model_initial_guesser<Two_scale>::NMSSMtower_standard_model_initial_guesser(
   NMSSMtower<Two_scale>* model_,
   standard_model::StandardModel<Two_scale>* eft_,
   const softsusy::QedQcd& qedqcd_,
   const standard_model::Standard_model_low_scale_constraint<Two_scale>& low_constraint_,
   const NMSSMtower_susy_scale_constraint<Two_scale>& susy_constraint_
)
   : Initial_guesser<Two_scale>()
   , model(model_)
   , eft(eft_)
   , qedqcd(qedqcd_)
   , mu_guess(0.)
   , mc_guess(0.)
   , mt_guess(0.)
   , md_guess(0.)
   , ms_guess(0.)
   , mb_guess(0.)
   , me_guess(0.)
   , mm_guess(0.)
   , mtau_guess(0.)
   , running_precision(1.0e-3)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
{
   assert(model && "NMSSMtower_initial_guesser: Error: pointer to model"
          " NMSSMtower<Two_scale> must not be zero");
}

NMSSMtower_standard_model_initial_guesser<Two_scale>::~NMSSMtower_standard_model_initial_guesser()
{
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::guess()
{
   guess_eft_parameters();
   guess_model_parameters();
}

/**
 * Guesses the effective parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) are ignored
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::guess_eft_parameters()
{
   using namespace softsusy;

   softsusy::QedQcd leAtMt(qedqcd);
   const double MZ = Electroweak_constants::MZ;
   const double MW = Electroweak_constants::MW;
   const double sinThetaW2 = 1.0 - Sqr(MW / MZ);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(mUp);
   mc_guess = leAtMt.displayMass(mCharm);
   mt_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(mTop) - 30.0 :
      leAtMt.displayPoleMt();
   md_guess = leAtMt.displayMass(mDown);
   ms_guess = leAtMt.displayMass(mStrange);
   mb_guess = leAtMt.displayMass(mBottom);
   me_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(mElectron) :
      leAtMt.displayPoleMel();
   mm_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(mMuon) :
      leAtMt.displayPoleMmuon();
   mtau_guess = leAtMt.displayMass(mTau);

   // guess gauge couplings at mt
   const DoubleVector alpha_sm(leAtMt.getGaugeMu(mtpole, sinThetaW2));

   eft->set_g1(sqrt(4.0 * M_PI * alpha_sm(1)));
   eft->set_g2(sqrt(4.0 * M_PI * alpha_sm(2)));
   eft->set_g3(sqrt(4.0 * M_PI * alpha_sm(3)));
   eft->set_scale(mtpole);

   eft->set_v(246.22);
   eft->set_Yu(ZEROMATRIX(3,3));
   eft->set_Yd(ZEROMATRIX(3,3));
   eft->set_Ye(ZEROMATRIX(3,3));

   eft->set_Yu(0, 0, -Sqrt(2.)* mu_guess/ eft->get_v());
   eft->set_Yu(1, 1, -Sqrt(2.)* mc_guess/ eft->get_v());
   eft->set_Yu(2, 2, -Sqrt(2.)* mt_guess/ eft->get_v());

   eft->set_Yd(0, 0, Sqrt(2.)* md_guess/ eft->get_v());
   eft->set_Yd(1, 1, Sqrt(2.)* ms_guess/ eft->get_v());
   eft->set_Yd(2, 2, Sqrt(2.)* mb_guess/ eft->get_v());

   eft->set_Ye(0, 0, Sqrt(2.)* me_guess/ eft->get_v());
   eft->set_Ye(1, 1, Sqrt(2.)* mm_guess/ eft->get_v());
   eft->set_Ye(2, 2, Sqrt(2.)* mtau_guess/ eft->get_v());

   eft->set_Lambdax(0.12604);
   eft->solve_ewsb_tree_level();

}

void NMSSMtower_standard_model_initial_guesser<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::calculate_Yu_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar(ZEROMATRIXCOMPLEX(3,3));
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;


}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::calculate_Yd_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar(ZEROMATRIXCOMPLEX(3,3));
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;


}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::calculate_Ye_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar(ZEROMATRIXCOMPLEX(3,3));
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;


}

/**
 * Guesses the full model parameters.  At first it runs to the
 * guess of the SUSY-scale (SUSYScaleFirstGuess) and imposes the
 * SUSY scale initial guess and the
 * SUSY-scale constraint (SUSYScaleInput).  Afterwards, it solves the
 * EWSB conditions at the tree-level.
 * Finally the DR-bar mass spectrum is calculated.
 */
void NMSSMtower_standard_model_initial_guesser<Two_scale>::guess_model_parameters()
{
   const double susy_scale_guess = susy_constraint.get_initial_scale_guess();

   model->set_scale(susy_scale_guess);

   // apply susy-scale first guess
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M3Input = INPUTPARAMETER(M3Input);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto AuInput = INPUTPARAMETER(AuInput);
   const auto AdInput = INPUTPARAMETER(AdInput);
   const auto AeInput = INPUTPARAMETER(AeInput);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto AKappaInput = INPUTPARAMETER(AKappaInput);
   const auto ALambdaInput = INPUTPARAMETER(ALambdaInput);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_vu(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vd(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_MassB(Re(M1Input));
   MODEL->set_MassWB(Re(M2Input));
   MODEL->set_MassG(Re(M3Input));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_TYu((AuInput*Yu).real());
   MODEL->set_TYd((AdInput*Yd).real());
   MODEL->set_TYe((AeInput*Ye).real());
   MODEL->set_vS(Re((1.4142135623730951*MuInput)/LambdaInput));
   MODEL->set_Kappa(Re(KappaInput));
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_TKappa(Re(AKappaInput*KappaInput));
   MODEL->set_TLambdax(Re(ALambdaInput*LambdaInput));


   eft->run_to(susy_scale_guess, running_precision);
   eft->calculate_DRbar_masses();

   //get gauge and Yukawa couplings from effective theory
   NMSSMtower_standard_model_Matching<Two_scale> matching;
   matching.set_models(eft, model);
   matching.set_constraint(&susy_constraint);
   matching.match_low_to_high_scale_model_tree_level();

   model->run_to(susy_scale_guess, running_precision);

   // apply susy-scale constraint
   susy_constraint.set_model(model);
   susy_constraint.apply();

   // apply EWSB constraint
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

} // namespace flexiblesusy
