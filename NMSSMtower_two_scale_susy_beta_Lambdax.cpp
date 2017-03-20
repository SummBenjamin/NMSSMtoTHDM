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

// File generated at Mon 20 Mar 2017 08:22:47

#include "NMSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(3*traceYdAdjYd*Lambdax +
      traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax + 2*AbsSqr(Kappa)*Lambdax -
      0.6*Lambdax*Sqr(g1) - 3*Lambdax*Sqr(g2) + 4*Conj(Lambdax)*Sqr(Lambdax)))
      ;


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.02*twoLoop*Lambdax*(-207*Power(g1,4) - 375*Power(
      g2,4) + 450*traceYdAdjYdYdAdjYd + 300*traceYdAdjYuYuAdjYd + 150*
      traceYeAdjYeYeAdjYe + 450*traceYuAdjYuYuAdjYu + 600*AbsSqr(Kappa)*AbsSqr(
      Lambdax) + 20*traceYdAdjYd*Sqr(g1) - 60*traceYeAdjYe*Sqr(g1) - 40*
      traceYuAdjYu*Sqr(g1) - 90*Sqr(g1)*Sqr(g2) - 30*AbsSqr(Lambdax)*(-15*
      traceYdAdjYd - 5*traceYeAdjYe - 15*traceYuAdjYu + 2*Sqr(g1) + 10*Sqr(g2))
      - 800*traceYdAdjYd*Sqr(g3) - 800*traceYuAdjYu*Sqr(g3) + 400*Sqr(Conj(
      Kappa))*Sqr(Kappa) + 500*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the three-loop beta function of Lambdax.
 *
 * @return three-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_Lambdax_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
