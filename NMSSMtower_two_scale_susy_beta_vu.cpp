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

// File generated at Mon 20 Mar 2017 08:22:51

#include "NMSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(-10*AbsSqr(Lambdax) + 3*(Sqr(g1) +
      5*(-2*traceYuAdjYu + Sqr(g2)))));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(207*Power(g1,4) + 10*Sqr(g1)*(34*
      traceYuAdjYu + 9*Sqr(g2)) + 20*AbsSqr(Lambdax)*(-30*traceYdAdjYd - 10*
      traceYeAdjYe - 20*AbsSqr(Kappa) + 3*Sqr(g1) + 15*Sqr(g2)) + 25*(11*Power(
      g2,4) - 24*traceYdAdjYuYuAdjYd - 72*traceYuAdjYuYuAdjYu + 36*traceYuAdjYu
      *Sqr(g2) + 128*traceYuAdjYu*Sqr(g3)) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax
      )));


   return beta_vu;
}

/**
 * Calculates the three-loop beta function of vu.
 *
 * @return three-loop beta function
 */
double NMSSMtower_susy_parameters::calc_beta_vu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

} // namespace flexiblesusy
