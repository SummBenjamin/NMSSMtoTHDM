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

// File generated at Mon 20 Mar 2017 08:24:11

#include "NMSSMtower_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 6*
      traceconjTYuTpTYu + 6*tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*
      traceYuAdjYu + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*
      AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*
      AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(0.04*twoLoop*(Conj(MassB)*Sqr(g1)*(621*MassB*Sqr(g1) +
      5*(-8*(traceAdjYuTYu - 2*MassB*traceYuAdjYu) + 9*(2*MassB + MassWB)*Sqr(
      g2))) + 5*(3*Conj(MassWB)*Sqr(g2)*(3*(MassB + 2*MassWB)*Sqr(g1) + 55*
      MassWB*Sqr(g2)) - 160*(traceAdjYuTYu - 2*MassG*traceYuAdjYu)*Conj(MassG)*
      Sqr(g3) + 2*(15*Power(g2,4)*Tr22 + 7.745966692414834*g1*Tr31 - 15*
      tracemd2YdAdjYuYuAdjYd - 15*tracemq2AdjYdYdAdjYuYu - 15*
      tracemq2AdjYuYuAdjYdYd - 90*tracemq2AdjYuYuAdjYuYu - 15*
      tracemu2YuAdjYdYdAdjYu - 90*tracemu2YuAdjYuYuAdjYu - 15*
      traceYdAdjTYuTYuAdjYd - 15*traceYdAdjYuTYuAdjTYd - 15*mHd2*
      traceYdAdjYuYuAdjYd - 15*mHu2*traceYdAdjYuYuAdjYd - 15*
      traceYuAdjTYdTYdAdjYu - 90*traceYuAdjTYuTYuAdjYu - 15*
      traceYuAdjYdTYdAdjTYu - 90*traceYuAdjYuTYuAdjTYu - 90*mHu2*
      traceYuAdjYuYuAdjYu + 3*Tr2U111*Sqr(g1) + 4*traceconjTYuTpTYu*Sqr(g1) - 4
      *MassB*traceconjTYuTpYu*Sqr(g1) + 4*tracemq2AdjYuYu*Sqr(g1) + 4*
      tracemu2YuAdjYu*Sqr(g1) + 4*mHu2*traceYuAdjYu*Sqr(g1) + 80*
      traceconjTYuTpTYu*Sqr(g3) - 80*MassG*traceconjTYuTpYu*Sqr(g3) + 80*
      tracemq2AdjYuYu*Sqr(g3) + 80*tracemu2YuAdjYu*Sqr(g3) + 80*mHu2*
      traceYuAdjYu*Sqr(g3) - 30*(mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 5*Conj(TLambdax)*((3*traceAdjYdTYd + traceAdjYeTYe)*Lambdax +
      (3*traceYdAdjYd + traceYeAdjYe)*TLambdax + 2*Conj(Kappa)*(Lambdax*TKappa
      + Kappa*TLambdax)) - 5*Conj(Lambdax)*(3*traceconjTYdTpTYd*Lambdax +
      traceconjTYeTpTYe*Lambdax + 3*tracemd2YdAdjYd*Lambdax + traceme2YeAdjYe*
      Lambdax + traceml2AdjYeYe*Lambdax + 3*tracemq2AdjYdYd*Lambdax + 6*mHd2*
      traceYdAdjYd*Lambdax + 3*mHu2*traceYdAdjYd*Lambdax + 3*ms2*traceYdAdjYd*
      Lambdax + 2*mHd2*traceYeAdjYe*Lambdax + mHu2*traceYeAdjYe*Lambdax + ms2*
      traceYeAdjYe*Lambdax + 2*(mHd2 + mHu2 + 4*ms2)*AbsSqr(Kappa)*Lambdax + 12
      *AbsSqr(TLambdax)*Lambdax + 3*traceconjTYdTpYd*TLambdax +
      traceconjTYeTpYe*TLambdax + 2*Conj(TKappa)*(Lambdax*TKappa + Kappa*
      TLambdax))))));


   return beta_mHu2;
}

/**
 * Calculates the three-loop beta function of mHu2.
 *
 * @return three-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_mHu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
