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

// File generated at Mon 20 Mar 2017 08:24:23

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
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_ms2;

   beta_ms2 = Re(4*oneOver16PiSqr*(3*ms2*AbsSqr(Kappa) + (mHd2 + mHu2 +
      ms2)*AbsSqr(Lambdax) + AbsSqr(TKappa) + AbsSqr(TLambdax)));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;


   double beta_ms2;

   beta_ms2 = Re(-0.8*twoLoop*(120*ms2*Sqr(Conj(Kappa))*Sqr(Kappa) + 20*(
      mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) + Conj(TLambdax)*(
      Lambdax*(3*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*
      traceAdjYuTYu + 3*MassWB*Sqr(g2))) + (-3*Sqr(g1) + 5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu - 3*Sqr(g2)))*TLambdax) + Conj(Lambdax)*(15
      *traceconjTYdTpTYd*Lambdax + 5*traceconjTYeTpTYe*Lambdax + 15*
      traceconjTYuTpTYu*Lambdax + 15*tracemd2YdAdjYd*Lambdax + 5*
      traceme2YeAdjYe*Lambdax + 5*traceml2AdjYeYe*Lambdax + 15*tracemq2AdjYdYd*
      Lambdax + 15*tracemq2AdjYuYu*Lambdax + 15*tracemu2YuAdjYu*Lambdax + 30*
      mHd2*traceYdAdjYd*Lambdax + 15*mHu2*traceYdAdjYd*Lambdax + 15*ms2*
      traceYdAdjYd*Lambdax + 10*mHd2*traceYeAdjYe*Lambdax + 5*mHu2*traceYeAdjYe
      *Lambdax + 5*ms2*traceYeAdjYe*Lambdax + 15*mHd2*traceYuAdjYu*Lambdax + 30
      *mHu2*traceYuAdjYu*Lambdax + 15*ms2*traceYuAdjYu*Lambdax + 20*AbsSqr(
      TKappa)*Lambdax + 40*AbsSqr(TLambdax)*Lambdax - 3*mHd2*Lambdax*Sqr(g1) -
      3*mHu2*Lambdax*Sqr(g1) - 3*ms2*Lambdax*Sqr(g1) - 15*mHd2*Lambdax*Sqr(g2)
      - 15*mHu2*Lambdax*Sqr(g2) - 15*ms2*Lambdax*Sqr(g2) + 15*traceconjTYdTpYd*
      TLambdax + 5*traceconjTYeTpYe*TLambdax + 15*traceconjTYuTpYu*TLambdax +
      20*Conj(TKappa)*Kappa*TLambdax + 3*Conj(MassB)*Sqr(g1)*(-2*MassB*Lambdax
      + TLambdax) + 15*Conj(MassWB)*Sqr(g2)*(-2*MassWB*Lambdax + TLambdax)) +
      20*Conj(Kappa)*((mHd2 + mHu2 + 4*ms2)*AbsSqr(Lambdax)*Kappa + 4*AbsSqr(
      TKappa)*Kappa + Conj(TLambdax)*(Lambdax*TKappa + Kappa*TLambdax))));


   return beta_ms2;
}

/**
 * Calculates the three-loop beta function of ms2.
 *
 * @return three-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_ms2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

} // namespace flexiblesusy
