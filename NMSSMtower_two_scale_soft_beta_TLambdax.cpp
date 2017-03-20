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

// File generated at Mon 20 Mar 2017 08:23:14

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
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*(6*traceAdjYdTYd*Lambdax + 2*
      traceAdjYeTYe*Lambdax + 6*traceAdjYuTYu*Lambdax + 1.2*MassB*Lambdax*Sqr(
      g1) + 6*MassWB*Lambdax*Sqr(g2) + (3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + 12*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2))*TLambdax + 2
      *Conj(Kappa)*(2*Lambdax*TKappa + Kappa*TLambdax)));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_TLambdax;

   beta_TLambdax = Re(twoLoop*(-16.56*Power(g1,4)*MassB*Lambdax - 30*
      Power(g2,4)*MassWB*Lambdax - 36*traceYdAdjYdTYdAdjYd*Lambdax - 12*
      traceYdAdjYuTYuAdjYd*Lambdax - 12*traceYeAdjYeTYeAdjYe*Lambdax - 12*
      traceYuAdjYdTYdAdjYu*Lambdax - 36*traceYuAdjYuTYuAdjYu*Lambdax - 0.8*
      traceAdjYdTYd*Lambdax*Sqr(g1) + 2.4*traceAdjYeTYe*Lambdax*Sqr(g1) + 1.6*
      traceAdjYuTYu*Lambdax*Sqr(g1) + 0.8*MassB*traceYdAdjYd*Lambdax*Sqr(g1) -
      2.4*MassB*traceYeAdjYe*Lambdax*Sqr(g1) - 1.6*MassB*traceYuAdjYu*Lambdax*
      Sqr(g1) - 3.6*MassB*Lambdax*Sqr(g1)*Sqr(g2) - 3.6*MassWB*Lambdax*Sqr(g1)*
      Sqr(g2) + 32*traceAdjYdTYd*Lambdax*Sqr(g3) + 32*traceAdjYuTYu*Lambdax*Sqr
      (g3) - 32*MassG*traceYdAdjYd*Lambdax*Sqr(g3) - 32*MassG*traceYuAdjYu*
      Lambdax*Sqr(g3) + 4.14*Power(g1,4)*TLambdax + 7.5*Power(g2,4)*TLambdax -
      9*traceYdAdjYdYdAdjYd*TLambdax - 6*traceYdAdjYuYuAdjYd*TLambdax - 3*
      traceYeAdjYeYeAdjYe*TLambdax - 9*traceYuAdjYuYuAdjYu*TLambdax - 0.4*
      traceYdAdjYd*Sqr(g1)*TLambdax + 1.2*traceYeAdjYe*Sqr(g1)*TLambdax + 0.8*
      traceYuAdjYu*Sqr(g1)*TLambdax + 1.8*Sqr(g1)*Sqr(g2)*TLambdax + 16*
      traceYdAdjYd*Sqr(g3)*TLambdax + 16*traceYuAdjYu*Sqr(g3)*TLambdax - 50*Sqr
      (Conj(Lambdax))*Sqr(Lambdax)*TLambdax - 8*Kappa*Sqr(Conj(Kappa))*(4*
      Lambdax*TKappa + Kappa*TLambdax) - 0.6*AbsSqr(Lambdax)*(2*Lambdax*(2*
      MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*traceAdjYuTYu + 2*
      MassWB*Sqr(g2))) + (-6*Sqr(g1) + 15*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 2*Sqr(g2)))*TLambdax + 20*Conj(Kappa)*(2*Lambdax*TKappa +
      3*Kappa*TLambdax))));


   return beta_TLambdax;
}

/**
 * Calculates the three-loop beta function of TLambdax.
 *
 * @return three-loop beta function
 */
double NMSSMtower_soft_parameters::calc_beta_TLambdax_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
