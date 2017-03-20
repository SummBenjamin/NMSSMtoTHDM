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

// File generated at Mon 20 Mar 2017 08:23:08

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
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSMtower_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*((3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3))*TYd + 0.13333333333333333*Yd*(45*traceAdjYdTYd + 15*traceAdjYeTYe +
      7*MassB*Sqr(g1) + 45*MassWB*Sqr(g2) + 80*MassG*Sqr(g3) + 15*Conj(Lambdax
      )*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*
      Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSMtower_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(-0.044444444444444446*Yd*(287*Power(g1,4)*MassB +
      5*(-32*Power(g3,4)*MassG + 135*Power(g2,4)*MassWB + 27*(6*
      traceYdAdjYdTYdAdjYd + traceYdAdjYuTYuAdjYd + 2*traceYeAdjYeTYeAdjYe +
      traceYuAdjYdTYdAdjYu) - 144*(traceAdjYdTYd - MassG*traceYdAdjYd)*Sqr(g3)
      + 72*(MassG + MassWB)*Sqr(g2)*Sqr(g3)) + Sqr(g1)*(45*(MassB + MassWB)*Sqr
      (g2) + 2*(9*(traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*
      MassB*traceYeAdjYe) + 20*(MassB + MassG)*Sqr(g3)))) + (3.188888888888889*
      Power(g1,4) + 7.5*Power(g2,4) - 1.7777777777777777*Power(g3,4) - 9*
      traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + Sqr
      (g1)*(-0.4*traceYdAdjYd + 1.2*traceYeAdjYe + Sqr(g2) + 0.8888888888888888
      *Sqr(g3)) + 16*traceYdAdjYd*Sqr(g3) + 8*Sqr(g2)*Sqr(g3))*TYd - 3*Lambdax*
      Sqr(Conj(Lambdax))*(Lambdax*TYd + 4*Yd*TLambdax) - Conj(Lambdax)*((3*
      traceYuAdjYu + 2*AbsSqr(Kappa))*Lambdax*TYd + 2*Yd*(3*(traceAdjYuTYu*
      Lambdax + traceYuAdjYu*TLambdax) + 2*Conj(Kappa)*(Lambdax*TKappa + Kappa*
      TLambdax))) - 0.4*(4*MassB*Sqr(g1) + 15*(3*traceAdjYdTYd + traceAdjYeTYe
      + 2*MassWB*Sqr(g2)) + 15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) + (
      -12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*
      Sqr(g2))*(Yd*Yd.adjoint()*TYd) + (-6*traceAdjYuTYu - 1.6*MassB*Sqr(g1) -
      2*Conj(Lambdax)*TLambdax)*(Yd*Yu.adjoint()*Yu) + (-6*traceYuAdjYu - 2*
      AbsSqr(Lambdax) + 1.6*Sqr(g1))*(Yd*Yu.adjoint()*TYu) + (-15*traceYdAdjYd
      - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 12*Sqr(g2))*(TYd*
      Yd.adjoint()*Yd) + (-3*traceYuAdjYu - AbsSqr(Lambdax) + 0.8*Sqr(g1))*(TYd
      *Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*
      Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*
      TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd
      ) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYd;
}

/**
 * Calculates the three-loop beta function of TYd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> NMSSMtower_soft_parameters::calc_beta_TYd_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
