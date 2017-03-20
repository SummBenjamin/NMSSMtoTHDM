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

// File generated at Mon 20 Mar 2017 08:52:43

/**
 * @file NMSSMtower_mass_eigenstates.cpp
 * @brief implementation of the NMSSMtower model class
 *
 * Contains the definition of the NMSSMtower model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 20 Mar 2017 08:52:43 with FlexibleSUSY
 * 1.7.3 (git commit: 4b72fcc4bd7cf8f9d865a5278b43566468f1243d) and SARAH 4.10.0 .
 */

#include "NMSSMtower_mass_eigenstates.hpp"
#include "eigen_utils.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "parallel.hpp"
#include "pv.hpp"
#include "functors.hpp"

#include "sfermions.hpp"
#include "mssm_twoloophiggs.hpp"
#include "nmssm_twoloophiggs.hpp"



#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

#define CLASSNAME NMSSMtower_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS     two_loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     two_loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     two_loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU two_loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION          two_loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS  1

CLASSNAME::NMSSMtower_mass_eigenstates(const NMSSMtower_input_parameters& input_)
   : NMSSMtower_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , calculate_bsm_pole_masses(true)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(NMSSMtower_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSd(Eigen::Array<
      double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>::Zero()), MSu(
      Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()),
      Mhh(Eigen::Array<double,3,1>::Zero()), MAh(Eigen::Array<double,3,1>::Zero())
      , MHpm(Eigen::Array<double,2,1>::Zero()), MChi(Eigen::Array<double,5,1>
      ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,3
      ,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double
      ,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,5,5>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,
      2>::Zero()), UP(Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~NMSSMtower_mass_eigenstates()
{
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::do_calculate_bsm_pole_masses(bool flag)
{
   calculate_bsm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_bsm_pole_masses() const
{
   return calculate_bsm_pole_masses;
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_two_loop_corrections(const Two_loop_corrections& two_loop_corrections_)
{
   two_loop_corrections = two_loop_corrections_;
}

const Two_loop_corrections& CLASSNAME::get_two_loop_corrections() const
{
   return two_loop_corrections;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_ewsb_iterations() const
{
   return number_of_ewsb_iterations;
}

void CLASSNAME::set_number_of_mass_iterations(std::size_t iterations)
{
   number_of_mass_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_mass_iterations() const
{
   return number_of_mass_iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void CLASSNAME::set_pole_mass_loop_order(unsigned loop_order)
{
   pole_mass_loop_order = loop_order;
}

unsigned CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const NMSSMtower_physical& CLASSNAME::get_physical() const
{
   return physical;
}

NMSSMtower_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const NMSSMtower_physical& physical_)
{
   physical = physical_;
}

const Problems<NMSSMtower_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<NMSSMtower_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));
      tadpole[2] -= Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {
         const auto tadpole_2l(tadpole_hh_2loop());
         tadpole[0] -= tadpole_2l(0);
         tadpole[1] -= tadpole_2l(1);
         tadpole[2] -= tadpole_2l(2);

      }
   }

   return tadpole;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   const auto tadpole_(tadpole_equations());
   std::copy(tadpole_.data(), tadpole_.data() + number_of_ewsb_equations, tadpole);
}

/**
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::EWSB_args struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::EWSB_args struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   NMSSMtower_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mHd2(gsl_vector_get(x, 0));
   model->set_mHu2(gsl_vector_get(x, 1));
   model->set_ms2(gsl_vector_get(x, 2));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   const auto tadpole(model->tadpole_equations());

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return IsFinite(tadpole) ? GSL_SUCCESS : GSL_EDOM;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, fixed_point_iterator::Convergence_tester_relative(ewsb_iteration_precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden))
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   const auto x_init(ewsb_initial_guess());

   VERBOSE_MSG("Solving EWSB equations ...");
   VERBOSE_MSG("\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i].get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const Eigen::Matrix<double, number_of_ewsb_equations, 1>& x_init
)
{
   const int status = solver->solve(&x_init[0]);

   mHd2 = solver->get_solution(0);
   mHu2 = solver->get_solution(1);
   ms2 = solver->get_solution(2);


   return status;
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Power(vd,3)*
      Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(Lambdax)*Sqr(vS) + 10*vu*Conj
      (Lambdax)*Kappa*Sqr(vS) + 10*vu*Conj(Kappa)*Lambdax*Sqr(vS) - 20*vd*AbsSqr(
      Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) +
      14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Power(vu,3)*
      Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(
      g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) + 10*vd*
      Conj(Lambdax)*Kappa*Sqr(vS) + 10*vd*Conj(Kappa)*Lambdax*Sqr(vS) +
      14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(-4*Power(vS,3)*AbsSqr(Kappa) + 1.4142135623730951*vd*vu*Conj
      (TLambdax) + 2*vd*vS*vu*Conj(Lambdax)*Kappa + 2*vd*vS*vu*Conj(Kappa)*Lambdax
      - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(TKappa)*Sqr(vS) -
      2*vS*AbsSqr(Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*TKappa +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      ms2 = old_ms2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = 0;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Power(vd,3)*
      Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(Lambdax)*Sqr(vS) + 10*vu*Conj
      (Lambdax)*Kappa*Sqr(vS) + 10*vu*Conj(Kappa)*Lambdax*Sqr(vS) - 20*vd*AbsSqr(
      Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) +
      14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Power(vu,3)*
      Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(
      g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) + 10*vd*
      Conj(Lambdax)*Kappa*Sqr(vS) + 10*vd*Conj(Kappa)*Lambdax*Sqr(vS) +
      14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(-4*Power(vS,3)*AbsSqr(Kappa) + 1.4142135623730951*vd*vu*Conj
      (TLambdax) + 2*vd*vS*vu*Conj(Lambdax)*Kappa + 2*vd*vS*vu*Conj(Kappa)*Lambdax
      - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(TKappa)*Sqr(vS) -
      2*vS*AbsSqr(Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*TKappa +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      ms2 = old_ms2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_initial_guess()
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> x_init(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   x_init[0] = mHd2;
   x_init[1] = mHu2;
   x_init[2] = ms2;


   return x_init;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_step() const
{
   double tadpole[number_of_ewsb_equations] = { 0. };
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh(0));
      tadpole[1] += Re(tadpole_hh(1));
      tadpole[2] += Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {
         const auto tadpole_2l(tadpole_hh_2loop());
         tadpole[0] += tadpole_2l(0);
         tadpole[1] += tadpole_2l(1);
         tadpole[2] += tadpole_2l(2);

      }
   }

   double mHd2;
   double mHu2;
   double ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) + 40*tadpole[0] -
      3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(Lambdax)*Sqr(vS
      ) + 10*vu*Conj(Lambdax)*Kappa*Sqr(vS) + 10*vu*Conj(Kappa)*Lambdax*Sqr(vS) -
      20*vd*AbsSqr(Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu)
      + 14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) + 40*tadpole[1] -
      3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*vu*AbsSqr(Lambdax)*Sqr(vd
      ) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(
      vS) + 10*vd*Conj(Lambdax)*Kappa*Sqr(vS) + 10*vd*Conj(Kappa)*Lambdax*Sqr(vS)
      + 14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(-4*Power(vS,3)*AbsSqr(Kappa) + 1.4142135623730951*vd*vu*Conj
      (TLambdax) + 2*vd*vS*vu*Conj(Lambdax)*Kappa + 2*vd*vS*vu*Conj(Kappa)*Lambdax
      + 4*tadpole[2] - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 1.4142135623730951*Conj(
      TKappa)*Sqr(vS) - 2*vS*AbsSqr(Lambdax)*Sqr(vu) - 1.4142135623730951*Sqr(vS)*
      TKappa + 1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mHd2;
   ewsb_parameters[1] = mHu2;
   ewsb_parameters[2] = ms2;


   return ewsb_parameters;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param x old EWSB output parameters
 * @param params further function parameters
 * @param f new EWSB output parameters
 *
 * @return Returns status of CLASSNAME::ewsb_step
 */
int CLASSNAME::ewsb_step(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   NMSSMtower_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mHd2 = gsl_vector_get(x, 0);
   const double mHu2 = gsl_vector_get(x, 1);
   const double ms2 = gsl_vector_get(x, 2);

   model->set_mHd2(mHd2);
   model->set_mHu2(mHu2);
   model->set_ms2(ms2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = mHd2;
   ewsb_parameters[1] = mHu2;
   ewsb_parameters[2] = ms2;


   int status = GSL_SUCCESS;

   try {
      ewsb_parameters = model->ewsb_step();
      status = GSL_SUCCESS;
   } catch (...) {
      status = GSL_EDOM;
   }

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "NMSSMtower\n"
           "========================================\n";
   NMSSMtower_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 */

double CLASSNAME::A0(double m) const
{
   return passarino_veltman::ReA0(m*m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const
{
   return passarino_veltman::ReB0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const
{
   return passarino_veltman::ReB1(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const
{
   return passarino_veltman::ReB00(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const
{
   return passarino_veltman::ReB22(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const
{
   return passarino_veltman::ReH0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const
{
   return passarino_veltman::ReF0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const
{
   return passarino_veltman::ReG0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto old_mHd2 = mHd2;
   const auto old_mHu2 = mHu2;
   const auto old_ms2 = ms2;

   solve_ewsb_tree_level_custom();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
   calculate_MFv();
   calculate_MGlu();
   calculate_MVG();

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;
   ms2 = old_ms2;

}

/**
 * Backward compatibility routine which finds the DRbar mass
 * eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_parameters()
{
   calculate_DRbar_masses();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   auto obj_ptr = this;

   std::future<void> fut_MVG;
   std::future<void> fut_MGlu;
   std::future<void> fut_MFv;
   std::future<void> fut_MVP;
   std::future<void> fut_MVZ;
   std::future<void> fut_MSd;
   std::future<void> fut_MSv;
   std::future<void> fut_MSu;
   std::future<void> fut_MSe;
   std::future<void> fut_Mhh;
   std::future<void> fut_MAh;
   std::future<void> fut_MHpm;
   std::future<void> fut_MChi;
   std::future<void> fut_MCha;
   std::future<void> fut_MFe;
   std::future<void> fut_MFd;
   std::future<void> fut_MFu;
   std::future<void> fut_MVWm;

   if (calculate_bsm_pole_masses) {
      fut_MAh = run_async([obj_ptr] () { obj_ptr->calculate_MAh_pole(); });
      fut_MCha = run_async([obj_ptr] () { obj_ptr->calculate_MCha_pole(); });
      fut_MChi = run_async([obj_ptr] () { obj_ptr->calculate_MChi_pole(); });
      fut_MGlu = run_async([obj_ptr] () { obj_ptr->calculate_MGlu_pole(); });
      fut_Mhh = run_async([obj_ptr] () { obj_ptr->calculate_Mhh_pole(); });
      fut_MHpm = run_async([obj_ptr] () { obj_ptr->calculate_MHpm_pole(); });
      fut_MSd = run_async([obj_ptr] () { obj_ptr->calculate_MSd_pole(); });
      fut_MSe = run_async([obj_ptr] () { obj_ptr->calculate_MSe_pole(); });
      fut_MSu = run_async([obj_ptr] () { obj_ptr->calculate_MSu_pole(); });
      fut_MSv = run_async([obj_ptr] () { obj_ptr->calculate_MSv_pole(); });
   }

   if (calculate_sm_pole_masses) {
      fut_MVG = run_async([obj_ptr] () { obj_ptr->calculate_MVG_pole(); });
      fut_MFv = run_async([obj_ptr] () { obj_ptr->calculate_MFv_pole(); });
      fut_MVP = run_async([obj_ptr] () { obj_ptr->calculate_MVP_pole(); });
      fut_MVZ = run_async([obj_ptr] () { obj_ptr->calculate_MVZ_pole(); });
      fut_MFe = run_async([obj_ptr] () { obj_ptr->calculate_MFe_pole(); });
      fut_MFd = run_async([obj_ptr] () { obj_ptr->calculate_MFd_pole(); });
      fut_MFu = run_async([obj_ptr] () { obj_ptr->calculate_MFu_pole(); });
      fut_MVWm = run_async([obj_ptr] () { obj_ptr->calculate_MVWm_pole(); });
   }

   if (fut_MAh.valid()) fut_MAh.get();
   if (fut_MCha.valid()) fut_MCha.get();
   if (fut_MChi.valid()) fut_MChi.get();
   if (fut_MGlu.valid()) fut_MGlu.get();
   if (fut_Mhh.valid()) fut_Mhh.get();
   if (fut_MHpm.valid()) fut_MHpm.get();
   if (fut_MSd.valid()) fut_MSd.get();
   if (fut_MSe.valid()) fut_MSe.get();
   if (fut_MSu.valid()) fut_MSu.get();
   if (fut_MSv.valid()) fut_MSv.get();
   if (fut_MVG.valid()) fut_MVG.get();
   if (fut_MFv.valid()) fut_MFv.get();
   if (fut_MVP.valid()) fut_MVP.get();
   if (fut_MVZ.valid()) fut_MVZ.get();
   if (fut_MFe.valid()) fut_MFe.get();
   if (fut_MFd.valid()) fut_MFd.get();
   if (fut_MFu.valid()) fut_MFu.get();
   if (fut_MVWm.valid()) fut_MVWm.get();

#else
   if (calculate_bsm_pole_masses) {
      calculate_MAh_pole();
      calculate_MCha_pole();
      calculate_MChi_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHpm_pole();
      calculate_MSd_pole();
      calculate_MSe_pole();
      calculate_MSu_pole();
      calculate_MSv_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFe_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MVWm_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSv) = MSv;
   PHYSICAL(ZV) = ZV;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(ZEL) = ZEL;
   PHYSICAL(ZER) = ZER;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(ZDL) = ZDL;
   PHYSICAL(ZDR) = ZDR;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(ZUL) = ZUL;
   PHYSICAL(ZUR) = ZUR;
   PHYSICAL(MVWm) = MVWm;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHpm, ZP);

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Sd);
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Sv);
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Su);
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Se);
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::hh);
   if (PHYSICAL(MAh).tail<2>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Ah);
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) problems.flag_tachyon(NMSSMtower_info::Hpm);

}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,5,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,5,5>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   NMSSMtower_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MSd(0) = pars(5);
   MSd(1) = pars(6);
   MSd(2) = pars(7);
   MSd(3) = pars(8);
   MSd(4) = pars(9);
   MSd(5) = pars(10);
   MSv(0) = pars(11);
   MSv(1) = pars(12);
   MSv(2) = pars(13);
   MSu(0) = pars(14);
   MSu(1) = pars(15);
   MSu(2) = pars(16);
   MSu(3) = pars(17);
   MSu(4) = pars(18);
   MSu(5) = pars(19);
   MSe(0) = pars(20);
   MSe(1) = pars(21);
   MSe(2) = pars(22);
   MSe(3) = pars(23);
   MSe(4) = pars(24);
   MSe(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   Mhh(2) = pars(28);
   MAh(0) = pars(29);
   MAh(1) = pars(30);
   MAh(2) = pars(31);
   MHpm(0) = pars(32);
   MHpm(1) = pars(33);
   MChi(0) = pars(34);
   MChi(1) = pars(35);
   MChi(2) = pars(36);
   MChi(3) = pars(37);
   MChi(4) = pars(38);
   MCha(0) = pars(39);
   MCha(1) = pars(40);
   MFe(0) = pars(41);
   MFe(1) = pars(42);
   MFe(2) = pars(43);
   MFd(0) = pars(44);
   MFd(1) = pars(45);
   MFd(2) = pars(46);
   MFu(0) = pars(47);
   MFu(1) = pars(48);
   MFu(2) = pars(49);
   MVWm = pars(50);
   MVP = pars(51);
   MVZ = pars(52);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(53);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSd(0);
   pars(6) = MSd(1);
   pars(7) = MSd(2);
   pars(8) = MSd(3);
   pars(9) = MSd(4);
   pars(10) = MSd(5);
   pars(11) = MSv(0);
   pars(12) = MSv(1);
   pars(13) = MSv(2);
   pars(14) = MSu(0);
   pars(15) = MSu(1);
   pars(16) = MSu(2);
   pars(17) = MSu(3);
   pars(18) = MSu(4);
   pars(19) = MSu(5);
   pars(20) = MSe(0);
   pars(21) = MSe(1);
   pars(22) = MSe(2);
   pars(23) = MSe(3);
   pars(24) = MSe(4);
   pars(25) = MSe(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = Mhh(2);
   pars(29) = MAh(0);
   pars(30) = MAh(1);
   pars(31) = MAh(2);
   pars(32) = MHpm(0);
   pars(33) = MHpm(1);
   pars(34) = MChi(0);
   pars(35) = MChi(1);
   pars(36) = MChi(2);
   pars(37) = MChi(3);
   pars(38) = MChi(4);
   pars(39) = MCha(0);
   pars(40) = MCha(1);
   pars(41) = MFe(0);
   pars(42) = MFe(1);
   pars(43) = MFe(2);
   pars(44) = MFd(0);
   pars(45) = MFd(1);
   pars(46) = MFd(2);
   pars(47) = MFu(0);
   pars(48) = MFu(1);
   pars(49) = MFu(2);
   pars(50) = MVWm;
   pars(51) = MVP;
   pars(52) = MVZ;

   return pars;
}

std::string CLASSNAME::name() const
{
   return "NMSSMtower";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   NMSSMtower_soft_parameters::run_to(scale, eps);
}


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_goldstone;
   MHpm_goldstone(0) = MVWm;

   return remove_if_equal(MHpm, MHpm_goldstone);
}

Eigen::Array<double,2,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;

   return remove_if_equal(MAh, MAh_goldstone);
}


/**
 * @brief finds the LSP and returns it's mass
 *
 * This function finds the lightest supersymmetric particle (LSP) and
 * returns it's mass.  The corresponding particle type is retured in
 * the reference parameter.  The list of potential LSPs is set in the
 * model file varible PotentialLSPParticles.  For this model it is set
 * to:
 * {Chi, Sv, Su, Sd, Se, Cha, Glu}
 *
 * @param particle_type particle type
 * @return mass of LSP
 */
double CLASSNAME::get_lsp(NMSSMtower_info::Particles& particle_type) const
{
   double lsp_mass = std::numeric_limits<double>::max();
   double tmp_mass;
   particle_type = NMSSMtower_info::NUMBER_OF_PARTICLES;

   tmp_mass = Abs(PHYSICAL(MChi(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Chi;
   }

   tmp_mass = Abs(PHYSICAL(MSv(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Sv;
   }

   tmp_mass = Abs(PHYSICAL(MSu(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Su;
   }

   tmp_mass = Abs(PHYSICAL(MSd(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Sd;
   }

   tmp_mass = Abs(PHYSICAL(MSe(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Se;
   }

   tmp_mass = Abs(PHYSICAL(MCha(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Cha;
   }

   tmp_mass = Abs(PHYSICAL(MGlu));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = NMSSMtower_info::Glu;
   }

   return lsp_mass;
}


double CLASSNAME::get_mass_matrix_VG() const
{
   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{
   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

double CLASSNAME::get_mass_matrix_Glu() const
{
   const double mass_matrix_Glu = Re(MassG);

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{
   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   MGlu = calculate_majorana_singlet_mass(mass_matrix_Glu, PhaseGlu);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   MFv.setConstant(0);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)
      ) + AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) +
      Conj(Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) +
      Conj(Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vS*vu
      *Conj(Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) - 0.5*vS*vu
      *Conj(Yd(1,0))*Lambdax;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) - 0.5*vS*vu
      *Conj(Yd(2,0))*Lambdax;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) +
      Conj(Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) - 0.5*vS*vu
      *Conj(Yd(0,1))*Lambdax;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vS*vu
      *Conj(Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) - 0.5*vS*vu
      *Conj(Yd(2,1))*Lambdax;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) - 0.5*vS*vu
      *Conj(Yd(0,2))*Lambdax;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) - 0.5*vS*vu
      *Conj(Yd(1,2))*Lambdax;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vS*vu
      *Conj(Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)
      ) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Sd, eigenvalue_error >
      precision * Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif


   if (MSd.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Sv, eigenvalue_error >
      precision * Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif


   if (MSv.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) +
      Conj(Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) +
      Conj(Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vS
      *Conj(Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) - 0.5*vd*vS
      *Conj(Yu(1,0))*Lambdax;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) - 0.5*vd*vS
      *Conj(Yu(2,0))*Lambdax;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) +
      Conj(Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) - 0.5*vd*vS
      *Conj(Yu(0,1))*Lambdax;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vS
      *Conj(Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) - 0.5*vd*vS
      *Conj(Yu(2,1))*Lambdax;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) - 0.5*vd*vS
      *Conj(Yu(0,2))*Lambdax;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) - 0.5*vd*vS
      *Conj(Yu(1,2))*Lambdax;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vS
      *Conj(Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      0,0)) + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      1,0)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      2,0)) + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Su, eigenvalue_error >
      precision * Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif


   if (MSu.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)
      ) + AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) +
      Conj(Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) +
      Conj(Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vS*vu
      *Conj(Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) - 0.5*vS*vu
      *Conj(Ye(1,0))*Lambdax;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) - 0.5*vS*vu
      *Conj(Ye(2,0))*Lambdax;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) +
      Conj(Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) - 0.5*vS*vu
      *Conj(Ye(0,1))*Lambdax;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vS*vu
      *Conj(Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) - 0.5*vS*vu
      *Conj(Ye(2,1))*Lambdax;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) - 0.5*vS*vu
      *Conj(Ye(0,2))*Lambdax;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) - 0.5*vS*vu
      *Conj(Ye(1,2))*Lambdax;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vS*vu
      *Conj(Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)
      ) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Se, eigenvalue_error >
      precision * Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif


   if (MSe.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr
      (vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*
      Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vS*
      Conj(TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) - 0.25*Conj(
      Lambdax)*Kappa*Sqr(vS) - 0.25*Conj(Kappa)*Lambdax*Sqr(vS) -
      0.35355339059327373*vS*TLambdax;
   mass_matrix_hh(0,2) = vd*vS*AbsSqr(Lambdax) - 0.35355339059327373*vu*
      Conj(TLambdax) - 0.5*vS*vu*Conj(Lambdax)*Kappa - 0.5*vS*vu*Conj(Kappa)*
      Lambdax - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(
      g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.225
      *Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = vS*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*
      Conj(TLambdax) - 0.5*vd*vS*Conj(Lambdax)*Kappa - 0.5*vd*vS*Conj(Kappa)*
      Lambdax - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.7071067811865475*vS*Conj(TKappa) - 0.5*
      vd*vu*Conj(Lambdax)*Kappa - 0.5*vd*vu*Conj(Kappa)*Lambdax + 0.5*AbsSqr(
      Lambdax)*Sqr(vd) + 3*AbsSqr(Kappa)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu)
      + 0.7071067811865475*vS*TKappa;

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::hh, eigenvalue_error >
      precision * Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif


   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.3872983346207417*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) +
      0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)
      *Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW())
      ) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vS*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) + 0.25*Conj(
      Lambdax)*Kappa*Sqr(vS) + 0.25*Conj(Kappa)*Lambdax*Sqr(vS) - 0.25*vd*vu*
      Sqr(g2)*Sqr(Cos(ThetaW())) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW())) +
      0.35355339059327373*vS*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) - 0.5*vS*
      vu*Conj(Lambdax)*Kappa - 0.5*vS*vu*Conj(Kappa)*Lambdax +
      0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(
      g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW
      ())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) - 0.5*vd*
      vS*Conj(Lambdax)*Kappa - 0.5*vd*vS*Conj(Kappa)*Lambdax +
      0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 - 0.7071067811865475*vS*Conj(TKappa) + 0.5*
      vd*vu*Conj(Lambdax)*Kappa + 0.5*vd*vu*Conj(Kappa)*Lambdax + 0.5*AbsSqr(
      Lambdax)*Sqr(vd) + AbsSqr(Kappa)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) -
      0.7071067811865475*vS*TKappa;

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Ah, eigenvalue_error >
      precision * Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif


   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*
      Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr
      (g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475
      *vS*Conj(TLambdax) + 0.5*Conj(Lambdax)*Kappa*Sqr(vS);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr
      (g2)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Hpm, eigenvalue_error >
      precision * Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif


   if (MHpm.minCoeff() < 0.) {
      problems.flag_tachyon(NMSSMtower_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,5,5> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,5,5> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(0,4) = 0;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.7071067811865475*vS*Lambdax;
   mass_matrix_Chi(2,4) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(4,4) = 1.4142135623730951*vS*Kappa;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Chi, eigenvalue_error >
      precision * Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = 0.7071067811865475*vS*Lambdax;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Cha, eigenvalue_error >
      precision * Abs(MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*vd*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*vd*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*vd*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*vd*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*vd*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*vd*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*vd*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*vd*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*vd*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Fe, eigenvalue_error >
      precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*vd*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*vd*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*vd*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*vd*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*vd*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*vd*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*vd*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*vd*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*vd*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*vu*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*vu*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*vu*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*vu*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*vu*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*vu*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*vu*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*vu*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*vu*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR, eigenvalue_error);
   problems.flag_bad_mass(NMSSMtower_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + Sqr(vu)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_tachyon(NMSSMtower_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void CLASSNAME::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> MVPVZ;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error
      );
   ZZ.transposeInPlace();
#else
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
   ZZ.transposeInPlace();
#endif


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = Re(mHd2*vd - 0.35355339059327373*vS*vu*Conj(TLambdax) +
      0.075*Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) + 0.5*vd*AbsSqr(
      Lambdax)*Sqr(vS) - 0.25*vu*Conj(Lambdax)*Kappa*Sqr(vS) - 0.25*vu*Conj(Kappa)
      *Lambdax*Sqr(vS) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu)
      - 0.125*vd*Sqr(g2)*Sqr(vu) - 0.35355339059327373*vS*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = Re(mHu2*vu - 0.35355339059327373*vd*vS*Conj(TLambdax) +
      0.075*Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) + 0.5*vu*AbsSqr(
      Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5
      *vu*AbsSqr(Lambdax)*Sqr(vS) - 0.25*vd*Conj(Lambdax)*Kappa*Sqr(vS) - 0.25*vd*
      Conj(Kappa)*Lambdax*Sqr(vS) - 0.35355339059327373*vd*vS*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   double result = Re(ms2*vS + Power(vS,3)*AbsSqr(Kappa) - 0.35355339059327373*
      vd*vu*Conj(TLambdax) - 0.5*vd*vS*vu*Conj(Lambdax)*Kappa - 0.5*vd*vS*vu*Conj(
      Kappa)*Lambdax + 0.5*vS*AbsSqr(Lambdax)*Sqr(vd) + 0.35355339059327373*Conj(
      TKappa)*Sqr(vS) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vu) + 0.35355339059327373*Sqr(
      vS)*TKappa - 0.35355339059327373*vd*vu*TLambdax);

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   std::complex<double> tmp_1;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_0 += tmp_1;
   result += (0.13333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_0;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += 0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()));
   }

   return result;
}

double CLASSNAME::CpUSdconjUSdconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2;
   std::complex<double> tmp_3;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_2 += tmp_3;
   result += (0.1*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_2;
   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_6;
      std::complex<double> tmp_7;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_8;
         std::complex<double> tmp_9;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_9 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_8 += tmp_9;
         tmp_7 += (KroneckerDelta(gO2,3 + j2)) * tmp_8;
      }
      tmp_6 += tmp_7;
      tmp_5 += (KroneckerDelta(gO1,3 + j3)) * tmp_6;
   }
   tmp_4 += tmp_5;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_4;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_10;
   std::complex<double> tmp_11;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_11 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_10 += tmp_11;
   result += (-0.1*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_10;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_12;
      std::complex<double> tmp_13;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_13 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_12 += tmp_13;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_12;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_14;
   std::complex<double> tmp_15;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_15 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_14 += tmp_15;
   result += (0.1*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_14;
   std::complex<double> tmp_16;
   std::complex<double> tmp_17;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_18;
      std::complex<double> tmp_19;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_20;
         std::complex<double> tmp_21;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_21 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_20 += tmp_21;
         tmp_19 += (KroneckerDelta(gO2,3 + j2)) * tmp_20;
      }
      tmp_18 += tmp_19;
      tmp_17 += (KroneckerDelta(gO1,3 + j3)) * tmp_18;
   }
   tmp_16 += tmp_17;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_16;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_22;
      std::complex<double> tmp_23;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_23 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_22 += tmp_23;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_22;
   }
   std::complex<double> tmp_24;
   std::complex<double> tmp_25;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_25 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_24 += tmp_25;
   result += (-0.1*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_24;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_26;
      std::complex<double> tmp_27;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_27 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_26 += tmp_27;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,1)) * tmp_26;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_28;
      std::complex<double> tmp_29;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_29 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_28 += tmp_29;
      result += (-0.5*Lambdax*ZA(gI1,2)*ZA(gI2,1)) * tmp_28;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_30;
      std::complex<double> tmp_31;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_31 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_30 += tmp_31;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,1)*ZA(gI2,2)) * tmp_30;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_32;
      std::complex<double> tmp_33;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_33 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_32 += tmp_33;
      result += (-0.5*Lambdax*ZA(gI1,1)*ZA(gI2,2)) * tmp_32;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_34;
   std::complex<double> tmp_35;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_35 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_34 += tmp_35;
   result += (0.1*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_34;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_36;
   std::complex<double> tmp_37;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_37 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_36 += tmp_37;
   result += (0.1*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_36;
   std::complex<double> tmp_38;
   std::complex<double> tmp_39;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_40;
      std::complex<double> tmp_41;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_42;
         std::complex<double> tmp_43;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_43 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_42 += tmp_43;
         tmp_41 += (KroneckerDelta(gO2,3 + j2)) * tmp_42;
      }
      tmp_40 += tmp_41;
      tmp_39 += (KroneckerDelta(gO1,3 + j3)) * tmp_40;
   }
   tmp_38 += tmp_39;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_38;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_44;
      std::complex<double> tmp_45;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_45 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_44 += tmp_45;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_44;
   }
   std::complex<double> tmp_46;
   std::complex<double> tmp_47;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_47 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_46 += tmp_47;
   result += (-0.1*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_46;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_48;
      std::complex<double> tmp_49;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_49 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_48 += tmp_49;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) * tmp_48;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_50;
      std::complex<double> tmp_51;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_51 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_50 += tmp_51;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_50;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_52;
      std::complex<double> tmp_53;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_53 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_52 += tmp_53;
      result += (0.5*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) * tmp_52;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_54;
      std::complex<double> tmp_55;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_55 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_54 += tmp_55;
      result += (0.5*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_54;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_56;
      std::complex<double> tmp_57;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_57 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_56 += tmp_57;
      result += (UP(gI2,1)) * tmp_56;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_60;
      std::complex<double> tmp_61;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_61 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_60 += tmp_61;
      tmp_59 += (Conj(ZUL(gI1,j2))) * tmp_60;
   }
   tmp_58 += tmp_59;
   result += (Conj(UM(gI2,1))) * tmp_58;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_62;
   std::complex<double> tmp_63;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_63 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_62 += tmp_63;
   result += (-0.3651483716701107*g1*ZN(gI2,0)) * tmp_62;
   if (gO2 < 3) {
      std::complex<double> tmp_64;
      std::complex<double> tmp_65;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_65 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_64 += tmp_65;
      result += (-ZN(gI2,2)) * tmp_64;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_66;
   std::complex<double> tmp_67;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_69 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_68 += tmp_69;
      tmp_67 += (Conj(ZDL(gI1,j2))) * tmp_68;
   }
   tmp_66 += tmp_67;
   result += (-Conj(ZN(gI2,2))) * tmp_66;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_70;
   std::complex<double> tmp_72;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_72 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_70 += tmp_72;
   std::complex<double> tmp_71;
   std::complex<double> tmp_73;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_73 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_71 += tmp_73;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_70 * tmp_71;
   std::complex<double> tmp_74;
   std::complex<double> tmp_76;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_76 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_74 += tmp_76;
   std::complex<double> tmp_75;
   std::complex<double> tmp_77;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_77 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_75 += tmp_77;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_74 * tmp_75;
   std::complex<double> tmp_78;
   std::complex<double> tmp_80;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_80 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_78 += tmp_80;
   std::complex<double> tmp_79;
   std::complex<double> tmp_81;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_81 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_79 += tmp_81;
   result += (-0.05*Sqr(g1)) * tmp_78 * tmp_79;
   std::complex<double> tmp_82;
   std::complex<double> tmp_84;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_84 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_82 += tmp_84;
   std::complex<double> tmp_83;
   std::complex<double> tmp_85;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_85 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_83 += tmp_85;
   result += (-0.1*Sqr(g1)) * tmp_82 * tmp_83;
   std::complex<double> tmp_86;
   std::complex<double> tmp_88;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_88 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_86 += tmp_88;
   std::complex<double> tmp_87;
   std::complex<double> tmp_89;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_89 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_87 += tmp_89;
   result += (-0.05*Sqr(g1)) * tmp_86 * tmp_87;
   std::complex<double> tmp_90;
   std::complex<double> tmp_92;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_92 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_90 += tmp_92;
   std::complex<double> tmp_91;
   std::complex<double> tmp_93;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_93 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_91 += tmp_93;
   result += (-0.1*Sqr(g1)) * tmp_90 * tmp_91;
   std::complex<double> tmp_94;
   std::complex<double> tmp_96;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_96 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_94 += tmp_96;
   std::complex<double> tmp_95;
   std::complex<double> tmp_97;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_97 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_95 += tmp_97;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_94 * tmp_95;
   std::complex<double> tmp_98;
   std::complex<double> tmp_100;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_100 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_98 += tmp_100;
   std::complex<double> tmp_99;
   std::complex<double> tmp_101;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_101 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_99 += tmp_101;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_98 * tmp_99;
   std::complex<double> tmp_102;
   std::complex<double> tmp_104;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_105;
      std::complex<double> tmp_106;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_106 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_105 += tmp_106;
      tmp_104 += (Conj(ZD(gI2,j2))) * tmp_105;
   }
   tmp_102 += tmp_104;
   std::complex<double> tmp_103;
   std::complex<double> tmp_107;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_108;
      std::complex<double> tmp_109;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_109 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_108 += tmp_109;
      tmp_107 += (ZD(gI1,j4)) * tmp_108;
   }
   tmp_103 += tmp_107;
   result += (-1) * tmp_102 * tmp_103;
   if (gO1 < 3) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_111 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_110 += tmp_111;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_110;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_113 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_112 += tmp_113;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_112;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_114;
      std::complex<double> tmp_115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_115 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_114 += tmp_115;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_114;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_117 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_116 += tmp_117;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_116;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_118;
      std::complex<double> tmp_119;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_119 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_118 += tmp_119;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_118;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_121 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_120 += tmp_121;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_120;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_124;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_124 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_122 += tmp_124;
      std::complex<double> tmp_123;
      std::complex<double> tmp_125;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_126;
         std::complex<double> tmp_127;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_127 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_126 += tmp_127;
         tmp_125 += (ZD(gI1,j4)) * tmp_126;
      }
      tmp_123 += tmp_125;
      result += (-3) * tmp_122 * tmp_123;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_129;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_129 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_128 += tmp_129;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_128;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_130;
      std::complex<double> tmp_131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_131 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_130 += tmp_131;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_130;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_133 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_132 += tmp_133;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_132;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_135 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_134 += tmp_135;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_134;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_136;
      std::complex<double> tmp_138;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_139;
         std::complex<double> tmp_140;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_140 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_139 += tmp_140;
         tmp_138 += (Conj(ZD(gI2,j2))) * tmp_139;
      }
      tmp_136 += tmp_138;
      std::complex<double> tmp_137;
      std::complex<double> tmp_141;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_141 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_137 += tmp_141;
      result += (-3) * tmp_136 * tmp_137;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_142;
      std::complex<double> tmp_144;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_144 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_142 += tmp_144;
      std::complex<double> tmp_143;
      std::complex<double> tmp_145;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_145 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_143 += tmp_145;
      result += (-1) * tmp_142 * tmp_143;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_146;
      std::complex<double> tmp_147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_147 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_146 += tmp_147;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_146;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_148;
      std::complex<double> tmp_149;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_149 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_148 += tmp_149;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_148;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_150;
      std::complex<double> tmp_151;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_151 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_150 += tmp_151;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_150;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_153 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_152 += tmp_153;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_152;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*ZD(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -1.3333333333333333*Conj(ZD(gI2,gO2))*Sqr(g3)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_154;
   std::complex<double> tmp_156;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_156 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_154 += tmp_156;
   std::complex<double> tmp_155;
   std::complex<double> tmp_157;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_157 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_155 += tmp_157;
   result += (0.05*Sqr(g1)) * tmp_154 * tmp_155;
   std::complex<double> tmp_158;
   std::complex<double> tmp_160;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_160 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_158 += tmp_160;
   std::complex<double> tmp_159;
   std::complex<double> tmp_161;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_161 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_159 += tmp_161;
   result += (-0.1*Sqr(g1)) * tmp_158 * tmp_159;
   std::complex<double> tmp_162;
   std::complex<double> tmp_164;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_164 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_162 += tmp_164;
   std::complex<double> tmp_163;
   std::complex<double> tmp_165;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_165 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_163 += tmp_165;
   result += (0.05*Sqr(g1)) * tmp_162 * tmp_163;
   std::complex<double> tmp_166;
   std::complex<double> tmp_168;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_168 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_166 += tmp_168;
   std::complex<double> tmp_167;
   std::complex<double> tmp_169;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_169 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_167 += tmp_169;
   result += (-0.1*Sqr(g1)) * tmp_166 * tmp_167;
   if (gO1 < 3) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_171;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_171 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_170 += tmp_171;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_170;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_172;
      std::complex<double> tmp_173;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_173 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_172 += tmp_173;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_172;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_175;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_175 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_174 += tmp_175;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_174;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_176;
      std::complex<double> tmp_177;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_177 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_176 += tmp_177;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_176;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_178;
      std::complex<double> tmp_179;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_179 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_178 += tmp_179;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_181;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_181 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_180 += tmp_181;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_180;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_182;
      std::complex<double> tmp_184;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_184 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_182 += tmp_184;
      std::complex<double> tmp_183;
      std::complex<double> tmp_185;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_186;
         std::complex<double> tmp_187;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_187 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_186 += tmp_187;
         tmp_185 += (ZE(gI1,j4)) * tmp_186;
      }
      tmp_183 += tmp_185;
      result += (-1) * tmp_182 * tmp_183;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_188;
      std::complex<double> tmp_190;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_191;
         std::complex<double> tmp_192;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_192 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_191 += tmp_192;
         tmp_190 += (Conj(ZE(gI2,j2))) * tmp_191;
      }
      tmp_188 += tmp_190;
      std::complex<double> tmp_189;
      std::complex<double> tmp_193;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_193 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_189 += tmp_193;
      result += (-1) * tmp_188 * tmp_189;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_194;
   std::complex<double> tmp_196;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_196 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_194 += tmp_196;
   std::complex<double> tmp_195;
   std::complex<double> tmp_197;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_197 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_195 += tmp_197;
   result += (-0.05*Sqr(g1)) * tmp_194 * tmp_195;
   std::complex<double> tmp_198;
   std::complex<double> tmp_200;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_200 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_198 += tmp_200;
   std::complex<double> tmp_199;
   std::complex<double> tmp_201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_201 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_199 += tmp_201;
   result += (0.2*Sqr(g1)) * tmp_198 * tmp_199;
   std::complex<double> tmp_202;
   std::complex<double> tmp_204;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_204 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_202 += tmp_204;
   std::complex<double> tmp_203;
   std::complex<double> tmp_205;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_205 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_203 += tmp_205;
   result += (-0.05*Sqr(g1)) * tmp_202 * tmp_203;
   std::complex<double> tmp_206;
   std::complex<double> tmp_208;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_208 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_206 += tmp_208;
   std::complex<double> tmp_207;
   std::complex<double> tmp_209;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_209 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_207 += tmp_209;
   result += (0.2*Sqr(g1)) * tmp_206 * tmp_207;
   std::complex<double> tmp_210;
   std::complex<double> tmp_212;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_213;
      std::complex<double> tmp_214;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_214 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_213 += tmp_214;
      tmp_212 += (Conj(ZU(gI2,j2))) * tmp_213;
   }
   tmp_210 += tmp_212;
   std::complex<double> tmp_211;
   std::complex<double> tmp_215;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_216;
      std::complex<double> tmp_217;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_217 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_216 += tmp_217;
      tmp_215 += (ZU(gI1,j4)) * tmp_216;
   }
   tmp_211 += tmp_215;
   result += (-1) * tmp_210 * tmp_211;
   if (gO1 < 3) {
      std::complex<double> tmp_218;
      std::complex<double> tmp_219;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_219 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_218 += tmp_219;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_218;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_220;
      std::complex<double> tmp_221;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_221 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_220 += tmp_221;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_220;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_222;
      std::complex<double> tmp_223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_223 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_222 += tmp_223;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_222;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_224;
      std::complex<double> tmp_225;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_225 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_224 += tmp_225;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_224;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_226;
      std::complex<double> tmp_227;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_227 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_226 += tmp_227;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_226;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_228;
      std::complex<double> tmp_229;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_229 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_228 += tmp_229;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_228;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_230;
      std::complex<double> tmp_232;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_232 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_230 += tmp_232;
      std::complex<double> tmp_231;
      std::complex<double> tmp_233;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_233 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_231 += tmp_233;
      result += (-1) * tmp_230 * tmp_231;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_234;
   std::complex<double> tmp_235;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_236;
      std::complex<double> tmp_237;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_237 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_236 += tmp_237;
      tmp_235 += (Conj(ZU(gI1,j2))) * tmp_236;
   }
   tmp_234 += tmp_235;
   result += (ZP(gI2,0)) * tmp_234;
   std::complex<double> tmp_238;
   std::complex<double> tmp_239;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_240;
      std::complex<double> tmp_241;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_242;
         std::complex<double> tmp_243;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_243 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_242 += tmp_243;
         tmp_241 += (KroneckerDelta(gO2,3 + j2)) * tmp_242;
      }
      tmp_240 += tmp_241;
      tmp_239 += (Conj(ZU(gI1,3 + j3))) * tmp_240;
   }
   tmp_238 += tmp_239;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_238;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_244;
      std::complex<double> tmp_245;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_245 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_244 += tmp_245;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI2,0)) * tmp_244;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_246;
      std::complex<double> tmp_247;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_248;
         std::complex<double> tmp_249;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_249 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_248 += tmp_249;
         tmp_247 += (Conj(ZU(gI1,j2))) * tmp_248;
      }
      tmp_246 += tmp_247;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_246;
   }
   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_253 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_252 += tmp_253;
      tmp_251 += (Conj(ZU(gI1,j2))) * tmp_252;
   }
   tmp_250 += tmp_251;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_250;
   std::complex<double> tmp_254;
   std::complex<double> tmp_255;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_258;
         std::complex<double> tmp_259;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_259 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_258 += tmp_259;
         tmp_257 += (KroneckerDelta(gO2,3 + j2)) * tmp_258;
      }
      tmp_256 += tmp_257;
      tmp_255 += (Conj(ZU(gI1,3 + j3))) * tmp_256;
   }
   tmp_254 += tmp_255;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_254;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_260;
      std::complex<double> tmp_261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_261 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_260 += tmp_261;
      result += (ZP(gI2,1)) * tmp_260;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_262;
      std::complex<double> tmp_263;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_264;
         std::complex<double> tmp_265;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_265 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_264 += tmp_265;
         tmp_263 += (Conj(ZU(gI1,j2))) * tmp_264;
      }
      tmp_262 += tmp_263;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_262;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_269 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_268 += tmp_269;
      tmp_267 += (Conj(ZD(gI1,j2))) * tmp_268;
   }
   tmp_266 += tmp_267;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_266
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_270;
      std::complex<double> tmp_271;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_271 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_270 += tmp_271;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_270;
   }
   std::complex<double> tmp_272;
   std::complex<double> tmp_273;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_274;
      std::complex<double> tmp_275;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_275 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_274 += tmp_275;
      tmp_273 += (Conj(ZD(gI1,j2))) * tmp_274;
   }
   tmp_272 += tmp_273;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*ZA(gI2,1)) *
      tmp_272;
   if (gO2 < 3) {
      std::complex<double> tmp_276;
      std::complex<double> tmp_277;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_277 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_276 += tmp_277;
      result += (std::complex<double>(0,0.5)*vS*Lambdax*ZA(gI2,1)) * tmp_276
         ;
   }
   std::complex<double> tmp_278;
   std::complex<double> tmp_279;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_280;
      std::complex<double> tmp_281;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_281 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_280 += tmp_281;
      tmp_279 += (Conj(ZD(gI1,j2))) * tmp_280;
   }
   tmp_278 += tmp_279;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*ZA(gI2,2)) *
      tmp_278;
   if (gO2 < 3) {
      std::complex<double> tmp_282;
      std::complex<double> tmp_283;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_283 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_282 += tmp_283;
      result += (std::complex<double>(0,0.5)*vu*Lambdax*ZA(gI2,2)) * tmp_282
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_284;
   std::complex<double> tmp_285;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_285 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_284 += tmp_285;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_284;
   std::complex<double> tmp_286;
   std::complex<double> tmp_287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_288;
      std::complex<double> tmp_289;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_289 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_288 += tmp_289;
      tmp_287 += (Conj(ZD(gI1,j2))) * tmp_288;
   }
   tmp_286 += tmp_287;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_286;
   std::complex<double> tmp_290;
   std::complex<double> tmp_291;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_292;
      std::complex<double> tmp_293;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_294;
         std::complex<double> tmp_295;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_295 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_294 += tmp_295;
         tmp_293 += (KroneckerDelta(gO2,3 + j2)) * tmp_294;
      }
      tmp_292 += tmp_293;
      tmp_291 += (Conj(ZD(gI1,3 + j3))) * tmp_292;
   }
   tmp_290 += tmp_291;
   result += (-(vd*ZH(gI2,0))) * tmp_290;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_296;
      std::complex<double> tmp_297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_297 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_296 += tmp_297;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_296;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_298;
      std::complex<double> tmp_299;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_300;
         std::complex<double> tmp_301;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_301 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_300 += tmp_301;
         tmp_299 += (Conj(ZD(gI1,j2))) * tmp_300;
      }
      tmp_298 += tmp_299;
      result += (-(vd*ZH(gI2,0))) * tmp_298;
   }
   std::complex<double> tmp_302;
   std::complex<double> tmp_303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_303 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_302 += tmp_303;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_302;
   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_306;
      std::complex<double> tmp_307;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_307 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_306 += tmp_307;
      tmp_305 += (Conj(ZD(gI1,j2))) * tmp_306;
   }
   tmp_304 += tmp_305;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,1)) * tmp_304;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_308;
      std::complex<double> tmp_309;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_309 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_308 += tmp_309;
      result += (0.5*vS*Lambdax*ZH(gI2,1)) * tmp_308;
   }
   std::complex<double> tmp_310;
   std::complex<double> tmp_311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_312;
      std::complex<double> tmp_313;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_313 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_312 += tmp_313;
      tmp_311 += (Conj(ZD(gI1,j2))) * tmp_312;
   }
   tmp_310 += tmp_311;
   result += (0.5*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_310;
   if (gO2 < 3) {
      std::complex<double> tmp_314;
      std::complex<double> tmp_315;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_315 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_314 += tmp_315;
      result += (0.5*vu*Lambdax*ZH(gI2,2)) * tmp_314;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_316;
   std::complex<double> tmp_317;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_317 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_316 += tmp_317;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_316;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*PhaseGlu*Conj(ZDL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVGSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 6) {
      result += g3*Conj(ZD(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVPSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_318;
   std::complex<double> tmp_319;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_319 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_318 += tmp_319;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_318;
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZD(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVZSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_320;
   std::complex<double> tmp_321;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_321 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_320 += tmp_321;
   result += (0.2581988897471611*g1*Sin(ThetaW())) * tmp_320;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVWmSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZU(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.1*KroneckerDelta(gO1,gO2)*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_322;
      std::complex<double> tmp_323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_323 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_322 += tmp_323;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_322;
   }
   result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += Conj(Ye(j1,gO2))*ZER(gI2,j1);
      }
      tmp_324 += tmp_325;
      result += (UM(gI1,1)) * tmp_324;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZEL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_326;
      std::complex<double> tmp_327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_327 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_326 += tmp_327;
      result += (ZP(gI1,0)) * tmp_326;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_328;
      std::complex<double> tmp_329;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_330;
         std::complex<double> tmp_331;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_331 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_330 += tmp_331;
         tmp_329 += (Conj(ZE(gI2,j2))) * tmp_330;
      }
      tmp_328 += tmp_329;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_328;
   }
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_332;
      std::complex<double> tmp_333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_333 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_332 += tmp_333;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI1,1)) * tmp_332;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,0)*ZA
      (gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1);
   result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2);
   if (gI1 < 3 && gI2 < 3) {
      result += -0.15*Conj(ZV(gI2,gO2))*Sqr(g1)*ZV(gI1,gO1);
   }
   if (gI1 < 3 && gI2 < 3) {
      result += -0.25*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,0)*ZH
      (gI2,0) - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -0.15*vd*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += -0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += 0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gI1 < 3) {
      result += 0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }

   return result;
}

double CLASSNAME::CpconjUSvFvChiPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZN(gI2,0))*KroneckerDelta(gI1,gO1
         );
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*KroneckerDelta(gI1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_334;
   std::complex<double> tmp_335;
   std::complex<double> tmp_336;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_336 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_335 += tmp_336;
   tmp_334 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_335;
   std::complex<double> tmp_337;
   std::complex<double> tmp_338;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_338 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_337 += tmp_338;
   tmp_334 += (std::complex<double>(0,0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_337;
   std::complex<double> tmp_339;
   std::complex<double> tmp_340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_340 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_339 += tmp_340;
   tmp_334 += (std::complex<double>(0,0.1)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_339;
   result += (std::complex<double>(0,-1)) * tmp_334;

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_341;
   std::complex<double> tmp_342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_342 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_341 += tmp_342;
   result += (-0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_341;
   std::complex<double> tmp_343;
   std::complex<double> tmp_344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_344 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_343 += tmp_344;
   result += (0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_343;
   std::complex<double> tmp_345;
   std::complex<double> tmp_346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_346 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_345 += tmp_346;
   result += (0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_345;
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_347;
      std::complex<double> tmp_349;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_349 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_347 += tmp_349;
      std::complex<double> tmp_348;
      std::complex<double> tmp_350;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_350 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_348 += tmp_350;
      result += (-1) * tmp_347 * tmp_348;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_351;
   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_353 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_352 += tmp_353;
   tmp_351 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_352;
   std::complex<double> tmp_354;
   std::complex<double> tmp_355;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_355 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_354 += tmp_355;
   tmp_351 += (std::complex<double>(0,-0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_354;
   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_357 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_356 += tmp_357;
   tmp_351 += (std::complex<double>(0,-0.2)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_356;
   result += (std::complex<double>(0,-1)) * tmp_351;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjVWmSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_359 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_358 += tmp_359;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_358;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += -0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()));
   }

   return result;
}

double CLASSNAME::CpUSuconjUSuconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_361 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_360 += tmp_361;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_360;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_362;
      std::complex<double> tmp_363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_363 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_362 += tmp_363;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_362;
   }
   std::complex<double> tmp_364;
   std::complex<double> tmp_365;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_365 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_364 += tmp_365;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_364;
   std::complex<double> tmp_366;
   std::complex<double> tmp_367;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_368;
      std::complex<double> tmp_369;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_370;
         std::complex<double> tmp_371;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_371 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_370 += tmp_371;
         tmp_369 += (KroneckerDelta(gO2,3 + j2)) * tmp_370;
      }
      tmp_368 += tmp_369;
      tmp_367 += (KroneckerDelta(gO1,3 + j3)) * tmp_368;
   }
   tmp_366 += tmp_367;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_366;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_372;
      std::complex<double> tmp_373;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_373 += Conj(Yd(j1,gO2))*ZDR(gI2,j1);
      }
      tmp_372 += tmp_373;
      result += (UM(gI1,1)) * tmp_372;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_374;
   std::complex<double> tmp_375;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_376;
      std::complex<double> tmp_377;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_377 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_376 += tmp_377;
      tmp_375 += (Conj(ZDL(gI2,j2))) * tmp_376;
   }
   tmp_374 += tmp_375;
   result += (Conj(UP(gI1,1))) * tmp_374;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_378;
   std::complex<double> tmp_379;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_380;
      std::complex<double> tmp_381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_381 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_380 += tmp_381;
      tmp_379 += (Conj(ZD(gI2,j2))) * tmp_380;
   }
   tmp_378 += tmp_379;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI1,0)) * tmp_378;
   std::complex<double> tmp_382;
   std::complex<double> tmp_383;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_384;
      std::complex<double> tmp_385;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_386;
         std::complex<double> tmp_387;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_387 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_386 += tmp_387;
         tmp_385 += (KroneckerDelta(gO2,3 + j2)) * tmp_386;
      }
      tmp_384 += tmp_385;
      tmp_383 += (Conj(ZD(gI2,3 + j3))) * tmp_384;
   }
   tmp_382 += tmp_383;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_382;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_388;
      std::complex<double> tmp_389;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_389 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_388 += tmp_389;
      result += (ZP(gI1,0)) * tmp_388;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_390;
      std::complex<double> tmp_391;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_392;
         std::complex<double> tmp_393;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_393 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_392 += tmp_393;
         tmp_391 += (Conj(ZD(gI2,j2))) * tmp_392;
      }
      tmp_390 += tmp_391;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_390;
   }
   std::complex<double> tmp_394;
   std::complex<double> tmp_395;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_396;
      std::complex<double> tmp_397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_397 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_396 += tmp_397;
      tmp_395 += (Conj(ZD(gI2,j2))) * tmp_396;
   }
   tmp_394 += tmp_395;
   result += (ZP(gI1,1)) * tmp_394;
   std::complex<double> tmp_398;
   std::complex<double> tmp_399;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_400;
      std::complex<double> tmp_401;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_402;
         std::complex<double> tmp_403;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_403 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_402 += tmp_403;
         tmp_401 += (KroneckerDelta(gO2,3 + j2)) * tmp_402;
      }
      tmp_400 += tmp_401;
      tmp_399 += (Conj(ZD(gI2,3 + j3))) * tmp_400;
   }
   tmp_398 += tmp_399;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_398;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_404;
      std::complex<double> tmp_405;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_405 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_404 += tmp_405;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI1,1)) * tmp_404;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_406;
      std::complex<double> tmp_407;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_408;
         std::complex<double> tmp_409;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_409 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_408 += tmp_409;
         tmp_407 += (Conj(ZD(gI2,j2))) * tmp_408;
      }
      tmp_406 += tmp_407;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_406;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_410;
   std::complex<double> tmp_411;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_411 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_410 += tmp_411;
   result += (-0.2*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_410;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_412;
      std::complex<double> tmp_413;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_413 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_412 += tmp_413;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,0)) * tmp_412;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_414;
      std::complex<double> tmp_415;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_415 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_414 += tmp_415;
      result += (-0.5*Lambdax*ZA(gI1,2)*ZA(gI2,0)) * tmp_414;
   }
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_417 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_416 += tmp_417;
   result += (0.2*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_416;
   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_420;
      std::complex<double> tmp_421;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_422;
         std::complex<double> tmp_423;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_423 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_422 += tmp_423;
         tmp_421 += (KroneckerDelta(gO2,3 + j2)) * tmp_422;
      }
      tmp_420 += tmp_421;
      tmp_419 += (KroneckerDelta(gO1,3 + j3)) * tmp_420;
   }
   tmp_418 += tmp_419;
   result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_418;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_424;
      std::complex<double> tmp_425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_425 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_424 += tmp_425;
      result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_424;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_426;
      std::complex<double> tmp_427;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_427 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_426 += tmp_427;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,0)*ZA(gI2,2)) * tmp_426;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_429 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_428 += tmp_429;
      result += (-0.5*Lambdax*ZA(gI1,0)*ZA(gI2,2)) * tmp_428;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_430;
   std::complex<double> tmp_431;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_431 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_430 += tmp_431;
   result += (-0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_430;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_433 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_432 += tmp_433;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_432;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_434;
      std::complex<double> tmp_435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_435 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_434 += tmp_435;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,0)) * tmp_434;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_436;
      std::complex<double> tmp_437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_437 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_436 += tmp_437;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,0)) * tmp_436;
   }
   std::complex<double> tmp_438;
   std::complex<double> tmp_439;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_439 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_438 += tmp_439;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_438;
   std::complex<double> tmp_440;
   std::complex<double> tmp_441;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_442;
      std::complex<double> tmp_443;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_444;
         std::complex<double> tmp_445;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_445 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_444 += tmp_445;
         tmp_443 += (KroneckerDelta(gO2,3 + j2)) * tmp_444;
      }
      tmp_442 += tmp_443;
      tmp_441 += (KroneckerDelta(gO1,3 + j3)) * tmp_442;
   }
   tmp_440 += tmp_441;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_440;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_446;
      std::complex<double> tmp_447;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_447 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_446 += tmp_447;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_446;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_448;
      std::complex<double> tmp_449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_449 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_448 += tmp_449;
      result += (0.5*Conj(Lambdax)*ZH(gI1,0)*ZH(gI2,2)) * tmp_448;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_450;
      std::complex<double> tmp_451;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_451 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_450 += tmp_451;
      result += (0.5*Lambdax*ZH(gI1,0)*ZH(gI2,2)) * tmp_450;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_452;
   std::complex<double> tmp_453;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_453 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_452 += tmp_453;
   result += (0.7302967433402214*g1*ZN(gI2,0)) * tmp_452;
   if (gO2 < 3) {
      std::complex<double> tmp_454;
      std::complex<double> tmp_455;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_455 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_454 += tmp_455;
      result += (-ZN(gI2,3)) * tmp_454;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_456;
   std::complex<double> tmp_457;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_458;
      std::complex<double> tmp_459;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_459 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_458 += tmp_459;
      tmp_457 += (Conj(ZUL(gI1,j2))) * tmp_458;
   }
   tmp_456 += tmp_457;
   result += (-Conj(ZN(gI2,3))) * tmp_456;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZN(gI2,0))*Conj(ZUL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*Conj(ZUL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_460;
   std::complex<double> tmp_462;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_462 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_460 += tmp_462;
   std::complex<double> tmp_461;
   std::complex<double> tmp_463;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_463 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_461 += tmp_463;
   result += (0.1*Sqr(g1)) * tmp_460 * tmp_461;
   std::complex<double> tmp_464;
   std::complex<double> tmp_466;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_466 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_464 += tmp_466;
   std::complex<double> tmp_465;
   std::complex<double> tmp_467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_467 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_465 += tmp_467;
   result += (0.2*Sqr(g1)) * tmp_464 * tmp_465;
   std::complex<double> tmp_468;
   std::complex<double> tmp_470;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_470 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_468 += tmp_470;
   std::complex<double> tmp_469;
   std::complex<double> tmp_471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_471 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_469 += tmp_471;
   result += (0.1*Sqr(g1)) * tmp_468 * tmp_469;
   std::complex<double> tmp_472;
   std::complex<double> tmp_474;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_474 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_472 += tmp_474;
   std::complex<double> tmp_473;
   std::complex<double> tmp_475;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_475 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_473 += tmp_475;
   result += (0.2*Sqr(g1)) * tmp_472 * tmp_473;
   std::complex<double> tmp_476;
   std::complex<double> tmp_478;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_479;
      std::complex<double> tmp_480;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_480 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_479 += tmp_480;
      tmp_478 += (Conj(ZD(gI2,j2))) * tmp_479;
   }
   tmp_476 += tmp_478;
   std::complex<double> tmp_477;
   std::complex<double> tmp_481;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_482;
      std::complex<double> tmp_483;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_483 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_482 += tmp_483;
      tmp_481 += (ZD(gI1,j4)) * tmp_482;
   }
   tmp_477 += tmp_481;
   result += (-1) * tmp_476 * tmp_477;
   if (gO1 < 3) {
      std::complex<double> tmp_484;
      std::complex<double> tmp_485;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_485 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_484 += tmp_485;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_484;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_486;
      std::complex<double> tmp_487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_487 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_486 += tmp_487;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_486;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_488;
      std::complex<double> tmp_489;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_489 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_488 += tmp_489;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_488;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_490;
      std::complex<double> tmp_491;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_491 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_490 += tmp_491;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_490;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_492;
      std::complex<double> tmp_493;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_493 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_492 += tmp_493;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_492;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_494;
      std::complex<double> tmp_495;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_495 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_494 += tmp_495;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_494;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_496;
      std::complex<double> tmp_498;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_498 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_496 += tmp_498;
      std::complex<double> tmp_497;
      std::complex<double> tmp_499;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_499 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_497 += tmp_499;
      result += (-1) * tmp_496 * tmp_497;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_500;
   std::complex<double> tmp_502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_502 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_500 += tmp_502;
   std::complex<double> tmp_501;
   std::complex<double> tmp_503;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_503 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_501 += tmp_503;
   result += (-0.1*Sqr(g1)) * tmp_500 * tmp_501;
   std::complex<double> tmp_504;
   std::complex<double> tmp_506;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_506 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_504 += tmp_506;
   std::complex<double> tmp_505;
   std::complex<double> tmp_507;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_507 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_505 += tmp_507;
   result += (0.2*Sqr(g1)) * tmp_504 * tmp_505;
   std::complex<double> tmp_508;
   std::complex<double> tmp_510;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_510 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_508 += tmp_510;
   std::complex<double> tmp_509;
   std::complex<double> tmp_511;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_511 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_509 += tmp_511;
   result += (-0.1*Sqr(g1)) * tmp_508 * tmp_509;
   std::complex<double> tmp_512;
   std::complex<double> tmp_514;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_514 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_512 += tmp_514;
   std::complex<double> tmp_513;
   std::complex<double> tmp_515;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_515 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_513 += tmp_515;
   result += (0.2*Sqr(g1)) * tmp_512 * tmp_513;
   if (gO1 < 3) {
      std::complex<double> tmp_516;
      std::complex<double> tmp_517;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_517 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_516 += tmp_517;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_516;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_518;
      std::complex<double> tmp_519;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_519 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_518 += tmp_519;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_518;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_520;
      std::complex<double> tmp_521;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_521 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_520 += tmp_521;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_520;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_522;
      std::complex<double> tmp_523;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_523 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_522 += tmp_523;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_522;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_525 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_524 += tmp_525;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_524;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_526;
      std::complex<double> tmp_527;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_527 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_526 += tmp_527;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_526;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_528;
   std::complex<double> tmp_530;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_530 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_528 += tmp_530;
   std::complex<double> tmp_529;
   std::complex<double> tmp_531;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_531 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_529 += tmp_531;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_528 * tmp_529;
   std::complex<double> tmp_532;
   std::complex<double> tmp_534;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_534 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_532 += tmp_534;
   std::complex<double> tmp_533;
   std::complex<double> tmp_535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_535 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_533 += tmp_535;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_532 * tmp_533;
   std::complex<double> tmp_536;
   std::complex<double> tmp_538;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_538 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_536 += tmp_538;
   std::complex<double> tmp_537;
   std::complex<double> tmp_539;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_539 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_537 += tmp_539;
   result += (0.1*Sqr(g1)) * tmp_536 * tmp_537;
   std::complex<double> tmp_540;
   std::complex<double> tmp_542;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_542 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_540 += tmp_542;
   std::complex<double> tmp_541;
   std::complex<double> tmp_543;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_543 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_541 += tmp_543;
   result += (-0.4*Sqr(g1)) * tmp_540 * tmp_541;
   std::complex<double> tmp_544;
   std::complex<double> tmp_546;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_546 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_544 += tmp_546;
   std::complex<double> tmp_545;
   std::complex<double> tmp_547;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_547 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_545 += tmp_547;
   result += (0.1*Sqr(g1)) * tmp_544 * tmp_545;
   std::complex<double> tmp_548;
   std::complex<double> tmp_550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_550 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_548 += tmp_550;
   std::complex<double> tmp_549;
   std::complex<double> tmp_551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_551 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_549 += tmp_551;
   result += (-0.4*Sqr(g1)) * tmp_548 * tmp_549;
   std::complex<double> tmp_552;
   std::complex<double> tmp_554;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_554 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_552 += tmp_554;
   std::complex<double> tmp_553;
   std::complex<double> tmp_555;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_555 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_553 += tmp_555;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_552 * tmp_553;
   std::complex<double> tmp_556;
   std::complex<double> tmp_558;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_558 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_556 += tmp_558;
   std::complex<double> tmp_557;
   std::complex<double> tmp_559;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_559 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_557 += tmp_559;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_556 * tmp_557;
   std::complex<double> tmp_560;
   std::complex<double> tmp_562;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_563;
      std::complex<double> tmp_564;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_564 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_563 += tmp_564;
      tmp_562 += (Conj(ZU(gI2,j2))) * tmp_563;
   }
   tmp_560 += tmp_562;
   std::complex<double> tmp_561;
   std::complex<double> tmp_565;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_566;
      std::complex<double> tmp_567;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_567 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_566 += tmp_567;
      tmp_565 += (ZU(gI1,j4)) * tmp_566;
   }
   tmp_561 += tmp_565;
   result += (-1) * tmp_560 * tmp_561;
   if (gO1 < 3) {
      std::complex<double> tmp_568;
      std::complex<double> tmp_569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_569 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_568 += tmp_569;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_568;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_570;
      std::complex<double> tmp_571;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_571 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_570 += tmp_571;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_570;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_572;
      std::complex<double> tmp_573;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_573 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_572 += tmp_573;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_572;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_574;
      std::complex<double> tmp_575;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_575 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_574 += tmp_575;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_574;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_576;
      std::complex<double> tmp_577;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_577 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_576 += tmp_577;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_576;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_578;
      std::complex<double> tmp_579;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_579 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_578 += tmp_579;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_578;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_580;
      std::complex<double> tmp_582;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_582 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_580 += tmp_582;
      std::complex<double> tmp_581;
      std::complex<double> tmp_583;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_584;
         std::complex<double> tmp_585;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_585 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_584 += tmp_585;
         tmp_583 += (ZU(gI1,j4)) * tmp_584;
      }
      tmp_581 += tmp_583;
      result += (-3) * tmp_580 * tmp_581;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_586;
      std::complex<double> tmp_587;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_587 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_586 += tmp_587;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_586;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_588;
      std::complex<double> tmp_589;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_589 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_588 += tmp_589;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_588;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_590;
      std::complex<double> tmp_591;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_591 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_590 += tmp_591;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_590;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_592;
      std::complex<double> tmp_593;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_593 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_592 += tmp_593;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_592;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_594;
      std::complex<double> tmp_596;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_597;
         std::complex<double> tmp_598;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_598 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_597 += tmp_598;
         tmp_596 += (Conj(ZU(gI2,j2))) * tmp_597;
      }
      tmp_594 += tmp_596;
      std::complex<double> tmp_595;
      std::complex<double> tmp_599;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_599 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_595 += tmp_599;
      result += (-3) * tmp_594 * tmp_595;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_600;
      std::complex<double> tmp_602;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_602 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_600 += tmp_602;
      std::complex<double> tmp_601;
      std::complex<double> tmp_603;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_603 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_601 += tmp_603;
      result += (-1) * tmp_600 * tmp_601;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_604;
      std::complex<double> tmp_605;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_605 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_604 += tmp_605;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_604;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_606;
      std::complex<double> tmp_607;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_607 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_606 += tmp_607;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_606;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_608;
      std::complex<double> tmp_609;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_609 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_608 += tmp_609;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_608;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_610;
      std::complex<double> tmp_611;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_611 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_610 += tmp_611;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_610;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.016666666666666666*Conj(ZU(gI2,gO2))*Sqr(g1)*ZU(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -1.3333333333333333*Conj(ZU(gI2,gO2))*Sqr(g3)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_612;
   std::complex<double> tmp_613;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_615 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_614 += tmp_615;
      tmp_613 += (Conj(ZU(gI1,j2))) * tmp_614;
   }
   tmp_612 += tmp_613;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*ZA(gI2,0)) *
      tmp_612;
   if (gO2 < 3) {
      std::complex<double> tmp_616;
      std::complex<double> tmp_617;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_617 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_616 += tmp_617;
      result += (std::complex<double>(0,0.5)*vS*Lambdax*ZA(gI2,0)) * tmp_616
         ;
   }
   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_620;
      std::complex<double> tmp_621;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_621 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_620 += tmp_621;
      tmp_619 += (Conj(ZU(gI1,j2))) * tmp_620;
   }
   tmp_618 += tmp_619;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) * tmp_618
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_622;
      std::complex<double> tmp_623;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_623 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_622 += tmp_623;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_622;
   }
   std::complex<double> tmp_624;
   std::complex<double> tmp_625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_626 += tmp_627;
      tmp_625 += (Conj(ZU(gI1,j2))) * tmp_626;
   }
   tmp_624 += tmp_625;
   result += (std::complex<double>(0,-0.5)*vd*Conj(Lambdax)*ZA(gI2,2)) *
      tmp_624;
   if (gO2 < 3) {
      std::complex<double> tmp_628;
      std::complex<double> tmp_629;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_629 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_628 += tmp_629;
      result += (std::complex<double>(0,0.5)*vd*Lambdax*ZA(gI2,2)) * tmp_628
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_630;
   std::complex<double> tmp_631;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_631 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_630 += tmp_631;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_630;
   std::complex<double> tmp_632;
   std::complex<double> tmp_633;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_634;
      std::complex<double> tmp_635;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_635 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_634 += tmp_635;
      tmp_633 += (Conj(ZU(gI1,j2))) * tmp_634;
   }
   tmp_632 += tmp_633;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,0)) * tmp_632;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_636;
      std::complex<double> tmp_637;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_637 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_636 += tmp_637;
      result += (0.5*vS*Lambdax*ZH(gI2,0)) * tmp_636;
   }
   std::complex<double> tmp_638;
   std::complex<double> tmp_639;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_639 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_638 += tmp_639;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_638;
   std::complex<double> tmp_640;
   std::complex<double> tmp_641;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_642;
      std::complex<double> tmp_643;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_643 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_642 += tmp_643;
      tmp_641 += (Conj(ZU(gI1,j2))) * tmp_642;
   }
   tmp_640 += tmp_641;
   result += (-0.7071067811865475*ZH(gI2,1)) * tmp_640;
   std::complex<double> tmp_644;
   std::complex<double> tmp_645;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_646;
      std::complex<double> tmp_647;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_648;
         std::complex<double> tmp_649;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_649 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_648 += tmp_649;
         tmp_647 += (KroneckerDelta(gO2,3 + j2)) * tmp_648;
      }
      tmp_646 += tmp_647;
      tmp_645 += (Conj(ZU(gI1,3 + j3))) * tmp_646;
   }
   tmp_644 += tmp_645;
   result += (-(vu*ZH(gI2,1))) * tmp_644;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_650;
      std::complex<double> tmp_651;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_651 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_650 += tmp_651;
      result += (-0.7071067811865475*ZH(gI2,1)) * tmp_650;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_652;
      std::complex<double> tmp_653;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_654;
         std::complex<double> tmp_655;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_655 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_654 += tmp_655;
         tmp_653 += (Conj(ZU(gI1,j2))) * tmp_654;
      }
      tmp_652 += tmp_653;
      result += (-(vu*ZH(gI2,1))) * tmp_652;
   }
   std::complex<double> tmp_656;
   std::complex<double> tmp_657;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_658;
      std::complex<double> tmp_659;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_659 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_658 += tmp_659;
      tmp_657 += (Conj(ZU(gI1,j2))) * tmp_658;
   }
   tmp_656 += tmp_657;
   result += (0.5*vd*Conj(Lambdax)*ZH(gI2,2)) * tmp_656;
   if (gO2 < 3) {
      std::complex<double> tmp_660;
      std::complex<double> tmp_661;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_661 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_660 += tmp_661;
      result += (0.5*vd*Lambdax*ZH(gI2,2)) * tmp_660;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_662;
   std::complex<double> tmp_663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_663 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_662 += tmp_663;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_662;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*PhaseGlu*Conj(ZUL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjVWmSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVGSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 6) {
      result += g3*Conj(ZU(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVPSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_664;
   std::complex<double> tmp_665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_665 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_664 += tmp_665;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_664;
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZU(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVZSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_666;
   std::complex<double> tmp_667;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_667 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_666 += tmp_667;
   result += (-0.5163977794943222*g1*Sin(ThetaW())) * tmp_666;
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZU(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_668;
   std::complex<double> tmp_669;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_669 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_668 += tmp_669;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_668;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += -0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()));
   }

   return result;
}

double CLASSNAME::CpUSeconjUSeconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_670;
   std::complex<double> tmp_671;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_671 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_670 += tmp_671;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_670;
   std::complex<double> tmp_672;
   std::complex<double> tmp_673;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_674;
      std::complex<double> tmp_675;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_676;
         std::complex<double> tmp_677;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_677 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_676 += tmp_677;
         tmp_675 += (KroneckerDelta(gO2,3 + j2)) * tmp_676;
      }
      tmp_674 += tmp_675;
      tmp_673 += (KroneckerDelta(gO1,3 + j3)) * tmp_674;
   }
   tmp_672 += tmp_673;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_672;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_678;
   std::complex<double> tmp_679;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_679 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_678 += tmp_679;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_678;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_680;
   std::complex<double> tmp_681;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_681 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_680 += tmp_681;
   result += (0.3*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_680;
   std::complex<double> tmp_682;
   std::complex<double> tmp_683;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_684;
      std::complex<double> tmp_685;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_686;
         std::complex<double> tmp_687;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_687 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_686 += tmp_687;
         tmp_685 += (KroneckerDelta(gO2,3 + j2)) * tmp_686;
      }
      tmp_684 += tmp_685;
      tmp_683 += (KroneckerDelta(gO1,3 + j3)) * tmp_684;
   }
   tmp_682 += tmp_683;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_682;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_689 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_688 += tmp_689;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_688;
   }
   std::complex<double> tmp_690;
   std::complex<double> tmp_691;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_691 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_690 += tmp_691;
   result += (-0.3*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_690;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_692;
      std::complex<double> tmp_693;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_693 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_692 += tmp_693;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,1)) * tmp_692;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_694;
      std::complex<double> tmp_695;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_695 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_694 += tmp_695;
      result += (-0.5*Lambdax*ZA(gI1,2)*ZA(gI2,1)) * tmp_694;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_696;
      std::complex<double> tmp_697;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_697 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_696 += tmp_697;
      result += (-0.5*Conj(Lambdax)*ZA(gI1,1)*ZA(gI2,2)) * tmp_696;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_698;
      std::complex<double> tmp_699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_699 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_698 += tmp_699;
      result += (-0.5*Lambdax*ZA(gI1,1)*ZA(gI2,2)) * tmp_698;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_700;
   std::complex<double> tmp_701;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_701 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_700 += tmp_701;
   result += (0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_700;
   std::complex<double> tmp_702;
   std::complex<double> tmp_704;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_705;
      std::complex<double> tmp_706;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_706 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_705 += tmp_706;
      tmp_704 += (Conj(ZV(gI2,j2))) * tmp_705;
   }
   tmp_702 += tmp_704;
   std::complex<double> tmp_703;
   std::complex<double> tmp_707;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_708;
      std::complex<double> tmp_709;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_709 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_708 += tmp_709;
      tmp_707 += (ZV(gI1,j4)) * tmp_708;
   }
   tmp_703 += tmp_707;
   result += (-1) * tmp_702 * tmp_703;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1
         );
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_710;
   std::complex<double> tmp_711;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_711 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_710 += tmp_711;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_710;
   std::complex<double> tmp_712;
   std::complex<double> tmp_713;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_714;
      std::complex<double> tmp_715;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_716;
         std::complex<double> tmp_717;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_717 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_716 += tmp_717;
         tmp_715 += (KroneckerDelta(gO2,3 + j2)) * tmp_716;
      }
      tmp_714 += tmp_715;
      tmp_713 += (KroneckerDelta(gO1,3 + j3)) * tmp_714;
   }
   tmp_712 += tmp_713;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_712;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_718;
      std::complex<double> tmp_719;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_719 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_718 += tmp_719;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_718;
   }
   std::complex<double> tmp_720;
   std::complex<double> tmp_721;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_721 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_720 += tmp_721;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_720;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_722;
      std::complex<double> tmp_723;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_723 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_722 += tmp_723;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) * tmp_722;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_724;
      std::complex<double> tmp_725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_725 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_724 += tmp_725;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_724;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_726;
      std::complex<double> tmp_727;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_727 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_726 += tmp_727;
      result += (0.5*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) * tmp_726;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_728;
      std::complex<double> tmp_729;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_729 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_728 += tmp_729;
      result += (0.5*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_728;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_730;
   std::complex<double> tmp_731;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_732;
      std::complex<double> tmp_733;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_733 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_732 += tmp_733;
      tmp_731 += (Conj(ZV(gI1,j2))) * tmp_732;
   }
   tmp_730 += tmp_731;
   result += (ZP(gI2,0)) * tmp_730;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_734;
      std::complex<double> tmp_735;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_736;
         std::complex<double> tmp_737;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_737 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_736 += tmp_737;
         tmp_735 += (Conj(ZV(gI1,j2))) * tmp_736;
      }
      tmp_734 += tmp_735;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_734;
   }
   std::complex<double> tmp_738;
   std::complex<double> tmp_739;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_740;
      std::complex<double> tmp_741;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_741 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_740 += tmp_741;
      tmp_739 += (Conj(ZV(gI1,j2))) * tmp_740;
   }
   tmp_738 += tmp_739;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_738;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }

   return result;
}

double CLASSNAME::CpconjUSeFvChaPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_742;
   std::complex<double> tmp_743;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_743 += KroneckerDelta(gO1,3 + j1)*Ye(j1,gI1);
   }
   tmp_742 += tmp_743;
   result += (Conj(UM(gI2,1))) * tmp_742;
   if (gI1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*KroneckerDelta(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_744;
   std::complex<double> tmp_745;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_745 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_744 += tmp_745;
   result += (-1.0954451150103321*g1*ZN(gI2,0)) * tmp_744;
   if (gO2 < 3) {
      std::complex<double> tmp_746;
      std::complex<double> tmp_747;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_747 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_746 += tmp_747;
      result += (-ZN(gI2,2)) * tmp_746;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_748;
   std::complex<double> tmp_749;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_750;
      std::complex<double> tmp_751;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_751 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_750 += tmp_751;
      tmp_749 += (Conj(ZEL(gI1,j2))) * tmp_750;
   }
   tmp_748 += tmp_749;
   result += (-Conj(ZN(gI2,2))) * tmp_748;
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_752;
   std::complex<double> tmp_754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_754 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_752 += tmp_754;
   std::complex<double> tmp_753;
   std::complex<double> tmp_755;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_755 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_753 += tmp_755;
   result += (-0.05*Sqr(g1)) * tmp_752 * tmp_753;
   std::complex<double> tmp_756;
   std::complex<double> tmp_758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_758 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_756 += tmp_758;
   std::complex<double> tmp_757;
   std::complex<double> tmp_759;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_759 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_757 += tmp_759;
   result += (-0.1*Sqr(g1)) * tmp_756 * tmp_757;
   std::complex<double> tmp_760;
   std::complex<double> tmp_762;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_762 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_760 += tmp_762;
   std::complex<double> tmp_761;
   std::complex<double> tmp_763;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_763 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_761 += tmp_763;
   result += (-0.05*Sqr(g1)) * tmp_760 * tmp_761;
   std::complex<double> tmp_764;
   std::complex<double> tmp_766;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_766 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_764 += tmp_766;
   std::complex<double> tmp_765;
   std::complex<double> tmp_767;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_767 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_765 += tmp_767;
   result += (-0.1*Sqr(g1)) * tmp_764 * tmp_765;
   if (gO1 < 3) {
      std::complex<double> tmp_768;
      std::complex<double> tmp_769;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_769 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_768 += tmp_769;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_768;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_770;
      std::complex<double> tmp_771;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_771 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_770 += tmp_771;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_770;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_772;
      std::complex<double> tmp_773;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_773 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_772 += tmp_773;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_772;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_774;
      std::complex<double> tmp_775;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_775 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_774 += tmp_775;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_774;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_776;
      std::complex<double> tmp_777;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_777 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_776 += tmp_777;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_776;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_778;
      std::complex<double> tmp_779;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_779 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_778 += tmp_779;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_778;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_780;
      std::complex<double> tmp_782;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_782 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_780 += tmp_782;
      std::complex<double> tmp_781;
      std::complex<double> tmp_783;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_784;
         std::complex<double> tmp_785;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_785 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_784 += tmp_785;
         tmp_783 += (ZD(gI1,j4)) * tmp_784;
      }
      tmp_781 += tmp_783;
      result += (-1) * tmp_780 * tmp_781;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_786;
      std::complex<double> tmp_788;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_789;
         std::complex<double> tmp_790;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_790 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_789 += tmp_790;
         tmp_788 += (Conj(ZD(gI2,j2))) * tmp_789;
      }
      tmp_786 += tmp_788;
      std::complex<double> tmp_787;
      std::complex<double> tmp_791;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_791 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_787 += tmp_791;
      result += (-1) * tmp_786 * tmp_787;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_792;
   std::complex<double> tmp_794;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_794 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_792 += tmp_794;
   std::complex<double> tmp_793;
   std::complex<double> tmp_795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_795 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_793 += tmp_795;
   result += (-0.3*Sqr(g1)) * tmp_792 * tmp_793;
   std::complex<double> tmp_796;
   std::complex<double> tmp_798;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_798 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_796 += tmp_798;
   std::complex<double> tmp_797;
   std::complex<double> tmp_799;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_799 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_797 += tmp_799;
   result += (0.15*Sqr(g1)) * tmp_796 * tmp_797;
   std::complex<double> tmp_800;
   std::complex<double> tmp_802;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_802 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_800 += tmp_802;
   std::complex<double> tmp_801;
   std::complex<double> tmp_803;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_803 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_801 += tmp_803;
   result += (-0.3*Sqr(g1)) * tmp_800 * tmp_801;
   std::complex<double> tmp_804;
   std::complex<double> tmp_806;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_806 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_804 += tmp_806;
   std::complex<double> tmp_805;
   std::complex<double> tmp_807;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_807 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_805 += tmp_807;
   result += (0.15*Sqr(g1)) * tmp_804 * tmp_805;
   std::complex<double> tmp_808;
   std::complex<double> tmp_810;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_810 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_808 += tmp_810;
   std::complex<double> tmp_809;
   std::complex<double> tmp_811;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_811 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_809 += tmp_811;
   result += (-0.3*Sqr(g1)) * tmp_808 * tmp_809;
   std::complex<double> tmp_812;
   std::complex<double> tmp_814;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_814 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_812 += tmp_814;
   std::complex<double> tmp_813;
   std::complex<double> tmp_815;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_815 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_813 += tmp_815;
   result += (-0.3*Sqr(g1)) * tmp_812 * tmp_813;
   std::complex<double> tmp_816;
   std::complex<double> tmp_818;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_819;
      std::complex<double> tmp_820;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_820 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_819 += tmp_820;
      tmp_818 += (Conj(ZE(gI2,j2))) * tmp_819;
   }
   tmp_816 += tmp_818;
   std::complex<double> tmp_817;
   std::complex<double> tmp_821;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_823;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_823 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_822 += tmp_823;
      tmp_821 += (ZE(gI1,j4)) * tmp_822;
   }
   tmp_817 += tmp_821;
   result += (-1) * tmp_816 * tmp_817;
   if (gO1 < 3) {
      std::complex<double> tmp_824;
      std::complex<double> tmp_825;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_825 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_824 += tmp_825;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_824;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_826;
      std::complex<double> tmp_827;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_827 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_826 += tmp_827;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_826;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_828;
      std::complex<double> tmp_829;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_829 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_828 += tmp_829;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_828;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_830;
      std::complex<double> tmp_831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_831 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_830 += tmp_831;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_830;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_832;
      std::complex<double> tmp_833;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_833 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_832 += tmp_833;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_832;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_834;
      std::complex<double> tmp_835;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_835 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_834 += tmp_835;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_834;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_836;
      std::complex<double> tmp_838;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_838 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_836 += tmp_838;
      std::complex<double> tmp_837;
      std::complex<double> tmp_839;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_840;
         std::complex<double> tmp_841;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_841 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_840 += tmp_841;
         tmp_839 += (ZE(gI1,j4)) * tmp_840;
      }
      tmp_837 += tmp_839;
      result += (-1) * tmp_836 * tmp_837;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_842;
      std::complex<double> tmp_843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_843 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_842 += tmp_843;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_842;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_844;
      std::complex<double> tmp_845;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_845 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_844 += tmp_845;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_844;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_846;
      std::complex<double> tmp_848;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_849;
         std::complex<double> tmp_850;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_850 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_849 += tmp_850;
         tmp_848 += (Conj(ZE(gI2,j2))) * tmp_849;
      }
      tmp_846 += tmp_848;
      std::complex<double> tmp_847;
      std::complex<double> tmp_851;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_851 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_847 += tmp_851;
      result += (-1) * tmp_846 * tmp_847;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_852;
      std::complex<double> tmp_854;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_854 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_852 += tmp_854;
      std::complex<double> tmp_853;
      std::complex<double> tmp_855;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_855 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_853 += tmp_855;
      result += (-1) * tmp_852 * tmp_853;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_856;
      std::complex<double> tmp_857;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_857 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_856 += tmp_857;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_856;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_858;
      std::complex<double> tmp_859;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_859 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_858 += tmp_859;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_858;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.15*Conj(ZE(gI2,gO2))*Sqr(g1)*ZE(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_860;
   std::complex<double> tmp_862;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_862 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_860 += tmp_862;
   std::complex<double> tmp_861;
   std::complex<double> tmp_863;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_863 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_861 += tmp_863;
   result += (-0.05*Sqr(g1)) * tmp_860 * tmp_861;
   std::complex<double> tmp_864;
   std::complex<double> tmp_866;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_866 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_864 += tmp_866;
   std::complex<double> tmp_865;
   std::complex<double> tmp_867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_867 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_865 += tmp_867;
   result += (0.2*Sqr(g1)) * tmp_864 * tmp_865;
   std::complex<double> tmp_868;
   std::complex<double> tmp_870;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_870 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_868 += tmp_870;
   std::complex<double> tmp_869;
   std::complex<double> tmp_871;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_871 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_869 += tmp_871;
   result += (-0.05*Sqr(g1)) * tmp_868 * tmp_869;
   std::complex<double> tmp_872;
   std::complex<double> tmp_874;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_874 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_872 += tmp_874;
   std::complex<double> tmp_873;
   std::complex<double> tmp_875;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_875 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_873 += tmp_875;
   result += (0.2*Sqr(g1)) * tmp_872 * tmp_873;
   if (gO1 < 3) {
      std::complex<double> tmp_876;
      std::complex<double> tmp_877;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_877 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_876 += tmp_877;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_876;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_879 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_878 += tmp_879;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_878;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_880;
      std::complex<double> tmp_881;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_881 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_880 += tmp_881;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_880;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_882;
      std::complex<double> tmp_883;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_883 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_882 += tmp_883;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_882;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_884;
      std::complex<double> tmp_885;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_885 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_884 += tmp_885;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_884;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_886;
      std::complex<double> tmp_887;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_887 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_886 += tmp_887;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_886;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_888;
   std::complex<double> tmp_889;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_890;
      std::complex<double> tmp_891;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_891 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_890 += tmp_891;
      tmp_889 += (Conj(ZE(gI1,j2))) * tmp_890;
   }
   tmp_888 += tmp_889;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_888
      ;
   if (gO2 < 3) {
      std::complex<double> tmp_892;
      std::complex<double> tmp_893;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_893 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_892 += tmp_893;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_892;
   }
   std::complex<double> tmp_894;
   std::complex<double> tmp_895;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_896;
      std::complex<double> tmp_897;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_897 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_896 += tmp_897;
      tmp_895 += (Conj(ZE(gI1,j2))) * tmp_896;
   }
   tmp_894 += tmp_895;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*ZA(gI2,1)) *
      tmp_894;
   if (gO2 < 3) {
      std::complex<double> tmp_898;
      std::complex<double> tmp_899;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_899 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_898 += tmp_899;
      result += (std::complex<double>(0,0.5)*vS*Lambdax*ZA(gI2,1)) * tmp_898
         ;
   }
   std::complex<double> tmp_900;
   std::complex<double> tmp_901;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_902;
      std::complex<double> tmp_903;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_903 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_902 += tmp_903;
      tmp_901 += (Conj(ZE(gI1,j2))) * tmp_902;
   }
   tmp_900 += tmp_901;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*ZA(gI2,2)) *
      tmp_900;
   if (gO2 < 3) {
      std::complex<double> tmp_904;
      std::complex<double> tmp_905;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_905 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_904 += tmp_905;
      result += (std::complex<double>(0,0.5)*vu*Lambdax*ZA(gI2,2)) * tmp_904
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_906;
   std::complex<double> tmp_907;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_907 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_906 += tmp_907;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_906;
   std::complex<double> tmp_908;
   std::complex<double> tmp_909;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_910;
      std::complex<double> tmp_911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_911 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_910 += tmp_911;
      tmp_909 += (Conj(ZE(gI1,j2))) * tmp_910;
   }
   tmp_908 += tmp_909;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_908;
   std::complex<double> tmp_912;
   std::complex<double> tmp_913;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_914;
      std::complex<double> tmp_915;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_916;
         std::complex<double> tmp_917;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_917 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_916 += tmp_917;
         tmp_915 += (KroneckerDelta(gO2,3 + j2)) * tmp_916;
      }
      tmp_914 += tmp_915;
      tmp_913 += (Conj(ZE(gI1,3 + j3))) * tmp_914;
   }
   tmp_912 += tmp_913;
   result += (-(vd*ZH(gI2,0))) * tmp_912;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_918;
      std::complex<double> tmp_919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_919 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_918 += tmp_919;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_918;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_920;
      std::complex<double> tmp_921;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_922;
         std::complex<double> tmp_923;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_923 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_922 += tmp_923;
         tmp_921 += (Conj(ZE(gI1,j2))) * tmp_922;
      }
      tmp_920 += tmp_921;
      result += (-(vd*ZH(gI2,0))) * tmp_920;
   }
   std::complex<double> tmp_924;
   std::complex<double> tmp_925;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_925 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_924 += tmp_925;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_924;
   std::complex<double> tmp_926;
   std::complex<double> tmp_927;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_928;
      std::complex<double> tmp_929;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_929 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_928 += tmp_929;
      tmp_927 += (Conj(ZE(gI1,j2))) * tmp_928;
   }
   tmp_926 += tmp_927;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,1)) * tmp_926;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_930;
      std::complex<double> tmp_931;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_931 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_930 += tmp_931;
      result += (0.5*vS*Lambdax*ZH(gI2,1)) * tmp_930;
   }
   std::complex<double> tmp_932;
   std::complex<double> tmp_933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_934;
      std::complex<double> tmp_935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_935 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_934 += tmp_935;
      tmp_933 += (Conj(ZE(gI1,j2))) * tmp_934;
   }
   tmp_932 += tmp_933;
   result += (0.5*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_932;
   if (gO2 < 3) {
      std::complex<double> tmp_936;
      std::complex<double> tmp_937;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_937 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_936 += tmp_937;
      result += (0.5*vu*Lambdax*ZH(gI2,2)) * tmp_936;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVWmSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZV(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVPSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_938;
   std::complex<double> tmp_939;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_939 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_938 += tmp_939;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_938;
   if (gO2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVZSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_941 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_940 += tmp_941;
   result += (0.7745966692414834*g1*Sin(ThetaW())) * tmp_940;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(vd*KroneckerDelta(0,gO2) + vu*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(vd*KroneckerDelta(0,gO2) + vu*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.025*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1)))) + KroneckerDelta(1,gO1)*(-5*
      KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) +
      ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1
      ,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(Conj(Kappa)*Lambdax*ZP(gI1,1)*
      ZP(gI2,0) + Conj(Lambdax)*(Lambdax*ZP(gI1,1)*ZP(gI2,1) + ZP(gI1,0)*(Lambdax*
      ZP(gI2,0) + Kappa*ZP(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(ZP(gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2)
      )*ZP(gI2,0) + 5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,1)) + ZP(gI1,1)*(5*
      vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP
      (gI2,1)))) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*(vu*(3*Sqr(g1) - 5*Sqr(g2))*ZP
      (gI2,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,1)) - ZP(gI1,1)*(5*vd*(
      -2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,
      1))) - 10*KroneckerDelta(2,gO2)*(2*vS*Conj(Kappa)*Lambdax*ZP(gI1,1)*ZP(gI2,0
      ) + 1.4142135623730951*(TLambdax*ZP(gI1,1)*ZP(gI2,0) + Conj(TLambdax)*ZP(gI1
      ,0)*ZP(gI2,1)) + 2*vS*Conj(Lambdax)*(Lambdax*ZP(gI1,1)*ZP(gI2,1) + ZP(gI1,0)
      *(Lambdax*ZP(gI2,0) + Kappa*ZP(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*KroneckerDelta(0,gO2)*UM(gI1,1)*UP(gI2,0) +
      (g2*KroneckerDelta(1,gO2)*UM(gI1,0) + Conj(Lambdax)*KroneckerDelta(2,gO2)*
      UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*Conj(UP(gI1,1))*
      KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*KroneckerDelta(0
      ,gO1) + Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(
      g2))*ZA(gI1,1)*ZA(gI2,1) - 20*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + 10*(
      Conj(Lambdax)*Kappa + Conj(Kappa)*Lambdax)*(-(KroneckerDelta(0,gO2)*ZA(gI1,2
      )*ZA(gI2,2)) + KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2
      ,2)))) + KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(
      g2))*ZA(gI1,0)*ZA(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZA(
      gI1,1)*ZA(gI2,1) + 20*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))) + 10*(Conj(
      Lambdax)*Kappa + Conj(Kappa)*Lambdax)*(-(KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(
      gI2,2)) + KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)))
      ) - 10*KroneckerDelta(2,gO1)*(Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,0
      )*(2*Lambdax*ZA(gI2,0) + Kappa*ZA(gI2,1)) + ZA(gI1,1)*(Kappa*ZA(gI2,0) + 2*
      Lambdax*ZA(gI2,1))) - Kappa*(KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*
      ZA(gI2,2)))) + Conj(Kappa)*(KroneckerDelta(2,gO2)*(Lambdax*ZA(gI1,1)*ZA(gI2,
      0) + Lambdax*ZA(gI1,0)*ZA(gI2,1) + 4*Kappa*ZA(gI1,2)*ZA(gI2,2)) - Lambdax*(
      KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) +
      KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*(3*(3*Sqr(g1)
      + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(
      g2))*ZH(gI1,1)*ZH(gI2,1) + 20*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2))) +
      KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,
      1)*ZH(gI2,0) + (-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(
      gI2,1) + 10*(Conj(Lambdax)*Kappa + Conj(Kappa)*Lambdax)*ZH(gI1,2)*ZH(gI2,2))
      + 10*KroneckerDelta(2,gO2)*(Conj(Kappa)*Lambdax*(ZH(gI1,2)*ZH(gI2,1) + ZH(
      gI1,1)*ZH(gI2,2)) + Conj(Lambdax)*(ZH(gI1,2)*(-2*Lambdax*ZH(gI2,0) + Kappa*
      ZH(gI2,1)) + (-2*Lambdax*ZH(gI1,0) + Kappa*ZH(gI1,1))*ZH(gI2,2)))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(
      g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*
      ZH(gI2,1) - 20*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*
      ((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,0) + (-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,1) + 10*(Conj(
      Lambdax)*Kappa + Conj(Kappa)*Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + 10*
      KroneckerDelta(2,gO2)*(Conj(Kappa)*Lambdax*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*
      ZH(gI2,2)) + Conj(Lambdax)*(ZH(gI1,2)*(Kappa*ZH(gI2,0) - 2*Lambdax*ZH(gI2,1)
      ) + (Kappa*ZH(gI1,0) - 2*Lambdax*ZH(gI1,1))*ZH(gI2,2)))) - 10*KroneckerDelta
      (2,gO1)*(Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZH(gI1,0)*(2*Lambdax*ZH(gI2,0
      ) - Kappa*ZH(gI2,1)) + ZH(gI1,1)*(-(Kappa*ZH(gI2,0)) + 2*Lambdax*ZH(gI2,1)))
      + KroneckerDelta(0,gO2)*(ZH(gI1,2)*(2*Lambdax*ZH(gI2,0) - Kappa*ZH(gI2,1))
      + (2*Lambdax*ZH(gI1,0) - Kappa*ZH(gI1,1))*ZH(gI2,2)) - KroneckerDelta(1,gO2)
      *(ZH(gI1,2)*(Kappa*ZH(gI2,0) - 2*Lambdax*ZH(gI2,1)) + (Kappa*ZH(gI1,0) - 2*
      Lambdax*ZH(gI1,1))*ZH(gI2,2))) - Conj(Kappa)*(KroneckerDelta(2,gO2)*(Lambdax
      *ZH(gI1,1)*ZH(gI2,0) + Lambdax*ZH(gI1,0)*ZH(gI2,1) - 12*Kappa*ZH(gI1,2)*ZH(
      gI2,2)) + Lambdax*(KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH
      (gI2,2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO2)*(-7.0710678118654755*Conj(TLambdax)*ZA(
      gI1,2)*ZA(gI2,0) + 10*vS*Conj(Lambdax)*Kappa*ZA(gI1,2)*ZA(gI2,0) + 10*vS*
      Conj(Kappa)*Lambdax*ZA(gI1,2)*ZA(gI2,0) - 7.0710678118654755*TLambdax*ZA(gI1
      ,2)*ZA(gI2,0) - 3*vu*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1) - 5*vu*Sqr(g2)*ZA(gI1,1)*ZA
      (gI2,1) - 20*vu*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2) - 10*vd*Conj(Lambdax)*
      Kappa*ZA(gI1,2)*ZA(gI2,2) - 10*vd*Conj(Kappa)*Lambdax*ZA(gI1,2)*ZA(gI2,2) +
      ZA(gI1,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI2,0) + 5*(2
      *vS*Conj(Lambdax)*Kappa + 2*vS*Conj(Kappa)*Lambdax - 1.4142135623730951*(
      Conj(TLambdax) + TLambdax))*ZA(gI2,2))) - 5*KroneckerDelta(2,gO2)*(
      1.4142135623730951*(Conj(TLambdax)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1
      )) + TLambdax*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) - 2*(Conj(TKappa)
      + TKappa)*ZA(gI1,2)*ZA(gI2,2)) + 2*Conj(Lambdax)*(-(Kappa*ZA(gI1,2)*(vu*ZA(
      gI2,0) + vd*ZA(gI2,1))) + ZA(gI1,1)*(vS*Kappa*ZA(gI2,0) + 2*vS*Lambdax*ZA(
      gI2,1) - vd*Kappa*ZA(gI2,2)) + ZA(gI1,0)*(2*vS*Lambdax*ZA(gI2,0) + vS*Kappa*
      ZA(gI2,1) - vu*Kappa*ZA(gI2,2))) - 2*Conj(Kappa)*(Lambdax*ZA(gI1,0)*(-(vS*ZA
      (gI2,1)) + vu*ZA(gI2,2)) + ZA(gI1,2)*(vu*Lambdax*ZA(gI2,0) + vd*Lambdax*ZA(
      gI2,1) - 4*vS*Kappa*ZA(gI2,2)) + ZA(gI1,1)*(-(vS*Lambdax*ZA(gI2,0)) + vd*
      Lambdax*ZA(gI2,2)))) - KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZA(
      gI1,0)*ZA(gI2,0) + ZA(gI1,1)*(vd*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2)
      )*ZA(gI2,1) + 5*(-2*vS*Conj(Lambdax)*Kappa - 2*vS*Conj(Kappa)*Lambdax +
      1.4142135623730951*(Conj(TLambdax) + TLambdax))*ZA(gI2,2)) + 5*ZA(gI1,2)*(
      1.4142135623730951*(Conj(TLambdax) + TLambdax)*ZA(gI2,1) + Conj(Kappa)*(-2*
      vS*Lambdax*ZA(gI2,1) + 2*vu*Lambdax*ZA(gI2,2)) + Conj(Lambdax)*(-2*vS*Kappa*
      ZA(gI2,1) + 2*(vu*Kappa + 2*vd*Lambdax)*ZA(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(-2*Conj(Lambdax)*Kappa*(-(
      KroneckerDelta(1,gO2)*(-(vS*ZA(gI2,0)*ZH(gI1,2)) + ZA(gI2,2)*(vS*ZH(gI1,0) +
      vd*ZH(gI1,2)))) + KroneckerDelta(2,gO2)*(-(ZA(gI2,2)*(vu*ZH(gI1,0) + vd*ZH(
      gI1,1))) + ZA(gI2,1)*(vS*ZH(gI1,0) + vd*ZH(gI1,2)) + ZA(gI2,0)*(vS*ZH(gI1,1)
      + vu*ZH(gI1,2))) - KroneckerDelta(0,gO2)*(-(vS*ZA(gI2,1)*ZH(gI1,2)) + ZA(
      gI2,2)*(vS*ZH(gI1,1) + vu*ZH(gI1,2)))) + 2*Conj(Kappa)*Lambdax*(-(
      KroneckerDelta(1,gO2)*(-(vS*ZA(gI2,0)*ZH(gI1,2)) + ZA(gI2,2)*(vS*ZH(gI1,0) +
      vd*ZH(gI1,2)))) + KroneckerDelta(2,gO2)*(-(ZA(gI2,2)*(vu*ZH(gI1,0) + vd*ZH(
      gI1,1))) + ZA(gI2,1)*(vS*ZH(gI1,0) + vd*ZH(gI1,2)) + ZA(gI2,0)*(vS*ZH(gI1,1)
      + vu*ZH(gI1,2))) - KroneckerDelta(0,gO2)*(-(vS*ZA(gI2,1)*ZH(gI1,2)) + ZA(
      gI2,2)*(vS*ZH(gI1,1) + vu*ZH(gI1,2)))) + 1.4142135623730951*(KroneckerDelta(
      2,gO2)*(TLambdax*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + 2*(Conj(
      TKappa) - TKappa)*ZA(gI2,2)*ZH(gI1,2)) + TLambdax*(KroneckerDelta(1,gO2)*(ZA
      (gI2,2)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)) + KroneckerDelta(0,gO2)*(ZA(gI2,2)*
      ZH(gI1,1) + ZA(gI2,1)*ZH(gI1,2))) - Conj(TLambdax)*(KroneckerDelta(2,gO2)*(
      ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + KroneckerDelta(1,gO2)*(ZA(gI2,2
      )*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)) + KroneckerDelta(0,gO2)*(ZA(gI2,2)*ZH(gI1
      ,1) + ZA(gI2,1)*ZH(gI1,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO2)*(-(ZH(gI1,0)*(3*vd*(3*Sqr(g1) + 5*Sqr(
      g2))*ZH(gI2,0) + vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI2,1) +
      20*vS*AbsSqr(Lambdax)*ZH(gI2,2))) + ZH(gI1,1)*(vu*(-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr
      (g2))*ZH(gI2,1) + 5*(2*vS*Conj(Lambdax)*Kappa + 2*vS*Conj(Kappa)*Lambdax +
      1.4142135623730951*(Conj(TLambdax) + TLambdax))*ZH(gI2,2)) + 5*ZH(gI1,2)*(
      1.4142135623730951*(Conj(TLambdax) + TLambdax)*ZH(gI2,1) + 2*Conj(Kappa)*
      Lambdax*(vS*ZH(gI2,1) + vu*ZH(gI2,2)) + 2*Conj(Lambdax)*(-2*vS*Lambdax*ZH(
      gI2,0) + vS*Kappa*ZH(gI2,1) + (vu*Kappa - 2*vd*Lambdax)*ZH(gI2,2)))) +
      KroneckerDelta(1,gO2)*(ZH(gI1,1)*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2))*ZH(gI2,0) - 3*vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) - 20*vS*AbsSqr(
      Lambdax)*ZH(gI2,2)) + ZH(gI1,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr
      (g2))*ZH(gI2,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1)
      + 5*(2*vS*Conj(Lambdax)*Kappa + 2*vS*Conj(Kappa)*Lambdax +
      1.4142135623730951*(Conj(TLambdax) + TLambdax))*ZH(gI2,2)) + 5*ZH(gI1,2)*(
      1.4142135623730951*(Conj(TLambdax) + TLambdax)*ZH(gI2,0) + 2*Conj(Kappa)*
      Lambdax*(vS*ZH(gI2,0) + vd*ZH(gI2,2)) + 2*Conj(Lambdax)*(vS*Kappa*ZH(gI2,0)
      - 2*vS*Lambdax*ZH(gI2,1) + (vd*Kappa - 2*vu*Lambdax)*ZH(gI2,2)))) - 5*
      KroneckerDelta(2,gO2)*(-1.4142135623730951*(Conj(TLambdax)*(ZH(gI1,1)*ZH(gI2
      ,0) + ZH(gI1,0)*ZH(gI2,1)) + TLambdax*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(
      gI2,1)) - 2*(Conj(TKappa) + TKappa)*ZH(gI1,2)*ZH(gI2,2)) - 2*Conj(Kappa)*(
      Lambdax*ZH(gI1,1)*(vS*ZH(gI2,0) + vd*ZH(gI2,2)) + Lambdax*ZH(gI1,0)*(vS*ZH(
      gI2,1) + vu*ZH(gI2,2)) + ZH(gI1,2)*(vu*Lambdax*ZH(gI2,0) + vd*Lambdax*ZH(gI2
      ,1) - 12*vS*Kappa*ZH(gI2,2))) - 2*Conj(Lambdax)*(ZH(gI1,2)*((vu*Kappa - 2*vd
      *Lambdax)*ZH(gI2,0) + (vd*Kappa - 2*vu*Lambdax)*ZH(gI2,1)) + ZH(gI1,0)*(-2*
      vS*Lambdax*ZH(gI2,0) + vS*Kappa*ZH(gI2,1) + (vu*Kappa - 2*vd*Lambdax)*ZH(gI2
      ,2)) + ZH(gI1,1)*(vS*Kappa*ZH(gI2,0) - 2*vS*Lambdax*ZH(gI2,1) + (vd*Kappa -
      2*vu*Lambdax)*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_942;
   std::complex<double> tmp_943;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_944;
      std::complex<double> tmp_945;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_945 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_944 += tmp_945;
      tmp_943 += (ZDL(gI1,j2)) * tmp_944;
   }
   tmp_942 += tmp_943;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_942;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_946;
   std::complex<double> tmp_947;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_948;
      std::complex<double> tmp_949;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_949 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_948 += tmp_949;
      tmp_947 += (Conj(ZDL(gI2,j2))) * tmp_948;
   }
   tmp_946 += tmp_947;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_946;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_950;
   std::complex<double> tmp_951;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_952;
      std::complex<double> tmp_953;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_953 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_952 += tmp_953;
      tmp_951 += (ZEL(gI1,j2)) * tmp_952;
   }
   tmp_950 += tmp_951;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_950;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_954;
   std::complex<double> tmp_955;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_956;
      std::complex<double> tmp_957;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_957 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_956 += tmp_957;
      tmp_955 += (Conj(ZEL(gI2,j2))) * tmp_956;
   }
   tmp_954 += tmp_955;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_954;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_958;
   std::complex<double> tmp_959;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_960;
      std::complex<double> tmp_961;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_961 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_960 += tmp_961;
      tmp_959 += (ZUL(gI1,j2)) * tmp_960;
   }
   tmp_958 += tmp_959;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_958;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_962;
   std::complex<double> tmp_963;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_964;
      std::complex<double> tmp_965;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_965 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_964 += tmp_965;
      tmp_963 += (Conj(ZUL(gI2,j2))) * tmp_964;
   }
   tmp_962 += tmp_963;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_962;

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO2)*(ZN(gI1,2)*(3.872983346207417*g1*ZN(gI2,
      0) - 5*g2*ZN(gI2,1)) + 3.872983346207417*g1*ZN(gI1,0)*ZN(gI2,2) - 5*g2*ZN(
      gI1,1)*ZN(gI2,2) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,3) +
      7.0710678118654755*Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4)) + 7.0710678118654755*
      KroneckerDelta(2,gO2)*(Conj(Lambdax)*(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2
      ,3)) - 2*Conj(Kappa)*ZN(gI1,4)*ZN(gI2,4)) + KroneckerDelta(1,gO2)*(ZN(gI1,3)
      *(-3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(gI2,1)) + (-3.872983346207417*g1
      *ZN(gI1,0) + 5*g2*ZN(gI1,1))*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*(
      ZN(gI1,4)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) -
      3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) +
      5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(
      gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN
      (gI1,0))*(Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) - Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1)) - 14.142135623730951*Conj(ZN(gI1,4))*Conj(ZN(gI2,4))*
      KroneckerDelta(2,gO1)*Kappa + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(gI2
      ,3))*KroneckerDelta(0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,3))*Conj
      (ZN(gI2,4))*KroneckerDelta(0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,4
      ))*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1)*Lambdax + 7.0710678118654755*Conj(
      ZN(gI1,3))*Conj(ZN(gI2,2))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(
      3.872983346207417*g1*Conj(ZN(gI2,0))*KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(
      gI2,1))*KroneckerDelta(0,gO1) + 7.0710678118654755*(Conj(ZN(gI2,4))*
      KroneckerDelta(1,gO1) + Conj(ZN(gI2,3))*KroneckerDelta(2,gO1))*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_966;
   std::complex<double> tmp_967;
   std::complex<double> tmp_968;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_968 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_967 += tmp_968;
   tmp_966 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_967;
   std::complex<double> tmp_969;
   std::complex<double> tmp_970;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_970 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_969 += tmp_970;
   tmp_966 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_969;
   std::complex<double> tmp_971;
   std::complex<double> tmp_972;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_972 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_971 += tmp_972;
   tmp_966 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_971;
   std::complex<double> tmp_973;
   std::complex<double> tmp_974;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_974 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_973 += tmp_974;
   tmp_966 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_973;
   std::complex<double> tmp_975;
   std::complex<double> tmp_976;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_976 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_975 += tmp_976;
   tmp_966 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*Sqr(g1)) * tmp_975;
   std::complex<double> tmp_977;
   std::complex<double> tmp_978;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_978 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_977 += tmp_978;
   tmp_966 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_977;
   std::complex<double> tmp_979;
   std::complex<double> tmp_980;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_981;
      std::complex<double> tmp_982;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_982 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_981 += tmp_982;
      tmp_980 += (Conj(ZD(gI2,j2))) * tmp_981;
   }
   tmp_979 += tmp_980;
   tmp_966 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)) * tmp_979;
   std::complex<double> tmp_983;
   std::complex<double> tmp_984;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_985;
      std::complex<double> tmp_986;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_986 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_985 += tmp_986;
      tmp_984 += (Conj(ZD(gI2,j2))) * tmp_985;
   }
   tmp_983 += tmp_984;
   tmp_966 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)) * tmp_983;
   std::complex<double> tmp_987;
   std::complex<double> tmp_988;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_989;
      std::complex<double> tmp_990;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_990 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_989 += tmp_990;
      tmp_988 += (ZD(gI1,j2)) * tmp_989;
   }
   tmp_987 += tmp_988;
   tmp_966 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*KroneckerDelta
      (2,gO1)*Lambdax) * tmp_987;
   std::complex<double> tmp_991;
   std::complex<double> tmp_992;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_993;
      std::complex<double> tmp_994;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_994 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_993 += tmp_994;
      tmp_992 += (ZD(gI1,j2)) * tmp_993;
   }
   tmp_991 += tmp_992;
   tmp_966 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*KroneckerDelta
      (2,gO2)*Lambdax) * tmp_991;
   std::complex<double> tmp_995;
   std::complex<double> tmp_996;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_997;
      std::complex<double> tmp_998;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_999;
         std::complex<double> tmp_1000;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1000 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_999 += tmp_1000;
         tmp_998 += (ZD(gI1,3 + j2)) * tmp_999;
      }
      tmp_997 += tmp_998;
      tmp_996 += (Conj(ZD(gI2,3 + j3))) * tmp_997;
   }
   tmp_995 += tmp_996;
   tmp_966 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_995;
   std::complex<double> tmp_1001;
   std::complex<double> tmp_1002;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1003;
      std::complex<double> tmp_1004;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1005;
         std::complex<double> tmp_1006;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1006 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1005 += tmp_1006;
         tmp_1004 += (Conj(ZD(gI2,j2))) * tmp_1005;
      }
      tmp_1003 += tmp_1004;
      tmp_1002 += (ZD(gI1,j3)) * tmp_1003;
   }
   tmp_1001 += tmp_1002;
   tmp_966 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_1001;
   result += (std::complex<double>(0,-1)) * tmp_966;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1007;
   std::complex<double> tmp_1008;
   std::complex<double> tmp_1009;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1009 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1008 += tmp_1009;
   tmp_1007 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1008;
   std::complex<double> tmp_1010;
   std::complex<double> tmp_1011;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1011 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1010 += tmp_1011;
   tmp_1007 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1010;
   std::complex<double> tmp_1012;
   std::complex<double> tmp_1013;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1013 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1012 += tmp_1013;
   tmp_1007 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1012;
   std::complex<double> tmp_1014;
   std::complex<double> tmp_1015;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1015 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1014 += tmp_1015;
   tmp_1007 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1014;
   std::complex<double> tmp_1016;
   std::complex<double> tmp_1017;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1017 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1016 += tmp_1017;
   tmp_1007 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1016;
   std::complex<double> tmp_1018;
   std::complex<double> tmp_1019;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1019 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1018 += tmp_1019;
   tmp_1007 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1018;
   std::complex<double> tmp_1020;
   std::complex<double> tmp_1021;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1022;
      std::complex<double> tmp_1023;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1023 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1022 += tmp_1023;
      tmp_1021 += (Conj(ZE(gI2,j2))) * tmp_1022;
   }
   tmp_1020 += tmp_1021;
   tmp_1007 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1)) * tmp_1020;
   std::complex<double> tmp_1024;
   std::complex<double> tmp_1025;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1026;
      std::complex<double> tmp_1027;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1027 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1026 += tmp_1027;
      tmp_1025 += (Conj(ZE(gI2,j2))) * tmp_1026;
   }
   tmp_1024 += tmp_1025;
   tmp_1007 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)
      *KroneckerDelta(2,gO2)) * tmp_1024;
   std::complex<double> tmp_1028;
   std::complex<double> tmp_1029;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1030;
      std::complex<double> tmp_1031;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1031 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1030 += tmp_1031;
      tmp_1029 += (ZE(gI1,j2)) * tmp_1030;
   }
   tmp_1028 += tmp_1029;
   tmp_1007 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1028;
   std::complex<double> tmp_1032;
   std::complex<double> tmp_1033;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1034;
      std::complex<double> tmp_1035;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1035 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1034 += tmp_1035;
      tmp_1033 += (ZE(gI1,j2)) * tmp_1034;
   }
   tmp_1032 += tmp_1033;
   tmp_1007 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1032;
   std::complex<double> tmp_1036;
   std::complex<double> tmp_1037;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1038;
      std::complex<double> tmp_1039;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1040;
         std::complex<double> tmp_1041;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1041 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1040 += tmp_1041;
         tmp_1039 += (ZE(gI1,3 + j2)) * tmp_1040;
      }
      tmp_1038 += tmp_1039;
      tmp_1037 += (Conj(ZE(gI2,3 + j3))) * tmp_1038;
   }
   tmp_1036 += tmp_1037;
   tmp_1007 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1036;
   std::complex<double> tmp_1042;
   std::complex<double> tmp_1043;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1044;
      std::complex<double> tmp_1045;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1046;
         std::complex<double> tmp_1047;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1047 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1046 += tmp_1047;
         tmp_1045 += (Conj(ZE(gI2,j2))) * tmp_1046;
      }
      tmp_1044 += tmp_1045;
      tmp_1043 += (ZE(gI1,j3)) * tmp_1044;
   }
   tmp_1042 += tmp_1043;
   tmp_1007 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1042;
   result += (std::complex<double>(0,-1)) * tmp_1007;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1048;
   std::complex<double> tmp_1049;
   std::complex<double> tmp_1050;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1050 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1049 += tmp_1050;
   tmp_1048 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1049;
   std::complex<double> tmp_1051;
   std::complex<double> tmp_1052;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1052 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1051 += tmp_1052;
   tmp_1048 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1051;
   std::complex<double> tmp_1053;
   std::complex<double> tmp_1054;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1054 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1053 += tmp_1054;
   tmp_1048 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1053;
   std::complex<double> tmp_1055;
   std::complex<double> tmp_1056;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1056 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1055 += tmp_1056;
   tmp_1048 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1055;
   std::complex<double> tmp_1057;
   std::complex<double> tmp_1058;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1058 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1057 += tmp_1058;
   tmp_1048 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1057;
   std::complex<double> tmp_1059;
   std::complex<double> tmp_1060;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1060 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1059 += tmp_1060;
   tmp_1048 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1059;
   std::complex<double> tmp_1061;
   std::complex<double> tmp_1062;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1063;
      std::complex<double> tmp_1064;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1064 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1063 += tmp_1064;
      tmp_1062 += (Conj(ZU(gI2,j2))) * tmp_1063;
   }
   tmp_1061 += tmp_1062;
   tmp_1048 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2)
      *KroneckerDelta(2,gO1)) * tmp_1061;
   std::complex<double> tmp_1065;
   std::complex<double> tmp_1066;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1067;
      std::complex<double> tmp_1068;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1068 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1067 += tmp_1068;
      tmp_1066 += (Conj(ZU(gI2,j2))) * tmp_1067;
   }
   tmp_1065 += tmp_1066;
   tmp_1048 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1)
      *KroneckerDelta(2,gO2)) * tmp_1065;
   std::complex<double> tmp_1069;
   std::complex<double> tmp_1070;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1071;
      std::complex<double> tmp_1072;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1072 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1071 += tmp_1072;
      tmp_1070 += (ZU(gI1,j2)) * tmp_1071;
   }
   tmp_1069 += tmp_1070;
   tmp_1048 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1069;
   std::complex<double> tmp_1073;
   std::complex<double> tmp_1074;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1075;
      std::complex<double> tmp_1076;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1076 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1075 += tmp_1076;
      tmp_1074 += (ZU(gI1,j2)) * tmp_1075;
   }
   tmp_1073 += tmp_1074;
   tmp_1048 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1073;
   std::complex<double> tmp_1077;
   std::complex<double> tmp_1078;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1079;
      std::complex<double> tmp_1080;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1081;
         std::complex<double> tmp_1082;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1082 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1081 += tmp_1082;
         tmp_1080 += (ZU(gI1,3 + j2)) * tmp_1081;
      }
      tmp_1079 += tmp_1080;
      tmp_1078 += (Conj(ZU(gI2,3 + j3))) * tmp_1079;
   }
   tmp_1077 += tmp_1078;
   tmp_1048 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1077;
   std::complex<double> tmp_1083;
   std::complex<double> tmp_1084;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1085;
      std::complex<double> tmp_1086;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1087;
         std::complex<double> tmp_1088;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1088 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1087 += tmp_1088;
         tmp_1086 += (Conj(ZU(gI2,j2))) * tmp_1087;
      }
      tmp_1085 += tmp_1086;
      tmp_1084 += (ZU(gI1,j3)) * tmp_1085;
   }
   tmp_1083 += tmp_1084;
   tmp_1048 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1083;
   result += (std::complex<double>(0,-1)) * tmp_1048;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1089;
   std::complex<double> tmp_1090;
   std::complex<double> tmp_1091;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1091 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1090 += tmp_1091;
   tmp_1089 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1090;
   std::complex<double> tmp_1092;
   std::complex<double> tmp_1093;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1093 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1092 += tmp_1093;
   tmp_1089 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1092;
   std::complex<double> tmp_1094;
   std::complex<double> tmp_1095;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1095 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1094 += tmp_1095;
   tmp_1089 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1094;
   std::complex<double> tmp_1096;
   std::complex<double> tmp_1097;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1097 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1096 += tmp_1097;
   tmp_1089 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1096;
   std::complex<double> tmp_1098;
   std::complex<double> tmp_1099;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1099 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1098 += tmp_1099;
   tmp_1089 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1098;
   std::complex<double> tmp_1100;
   std::complex<double> tmp_1101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1101 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1100 += tmp_1101;
   tmp_1089 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1100;
   std::complex<double> tmp_1102;
   std::complex<double> tmp_1103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1104;
      std::complex<double> tmp_1105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1105 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1104 += tmp_1105;
      tmp_1103 += (Conj(ZD(gI2,j2))) * tmp_1104;
   }
   tmp_1102 += tmp_1103;
   tmp_1089 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_1102;
   std::complex<double> tmp_1106;
   std::complex<double> tmp_1107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1108;
      std::complex<double> tmp_1109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1109 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1108 += tmp_1109;
      tmp_1107 += (Conj(ZD(gI2,j2))) * tmp_1108;
   }
   tmp_1106 += tmp_1107;
   tmp_1089 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1106;
   std::complex<double> tmp_1110;
   std::complex<double> tmp_1111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1112;
      std::complex<double> tmp_1113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1113 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1112 += tmp_1113;
      tmp_1111 += (Conj(ZD(gI2,j2))) * tmp_1112;
   }
   tmp_1110 += tmp_1111;
   tmp_1089 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1110;
   std::complex<double> tmp_1114;
   std::complex<double> tmp_1115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1116;
      std::complex<double> tmp_1117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1117 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1116 += tmp_1117;
      tmp_1115 += (ZD(gI1,j2)) * tmp_1116;
   }
   tmp_1114 += tmp_1115;
   tmp_1089 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1114;
   std::complex<double> tmp_1118;
   std::complex<double> tmp_1119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1120;
      std::complex<double> tmp_1121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1121 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1120 += tmp_1121;
      tmp_1119 += (ZD(gI1,j2)) * tmp_1120;
   }
   tmp_1118 += tmp_1119;
   tmp_1089 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1118;
   std::complex<double> tmp_1122;
   std::complex<double> tmp_1123;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1124;
      std::complex<double> tmp_1125;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1125 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1124 += tmp_1125;
      tmp_1123 += (ZD(gI1,j2)) * tmp_1124;
   }
   tmp_1122 += tmp_1123;
   tmp_1089 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1122;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1127;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1128;
      std::complex<double> tmp_1129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1130;
         std::complex<double> tmp_1131;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1131 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1130 += tmp_1131;
         tmp_1129 += (ZD(gI1,3 + j2)) * tmp_1130;
      }
      tmp_1128 += tmp_1129;
      tmp_1127 += (Conj(ZD(gI2,3 + j3))) * tmp_1128;
   }
   tmp_1126 += tmp_1127;
   tmp_1089 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1126
      ;
   std::complex<double> tmp_1132;
   std::complex<double> tmp_1133;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1134;
      std::complex<double> tmp_1135;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1136;
         std::complex<double> tmp_1137;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1137 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1136 += tmp_1137;
         tmp_1135 += (Conj(ZD(gI2,j2))) * tmp_1136;
      }
      tmp_1134 += tmp_1135;
      tmp_1133 += (ZD(gI1,j3)) * tmp_1134;
   }
   tmp_1132 += tmp_1133;
   tmp_1089 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1132
      ;
   result += (std::complex<double>(0,-1)) * tmp_1089;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1138;
   std::complex<double> tmp_1139;
   std::complex<double> tmp_1140;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1140 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1139 += tmp_1140;
   tmp_1138 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1139;
   std::complex<double> tmp_1141;
   std::complex<double> tmp_1142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1142 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1141 += tmp_1142;
   tmp_1138 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1141;
   std::complex<double> tmp_1143;
   std::complex<double> tmp_1144;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1144 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1143 += tmp_1144;
   tmp_1138 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1143;
   std::complex<double> tmp_1145;
   std::complex<double> tmp_1146;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1146 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1145 += tmp_1146;
   tmp_1138 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1145;
   std::complex<double> tmp_1147;
   std::complex<double> tmp_1148;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1148 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1147 += tmp_1148;
   tmp_1138 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1147;
   std::complex<double> tmp_1149;
   std::complex<double> tmp_1150;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1150 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1149 += tmp_1150;
   tmp_1138 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1149;
   std::complex<double> tmp_1151;
   std::complex<double> tmp_1152;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1153;
      std::complex<double> tmp_1154;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1154 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1153 += tmp_1154;
      tmp_1152 += (Conj(ZE(gI2,j2))) * tmp_1153;
   }
   tmp_1151 += tmp_1152;
   tmp_1138 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_1151;
   std::complex<double> tmp_1155;
   std::complex<double> tmp_1156;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1157;
      std::complex<double> tmp_1158;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1158 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1157 += tmp_1158;
      tmp_1156 += (Conj(ZE(gI2,j2))) * tmp_1157;
   }
   tmp_1155 += tmp_1156;
   tmp_1138 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1155;
   std::complex<double> tmp_1159;
   std::complex<double> tmp_1160;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1161;
      std::complex<double> tmp_1162;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1162 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1161 += tmp_1162;
      tmp_1160 += (Conj(ZE(gI2,j2))) * tmp_1161;
   }
   tmp_1159 += tmp_1160;
   tmp_1138 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1159;
   std::complex<double> tmp_1163;
   std::complex<double> tmp_1164;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1165;
      std::complex<double> tmp_1166;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1166 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1165 += tmp_1166;
      tmp_1164 += (ZE(gI1,j2)) * tmp_1165;
   }
   tmp_1163 += tmp_1164;
   tmp_1138 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1163;
   std::complex<double> tmp_1167;
   std::complex<double> tmp_1168;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1169;
      std::complex<double> tmp_1170;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1170 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1169 += tmp_1170;
      tmp_1168 += (ZE(gI1,j2)) * tmp_1169;
   }
   tmp_1167 += tmp_1168;
   tmp_1138 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1167;
   std::complex<double> tmp_1171;
   std::complex<double> tmp_1172;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1173;
      std::complex<double> tmp_1174;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1174 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1173 += tmp_1174;
      tmp_1172 += (ZE(gI1,j2)) * tmp_1173;
   }
   tmp_1171 += tmp_1172;
   tmp_1138 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1171;
   std::complex<double> tmp_1175;
   std::complex<double> tmp_1176;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1177;
      std::complex<double> tmp_1178;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1179;
         std::complex<double> tmp_1180;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1180 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1179 += tmp_1180;
         tmp_1178 += (ZE(gI1,3 + j2)) * tmp_1179;
      }
      tmp_1177 += tmp_1178;
      tmp_1176 += (Conj(ZE(gI2,3 + j3))) * tmp_1177;
   }
   tmp_1175 += tmp_1176;
   tmp_1138 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1175
      ;
   std::complex<double> tmp_1181;
   std::complex<double> tmp_1182;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1183;
      std::complex<double> tmp_1184;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1185;
         std::complex<double> tmp_1186;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1186 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1185 += tmp_1186;
         tmp_1184 += (Conj(ZE(gI2,j2))) * tmp_1185;
      }
      tmp_1183 += tmp_1184;
      tmp_1182 += (ZE(gI1,j3)) * tmp_1183;
   }
   tmp_1181 += tmp_1182;
   tmp_1138 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1181
      ;
   result += (std::complex<double>(0,-1)) * tmp_1138;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1187;
   std::complex<double> tmp_1188;
   std::complex<double> tmp_1189;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1189 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1188 += tmp_1189;
   tmp_1187 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1188;
   std::complex<double> tmp_1190;
   std::complex<double> tmp_1191;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1191 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1190 += tmp_1191;
   tmp_1187 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1190;
   std::complex<double> tmp_1192;
   std::complex<double> tmp_1193;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1193 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1192 += tmp_1193;
   tmp_1187 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1192;
   std::complex<double> tmp_1194;
   std::complex<double> tmp_1195;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1195 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1194 += tmp_1195;
   tmp_1187 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1194;
   std::complex<double> tmp_1196;
   std::complex<double> tmp_1197;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1197 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1196 += tmp_1197;
   tmp_1187 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1196;
   std::complex<double> tmp_1198;
   std::complex<double> tmp_1199;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1199 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1198 += tmp_1199;
   tmp_1187 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1198;
   std::complex<double> tmp_1200;
   std::complex<double> tmp_1201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1202;
      std::complex<double> tmp_1203;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1203 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1202 += tmp_1203;
      tmp_1201 += (Conj(ZU(gI2,j2))) * tmp_1202;
   }
   tmp_1200 += tmp_1201;
   tmp_1187 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(0,
      gO2)) * tmp_1200;
   std::complex<double> tmp_1204;
   std::complex<double> tmp_1205;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1206;
      std::complex<double> tmp_1207;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1207 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1206 += tmp_1207;
      tmp_1205 += (Conj(ZU(gI2,j2))) * tmp_1206;
   }
   tmp_1204 += tmp_1205;
   tmp_1187 += (std::complex<double>(0,0.5)*vd*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1204;
   std::complex<double> tmp_1208;
   std::complex<double> tmp_1209;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1210;
      std::complex<double> tmp_1211;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1211 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1210 += tmp_1211;
      tmp_1209 += (Conj(ZU(gI2,j2))) * tmp_1210;
   }
   tmp_1208 += tmp_1209;
   tmp_1187 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1208;
   std::complex<double> tmp_1212;
   std::complex<double> tmp_1213;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1214;
      std::complex<double> tmp_1215;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1215 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1214 += tmp_1215;
      tmp_1213 += (ZU(gI1,j2)) * tmp_1214;
   }
   tmp_1212 += tmp_1213;
   tmp_1187 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(0,gO2)*Lambdax) *
      tmp_1212;
   std::complex<double> tmp_1216;
   std::complex<double> tmp_1217;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1218;
      std::complex<double> tmp_1219;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1219 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1218 += tmp_1219;
      tmp_1217 += (ZU(gI1,j2)) * tmp_1218;
   }
   tmp_1216 += tmp_1217;
   tmp_1187 += (std::complex<double>(0,0.5)*vd*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1216;
   std::complex<double> tmp_1220;
   std::complex<double> tmp_1221;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1222;
      std::complex<double> tmp_1223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1223 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1222 += tmp_1223;
      tmp_1221 += (ZU(gI1,j2)) * tmp_1222;
   }
   tmp_1220 += tmp_1221;
   tmp_1187 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1220;
   std::complex<double> tmp_1224;
   std::complex<double> tmp_1225;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1226;
      std::complex<double> tmp_1227;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1228;
         std::complex<double> tmp_1229;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1229 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1228 += tmp_1229;
         tmp_1227 += (ZU(gI1,3 + j2)) * tmp_1228;
      }
      tmp_1226 += tmp_1227;
      tmp_1225 += (Conj(ZU(gI2,3 + j3))) * tmp_1226;
   }
   tmp_1224 += tmp_1225;
   tmp_1187 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1224
      ;
   std::complex<double> tmp_1230;
   std::complex<double> tmp_1231;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1232;
      std::complex<double> tmp_1233;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1234;
         std::complex<double> tmp_1235;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1235 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1234 += tmp_1235;
         tmp_1233 += (Conj(ZU(gI2,j2))) * tmp_1234;
      }
      tmp_1232 += tmp_1233;
      tmp_1231 += (ZU(gI1,j3)) * tmp_1232;
   }
   tmp_1230 += tmp_1231;
   tmp_1187 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1230
      ;
   result += (std::complex<double>(0,-1)) * tmp_1187;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) - KroneckerDelta(1,gO2)*ZP(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZA(gI2,0) -
      KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(vd*KroneckerDelta(0,gO1) - vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vd*KroneckerDelta(0,gO1) - vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(-5*KroneckerDelta(1,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1)))) + KroneckerDelta(1,gO1)*(5*
      KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) +
      ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1
      ,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(-(Conj(Kappa)*Lambdax*ZP(gI1,1)
      *ZP(gI2,0)) + Conj(Lambdax)*(Lambdax*ZP(gI1,1)*ZP(gI2,1) + ZP(gI1,0)*(
      Lambdax*ZP(gI2,0) - Kappa*ZP(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vu*KroneckerDelta(0,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1)) + vd*
      KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) -
      ZP(gI1,0)*ZP(gI2,1)) + 2*KroneckerDelta(2,gO2)*(2*vS*Conj(Kappa)*Lambdax*ZP(
      gI1,1)*ZP(gI2,0) - 1.4142135623730951*TLambdax*ZP(gI1,1)*ZP(gI2,0) + (
      1.4142135623730951*Conj(TLambdax) - 2*vS*Conj(Lambdax)*Kappa)*ZP(gI1,0)*ZP(
      gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2*KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) + (g2*KroneckerDelta(1,gO2)*UM(gI1,0) - Conj(
      Lambdax)*KroneckerDelta(2,gO2)*UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gI2,0))*
      Conj(UP(gI1,1))*KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) - Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*(3*(3*Sqr(g1)
      + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(
      g2))*ZA(gI1,1)*ZA(gI2,1) + 20*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))) +
      KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,
      1)*ZA(gI2,0) + (-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(
      gI2,1) + 10*(Conj(Lambdax)*Kappa + Conj(Kappa)*Lambdax)*ZA(gI1,2)*ZA(gI2,2))
      + 10*KroneckerDelta(2,gO2)*(Conj(Kappa)*Lambdax*(ZA(gI1,2)*ZA(gI2,1) + ZA(
      gI1,1)*ZA(gI2,2)) + Conj(Lambdax)*(ZA(gI1,2)*(-2*Lambdax*ZA(gI2,0) + Kappa*
      ZA(gI2,1)) + (-2*Lambdax*ZA(gI1,0) + Kappa*ZA(gI1,1))*ZA(gI2,2)))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(
      g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*
      ZA(gI2,1) - 20*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*
      ((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,0) + (-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,1) + 10*(Conj(
      Lambdax)*Kappa + Conj(Kappa)*Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + 10*
      KroneckerDelta(2,gO2)*(Conj(Kappa)*Lambdax*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*
      ZA(gI2,2)) + Conj(Lambdax)*(ZA(gI1,2)*(Kappa*ZA(gI2,0) - 2*Lambdax*ZA(gI2,1)
      ) + (Kappa*ZA(gI1,0) - 2*Lambdax*ZA(gI1,1))*ZA(gI2,2)))) - 10*KroneckerDelta
      (2,gO1)*(Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*(2*Lambdax*ZA(gI2,0
      ) - Kappa*ZA(gI2,1)) + ZA(gI1,1)*(-(Kappa*ZA(gI2,0)) + 2*Lambdax*ZA(gI2,1)))
      + KroneckerDelta(0,gO2)*(ZA(gI1,2)*(2*Lambdax*ZA(gI2,0) - Kappa*ZA(gI2,1))
      + (2*Lambdax*ZA(gI1,0) - Kappa*ZA(gI1,1))*ZA(gI2,2)) - KroneckerDelta(1,gO2)
      *(ZA(gI1,2)*(Kappa*ZA(gI2,0) - 2*Lambdax*ZA(gI2,1)) + (Kappa*ZA(gI1,0) - 2*
      Lambdax*ZA(gI1,1))*ZA(gI2,2))) - Conj(Kappa)*(KroneckerDelta(2,gO2)*(Lambdax
      *ZA(gI1,1)*ZA(gI2,0) + Lambdax*ZA(gI1,0)*ZA(gI2,1) - 12*Kappa*ZA(gI1,2)*ZA(
      gI2,2)) + Lambdax*(KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA
      (gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(
      g2))*ZH(gI1,1)*ZH(gI2,1) - 20*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + 10*(
      Conj(Lambdax)*Kappa + Conj(Kappa)*Lambdax)*(-(KroneckerDelta(0,gO2)*ZH(gI1,2
      )*ZH(gI2,2)) + KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2
      ,2)))) + KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(
      g2))*ZH(gI1,0)*ZH(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(
      gI1,1)*ZH(gI2,1) + 20*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2))) + 10*(Conj(
      Lambdax)*Kappa + Conj(Kappa)*Lambdax)*(-(KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(
      gI2,2)) + KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))
      ) - 10*KroneckerDelta(2,gO1)*(Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZH(gI1,0
      )*(2*Lambdax*ZH(gI2,0) + Kappa*ZH(gI2,1)) + ZH(gI1,1)*(Kappa*ZH(gI2,0) + 2*
      Lambdax*ZH(gI2,1))) - Kappa*(KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH
      (gI1,0)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*
      ZH(gI2,2)))) + Conj(Kappa)*(KroneckerDelta(2,gO2)*(Lambdax*ZH(gI1,1)*ZH(gI2,
      0) + Lambdax*ZH(gI1,0)*ZH(gI2,1) + 4*Kappa*ZH(gI1,2)*ZH(gI2,2)) - Lambdax*(
      KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(-2*Conj(Lambdax)*Kappa*(
      KroneckerDelta(1,gO2)*(vS*ZA(gI1,0)*ZA(gI2,2) + ZA(gI1,2)*(vS*ZA(gI2,0) - vd
      *ZA(gI2,2))) + KroneckerDelta(0,gO2)*(vS*ZA(gI1,1)*ZA(gI2,2) + ZA(gI1,2)*(vS
      *ZA(gI2,1) - vu*ZA(gI2,2))) - KroneckerDelta(2,gO2)*(ZA(gI1,2)*(vu*ZA(gI2,0)
      + vd*ZA(gI2,1)) + ZA(gI1,1)*(-(vS*ZA(gI2,0)) + vd*ZA(gI2,2)) + ZA(gI1,0)*(-
      (vS*ZA(gI2,1)) + vu*ZA(gI2,2)))) + 2*Conj(Kappa)*Lambdax*(KroneckerDelta(1,
      gO2)*(vS*ZA(gI1,0)*ZA(gI2,2) + ZA(gI1,2)*(vS*ZA(gI2,0) - vd*ZA(gI2,2))) +
      KroneckerDelta(0,gO2)*(vS*ZA(gI1,1)*ZA(gI2,2) + ZA(gI1,2)*(vS*ZA(gI2,1) - vu
      *ZA(gI2,2))) - KroneckerDelta(2,gO2)*(ZA(gI1,2)*(vu*ZA(gI2,0) + vd*ZA(gI2,1)
      ) + ZA(gI1,1)*(-(vS*ZA(gI2,0)) + vd*ZA(gI2,2)) + ZA(gI1,0)*(-(vS*ZA(gI2,1))
      + vu*ZA(gI2,2)))) + 1.4142135623730951*(-(KroneckerDelta(2,gO2)*(TLambdax*(
      ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + 2*(Conj(TKappa) - TKappa)*ZA(
      gI1,2)*ZA(gI2,2))) - TLambdax*(KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1
      )*ZA(gI2,2))) + Conj(TLambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,1)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,
      0)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(
      gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(5*KroneckerDelta(2,gO2)*(-1.4142135623730951*(Conj(TLambdax)*
      (ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + TLambdax*(ZA(gI2,1)*ZH(gI1,0)
      + ZA(gI2,0)*ZH(gI1,1)) - 2*(Conj(TKappa) + TKappa)*ZA(gI2,2)*ZH(gI1,2)) + 2*
      Conj(Lambdax)*(-(ZA(gI2,2)*((vu*Kappa + 2*vd*Lambdax)*ZH(gI1,0) + (vd*Kappa
      + 2*vu*Lambdax)*ZH(gI1,1))) + Kappa*ZA(gI2,1)*(vS*ZH(gI1,0) + vd*ZH(gI1,2))
      + Kappa*ZA(gI2,0)*(vS*ZH(gI1,1) + vu*ZH(gI1,2))) + 2*Conj(Kappa)*(Lambdax*ZA
      (gI2,1)*(vS*ZH(gI1,0) + vd*ZH(gI1,2)) + Lambdax*ZA(gI2,0)*(vS*ZH(gI1,1) + vu
      *ZH(gI1,2)) - ZA(gI2,2)*(vu*Lambdax*ZH(gI1,0) + vd*Lambdax*ZH(gI1,1) + 4*vS*
      Kappa*ZH(gI1,2)))) + KroneckerDelta(1,gO2)*(ZA(gI2,1)*(vd*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(
      gI1,1) - 20*vS*AbsSqr(Lambdax)*ZH(gI1,2)) + 5*(-1.4142135623730951*(Conj(
      TLambdax) + TLambdax)*(ZA(gI2,2)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)) + 2*Conj(
      Lambdax)*Kappa*(-(vS*ZA(gI2,0)*ZH(gI1,2)) + ZA(gI2,2)*(vS*ZH(gI1,0) + vd*ZH(
      gI1,2))) + 2*Conj(Kappa)*Lambdax*(-(vS*ZA(gI2,0)*ZH(gI1,2)) + ZA(gI2,2)*(vS*
      ZH(gI1,0) + vd*ZH(gI1,2))))) - KroneckerDelta(0,gO2)*(ZA(gI2,0)*(vd*(3*Sqr(
      g1) + 5*Sqr(g2))*ZH(gI1,0) + vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))
      *ZH(gI1,1) + 20*vS*AbsSqr(Lambdax)*ZH(gI1,2)) - 5*(-1.4142135623730951*(Conj
      (TLambdax) + TLambdax)*(ZA(gI2,2)*ZH(gI1,1) + ZA(gI2,1)*ZH(gI1,2)) + 2*Conj(
      Lambdax)*Kappa*(-(vS*ZA(gI2,1)*ZH(gI1,2)) + ZA(gI2,2)*(vS*ZH(gI1,1) + vu*ZH(
      gI1,2))) + 2*Conj(Kappa)*Lambdax*(-(vS*ZA(gI2,1)*ZH(gI1,2)) + ZA(gI2,2)*(vS*
      ZH(gI1,1) + vu*ZH(gI1,2))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(2*Conj(Lambdax)*Kappa*(-(
      KroneckerDelta(1,gO2)*(vS*ZH(gI1,0)*ZH(gI2,2) + ZH(gI1,2)*(vS*ZH(gI2,0) + vd
      *ZH(gI2,2)))) + KroneckerDelta(2,gO2)*(ZH(gI1,2)*(vu*ZH(gI2,0) + vd*ZH(gI2,1
      )) + ZH(gI1,1)*(vS*ZH(gI2,0) + vd*ZH(gI2,2)) + ZH(gI1,0)*(vS*ZH(gI2,1) + vu*
      ZH(gI2,2))) - KroneckerDelta(0,gO2)*(vS*ZH(gI1,1)*ZH(gI2,2) + ZH(gI1,2)*(vS*
      ZH(gI2,1) + vu*ZH(gI2,2)))) - 2*Conj(Kappa)*Lambdax*(-(KroneckerDelta(1,gO2)
      *(vS*ZH(gI1,0)*ZH(gI2,2) + ZH(gI1,2)*(vS*ZH(gI2,0) + vd*ZH(gI2,2)))) +
      KroneckerDelta(2,gO2)*(ZH(gI1,2)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + ZH(gI1,1)*(
      vS*ZH(gI2,0) + vd*ZH(gI2,2)) + ZH(gI1,0)*(vS*ZH(gI2,1) + vu*ZH(gI2,2))) -
      KroneckerDelta(0,gO2)*(vS*ZH(gI1,1)*ZH(gI2,2) + ZH(gI1,2)*(vS*ZH(gI2,1) + vu
      *ZH(gI2,2)))) + 1.4142135623730951*(KroneckerDelta(2,gO2)*(TLambdax*(ZH(gI1,
      1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + 2*(Conj(TKappa) - TKappa)*ZH(gI1,2)*ZH
      (gI2,2)) + TLambdax*(KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*
      ZH(gI2,2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2
      ))) - Conj(TLambdax)*(KroneckerDelta(2,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)
      *ZH(gI2,1)) + KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,
      2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1236;
   std::complex<double> tmp_1237;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1238;
      std::complex<double> tmp_1239;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1239 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1238 += tmp_1239;
      tmp_1237 += (ZDL(gI1,j2)) * tmp_1238;
   }
   tmp_1236 += tmp_1237;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1236;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1240;
   std::complex<double> tmp_1241;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1242;
      std::complex<double> tmp_1243;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1243 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1242 += tmp_1243;
      tmp_1241 += (Conj(ZDL(gI2,j2))) * tmp_1242;
   }
   tmp_1240 += tmp_1241;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1240;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1244;
   std::complex<double> tmp_1245;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1246;
      std::complex<double> tmp_1247;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1247 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1246 += tmp_1247;
      tmp_1245 += (ZEL(gI1,j2)) * tmp_1246;
   }
   tmp_1244 += tmp_1245;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_1244;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1248;
   std::complex<double> tmp_1249;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1250;
      std::complex<double> tmp_1251;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1251 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1250 += tmp_1251;
      tmp_1249 += (Conj(ZEL(gI2,j2))) * tmp_1250;
   }
   tmp_1248 += tmp_1249;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_1248;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1252;
   std::complex<double> tmp_1253;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1254;
      std::complex<double> tmp_1255;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1255 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1254 += tmp_1255;
      tmp_1253 += (ZUL(gI1,j2)) * tmp_1254;
   }
   tmp_1252 += tmp_1253;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_1252;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1256;
   std::complex<double> tmp_1257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1258;
      std::complex<double> tmp_1259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1259 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1258 += tmp_1259;
      tmp_1257 += (Conj(ZUL(gI2,j2))) * tmp_1258;
   }
   tmp_1256 += tmp_1257;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_1256;

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(-7.0710678118654755*KroneckerDelta(2,
      gO2)*(Conj(Lambdax)*(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,3)) - 2*Conj(
      Kappa)*ZN(gI1,4)*ZN(gI2,4)) - KroneckerDelta(1,gO2)*(ZN(gI1,3)*(
      3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (3.872983346207417*g1*ZN(
      gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*(ZN(
      gI1,4)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,4))) + KroneckerDelta(0,gO2)*(ZN(gI1,2)*
      (3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + 3.872983346207417*g1*ZN(
      gI1,0)*ZN(gI2,2) - 5*(g2*ZN(gI1,1)*ZN(gI2,2) + 1.4142135623730951*Conj(
      Lambdax)*(ZN(gI1,4)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,4)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*
      KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))
      *KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta
      (1,gO1) - 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) +
      3.872983346207417*g1*Conj(ZN(gI1,0))*(-(Conj(ZN(gI2,2))*KroneckerDelta(0,gO1
      )) + Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) - 14.142135623730951*Conj(ZN(gI1
      ,4))*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)*Kappa + 7.0710678118654755*Conj(
      ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1)*Lambdax +
      7.0710678118654755*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*KroneckerDelta(0,gO1)*
      Lambdax + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*KroneckerDelta(
      1,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,3))*Conj(ZN(gI2,2))*
      KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(-3.872983346207417*g1*Conj(
      ZN(gI2,0))*KroneckerDelta(0,gO1) + 5*(g2*Conj(ZN(gI2,1))*KroneckerDelta(0,
      gO1) + 1.4142135623730951*(Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) + Conj(ZN(
      gI2,3))*KroneckerDelta(2,gO1))*Lambdax)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1260;
   std::complex<double> tmp_1261;
   std::complex<double> tmp_1262;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1262 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1261 += tmp_1262;
   tmp_1260 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1261;
   std::complex<double> tmp_1263;
   std::complex<double> tmp_1264;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1264 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1263 += tmp_1264;
   tmp_1260 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1263;
   std::complex<double> tmp_1265;
   std::complex<double> tmp_1266;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1266 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1265 += tmp_1266;
   tmp_1260 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1265;
   std::complex<double> tmp_1267;
   std::complex<double> tmp_1268;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1268 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1267 += tmp_1268;
   tmp_1260 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1267;
   std::complex<double> tmp_1269;
   std::complex<double> tmp_1270;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1270 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1269 += tmp_1270;
   tmp_1260 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1269;
   std::complex<double> tmp_1271;
   std::complex<double> tmp_1272;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1272 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1271 += tmp_1272;
   tmp_1260 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1271;
   std::complex<double> tmp_1273;
   std::complex<double> tmp_1274;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1275;
      std::complex<double> tmp_1276;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1276 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1275 += tmp_1276;
      tmp_1274 += (Conj(ZD(gI2,j2))) * tmp_1275;
   }
   tmp_1273 += tmp_1274;
   tmp_1260 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1273;
   std::complex<double> tmp_1277;
   std::complex<double> tmp_1278;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1279;
      std::complex<double> tmp_1280;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1280 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1279 += tmp_1280;
      tmp_1278 += (Conj(ZD(gI2,j2))) * tmp_1279;
   }
   tmp_1277 += tmp_1278;
   tmp_1260 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1277;
   std::complex<double> tmp_1281;
   std::complex<double> tmp_1282;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1283;
      std::complex<double> tmp_1284;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1284 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1283 += tmp_1284;
      tmp_1282 += (ZD(gI1,j2)) * tmp_1283;
   }
   tmp_1281 += tmp_1282;
   tmp_1260 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1281;
   std::complex<double> tmp_1285;
   std::complex<double> tmp_1286;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1287;
      std::complex<double> tmp_1288;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1288 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1287 += tmp_1288;
      tmp_1286 += (ZD(gI1,j2)) * tmp_1287;
   }
   tmp_1285 += tmp_1286;
   tmp_1260 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1285;
   std::complex<double> tmp_1289;
   std::complex<double> tmp_1290;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1291;
      std::complex<double> tmp_1292;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1293;
         std::complex<double> tmp_1294;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1294 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1293 += tmp_1294;
         tmp_1292 += (ZD(gI1,3 + j2)) * tmp_1293;
      }
      tmp_1291 += tmp_1292;
      tmp_1290 += (Conj(ZD(gI2,3 + j3))) * tmp_1291;
   }
   tmp_1289 += tmp_1290;
   tmp_1260 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1289;
   std::complex<double> tmp_1295;
   std::complex<double> tmp_1296;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1297;
      std::complex<double> tmp_1298;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1299;
         std::complex<double> tmp_1300;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1300 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1299 += tmp_1300;
         tmp_1298 += (Conj(ZD(gI2,j2))) * tmp_1299;
      }
      tmp_1297 += tmp_1298;
      tmp_1296 += (ZD(gI1,j3)) * tmp_1297;
   }
   tmp_1295 += tmp_1296;
   tmp_1260 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1295;
   result += (std::complex<double>(0,-1)) * tmp_1260;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1301;
   std::complex<double> tmp_1302;
   std::complex<double> tmp_1303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1303 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1302 += tmp_1303;
   tmp_1301 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1302;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1305 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1304 += tmp_1305;
   tmp_1301 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1304;
   std::complex<double> tmp_1306;
   std::complex<double> tmp_1307;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1307 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1306 += tmp_1307;
   tmp_1301 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1306;
   std::complex<double> tmp_1308;
   std::complex<double> tmp_1309;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1309 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1308 += tmp_1309;
   tmp_1301 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1308;
   std::complex<double> tmp_1310;
   std::complex<double> tmp_1311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1311 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1310 += tmp_1311;
   tmp_1301 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1310;
   std::complex<double> tmp_1312;
   std::complex<double> tmp_1313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1313 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1312 += tmp_1313;
   tmp_1301 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1312;
   std::complex<double> tmp_1314;
   std::complex<double> tmp_1315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1316;
      std::complex<double> tmp_1317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1317 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1316 += tmp_1317;
      tmp_1315 += (Conj(ZE(gI2,j2))) * tmp_1316;
   }
   tmp_1314 += tmp_1315;
   tmp_1301 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1314;
   std::complex<double> tmp_1318;
   std::complex<double> tmp_1319;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1320;
      std::complex<double> tmp_1321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1321 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1320 += tmp_1321;
      tmp_1319 += (Conj(ZE(gI2,j2))) * tmp_1320;
   }
   tmp_1318 += tmp_1319;
   tmp_1301 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1318;
   std::complex<double> tmp_1322;
   std::complex<double> tmp_1323;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1324;
      std::complex<double> tmp_1325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1325 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1324 += tmp_1325;
      tmp_1323 += (ZE(gI1,j2)) * tmp_1324;
   }
   tmp_1322 += tmp_1323;
   tmp_1301 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1322;
   std::complex<double> tmp_1326;
   std::complex<double> tmp_1327;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1328;
      std::complex<double> tmp_1329;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1329 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1328 += tmp_1329;
      tmp_1327 += (ZE(gI1,j2)) * tmp_1328;
   }
   tmp_1326 += tmp_1327;
   tmp_1301 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1326;
   std::complex<double> tmp_1330;
   std::complex<double> tmp_1331;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1332;
      std::complex<double> tmp_1333;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1334;
         std::complex<double> tmp_1335;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1335 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1334 += tmp_1335;
         tmp_1333 += (ZE(gI1,3 + j2)) * tmp_1334;
      }
      tmp_1332 += tmp_1333;
      tmp_1331 += (Conj(ZE(gI2,3 + j3))) * tmp_1332;
   }
   tmp_1330 += tmp_1331;
   tmp_1301 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1330;
   std::complex<double> tmp_1336;
   std::complex<double> tmp_1337;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1338;
      std::complex<double> tmp_1339;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1340;
         std::complex<double> tmp_1341;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1341 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1340 += tmp_1341;
         tmp_1339 += (Conj(ZE(gI2,j2))) * tmp_1340;
      }
      tmp_1338 += tmp_1339;
      tmp_1337 += (ZE(gI1,j3)) * tmp_1338;
   }
   tmp_1336 += tmp_1337;
   tmp_1301 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1336;
   result += (std::complex<double>(0,-1)) * tmp_1301;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1342;
   std::complex<double> tmp_1343;
   std::complex<double> tmp_1344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1344 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1343 += tmp_1344;
   tmp_1342 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1343;
   std::complex<double> tmp_1345;
   std::complex<double> tmp_1346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1346 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1345 += tmp_1346;
   tmp_1342 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1345;
   std::complex<double> tmp_1347;
   std::complex<double> tmp_1348;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1348 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1347 += tmp_1348;
   tmp_1342 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1347;
   std::complex<double> tmp_1349;
   std::complex<double> tmp_1350;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1350 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1349 += tmp_1350;
   tmp_1342 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1349;
   std::complex<double> tmp_1351;
   std::complex<double> tmp_1352;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1352 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1351 += tmp_1352;
   tmp_1342 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1351;
   std::complex<double> tmp_1353;
   std::complex<double> tmp_1354;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1354 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1353 += tmp_1354;
   tmp_1342 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1353;
   std::complex<double> tmp_1355;
   std::complex<double> tmp_1356;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1357;
      std::complex<double> tmp_1358;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1358 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1357 += tmp_1358;
      tmp_1356 += (Conj(ZU(gI2,j2))) * tmp_1357;
   }
   tmp_1355 += tmp_1356;
   tmp_1342 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1355;
   std::complex<double> tmp_1359;
   std::complex<double> tmp_1360;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1361;
      std::complex<double> tmp_1362;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1362 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1361 += tmp_1362;
      tmp_1360 += (Conj(ZU(gI2,j2))) * tmp_1361;
   }
   tmp_1359 += tmp_1360;
   tmp_1342 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1359;
   std::complex<double> tmp_1363;
   std::complex<double> tmp_1364;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1365;
      std::complex<double> tmp_1366;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1366 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1365 += tmp_1366;
      tmp_1364 += (ZU(gI1,j2)) * tmp_1365;
   }
   tmp_1363 += tmp_1364;
   tmp_1342 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1363;
   std::complex<double> tmp_1367;
   std::complex<double> tmp_1368;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1369;
      std::complex<double> tmp_1370;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1370 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1369 += tmp_1370;
      tmp_1368 += (ZU(gI1,j2)) * tmp_1369;
   }
   tmp_1367 += tmp_1368;
   tmp_1342 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1367;
   std::complex<double> tmp_1371;
   std::complex<double> tmp_1372;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1373;
      std::complex<double> tmp_1374;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1375;
         std::complex<double> tmp_1376;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1376 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1375 += tmp_1376;
         tmp_1374 += (ZU(gI1,3 + j2)) * tmp_1375;
      }
      tmp_1373 += tmp_1374;
      tmp_1372 += (Conj(ZU(gI2,3 + j3))) * tmp_1373;
   }
   tmp_1371 += tmp_1372;
   tmp_1342 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1371;
   std::complex<double> tmp_1377;
   std::complex<double> tmp_1378;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1379;
      std::complex<double> tmp_1380;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1381;
         std::complex<double> tmp_1382;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1382 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1381 += tmp_1382;
         tmp_1380 += (Conj(ZU(gI2,j2))) * tmp_1381;
      }
      tmp_1379 += tmp_1380;
      tmp_1378 += (ZU(gI1,j3)) * tmp_1379;
   }
   tmp_1377 += tmp_1378;
   tmp_1342 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1377;
   result += (std::complex<double>(0,-1)) * tmp_1342;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1383;
   std::complex<double> tmp_1384;
   std::complex<double> tmp_1385;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1386;
      std::complex<double> tmp_1387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1387 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1386 += tmp_1387;
      tmp_1385 += (Conj(ZD(gI2,j2))) * tmp_1386;
   }
   tmp_1384 += tmp_1385;
   tmp_1383 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_1384;
   std::complex<double> tmp_1388;
   std::complex<double> tmp_1389;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1390;
      std::complex<double> tmp_1391;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1391 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1390 += tmp_1391;
      tmp_1389 += (Conj(ZD(gI2,j2))) * tmp_1390;
   }
   tmp_1388 += tmp_1389;
   tmp_1383 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1388;
   std::complex<double> tmp_1392;
   std::complex<double> tmp_1393;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1394;
      std::complex<double> tmp_1395;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1395 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1394 += tmp_1395;
      tmp_1393 += (Conj(ZD(gI2,j2))) * tmp_1394;
   }
   tmp_1392 += tmp_1393;
   tmp_1383 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1392;
   std::complex<double> tmp_1396;
   std::complex<double> tmp_1397;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1398;
      std::complex<double> tmp_1399;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1399 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1398 += tmp_1399;
      tmp_1397 += (ZD(gI1,j2)) * tmp_1398;
   }
   tmp_1396 += tmp_1397;
   tmp_1383 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_1396;
   std::complex<double> tmp_1400;
   std::complex<double> tmp_1401;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1402;
      std::complex<double> tmp_1403;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1403 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1402 += tmp_1403;
      tmp_1401 += (ZD(gI1,j2)) * tmp_1402;
   }
   tmp_1400 += tmp_1401;
   tmp_1383 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_1400;
   std::complex<double> tmp_1404;
   std::complex<double> tmp_1405;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1406;
      std::complex<double> tmp_1407;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1407 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1406 += tmp_1407;
      tmp_1405 += (ZD(gI1,j2)) * tmp_1406;
   }
   tmp_1404 += tmp_1405;
   tmp_1383 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1404;
   result += (std::complex<double>(0,-1)) * tmp_1383;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1408;
   std::complex<double> tmp_1409;
   std::complex<double> tmp_1410;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1411;
      std::complex<double> tmp_1412;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1412 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1411 += tmp_1412;
      tmp_1410 += (Conj(ZE(gI2,j2))) * tmp_1411;
   }
   tmp_1409 += tmp_1410;
   tmp_1408 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_1409;
   std::complex<double> tmp_1413;
   std::complex<double> tmp_1414;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1415;
      std::complex<double> tmp_1416;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1416 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1415 += tmp_1416;
      tmp_1414 += (Conj(ZE(gI2,j2))) * tmp_1415;
   }
   tmp_1413 += tmp_1414;
   tmp_1408 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1413;
   std::complex<double> tmp_1417;
   std::complex<double> tmp_1418;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1419;
      std::complex<double> tmp_1420;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1420 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1419 += tmp_1420;
      tmp_1418 += (Conj(ZE(gI2,j2))) * tmp_1419;
   }
   tmp_1417 += tmp_1418;
   tmp_1408 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1417;
   std::complex<double> tmp_1421;
   std::complex<double> tmp_1422;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1423;
      std::complex<double> tmp_1424;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1424 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1423 += tmp_1424;
      tmp_1422 += (ZE(gI1,j2)) * tmp_1423;
   }
   tmp_1421 += tmp_1422;
   tmp_1408 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_1421;
   std::complex<double> tmp_1425;
   std::complex<double> tmp_1426;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1427;
      std::complex<double> tmp_1428;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1428 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1427 += tmp_1428;
      tmp_1426 += (ZE(gI1,j2)) * tmp_1427;
   }
   tmp_1425 += tmp_1426;
   tmp_1408 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_1425;
   std::complex<double> tmp_1429;
   std::complex<double> tmp_1430;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1431;
      std::complex<double> tmp_1432;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1432 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1431 += tmp_1432;
      tmp_1430 += (ZE(gI1,j2)) * tmp_1431;
   }
   tmp_1429 += tmp_1430;
   tmp_1408 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1429;
   result += (std::complex<double>(0,-1)) * tmp_1408;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1433;
   std::complex<double> tmp_1434;
   std::complex<double> tmp_1435;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1436;
      std::complex<double> tmp_1437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1437 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1436 += tmp_1437;
      tmp_1435 += (Conj(ZU(gI2,j2))) * tmp_1436;
   }
   tmp_1434 += tmp_1435;
   tmp_1433 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(0,gO2)) * tmp_1434;
   std::complex<double> tmp_1438;
   std::complex<double> tmp_1439;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1440;
      std::complex<double> tmp_1441;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1441 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1440 += tmp_1441;
      tmp_1439 += (Conj(ZU(gI2,j2))) * tmp_1440;
   }
   tmp_1438 += tmp_1439;
   tmp_1433 += (0.5*vd*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1438;
   std::complex<double> tmp_1442;
   std::complex<double> tmp_1443;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1444;
      std::complex<double> tmp_1445;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1445 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1444 += tmp_1445;
      tmp_1443 += (Conj(ZU(gI2,j2))) * tmp_1444;
   }
   tmp_1442 += tmp_1443;
   tmp_1433 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1442;
   std::complex<double> tmp_1446;
   std::complex<double> tmp_1447;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1448;
      std::complex<double> tmp_1449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1449 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1448 += tmp_1449;
      tmp_1447 += (ZU(gI1,j2)) * tmp_1448;
   }
   tmp_1446 += tmp_1447;
   tmp_1433 += (-0.5*vS*KroneckerDelta(0,gO2)*Lambdax) * tmp_1446;
   std::complex<double> tmp_1450;
   std::complex<double> tmp_1451;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1452;
      std::complex<double> tmp_1453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1453 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1452 += tmp_1453;
      tmp_1451 += (ZU(gI1,j2)) * tmp_1452;
   }
   tmp_1450 += tmp_1451;
   tmp_1433 += (-0.5*vd*KroneckerDelta(2,gO2)*Lambdax) * tmp_1450;
   std::complex<double> tmp_1454;
   std::complex<double> tmp_1455;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1456;
      std::complex<double> tmp_1457;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1457 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1456 += tmp_1457;
      tmp_1455 += (ZU(gI1,j2)) * tmp_1456;
   }
   tmp_1454 += tmp_1455;
   tmp_1433 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1454;
   result += (std::complex<double>(0,-1)) * tmp_1433;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZH(gI2,0) -
      KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVP(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*KroneckerDelta(0,gO2) -
      vu*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*(vd*KroneckerDelta(0,gO2) - vu*
      KroneckerDelta(1,gO2))*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargWmCgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*(5*g2
      *Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgWmCbargZ(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*(5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargZgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*(5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgZbargWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*(5*g2
      *Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,0) + KroneckerDelta(1,gO2
      )*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - 2*(3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(1,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0
      )*ZP(gI2,1) + KroneckerDelta(0,gO2)*(-2*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP
      (gI2,0) + (-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(KroneckerDelta(1,gO2)*(vu*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZA(gI2,0) + vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,1)
      + 2*(2*vS*Conj(Kappa)*Lambdax - 1.4142135623730951*TLambdax)*ZA(gI2,2))*ZP(
      gI1,0) - KroneckerDelta(0,gO2)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,0)
      + vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,1) + 2*(-1.4142135623730951*Conj(
      TLambdax) + 2*vS*Conj(Lambdax)*Kappa)*ZA(gI2,2))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(1,gO2)*(10*ZH(gI2,2)*(2*vS*Conj(Kappa)*
      Lambdax*ZP(gI1,0) + 1.4142135623730951*TLambdax*ZP(gI1,0) + 2*vS*AbsSqr(
      Lambdax)*ZP(gI1,1)) + ZH(gI2,0)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,
      0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)) + ZH(gI2,1)*(5*vd*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)))) -
      KroneckerDelta(0,gO2)*(ZH(gI2,1)*((-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZP(gI1,0) +
      5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) + ZH(gI2,0)*(vd*(3*Sqr(g1) +
      5*Sqr(g2))*ZP(gI1,0) + 5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) + 10*
      ZH(gI2,2)*(1.4142135623730951*Conj(TLambdax)*ZP(gI1,1) + 2*vS*Conj(Lambdax)*
      (Lambdax*ZP(gI1,0) + Kappa*ZP(gI1,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*
      Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) -
      20*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + 5*KroneckerDelta(0,gO2)*((-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZA(gI1,1)*ZA(gI2,0) + (-2*AbsSqr(Lambdax) + Sqr(
      g2))*ZA(gI1,0)*ZA(gI2,1) + 4*Conj(Lambdax)*Kappa*ZA(gI1,2)*ZA(gI2,2))) -
      KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1
      ,0)*ZA(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 20*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2)) - 5*KroneckerDelta(1,gO2)*((-2*AbsSqr(Lambdax)
      + Sqr(g2))*ZA(gI1,1)*ZA(gI2,0) + (-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI1,0)*
      ZA(gI2,1) + 4*Conj(Kappa)*Lambdax*ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1458;
   tmp_1458 += std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1458 += std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1458 += std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1458 += std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   std::complex<double> tmp_1459;
   std::complex<double> tmp_1460;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1461;
      std::complex<double> tmp_1462;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1463;
         std::complex<double> tmp_1464;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1464 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1463 += tmp_1464;
         tmp_1462 += (Conj(ZV(gI2,j2))) * tmp_1463;
      }
      tmp_1461 += tmp_1462;
      tmp_1460 += (ZV(gI1,j3)) * tmp_1461;
   }
   tmp_1459 += tmp_1460;
   tmp_1458 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1459;
   result += (std::complex<double>(0,-1)) * tmp_1458;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*
      Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) -
      20*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) - 5*KroneckerDelta(0,gO2)*((-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZH(gI1,1)*ZH(gI2,0) + (-2*AbsSqr(Lambdax) + Sqr(
      g2))*ZH(gI1,0)*ZH(gI2,1) + 4*Conj(Lambdax)*Kappa*ZH(gI1,2)*ZH(gI2,2))) -
      KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1
      ,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 20*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + 5*KroneckerDelta(1,gO2)*((-2*AbsSqr(Lambdax)
      + Sqr(g2))*ZH(gI1,1)*ZH(gI2,0) + (-2*AbsSqr(Lambdax) + Sqr(g2))*ZH(gI1,0)*
      ZH(gI2,1) + 4*Conj(Kappa)*Lambdax*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1465;
   std::complex<double> tmp_1466;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1467;
      std::complex<double> tmp_1468;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1468 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1467 += tmp_1468;
      tmp_1466 += (ZUL(gI1,j2)) * tmp_1467;
   }
   tmp_1465 += tmp_1466;
   result += (KroneckerDelta(0,gO2)) * tmp_1465;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1469;
   std::complex<double> tmp_1470;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1471;
      std::complex<double> tmp_1472;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1472 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1471 += tmp_1472;
      tmp_1470 += (Conj(ZDL(gI2,j2))) * tmp_1471;
   }
   tmp_1469 += tmp_1470;
   result += (KroneckerDelta(1,gO1)) * tmp_1469;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1473;
   std::complex<double> tmp_1474;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1474 += Conj(Ye(j1,gI1))*ZER(gI2,j1);
   }
   tmp_1473 += tmp_1474;
   result += (KroneckerDelta(0,gO2)) * tmp_1473;

   return result;
}

double CLASSNAME::CpconjUHpmbarFvFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1475;
   std::complex<double> tmp_1476;
   std::complex<double> tmp_1477;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1477 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1476 += tmp_1477;
   tmp_1475 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1476;
   std::complex<double> tmp_1478;
   std::complex<double> tmp_1479;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1479 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1478 += tmp_1479;
   tmp_1475 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1478;
   std::complex<double> tmp_1480;
   std::complex<double> tmp_1481;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1482;
      std::complex<double> tmp_1483;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1483 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1482 += tmp_1483;
      tmp_1481 += (ZV(gI1,j2)) * tmp_1482;
   }
   tmp_1480 += tmp_1481;
   tmp_1475 += (std::complex<double>(0.,0.7071067811865475)*vS*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_1480;
   std::complex<double> tmp_1484;
   std::complex<double> tmp_1485;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1486;
      std::complex<double> tmp_1487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1487 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1486 += tmp_1487;
      tmp_1485 += (ZV(gI1,j2)) * tmp_1486;
   }
   tmp_1484 += tmp_1485;
   tmp_1475 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1484;
   std::complex<double> tmp_1488;
   std::complex<double> tmp_1489;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1490;
      std::complex<double> tmp_1491;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1492;
         std::complex<double> tmp_1493;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1493 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1492 += tmp_1493;
         tmp_1491 += (Conj(ZE(gI2,j2))) * tmp_1492;
      }
      tmp_1490 += tmp_1491;
      tmp_1489 += (ZV(gI1,j3)) * tmp_1490;
   }
   tmp_1488 += tmp_1489;
   tmp_1475 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1488;
   result += (std::complex<double>(0,-1)) * tmp_1475;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*KroneckerDelta(1,gO2)*(UP(gI2,1)*(5.477225575051661*g1*ZN(gI1,
      0) + 7.0710678118654755*g2*ZN(gI1,1)) + 10*g2*UP(gI2,0)*ZN(gI1,3)) - Conj(
      Lambdax)*KroneckerDelta(0,gO2)*UP(gI2,1)*ZN(gI1,4);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*Conj(ZN(gI1,2))*KroneckerDelta(0,gO1)) + Conj(
      UM(gI2,1))*(0.5477225575051661*g1*Conj(ZN(gI1,0))*KroneckerDelta(0,gO1) +
      0.7071067811865475*g2*Conj(ZN(gI1,1))*KroneckerDelta(0,gO1) - Conj(ZN(gI1,4)
      )*KroneckerDelta(1,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1494;
   std::complex<double> tmp_1495;
   std::complex<double> tmp_1496;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1496 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1495 += tmp_1496;
   tmp_1494 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1495;
   std::complex<double> tmp_1497;
   std::complex<double> tmp_1498;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1498 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1497 += tmp_1498;
   tmp_1494 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1497;
   std::complex<double> tmp_1499;
   std::complex<double> tmp_1500;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1500 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1499 += tmp_1500;
   tmp_1494 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1499;
   std::complex<double> tmp_1501;
   std::complex<double> tmp_1502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1502 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1501 += tmp_1502;
   tmp_1494 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1501;
   std::complex<double> tmp_1503;
   std::complex<double> tmp_1504;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1504 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1503 += tmp_1504;
   tmp_1494 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1503;
   std::complex<double> tmp_1505;
   std::complex<double> tmp_1506;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1506 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1505 += tmp_1506;
   tmp_1494 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1505;
   std::complex<double> tmp_1507;
   std::complex<double> tmp_1508;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1509;
      std::complex<double> tmp_1510;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1511;
         std::complex<double> tmp_1512;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1512 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1511 += tmp_1512;
         tmp_1510 += (ZD(gI1,3 + j2)) * tmp_1511;
      }
      tmp_1509 += tmp_1510;
      tmp_1508 += (Conj(ZD(gI2,3 + j3))) * tmp_1509;
   }
   tmp_1507 += tmp_1508;
   tmp_1494 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1507;
   std::complex<double> tmp_1513;
   std::complex<double> tmp_1514;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1515;
      std::complex<double> tmp_1516;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1517;
         std::complex<double> tmp_1518;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1518 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1517 += tmp_1518;
         tmp_1516 += (Conj(ZD(gI2,j2))) * tmp_1517;
      }
      tmp_1515 += tmp_1516;
      tmp_1514 += (ZD(gI1,j3)) * tmp_1515;
   }
   tmp_1513 += tmp_1514;
   tmp_1494 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1513;
   result += (std::complex<double>(0,-1)) * tmp_1494;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1519;
   std::complex<double> tmp_1520;
   std::complex<double> tmp_1521;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1521 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1520 += tmp_1521;
   tmp_1519 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1520;
   std::complex<double> tmp_1522;
   std::complex<double> tmp_1523;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1523 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1522 += tmp_1523;
   tmp_1519 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1522;
   std::complex<double> tmp_1524;
   std::complex<double> tmp_1525;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1525 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1524 += tmp_1525;
   tmp_1519 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1524;
   std::complex<double> tmp_1526;
   std::complex<double> tmp_1527;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1527 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1526 += tmp_1527;
   tmp_1519 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1526;
   std::complex<double> tmp_1528;
   std::complex<double> tmp_1529;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1529 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1528 += tmp_1529;
   tmp_1519 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1528;
   std::complex<double> tmp_1530;
   std::complex<double> tmp_1531;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1531 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1530 += tmp_1531;
   tmp_1519 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1530;
   std::complex<double> tmp_1532;
   std::complex<double> tmp_1533;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1534;
      std::complex<double> tmp_1535;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1536;
         std::complex<double> tmp_1537;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1537 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1536 += tmp_1537;
         tmp_1535 += (ZE(gI1,3 + j2)) * tmp_1536;
      }
      tmp_1534 += tmp_1535;
      tmp_1533 += (Conj(ZE(gI2,3 + j3))) * tmp_1534;
   }
   tmp_1532 += tmp_1533;
   tmp_1519 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1532;
   result += (std::complex<double>(0,-1)) * tmp_1519;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1538;
   std::complex<double> tmp_1539;
   std::complex<double> tmp_1540;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1540 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1539 += tmp_1540;
   tmp_1538 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1539;
   std::complex<double> tmp_1541;
   std::complex<double> tmp_1542;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1542 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1541 += tmp_1542;
   tmp_1538 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1541;
   std::complex<double> tmp_1543;
   std::complex<double> tmp_1544;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1544 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1543 += tmp_1544;
   tmp_1538 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1543;
   std::complex<double> tmp_1545;
   std::complex<double> tmp_1546;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1546 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1545 += tmp_1546;
   tmp_1538 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1545;
   std::complex<double> tmp_1547;
   std::complex<double> tmp_1548;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1548 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1547 += tmp_1548;
   tmp_1538 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1547;
   std::complex<double> tmp_1549;
   std::complex<double> tmp_1550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1550 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1549 += tmp_1550;
   tmp_1538 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1549;
   std::complex<double> tmp_1551;
   std::complex<double> tmp_1552;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1553;
      std::complex<double> tmp_1554;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1555;
         std::complex<double> tmp_1556;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1556 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1555 += tmp_1556;
         tmp_1554 += (ZU(gI1,3 + j2)) * tmp_1555;
      }
      tmp_1553 += tmp_1554;
      tmp_1552 += (Conj(ZU(gI2,3 + j3))) * tmp_1553;
   }
   tmp_1551 += tmp_1552;
   tmp_1538 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1551;
   std::complex<double> tmp_1557;
   std::complex<double> tmp_1558;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1559;
      std::complex<double> tmp_1560;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1561;
         std::complex<double> tmp_1562;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1562 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1561 += tmp_1562;
         tmp_1560 += (Conj(ZU(gI2,j2))) * tmp_1561;
      }
      tmp_1559 += tmp_1560;
      tmp_1558 += (ZU(gI1,j3)) * tmp_1559;
   }
   tmp_1557 += tmp_1558;
   tmp_1538 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1557;
   result += (std::complex<double>(0,-1)) * tmp_1538;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1563;
   std::complex<double> tmp_1564;
   std::complex<double> tmp_1565;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1565 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1564 += tmp_1565;
   tmp_1563 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_1564;
   std::complex<double> tmp_1566;
   std::complex<double> tmp_1567;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1567 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1566 += tmp_1567;
   tmp_1563 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_1566;
   std::complex<double> tmp_1568;
   std::complex<double> tmp_1569;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1570;
      std::complex<double> tmp_1571;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1571 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1570 += tmp_1571;
      tmp_1569 += (Conj(ZD(gI2,j2))) * tmp_1570;
   }
   tmp_1568 += tmp_1569;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vS*Conj(Lambdax)*
      KroneckerDelta(0,gO2)) * tmp_1568;
   std::complex<double> tmp_1572;
   std::complex<double> tmp_1573;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1574;
      std::complex<double> tmp_1575;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1575 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1574 += tmp_1575;
      tmp_1573 += (Conj(ZD(gI2,j2))) * tmp_1574;
   }
   tmp_1572 += tmp_1573;
   tmp_1563 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1572;
   std::complex<double> tmp_1576;
   std::complex<double> tmp_1577;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1578;
      std::complex<double> tmp_1579;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1579 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1578 += tmp_1579;
      tmp_1577 += (ZU(gI1,j2)) * tmp_1578;
   }
   tmp_1576 += tmp_1577;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vS*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_1576;
   std::complex<double> tmp_1580;
   std::complex<double> tmp_1581;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1582;
      std::complex<double> tmp_1583;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1583 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1582 += tmp_1583;
      tmp_1581 += (ZU(gI1,j2)) * tmp_1582;
   }
   tmp_1580 += tmp_1581;
   tmp_1563 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1580;
   std::complex<double> tmp_1584;
   std::complex<double> tmp_1585;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1586;
      std::complex<double> tmp_1587;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1588;
         std::complex<double> tmp_1589;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1589 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1588 += tmp_1589;
         tmp_1587 += (ZU(gI1,3 + j2)) * tmp_1588;
      }
      tmp_1586 += tmp_1587;
      tmp_1585 += (Conj(ZD(gI2,3 + j3))) * tmp_1586;
   }
   tmp_1584 += tmp_1585;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_1584;
   std::complex<double> tmp_1590;
   std::complex<double> tmp_1591;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1592;
      std::complex<double> tmp_1593;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1594;
         std::complex<double> tmp_1595;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1595 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1594 += tmp_1595;
         tmp_1593 += (ZU(gI1,3 + j2)) * tmp_1594;
      }
      tmp_1592 += tmp_1593;
      tmp_1591 += (Conj(ZD(gI2,3 + j3))) * tmp_1592;
   }
   tmp_1590 += tmp_1591;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_1590;
   std::complex<double> tmp_1596;
   std::complex<double> tmp_1597;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1598;
      std::complex<double> tmp_1599;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1600;
         std::complex<double> tmp_1601;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1601 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1600 += tmp_1601;
         tmp_1599 += (Conj(ZD(gI2,j2))) * tmp_1600;
      }
      tmp_1598 += tmp_1599;
      tmp_1597 += (ZU(gI1,j3)) * tmp_1598;
   }
   tmp_1596 += tmp_1597;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_1596;
   std::complex<double> tmp_1602;
   std::complex<double> tmp_1603;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1604;
      std::complex<double> tmp_1605;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1606;
         std::complex<double> tmp_1607;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1607 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1606 += tmp_1607;
         tmp_1605 += (Conj(ZD(gI2,j2))) * tmp_1606;
      }
      tmp_1604 += tmp_1605;
      tmp_1603 += (ZU(gI1,j3)) * tmp_1604;
   }
   tmp_1602 += tmp_1603;
   tmp_1563 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_1602;
   result += (std::complex<double>(0,-1)) * tmp_1563;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVPHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.3872983346207417*g1*Cos(ThetaW())*ZP(gI2,gO2);
   }
   if (gI2 < 2) {
      result += -0.5*g2*Sin(ThetaW())*ZP(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.5*g2*Cos(ThetaW())*ZP(gI2,gO2);
   }
   if (gI2 < 2) {
      result += 0.3872983346207417*g1*Sin(ThetaW())*ZP(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZA(gI2,0) +
      KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) - KroneckerDelta(1,gO2)*ZH(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpVGbargGgG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpVGbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGVGconjSdSd(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGconjSuSu(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGconjSdSd(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVGconjSuSu(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVGGluGluPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGGluGluPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   double result = 0.0;

   result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   double result = 0.0;

   result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVPbargWmgWm() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbargWmCgWmC() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + Cos(2*
      ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1
      ,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpVPconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Sin(ThetaW())*UM(gI1,0) + Conj(UM(gI2,1))
      *(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Sin(ThetaW())*UP(gI2,0) + Conj(UP(gI1,1))
      *(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UP(gI2,1));

   return result;
}

double CLASSNAME::CpVPbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1608;
   std::complex<double> tmp_1609;
   std::complex<double> tmp_1610;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1610 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1609 += tmp_1610;
   tmp_1608 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1609;
   std::complex<double> tmp_1611;
   std::complex<double> tmp_1612;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1612 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1611 += tmp_1612;
   tmp_1608 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1611;
   std::complex<double> tmp_1613;
   std::complex<double> tmp_1614;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1614 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1613 += tmp_1614;
   tmp_1608 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1613;
   std::complex<double> tmp_1615;
   std::complex<double> tmp_1616;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1616 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1615 += tmp_1616;
   tmp_1608 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1615;
   result += (std::complex<double>(0,-1)) * tmp_1608;

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1617;
   std::complex<double> tmp_1618;
   std::complex<double> tmp_1619;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1619 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1618 += tmp_1619;
   tmp_1617 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Cos(ThetaW()))) *
      tmp_1618;
   std::complex<double> tmp_1620;
   std::complex<double> tmp_1621;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1621 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1620 += tmp_1621;
   tmp_1617 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Cos(ThetaW()))) *
      tmp_1620;
   std::complex<double> tmp_1622;
   std::complex<double> tmp_1623;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1623 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1622 += tmp_1623;
   tmp_1617 += (std::complex<double>(0.,0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1622;
   std::complex<double> tmp_1624;
   std::complex<double> tmp_1625;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1625 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1624 += tmp_1625;
   tmp_1617 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1624;
   result += (std::complex<double>(0,-1)) * tmp_1617;

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1626;
   std::complex<double> tmp_1627;
   std::complex<double> tmp_1628;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1628 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1627 += tmp_1628;
   tmp_1626 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1627;
   std::complex<double> tmp_1629;
   std::complex<double> tmp_1630;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1630 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1629 += tmp_1630;
   tmp_1626 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaW()))) * tmp_1629;
   std::complex<double> tmp_1631;
   std::complex<double> tmp_1632;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1632 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1631 += tmp_1632;
   tmp_1626 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1631;
   std::complex<double> tmp_1633;
   std::complex<double> tmp_1634;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1634 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1633 += tmp_1634;
   tmp_1626 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Sin(ThetaW()))) *
      tmp_1633;
   result += (std::complex<double>(0,-1)) * tmp_1626;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1635;
   std::complex<double> tmp_1636;
   std::complex<double> tmp_1637;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1637 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1636 += tmp_1637;
   tmp_1635 += (-1.5491933384829668*g1*Cos(ThetaW())) * tmp_1636;
   std::complex<double> tmp_1638;
   std::complex<double> tmp_1639;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1639 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1638 += tmp_1639;
   tmp_1635 += (0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW())) *
      tmp_1638;
   result += (0.16666666666666666) * tmp_1635;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1640;
   std::complex<double> tmp_1641;
   std::complex<double> tmp_1642;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1642 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1641 += tmp_1642;
   tmp_1640 += (-1.5491933384829668*g1*Cos(ThetaW())) * tmp_1641;
   std::complex<double> tmp_1643;
   std::complex<double> tmp_1644;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1644 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1643 += tmp_1644;
   tmp_1640 += (-0.7745966692414834*g1*Cos(ThetaW()) - g2*Sin(ThetaW())) *
      tmp_1643;
   result += (0.5) * tmp_1640;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1645;
   std::complex<double> tmp_1646;
   std::complex<double> tmp_1647;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1647 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1646 += tmp_1647;
   tmp_1645 += (3.0983866769659336*g1*Cos(ThetaW())) * tmp_1646;
   std::complex<double> tmp_1648;
   std::complex<double> tmp_1649;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1649 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1648 += tmp_1649;
   tmp_1645 += (0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW())) *
      tmp_1648;
   result += (0.16666666666666666) * tmp_1645;

   return result;
}

std::complex<double> CLASSNAME::CpVPconjVWmHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1
      ));

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWmgWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWmCgWmC() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZP(gI1,0)*ZP(gI2,0) + ZP(
      gI1,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.1*KroneckerDelta(gI1,gI2)*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Cos(ThetaW())*UM(gI1,0) + Conj(UM(gI2,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Cos(ThetaW())*UP(gI2,0) + Conj(UP(gI1,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZA(gI1,0)*ZA(gI2,0) + ZA(
      gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*KroneckerDelta(gI1,gI2)*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZH(gI1,0)*ZH(gI2,0) + ZH(
      gI1,1)*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(gI2,1)*ZH(gI1,1));

   return result;
}

double CLASSNAME::CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj
      (ZN(gI2,2))*ZN(gI1,2) - Conj(ZN(gI2,3))*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj(
      ZN(gI1,2))*ZN(gI2,2) - Conj(ZN(gI1,3))*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1650;
   std::complex<double> tmp_1651;
   std::complex<double> tmp_1652;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1652 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1651 += tmp_1652;
   tmp_1650 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1651;
   std::complex<double> tmp_1653;
   std::complex<double> tmp_1654;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1654 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1653 += tmp_1654;
   tmp_1650 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1653;
   std::complex<double> tmp_1655;
   std::complex<double> tmp_1656;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1656 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1655 += tmp_1656;
   tmp_1650 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1655;
   std::complex<double> tmp_1657;
   std::complex<double> tmp_1658;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1658 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1657 += tmp_1658;
   tmp_1650 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1657;
   result += (std::complex<double>(0,-1)) * tmp_1650;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1659;
   std::complex<double> tmp_1660;
   std::complex<double> tmp_1661;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1661 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1660 += tmp_1661;
   tmp_1659 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1660;
   std::complex<double> tmp_1662;
   std::complex<double> tmp_1663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1663 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1662 += tmp_1663;
   tmp_1659 += (std::complex<double>(0.,-0.7745966692414834)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1662;
   std::complex<double> tmp_1664;
   std::complex<double> tmp_1665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1665 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1664 += tmp_1665;
   tmp_1659 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1664;
   std::complex<double> tmp_1666;
   std::complex<double> tmp_1667;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1667 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1666 += tmp_1667;
   tmp_1659 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1666;
   result += (std::complex<double>(0,-1)) * tmp_1659;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1668;
   std::complex<double> tmp_1669;
   std::complex<double> tmp_1670;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1670 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1669 += tmp_1670;
   tmp_1668 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1669;
   std::complex<double> tmp_1671;
   std::complex<double> tmp_1672;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1672 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1671 += tmp_1672;
   tmp_1668 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())) * tmp_1671;
   std::complex<double> tmp_1673;
   std::complex<double> tmp_1674;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1674 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1673 += tmp_1674;
   tmp_1668 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1673;
   std::complex<double> tmp_1675;
   std::complex<double> tmp_1676;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1676 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1675 += tmp_1676;
   tmp_1668 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1675;
   result += (std::complex<double>(0,-1)) * tmp_1668;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1677;
   std::complex<double> tmp_1678;
   std::complex<double> tmp_1679;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1679 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1678 += tmp_1679;
   tmp_1677 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1678;
   std::complex<double> tmp_1680;
   std::complex<double> tmp_1681;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1681 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1680 += tmp_1681;
   tmp_1677 += (-3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1680;
   result += (0.16666666666666666) * tmp_1677;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1682;
   std::complex<double> tmp_1683;
   std::complex<double> tmp_1684;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1684 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1683 += tmp_1684;
   tmp_1682 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1683;
   std::complex<double> tmp_1685;
   std::complex<double> tmp_1686;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1686 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1685 += tmp_1686;
   tmp_1682 += (-(g2*Cos(ThetaW())) + 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1685;
   result += (0.5) * tmp_1682;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1687;
   std::complex<double> tmp_1688;
   std::complex<double> tmp_1689;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1689 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1688 += tmp_1689;
   tmp_1687 += (-3.0983866769659336*g1*Sin(ThetaW())) * tmp_1688;
   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1691 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1690 += tmp_1691;
   tmp_1687 += (3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1690;
   result += (0.16666666666666666) * tmp_1687;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjVWmHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Sin(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(vd
      *ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargPgWm() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmbargWmCgP() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargWmCgZ() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargZgWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVWmVP() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVZVWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(ZA(gI2,0)*ZP(gI1,0) + ZA(gI2,1)*ZP
      (gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1));

   return result;
}

double CLASSNAME::CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1692;
   std::complex<double> tmp_1693;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1693 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1692 += tmp_1693;
   result += (-0.7071067811865475*g2) * tmp_1692;

   return result;
}

double CLASSNAME::CpconjVWmbarFuFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFvFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gI1));
   }

   return result;
}

double CLASSNAME::CpconjVWmbarFvFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1694;
   std::complex<double> tmp_1695;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1695 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1694 += tmp_1695;
   result += (0.7071067811865475*g2) * tmp_1694;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,1) + 1.4142135623730951*Conj(UM(
      gI2,1))*ZN(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(ZN(gI1,1))*UP(gI2,0)) + 0.7071067811865475*g2*Conj(ZN(gI1
      ,3))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1696;
   std::complex<double> tmp_1697;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1697 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1696 += tmp_1697;
   result += (0.5*Sqr(g2)) * tmp_1696;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1698;
   std::complex<double> tmp_1699;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1699 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1698 += tmp_1699;
   result += (0.5*Sqr(g2)) * tmp_1698;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1700;
   std::complex<double> tmp_1701;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1701 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1700 += tmp_1701;
   result += (0.5*Sqr(g2)) * tmp_1700;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1702;
   std::complex<double> tmp_1703;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1703 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1702 += tmp_1703;
   result += (0.7071067811865475*g2) * tmp_1702;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVPHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Sin(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm1() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm2() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm3() const
{
   double result = 0.0;

   result = 2*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaHpmPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(3,gO2)*ZP(gI2,1)) - 0.1*Conj(UP
      (gI1,1))*(10*KroneckerDelta(4,gO2)*Lambdax*ZP(gI2,0) + 1.4142135623730951*(
      3.872983346207417*g1*KroneckerDelta(0,gO2) + 5*g2*KroneckerDelta(1,gO2))*ZP(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaHpmPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(2,gO1)*UM(gI1,0)*ZP(gI2,0)) + 0.1*UM(gI1,1)*(
      5.477225575051661*g1*KroneckerDelta(0,gO1)*ZP(gI2,0) + 7.0710678118654755*g2
      *KroneckerDelta(1,gO1)*ZP(gI2,0) - 10*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZP
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*KroneckerDelta(2,gO2)*ZP(gI1,0)) + Conj(UM(gI2
      ,1))*(0.5477225575051661*g1*KroneckerDelta(0,gO2)*ZP(gI1,0) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZP(gI1,0) - KroneckerDelta(4,gO2
      )*Lambdax*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(Lambdax)*KroneckerDelta(4,gO1)*UP(gI2,1)*ZP(gI1,0)) - 0.1*(
      10*g2*KroneckerDelta(3,gO1)*UP(gI2,0) + 1.4142135623730951*(
      3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*KroneckerDelta(1,gO1))*UP(
      gI2,1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaVWmPR(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(1,gO2)*UM(gI1,0) + 0.7071067811865475*g2*
      KroneckerDelta(2,gO2)*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaVWmPL(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = g2*Conj(UP(gI1,0))*KroneckerDelta(1,gO1) - 0.7071067811865475*g2*
      Conj(UP(gI1,1))*KroneckerDelta(3,gO1);

   return result;
}

double CLASSNAME::CpUChibarFvSvPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFvSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZV(gI2,gI1))*KroneckerDelta(0,gO1
         );
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZV(gI2,gI1))*KroneckerDelta(1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5477225575051661*g1*KroneckerDelta(0,gO2)*ZV(gI1,gI2);
   }
   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZV(gI1,gI2);
   }

   return result;
}

double CLASSNAME::CpUChiconjSvFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-5*g2*Conj(ZN(gI2,1))*KroneckerDelta(2,gO2)*ZH(gI1,0) +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,0) +
      7.0710678118654755*Conj(ZN(gI2,3))*KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,0)
      - 3.872983346207417*g1*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2)*ZH(gI1,1) + 5*
      g2*Conj(ZN(gI2,3))*KroneckerDelta(1,gO2)*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,1))*
      KroneckerDelta(3,gO2)*ZH(gI1,1) + 7.0710678118654755*Conj(ZN(gI2,4))*
      KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1) + 3.872983346207417*g1*Conj(ZN(gI2,0
      ))*(KroneckerDelta(2,gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1)) -
      14.142135623730951*Conj(ZN(gI2,4))*KroneckerDelta(4,gO2)*Kappa*ZH(gI1,2) +
      7.0710678118654755*Conj(ZN(gI2,3))*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,2) +
      Conj(ZN(gI2,2))*(3.872983346207417*g1*KroneckerDelta(0,gO2)*ZH(gI1,0) - 5*
      g2*KroneckerDelta(1,gO2)*ZH(gI1,0) + 7.0710678118654755*Lambdax*(
      KroneckerDelta(4,gO2)*ZH(gI1,1) + KroneckerDelta(3,gO2)*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN(gI2,2)
      - 5*g2*KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN(gI2,2) + 7.0710678118654755*Conj(
      Lambdax)*KroneckerDelta(4,gO1)*ZH(gI1,1)*ZN(gI2,2) + 7.0710678118654755*Conj
      (Lambdax)*KroneckerDelta(4,gO1)*ZH(gI1,0)*ZN(gI2,3) - 3.872983346207417*g1*
      KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN(gI2,3) + 5*g2*KroneckerDelta(1,gO1)*ZH(
      gI1,1)*ZN(gI2,3) - 14.142135623730951*Conj(Kappa)*KroneckerDelta(4,gO1)*ZH(
      gI1,2)*ZN(gI2,4) + KroneckerDelta(3,gO1)*(ZH(gI1,1)*(-3.872983346207417*g1*
      ZN(gI2,0) + 5*g2*ZN(gI2,1)) + 7.0710678118654755*Conj(Lambdax)*(ZH(gI1,2)*ZN
      (gI2,2) + ZH(gI1,0)*ZN(gI2,4))) + KroneckerDelta(2,gO1)*(ZH(gI1,0)*(
      3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + 7.0710678118654755*Conj(
      Lambdax)*(ZH(gI1,2)*ZN(gI2,3) + ZH(gI1,1)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFdSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1704;
   std::complex<double> tmp_1705;
   std::complex<double> tmp_1706;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1706 += Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1));
   }
   tmp_1705 += tmp_1706;
   tmp_1704 += (std::complex<double>(0.,-0.3651483716701107)*g1*KroneckerDelta(
      0,gO2)) * tmp_1705;
   std::complex<double> tmp_1707;
   std::complex<double> tmp_1708;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1709;
      std::complex<double> tmp_1710;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1710 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1709 += tmp_1710;
      tmp_1708 += (Conj(ZD(gI2,j2))) * tmp_1709;
   }
   tmp_1707 += tmp_1708;
   tmp_1704 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1707;
   result += (std::complex<double>(0,-1)) * tmp_1704;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFdSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1711;
   std::complex<double> tmp_1712;
   std::complex<double> tmp_1713;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1713 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1712 += tmp_1713;
   tmp_1711 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO1)) * tmp_1712;
   std::complex<double> tmp_1714;
   std::complex<double> tmp_1715;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1715 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1714 += tmp_1715;
   tmp_1711 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO1)) * tmp_1714;
   std::complex<double> tmp_1716;
   std::complex<double> tmp_1717;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1718;
      std::complex<double> tmp_1719;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1719 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1718 += tmp_1719;
      tmp_1717 += (ZDL(gI1,j2)) * tmp_1718;
   }
   tmp_1716 += tmp_1717;
   tmp_1711 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1716;
   result += (std::complex<double>(0,-1)) * tmp_1711;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFeSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1720;
   std::complex<double> tmp_1721;
   std::complex<double> tmp_1722;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1722 += Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1));
   }
   tmp_1721 += tmp_1722;
   tmp_1720 += (std::complex<double>(0.,-1.0954451150103321)*g1*KroneckerDelta(
      0,gO2)) * tmp_1721;
   std::complex<double> tmp_1723;
   std::complex<double> tmp_1724;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1725;
      std::complex<double> tmp_1726;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1726 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1725 += tmp_1726;
      tmp_1724 += (Conj(ZE(gI2,j2))) * tmp_1725;
   }
   tmp_1723 += tmp_1724;
   tmp_1720 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1723;
   result += (std::complex<double>(0,-1)) * tmp_1720;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFeSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1727;
   std::complex<double> tmp_1728;
   std::complex<double> tmp_1729;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1729 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1728 += tmp_1729;
   tmp_1727 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(0
      ,gO1)) * tmp_1728;
   std::complex<double> tmp_1730;
   std::complex<double> tmp_1731;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1731 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1730 += tmp_1731;
   tmp_1727 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO1)) * tmp_1730;
   std::complex<double> tmp_1732;
   std::complex<double> tmp_1733;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1734;
      std::complex<double> tmp_1735;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1735 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1734 += tmp_1735;
      tmp_1733 += (ZEL(gI1,j2)) * tmp_1734;
   }
   tmp_1732 += tmp_1733;
   tmp_1727 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1732;
   result += (std::complex<double>(0,-1)) * tmp_1727;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFuSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1736;
   std::complex<double> tmp_1737;
   std::complex<double> tmp_1738;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1738 += Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1));
   }
   tmp_1737 += tmp_1738;
   tmp_1736 += (std::complex<double>(0.,0.7302967433402214)*g1*KroneckerDelta(0
      ,gO2)) * tmp_1737;
   std::complex<double> tmp_1739;
   std::complex<double> tmp_1740;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1741;
      std::complex<double> tmp_1742;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1742 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1741 += tmp_1742;
      tmp_1740 += (Conj(ZU(gI2,j2))) * tmp_1741;
   }
   tmp_1739 += tmp_1740;
   tmp_1736 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_1739;
   result += (std::complex<double>(0,-1)) * tmp_1736;

   return result;
}

std::complex<double> CLASSNAME::CpUChibarFuSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   std::complex<double> tmp_1745;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1745 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1744 += tmp_1745;
   tmp_1743 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO1)) * tmp_1744;
   std::complex<double> tmp_1746;
   std::complex<double> tmp_1747;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1747 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1746 += tmp_1747;
   tmp_1743 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      1,gO1)) * tmp_1746;
   std::complex<double> tmp_1748;
   std::complex<double> tmp_1749;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1750;
      std::complex<double> tmp_1751;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1751 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1750 += tmp_1751;
      tmp_1749 += (ZUL(gI1,j2)) * tmp_1750;
   }
   tmp_1748 += tmp_1749;
   tmp_1743 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_1748;
   result += (std::complex<double>(0,-1)) * tmp_1743;

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(gI1,1))*KroneckerDelta(2,
      gO2)*ZA(gI2,0) + 7.0710678118654755*Conj(ZN(gI1,4))*KroneckerDelta(3,gO2)*
      Lambdax*ZA(gI2,0) + 7.0710678118654755*Conj(ZN(gI1,3))*KroneckerDelta(4,gO2)
      *Lambdax*ZA(gI2,0) + 3.872983346207417*g1*Conj(ZN(gI1,3))*KroneckerDelta(0,
      gO2)*ZA(gI2,1) - 5*g2*Conj(ZN(gI1,3))*KroneckerDelta(1,gO2)*ZA(gI2,1) - 5*g2
      *Conj(ZN(gI1,1))*KroneckerDelta(3,gO2)*ZA(gI2,1) + 7.0710678118654755*Conj(
      ZN(gI1,4))*KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,1) + 3.872983346207417*g1*
      Conj(ZN(gI1,0))*(-(KroneckerDelta(2,gO2)*ZA(gI2,0)) + KroneckerDelta(3,gO2)*
      ZA(gI2,1)) - 14.142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(4,gO2)*Kappa*
      ZA(gI2,2) + 7.0710678118654755*Conj(ZN(gI1,3))*KroneckerDelta(2,gO2)*Lambdax
      *ZA(gI2,2) + Conj(ZN(gI1,2))*(-3.872983346207417*g1*KroneckerDelta(0,gO2)*ZA
      (gI2,0) + 5*(g2*KroneckerDelta(1,gO2)*ZA(gI2,0) + 1.4142135623730951*Lambdax
      *(KroneckerDelta(4,gO2)*ZA(gI2,1) + KroneckerDelta(3,gO2)*ZA(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(3.872983346207417*g1*KroneckerDelta(0,
      gO1)*ZA(gI2,0)*ZN(gI1,2) - 5*g2*KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN(gI1,2) -
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZA(gI2,1)*ZN(gI1,2) -
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZA(gI2,0)*ZN(gI1,3)
      - 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZA(gI2,1)*ZN(gI1,3) + 5*g2*
      KroneckerDelta(1,gO1)*ZA(gI2,1)*ZN(gI1,3) + 14.142135623730951*Conj(Kappa)*
      KroneckerDelta(4,gO1)*ZA(gI2,2)*ZN(gI1,4) - KroneckerDelta(3,gO1)*(ZA(gI2,1)
      *(3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1)) + 7.0710678118654755*Conj
      (Lambdax)*(ZA(gI2,2)*ZN(gI1,2) + ZA(gI2,0)*ZN(gI1,4))) + KroneckerDelta(2,
      gO1)*(ZA(gI2,0)*(3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1)) -
      7.0710678118654755*Conj(Lambdax)*(ZA(gI2,2)*ZN(gI1,3) + ZA(gI2,1)*ZN(gI1,4))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1752;
   std::complex<double> tmp_1753;
   std::complex<double> tmp_1754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1754 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1753 += tmp_1754;
   tmp_1752 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO2)) * tmp_1753;
   std::complex<double> tmp_1755;
   std::complex<double> tmp_1756;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1756 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1755 += tmp_1756;
   tmp_1752 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO2)) * tmp_1755;
   std::complex<double> tmp_1757;
   std::complex<double> tmp_1758;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1759;
      std::complex<double> tmp_1760;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1760 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1759 += tmp_1760;
      tmp_1758 += (Conj(ZDL(gI2,j2))) * tmp_1759;
   }
   tmp_1757 += tmp_1758;
   tmp_1752 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1757;
   result += (std::complex<double>(0,-1)) * tmp_1752;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1761;
   std::complex<double> tmp_1762;
   std::complex<double> tmp_1763;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1763 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1762 += tmp_1763;
   tmp_1761 += (std::complex<double>(0.,-0.3651483716701107)*g1*KroneckerDelta(
      0,gO1)) * tmp_1762;
   std::complex<double> tmp_1764;
   std::complex<double> tmp_1765;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1766;
      std::complex<double> tmp_1767;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1767 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1766 += tmp_1767;
      tmp_1765 += (ZD(gI1,j2)) * tmp_1766;
   }
   tmp_1764 += tmp_1765;
   tmp_1761 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1764;
   result += (std::complex<double>(0,-1)) * tmp_1761;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1768;
   std::complex<double> tmp_1769;
   std::complex<double> tmp_1770;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1770 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1769 += tmp_1770;
   tmp_1768 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(0
      ,gO2)) * tmp_1769;
   std::complex<double> tmp_1771;
   std::complex<double> tmp_1772;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1772 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1771 += tmp_1772;
   tmp_1768 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(1
      ,gO2)) * tmp_1771;
   std::complex<double> tmp_1773;
   std::complex<double> tmp_1774;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1775;
      std::complex<double> tmp_1776;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1776 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1775 += tmp_1776;
      tmp_1774 += (Conj(ZEL(gI2,j2))) * tmp_1775;
   }
   tmp_1773 += tmp_1774;
   tmp_1768 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1773;
   result += (std::complex<double>(0,-1)) * tmp_1768;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1777;
   std::complex<double> tmp_1778;
   std::complex<double> tmp_1779;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1779 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_1778 += tmp_1779;
   tmp_1777 += (std::complex<double>(0.,-1.0954451150103321)*g1*KroneckerDelta(
      0,gO1)) * tmp_1778;
   std::complex<double> tmp_1780;
   std::complex<double> tmp_1781;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1782;
      std::complex<double> tmp_1783;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1783 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1782 += tmp_1783;
      tmp_1781 += (ZE(gI1,j2)) * tmp_1782;
   }
   tmp_1780 += tmp_1781;
   tmp_1777 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1780;
   result += (std::complex<double>(0,-1)) * tmp_1777;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1784;
   std::complex<double> tmp_1785;
   std::complex<double> tmp_1786;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1786 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1785 += tmp_1786;
   tmp_1784 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (0,gO2)) * tmp_1785;
   std::complex<double> tmp_1787;
   std::complex<double> tmp_1788;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1788 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1787 += tmp_1788;
   tmp_1784 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      1,gO2)) * tmp_1787;
   std::complex<double> tmp_1789;
   std::complex<double> tmp_1790;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1791;
      std::complex<double> tmp_1792;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1792 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1791 += tmp_1792;
      tmp_1790 += (Conj(ZUL(gI2,j2))) * tmp_1791;
   }
   tmp_1789 += tmp_1790;
   tmp_1784 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_1789;
   result += (std::complex<double>(0,-1)) * tmp_1784;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1793;
   std::complex<double> tmp_1794;
   std::complex<double> tmp_1795;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1795 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1794 += tmp_1795;
   tmp_1793 += (std::complex<double>(0.,0.7302967433402214)*g1*KroneckerDelta(0
      ,gO1)) * tmp_1794;
   std::complex<double> tmp_1796;
   std::complex<double> tmp_1797;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1798;
      std::complex<double> tmp_1799;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1799 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1798 += tmp_1799;
      tmp_1797 += (ZU(gI1,j2)) * tmp_1798;
   }
   tmp_1796 += tmp_1797;
   tmp_1793 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_1796;
   result += (std::complex<double>(0,-1)) * tmp_1793;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(3,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1,gO1) +
      1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(
      KroneckerDelta(2,gO2)*ZN(gI2,2) - KroneckerDelta(3,gO2)*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(Conj(ZN(gI2,2))*KroneckerDelta(2,gO1) - Conj(ZN(gI2,3))*
      KroneckerDelta(3,gO1))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gI1,0))*
      KroneckerDelta(1,gO2)*ZA(gI2,1) + Conj(UM(gI1,1))*(g2*KroneckerDelta(0,gO2)*
      ZA(gI2,0) - KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2*KroneckerDelta(0,
      gO1)*UP(gI1,1)*ZA(gI2,1) + KroneckerDelta(1,gO1)*(g2*UP(gI1,0)*ZA(gI2,0) -
      Conj(Lambdax)*UP(gI1,1)*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(ZN(gI2,4))*KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,0)) - 0.1*(
      10*g2*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2) + 1.4142135623730951*(
      3.872983346207417*g1*Conj(ZN(gI2,0)) + 5*g2*Conj(ZN(gI2,1)))*KroneckerDelta(
      1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO1)*ZN(gI2,2)*ZP(gI1,0)) + KroneckerDelta(1,
      gO1)*(0.5477225575051661*g1*ZN(gI2,0)*ZP(gI1,0) + 0.7071067811865475*g2*ZN(
      gI2,1)*ZP(gI1,0) - Conj(Lambdax)*ZN(gI2,4)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*KroneckerDelta(1,gO2)*ZH(
      gI1,1) + Conj(UM(gI2,1))*(g2*KroneckerDelta(0,gO2)*ZH(gI1,0) +
      KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*KroneckerDelta(0,gO1)*UP(gI2,1)*ZH(gI1,1) +
      KroneckerDelta(1,gO1)*(g2*UP(gI2,0)*ZH(gI1,0) + Conj(Lambdax)*UP(gI2,1)*ZH(
      gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1800;
   std::complex<double> tmp_1801;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1801 += Conj(ZEL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1800 += tmp_1801;
   result += (-(g2*KroneckerDelta(0,gO2))) * tmp_1800;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1802;
   std::complex<double> tmp_1803;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1804;
      std::complex<double> tmp_1805;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1805 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1804 += tmp_1805;
      tmp_1803 += (ZV(gI1,j2)) * tmp_1804;
   }
   tmp_1802 += tmp_1803;
   result += (KroneckerDelta(1,gO1)) * tmp_1802;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1806;
   std::complex<double> tmp_1807;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1808;
      std::complex<double> tmp_1809;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1809 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1808 += tmp_1809;
      tmp_1807 += (Conj(ZD(gI2,j2))) * tmp_1808;
   }
   tmp_1806 += tmp_1807;
   result += (KroneckerDelta(1,gO2)) * tmp_1806;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1810;
   std::complex<double> tmp_1811;
   std::complex<double> tmp_1812;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1812 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1811 += tmp_1812;
   tmp_1810 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_1811
      ;
   std::complex<double> tmp_1813;
   std::complex<double> tmp_1814;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1815;
      std::complex<double> tmp_1816;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1816 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1815 += tmp_1816;
      tmp_1814 += (ZUL(gI1,j2)) * tmp_1815;
   }
   tmp_1813 += tmp_1814;
   tmp_1810 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_1813;
   result += (std::complex<double>(0,-1)) * tmp_1810;

   return result;
}

double CLASSNAME::CpbarUChabarFvSePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1817;
   std::complex<double> tmp_1818;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1818 += Conj(Ye(j1,gI1))*Conj(ZE(gI2,3 + j1));
   }
   tmp_1817 += tmp_1818;
   result += (KroneckerDelta(1,gO1)) * tmp_1817;
   if (gI1 < 3) {
      result += -(g2*Conj(ZE(gI2,gI1))*KroneckerDelta(0,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1819;
   std::complex<double> tmp_1820;
   std::complex<double> tmp_1821;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1821 += Conj(ZDL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1820 += tmp_1821;
   tmp_1819 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_1820
      ;
   std::complex<double> tmp_1822;
   std::complex<double> tmp_1823;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1824;
      std::complex<double> tmp_1825;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1825 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1824 += tmp_1825;
      tmp_1823 += (Conj(ZDL(gI2,j2))) * tmp_1824;
   }
   tmp_1822 += tmp_1823;
   tmp_1819 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1822;
   result += (std::complex<double>(0,-1)) * tmp_1819;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1826;
   std::complex<double> tmp_1827;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1828;
      std::complex<double> tmp_1829;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1829 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1828 += tmp_1829;
      tmp_1827 += (ZU(gI1,j2)) * tmp_1828;
   }
   tmp_1826 += tmp_1827;
   result += (KroneckerDelta(1,gO1)) * tmp_1826;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*KroneckerDelta(0,gO1)*Sin(ThetaW()) + 0.1*Conj(
      UM(gI2,1))*KroneckerDelta(1,gO1)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      )))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*KroneckerDelta(0,gO1) + 0.1*Conj(
      UM(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN(gI2,1))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += Ye(gO2,gI2)*ZP(gI1,0);
   }

   return result;
}

double CLASSNAME::CpbarUFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1830;
      std::complex<double> tmp_1831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1831 += Conj(ZV(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1830 += tmp_1831;
      result += (Conj(UM(gI2,1))) * tmp_1830;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZV(gI1,gO1))*UP(gI2,0));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1832;
      std::complex<double> tmp_1833;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1833 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1832 += tmp_1833;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1832;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1834;
      std::complex<double> tmp_1835;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1835 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_1834 += tmp_1835;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1834;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1836;
      std::complex<double> tmp_1837;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1837 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1836 += tmp_1837;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1836;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1838;
      std::complex<double> tmp_1839;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1839 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_1838 += tmp_1839;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1838;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.0954451150103321*g1*Conj(ZE(gI1,3 + gO2))*Conj(ZN(gI2,0))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1840;
      std::complex<double> tmp_1841;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1841 += Conj(ZE(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1840 += tmp_1841;
      result += (-Conj(ZN(gI2,2))) * tmp_1840;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZE(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1842;
      std::complex<double> tmp_1843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1843 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1842 += tmp_1843;
      result += (-ZN(gI2,2)) * tmp_1842;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*ZER(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZEL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFeVWmFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(gI2,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Sin(ThetaW())*ZER(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZEL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1844;
      std::complex<double> tmp_1845;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1845 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1844 += tmp_1845;
      result += (ZP(gI1,0)) * tmp_1844;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1846;
      std::complex<double> tmp_1847;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1847 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1846 += tmp_1847;
      result += (ZP(gI1,1)) * tmp_1846;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1848;
      std::complex<double> tmp_1849;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1849 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1848 += tmp_1849;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1848;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1850;
      std::complex<double> tmp_1851;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1851 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_1850 += tmp_1851;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1850;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1852;
      std::complex<double> tmp_1853;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1853 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1852 += tmp_1853;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1852;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1854;
      std::complex<double> tmp_1855;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1855 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1854 += tmp_1855;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1854;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1856;
      std::complex<double> tmp_1857;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1857 += Conj(ZU(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1856 += tmp_1857;
      result += (Conj(UM(gI2,1))) * tmp_1856;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZU(gI1,gO1))*UP(gI2,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1858;
      std::complex<double> tmp_1859;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1859 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1858 += tmp_1859;
      result += (UP(gI2,1)) * tmp_1858;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.3651483716701107*g1*Conj(ZD(gI1,3 + gO2))*Conj(ZN(gI2,0))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1860;
      std::complex<double> tmp_1861;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1861 += Conj(ZD(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1860 += tmp_1861;
      result += (-Conj(ZN(gI2,2))) * tmp_1860;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZD(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1862;
      std::complex<double> tmp_1863;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1863 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1862 += tmp_1863;
      result += (-ZN(gI2,2)) * tmp_1862;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*PhaseGlu*Conj(ZD(gI1,3 + gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(PhaseGlu)*Conj(ZD(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*ZDR(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Cos(ThetaW())*ZDR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZDL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFdVWmFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVWmFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZUL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Sin(ThetaW())*ZDR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZDL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1864;
      std::complex<double> tmp_1865;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1865 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1864 += tmp_1865;
      result += (ZP(gI1,1)) * tmp_1864;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1866;
      std::complex<double> tmp_1867;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1867 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1866 += tmp_1867;
      result += (ZP(gI1,0)) * tmp_1866;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1868;
      std::complex<double> tmp_1869;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1869 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1868 += tmp_1869;
      result += (Conj(UP(gI1,1))) * tmp_1868;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZD(gI2,gO1))*UM(gI1,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1870;
      std::complex<double> tmp_1871;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1871 += Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1870 += tmp_1871;
      result += (UM(gI1,1)) * tmp_1870;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1872;
      std::complex<double> tmp_1873;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1873 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1872 += tmp_1873;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_1872;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1874;
      std::complex<double> tmp_1875;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1875 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_1874 += tmp_1875;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_1874;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1876;
      std::complex<double> tmp_1877;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1877 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1876 += tmp_1877;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1876;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1878;
      std::complex<double> tmp_1879;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1879 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1878 += tmp_1879;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1878;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7302967433402214*g1*Conj(ZN(gI2,0))*Conj(ZU(gI1,3 + gO2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1880;
      std::complex<double> tmp_1881;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1881 += Conj(ZU(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1880 += tmp_1881;
      result += (-Conj(ZN(gI2,3))) * tmp_1880;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZU(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZU(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1882;
      std::complex<double> tmp_1883;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1883 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1882 += tmp_1883;
      result += (-ZN(gI2,3)) * tmp_1882;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*PhaseGlu*Conj(ZU(gI1,3 + gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(PhaseGlu)*Conj(ZU(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*ZUR(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(ZUL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Cos(ThetaW())*ZUR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZUL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Sin(ThetaW())*ZUR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZUL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFuconjVWmFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjVWmFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZDL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpGlubarFdSdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1884;
   std::complex<double> tmp_1885;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1885 += Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1));
   }
   tmp_1884 += tmp_1885;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_1884;

   return result;
}

std::complex<double> CLASSNAME::CpGlubarFdSdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1886;
   std::complex<double> tmp_1887;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1887 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1886 += tmp_1887;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1886;

   return result;
}

std::complex<double> CLASSNAME::CpGlubarFuSuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1888;
   std::complex<double> tmp_1889;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1889 += Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1));
   }
   tmp_1888 += tmp_1889;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_1888;

   return result;
}

std::complex<double> CLASSNAME::CpGlubarFuSuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1890;
   std::complex<double> tmp_1891;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1891 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1890 += tmp_1891;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1890;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1892;
   std::complex<double> tmp_1893;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1893 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1892 += tmp_1893;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_1892;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1894;
   std::complex<double> tmp_1895;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1895 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1894 += tmp_1895;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1894;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1896;
   std::complex<double> tmp_1897;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1897 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1896 += tmp_1897;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_1896;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1898;
   std::complex<double> tmp_1899;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1899 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1898 += tmp_1899;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1898;

   return result;
}

std::complex<double> CLASSNAME::CpGluVGGluPR() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3*AbsSqr(PhaseGlu);

   return result;
}

std::complex<double> CLASSNAME::CpGluVGGluPL() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3*AbsSqr(PhaseGlu);

   return result;
}

double CLASSNAME::CpbarFvconjHpmFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvconjHpmFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1900;
   std::complex<double> tmp_1901;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1901 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
   }
   tmp_1900 += tmp_1901;
   result += (ZP(gI1,0)) * tmp_1900;

   return result;
}

double CLASSNAME::CpbarFvbarChaSePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvbarChaSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZE(gI2,gO1))*UM(gI1,0));
   }
   std::complex<double> tmp_1902;
   std::complex<double> tmp_1903;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1903 += Conj(Ye(j1,gO1))*Conj(ZE(gI2,3 + j1));
   }
   tmp_1902 += tmp_1903;
   result += (UM(gI1,1)) * tmp_1902;

   return result;
}

double CLASSNAME::CpbarFvSvChiPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvSvChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZV(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZV(gI1,gO1))*ZN(gI2,1);
   }

   return result;
}

double CLASSNAME::CpbarFvVZFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFvVZFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFvconjVWmFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvconjVWmFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1904;
   std::complex<double> tmp_1905;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1905 += Conj(ZER(gO2,j1))*Ye(j1,gI2);
   }
   tmp_1904 += tmp_1905;
   result += (ZP(gI1,0)) * tmp_1904;

   return result;
}

double CLASSNAME::CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1906;
   std::complex<double> tmp_1907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1908 += tmp_1909;
      tmp_1907 += (Conj(ZV(gI1,j2))) * tmp_1908;
   }
   tmp_1906 += tmp_1907;
   result += (Conj(UM(gI2,1))) * tmp_1906;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1910;
   std::complex<double> tmp_1911;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1911 += Conj(ZV(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_1910 += tmp_1911;
   result += (-(g2*UP(gI2,0))) * tmp_1910;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1912;
   std::complex<double> tmp_1913;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1914;
      std::complex<double> tmp_1915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1915 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1914 += tmp_1915;
      tmp_1913 += (Conj(ZEL(gI1,j2))) * tmp_1914;
   }
   tmp_1912 += tmp_1913;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_1912;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1916;
   std::complex<double> tmp_1917;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1918;
      std::complex<double> tmp_1919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1919 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_1918 += tmp_1919;
      tmp_1917 += (ZEL(gO1,j2)) * tmp_1918;
   }
   tmp_1916 += tmp_1917;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_1916
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1920;
   std::complex<double> tmp_1921;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1922;
      std::complex<double> tmp_1923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1923 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1922 += tmp_1923;
      tmp_1921 += (Conj(ZEL(gI2,j2))) * tmp_1922;
   }
   tmp_1920 += tmp_1921;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1920;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1924;
   std::complex<double> tmp_1925;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1926;
      std::complex<double> tmp_1927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1927 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1926 += tmp_1927;
      tmp_1925 += (ZEL(gO1,j2)) * tmp_1926;
   }
   tmp_1924 += tmp_1925;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1924;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1928;
   std::complex<double> tmp_1929;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1929 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_1928 += tmp_1929;
   result += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_1928;
   std::complex<double> tmp_1930;
   std::complex<double> tmp_1931;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1932;
      std::complex<double> tmp_1933;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1933 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1932 += tmp_1933;
      tmp_1931 += (Conj(ZE(gI1,j2))) * tmp_1932;
   }
   tmp_1930 += tmp_1931;
   result += (-Conj(ZN(gI2,2))) * tmp_1930;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1934;
   std::complex<double> tmp_1935;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1935 += Conj(ZE(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_1934 += tmp_1935;
   result += (0.7071067811865475*(0.7745966692414834*g1*ZN(gI2,0) + g2*ZN(gI2,1
      ))) * tmp_1934;
   std::complex<double> tmp_1936;
   std::complex<double> tmp_1937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1938;
      std::complex<double> tmp_1939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1939 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1938 += tmp_1939;
      tmp_1937 += (ZEL(gO1,j2)) * tmp_1938;
   }
   tmp_1936 += tmp_1937;
   result += (-ZN(gI2,2)) * tmp_1936;

   return result;
}

double CLASSNAME::CpbarFeVWmFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*ZEL(gO1,gI2);
   }

   return result;
}

double CLASSNAME::CpbarFeVZFePR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFeVZFePL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1940;
   std::complex<double> tmp_1941;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1942;
      std::complex<double> tmp_1943;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1943 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1942 += tmp_1943;
      tmp_1941 += (Conj(ZUL(gI2,j2))) * tmp_1942;
   }
   tmp_1940 += tmp_1941;
   result += (ZP(gI1,0)) * tmp_1940;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1944;
   std::complex<double> tmp_1945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1946 += tmp_1947;
      tmp_1945 += (ZDL(gO1,j2)) * tmp_1946;
   }
   tmp_1944 += tmp_1945;
   result += (ZP(gI1,1)) * tmp_1944;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1948;
   std::complex<double> tmp_1949;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1950;
      std::complex<double> tmp_1951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1951 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1950 += tmp_1951;
      tmp_1949 += (Conj(ZDL(gI1,j2))) * tmp_1950;
   }
   tmp_1948 += tmp_1949;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
      tmp_1948;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1952;
   std::complex<double> tmp_1953;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1954;
      std::complex<double> tmp_1955;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1955 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_1954 += tmp_1955;
      tmp_1953 += (ZDL(gO1,j2)) * tmp_1954;
   }
   tmp_1952 += tmp_1953;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_1952
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1956;
   std::complex<double> tmp_1957;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1958;
      std::complex<double> tmp_1959;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1959 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1958 += tmp_1959;
      tmp_1957 += (Conj(ZDL(gI2,j2))) * tmp_1958;
   }
   tmp_1956 += tmp_1957;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1956;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1960;
   std::complex<double> tmp_1961;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1962;
      std::complex<double> tmp_1963;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1963 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1962 += tmp_1963;
      tmp_1961 += (ZDL(gO1,j2)) * tmp_1962;
   }
   tmp_1960 += tmp_1961;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1960;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1964;
   std::complex<double> tmp_1965;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1966;
      std::complex<double> tmp_1967;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1967 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1966 += tmp_1967;
      tmp_1965 += (Conj(ZU(gI1,j2))) * tmp_1966;
   }
   tmp_1964 += tmp_1965;
   result += (Conj(UM(gI2,1))) * tmp_1964;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1968;
   std::complex<double> tmp_1969;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1969 += Conj(ZU(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_1968 += tmp_1969;
   result += (-(g2*UP(gI2,0))) * tmp_1968;
   std::complex<double> tmp_1970;
   std::complex<double> tmp_1971;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1972;
      std::complex<double> tmp_1973;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1973 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1972 += tmp_1973;
      tmp_1971 += (ZDL(gO1,j2)) * tmp_1972;
   }
   tmp_1970 += tmp_1971;
   result += (UP(gI2,1)) * tmp_1970;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1974;
   std::complex<double> tmp_1975;
   std::complex<double> tmp_1976;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1976 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_1975 += tmp_1976;
   tmp_1974 += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_1975;
   std::complex<double> tmp_1977;
   std::complex<double> tmp_1978;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1979;
      std::complex<double> tmp_1980;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1980 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1979 += tmp_1980;
      tmp_1978 += (Conj(ZD(gI1,j2))) * tmp_1979;
   }
   tmp_1977 += tmp_1978;
   tmp_1974 += (-3*Conj(ZN(gI2,2))) * tmp_1977;
   result += (0.3333333333333333) * tmp_1974;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1981;
   std::complex<double> tmp_1982;
   std::complex<double> tmp_1983;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1983 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_1982 += tmp_1983;
   tmp_1981 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) - 3*g2*ZN(
      gI2,1))) * tmp_1982;
   std::complex<double> tmp_1984;
   std::complex<double> tmp_1985;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1986;
      std::complex<double> tmp_1987;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1987 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1986 += tmp_1987;
      tmp_1985 += (ZDL(gO1,j2)) * tmp_1986;
   }
   tmp_1984 += tmp_1985;
   tmp_1981 += (-6*ZN(gI2,2)) * tmp_1984;
   result += (0.16666666666666666) * tmp_1981;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_1988;
   std::complex<double> tmp_1989;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1989 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_1988 += tmp_1989;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_1988;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_1990;
   std::complex<double> tmp_1991;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1991 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_1990 += tmp_1991;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1990;

   return result;
}

double CLASSNAME::CpbarFdVWmFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdVWmFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1992;
   std::complex<double> tmp_1993;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1993 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_1992 += tmp_1993;
   result += (-0.7071067811865475*g2) * tmp_1992;

   return result;
}

double CLASSNAME::CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1994;
   std::complex<double> tmp_1995;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1996;
      std::complex<double> tmp_1997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1997 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_1996 += tmp_1997;
      tmp_1995 += (Conj(ZDL(gI2,j2))) * tmp_1996;
   }
   tmp_1994 += tmp_1995;
   result += (ZP(gI1,1)) * tmp_1994;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1998;
   std::complex<double> tmp_1999;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2000;
      std::complex<double> tmp_2001;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2001 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2000 += tmp_2001;
      tmp_1999 += (ZUL(gO1,j2)) * tmp_2000;
   }
   tmp_1998 += tmp_1999;
   result += (ZP(gI1,0)) * tmp_1998;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2002;
   std::complex<double> tmp_2003;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2004;
      std::complex<double> tmp_2005;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2005 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2004 += tmp_2005;
      tmp_2003 += (Conj(ZD(gI2,j2))) * tmp_2004;
   }
   tmp_2002 += tmp_2003;
   result += (Conj(UP(gI1,1))) * tmp_2002;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2006;
   std::complex<double> tmp_2007;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2007 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2006 += tmp_2007;
   result += (-(g2*UM(gI1,0))) * tmp_2006;
   std::complex<double> tmp_2008;
   std::complex<double> tmp_2009;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2010;
      std::complex<double> tmp_2011;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2011 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2010 += tmp_2011;
      tmp_2009 += (ZUL(gO1,j2)) * tmp_2010;
   }
   tmp_2008 += tmp_2009;
   result += (UM(gI1,1)) * tmp_2008;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2012;
   std::complex<double> tmp_2013;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2014;
      std::complex<double> tmp_2015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2015 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2014 += tmp_2015;
      tmp_2013 += (Conj(ZUL(gI1,j2))) * tmp_2014;
   }
   tmp_2012 += tmp_2013;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
      tmp_2012;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2016;
   std::complex<double> tmp_2017;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2018;
      std::complex<double> tmp_2019;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2019 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2018 += tmp_2019;
      tmp_2017 += (ZUL(gO1,j2)) * tmp_2018;
   }
   tmp_2016 += tmp_2017;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_2016
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2020;
   std::complex<double> tmp_2021;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2022;
      std::complex<double> tmp_2023;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2023 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2022 += tmp_2023;
      tmp_2021 += (Conj(ZUL(gI2,j2))) * tmp_2022;
   }
   tmp_2020 += tmp_2021;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2020;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2024;
   std::complex<double> tmp_2025;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2026;
      std::complex<double> tmp_2027;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2027 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2026 += tmp_2027;
      tmp_2025 += (ZUL(gO1,j2)) * tmp_2026;
   }
   tmp_2024 += tmp_2025;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2024;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2028;
   std::complex<double> tmp_2029;
   std::complex<double> tmp_2030;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2030 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2029 += tmp_2030;
   tmp_2028 += (2.1908902300206643*g1*Conj(ZN(gI2,0))) * tmp_2029;
   std::complex<double> tmp_2031;
   std::complex<double> tmp_2032;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2033;
      std::complex<double> tmp_2034;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2034 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2033 += tmp_2034;
      tmp_2032 += (Conj(ZU(gI1,j2))) * tmp_2033;
   }
   tmp_2031 += tmp_2032;
   tmp_2028 += (-3*Conj(ZN(gI2,3))) * tmp_2031;
   result += (0.3333333333333333) * tmp_2028;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2035;
   std::complex<double> tmp_2036;
   std::complex<double> tmp_2037;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2037 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2036 += tmp_2037;
   tmp_2035 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) + 3*g2*ZN(
      gI2,1))) * tmp_2036;
   std::complex<double> tmp_2038;
   std::complex<double> tmp_2039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2040;
      std::complex<double> tmp_2041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2041 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2040 += tmp_2041;
      tmp_2039 += (ZUL(gO1,j2)) * tmp_2040;
   }
   tmp_2038 += tmp_2039;
   tmp_2035 += (-6*ZN(gI2,3)) * tmp_2038;
   result += (0.16666666666666666) * tmp_2035;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2042;
   std::complex<double> tmp_2043;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2043 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2042 += tmp_2043;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2042;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2044;
   std::complex<double> tmp_2045;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2045 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2044 += tmp_2045;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2044;

   return result;
}

double CLASSNAME::CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI2,gO2);

   return result;
}

double CLASSNAME::CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuconjVWmFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjVWmFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2046;
   std::complex<double> tmp_2047;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2047 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2046 += tmp_2047;
   result += (-0.7071067811865475*g2) * tmp_2046;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2048;
   std::complex<double> tmp_2049;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2049 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2048 += tmp_2049;
   result += (-1) * tmp_2048;
   std::complex<double> tmp_2050;
   std::complex<double> tmp_2051;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2051 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2050 += tmp_2051;
   result += (-0.5) * tmp_2050;
   std::complex<double> tmp_2052;
   std::complex<double> tmp_2053;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2053 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2052 += tmp_2053;
   result += (-1) * tmp_2052;
   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2055 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2054 += tmp_2055;
   result += (-0.5) * tmp_2054;
   std::complex<double> tmp_2056;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2057;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2057 += (Conj(CpconjUSdFuChaPL(gO2,gI1,gI2))*
            CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPR(gO2,gI1,gI2))*
            CpconjUSdFuChaPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha(gI2));
      }
      tmp_2056 += tmp_2057;
   }
   result += tmp_2056;
   std::complex<double> tmp_2058;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2059;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2059 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2058 += tmp_2059;
   }
   result += tmp_2058;
   std::complex<double> tmp_2060;
   std::complex<double> tmp_2061;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2062;
      std::complex<double> tmp_2063;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2063 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2062 += tmp_2063;
      tmp_2061 += (MFd(gI1)) * tmp_2062;
   }
   tmp_2060 += tmp_2061;
   result += (-2) * tmp_2060;
   std::complex<double> tmp_2064;
   std::complex<double> tmp_2065;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2066;
      std::complex<double> tmp_2067;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2067 += B0(p,MFu(gI1),MCha(gI2))*(Conj(CpconjUSdFuChaPR(gO2,
            gI1,gI2))*CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPL(gO2,
            gI1,gI2))*CpconjUSdFuChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2066 += tmp_2067;
      tmp_2065 += (MFu(gI1)) * tmp_2066;
   }
   tmp_2064 += tmp_2065;
   result += (-2) * tmp_2064;
   std::complex<double> tmp_2068;
   std::complex<double> tmp_2069;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2069 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2068 += tmp_2069;
   result += (-1) * tmp_2068;
   std::complex<double> tmp_2070;
   std::complex<double> tmp_2071;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2071 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2070 += tmp_2071;
   result += (-1) * tmp_2070;
   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2073 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2072 += tmp_2073;
   result += (-1) * tmp_2072;
   std::complex<double> tmp_2074;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2075;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2075 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2074 += tmp_2075;
   }
   result += tmp_2074;
   std::complex<double> tmp_2076;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2077;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2077 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2076 += tmp_2077;
   }
   result += tmp_2076;
   std::complex<double> tmp_2078;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2079;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2079 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2078 += tmp_2079;
   }
   result += tmp_2078;
   std::complex<double> tmp_2080;
   std::complex<double> tmp_2081;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2081 += (Conj(CpconjUSdGluFdPL(gO2,0,gI2))*CpconjUSdGluFdPL(gO1,0,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,0,gI2))*CpconjUSdGluFdPR(gO1,0,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2080 += tmp_2081;
   result += (1.3333333333333333) * tmp_2080;
   std::complex<double> tmp_2082;
   std::complex<double> tmp_2083;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2083 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2082 += tmp_2083;
   result += (1.3333333333333333) * tmp_2082;
   std::complex<double> tmp_2084;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2084 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2084;
   std::complex<double> tmp_2085;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2085 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2085;
   std::complex<double> tmp_2086;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2086 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2086;
   std::complex<double> tmp_2087;
   std::complex<double> tmp_2088;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2088 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,0,gI2))*
         CpconjUSdGluFdPL(gO1,0,gI2) + Conj(CpconjUSdGluFdPL(gO2,0,gI2))*
         CpconjUSdGluFdPR(gO1,0,gI2))*MFd(gI2);
   }
   tmp_2087 += tmp_2088;
   result += (-2.6666666666666665*MGlu) * tmp_2087;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2089;
   std::complex<double> tmp_2090;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2090 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2089 += tmp_2090;
   result += (-1) * tmp_2089;
   std::complex<double> tmp_2091;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2092;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2092 += (Conj(CpconjUSvbarChaFePL(gO2,gI1,gI2))*
            CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(CpconjUSvbarChaFePR(gO2,gI1,
            gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFe(gI2));
      }
      tmp_2091 += tmp_2092;
   }
   result += tmp_2091;
   std::complex<double> tmp_2093;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2094;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2094 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2093 += tmp_2094;
   }
   result += tmp_2093;
   std::complex<double> tmp_2095;
   std::complex<double> tmp_2096;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2097;
      std::complex<double> tmp_2098;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2098 += B0(p,MCha(gI1),MFe(gI2))*(Conj(CpconjUSvbarChaFePR(
            gO2,gI1,gI2))*CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(
            CpconjUSvbarChaFePL(gO2,gI1,gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2097 += tmp_2098;
      tmp_2096 += (MCha(gI1)) * tmp_2097;
   }
   tmp_2095 += tmp_2096;
   result += (-2) * tmp_2095;
   std::complex<double> tmp_2099;
   std::complex<double> tmp_2100;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2100 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2099 += tmp_2100;
   result += (-0.5) * tmp_2099;
   std::complex<double> tmp_2101;
   std::complex<double> tmp_2102;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2102 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2101 += tmp_2102;
   result += (-1) * tmp_2101;
   std::complex<double> tmp_2103;
   std::complex<double> tmp_2104;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2104 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2103 += tmp_2104;
   result += (-0.5) * tmp_2103;
   std::complex<double> tmp_2105;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2106;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2106 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2105 += tmp_2106;
   }
   result += tmp_2105;
   std::complex<double> tmp_2107;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2108;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2108 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2107 += tmp_2108;
   }
   result += tmp_2107;
   std::complex<double> tmp_2109;
   std::complex<double> tmp_2110;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2111;
      std::complex<double> tmp_2112;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2112 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2111 += tmp_2112;
      tmp_2110 += (MFv(gI1)) * tmp_2111;
   }
   tmp_2109 += tmp_2110;
   result += (-2) * tmp_2109;
   std::complex<double> tmp_2113;
   std::complex<double> tmp_2114;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2114 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2113 += tmp_2114;
   result += (-3) * tmp_2113;
   std::complex<double> tmp_2115;
   std::complex<double> tmp_2116;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2116 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2115 += tmp_2116;
   result += (-1) * tmp_2115;
   std::complex<double> tmp_2117;
   std::complex<double> tmp_2118;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2118 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2117 += tmp_2118;
   result += (-3) * tmp_2117;
   std::complex<double> tmp_2119;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2119 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2119;
   std::complex<double> tmp_2120;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2120 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2120;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2121;
   std::complex<double> tmp_2122;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2122 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2121 += tmp_2122;
   result += (-1) * tmp_2121;
   std::complex<double> tmp_2123;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2124;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2124 += (Conj(CpconjUSubarChaFdPL(gO2,gI1,gI2))*
            CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(CpconjUSubarChaFdPR(gO2,gI1,
            gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFd(gI2));
      }
      tmp_2123 += tmp_2124;
   }
   result += tmp_2123;
   std::complex<double> tmp_2125;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2126;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2126 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2125 += tmp_2126;
   }
   result += tmp_2125;
   std::complex<double> tmp_2127;
   std::complex<double> tmp_2128;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2129;
      std::complex<double> tmp_2130;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2130 += B0(p,MCha(gI1),MFd(gI2))*(Conj(CpconjUSubarChaFdPR(
            gO2,gI1,gI2))*CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChaFdPL(gO2,gI1,gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2129 += tmp_2130;
      tmp_2128 += (MCha(gI1)) * tmp_2129;
   }
   tmp_2127 += tmp_2128;
   result += (-2) * tmp_2127;
   std::complex<double> tmp_2131;
   std::complex<double> tmp_2132;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2132 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2131 += tmp_2132;
   result += (-0.5) * tmp_2131;
   std::complex<double> tmp_2133;
   std::complex<double> tmp_2134;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2134 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2133 += tmp_2134;
   result += (-1) * tmp_2133;
   std::complex<double> tmp_2135;
   std::complex<double> tmp_2136;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2136 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2135 += tmp_2136;
   result += (-0.5) * tmp_2135;
   std::complex<double> tmp_2137;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2138;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2138 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2137 += tmp_2138;
   }
   result += tmp_2137;
   std::complex<double> tmp_2139;
   std::complex<double> tmp_2140;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2141;
      std::complex<double> tmp_2142;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2142 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2141 += tmp_2142;
      tmp_2140 += (MFu(gI1)) * tmp_2141;
   }
   tmp_2139 += tmp_2140;
   result += (-2) * tmp_2139;
   std::complex<double> tmp_2143;
   std::complex<double> tmp_2144;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2144 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2143 += tmp_2144;
   result += (-1) * tmp_2143;
   std::complex<double> tmp_2145;
   std::complex<double> tmp_2146;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2146 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2145 += tmp_2146;
   result += (-1) * tmp_2145;
   std::complex<double> tmp_2147;
   std::complex<double> tmp_2148;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2148 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2147 += tmp_2148;
   result += (-1) * tmp_2147;
   std::complex<double> tmp_2149;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2150;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2150 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2149 += tmp_2150;
   }
   result += tmp_2149;
   std::complex<double> tmp_2151;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2152;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2152 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2151 += tmp_2152;
   }
   result += tmp_2151;
   std::complex<double> tmp_2153;
   std::complex<double> tmp_2154;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2154 += (Conj(CpconjUSuGluFuPL(gO2,0,gI2))*CpconjUSuGluFuPL(gO1,0,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,0,gI2))*CpconjUSuGluFuPR(gO1,0,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2153 += tmp_2154;
   result += (1.3333333333333333) * tmp_2153;
   std::complex<double> tmp_2155;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2155 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2155;
   std::complex<double> tmp_2156;
   std::complex<double> tmp_2157;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2157 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2156 += tmp_2157;
   result += (1.3333333333333333) * tmp_2156;
   std::complex<double> tmp_2158;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2158 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2158;
   std::complex<double> tmp_2159;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2159 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2159;
   std::complex<double> tmp_2160;
   std::complex<double> tmp_2161;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2161 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,0,gI2))*
         CpconjUSuGluFuPL(gO1,0,gI2) + Conj(CpconjUSuGluFuPL(gO2,0,gI2))*
         CpconjUSuGluFuPR(gO1,0,gI2))*MFu(gI2);
   }
   tmp_2160 += tmp_2161;
   result += (-2.6666666666666665*MGlu) * tmp_2160;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2162;
   std::complex<double> tmp_2163;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2163 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2162 += tmp_2163;
   result += (-1) * tmp_2162;
   std::complex<double> tmp_2164;
   std::complex<double> tmp_2165;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2165 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2164 += tmp_2165;
   result += (-0.5) * tmp_2164;
   std::complex<double> tmp_2166;
   std::complex<double> tmp_2167;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2167 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2166 += tmp_2167;
   result += (-1) * tmp_2166;
   std::complex<double> tmp_2168;
   std::complex<double> tmp_2169;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2169 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2168 += tmp_2169;
   result += (-0.5) * tmp_2168;
   std::complex<double> tmp_2170;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2171;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2171 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_2170 += tmp_2171;
   }
   result += tmp_2170;
   std::complex<double> tmp_2172;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2173;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2173 += (Conj(CpconjUSeFvChaPL(gO2,gI1,gI2))*
            CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPR(gO2,gI1,gI2))*
            CpconjUSeFvChaPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha(gI2));
      }
      tmp_2172 += tmp_2173;
   }
   result += tmp_2172;
   std::complex<double> tmp_2174;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2175;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2175 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_2174 += tmp_2175;
   }
   result += tmp_2174;
   std::complex<double> tmp_2176;
   std::complex<double> tmp_2177;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2178;
      std::complex<double> tmp_2179;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2179 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2178 += tmp_2179;
      tmp_2177 += (MFe(gI1)) * tmp_2178;
   }
   tmp_2176 += tmp_2177;
   result += (-2) * tmp_2176;
   std::complex<double> tmp_2180;
   std::complex<double> tmp_2181;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2182;
      std::complex<double> tmp_2183;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2183 += B0(p,MFv(gI1),MCha(gI2))*(Conj(CpconjUSeFvChaPR(gO2,
            gI1,gI2))*CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPL(gO2,
            gI1,gI2))*CpconjUSeFvChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2182 += tmp_2183;
      tmp_2181 += (MFv(gI1)) * tmp_2182;
   }
   tmp_2180 += tmp_2181;
   result += (-2) * tmp_2180;
   std::complex<double> tmp_2184;
   std::complex<double> tmp_2185;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2185 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2184 += tmp_2185;
   result += (-3) * tmp_2184;
   std::complex<double> tmp_2186;
   std::complex<double> tmp_2187;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2187 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2186 += tmp_2187;
   result += (-1) * tmp_2186;
   std::complex<double> tmp_2188;
   std::complex<double> tmp_2189;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2189 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2188 += tmp_2189;
   result += (-3) * tmp_2188;
   std::complex<double> tmp_2190;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2191;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2191 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2190 += tmp_2191;
   }
   result += tmp_2190;
   std::complex<double> tmp_2192;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2193;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2193 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2192 += tmp_2193;
   }
   result += tmp_2192;
   std::complex<double> tmp_2194;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2194 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2194;
   std::complex<double> tmp_2195;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2195 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2195;
   std::complex<double> tmp_2196;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2196 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2196;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmCgWmC(gO1)*CpUhhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmgWm(gO1)*CpUhhbargWmgWm(gO2));
   result += -(B0(p,MVZ,MVZ)*CpUhhbargZgZ(gO1)*CpUhhbargZgZ(gO2));
   result += 4*B0(p,MVWm,MVWm)*Conj(CpUhhconjVWmVWm(gO2))*CpUhhconjVWmVWm(gO1);
   result += 4*A0(MVWm)*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(p,MVZ,MVZ)*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   std::complex<double> tmp_2197;
   std::complex<double> tmp_2198;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2198 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2197 += tmp_2198;
   result += (-1) * tmp_2197;
   std::complex<double> tmp_2199;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2200;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2200 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2199 += tmp_2200;
   }
   result += tmp_2199;
   std::complex<double> tmp_2201;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2202;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2202 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*
            CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*
            CpUhhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2201 += tmp_2202;
   }
   result += tmp_2201;
   std::complex<double> tmp_2203;
   std::complex<double> tmp_2204;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2205;
      std::complex<double> tmp_2206;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2206 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2
            ,gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2205 += tmp_2206;
      tmp_2204 += (MCha(gI1)) * tmp_2205;
   }
   tmp_2203 += tmp_2204;
   result += (-2) * tmp_2203;
   std::complex<double> tmp_2207;
   std::complex<double> tmp_2208;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2208 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2207 += tmp_2208;
   result += (-0.5) * tmp_2207;
   std::complex<double> tmp_2209;
   std::complex<double> tmp_2210;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2210 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2209 += tmp_2210;
   result += (-1) * tmp_2209;
   std::complex<double> tmp_2211;
   std::complex<double> tmp_2212;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2212 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2211 += tmp_2212;
   result += (-0.5) * tmp_2211;
   std::complex<double> tmp_2213;
   std::complex<double> tmp_2214;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2215;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2215 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2214 += tmp_2215;
   }
   tmp_2213 += tmp_2214;
   result += (0.5) * tmp_2213;
   std::complex<double> tmp_2216;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2217;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2217 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2216 += tmp_2217;
   }
   result += tmp_2216;
   std::complex<double> tmp_2218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2219;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2219 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_2218 += tmp_2219;
   }
   result += tmp_2218;
   std::complex<double> tmp_2220;
   std::complex<double> tmp_2221;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2222;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2222 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2221 += tmp_2222;
   }
   tmp_2220 += tmp_2221;
   result += (0.5) * tmp_2220;
   std::complex<double> tmp_2223;
   std::complex<double> tmp_2224;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2225;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2225 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2224 += tmp_2225;
   }
   tmp_2223 += tmp_2224;
   result += (3) * tmp_2223;
   std::complex<double> tmp_2226;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2227;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2227 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2226 += tmp_2227;
   }
   result += tmp_2226;
   std::complex<double> tmp_2228;
   std::complex<double> tmp_2229;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2230;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2230 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2229 += tmp_2230;
   }
   tmp_2228 += tmp_2229;
   result += (3) * tmp_2228;
   std::complex<double> tmp_2231;
   std::complex<double> tmp_2232;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2233;
      std::complex<double> tmp_2234;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2234 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2233 += tmp_2234;
      tmp_2232 += (MFd(gI1)) * tmp_2233;
   }
   tmp_2231 += tmp_2232;
   result += (-6) * tmp_2231;
   std::complex<double> tmp_2235;
   std::complex<double> tmp_2236;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2237;
      std::complex<double> tmp_2238;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2238 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2237 += tmp_2238;
      tmp_2236 += (MFe(gI1)) * tmp_2237;
   }
   tmp_2235 += tmp_2236;
   result += (-2) * tmp_2235;
   std::complex<double> tmp_2239;
   std::complex<double> tmp_2240;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2241;
      std::complex<double> tmp_2242;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2242 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2241 += tmp_2242;
      tmp_2240 += (MFu(gI1)) * tmp_2241;
   }
   tmp_2239 += tmp_2240;
   result += (-6) * tmp_2239;
   std::complex<double> tmp_2243;
   std::complex<double> tmp_2244;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2245;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2245 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2244 += tmp_2245;
   }
   tmp_2243 += tmp_2244;
   result += (0.5) * tmp_2243;
   std::complex<double> tmp_2246;
   std::complex<double> tmp_2247;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2248;
      std::complex<double> tmp_2249;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2249 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,
            gI1,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2)
            )*CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2248 += tmp_2249;
      tmp_2247 += (MChi(gI1)) * tmp_2248;
   }
   tmp_2246 += tmp_2247;
   result += (-1) * tmp_2246;
   std::complex<double> tmp_2250;
   std::complex<double> tmp_2251;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2251 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2250 += tmp_2251;
   result += (-3) * tmp_2250;
   std::complex<double> tmp_2252;
   std::complex<double> tmp_2253;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2253 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2252 += tmp_2253;
   result += (-1) * tmp_2252;
   std::complex<double> tmp_2254;
   std::complex<double> tmp_2255;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2255 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2254 += tmp_2255;
   result += (-3) * tmp_2254;
   std::complex<double> tmp_2256;
   std::complex<double> tmp_2257;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2258;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2258 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2257 += tmp_2258;
   }
   tmp_2256 += tmp_2257;
   result += (3) * tmp_2256;
   std::complex<double> tmp_2259;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2260;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2260 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2259 += tmp_2260;
   }
   result += tmp_2259;
   std::complex<double> tmp_2261;
   std::complex<double> tmp_2262;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2263;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2263 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2262 += tmp_2263;
   }
   tmp_2261 += tmp_2262;
   result += (3) * tmp_2261;
   std::complex<double> tmp_2264;
   std::complex<double> tmp_2265;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2265 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2264 += tmp_2265;
   result += (2) * tmp_2264;
   std::complex<double> tmp_2266;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2266 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2266;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2267;
   std::complex<double> tmp_2268;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2268 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2267 += tmp_2268;
   result += (-1) * tmp_2267;
   std::complex<double> tmp_2269;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2270;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2270 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2269 += tmp_2270;
   }
   result += tmp_2269;
   std::complex<double> tmp_2271;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2272;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2272 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*
            CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*
            CpUAhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2271 += tmp_2272;
   }
   result += tmp_2271;
   std::complex<double> tmp_2273;
   std::complex<double> tmp_2274;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2275;
      std::complex<double> tmp_2276;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2276 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2
            ,gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2275 += tmp_2276;
      tmp_2274 += (MCha(gI1)) * tmp_2275;
   }
   tmp_2273 += tmp_2274;
   result += (-2) * tmp_2273;
   std::complex<double> tmp_2277;
   std::complex<double> tmp_2278;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2278 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2277 += tmp_2278;
   result += (-0.5) * tmp_2277;
   std::complex<double> tmp_2279;
   std::complex<double> tmp_2280;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2280 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2279 += tmp_2280;
   result += (-1) * tmp_2279;
   std::complex<double> tmp_2281;
   std::complex<double> tmp_2282;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2282 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2281 += tmp_2282;
   result += (-0.5) * tmp_2281;
   std::complex<double> tmp_2283;
   std::complex<double> tmp_2284;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2285;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2285 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_2284 += tmp_2285;
   }
   tmp_2283 += tmp_2284;
   result += (0.5) * tmp_2283;
   std::complex<double> tmp_2286;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2287;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2287 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2286 += tmp_2287;
   }
   result += tmp_2286;
   std::complex<double> tmp_2288;
   std::complex<double> tmp_2289;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2290;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2290 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_2289 += tmp_2290;
   }
   tmp_2288 += tmp_2289;
   result += (0.5) * tmp_2288;
   std::complex<double> tmp_2291;
   std::complex<double> tmp_2292;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2293;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2293 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2292 += tmp_2293;
   }
   tmp_2291 += tmp_2292;
   result += (3) * tmp_2291;
   std::complex<double> tmp_2294;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2295;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2295 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2294 += tmp_2295;
   }
   result += tmp_2294;
   std::complex<double> tmp_2296;
   std::complex<double> tmp_2297;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2298;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2298 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2297 += tmp_2298;
   }
   tmp_2296 += tmp_2297;
   result += (3) * tmp_2296;
   std::complex<double> tmp_2299;
   std::complex<double> tmp_2300;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2301;
      std::complex<double> tmp_2302;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2302 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2301 += tmp_2302;
      tmp_2300 += (MFd(gI1)) * tmp_2301;
   }
   tmp_2299 += tmp_2300;
   result += (-6) * tmp_2299;
   std::complex<double> tmp_2303;
   std::complex<double> tmp_2304;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2305;
      std::complex<double> tmp_2306;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2306 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2305 += tmp_2306;
      tmp_2304 += (MFe(gI1)) * tmp_2305;
   }
   tmp_2303 += tmp_2304;
   result += (-2) * tmp_2303;
   std::complex<double> tmp_2307;
   std::complex<double> tmp_2308;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2309;
      std::complex<double> tmp_2310;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2310 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2309 += tmp_2310;
      tmp_2308 += (MFu(gI1)) * tmp_2309;
   }
   tmp_2307 += tmp_2308;
   result += (-6) * tmp_2307;
   std::complex<double> tmp_2311;
   std::complex<double> tmp_2312;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2313;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2313 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2312 += tmp_2313;
   }
   tmp_2311 += tmp_2312;
   result += (0.5) * tmp_2311;
   std::complex<double> tmp_2314;
   std::complex<double> tmp_2315;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2316;
      std::complex<double> tmp_2317;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2317 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,
            gI1,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2)
            )*CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2316 += tmp_2317;
      tmp_2315 += (MChi(gI1)) * tmp_2316;
   }
   tmp_2314 += tmp_2315;
   result += (-1) * tmp_2314;
   std::complex<double> tmp_2318;
   std::complex<double> tmp_2319;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2319 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2318 += tmp_2319;
   result += (-3) * tmp_2318;
   std::complex<double> tmp_2320;
   std::complex<double> tmp_2321;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2321 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2320 += tmp_2321;
   result += (-1) * tmp_2320;
   std::complex<double> tmp_2322;
   std::complex<double> tmp_2323;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2323 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2322 += tmp_2323;
   result += (-3) * tmp_2322;
   std::complex<double> tmp_2324;
   std::complex<double> tmp_2325;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2326;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2326 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2325 += tmp_2326;
   }
   tmp_2324 += tmp_2325;
   result += (3) * tmp_2324;
   std::complex<double> tmp_2327;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2328;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2328 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2327 += tmp_2328;
   }
   result += tmp_2327;
   std::complex<double> tmp_2329;
   std::complex<double> tmp_2330;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2331;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2331 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2330 += tmp_2331;
   }
   tmp_2329 += tmp_2330;
   result += (3) * tmp_2329;
   std::complex<double> tmp_2332;
   std::complex<double> tmp_2333;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2333 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2332 += tmp_2333;
   result += (2) * tmp_2332;
   std::complex<double> tmp_2334;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2334 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2334;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*B0(p,0,MVWm)*Conj(CpconjUHpmVWmVP(gO2))*CpconjUHpmVWmVP(gO1);
   result += 4*B0(p,MVWm,MVZ)*Conj(CpconjUHpmVZVWm(gO2))*CpconjUHpmVZVWm(gO1);
   result += 4*A0(MVWm)*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += -(B0(p,MVZ,MVWm)*CpconjUHpmbargWmCgZ(gO1)*CpUHpmgWmCbargZ(gO2));
   result += -(B0(p,MVWm,MVZ)*CpconjUHpmbargZgWm(gO1)*CpUHpmgZbargWm(gO2));
   std::complex<double> tmp_2335;
   std::complex<double> tmp_2336;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2336 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2335 += tmp_2336;
   result += (-1) * tmp_2335;
   std::complex<double> tmp_2337;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2338;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2338 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2337 += tmp_2338;
   }
   result += tmp_2337;
   std::complex<double> tmp_2339;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2340;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2340 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2339 += tmp_2340;
   }
   result += tmp_2339;
   std::complex<double> tmp_2341;
   std::complex<double> tmp_2342;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2342 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2341 += tmp_2342;
   result += (-0.5) * tmp_2341;
   std::complex<double> tmp_2343;
   std::complex<double> tmp_2344;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2344 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2343 += tmp_2344;
   result += (-1) * tmp_2343;
   std::complex<double> tmp_2345;
   std::complex<double> tmp_2346;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2346 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2345 += tmp_2346;
   result += (-0.5) * tmp_2345;
   std::complex<double> tmp_2347;
   std::complex<double> tmp_2348;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2349;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2349 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2348 += tmp_2349;
   }
   tmp_2347 += tmp_2348;
   result += (3) * tmp_2347;
   std::complex<double> tmp_2350;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2351;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2351 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2350 += tmp_2351;
   }
   result += tmp_2350;
   std::complex<double> tmp_2352;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2353;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2353 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2352 += tmp_2353;
   }
   result += tmp_2352;
   std::complex<double> tmp_2354;
   std::complex<double> tmp_2355;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2356;
      std::complex<double> tmp_2357;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2357 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2356 += tmp_2357;
      tmp_2355 += (MFu(gI1)) * tmp_2356;
   }
   tmp_2354 += tmp_2355;
   result += (-6) * tmp_2354;
   std::complex<double> tmp_2358;
   std::complex<double> tmp_2359;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2360;
      std::complex<double> tmp_2361;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2361 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2360 += tmp_2361;
      tmp_2359 += (MFv(gI1)) * tmp_2360;
   }
   tmp_2358 += tmp_2359;
   result += (-2) * tmp_2358;
   std::complex<double> tmp_2362;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2363;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2363 += (Conj(CpconjUHpmChiChaPL(gO2,gI1,gI2))*
            CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHpmChiChaPR(gO2,gI1,gI2)
            )*CpconjUHpmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_2362 += tmp_2363;
   }
   result += tmp_2362;
   std::complex<double> tmp_2364;
   std::complex<double> tmp_2365;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2366;
      std::complex<double> tmp_2367;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2367 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHpmChiChaPR(
            gO2,gI1,gI2))*CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmChiChaPL(gO2,gI1,gI2))*CpconjUHpmChiChaPR(gO1,gI1,gI2))*MCha
            (gI2);
      }
      tmp_2366 += tmp_2367;
      tmp_2365 += (MChi(gI1)) * tmp_2366;
   }
   tmp_2364 += tmp_2365;
   result += (-2) * tmp_2364;
   std::complex<double> tmp_2368;
   std::complex<double> tmp_2369;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2369 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2368 += tmp_2369;
   result += (-3) * tmp_2368;
   std::complex<double> tmp_2370;
   std::complex<double> tmp_2371;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2371 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2370 += tmp_2371;
   result += (-1) * tmp_2370;
   std::complex<double> tmp_2372;
   std::complex<double> tmp_2373;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2373 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2372 += tmp_2373;
   result += (-3) * tmp_2372;
   std::complex<double> tmp_2374;
   std::complex<double> tmp_2375;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2376;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2376 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2375 += tmp_2376;
   }
   tmp_2374 += tmp_2375;
   result += (3) * tmp_2374;
   std::complex<double> tmp_2377;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2377 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2377;
   std::complex<double> tmp_2378;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2378 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2378;
   std::complex<double> tmp_2379;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2379 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2379;
   std::complex<double> tmp_2380;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2380 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2380;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VG(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpVGbargGgG())*B00(p,MVG,MVG);
   result += -1.5*AbsSqr(CpVGVGVG())*(10*B00(p,0,0) + 4*B0(p,0,0)*Sqr(p));
   result += 0;
   std::complex<double> tmp_2381;
   std::complex<double> tmp_2382;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2383;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2383 += (AbsSqr(CpVGbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVGbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2383 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVGbarFdFdPL(gI1,gI2))*CpVGbarFdFdPR(gI1,gI2));
      }
      tmp_2382 += tmp_2383;
   }
   tmp_2381 += tmp_2382;
   result += (0.5) * tmp_2381;
   std::complex<double> tmp_2384;
   std::complex<double> tmp_2385;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2386;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2386 += (AbsSqr(CpVGbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVGbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2386 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVGbarFuFuPL(gI1,gI2))*CpVGbarFuFuPR(gI1,gI2));
      }
      tmp_2385 += tmp_2386;
   }
   tmp_2384 += tmp_2385;
   result += (0.5) * tmp_2384;
   std::complex<double> tmp_2387;
   std::complex<double> tmp_2388;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2388 += A0(MSd(gI1))*CpVGVGconjSdSd(gI1,gI1);
   }
   tmp_2387 += tmp_2388;
   result += (999) * tmp_2387;
   std::complex<double> tmp_2389;
   std::complex<double> tmp_2390;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2390 += A0(MSu(gI1))*CpVGVGconjSuSu(gI1,gI1);
   }
   tmp_2389 += tmp_2390;
   result += (999) * tmp_2389;
   std::complex<double> tmp_2391;
   std::complex<double> tmp_2392;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2393;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2393 += AbsSqr(CpVGconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2392 += tmp_2393;
   }
   tmp_2391 += tmp_2392;
   result += (-2) * tmp_2391;
   std::complex<double> tmp_2394;
   std::complex<double> tmp_2395;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2396;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2396 += AbsSqr(CpVGconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2395 += tmp_2396;
   }
   tmp_2394 += tmp_2395;
   result += (-2) * tmp_2394;
   result += 1.5*((AbsSqr(CpVGGluGluPL(0,0)) + AbsSqr(CpVGGluGluPR(0,0)))*H0(p,
      MGlu,MGlu) + 4*B0(p,MGlu,MGlu)*Re(Conj(CpVGGluGluPL(0,0))*CpVGGluGluPR(0,0))
      *Sqr(MGlu));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VP(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVPbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVPbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVPVPconjVWmVWm1() + CpVPVPconjVWmVWm2() +
      CpVPVPconjVWmVWm3()));
   std::complex<double> tmp_2397;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2397 += A0(MHpm(gI1))*CpVPVPconjHpmHpm(gI1,gI1);
   }
   result += tmp_2397;
   std::complex<double> tmp_2398;
   std::complex<double> tmp_2399;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2400;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2400 += AbsSqr(CpVPconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2399 += tmp_2400;
   }
   tmp_2398 += tmp_2399;
   result += (-4) * tmp_2398;
   std::complex<double> tmp_2401;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2402;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2402 += (AbsSqr(CpVPbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVPbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2402 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVPbarChaChaPL(gI1,gI2))*CpVPbarChaChaPR(gI1,gI2));
      }
      tmp_2401 += tmp_2402;
   }
   result += tmp_2401;
   std::complex<double> tmp_2403;
   std::complex<double> tmp_2404;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2405;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2405 += (AbsSqr(CpVPbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVPbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2405 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVPbarFdFdPL(gI1,gI2))*CpVPbarFdFdPR(gI1,gI2));
      }
      tmp_2404 += tmp_2405;
   }
   tmp_2403 += tmp_2404;
   result += (3) * tmp_2403;
   std::complex<double> tmp_2406;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2407;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2407 += (AbsSqr(CpVPbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVPbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2407 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVPbarFeFePL(gI1,gI2))*CpVPbarFeFePR(gI1,gI2));
      }
      tmp_2406 += tmp_2407;
   }
   result += tmp_2406;
   std::complex<double> tmp_2408;
   std::complex<double> tmp_2409;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2410;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2410 += (AbsSqr(CpVPbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVPbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2410 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVPbarFuFuPL(gI1,gI2))*CpVPbarFuFuPR(gI1,gI2));
      }
      tmp_2409 += tmp_2410;
   }
   tmp_2408 += tmp_2409;
   result += (3) * tmp_2408;
   std::complex<double> tmp_2411;
   std::complex<double> tmp_2412;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2412 += A0(MSd(gI1))*CpVPVPconjSdSd(gI1,gI1);
   }
   tmp_2411 += tmp_2412;
   result += (3) * tmp_2411;
   std::complex<double> tmp_2413;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2413 += A0(MSe(gI1))*CpVPVPconjSeSe(gI1,gI1);
   }
   result += tmp_2413;
   std::complex<double> tmp_2414;
   std::complex<double> tmp_2415;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2415 += A0(MSu(gI1))*CpVPVPconjSuSu(gI1,gI1);
   }
   tmp_2414 += tmp_2415;
   result += (3) * tmp_2414;
   std::complex<double> tmp_2416;
   std::complex<double> tmp_2417;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2418;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2418 += AbsSqr(CpVPconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2417 += tmp_2418;
   }
   tmp_2416 += tmp_2417;
   result += (-12) * tmp_2416;
   std::complex<double> tmp_2419;
   std::complex<double> tmp_2420;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2421;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2421 += AbsSqr(CpVPconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2420 += tmp_2421;
   }
   tmp_2419 += tmp_2420;
   result += (-4) * tmp_2419;
   std::complex<double> tmp_2422;
   std::complex<double> tmp_2423;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2424;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2424 += AbsSqr(CpVPconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2423 += tmp_2424;
   }
   tmp_2422 += tmp_2423;
   result += (-12) * tmp_2422;
   std::complex<double> tmp_2425;
   std::complex<double> tmp_2426;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2426 += AbsSqr(CpVPconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2425 += tmp_2426;
   result += (2) * tmp_2425;
   result += -(AbsSqr(CpVPconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0(
      p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2427;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2427 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2427;
   std::complex<double> tmp_2428;
   std::complex<double> tmp_2429;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2430;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2430 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2429 += tmp_2430;
   }
   tmp_2428 += tmp_2429;
   result += (-4) * tmp_2428;
   std::complex<double> tmp_2431;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2432;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2432 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2432 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2431 += tmp_2432;
   }
   result += tmp_2431;
   std::complex<double> tmp_2433;
   std::complex<double> tmp_2434;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2434 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2433 += tmp_2434;
   result += (0.5) * tmp_2433;
   std::complex<double> tmp_2435;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2435 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2435;
   std::complex<double> tmp_2436;
   std::complex<double> tmp_2437;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2437 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2436 += tmp_2437;
   result += (0.5) * tmp_2436;
   std::complex<double> tmp_2438;
   std::complex<double> tmp_2439;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2440;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2440 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2439 += tmp_2440;
   }
   tmp_2438 += tmp_2439;
   result += (-4) * tmp_2438;
   std::complex<double> tmp_2441;
   std::complex<double> tmp_2442;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2443;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2443 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2442 += tmp_2443;
   }
   tmp_2441 += tmp_2442;
   result += (-4) * tmp_2441;
   std::complex<double> tmp_2444;
   std::complex<double> tmp_2445;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2446;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2446 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2446 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2445 += tmp_2446;
   }
   tmp_2444 += tmp_2445;
   result += (3) * tmp_2444;
   std::complex<double> tmp_2447;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2448;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2448 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2448 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2447 += tmp_2448;
   }
   result += tmp_2447;
   std::complex<double> tmp_2449;
   std::complex<double> tmp_2450;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2451;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2451 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2451 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2450 += tmp_2451;
   }
   tmp_2449 += tmp_2450;
   result += (3) * tmp_2449;
   std::complex<double> tmp_2452;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2453;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2453 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2453 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_2452 += tmp_2453;
   }
   result += tmp_2452;
   std::complex<double> tmp_2454;
   std::complex<double> tmp_2455;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2456;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2456 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2456 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2455 += tmp_2456;
   }
   tmp_2454 += tmp_2455;
   result += (0.5) * tmp_2454;
   std::complex<double> tmp_2457;
   std::complex<double> tmp_2458;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2458 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2457 += tmp_2458;
   result += (3) * tmp_2457;
   std::complex<double> tmp_2459;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2459 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2459;
   std::complex<double> tmp_2460;
   std::complex<double> tmp_2461;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2461 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2460 += tmp_2461;
   result += (3) * tmp_2460;
   std::complex<double> tmp_2462;
   std::complex<double> tmp_2463;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2464;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2464 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2463 += tmp_2464;
   }
   tmp_2462 += tmp_2463;
   result += (-12) * tmp_2462;
   std::complex<double> tmp_2465;
   std::complex<double> tmp_2466;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2467;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2467 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2466 += tmp_2467;
   }
   tmp_2465 += tmp_2466;
   result += (-4) * tmp_2465;
   std::complex<double> tmp_2468;
   std::complex<double> tmp_2469;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2470;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2470 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2469 += tmp_2470;
   }
   tmp_2468 += tmp_2469;
   result += (-12) * tmp_2468;
   std::complex<double> tmp_2471;
   std::complex<double> tmp_2472;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2472 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2471 += tmp_2472;
   result += (2) * tmp_2471;
   std::complex<double> tmp_2473;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2473 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2473;
   result += -(AbsSqr(CpVZconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0(
      p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWmbargPgWm())*B00(p,MVWm,MVP);
   result += AbsSqr(CpconjVWmbargWmCgP())*B00(p,MVP,MVWm);
   result += AbsSqr(CpconjVWmbargWmCgZ())*B00(p,MVZ,MVWm);
   result += AbsSqr(CpconjVWmbargZgWm())*B00(p,MVWm,MVZ);
   result += -(A0(MVWm)*(4*CpVWmconjVWmconjVWmVWm1() + CpVWmconjVWmconjVWmVWm2(
      ) + CpVWmconjVWmconjVWmVWm3()));
   result += 0;
   result += -0.5*A0(MVZ)*(4*CpVWmconjVWmVZVZ1() + CpVWmconjVWmVZVZ2() +
      CpVWmconjVWmVZVZ3());
   std::complex<double> tmp_2474;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2474 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_2474;
   std::complex<double> tmp_2475;
   std::complex<double> tmp_2476;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2477;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2477 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_2476 += tmp_2477;
   }
   tmp_2475 += tmp_2476;
   result += (-4) * tmp_2475;
   std::complex<double> tmp_2478;
   std::complex<double> tmp_2479;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2480;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2480 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_2479 += tmp_2480;
   }
   tmp_2478 += tmp_2479;
   result += (-4) * tmp_2478;
   std::complex<double> tmp_2481;
   std::complex<double> tmp_2482;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2482 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_2481 += tmp_2482;
   result += (0.5) * tmp_2481;
   std::complex<double> tmp_2483;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2483 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2483;
   std::complex<double> tmp_2484;
   std::complex<double> tmp_2485;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2485 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2484 += tmp_2485;
   result += (0.5) * tmp_2484;
   std::complex<double> tmp_2486;
   std::complex<double> tmp_2487;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2488;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2488 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_2488 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_2487 += tmp_2488;
   }
   tmp_2486 += tmp_2487;
   result += (3) * tmp_2486;
   std::complex<double> tmp_2489;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2490;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2490 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_2490 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_2489 += tmp_2490;
   }
   result += tmp_2489;
   std::complex<double> tmp_2491;
   std::complex<double> tmp_2492;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2493;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2493 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2492 += tmp_2493;
   }
   tmp_2491 += tmp_2492;
   result += (-4) * tmp_2491;
   std::complex<double> tmp_2494;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2495;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2495 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_2495 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_2494 += tmp_2495;
   }
   result += tmp_2494;
   std::complex<double> tmp_2496;
   std::complex<double> tmp_2497;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2497 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2496 += tmp_2497;
   result += (3) * tmp_2496;
   std::complex<double> tmp_2498;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2498 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2498;
   std::complex<double> tmp_2499;
   std::complex<double> tmp_2500;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2500 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2499 += tmp_2500;
   result += (3) * tmp_2499;
   std::complex<double> tmp_2501;
   std::complex<double> tmp_2502;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2503;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2503 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2502 += tmp_2503;
   }
   tmp_2501 += tmp_2502;
   result += (-12) * tmp_2501;
   std::complex<double> tmp_2504;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2504 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_2504;
   std::complex<double> tmp_2505;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2505 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_2505;
   std::complex<double> tmp_2506;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2506 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2506;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2507;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2508;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2508 += B0(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2507 += tmp_2508;
   }
   result += tmp_2507;
   std::complex<double> tmp_2509;
   std::complex<double> tmp_2510;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2510 += B0(p,MCha(gI1),MVWm)*Conj(CpUChibarChaVWmPR(gO2,gI1))*
         CpUChibarChaVWmPL(gO1,gI1)*MCha(gI1);
   }
   tmp_2509 += tmp_2510;
   result += (-4) * tmp_2509;
   std::complex<double> tmp_2511;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2512;
      std::complex<double> tmp_2513;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2513 += B0(p,MCha(gI1),MHpm(gI2))*Conj(CpUChibarChaHpmPL(gO2
            ,gI1,gI2))*CpUChibarChaHpmPR(gO1,gI1,gI2);
      }
      tmp_2512 += tmp_2513;
      tmp_2511 += (MCha(gI1)) * tmp_2512;
   }
   result += tmp_2511;
   std::complex<double> tmp_2514;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2515;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2515 += B0(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2514 += tmp_2515;
   }
   result += tmp_2514;
   std::complex<double> tmp_2516;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2517;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2517 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2516 += tmp_2517;
   }
   result += tmp_2516;
   std::complex<double> tmp_2518;
   std::complex<double> tmp_2519;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2520;
      std::complex<double> tmp_2521;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2521 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpUChibarFdSdPL(gO2,gI1
            ,gI2))*CpUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2520 += tmp_2521;
      tmp_2519 += (MFd(gI1)) * tmp_2520;
   }
   tmp_2518 += tmp_2519;
   result += (3) * tmp_2518;
   std::complex<double> tmp_2522;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2523;
      std::complex<double> tmp_2524;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2524 += B0(p,MFe(gI1),MSe(gI2))*Conj(CpUChibarFeSePL(gO2,gI1
            ,gI2))*CpUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2523 += tmp_2524;
      tmp_2522 += (MFe(gI1)) * tmp_2523;
   }
   result += tmp_2522;
   std::complex<double> tmp_2525;
   std::complex<double> tmp_2526;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2527;
      std::complex<double> tmp_2528;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2528 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpUChibarFuSuPL(gO2,gI1
            ,gI2))*CpUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2527 += tmp_2528;
      tmp_2526 += (MFu(gI1)) * tmp_2527;
   }
   tmp_2525 += tmp_2526;
   result += (3) * tmp_2525;
   std::complex<double> tmp_2529;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2530;
      std::complex<double> tmp_2531;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2531 += B0(p,MFv(gI1),MSv(gI2))*Conj(CpUChibarFvSvPL(gO2,gI1
            ,gI2))*CpUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2530 += tmp_2531;
      tmp_2529 += (MFv(gI1)) * tmp_2530;
   }
   result += tmp_2529;
   std::complex<double> tmp_2532;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2533;
      std::complex<double> tmp_2534;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2534 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2533 += tmp_2534;
      tmp_2532 += (MChi(gI1)) * tmp_2533;
   }
   result += tmp_2532;
   std::complex<double> tmp_2535;
   std::complex<double> tmp_2536;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2537;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2537 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2536 += tmp_2537;
   }
   tmp_2535 += tmp_2536;
   result += (3) * tmp_2535;
   std::complex<double> tmp_2538;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2539;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2539 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2538 += tmp_2539;
   }
   result += tmp_2538;
   std::complex<double> tmp_2540;
   std::complex<double> tmp_2541;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2542;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2542 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2541 += tmp_2542;
   }
   tmp_2540 += tmp_2541;
   result += (3) * tmp_2540;
   std::complex<double> tmp_2543;
   std::complex<double> tmp_2544;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2544 += B0(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2543 += tmp_2544;
   result += (-4) * tmp_2543;
   std::complex<double> tmp_2545;
   std::complex<double> tmp_2546;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2546 += B0(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2545 += tmp_2546;
   result += (-4) * tmp_2545;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2547;
   std::complex<double> tmp_2548;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2548 += B1(p,MCha(gI1),MVWm)*Conj(CpUChibarChaVWmPL(gO2,gI1))*
         CpUChibarChaVWmPL(gO1,gI1);
   }
   tmp_2547 += tmp_2548;
   result += (-1) * tmp_2547;
   std::complex<double> tmp_2549;
   std::complex<double> tmp_2550;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2551;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2551 += B1(p,MCha(gI1),MHpm(gI2))*Conj(CpUChibarChaHpmPR(gO2
            ,gI1,gI2))*CpUChibarChaHpmPR(gO1,gI1,gI2);
      }
      tmp_2550 += tmp_2551;
   }
   tmp_2549 += tmp_2550;
   result += (-0.5) * tmp_2549;
   std::complex<double> tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2554;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2554 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPR(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2);
      }
      tmp_2553 += tmp_2554;
   }
   tmp_2552 += tmp_2553;
   result += (-0.5) * tmp_2552;
   std::complex<double> tmp_2555;
   std::complex<double> tmp_2556;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2557;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2557 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpUChibarFvSvPR(gO2,gI1
            ,gI2))*CpUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2556 += tmp_2557;
   }
   tmp_2555 += tmp_2556;
   result += (-0.5) * tmp_2555;
   std::complex<double> tmp_2558;
   std::complex<double> tmp_2559;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2560;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2560 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPR(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2);
      }
      tmp_2559 += tmp_2560;
   }
   tmp_2558 += tmp_2559;
   result += (-0.5) * tmp_2558;
   std::complex<double> tmp_2561;
   std::complex<double> tmp_2562;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2563;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2563 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_2562 += tmp_2563;
   }
   tmp_2561 += tmp_2562;
   result += (-0.5) * tmp_2561;
   std::complex<double> tmp_2564;
   std::complex<double> tmp_2565;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2566;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2566 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpUChibarFdSdPR(gO2,gI1
            ,gI2))*CpUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2565 += tmp_2566;
   }
   tmp_2564 += tmp_2565;
   result += (-1.5) * tmp_2564;
   std::complex<double> tmp_2567;
   std::complex<double> tmp_2568;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2569;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2569 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpUChibarFeSePR(gO2,gI1
            ,gI2))*CpUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2568 += tmp_2569;
   }
   tmp_2567 += tmp_2568;
   result += (-0.5) * tmp_2567;
   std::complex<double> tmp_2570;
   std::complex<double> tmp_2571;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2572;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2572 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpUChibarFuSuPR(gO2,gI1
            ,gI2))*CpUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2571 += tmp_2572;
   }
   tmp_2570 += tmp_2571;
   result += (-1.5) * tmp_2570;
   std::complex<double> tmp_2573;
   std::complex<double> tmp_2574;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2575;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2575 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2574 += tmp_2575;
   }
   tmp_2573 += tmp_2574;
   result += (-0.5) * tmp_2573;
   std::complex<double> tmp_2576;
   std::complex<double> tmp_2577;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2578;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2578 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPR(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_2577 += tmp_2578;
   }
   tmp_2576 += tmp_2577;
   result += (-1.5) * tmp_2576;
   std::complex<double> tmp_2579;
   std::complex<double> tmp_2580;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2581;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2581 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePR(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_2580 += tmp_2581;
   }
   tmp_2579 += tmp_2580;
   result += (-0.5) * tmp_2579;
   std::complex<double> tmp_2582;
   std::complex<double> tmp_2583;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2584;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2584 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPR(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_2583 += tmp_2584;
   }
   tmp_2582 += tmp_2583;
   result += (-1.5) * tmp_2582;
   std::complex<double> tmp_2585;
   std::complex<double> tmp_2586;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2586 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPL(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_2585 += tmp_2586;
   result += (-1) * tmp_2585;
   std::complex<double> tmp_2587;
   std::complex<double> tmp_2588;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2588 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_2587 += tmp_2588;
   result += (-1) * tmp_2587;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2589;
   std::complex<double> tmp_2590;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2590 += B1(p,MCha(gI1),MVWm)*Conj(CpUChibarChaVWmPR(gO2,gI1))*
         CpUChibarChaVWmPR(gO1,gI1);
   }
   tmp_2589 += tmp_2590;
   result += (-1) * tmp_2589;
   std::complex<double> tmp_2591;
   std::complex<double> tmp_2592;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2593;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2593 += B1(p,MCha(gI1),MHpm(gI2))*Conj(CpUChibarChaHpmPL(gO2
            ,gI1,gI2))*CpUChibarChaHpmPL(gO1,gI1,gI2);
      }
      tmp_2592 += tmp_2593;
   }
   tmp_2591 += tmp_2592;
   result += (-0.5) * tmp_2591;
   std::complex<double> tmp_2594;
   std::complex<double> tmp_2595;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2596;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2596 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPL(gO1,gI1,gI2);
      }
      tmp_2595 += tmp_2596;
   }
   tmp_2594 += tmp_2595;
   result += (-0.5) * tmp_2594;
   std::complex<double> tmp_2597;
   std::complex<double> tmp_2598;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2599;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2599 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpUChibarFvSvPL(gO2,gI1
            ,gI2))*CpUChibarFvSvPL(gO1,gI1,gI2);
      }
      tmp_2598 += tmp_2599;
   }
   tmp_2597 += tmp_2598;
   result += (-0.5) * tmp_2597;
   std::complex<double> tmp_2600;
   std::complex<double> tmp_2601;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2602;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2602 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPL(gO1,gI1,gI2);
      }
      tmp_2601 += tmp_2602;
   }
   tmp_2600 += tmp_2601;
   result += (-0.5) * tmp_2600;
   std::complex<double> tmp_2603;
   std::complex<double> tmp_2604;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2605;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2605 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_2604 += tmp_2605;
   }
   tmp_2603 += tmp_2604;
   result += (-0.5) * tmp_2603;
   std::complex<double> tmp_2606;
   std::complex<double> tmp_2607;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2608;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2608 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpUChibarFdSdPL(gO2,gI1
            ,gI2))*CpUChibarFdSdPL(gO1,gI1,gI2);
      }
      tmp_2607 += tmp_2608;
   }
   tmp_2606 += tmp_2607;
   result += (-1.5) * tmp_2606;
   std::complex<double> tmp_2609;
   std::complex<double> tmp_2610;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2611;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2611 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpUChibarFeSePL(gO2,gI1
            ,gI2))*CpUChibarFeSePL(gO1,gI1,gI2);
      }
      tmp_2610 += tmp_2611;
   }
   tmp_2609 += tmp_2610;
   result += (-0.5) * tmp_2609;
   std::complex<double> tmp_2612;
   std::complex<double> tmp_2613;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2614;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2614 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpUChibarFuSuPL(gO2,gI1
            ,gI2))*CpUChibarFuSuPL(gO1,gI1,gI2);
      }
      tmp_2613 += tmp_2614;
   }
   tmp_2612 += tmp_2613;
   result += (-1.5) * tmp_2612;
   std::complex<double> tmp_2615;
   std::complex<double> tmp_2616;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2617;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2617 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_2616 += tmp_2617;
   }
   tmp_2615 += tmp_2616;
   result += (-0.5) * tmp_2615;
   std::complex<double> tmp_2618;
   std::complex<double> tmp_2619;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2620;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2620 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_2619 += tmp_2620;
   }
   tmp_2618 += tmp_2619;
   result += (-1.5) * tmp_2618;
   std::complex<double> tmp_2621;
   std::complex<double> tmp_2622;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2623;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2623 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_2622 += tmp_2623;
   }
   tmp_2621 += tmp_2622;
   result += (-0.5) * tmp_2621;
   std::complex<double> tmp_2624;
   std::complex<double> tmp_2625;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2626;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2626 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_2625 += tmp_2626;
   }
   tmp_2624 += tmp_2625;
   result += (-1.5) * tmp_2624;
   std::complex<double> tmp_2627;
   std::complex<double> tmp_2628;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2628 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_2627 += tmp_2628;
   result += (-1) * tmp_2627;
   std::complex<double> tmp_2629;
   std::complex<double> tmp_2630;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2630 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_2629 += tmp_2630;
   result += (-1) * tmp_2629;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2631;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2632;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2632 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2631 += tmp_2632;
   }
   result += tmp_2631;
   std::complex<double> tmp_2633;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2634;
      std::complex<double> tmp_2635;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2635 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2634 += tmp_2635;
      tmp_2633 += (MCha(gI1)) * tmp_2634;
   }
   result += tmp_2633;
   std::complex<double> tmp_2636;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2637;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2637 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2636 += tmp_2637;
   }
   result += tmp_2636;
   std::complex<double> tmp_2638;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2639;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2639 += B0(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2638 += tmp_2639;
   }
   result += tmp_2638;
   std::complex<double> tmp_2640;
   std::complex<double> tmp_2641;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2642;
      std::complex<double> tmp_2643;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2643 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2642 += tmp_2643;
      tmp_2641 += (MFu(gI1)) * tmp_2642;
   }
   tmp_2640 += tmp_2641;
   result += (3) * tmp_2640;
   std::complex<double> tmp_2644;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2645;
      std::complex<double> tmp_2646;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2646 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_2645 += tmp_2646;
      tmp_2644 += (MFv(gI1)) * tmp_2645;
   }
   result += tmp_2644;
   std::complex<double> tmp_2647;
   std::complex<double> tmp_2648;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2649;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2649 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2648 += tmp_2649;
   }
   tmp_2647 += tmp_2648;
   result += (3) * tmp_2647;
   std::complex<double> tmp_2650;
   std::complex<double> tmp_2651;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2651 += B0(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2650 += tmp_2651;
   result += (-4) * tmp_2650;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2653 += B0(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2652 += tmp_2653;
   result += (-4) * tmp_2652;
   std::complex<double> tmp_2654;
   std::complex<double> tmp_2655;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2655 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2654 += tmp_2655;
   result += (-4) * tmp_2654;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2656;
   std::complex<double> tmp_2657;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2658;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2658 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2657 += tmp_2658;
   }
   tmp_2656 += tmp_2657;
   result += (-0.5) * tmp_2656;
   std::complex<double> tmp_2659;
   std::complex<double> tmp_2660;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2661;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2661 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPR(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2);
      }
      tmp_2660 += tmp_2661;
   }
   tmp_2659 += tmp_2660;
   result += (-0.5) * tmp_2659;
   std::complex<double> tmp_2662;
   std::complex<double> tmp_2663;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2664;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2664 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_2663 += tmp_2664;
   }
   tmp_2662 += tmp_2663;
   result += (-0.5) * tmp_2662;
   std::complex<double> tmp_2665;
   std::complex<double> tmp_2666;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2667;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2667 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePR(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2);
      }
      tmp_2666 += tmp_2667;
   }
   tmp_2665 += tmp_2666;
   result += (-0.5) * tmp_2665;
   std::complex<double> tmp_2668;
   std::complex<double> tmp_2669;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2670;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2670 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPR(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2669 += tmp_2670;
   }
   tmp_2668 += tmp_2669;
   result += (-1.5) * tmp_2668;
   std::complex<double> tmp_2671;
   std::complex<double> tmp_2672;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2673;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2673 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePR(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_2672 += tmp_2673;
   }
   tmp_2671 += tmp_2672;
   result += (-0.5) * tmp_2671;
   std::complex<double> tmp_2674;
   std::complex<double> tmp_2675;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2676;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2676 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPR(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2);
      }
      tmp_2675 += tmp_2676;
   }
   tmp_2674 += tmp_2675;
   result += (-1.5) * tmp_2674;
   std::complex<double> tmp_2677;
   std::complex<double> tmp_2678;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2678 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_2677 += tmp_2678;
   result += (-1) * tmp_2677;
   std::complex<double> tmp_2679;
   std::complex<double> tmp_2680;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2680 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPL(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_2679 += tmp_2680;
   result += (-1) * tmp_2679;
   std::complex<double> tmp_2681;
   std::complex<double> tmp_2682;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2682 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPL(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_2681 += tmp_2682;
   result += (-1) * tmp_2681;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2683;
   std::complex<double> tmp_2684;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2685;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2685 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_2684 += tmp_2685;
   }
   tmp_2683 += tmp_2684;
   result += (-0.5) * tmp_2683;
   std::complex<double> tmp_2686;
   std::complex<double> tmp_2687;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2688;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2688 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPL(gO1,gI1,gI2);
      }
      tmp_2687 += tmp_2688;
   }
   tmp_2686 += tmp_2687;
   result += (-0.5) * tmp_2686;
   std::complex<double> tmp_2689;
   std::complex<double> tmp_2690;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2691;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2691 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_2690 += tmp_2691;
   }
   tmp_2689 += tmp_2690;
   result += (-0.5) * tmp_2689;
   std::complex<double> tmp_2692;
   std::complex<double> tmp_2693;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2694;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2694 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePL(gO1,gI1,gI2);
      }
      tmp_2693 += tmp_2694;
   }
   tmp_2692 += tmp_2693;
   result += (-0.5) * tmp_2692;
   std::complex<double> tmp_2695;
   std::complex<double> tmp_2696;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2697;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2697 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2);
      }
      tmp_2696 += tmp_2697;
   }
   tmp_2695 += tmp_2696;
   result += (-1.5) * tmp_2695;
   std::complex<double> tmp_2698;
   std::complex<double> tmp_2699;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2700;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2700 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2);
      }
      tmp_2699 += tmp_2700;
   }
   tmp_2698 += tmp_2699;
   result += (-0.5) * tmp_2698;
   std::complex<double> tmp_2701;
   std::complex<double> tmp_2702;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2703;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2703 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPL(gO1,gI1,gI2);
      }
      tmp_2702 += tmp_2703;
   }
   tmp_2701 += tmp_2702;
   result += (-1.5) * tmp_2701;
   std::complex<double> tmp_2704;
   std::complex<double> tmp_2705;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2705 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_2704 += tmp_2705;
   result += (-1) * tmp_2704;
   std::complex<double> tmp_2706;
   std::complex<double> tmp_2707;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2707 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_2706 += tmp_2707;
   result += (-1) * tmp_2706;
   std::complex<double> tmp_2708;
   std::complex<double> tmp_2709;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2709 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_2708 += tmp_2709;
   result += (-1) * tmp_2708;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2710;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2711;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2711 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2710 += tmp_2711;
   }
   result += tmp_2710;
   std::complex<double> tmp_2712;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2713;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2713 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2712 += tmp_2713;
   }
   result += tmp_2712;
   std::complex<double> tmp_2714;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2715;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2715 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2714 += tmp_2715;
   }
   result += tmp_2714;
   std::complex<double> tmp_2716;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2717;
      std::complex<double> tmp_2718;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2718 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2717 += tmp_2718;
      tmp_2716 += (MFe(gI1)) * tmp_2717;
   }
   result += tmp_2716;
   std::complex<double> tmp_2719;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2720;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2720 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2719 += tmp_2720;
   }
   result += tmp_2719;
   std::complex<double> tmp_2721;
   std::complex<double> tmp_2722;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2722 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2721 += tmp_2722;
   result += (-4) * tmp_2721;
   std::complex<double> tmp_2723;
   std::complex<double> tmp_2724;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2724 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2723 += tmp_2724;
   result += (-4) * tmp_2723;
   std::complex<double> tmp_2725;
   std::complex<double> tmp_2726;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2726 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_2725 += tmp_2726;
   result += (-4) * tmp_2725;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2727;
   std::complex<double> tmp_2728;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2729;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2729 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_2728 += tmp_2729;
   }
   tmp_2727 += tmp_2728;
   result += (-0.5) * tmp_2727;
   std::complex<double> tmp_2730;
   std::complex<double> tmp_2731;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2732;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2732 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPR(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_2731 += tmp_2732;
   }
   tmp_2730 += tmp_2731;
   result += (-0.5) * tmp_2730;
   std::complex<double> tmp_2733;
   std::complex<double> tmp_2734;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2735;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2735 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2734 += tmp_2735;
   }
   tmp_2733 += tmp_2734;
   result += (-0.5) * tmp_2733;
   std::complex<double> tmp_2736;
   std::complex<double> tmp_2737;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2738;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2738 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_2737 += tmp_2738;
   }
   tmp_2736 += tmp_2737;
   result += (-0.5) * tmp_2736;
   std::complex<double> tmp_2739;
   std::complex<double> tmp_2740;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2741;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2741 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_2740 += tmp_2741;
   }
   tmp_2739 += tmp_2740;
   result += (-0.5) * tmp_2739;
   std::complex<double> tmp_2742;
   std::complex<double> tmp_2743;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2743 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_2742 += tmp_2743;
   result += (-1) * tmp_2742;
   std::complex<double> tmp_2744;
   std::complex<double> tmp_2745;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2745 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_2744 += tmp_2745;
   result += (-1) * tmp_2744;
   std::complex<double> tmp_2746;
   std::complex<double> tmp_2747;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2747 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_2746 += tmp_2747;
   result += (-1) * tmp_2746;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2748;
   std::complex<double> tmp_2749;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2750;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2750 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_2749 += tmp_2750;
   }
   tmp_2748 += tmp_2749;
   result += (-0.5) * tmp_2748;
   std::complex<double> tmp_2751;
   std::complex<double> tmp_2752;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2753;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2753 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_2752 += tmp_2753;
   }
   tmp_2751 += tmp_2752;
   result += (-0.5) * tmp_2751;
   std::complex<double> tmp_2754;
   std::complex<double> tmp_2755;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2756;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2756 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_2755 += tmp_2756;
   }
   tmp_2754 += tmp_2755;
   result += (-0.5) * tmp_2754;
   std::complex<double> tmp_2757;
   std::complex<double> tmp_2758;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2759;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2759 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_2758 += tmp_2759;
   }
   tmp_2757 += tmp_2758;
   result += (-0.5) * tmp_2757;
   std::complex<double> tmp_2760;
   std::complex<double> tmp_2761;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2762;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2762 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_2761 += tmp_2762;
   }
   tmp_2760 += tmp_2761;
   result += (-0.5) * tmp_2760;
   std::complex<double> tmp_2763;
   std::complex<double> tmp_2764;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2764 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_2763 += tmp_2764;
   result += (-1) * tmp_2763;
   std::complex<double> tmp_2765;
   std::complex<double> tmp_2766;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2766 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_2765 += tmp_2766;
   result += (-1) * tmp_2765;
   std::complex<double> tmp_2767;
   std::complex<double> tmp_2768;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2768 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_2767 += tmp_2768;
   result += (-1) * tmp_2767;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2769;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2770;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2770 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2769 += tmp_2770;
   }
   result += tmp_2769;
   std::complex<double> tmp_2771;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2772;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2772 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2771 += tmp_2772;
   }
   result += tmp_2771;
   std::complex<double> tmp_2773;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2774;
      std::complex<double> tmp_2775;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2775 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_2774 += tmp_2775;
      tmp_2773 += (MFd(gI1)) * tmp_2774;
   }
   result += tmp_2773;
   std::complex<double> tmp_2776;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2777;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2777 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2776 += tmp_2777;
   }
   result += tmp_2776;
   std::complex<double> tmp_2778;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2779;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2779 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2778 += tmp_2779;
   }
   result += tmp_2778;
   std::complex<double> tmp_2780;
   std::complex<double> tmp_2781;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2781 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2780 += tmp_2781;
   result += (-5.333333333333333) * tmp_2780;
   std::complex<double> tmp_2782;
   std::complex<double> tmp_2783;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2783 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2782 += tmp_2783;
   result += (-4) * tmp_2782;
   std::complex<double> tmp_2784;
   std::complex<double> tmp_2785;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2785 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2784 += tmp_2785;
   result += (-4) * tmp_2784;
   std::complex<double> tmp_2786;
   std::complex<double> tmp_2787;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2787 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2786 += tmp_2787;
   result += (-4) * tmp_2786;
   std::complex<double> tmp_2788;
   std::complex<double> tmp_2789;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2789 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,0))*
         CpbarUFdSdGluPR(gO1,gI1,0);
   }
   tmp_2788 += tmp_2789;
   result += (1.3333333333333333*MGlu) * tmp_2788;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2790;
   std::complex<double> tmp_2791;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2792;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2792 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_2791 += tmp_2792;
   }
   tmp_2790 += tmp_2791;
   result += (-0.5) * tmp_2790;
   std::complex<double> tmp_2793;
   std::complex<double> tmp_2794;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2795;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2795 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_2794 += tmp_2795;
   }
   tmp_2793 += tmp_2794;
   result += (-0.5) * tmp_2793;
   std::complex<double> tmp_2796;
   std::complex<double> tmp_2797;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2798;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2798 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_2797 += tmp_2798;
   }
   tmp_2796 += tmp_2797;
   result += (-0.5) * tmp_2796;
   std::complex<double> tmp_2799;
   std::complex<double> tmp_2800;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2800 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,0))*
         CpbarUFdSdGluPR(gO1,gI1,0);
   }
   tmp_2799 += tmp_2800;
   result += (-0.6666666666666666) * tmp_2799;
   std::complex<double> tmp_2801;
   std::complex<double> tmp_2802;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2803;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2803 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPR(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_2802 += tmp_2803;
   }
   tmp_2801 += tmp_2802;
   result += (-0.5) * tmp_2801;
   std::complex<double> tmp_2804;
   std::complex<double> tmp_2805;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2806;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2806 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_2805 += tmp_2806;
   }
   tmp_2804 += tmp_2805;
   result += (-0.5) * tmp_2804;
   std::complex<double> tmp_2807;
   std::complex<double> tmp_2808;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2808 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_2807 += tmp_2808;
   result += (-1.3333333333333333) * tmp_2807;
   std::complex<double> tmp_2809;
   std::complex<double> tmp_2810;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2810 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_2809 += tmp_2810;
   result += (-1) * tmp_2809;
   std::complex<double> tmp_2811;
   std::complex<double> tmp_2812;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2812 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_2811 += tmp_2812;
   result += (-1) * tmp_2811;
   std::complex<double> tmp_2813;
   std::complex<double> tmp_2814;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2814 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_2813 += tmp_2814;
   result += (-1) * tmp_2813;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2815;
   std::complex<double> tmp_2816;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2817;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2817 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_2816 += tmp_2817;
   }
   tmp_2815 += tmp_2816;
   result += (-0.5) * tmp_2815;
   std::complex<double> tmp_2818;
   std::complex<double> tmp_2819;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2820;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2820 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_2819 += tmp_2820;
   }
   tmp_2818 += tmp_2819;
   result += (-0.5) * tmp_2818;
   std::complex<double> tmp_2821;
   std::complex<double> tmp_2822;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2823;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2823 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_2822 += tmp_2823;
   }
   tmp_2821 += tmp_2822;
   result += (-0.5) * tmp_2821;
   std::complex<double> tmp_2824;
   std::complex<double> tmp_2825;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2825 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,0))*
         CpbarUFdSdGluPL(gO1,gI1,0);
   }
   tmp_2824 += tmp_2825;
   result += (-0.6666666666666666) * tmp_2824;
   std::complex<double> tmp_2826;
   std::complex<double> tmp_2827;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2828;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2828 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_2827 += tmp_2828;
   }
   tmp_2826 += tmp_2827;
   result += (-0.5) * tmp_2826;
   std::complex<double> tmp_2829;
   std::complex<double> tmp_2830;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2831;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2831 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_2830 += tmp_2831;
   }
   tmp_2829 += tmp_2830;
   result += (-0.5) * tmp_2829;
   std::complex<double> tmp_2832;
   std::complex<double> tmp_2833;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2833 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_2832 += tmp_2833;
   result += (-1.3333333333333333) * tmp_2832;
   std::complex<double> tmp_2834;
   std::complex<double> tmp_2835;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2835 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_2834 += tmp_2835;
   result += (-1) * tmp_2834;
   std::complex<double> tmp_2836;
   std::complex<double> tmp_2837;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2837 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_2836 += tmp_2837;
   result += (-1) * tmp_2836;
   std::complex<double> tmp_2838;
   std::complex<double> tmp_2839;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2839 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_2838 += tmp_2839;
   result += (-1) * tmp_2838;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2840;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2841;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2841 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2840 += tmp_2841;
   }
   result += tmp_2840;
   std::complex<double> tmp_2842;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2843;
      std::complex<double> tmp_2844;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2844 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_2843 += tmp_2844;
      tmp_2842 += (MCha(gI1)) * tmp_2843;
   }
   result += tmp_2842;
   std::complex<double> tmp_2845;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2846;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2846 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2845 += tmp_2846;
   }
   result += tmp_2845;
   std::complex<double> tmp_2847;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2848;
      std::complex<double> tmp_2849;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2849 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_2848 += tmp_2849;
      tmp_2847 += (MFu(gI1)) * tmp_2848;
   }
   result += tmp_2847;
   std::complex<double> tmp_2850;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2851;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2851 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2850 += tmp_2851;
   }
   result += tmp_2850;
   std::complex<double> tmp_2852;
   std::complex<double> tmp_2853;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2853 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2852 += tmp_2853;
   result += (-4) * tmp_2852;
   std::complex<double> tmp_2854;
   std::complex<double> tmp_2855;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2855 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2854 += tmp_2855;
   result += (-5.333333333333333) * tmp_2854;
   std::complex<double> tmp_2856;
   std::complex<double> tmp_2857;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2857 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2856 += tmp_2857;
   result += (-4) * tmp_2856;
   std::complex<double> tmp_2858;
   std::complex<double> tmp_2859;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2859 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2858 += tmp_2859;
   result += (-4) * tmp_2858;
   std::complex<double> tmp_2860;
   std::complex<double> tmp_2861;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2861 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_2860 += tmp_2861;
   result += (1.3333333333333333*MGlu) * tmp_2860;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2862;
   std::complex<double> tmp_2863;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2864;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2864 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_2863 += tmp_2864;
   }
   tmp_2862 += tmp_2863;
   result += (-0.5) * tmp_2862;
   std::complex<double> tmp_2865;
   std::complex<double> tmp_2866;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2867;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2867 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_2866 += tmp_2867;
   }
   tmp_2865 += tmp_2866;
   result += (-0.5) * tmp_2865;
   std::complex<double> tmp_2868;
   std::complex<double> tmp_2869;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2870;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2870 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_2869 += tmp_2870;
   }
   tmp_2868 += tmp_2869;
   result += (-0.5) * tmp_2868;
   std::complex<double> tmp_2871;
   std::complex<double> tmp_2872;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2873;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2873 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_2872 += tmp_2873;
   }
   tmp_2871 += tmp_2872;
   result += (-0.5) * tmp_2871;
   std::complex<double> tmp_2874;
   std::complex<double> tmp_2875;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2875 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_2874 += tmp_2875;
   result += (-0.6666666666666666) * tmp_2874;
   std::complex<double> tmp_2876;
   std::complex<double> tmp_2877;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2878;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2878 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_2877 += tmp_2878;
   }
   tmp_2876 += tmp_2877;
   result += (-0.5) * tmp_2876;
   std::complex<double> tmp_2879;
   std::complex<double> tmp_2880;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2880 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_2879 += tmp_2880;
   result += (-1) * tmp_2879;
   std::complex<double> tmp_2881;
   std::complex<double> tmp_2882;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2882 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_2881 += tmp_2882;
   result += (-1.3333333333333333) * tmp_2881;
   std::complex<double> tmp_2883;
   std::complex<double> tmp_2884;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2884 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_2883 += tmp_2884;
   result += (-1) * tmp_2883;
   std::complex<double> tmp_2885;
   std::complex<double> tmp_2886;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2886 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_2885 += tmp_2886;
   result += (-1) * tmp_2885;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2887;
   std::complex<double> tmp_2888;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2889;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2889 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_2888 += tmp_2889;
   }
   tmp_2887 += tmp_2888;
   result += (-0.5) * tmp_2887;
   std::complex<double> tmp_2890;
   std::complex<double> tmp_2891;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2892;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2892 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_2891 += tmp_2892;
   }
   tmp_2890 += tmp_2891;
   result += (-0.5) * tmp_2890;
   std::complex<double> tmp_2893;
   std::complex<double> tmp_2894;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2895;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2895 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_2894 += tmp_2895;
   }
   tmp_2893 += tmp_2894;
   result += (-0.5) * tmp_2893;
   std::complex<double> tmp_2896;
   std::complex<double> tmp_2897;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2898;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2898 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_2897 += tmp_2898;
   }
   tmp_2896 += tmp_2897;
   result += (-0.5) * tmp_2896;
   std::complex<double> tmp_2899;
   std::complex<double> tmp_2900;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2900 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPL(gO1,gI1,0);
   }
   tmp_2899 += tmp_2900;
   result += (-0.6666666666666666) * tmp_2899;
   std::complex<double> tmp_2901;
   std::complex<double> tmp_2902;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2903;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2903 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_2902 += tmp_2903;
   }
   tmp_2901 += tmp_2902;
   result += (-0.5) * tmp_2901;
   std::complex<double> tmp_2904;
   std::complex<double> tmp_2905;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2905 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_2904 += tmp_2905;
   result += (-1) * tmp_2904;
   std::complex<double> tmp_2906;
   std::complex<double> tmp_2907;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2907 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_2906 += tmp_2907;
   result += (-1.3333333333333333) * tmp_2906;
   std::complex<double> tmp_2908;
   std::complex<double> tmp_2909;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2909 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_2908 += tmp_2909;
   result += (-1) * tmp_2908;
   std::complex<double> tmp_2910;
   std::complex<double> tmp_2911;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2911 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_2910 += tmp_2911;
   result += (-1) * tmp_2910;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2912;
   std::complex<double> tmp_2913;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2914;
      std::complex<double> tmp_2915;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2915 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpGlubarFdSdPL(gI1,gI2)
            )*CpGlubarFdSdPR(gI1,gI2);
      }
      tmp_2914 += tmp_2915;
      tmp_2913 += (MFd(gI1)) * tmp_2914;
   }
   tmp_2912 += tmp_2913;
   result += (0.5) * tmp_2912;
   std::complex<double> tmp_2916;
   std::complex<double> tmp_2917;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2918;
      std::complex<double> tmp_2919;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2919 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpGlubarFuSuPL(gI1,gI2)
            )*CpGlubarFuSuPR(gI1,gI2);
      }
      tmp_2918 += tmp_2919;
      tmp_2917 += (MFu(gI1)) * tmp_2918;
   }
   tmp_2916 += tmp_2917;
   result += (0.5) * tmp_2916;
   std::complex<double> tmp_2920;
   std::complex<double> tmp_2921;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2922;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2922 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpGluconjSdFdPL(gI1,gI2
            ))*CpGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_2921 += tmp_2922;
   }
   tmp_2920 += tmp_2921;
   result += (0.5) * tmp_2920;
   std::complex<double> tmp_2923;
   std::complex<double> tmp_2924;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2925;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2925 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpGluconjSuFuPL(gI1,gI2
            ))*CpGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_2924 += tmp_2925;
   }
   tmp_2923 += tmp_2924;
   result += (0.5) * tmp_2923;
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpGluVGGluPR())*CpGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_2926;
   std::complex<double> tmp_2927;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2928;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2928 += AbsSqr(CpGlubarFdSdPR(gI1,gI2))*B1(p,MFd(gI1),MSd(
            gI2));
      }
      tmp_2927 += tmp_2928;
   }
   tmp_2926 += tmp_2927;
   result += (-0.25) * tmp_2926;
   std::complex<double> tmp_2929;
   std::complex<double> tmp_2930;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2931;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2931 += AbsSqr(CpGlubarFuSuPR(gI1,gI2))*B1(p,MFu(gI1),MSu(
            gI2));
      }
      tmp_2930 += tmp_2931;
   }
   tmp_2929 += tmp_2930;
   result += (-0.25) * tmp_2929;
   std::complex<double> tmp_2932;
   std::complex<double> tmp_2933;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2934;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2934 += AbsSqr(CpGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_2933 += tmp_2934;
   }
   tmp_2932 += tmp_2933;
   result += (-0.25) * tmp_2932;
   std::complex<double> tmp_2935;
   std::complex<double> tmp_2936;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2937;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2937 += AbsSqr(CpGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_2936 += tmp_2937;
   }
   tmp_2935 += tmp_2936;
   result += (-0.25) * tmp_2935;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_2938;
   std::complex<double> tmp_2939;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2940;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2940 += AbsSqr(CpGlubarFdSdPL(gI1,gI2))*B1(p,MFd(gI1),MSd(
            gI2));
      }
      tmp_2939 += tmp_2940;
   }
   tmp_2938 += tmp_2939;
   result += (-0.25) * tmp_2938;
   std::complex<double> tmp_2941;
   std::complex<double> tmp_2942;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2943;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2943 += AbsSqr(CpGlubarFuSuPL(gI1,gI2))*B1(p,MFu(gI1),MSu(
            gI2));
      }
      tmp_2942 += tmp_2943;
   }
   tmp_2941 += tmp_2942;
   result += (-0.25) * tmp_2941;
   std::complex<double> tmp_2944;
   std::complex<double> tmp_2945;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2946;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2946 += AbsSqr(CpGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_2945 += tmp_2946;
   }
   tmp_2944 += tmp_2945;
   result += (-0.25) * tmp_2944;
   std::complex<double> tmp_2947;
   std::complex<double> tmp_2948;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2949;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2949 += AbsSqr(CpGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_2948 += tmp_2949;
   }
   tmp_2947 += tmp_2948;
   result += (-0.25) * tmp_2947;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2950;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2951;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2951 += B0(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePL(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2950 += tmp_2951;
   }
   result += tmp_2950;
   std::complex<double> tmp_2952;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2953;
      std::complex<double> tmp_2954;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2954 += B0(p,MCha(gI1),MSe(gI2))*Conj(CpbarFvbarChaSePL(gO2,
            gI1,gI2))*CpbarFvbarChaSePR(gO1,gI1,gI2);
      }
      tmp_2953 += tmp_2954;
      tmp_2952 += (MCha(gI1)) * tmp_2953;
   }
   result += tmp_2952;
   std::complex<double> tmp_2955;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2956;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2956 += B0(p,MChi(gI2),MSv(gI1))*Conj(CpbarFvSvChiPL(gO2,gI1
            ,gI2))*CpbarFvSvChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2955 += tmp_2956;
   }
   result += tmp_2955;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2958 += B0(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePR(gO2,gI2))*
         CpbarFvconjVWmFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2957 += tmp_2958;
   result += (-4) * tmp_2957;
   std::complex<double> tmp_2959;
   std::complex<double> tmp_2960;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2960 += B0(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_2959 += tmp_2960;
   result += (-4) * tmp_2959;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2961;
   std::complex<double> tmp_2962;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2963;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2963 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePR(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePR(gO1,gI1,gI2);
      }
      tmp_2962 += tmp_2963;
   }
   tmp_2961 += tmp_2962;
   result += (-0.5) * tmp_2961;
   std::complex<double> tmp_2964;
   std::complex<double> tmp_2965;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2966;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2966 += B1(p,MCha(gI1),MSe(gI2))*Conj(CpbarFvbarChaSePR(gO2,
            gI1,gI2))*CpbarFvbarChaSePR(gO1,gI1,gI2);
      }
      tmp_2965 += tmp_2966;
   }
   tmp_2964 += tmp_2965;
   result += (-0.5) * tmp_2964;
   std::complex<double> tmp_2967;
   std::complex<double> tmp_2968;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2969;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2969 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpbarFvSvChiPR(gO2,gI1
            ,gI2))*CpbarFvSvChiPR(gO1,gI1,gI2);
      }
      tmp_2968 += tmp_2969;
   }
   tmp_2967 += tmp_2968;
   result += (-0.5) * tmp_2967;
   std::complex<double> tmp_2970;
   std::complex<double> tmp_2971;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2971 += B1(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePL(gO2,gI2))*
         CpbarFvconjVWmFePL(gO1,gI2);
   }
   tmp_2970 += tmp_2971;
   result += (-1) * tmp_2970;
   std::complex<double> tmp_2972;
   std::complex<double> tmp_2973;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2973 += B1(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPL(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2);
   }
   tmp_2972 += tmp_2973;
   result += (-1) * tmp_2972;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2974;
   std::complex<double> tmp_2975;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2976;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2976 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarFvconjHpmFePL(gO2
            ,gI1,gI2))*CpbarFvconjHpmFePL(gO1,gI1,gI2);
      }
      tmp_2975 += tmp_2976;
   }
   tmp_2974 += tmp_2975;
   result += (-0.5) * tmp_2974;
   std::complex<double> tmp_2977;
   std::complex<double> tmp_2978;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2979;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2979 += B1(p,MCha(gI1),MSe(gI2))*Conj(CpbarFvbarChaSePL(gO2,
            gI1,gI2))*CpbarFvbarChaSePL(gO1,gI1,gI2);
      }
      tmp_2978 += tmp_2979;
   }
   tmp_2977 += tmp_2978;
   result += (-0.5) * tmp_2977;
   std::complex<double> tmp_2980;
   std::complex<double> tmp_2981;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2982;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2982 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpbarFvSvChiPL(gO2,gI1
            ,gI2))*CpbarFvSvChiPL(gO1,gI1,gI2);
      }
      tmp_2981 += tmp_2982;
   }
   tmp_2980 += tmp_2981;
   result += (-0.5) * tmp_2980;
   std::complex<double> tmp_2983;
   std::complex<double> tmp_2984;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2984 += B1(p,MFe(gI2),MVWm)*Conj(CpbarFvconjVWmFePR(gO2,gI2))*
         CpbarFvconjVWmFePR(gO1,gI2);
   }
   tmp_2983 += tmp_2984;
   result += (-1) * tmp_2983;
   std::complex<double> tmp_2985;
   std::complex<double> tmp_2986;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2986 += B1(p,MFv(gI2),MVZ)*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPR(gO1,gI2);
   }
   tmp_2985 += tmp_2986;
   result += (-1) * tmp_2985;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2987;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2988;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2988 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2988 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2987 += tmp_2988;
   }
   result += tmp_2987;
   std::complex<double> tmp_2989;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2989 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2989;
   std::complex<double> tmp_2990;
   std::complex<double> tmp_2991;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2991 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2990 += tmp_2991;
   result += (0.5) * tmp_2990;
   std::complex<double> tmp_2992;
   std::complex<double> tmp_2993;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2994;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2994 += AbsSqr(CpVZhhAh(gI1,1 + gI2))*B00(p,MAh(1 + gI2),Mhh
            (gI1));
      }
      tmp_2993 += tmp_2994;
   }
   tmp_2992 += tmp_2993;
   result += (-4) * tmp_2992;
   std::complex<double> tmp_2995;
   std::complex<double> tmp_2996;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2997;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2997 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2996 += tmp_2997;
   }
   tmp_2995 += tmp_2996;
   result += (-4) * tmp_2995;
   std::complex<double> tmp_2998;
   std::complex<double> tmp_2999;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_3000;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3000 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3000 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2999 += tmp_3000;
   }
   tmp_2998 += tmp_2999;
   result += (0.5) * tmp_2998;
   std::complex<double> tmp_3001;
   std::complex<double> tmp_3002;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3002 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3001 += tmp_3002;
   result += (3) * tmp_3001;
   std::complex<double> tmp_3003;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3003 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3003;
   std::complex<double> tmp_3004;
   std::complex<double> tmp_3005;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3005 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3004 += tmp_3005;
   result += (3) * tmp_3004;
   std::complex<double> tmp_3006;
   std::complex<double> tmp_3007;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3008;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3008 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_3007 += tmp_3008;
   }
   tmp_3006 += tmp_3007;
   result += (-12) * tmp_3006;
   std::complex<double> tmp_3009;
   std::complex<double> tmp_3010;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3011;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3011 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_3010 += tmp_3011;
   }
   tmp_3009 += tmp_3010;
   result += (-4) * tmp_3009;
   std::complex<double> tmp_3012;
   std::complex<double> tmp_3013;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3014;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3014 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_3013 += tmp_3014;
   }
   tmp_3012 += tmp_3013;
   result += (-12) * tmp_3012;
   std::complex<double> tmp_3015;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3015 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3015;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3016;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3016 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3016;
   std::complex<double> tmp_3017;
   std::complex<double> tmp_3018;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3018 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3017 += tmp_3018;
   result += (0.5) * tmp_3017;
   std::complex<double> tmp_3019;
   std::complex<double> tmp_3020;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3021;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3021 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3020 += tmp_3021;
   }
   tmp_3019 += tmp_3020;
   result += (-4) * tmp_3019;
   std::complex<double> tmp_3022;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_3023;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3023 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_3023 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_3022 += tmp_3023;
   }
   result += tmp_3022;
   std::complex<double> tmp_3024;
   std::complex<double> tmp_3025;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3025 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3024 += tmp_3025;
   result += (3) * tmp_3024;
   std::complex<double> tmp_3026;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3026 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3026;
   std::complex<double> tmp_3027;
   std::complex<double> tmp_3028;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3028 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3027 += tmp_3028;
   result += (3) * tmp_3027;
   std::complex<double> tmp_3029;
   std::complex<double> tmp_3030;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3031;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3031 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3030 += tmp_3031;
   }
   tmp_3029 += tmp_3030;
   result += (-12) * tmp_3029;
   std::complex<double> tmp_3032;
   std::complex<double> tmp_3033;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3033 += AbsSqr(CpconjVWmHpmhh(1,gI2))*B00(p,Mhh(gI2),MHpm(1));
   }
   tmp_3032 += tmp_3033;
   result += (-4) * tmp_3032;
   std::complex<double> tmp_3034;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3034 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3034;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3035;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3036;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3036 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3035 += tmp_3036;
   }
   result += tmp_3035;
   std::complex<double> tmp_3037;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3038;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3038 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3037 += tmp_3038;
   }
   result += tmp_3037;
   std::complex<double> tmp_3039;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3040;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3040 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3039 += tmp_3040;
   }
   result += tmp_3039;
   std::complex<double> tmp_3041;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3042;
      std::complex<double> tmp_3043;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3043 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3042 += tmp_3043;
      tmp_3041 += (MFe(gI1)) * tmp_3042;
   }
   result += tmp_3041;
   std::complex<double> tmp_3044;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3045;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3045 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3044 += tmp_3045;
   }
   result += tmp_3044;
   std::complex<double> tmp_3046;
   std::complex<double> tmp_3047;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3047 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3046 += tmp_3047;
   result += (-4) * tmp_3046;
   std::complex<double> tmp_3048;
   std::complex<double> tmp_3049;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3049 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3048 += tmp_3049;
   result += (-4) * tmp_3048;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3050;
   std::complex<double> tmp_3051;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3052;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3052 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3051 += tmp_3052;
   }
   tmp_3050 += tmp_3051;
   result += (-0.5) * tmp_3050;
   std::complex<double> tmp_3053;
   std::complex<double> tmp_3054;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3055;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3055 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPR(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3054 += tmp_3055;
   }
   tmp_3053 += tmp_3054;
   result += (-0.5) * tmp_3053;
   std::complex<double> tmp_3056;
   std::complex<double> tmp_3057;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3058;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3058 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3057 += tmp_3058;
   }
   tmp_3056 += tmp_3057;
   result += (-0.5) * tmp_3056;
   std::complex<double> tmp_3059;
   std::complex<double> tmp_3060;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3061;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3061 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_3060 += tmp_3061;
   }
   tmp_3059 += tmp_3060;
   result += (-0.5) * tmp_3059;
   std::complex<double> tmp_3062;
   std::complex<double> tmp_3063;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3064;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3064 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3063 += tmp_3064;
   }
   tmp_3062 += tmp_3063;
   result += (-0.5) * tmp_3062;
   std::complex<double> tmp_3065;
   std::complex<double> tmp_3066;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3066 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3065 += tmp_3066;
   result += (-1) * tmp_3065;
   std::complex<double> tmp_3067;
   std::complex<double> tmp_3068;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3068 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3067 += tmp_3068;
   result += (-1) * tmp_3067;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3069;
   std::complex<double> tmp_3070;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3071;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3071 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3070 += tmp_3071;
   }
   tmp_3069 += tmp_3070;
   result += (-0.5) * tmp_3069;
   std::complex<double> tmp_3072;
   std::complex<double> tmp_3073;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3074;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3074 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3073 += tmp_3074;
   }
   tmp_3072 += tmp_3073;
   result += (-0.5) * tmp_3072;
   std::complex<double> tmp_3075;
   std::complex<double> tmp_3076;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3077;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3077 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3076 += tmp_3077;
   }
   tmp_3075 += tmp_3076;
   result += (-0.5) * tmp_3075;
   std::complex<double> tmp_3078;
   std::complex<double> tmp_3079;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3080;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3080 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3079 += tmp_3080;
   }
   tmp_3078 += tmp_3079;
   result += (-0.5) * tmp_3078;
   std::complex<double> tmp_3081;
   std::complex<double> tmp_3082;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3083;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3083 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3082 += tmp_3083;
   }
   tmp_3081 += tmp_3082;
   result += (-0.5) * tmp_3081;
   std::complex<double> tmp_3084;
   std::complex<double> tmp_3085;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3085 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3084 += tmp_3085;
   result += (-1) * tmp_3084;
   std::complex<double> tmp_3086;
   std::complex<double> tmp_3087;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3087 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3086 += tmp_3087;
   result += (-1) * tmp_3086;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3088;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3089;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3089 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3088 += tmp_3089;
   }
   result += tmp_3088;
   std::complex<double> tmp_3090;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3091;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3091 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3090 += tmp_3091;
   }
   result += tmp_3090;
   std::complex<double> tmp_3092;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3093;
      std::complex<double> tmp_3094;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3094 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3093 += tmp_3094;
      tmp_3092 += (MFd(gI1)) * tmp_3093;
   }
   result += tmp_3092;
   std::complex<double> tmp_3095;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3096;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3096 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3095 += tmp_3096;
   }
   result += tmp_3095;
   std::complex<double> tmp_3097;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3098;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3098 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3097 += tmp_3098;
   }
   result += tmp_3097;
   std::complex<double> tmp_3099;
   std::complex<double> tmp_3100;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3100 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3099 += tmp_3100;
   result += (-4) * tmp_3099;
   std::complex<double> tmp_3101;
   std::complex<double> tmp_3102;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3102 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3101 += tmp_3102;
   result += (-4) * tmp_3101;
   std::complex<double> tmp_3103;
   std::complex<double> tmp_3104;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3104 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,0))*
         CpbarFdSdGluPR(gO1,gI1,0);
   }
   tmp_3103 += tmp_3104;
   result += (1.3333333333333333*MGlu) * tmp_3103;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3105;
   std::complex<double> tmp_3106;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3107;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3107 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3106 += tmp_3107;
   }
   tmp_3105 += tmp_3106;
   result += (-0.5) * tmp_3105;
   std::complex<double> tmp_3108;
   std::complex<double> tmp_3109;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3110;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3110 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3109 += tmp_3110;
   }
   tmp_3108 += tmp_3109;
   result += (-0.5) * tmp_3108;
   std::complex<double> tmp_3111;
   std::complex<double> tmp_3112;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3113 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3112 += tmp_3113;
   }
   tmp_3111 += tmp_3112;
   result += (-0.5) * tmp_3111;
   std::complex<double> tmp_3114;
   std::complex<double> tmp_3115;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3115 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,0))*
         CpbarFdSdGluPR(gO1,gI1,0);
   }
   tmp_3114 += tmp_3115;
   result += (-0.6666666666666666) * tmp_3114;
   std::complex<double> tmp_3116;
   std::complex<double> tmp_3117;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3118;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3118 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPR(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3117 += tmp_3118;
   }
   tmp_3116 += tmp_3117;
   result += (-0.5) * tmp_3116;
   std::complex<double> tmp_3119;
   std::complex<double> tmp_3120;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3121;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3121 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3120 += tmp_3121;
   }
   tmp_3119 += tmp_3120;
   result += (-0.5) * tmp_3119;
   std::complex<double> tmp_3122;
   std::complex<double> tmp_3123;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3123 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3122 += tmp_3123;
   result += (-1) * tmp_3122;
   std::complex<double> tmp_3124;
   std::complex<double> tmp_3125;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3125 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3124 += tmp_3125;
   result += (-1) * tmp_3124;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3126;
   std::complex<double> tmp_3127;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3128;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3128 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3127 += tmp_3128;
   }
   tmp_3126 += tmp_3127;
   result += (-0.5) * tmp_3126;
   std::complex<double> tmp_3129;
   std::complex<double> tmp_3130;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3131;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3131 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3130 += tmp_3131;
   }
   tmp_3129 += tmp_3130;
   result += (-0.5) * tmp_3129;
   std::complex<double> tmp_3132;
   std::complex<double> tmp_3133;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3134;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3134 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3133 += tmp_3134;
   }
   tmp_3132 += tmp_3133;
   result += (-0.5) * tmp_3132;
   std::complex<double> tmp_3135;
   std::complex<double> tmp_3136;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3136 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,0))*
         CpbarFdSdGluPL(gO1,gI1,0);
   }
   tmp_3135 += tmp_3136;
   result += (-0.6666666666666666) * tmp_3135;
   std::complex<double> tmp_3137;
   std::complex<double> tmp_3138;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3139;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3139 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_3138 += tmp_3139;
   }
   tmp_3137 += tmp_3138;
   result += (-0.5) * tmp_3137;
   std::complex<double> tmp_3140;
   std::complex<double> tmp_3141;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3142;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3142 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3141 += tmp_3142;
   }
   tmp_3140 += tmp_3141;
   result += (-0.5) * tmp_3140;
   std::complex<double> tmp_3143;
   std::complex<double> tmp_3144;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3144 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3143 += tmp_3144;
   result += (-1) * tmp_3143;
   std::complex<double> tmp_3145;
   std::complex<double> tmp_3146;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3146 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3145 += tmp_3146;
   result += (-1) * tmp_3145;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3147;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3148;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3148 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3147 += tmp_3148;
   }
   result += tmp_3147;
   std::complex<double> tmp_3149;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3150;
      std::complex<double> tmp_3151;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3151 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3150 += tmp_3151;
      tmp_3149 += (MCha(gI1)) * tmp_3150;
   }
   result += tmp_3149;
   std::complex<double> tmp_3152;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3153;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3153 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3152 += tmp_3153;
   }
   result += tmp_3152;
   std::complex<double> tmp_3154;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3155;
      std::complex<double> tmp_3156;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3156 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3155 += tmp_3156;
      tmp_3154 += (MFu(gI1)) * tmp_3155;
   }
   result += tmp_3154;
   std::complex<double> tmp_3157;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3158;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3158 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3157 += tmp_3158;
   }
   result += tmp_3157;
   std::complex<double> tmp_3159;
   std::complex<double> tmp_3160;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3160 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3159 += tmp_3160;
   result += (-4) * tmp_3159;
   std::complex<double> tmp_3161;
   std::complex<double> tmp_3162;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3162 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3161 += tmp_3162;
   result += (-4) * tmp_3161;
   std::complex<double> tmp_3163;
   std::complex<double> tmp_3164;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3164 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3163 += tmp_3164;
   result += (-4) * tmp_3163;
   std::complex<double> tmp_3165;
   std::complex<double> tmp_3166;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3166 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,0))*
         CpbarFuSuGluPR(gO1,gI1,0);
   }
   tmp_3165 += tmp_3166;
   result += (1.3333333333333333*MGlu) * tmp_3165;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3167;
   std::complex<double> tmp_3168;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3169;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3169 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3168 += tmp_3169;
   }
   tmp_3167 += tmp_3168;
   result += (-0.5) * tmp_3167;
   std::complex<double> tmp_3170;
   std::complex<double> tmp_3171;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3172;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3172 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPR(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3171 += tmp_3172;
   }
   tmp_3170 += tmp_3171;
   result += (-0.5) * tmp_3170;
   std::complex<double> tmp_3173;
   std::complex<double> tmp_3174;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3175;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3175 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3174 += tmp_3175;
   }
   tmp_3173 += tmp_3174;
   result += (-0.5) * tmp_3173;
   std::complex<double> tmp_3176;
   std::complex<double> tmp_3177;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3178;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3178 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3177 += tmp_3178;
   }
   tmp_3176 += tmp_3177;
   result += (-0.5) * tmp_3176;
   std::complex<double> tmp_3179;
   std::complex<double> tmp_3180;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3180 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,0))*
         CpbarFuSuGluPR(gO1,gI1,0);
   }
   tmp_3179 += tmp_3180;
   result += (-0.6666666666666666) * tmp_3179;
   std::complex<double> tmp_3181;
   std::complex<double> tmp_3182;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3183;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3183 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3182 += tmp_3183;
   }
   tmp_3181 += tmp_3182;
   result += (-0.5) * tmp_3181;
   std::complex<double> tmp_3184;
   std::complex<double> tmp_3185;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3185 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3184 += tmp_3185;
   result += (-1) * tmp_3184;
   std::complex<double> tmp_3186;
   std::complex<double> tmp_3187;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3187 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3186 += tmp_3187;
   result += (-1) * tmp_3186;
   std::complex<double> tmp_3188;
   std::complex<double> tmp_3189;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3189 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3188 += tmp_3189;
   result += (-1) * tmp_3188;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3190;
   std::complex<double> tmp_3191;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3192;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3192 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3191 += tmp_3192;
   }
   tmp_3190 += tmp_3191;
   result += (-0.5) * tmp_3190;
   std::complex<double> tmp_3193;
   std::complex<double> tmp_3194;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3195;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3195 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3194 += tmp_3195;
   }
   tmp_3193 += tmp_3194;
   result += (-0.5) * tmp_3193;
   std::complex<double> tmp_3196;
   std::complex<double> tmp_3197;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3198;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3198 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3197 += tmp_3198;
   }
   tmp_3196 += tmp_3197;
   result += (-0.5) * tmp_3196;
   std::complex<double> tmp_3199;
   std::complex<double> tmp_3200;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3201;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3201 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3200 += tmp_3201;
   }
   tmp_3199 += tmp_3200;
   result += (-0.5) * tmp_3199;
   std::complex<double> tmp_3202;
   std::complex<double> tmp_3203;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3203 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,0))*
         CpbarFuSuGluPL(gO1,gI1,0);
   }
   tmp_3202 += tmp_3203;
   result += (-0.6666666666666666) * tmp_3202;
   std::complex<double> tmp_3204;
   std::complex<double> tmp_3205;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3206;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3206 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3205 += tmp_3206;
   }
   tmp_3204 += tmp_3205;
   result += (-0.5) * tmp_3204;
   std::complex<double> tmp_3207;
   std::complex<double> tmp_3208;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3208 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3207 += tmp_3208;
   result += (-1) * tmp_3207;
   std::complex<double> tmp_3209;
   std::complex<double> tmp_3210;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3210 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3209 += tmp_3210;
   result += (-1) * tmp_3209;
   std::complex<double> tmp_3211;
   std::complex<double> tmp_3212;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3212 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3211 += tmp_3212;
   result += (-1) * tmp_3211;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3213;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3214;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3214 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3213 += tmp_3214;
   }
   result += tmp_3213;
   std::complex<double> tmp_3215;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3216;
      std::complex<double> tmp_3217;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3217 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3216 += tmp_3217;
      tmp_3215 += (MCha(gI1)) * tmp_3216;
   }
   result += tmp_3215;
   std::complex<double> tmp_3218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3219;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3219 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3218 += tmp_3219;
   }
   result += tmp_3218;
   std::complex<double> tmp_3220;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3221;
      std::complex<double> tmp_3222;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3222 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3221 += tmp_3222;
      tmp_3220 += (MFu(gI1)) * tmp_3221;
   }
   result += tmp_3220;
   std::complex<double> tmp_3223;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3224;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3224 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3223 += tmp_3224;
   }
   result += tmp_3223;
   std::complex<double> tmp_3225;
   std::complex<double> tmp_3226;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3226 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3225 += tmp_3226;
   result += (-4) * tmp_3225;
   std::complex<double> tmp_3227;
   std::complex<double> tmp_3228;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3228 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3227 += tmp_3228;
   result += (-4) * tmp_3227;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3230 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3229 += tmp_3230;
   result += (-4) * tmp_3229;
   std::complex<double> tmp_3231;
   std::complex<double> tmp_3232;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3232 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3231 += tmp_3232;
   result += (1.3333333333333333*MGlu) * tmp_3231;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3233;
   std::complex<double> tmp_3234;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3235;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3235 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3234 += tmp_3235;
   }
   tmp_3233 += tmp_3234;
   result += (-0.5) * tmp_3233;
   std::complex<double> tmp_3236;
   std::complex<double> tmp_3237;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3238;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3238 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3237 += tmp_3238;
   }
   tmp_3236 += tmp_3237;
   result += (-0.5) * tmp_3236;
   std::complex<double> tmp_3239;
   std::complex<double> tmp_3240;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3241;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3241 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3240 += tmp_3241;
   }
   tmp_3239 += tmp_3240;
   result += (-0.5) * tmp_3239;
   std::complex<double> tmp_3242;
   std::complex<double> tmp_3243;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3244;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3244 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3243 += tmp_3244;
   }
   tmp_3242 += tmp_3243;
   result += (-0.5) * tmp_3242;
   std::complex<double> tmp_3245;
   std::complex<double> tmp_3246;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3246 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,0))*
         CpbarUFuSuGluPR(gO1,gI1,0);
   }
   tmp_3245 += tmp_3246;
   result += (-0.6666666666666666) * tmp_3245;
   std::complex<double> tmp_3247;
   std::complex<double> tmp_3248;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3249;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3249 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3248 += tmp_3249;
   }
   tmp_3247 += tmp_3248;
   result += (-0.5) * tmp_3247;
   std::complex<double> tmp_3250;
   std::complex<double> tmp_3251;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3251 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3250 += tmp_3251;
   result += (-1) * tmp_3250;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3253 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3252 += tmp_3253;
   result += (-1) * tmp_3252;
   std::complex<double> tmp_3254;
   std::complex<double> tmp_3255;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3255 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3254 += tmp_3255;
   result += (-1) * tmp_3254;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3256;
   std::complex<double> tmp_3257;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3258;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3258 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3257 += tmp_3258;
   }
   tmp_3256 += tmp_3257;
   result += (-0.5) * tmp_3256;
   std::complex<double> tmp_3259;
   std::complex<double> tmp_3260;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3261;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3261 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3260 += tmp_3261;
   }
   tmp_3259 += tmp_3260;
   result += (-0.5) * tmp_3259;
   std::complex<double> tmp_3262;
   std::complex<double> tmp_3263;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3264;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3264 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3263 += tmp_3264;
   }
   tmp_3262 += tmp_3263;
   result += (-0.5) * tmp_3262;
   std::complex<double> tmp_3265;
   std::complex<double> tmp_3266;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3267;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3267 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3266 += tmp_3267;
   }
   tmp_3265 += tmp_3266;
   result += (-0.5) * tmp_3265;
   std::complex<double> tmp_3268;
   std::complex<double> tmp_3269;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3269 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,0))*
         CpbarUFuSuGluPL(gO1,gI1,0);
   }
   tmp_3268 += tmp_3269;
   result += (-0.6666666666666666) * tmp_3268;
   std::complex<double> tmp_3270;
   std::complex<double> tmp_3271;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3272;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3272 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3271 += tmp_3272;
   }
   tmp_3270 += tmp_3271;
   result += (-0.5) * tmp_3270;
   std::complex<double> tmp_3273;
   std::complex<double> tmp_3274;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3274 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3273 += tmp_3274;
   result += (-1) * tmp_3273;
   std::complex<double> tmp_3275;
   std::complex<double> tmp_3276;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3276 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3275 += tmp_3276;
   result += (-1) * tmp_3275;
   std::complex<double> tmp_3277;
   std::complex<double> tmp_3278;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3278 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3277 += tmp_3278;
   result += (-1) * tmp_3277;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   result += 4*A0(MVWm)*CpUhhconjVWmVWm(gO1);
   result += 2*A0(MVZ)*CpUhhVZVZ(gO1);
   std::complex<double> tmp_3279;
   std::complex<double> tmp_3280;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3280 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3279 += tmp_3280;
   result += (-1) * tmp_3279;
   std::complex<double> tmp_3281;
   std::complex<double> tmp_3282;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3282 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_3281 += tmp_3282;
   result += (2) * tmp_3281;
   std::complex<double> tmp_3283;
   std::complex<double> tmp_3284;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3284 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3283 += tmp_3284;
   result += (-0.5) * tmp_3283;
   std::complex<double> tmp_3285;
   std::complex<double> tmp_3286;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3286 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3285 += tmp_3286;
   result += (-1) * tmp_3285;
   std::complex<double> tmp_3287;
   std::complex<double> tmp_3288;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3288 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3287 += tmp_3288;
   result += (-0.5) * tmp_3287;
   std::complex<double> tmp_3289;
   std::complex<double> tmp_3290;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3290 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3289 += tmp_3290;
   result += (6) * tmp_3289;
   std::complex<double> tmp_3291;
   std::complex<double> tmp_3292;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3292 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3291 += tmp_3292;
   result += (2) * tmp_3291;
   std::complex<double> tmp_3293;
   std::complex<double> tmp_3294;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3294 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3293 += tmp_3294;
   result += (6) * tmp_3293;
   std::complex<double> tmp_3295;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      tmp_3295 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_3295;
   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3297 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3296 += tmp_3297;
   result += (-3) * tmp_3296;
   std::complex<double> tmp_3298;
   std::complex<double> tmp_3299;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3299 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3298 += tmp_3299;
   result += (-1) * tmp_3298;
   std::complex<double> tmp_3300;
   std::complex<double> tmp_3301;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3301 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3300 += tmp_3301;
   result += (-3) * tmp_3300;

   return result * oneOver16PiSqr;

}


void CLASSNAME::calculate_MSu_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(2,2));
   sf_data.mr2 = Re(mu2(2,2));
   sf_data.yf  = Re(Yu(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYu(2,2));
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::up];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::up];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::up];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}

void CLASSNAME::calculate_MSd_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(2,2));
   sf_data.mr2 = Re(md2(2,2));
   sf_data.yf  = Re(Yd(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYd(2,2));
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::down];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::down];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::down];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}

void CLASSNAME::calculate_MSv_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(2,2));
   sf_data.mr2 = 0.;
   sf_data.yf  = 0.;
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = 0.;
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::neutrino];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::neutrino];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::neutrino];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}

void CLASSNAME::calculate_MSe_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(2,2));
   sf_data.mr2 = Re(me2(2,2));
   sf_data.yf  = Re(Ye(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYe(2,2));
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::electron];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::electron];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::electron];

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}


Eigen::Matrix<double,3,3> CLASSNAME::self_energy_hh_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double as = Sqr(gs) / (4.0 * Pi);
   const double rmt = MFu(2);
   const double rmtsq = Sqr(rmt);
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = ((Sqr(vd) + Sqr(vu))*(0.5*Kappa*Lambdax*Sqr(vS) + 0.7071067811865475*vS*TLambdax))/(vd*vu);
   const double cotb = 1.0 / tanb;
   const double rmb = MFd(2);
   const double rmbsq = Sqr(rmb);
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_2l += self_energy_higgs_2loop_at_as_nmssm(
         rmt, mg, mst1sq, mst2sq, sxt, cxt,
         scalesq, tanb, vev2, lam, svev, as, amu);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_2l += self_energy_higgs_2loop_ab_as_nmssm(
         rmb, mg, msb1sq, msb2sq, sxb, cxb,
         scalesq, cotb, vev2, lam, svev, as, amu);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq,
         msb2sq, sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   return self_energy_2l;
}

Eigen::Matrix<double,3,3> CLASSNAME::self_energy_Ah_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double as = Sqr(gs) / (4.0 * Pi);
   const double rmt = MFu(2);
   const double rmtsq = Sqr(rmt);
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = ((Sqr(vd) + Sqr(vu))*(0.5*Kappa*Lambdax*Sqr(vS) + 0.7071067811865475*vS*TLambdax))/(vd*vu);
   const double cotb = 1.0 / tanb;
   const double rmb = MFd(2);
   const double rmbsq = Sqr(rmb);
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_2l += self_energy_pseudoscalar_2loop_at_as_nmssm(
         rmt, mg, mst1sq, mst2sq, sxt, cxt,
         scalesq, tanb, vev2, lam, svev, as, amu);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_nmssm(
         rmb, mg, msb1sq, msb2sq, sxb, cxb,
         scalesq, cotb, vev2, lam, svev, as, amu);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
         sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   return self_energy_2l;
}



Eigen::Matrix<double,3,1> CLASSNAME::tadpole_hh_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double rmtsq = Sqr(MFu(2));
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = ((Sqr(vd) + Sqr(vu))*(0.5*Kappa*Lambdax*Sqr(vS) + 0.7071067811865475*vS*TLambdax))/(vd*vu);
   const double cotbeta = 1.0 / tanb;
   const double rmbsq = Sqr(MFd(2));
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,1> tadpole_2l(Eigen::Matrix<double,3,1>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      tadpole_2l += tadpole_higgs_2loop_at_as_nmssm(
         rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
         amu, tanb, vev2, gs, svev);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      tadpole_2l.head<2>() += tadpole_higgs_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
         sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      tadpole_2l += tadpole_higgs_2loop_ab_as_nmssm(
         rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
         amu, cotbeta, vev2, gs, svev);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      tadpole_2l.head<2>() += tadpole_higgs_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   tadpole_2l(0) *= vd;
   tadpole_2l(1) *= vu;
   tadpole_2l(2) *= vS;

   return tadpole_2l;
}




void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVG) old_MVG(MVG), new_MVG(MVG);

   do {
      PHYSICAL(MVG) = 0.;

      new_MVG = PHYSICAL(MVG);
      diff = MaxRelDiff(new_MVG, old_MVG);
      old_MVG = new_MVG;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::VG);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::VG);
}

void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MGlu) old_MGlu(MGlu), new_MGlu(MGlu);

   do {
      const double M_tree(MGlu);
      const double p = old_MGlu;
      const double self_energy_1  = Re(self_energy_Glu_1(p));
      const double self_energy_PL = Re(self_energy_Glu_PL(p));
      const double self_energy_PR = Re(self_energy_Glu_PR(p));
      const auto M_loop = M_tree - self_energy_1 - M_tree * (
         self_energy_PL + self_energy_PR);
      PHYSICAL(MGlu) = calculate_singlet_mass(M_loop);

      new_MGlu = PHYSICAL(MGlu);
      diff = MaxRelDiff(new_MGlu, old_MGlu);
      old_MGlu = new_MGlu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Glu);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Glu);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFv) old_MFv(MFv), new_MFv(MFv);

   do {
      PHYSICAL(MFv).setConstant(0.);

      new_MFv = PHYSICAL(MFv);
      diff = MaxRelDiff(new_MFv, old_MFv);
      old_MFv = new_MFv;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Fv);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Fv);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVP) old_MVP(MVP), new_MVP(MVP);

   do {
      PHYSICAL(MVP) = 0.;

      new_MVP = PHYSICAL(MVP);
      diff = MaxRelDiff(new_MVP, old_MVP);
      old_MVP = new_MVP;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::VP);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::VP);
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::VZ))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVZ) old_MVZ(MVZ), new_MVZ(MVZ);

   do {
      const double M_tree(Sqr(MVZ));
      const double p = old_MVZ;
      const double self_energy = Re(self_energy_VZ(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(NMSSMtower_info::VZ);

      PHYSICAL(MVZ) = AbsSqrt(mass_sqr);

      new_MVZ = PHYSICAL(MVZ);
      diff = MaxRelDiff(new_MVZ, old_MVZ);
      old_MVZ = new_MVZ;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::VZ);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::VZ);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Sd))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSd) old_MSd(MSd), new_MSd(MSd);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sd());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSd(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Sd(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZD;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZD, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Sd,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZD);
         #endif

         PHYSICAL(MSd(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZD) = mix_ZD;
      }

      new_MSd = PHYSICAL(MSd);
      diff = MaxRelDiff(new_MSd, old_MSd);
      old_MSd = new_MSd;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Sd);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Sd);
}

void CLASSNAME::calculate_MSv_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Sv))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSv) old_MSv(MSv), new_MSv(MSv);

   do {
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Sv());

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MSv(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Sv(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZV;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZV, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Sv,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZV);
         #endif

         PHYSICAL(MSv(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZV) = mix_ZV;
      }

      new_MSv = PHYSICAL(MSv);
      diff = MaxRelDiff(new_MSv, old_MSv);
      old_MSv = new_MSv;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Sv);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Sv);
}

void CLASSNAME::calculate_MSu_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Su))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSu) old_MSu(MSu), new_MSu(MSu);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Su());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSu(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Su(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZU;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZU, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Su,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZU);
         #endif

         PHYSICAL(MSu(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZU) = mix_ZU;
      }

      new_MSu = PHYSICAL(MSu);
      diff = MaxRelDiff(new_MSu, old_MSu);
      old_MSu = new_MSu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Su);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Su);
}

void CLASSNAME::calculate_MSe_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Se))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MSe) old_MSe(MSe), new_MSe(MSe);

   do {
      Eigen::Matrix<double,6,6> self_energy;
      const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Se());

      for (unsigned es = 0; es < 6; ++es) {
         const double p = Abs(old_MSe(es));
         for (unsigned i1 = 0; i1 < 6; ++i1) {
            for (unsigned i2 = i1; i2 < 6; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Se(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,6,6> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,6,1> eigen_values;
         Eigen::Matrix<double,6,6> mix_ZE;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZE, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Se,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZE);
         #endif

         PHYSICAL(MSe(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZE) = mix_ZE;
      }

      new_MSe = PHYSICAL(MSe);
      diff = MaxRelDiff(new_MSe, old_MSe);
      old_MSe = new_MSe;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Se);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Se);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::hh))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_hh());

      // two-loop Higgs self-energy contributions
      Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,
         3>::Zero());
      if (pole_mass_loop_order > 1) {
         self_energy_2l = self_energy_hh_2loop();
         for (unsigned i = 0; i < 3; i++) {
            for (unsigned k = 0; k < 3; k++) {
               if (!std::isfinite(self_energy_2l(i,k))) {
                  self_energy_2l(i,k) = 0.;
                  problems.flag_bad_mass(
                     NMSSMtower_info::hh);
               }
            }
         }
      }

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_Mhh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_hh(p,i1,i2
                  ));
            }
         }

         self_energy += self_energy_2l;

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::hh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH);
         #endif

         PHYSICAL(Mhh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZH) = mix_ZH;
      }

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Ah))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Ah());

      // two-loop Higgs self-energy contributions
      Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,
         3>::Zero());
      if (pole_mass_loop_order > 1) {
         self_energy_2l = self_energy_Ah_2loop();
         for (unsigned i = 0; i < 3; i++) {
            for (unsigned k = 0; k < 3; k++) {
               if (!std::isfinite(self_energy_2l(i,k))) {
                  self_energy_2l(i,k) = 0.;
                  problems.flag_bad_mass(
                     NMSSMtower_info::Ah);
               }
            }
         }
      }

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MAh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2
                  ));
            }
         }

         self_energy += self_energy_2l;

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA);
         #endif

         PHYSICAL(MAh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 1)
            PHYSICAL(ZA) = mix_ZA;
      }

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Ah);
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::Hpm))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MHpm) old_MHpm(MHpm), new_MHpm(MHpm);

   do {
      Eigen::Matrix<double,2,2> self_energy;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Hpm());

      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MHpm(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = i1; i2 < 2; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Hpm(p,i1,
                  i2));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,2,2> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Hpm,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP);
         #endif

         PHYSICAL(MHpm(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 1)
            PHYSICAL(ZP) = mix_ZP;
      }

      new_MHpm = PHYSICAL(MHpm);
      diff = MaxRelDiff(new_MHpm, old_MHpm);
      old_MHpm = new_MHpm;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Hpm);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MChi) old_MChi(MChi), new_MChi(MChi);

   do {
      Eigen::Matrix<double,5,5> self_energy_1;
      Eigen::Matrix<double,5,5> self_energy_PL;
      Eigen::Matrix<double,5,5> self_energy_PR;
      const Eigen::Matrix<double,5,5> M_tree(get_mass_matrix_Chi());
      for (unsigned es = 0; es < 5; ++es) {
         const double p = Abs(old_MChi(es));
         for (unsigned i1 = 0; i1 < 5; ++i1) {
            for (unsigned i2 = 0; i2 < 5; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Chi_1(p
                  ,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Chi_PL(
                  p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Chi_PR(
                  p,i1,i2));
            }
         }
         const Eigen::Matrix<double,5,5> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,5,5> M_loop(M_tree + 0.5 * (
            delta_M + delta_M.transpose()));
         Eigen::Array<double,5,1> eigen_values;
         decltype(ZN) mix_ZN;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_symmetric(M_loop, eigen_values,
               mix_ZN, eigenvalue_error);
            problems.flag_bad_mass(NMSSMtower_info::Chi,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_symmetric(M_loop, eigen_values,
               mix_ZN);
         #endif
         if (es == 0)
            PHYSICAL(ZN) = mix_ZN;
         PHYSICAL(MChi(es)) = Abs(eigen_values(es));
      }

      new_MChi = PHYSICAL(MChi);
      diff = MaxRelDiff(new_MChi, old_MChi);
      old_MChi = new_MChi;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Chi);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Chi);
}

void CLASSNAME::calculate_MCha_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MCha) old_MCha(MCha), new_MCha(MCha);

   do {
      Eigen::Matrix<double,2,2> self_energy_1;
      Eigen::Matrix<double,2,2> self_energy_PL;
      Eigen::Matrix<double,2,2> self_energy_PR;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha());
      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_MCha(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = 0; i2 < 2; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Cha_1(p
                  ,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Cha_PL(
                  p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Cha_PR(
                  p,i1,i2));
            }
         }
         const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
         Eigen::Array<double,2,1> eigen_values;
         decltype(UM) mix_UM;
         decltype(UP) mix_UP;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_UM, mix_UP,
            eigenvalue_error);
         problems.flag_bad_mass(NMSSMtower_info::Cha,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_UM, mix_UP);
      #endif
         if (es == 0) {
            PHYSICAL(UM) = mix_UM;
            PHYSICAL(UP) = mix_UP;
         }
         PHYSICAL(MCha(es)) = Abs(eigen_values(es));
      }

      new_MCha = PHYSICAL(MCha);
      diff = MaxRelDiff(new_MCha, old_MCha);
      old_MCha = new_MCha;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Cha);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Cha);
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFe) old_MFe(MFe), new_MFe(MFe);

   do {
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFe(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = 0; i2 < 3; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fe_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fe_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fe_PR(p
                  ,i1,i2));
            }
         }
         const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZEL) mix_ZEL;
         decltype(ZER) mix_ZER;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER,
            eigenvalue_error);
         problems.flag_bad_mass(NMSSMtower_info::Fe,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER);
      #endif
         if (es == 0) {
            PHYSICAL(ZEL) = mix_ZEL;
            PHYSICAL(ZER) = mix_ZER;
         }
         PHYSICAL(MFe(es)) = Abs(eigen_values(es));
      }

      new_MFe = PHYSICAL(MFe);
      diff = MaxRelDiff(new_MFe, old_MFe);
      old_MFe = new_MFe;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Fe);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Fe);
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFd) old_MFd(MFd), new_MFd(MFd);

   do {
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFd(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = 0; i2 < 3; ++i2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fd_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fd_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fd_PR(p
                  ,i1,i2));
            }
         }
         const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZDL) mix_ZDL;
         decltype(ZDR) mix_ZDR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR,
            eigenvalue_error);
         problems.flag_bad_mass(NMSSMtower_info::Fd,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR);
      #endif
         if (es == 0) {
            PHYSICAL(ZDL) = mix_ZDL;
            PHYSICAL(ZDR) = mix_ZDR;
         }
         PHYSICAL(MFd(es)) = Abs(eigen_values(es));
      }

      new_MFd = PHYSICAL(MFd);
      diff = MaxRelDiff(new_MFd, old_MFd);
      old_MFd = new_MFd;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Fd);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Fd);
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MFu) old_MFu(MFu), new_MFu(MFu);

   do {
      double qcd_1l = 0.;

      {
         const double currentScale = get_scale();
         qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))
            /Sqr(currentScale)))*Sqr(g3);
      }

      double qcd_2l = 0.;

      if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
         const double currentScale = get_scale();
         qcd_2l = -0.005191204615668296*Power(g3,4) -
            0.0032883224409535764*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(2))
            ) - 0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(currentScale,2)
            /Sqr(MFu(2))));
      }

      double qcd_3l = 0.;

      if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
         const double currentScale = get_scale();
         qcd_3l = 0;
      }

      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MFu(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = 0; i2 < 3; ++i2) {
               if (i1 == 2 && i2 == 2) {
                  self_energy_1(i1,i2)  = Re(
                     self_energy_Fu_1_heavy(p,i1,i2));
                  self_energy_PL(i1,i2) = Re(
                     self_energy_Fu_PL_heavy(p,i1,i2));
                  self_energy_PR(i1,i2) = Re(
                     self_energy_Fu_PR_heavy(p,i1,i2));
               } else {
                  self_energy_1(i1,i2)  = Re(
                     self_energy_Fu_1(p,i1,i2));
                  self_energy_PL(i1,i2) = Re(
                     self_energy_Fu_PL(p,i1,i2));
                  self_energy_PR(i1,i2) = Re(
                     self_energy_Fu_PR(p,i1,i2));
               }
            }
         }
         Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
            M_tree - M_tree * self_energy_PL - self_energy_1);
         delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
         const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
         Eigen::Array<double,3,1> eigen_values;
         decltype(ZUL) mix_ZUL;
         decltype(ZUR) mix_ZUR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR,
            eigenvalue_error);
         problems.flag_bad_mass(NMSSMtower_info::Fu,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR);
      #endif
         if (es == 0) {
            PHYSICAL(ZUL) = mix_ZUL;
            PHYSICAL(ZUR) = mix_ZUR;
         }
         PHYSICAL(MFu(es)) = Abs(eigen_values(es));
      }

      new_MFu = PHYSICAL(MFu);
      diff = MaxRelDiff(new_MFu, old_MFu);
      old_MFu = new_MFu;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::Fu);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::Fu);
}

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::VWm))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MVWm) old_MVWm(MVWm), new_MVWm(MVWm);

   do {
      const double M_tree(Sqr(MVWm));
      const double p = old_MVWm;
      const double self_energy = Re(self_energy_VWm(p));
      const double mass_sqr = M_tree - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(NMSSMtower_info::VWm);

      PHYSICAL(MVWm) = AbsSqrt(mass_sqr);

      new_MVWm = PHYSICAL(MVWm);
      diff = MaxRelDiff(new_MVWm, old_MVWm);
      old_MVWm = new_MVWm;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(NMSSMtower_info::VWm);
   else
      problems.unflag_no_pole_mass_convergence(NMSSMtower_info::VWm);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(NMSSMtower_info::VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(NMSSMtower_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(NMSSMtower_info::VZ);

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1 - 0.0023747152416172916*(0.6*Sqr(g1)
      - Sqr(g2));
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   const double qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(idx))
      /Sqr(currentScale)))*Sqr(g3);
   double qcd_2l = 0., qcd_3l = 0.;

   if (get_thresholds() > 1) {
      qcd_2l = -0.003408916029785599*Power(g3,4) -
         0.0011495761378943394*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(currentScale,2)
         /Sqr(MFu(idx))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1 - 0.00020496318737651018*Power(g3,4)
      + 0.0006860288475783287*Sqr(g1) + 0.0023747152416172916*Sqr(g2) -
      0.008443431970194815*Sqr(g3);
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(NMSSMtower_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(NMSSMtower_info::VWm);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::v() const
{
   return Sqrt(Sqr(vd) + Sqr(vu));
}

double CLASSNAME::Betax() const
{
   return ArcSin(Abs(ZP(0,1)));
}

double CLASSNAME::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}


std::ostream& operator<<(std::ostream& ostr, const NMSSMtower_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
