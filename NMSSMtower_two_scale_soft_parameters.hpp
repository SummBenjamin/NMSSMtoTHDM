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

// File generated at Mon 20 Mar 2017 08:23:03

#ifndef NMSSMtower_TWO_SCALE_soft_parameters_H
#define NMSSMtower_TWO_SCALE_soft_parameters_H

#include "NMSSMtower_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class NMSSMtower_soft_parameters : public NMSSMtower_susy_parameters {
public:
   explicit NMSSMtower_soft_parameters(const NMSSMtower_input_parameters& input_ = NMSSMtower_input_parameters());
   NMSSMtower_soft_parameters(const NMSSMtower_susy_parameters& , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, double TLambdax_, double TKappa_, const Eigen::Matrix<double,3,3>&
   TYu_, const Eigen::Matrix<double,3,3>& mq2_, const Eigen::Matrix<double,3,3>
   & ml2_, double mHd2_, double mHu2_, const Eigen::Matrix<double,3,3>& md2_,
   const Eigen::Matrix<double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_
   , double ms2_, double MassB_, double MassWB_, double MassG_
);
   virtual ~NMSSMtower_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   NMSSMtower_soft_parameters calc_beta() const;
   NMSSMtower_soft_parameters calc_beta(unsigned) const;
   virtual void clear();

   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) { TYd = TYd_; }
   void set_TYd(int i, int k, double value) { TYd(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) { TYe = TYe_; }
   void set_TYe(int i, int k, double value) { TYe(i,k) = value; }
   void set_TLambdax(double TLambdax_) { TLambdax = TLambdax_; }
   void set_TKappa(double TKappa_) { TKappa = TKappa_; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) { TYu = TYu_; }
   void set_TYu(int i, int k, double value) { TYu(i,k) = value; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) { mq2 = mq2_; }
   void set_mq2(int i, int k, double value) { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) { ml2 = ml2_; }
   void set_ml2(int i, int k, double value) { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) { md2 = md2_; }
   void set_md2(int i, int k, double value) { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) { mu2 = mu2_; }
   void set_mu2(int i, int k, double value) { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) { me2 = me2_; }
   void set_me2(int i, int k, double value) { me2(i,k) = value; }
   void set_ms2(double ms2_) { ms2 = ms2_; }
   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }

   const Eigen::Matrix<double,3,3>& get_TYd() const { return TYd; }
   double get_TYd(int i, int k) const { return TYd(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYe() const { return TYe; }
   double get_TYe(int i, int k) const { return TYe(i,k); }
   double get_TLambdax() const { return TLambdax; }
   double get_TKappa() const { return TKappa; }
   const Eigen::Matrix<double,3,3>& get_TYu() const { return TYu; }
   double get_TYu(int i, int k) const { return TYu(i,k); }
   const Eigen::Matrix<double,3,3>& get_mq2() const { return mq2; }
   double get_mq2(int i, int k) const { return mq2(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2() const { return ml2; }
   double get_ml2(int i, int k) const { return ml2(i,k); }
   double get_mHd2() const { return mHd2; }
   double get_mHu2() const { return mHu2; }
   const Eigen::Matrix<double,3,3>& get_md2() const { return md2; }
   double get_md2(int i, int k) const { return md2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2() const { return mu2; }
   double get_mu2(int i, int k) const { return mu2(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2() const { return me2; }
   double get_me2(int i, int k) const { return me2(i,k); }
   double get_ms2() const { return ms2; }
   double get_MassB() const { return MassB; }
   double get_MassWB() const { return MassWB; }
   double get_MassG() const { return MassG; }


protected:
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,3> TYe;
   double TLambdax;
   double TKappa;
   Eigen::Matrix<double,3,3> TYu;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double ms2;
   double MassB;
   double MassWB;
   double MassG;


private:
   static const int numberOfParameters = 115;

   struct Soft_traces {
      double traceAdjYdTYd;
      double traceAdjYeTYe;
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceAdjYuTYu;
      double traceYdAdjYdTYdAdjYd;
      double traceYdAdjYuTYuAdjYd;
      double traceYeAdjYeTYeAdjYe;
      double traceYuAdjYdTYdAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;
      double traceYuAdjYuTYuAdjYu;
      double traceconjTYdTpTYd;
      double traceconjTYeTpTYe;
      double tracemd2YdAdjYd;
      double traceme2YeAdjYe;
      double traceml2AdjYeYe;
      double tracemq2AdjYdYd;
      double traceconjTYdTpYd;
      double traceconjTYeTpYe;
      double traceconjTYuTpTYu;
      double tracemq2AdjYuYu;
      double tracemu2YuAdjYu;
      double traceconjTYuTpYu;
      double traceYdAdjYdTYdAdjTYd;
      double traceYdAdjYuTYuAdjTYd;
      double traceYdAdjTYdTYdAdjYd;
      double traceYdAdjTYuTYuAdjYd;
      double traceYeAdjYeTYeAdjTYe;
      double traceYeAdjTYeTYeAdjYe;
      double traceYuAdjYdTYdAdjTYu;
      double traceYuAdjTYdTYdAdjYu;
      double tracemd2YdAdjYdYdAdjYd;
      double tracemd2YdAdjYuYuAdjYd;
      double traceme2YeAdjYeYeAdjYe;
      double traceml2AdjYeYeAdjYeYe;
      double tracemq2AdjYdYdAdjYdYd;
      double tracemq2AdjYdYdAdjYuYu;
      double tracemq2AdjYuYuAdjYdYd;
      double tracemu2YuAdjYdYdAdjYu;
      double traceYuAdjYuTYuAdjTYu;
      double traceYuAdjTYuTYuAdjYu;
      double tracemq2AdjYuYuAdjYuYu;
      double tracemu2YuAdjYuYuAdjYu;
      double Tr11;
      double Tr2U111;
      double Tr31;
      double Tr22;
      double Tr23;

   };
   void calc_soft_traces(Soft_traces&) const;

   Eigen::Matrix<double,3,3> calc_beta_TYd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappa_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappa_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappa_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_three_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const NMSSMtower_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
