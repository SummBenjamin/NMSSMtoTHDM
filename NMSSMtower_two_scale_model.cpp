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
 * @file NMSSMtower_two_scale_model.cpp
 * @brief implementation of the NMSSMtower model class
 *
 * Contains the definition of the NMSSMtower model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 20 Mar 2017 08:52:43 with FlexibleSUSY
 * 1.7.3 (git commit: 4b72fcc4bd7cf8f9d865a5278b43566468f1243d) and SARAH 4.10.0 .
 */

#include "NMSSMtower_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME NMSSMtower<Two_scale>

CLASSNAME::NMSSMtower(const NMSSMtower_input_parameters& input_)
   : Two_scale_model()
   , NMSSMtower_mass_eigenstates(input_)
{
}

CLASSNAME::~NMSSMtower()
{
}

void CLASSNAME::calculate_spectrum()
{
   NMSSMtower_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   NMSSMtower_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return NMSSMtower_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   NMSSMtower_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   NMSSMtower_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   NMSSMtower_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const NMSSMtower<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
