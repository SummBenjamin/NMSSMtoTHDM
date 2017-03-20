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

// File generated at Mon 20 Mar 2017 08:24:27

#ifndef NMSSMtower_TWO_SCALE_CONVERGENCE_TESTER_H
#define NMSSMtower_TWO_SCALE_CONVERGENCE_TESTER_H

#include "NMSSMtower_convergence_tester.hpp"
#include "NMSSMtower_two_scale_model.hpp"
#include "two_scale_convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class NMSSMtower_convergence_tester<Two_scale> : public Convergence_tester_DRbar<NMSSMtower<Two_scale> > {
public:
   NMSSMtower_convergence_tester(NMSSMtower<Two_scale>*, double);
   virtual ~NMSSMtower_convergence_tester();

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
