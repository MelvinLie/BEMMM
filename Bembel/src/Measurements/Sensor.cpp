// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SENSOR_H_
#define BEMBEL_SENSOR_H_

#include <Eigen/Dense>

#include "../LinearOperator/DifferentialFormEnum.hpp"
#include "../LinearOperator/LinearOperatorTraits.hpp"
#include "../util/Macros.hpp"

namespace Bembel {
/**
 *    \ingroup Sensor
 *    \brief struct containing specifications on the functional
 *           has to be specialized or derived for any particular operator
 *           under consideration
 **/
template <typename Derived>
struct SensorTraits {
  enum { YOU_DID_NOT_SPECIFY_POTENTIAL_TRAITS = 1 };
};



}  // namespace Bembel
#endif
