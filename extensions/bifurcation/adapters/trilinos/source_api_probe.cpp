// Compile-time source/API probe for the pinned Trilinos tree. Keeping this in
// the adapter makes version drift fail at the dependency boundary instead of
// leaking NOX types into the mathematical provider contract.
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <LOCA_AnasaziOperator_AbstractStrategy.H>
#include <LOCA_MultiContinuation_AbstractGroup.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Tpetra_Vector.hpp>
#include <BelosLinearProblem.hpp>
#include <Teuchos_ParameterList.hpp>
#include <type_traits>

namespace bif::trilinos {
int sourceApiProbe() {
  Teuchos::ParameterList parameters;
  parameters.set("Linear Solver", "Belos GMRES");
  static_assert(std::is_polymorphic<NOX::Abstract::Group>::value, "NOX group API changed");
  static_assert(std::is_polymorphic<LOCA::MultiContinuation::AbstractGroup>::value, "LOCA group API changed");
  return parameters.isParameter("Linear Solver") ? 0 : 1;
}
} // namespace bif::trilinos
