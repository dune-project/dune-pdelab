#include "config.h"

#include <dune/pdelab/common/typetree/utility.hh>

namespace Dune {
  namespace PDELab {
    namespace TypeTree {

      namespace {
        static const shared_ptr<EmptyNode> _emptyNodePtr(make_shared<EmptyNode>());
      }

      const shared_ptr<EmptyNode>& emptyNodePtr()
      {
        return _emptyNodePtr;
      }

    } // namespace TypeTree
  } // namespace PDELab
} // namespace Dune
