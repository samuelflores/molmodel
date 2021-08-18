#ifndef _MOLMODEL_EXCEPTIONS_H
#define _MOLMODEL_EXCEPTIONS_H

#include "Linkage.h"
#include <stdexcept>

namespace SimTK {

class SimTK_MOLMODEL_EXPORT UnrecoverableMolmodelError : public std::runtime_error {
    public:
        using std::runtime_error::runtime_error;
	using std::runtime_error::what;
};

} // namespace SimTK

#endif // _MOLMODEL_EXCEPTIONS_H
