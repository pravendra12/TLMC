#include "LatticeSite.hpp"

namespace std {

    /// Stream operator.
    ostream& operator<<(ostream& os, const LatticeSite& lattice_site)
    {
        os << to_string(lattice_site.index()) + " :";
        for (size_t i = 0; i < 3; i++)
        {
            os << " " + to_string(lattice_site.unitcellOffset()[i]);
        }
        return os;
    }

}
