#ifndef slic3r_FillCustom_hpp_
#define slic3r_FillCustom_hpp_

#include <map>

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillCustom : public Fill
{
public:
    virtual ~FillCustom() {}

protected:
    virtual Fill* clone() const { return new FillCustom(*this); };
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction,
	    ExPolygon                       &expolygon,
	    Polylines*                      polylines_out
	);

	// Cache the pattern maths.
	struct CacheData
	{
        coord_t	distance;
        std::vector<std::vector<int>>  pattern;
        coord_t pattern_width;
        coord_t	pattern_height;
        coord_t x_offset;
        coord_t	y_offset;
        Point	hex_center;
    };
    typedef std::pair<float,coordf_t> CacheID;  // density, spacing
    typedef std::map<CacheID, CacheData> Cache;
	Cache cache;
};

} // namespace Slic3r

#endif // slic3r_FillCustom_hpp_
