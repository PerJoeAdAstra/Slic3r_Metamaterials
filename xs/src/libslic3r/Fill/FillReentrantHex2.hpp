#ifndef slic3r_FillReentrantHex2_hpp_
#define slic3r_FillReentrantHex2_hpp_

#include <map>

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillReentrantHex2 : public Fill
{
public:
    virtual ~FillReentrantHex2() {}

protected:
    virtual Fill* clone() const { return new FillReentrantHex2(*this); };
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction,
	    ExPolygon                       &expolygon,
	    Polylines*                      polylines_out
	);

	// Cache the hexagon math.
	struct CacheData
	{
        coord_t h;
        coord_t w;
        coordf_t theta;

        coord_t	distance;
        coord_t	pattern_height;
        coord_t y_short;
        coord_t x_offset;
        Point	hex_center;

        coord_t x_pattern_offset;
        coord_t y_offset;
    };
    typedef std::pair<float,coordf_t> CacheID;  // density, spacing
    typedef std::map<CacheID, CacheData> Cache;
	Cache cache;
};

} // namespace Slic3r

#endif // slic3r_FillReentrantHex_hpp_
