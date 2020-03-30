#ifndef slic3r_FillArrowhead_hpp_
#define slic3r_FillArrowhead_hpp_

#include <map>

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillArrowhead : public Fill
{
public:
    virtual ~FillArrowhead() {}

protected:
    virtual Fill* clone() const { return new FillArrowhead(*this); };
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
        coord_t hex_side;
        coord_t hex_width;
        coord_t	pattern_height;
        coord_t y_short;
        Point	hex_center;
    };
    typedef std::pair<float,coordf_t> CacheID;  // density, spacing
    typedef std::map<CacheID, CacheData> Cache;
	Cache cache;
};

} // namespace Slic3r

#endif // slic3r_FillArrowhead_hpp_
