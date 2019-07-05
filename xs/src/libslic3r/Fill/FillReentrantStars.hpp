#ifndef slic3r_FillReentrantStars_hpp_
#define slic3r_FillReentrantStars_hpp_

#include <map>

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillReentrantStars : public Fill
{
public:
    virtual ~FillReentrantStars() {}

protected:
    virtual Fill* clone() const { return new FillReentrantStars(*this); };
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction,
	    ExPolygon                       &expolygon,
	    Polylines*                      polylines_out
	);

	// Cache the hexagon math.
	struct CacheData
	{
		coord_t	distance;
        coord_t hex_side;
        coord_t hex_width;
        coord_t	pattern_height;
        coord_t y_short;
        coord_t x_offset;
        coord_t	y_offset;
        Point	hex_center;
    };
    typedef std::pair<float,coordf_t> CacheID;  // density, spacing
    typedef std::map<CacheID, CacheData> Cache;
	Cache cache;

    virtual float _layer_angle(size_t idx) const { return float(M_PI/3.) * (idx % 3); }
};

} // namespace Slic3r

#endif // slic3r_FillReentrantStars_hpp_
