#ifndef slic3r_FillBiholarBroke_hpp_
#define slic3r_FillBiholarBroke_hpp_

#include "../libslic3r.h"

#include "Fill.hpp"

namespace Slic3r {

class FillBiholar : public Fill
{
public:
    virtual Fill* clone() const { return new FillBiholar(*this); };
    virtual ~FillBiholar() {}
    virtual bool can_solid() const { return true; };

protected:
	virtual void _fill_surface_single(
	    unsigned int                     thickness_layers,
	    const direction_t               &direction,
	    ExPolygon                       &expolygon,
	    Polylines*                      polylines_out);

	void _fill_single_direction(ExPolygon expolygon, const direction_t &direction,
	    coord_t x_shift, Polylines* out);

      // Cache the polygon maths.
    	struct CacheData
    	{
            coord_t	distance;

            coord_t w;
            coord_t r1;
            coord_t r2;

            coord_t x_offset;
            coord_t y_offset;
        };
        typedef std::pair<float,coordf_t> CacheID;  // density, spacing
        typedef std::map<CacheID, CacheData> Cache;
    	Cache cache;
};

}; // namespace Slic3r

#endif // slic3r_FillBiholar_hpp_
