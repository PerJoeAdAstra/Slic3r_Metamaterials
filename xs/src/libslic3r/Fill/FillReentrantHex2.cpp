#include "FillReentrantHex2.hpp" 
#include "../ClipperUtils.hpp"
#include "../PolylineCollection.hpp"
#include "../Surface.hpp"


namespace Slic3r {

void
FillReentrantHex2::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    // cache hexagons math
    CacheID cache_id = std::make_pair(this->density, this->min_spacing);
    Cache::iterator it_m = this->cache.find(cache_id);
    if (it_m == this->cache.end()) {
        it_m = this->cache.insert(it_m, std::pair<CacheID,CacheData>(cache_id, CacheData()));
        CacheData &m = it_m->second;
        coord_t min_spacing = scale_(this->min_spacing);
        if(this->meta_isMM){
          m.distance          = 1; //used as scaling
          m.w                 = scale_(this->meta_l / 2);
          m.h                 = scale_(this->meta_h);
        }
        else
        {
          m.distance          = min_spacing / this->density; //used as scaling factor
          m.w                 = (this->meta_l * m.distance) / 2;
          m.h                 = this->meta_h * m.distance;
        }
        m.theta = this->meta_angle;
        m.y_short           = m.w * tan(m.theta);
        m.pattern_height    = m.h + (m.h - (2 * -m.y_short));
        m.x_offset          = 0;

        m.hex_center        = Point(m.w/2, m.h);

        m.x_pattern_offset          = scale_(this->x_offset);
        m.y_offset          = scale_(this->y_offset);
    }
    CacheData &m = it_m->second;

    Polylines polylines;
    {
        // adjust actual bounding box to the nearest multiple of our hex pattern
        // and align it so that it matches across layers

        BoundingBox bounding_box = expolygon.contour.bounding_box();
        {
            // rotate bounding box according to infill direction
            Polygon bb_polygon = bounding_box.polygon();
            bb_polygon.rotate(this->infill_angle, m.hex_center);
            bounding_box = bb_polygon.bounding_box();

            // extend bounding box so that our pattern will be aligned with other layers
            // $bounding_box->[X1] and [Y1] represent the displacement between new bounding box offset and old one
            // The infill is not aligned to the object bounding box, but to a world coordinate system. Supposedly good enough.
            bounding_box.min.align_to_grid(Point(m.w, m.pattern_height));
        }

        for (coord_t x = bounding_box.min.x - m.x_pattern_offset; x <= bounding_box.max.x; ) {
            coord_t ax[2] = { x + m.x_offset, x + m.w - m.x_offset };
            for (size_t i = 0; i < 2; ++ i) {
                Polyline p;
                Polylines polylines1;
                Polylines polylines2;
                for (coord_t y = bounding_box.min.y - m.y_offset; y <= bounding_box.max.y; y += - m.y_short + m.h - m.y_short + m.h) {
                  if(i == 0){
                    p.points.push_back(Point(ax[1], y));
                    p.points.push_back(Point(ax[0], y - m.y_short));
                    //New polyline
                    p.rotate(-this->infill_angle, m.hex_center);
                    polylines.push_back(p);
                    //
                    Polyline p1;
                    p = p1;

                    p.points.push_back(Point(ax[0], y - m.y_short + m.h));
                    p.points.push_back(Point(ax[1], y - m.y_short + m.h - m.y_short));
                    p.points.push_back(Point(ax[1], y - m.y_short + m.h - m.y_short + m.h));
                  }
                  else{
                    p.points.push_back(Point(ax[1], y));
                    p.points.push_back(Point(ax[0], y - m.y_short));
                    p.points.push_back(Point(ax[0], y - m.y_short + m.h));
                    p.points.push_back(Point(ax[1], y - m.y_short + m.h - m.y_short));
                    //New polyline
                    p.rotate(-this->infill_angle, m.hex_center);
                    polylines.push_back(p);
                    //
                    Polyline p1;
                    p = p1;
                    p.points.push_back(Point(ax[1], y - m.y_short + m.h - m.y_short + m.h));
                  }
                }
                p.rotate(-this->infill_angle, m.hex_center);
                polylines.push_back(p);

                ax[0] = ax[0] + m.w;
                ax[1] = ax[1] + m.w;
                std::swap(ax[0], ax[1]); // draw symmetrical pattern
                x += m.w;
            }
        }
    }

    if (false) {
        printf("Something went wrong?");
    } else {
        // consider polygons as polylines without re-appending the initial point:
        // this cuts the last segment on purpose, so that the jump to the next
        // path is more straight
        Polylines paths = intersection_pl(
            polylines,
            (Polygons)expolygon
        );

        // connect paths
        if (!paths.empty()) { // prevent calling leftmost_point() on empty collections
            Polylines chained = PolylineCollection::chained_path_from(
                STDMOVE(paths),
                PolylineCollection::leftmost_point(paths),
                false
            );
            printf("Reordered paths!\n");
            assert(paths.empty());
            printf("Post assert\n");
            paths.clear();

            for (Polylines::iterator it_path = chained.begin(); it_path != chained.end(); ++ it_path) {
                if (!paths.empty()) {
                    // distance between first point of this path and last point of last path
                    double distance = paths.back().last_point().distance_to(it_path->first_point());
                    if (distance <= m.w) {
                        paths.back().points.insert(paths.back().points.end(), it_path->points.begin(), it_path->points.end());
                        continue;
                    }
                }
                // Don't connect the paths.
                paths.push_back(*it_path);
            }
        }

        //clip paths again to prevent connection segments from crossing the expolygon boundaries
        paths = intersection_pl(paths, Polygons(expolygon));

        // Move the polylines to the output, avoid a deep copy.
        size_t j = polylines_out->size();
        polylines_out->resize(j + paths.size(), Polyline());
        for (size_t i = 0; i < paths.size(); ++ i)
            std::swap((*polylines_out)[j++], paths[i]);
    }
}

} // namespace Slic3r
