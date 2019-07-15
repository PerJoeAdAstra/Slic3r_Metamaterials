#include "FillReentrantStars.hpp"
#include "../ClipperUtils.hpp"
#include "../PolylineCollection.hpp"
#include "../Surface.hpp"


namespace Slic3r {

void
FillReentrantStars::_fill_surface_single(
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
        m.distance          = min_spacing / this->density;
        m.hex_side          = m.distance / (sqrt(3)/2);
        m.hex_width         = m.distance * 2; // $m->{hex_width} == $m->{hex_side} * sqrt(3);
        coord_t hex_height  = m.hex_side * 2;
        m.pattern_height    = hex_height + m.hex_side;
        m.y_short           = m.distance * sqrt(3)/3;
        m.x_offset          = 0;
        m.y_offset          = 0;
        m.hex_center        = Point(m.hex_width/2, m.hex_side);


        m.starOffset = m.distance * this->meta_l/8;
        m.starHeight = m.distance * this->meta_h/2;

        m.in_short = m.starOffset * sin(Geometry::deg2rad(30));
        m.in_long = m.starOffset * cos(Geometry::deg2rad(30));

        m.out_long = m.starHeight * cos(Geometry::deg2rad(30));
        m.out_short = m.starHeight * sin(Geometry::deg2rad(30));
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
            bb_polygon.rotate(0, m.hex_center); //bb_polygon.rotate(this->infill_angle, m.hex_center);
            bounding_box = bb_polygon.bounding_box();

            // extend bounding box so that our pattern will be aligned with other layers
            // $bounding_box->[X1] and [Y1] represent the displacement between new bounding box offset and old one
            // The infill is not aligned to the object bounding box, but to a world coordinate system. Supposedly good enough.
            bounding_box.min.align_to_grid(Point(m.hex_width, m.pattern_height));
        }
        for (coord_t x = bounding_box.min.x; x <= bounding_box.max.x; ) {
            coord_t ax[2] = { x + m.x_offset, x + m.distance - m.x_offset };
            for (coord_t i = 1; i > -2; i -= 2) {
                for (coord_t y = bounding_box.min.y; y <= bounding_box.max.y; y += m.y_short + m.hex_side + m.y_short + m.hex_side) {
                    //first star
                    Polyline p;
                    p.points.push_back(Point(ax[1], y - m.starOffset));
                    p.points.push_back(Point(ax[1] - (i * m.out_long), y - m.out_short));
                    p.points.push_back(Point(ax[1] - (i * m.in_long), y + m.in_short));
                    p.points.push_back(Point(ax[1], y + m.starHeight));
                    if(i==-1){
                      std::reverse(p.points.begin(), p.points.end());
                    }
                    polylines.push_back(p);

                    Polyline p1;
                    p1.points.push_back(Point(ax[1] - (i * m.in_long), y + m.in_short));
                    //second star
                    p1.points.push_back(Point(ax[0] + (i * m.in_long), y + m.y_short - m.in_short));
                    p1.points.push_back(Point(ax[0], y + m.y_short - m.starHeight));
                    if(i == -1){
                      std::reverse(p1.points.begin(), p1.points.end());
                    }
                    polylines.push_back(p1);

                    Polyline p2;
                    p2.points.push_back(Point(ax[0] + (i * m.in_long), y + m.y_short - m.in_short));
                    p2.points.push_back(Point(ax[0] + (i * m.out_long), y + m.y_short + m.out_short));
                    p2.points.push_back(Point(ax[0], y + m.y_short + m.starOffset));
                    //Third star
                    p2.points.push_back(Point(ax[0], y + m.y_short + m.hex_side - m.starOffset ));
                    p2.points.push_back(Point(ax[0] + (i * m.out_long), y + m.y_short + m.hex_side - m.out_short));
                    p2.points.push_back(Point(ax[0] + (i * m.in_long), y + m.y_short + m.hex_side + m.in_short));
                    p2.points.push_back(Point(ax[0], y + m.y_short + m.hex_side + m.starHeight));
                    if(i == -1){
                      std::reverse(p2.points.begin(), p2.points.end());
                    }
                    polylines.push_back(p2);

                    Polyline p3;
                    p3.points.push_back(Point(ax[0] + (i * m.in_long), y + m.y_short + m.hex_side + m.in_short));
                    p3.points.push_back(Point(ax[1] - (i * m.in_long), y + m.y_short + m.hex_side + m.y_short - m.in_long));
                    //Fourth star
                    p3.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short - m.starHeight));
                    if(i == -1){
                      std::reverse(p3.points.begin(), p3.points.end());
                    }
                    polylines.push_back(p3);

                    Polyline p4;
                    p4.points.push_back(Point(ax[1] - (i * m.in_long), y + m.y_short + m.hex_side + m.y_short - m.in_long));
                    p4.points.push_back(Point(ax[1] - (i * m.out_long), y + m.y_short + m.hex_side + m.y_short + m.out_short));
                    p4.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short + m.starOffset));
                    p4.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short + m.hex_side - m.starOffset));
                    if(i == -1){
                      std::reverse(p4.points.begin(), p4.points.end());
                    }
                    polylines.push_back(p4);
                }
                ax[0] = ax[0] + m.distance;
                ax[1] = ax[1] + m.distance;
                std::swap(ax[0], ax[1]); // draw symmetrical pattern
                x += m.distance;
            }
        }
    }

    if (false) {
        // we were requested to complete each loop;
        // in this case we don't try to make more continuous paths
        printf("Something went wrong!?\n");
    } else {
        // consider polygons as polylines without re-appending the initial point:
        // this cuts the last segment on purpose, so that the jump to the next
        // path is more straight
        Polylines paths = intersection_pl(
            polylines,
            (Polygons)expolygon
        );

        // clip paths again to prevent connection segments from crossing the expolygon boundaries
        // paths = intersection_pl(paths, to_polygons(offset_ex(expolygon, SCALED_EPSILON)));

        // Move the polylines to the output, avoid a deep copy.
        size_t j = polylines_out->size();
        polylines_out->resize(j + paths.size(), Polyline());
        for (size_t i = 0; i < paths.size(); ++ i)
            std::swap((*polylines_out)[j++], paths[i]);
    }
}

} // namespace Slic3r
