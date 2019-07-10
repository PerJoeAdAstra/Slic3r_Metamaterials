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


        m.starOffset = this->meta_l;
        m.starHeight = this->meta_h;

        m.in_short = this->meta_l * sin(Geometry::deg2rad(30));
        m.in_long = this->meta_l * cos(Geometry::deg2rad(30));

        m.out_long = this->meta_h * cos(Geometry::deg2rad(30));
        m.out_short = this->meta_h * sin(Geometry::deg2rad(30));

    }
    CacheData &m = it_m->second;

    Polygons polygons;
    {
        // adjust actual bounding box to the nearest multiple of our hex pattern
        // and align it so that it matches across layers

        BoundingBox bounding_box = expolygon.contour.bounding_box();
        {
            // rotate bounding box according to infill direction
            Polygon bb_polygon = bounding_box.polygon();
            bb_polygon.rotate(0, m.hex_center);
            bounding_box = bb_polygon.bounding_box();

            // extend bounding box so that our pattern will be aligned with other layers
            // $bounding_box->[X1] and [Y1] represent the displacement between new bounding box offset and old one
            // The infill is not aligned to the object bounding box, but to a world coordinate system. Supposedly good enough.
            bounding_box.min.align_to_grid(Point(m.hex_width, m.pattern_height));
        }

        for (coord_t x = bounding_box.min.x; x <= bounding_box.max.x; ) {
            coord_t ax[2] = { x + m.x_offset, x + m.distance - m.x_offset };
            // for (size_t i = 0; i < 2; ++ i) {
                // std::reverse(p.points.begin(), p.points.end()); // turn first half upside down
                for (coord_t y = bounding_box.min.y; y <= bounding_box.max.y; y += m.y_short + m.hex_side + m.y_short + m.hex_side) {

                    // Polyline polyline;
                    // polyline.points.reserve(this->points.size() + 1);
                    // for (Points::const_iterator it = this->points.begin() + index; it != this->points.end(); ++it)
                    //     polyline.points.push_back(*it);
                    // for (Points::const_iterator it = this->points.begin(); it != this->points.begin() + index + 1; ++it)
                    //     polyline.points.push_back(*it);
                    // return polyline;

                    // printf("loop!");
                    Polygon p1;
                    p1.points.push_back(Point(ax[0], y + m.y_short + m.hex_side + m.y_offset));
                    // p.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short - m.y_offset));
                    p1.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short - m.y_offset));

                    // p.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short + m.hex_side + m.y_offset));
                    p1.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short + m.hex_side + m.y_offset));
                    p1.rotate(0, m.hex_center);
                    polygons.push_back(p1);

                    Polygon p;

                    //Points p;
                    // Point point = Point(ax[1], y + m.y_offset)
                    p.points.push_back(Point(ax[1], y + m.y_offset));

                    p.points.push_back(Point(ax[0], y + m.y_short - m.y_offset));

                    p.points.push_back(Point(ax[0], y + m.y_short + m.hex_side + m.y_offset));

                    p.rotate(0, m.hex_center);
                    polygons.push_back(p); //<- a polygon of points


                    // Point point = Point(ax[1], y);
                    // //p.points.push_back(Point(ax[1], y + m.y_offset));
                    // p.points.push_back((point.x, point.y - m.starOffset));           //#1
                    // p.points.push_back((point.x - m.out_long, point.y - m.out_short));   //#2
                    // p.points.push_back((point.x - m.in_long, point.y + m.in_short));     //#3
                    // p.points.push_back((point.x, point.y + m.starHeight));           //#4
                    // p.rotate(0, m.hex_center);
                    // polygons.push_back(p);
                    // Polygon p1;                                                       //gap?
                    // p1.points.push_back((point.x - m.in_long, point.y + m.in_short));     //#5
                    //
                    // point = Point(ax[0], y);
                    // //p.points.push_back(Point(ax[0], y + m.y_short - m.y_offset));
                    // p1.points.push_back((point.x + m.in_long, point.y - m.in_short));
                    // p1.points.push_back((point.x, point.y - m.starHeight));
                    // p1.rotate(0, m.hex_center);
                    // polygons.push_back(p1);
                    // Polygon p2;
                    // p2.points.push_back((point.x + m.in_long, point.y - m.in_short));
                    // p2.points.push_back((point.x + m.out_long, point.y + m.out_short));
                    // p2.points.push_back((point.x, point.y + m.starOffset));
                    //
                    // point = Point(ax[0], y + m.y_short + m.hex_side + m.y_offset);
                    // //p.points.push_back(Point(ax[0], y + m.y_short + m.hex_side + m.y_offset));
                    // p2.points.push_back((point.x, point.y - m.starOffset));
                    // p2.points.push_back((point.x + m.out_long, point.y - m.out_short));
                    // p2.points.push_back((point.x + m.in_long, point.y + m.in_short));
                    // p2.points.push_back((point.x, point.y + m.starHeight));
                    // p2.rotate(0, m.hex_center);
                    // polygons.push_back(p2);
                    // Polygon p3;
                    // p3.points.push_back((point.x + m.in_long, point.y + m.in_short));
                    //
                    // point = Point(ax[1], y + m.y_short + m.hex_side + m.y_short - m.y_offset);
                    // //p.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short - m.y_offset));
                    // p3.points.push_back((point.x - m.in_long, point.y - m.in_short));
                    // p3.points.push_back((point.x, point.y - m.starHeight));
                    // p3.rotate(0, m.hex_center);
                    // polygons.push_back(p3);
                    // Polygon p4;
                    // p4.points.push_back((point.x - m.in_long, point.y - m.in_short));
                    // p4.points.push_back((point.x - m.out_long, point.y + m.out_short));
                    // p4.points.push_back((point.x, point.y + m.in_short));
                    //

                    // //p.push_back(Point[ax[1], y + m.y_short + m.hex_side + m.y_short + (m.hex_side)])
                    // //p.rotate(0, m.hex_center);
                    // p4.points.push_back(Point(ax[1], y + m.y_short + m.hex_side + m.y_short + m.hex_side + m.y_offset - m.in_short));
                    // p4.rotate(0, m.hex_center);

                    // //polygons.push_back(p); //<- a polygon of points
                    // polygons.push_back(p);
                }
                // ax[0] = ax[0] + m.distance;
                // ax[1] = ax[1] + m.distance;
                // std::swap(ax[0], ax[1]); // draw symmetrical pattern
                x += m.distance;
            // }
        }
    }

    if (true || this->complete) {
        // we were requested to complete each loop;
        // in this case we don't try to make more continuous paths
        Polygons polygons_trimmed = intersection((Polygons)expolygon, polygons);                    //clipping with the exterior
        for (Polygons::iterator it = polygons_trimmed.begin(); it != polygons_trimmed.end(); ++ it) //Adding trimmed polygons to polylines out
            // polylines_out->push_back(it->split_at_first_point());
            polylines_out->push_back(it->split_at_first_point_no_loop());

    } else {
        printf("How?");
        // consider polygons as polylines without re-appending the initial point:
        // this cuts the last segment on purpose, so that the jump to the next
        // path is more straight
        Polylines paths = intersection_pl(
            to_polylines(polygons),
            (Polygons)expolygon
        );

        // connect paths
        if (!paths.empty()) { // prevent calling leftmost_point() on empty collections
            Polylines chained = PolylineCollection::chained_path_from(
                STDMOVE(paths),
                PolylineCollection::leftmost_point(paths),
                false
            );
            assert(paths.empty());
            paths.clear();

            for (Polylines::iterator it_path = chained.begin(); it_path != chained.end(); ++ it_path) {
                if (!paths.empty()) {
                    // distance between first point of this path and last point of last path
                    double distance = paths.back().last_point().distance_to(it_path->first_point());
                    if (distance <= m.hex_width) {
                        paths.back().points.insert(paths.back().points.end(), it_path->points.begin(), it_path->points.end());
                        continue;
                    }
                }
                // Don't connect the paths.
                paths.push_back(*it_path);
            }
        }

        // clip paths again to prevent connection segments from crossing the expolygon boundaries
        paths = intersection_pl(paths, to_polygons(offset_ex(expolygon, SCALED_EPSILON)));

        // Move the polylines to the output, avoid a deep copy.
        size_t j = polylines_out->size();
        polylines_out->resize(j + paths.size(), Polyline());
        for (size_t i = 0; i < paths.size(); ++ i)
            std::swap((*polylines_out)[j++], paths[i]);
    }
}

} // namespace Slic3r
