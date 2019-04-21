#ifndef slic3r_SupportMaterial_hpp_
#define slic3r_SupportMaterial_hpp_

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include "libslic3r.h"

#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Fill/Fill.hpp"
#include "Flow.hpp"
#include "Geometry.hpp"
#include "Layer.hpp"
#include "Polygon.hpp"
#include "Print.hpp"
#include "PrintConfig.hpp"
#include "SVG.hpp"

using namespace std;

namespace Slic3r
{

// how much we extend support around the actual contact area
constexpr coordf_t SUPPORT_MATERIAL_MARGIN = 1.5;

constexpr coordf_t MARGIN_STEP = SUPPORT_MATERIAL_MARGIN / 3;

constexpr coordf_t PILLAR_SIZE = 2.5;

constexpr coordf_t PILLAR_SPACING = 10;

/// Struct for carrying the toolpaths parameters needed for each thread.
struct toolpaths_params
{
    int contact_loops;
    coordf_t circle_radius;
    coordf_t circle_distance;
    Polygon circle;
    SupportMaterialPattern pattern;
    vector<int> angles;
    double interface_angle{};
    double interface_spacing{};
    float interface_density{};
    double support_spacing{};
    double support_density{};

    toolpaths_params(int contact_loops = 0,
                     coordf_t circle_radius = 0,
                     coordf_t circle_distance = 0,
                     const Polygon &circle = Polygon(),
                     const SupportMaterialPattern &pattern = SupportMaterialPattern(),
                     const vector<int> &angles = vector<int>())
        : contact_loops(contact_loops),
          circle_radius(circle_radius),
          circle_distance(circle_distance),
          circle(circle),
          pattern(pattern),
          angles(angles)
    {}
};

class SupportMaterial
{
public:
    friend PrintObject;

    PrintConfig *config; ///< The print config
    PrintObjectConfig *object_config; ///< The object print config.
    Flow flow; ///< The intermediate layers print flow.
    Flow first_layer_flow; ///< The first (base) layers print flow.
    Flow interface_flow; ///< The interface layers print flow.

    /// Generate the extrusions paths for the support matterial generated for the given print object.
    void generate_toolpaths(PrintObject *object,
                            map<coordf_t, Polygons> overhang,
                            map<coordf_t, Polygons> contact,
                            map<int, Polygons> _interface,
                            map<int, Polygons> base);

    /// Generate support material for the given print object.
    void generate(PrintObject *object);

    /// Generate the support layers slicing z coordinates.
    vector<coordf_t> support_layers_z(vector<coordf_t> contact_z,
                                      vector<coordf_t> top_z,
                                      coordf_t max_object_layer_height);

    pair<map<coordf_t, Polygons>, map<coordf_t, Polygons>> contact_area(PrintObject *object);

    map<coordf_t, Polygons> object_top(PrintObject *object, map<coordf_t, Polygons> *contact);

    void generate_pillars_shape(const map<coordf_t, Polygons> &contact,
                                const vector<coordf_t> &support_z,
                                map<int, Polygons> &shape);

    map<int, Polygons> generate_base_layers(vector<coordf_t> support_z,
                                            map<coordf_t, Polygons> contact,
                                            map<int, Polygons> _interface,
                                            map<coordf_t, Polygons> top);

    map<int, Polygons> generate_interface_layers(vector<coordf_t> support_z,
                                                 map<coordf_t, Polygons> contact,
                                                 map<coordf_t, Polygons> top);

    void generate_bottom_interface_layers(const vector<coordf_t> &support_z,
                                          map<int, Polygons> &base,
                                          map<coordf_t, Polygons> &top,
                                          map<int, Polygons> &_interface);

    coordf_t contact_distance(coordf_t layer_height, coordf_t nozzle_diameter);

    /// This method returns the indices of the layers overlapping with the given one.
    vector<int> overlapping_layers(int layer_idx, const vector<coordf_t> &support_z);

    void clip_with_shape(map<int, Polygons> &support, map<int, Polygons> &shape);

    // This method removes object silhouette from support material
    // (it's used with interface and base only). It removes a bit more,
    // leaving a thin gap between object and support in the XY plane.
    void clip_with_object(map<int, Polygons> &support, vector<coordf_t> support_z, PrintObject &object);

    void process_layer(int layer_id,
                       toolpaths_params params,
                       const map<coord_t, Polygons>& in_overhang,
                       const map<coord_t, Polygons>& in_contact,
                       const map<int, Polygons>& in_interface,
                       const map<int, Polygons>& in_base);

private:
    /// SupportMaterial is generated by PrintObject.
    SupportMaterial(PrintConfig *print_config,
                    PrintObjectConfig *print_object_config,
                    Flow flow,
                    Flow first_layer_flow,
                    Flow interface_flow)
        : config(print_config),
          object_config(print_object_config),
          flow(flow),
          first_layer_flow(first_layer_flow),
          interface_flow(interface_flow),
          object(nullptr)
    {}

    // Get the maximum layer height given a print object.
    coordf_t get_max_layer_height(PrintObject *object);

    // (Deprecated) use append_to instead
    void append_polygons(Polygons &dst, Polygons &src);

    // Return polygon vector given a vector of surfaces.
    Polygons p(SurfacesPtr &surfaces);

    vector<coordf_t> get_keys_sorted(map<coordf_t, Polygons> _map);

    Polygon create_circle(coordf_t radius);

    // Used during generate_toolpaths function.
    PrintObject *object;
    map<coordf_t, Polygons> overhang;
    map<coordf_t, Polygons> contact;
    map<int, Polygons> _interface;
    map<int, Polygons> base;

};

}

#endif

