#include "FillCustom.hpp"
#include "../ClipperUtils.hpp"
#include "../PolylineCollection.hpp"
#include "../Surface.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


namespace Slic3r {

std::vector<int> parseLine(std::string line){
  std::istringstream iss(line); // e.g. data = "1,2"
  std::vector<int> point;
  int to_draw, x, y;
  char delimiter;
  char open_bracket;
  char close_bracket;
  iss >> open_bracket >> to_draw >> delimiter >> open_bracket >> x >> delimiter >> y >> close_bracket >> close_bracket;
  point.push_back(to_draw);
  point.push_back(x);
  point.push_back(y);
  return point;
}

//reads a file of points and returns a vector of vector points
std::vector<std::vector<int>> fileToPattern(std::string filename){
  std::vector<std::vector<int>> everyline;
  std::string line;
  std::ifstream myfile(filename);
  if(myfile.is_open()){
    while ( getline (myfile,line) )
    {
      everyline.push_back(parseLine(line));
    }
    myfile.close();
  }
  else std::cout << "Unable to open file" << std::endl;
  return everyline;
}

//takes in a vector of vector points and scales them by a scale factor
std::vector<std::vector<int>> scaleFromOrigin(std::vector<std::vector<int>> pattern, float scale){
  std::vector<std::vector<int>> scaled_pattern;
  for(int i = 0; i < pattern.size(); i++){
    scaled_pattern.push_back(std::vector<int>{pattern[i][0], int(pattern[i][1] * scale), int(pattern[i][2] * scale)});
  }
  return scaled_pattern;
}

//takes in a vector of vector points and scales them by a scale factor
std::vector<std::vector<int>> scaleFromOriginX(std::vector<std::vector<int>> pattern, int scale){
  std::vector<std::vector<int>> scaled_pattern;
  for(int i = 0; i < pattern.size(); i++){
    scaled_pattern.push_back(std::vector<int>{pattern[i][0], int(pattern[i][1] * scale), pattern[i][2]});
  }
  return scaled_pattern;
}

//takes in a vector of vector points and scales them by a scale factor
std::vector<std::vector<int>> scaleFromOriginY(std::vector<std::vector<int>> pattern, int scale){
  std::vector<std::vector<int>> scaled_pattern;
  for(int i = 0; i < pattern.size(); i++){
    scaled_pattern.push_back(std::vector<int>{pattern[i][0], pattern[i][1], int(pattern[i][2] * scale)});
  }
  return scaled_pattern;
}

//Adjusts the pattern so the bottom left is zero
std::vector<std::vector<int>> moveCentreToBottomLeft(std::vector<std::vector<int>> pattern){
  std::vector<std::vector<int>> scaled_pattern;
  int xmin = INT_MAX;
  int ymin = INT_MAX;
  for(int i = 0; i < pattern.size(); i++){
    if(pattern[i][1] < xmin) xmin = pattern[i][1];
    if(pattern[i][2] < ymin) ymin = pattern[i][2];
  }

  for(int i = 0; i < pattern.size(); i++){
    scaled_pattern.push_back(std::vector<int>{pattern[i][0], pattern[i][1] - xmin, pattern[i][2] - ymin});
  }
  return scaled_pattern;
}

int calculatePatternWidth(std::vector<std::vector<int>> pattern){
  int min = INT_MAX;
  int max = -INT_MAX;
  for(int i = 0; i < pattern.size(); i++){
    if(pattern[i][1] < min) min = pattern[i][1];
    if(pattern[i][1] > max) max = pattern[i][1];
  }
  return abs(max-min);
}

int calculatePatternHeight(std::vector<std::vector<int>> pattern){
  int min = INT_MAX;
  int max = -INT_MAX;
  for(int i = 0; i < pattern.size(); i++){
    if(pattern[i][2] < min) min = pattern[i][2];
    if(pattern[i][2] > max) max = pattern[i][2];
  }
  return abs(max-min);
}

void
FillCustom::_fill_surface_single(
    unsigned int                    thickness_layers,
    const direction_t               &direction,
    ExPolygon                       &expolygon,
    Polylines*                      polylines_out)
{
    // cache shape
    CacheID cache_id = std::make_pair(this->density, this->min_spacing);
    Cache::iterator it_m = this->cache.find(cache_id);
    if (it_m == this->cache.end()) {
        it_m = this->cache.insert(it_m, std::pair<CacheID,CacheData>(cache_id, CacheData()));
        CacheData &m = it_m->second;
        coord_t min_spacing = scale_(this->min_spacing);
        if(this->meta_isMM){
          m.distance          = 1; //used as scaling
          m.pattern           = moveCentreToBottomLeft(scaleFromOrigin(fileToPattern(this->filename), scale_(1)));
        }
        else{
          m.distance          = (min_spacing / this->density)/10.f; //used as scaling
          m.pattern           = moveCentreToBottomLeft(scaleFromOrigin(fileToPattern(this->filename), m.distance));
        }

        // m.distance          = (min_spacing / this->density); //used as scaling factor

        // caches the entire pattern
        m.pattern_width     = calculatePatternWidth(m.pattern);
        m.pattern_height    = calculatePatternHeight(m.pattern);


        m.x_offset = scale_(this->x_offset);
        m.y_offset = scale_(this->y_offset);

        m.pattern_center    = Point(m.pattern_width/2, m.pattern_height/2);
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
            bb_polygon.rotate(this->infill_angle, m.pattern_center);
            bounding_box = bb_polygon.bounding_box();

            // extend bounding box so that our pattern will be aligned with other layers
            // $bounding_box->[X1] and [Y1] represent the displacement between new bounding box offset and old one
            // The infill is not aligned to the object bounding box, but to a world coordinate system. Supposedly good enough.
            bounding_box.min.align_to_grid(Point(m.pattern_width, m.pattern_height));
        }

        //For x in range
        for (coord_t x = bounding_box.min.x - m.x_offset; x <= bounding_box.max.x; x += m.pattern_width) {
          //For y in range
          for (coord_t y = bounding_box.min.y - m.y_offset; y <= bounding_box.max.y; y += m.pattern_height) {
            bool addpoint = true;
            int i = 0;
            while(i < m.pattern.size())
            {
              Polyline polyline; //create a new polyline
              if (addpoint) {                                                                 //as we've incremented the pointer we need to add this
                polyline.points.push_back(Point(m.pattern[i][1] + x, m.pattern[i][2] + y));   // point if there is a new polyline as we skipped adding
                i++;                                                                          //it at the end of the last loop
              }

              while(i < m.pattern.size() && m.pattern[i][0] != 0)                             //read and add points until there is a "no draw" point
              {                                                                               //or its the end of the file
                polyline.points.push_back(Point(m.pattern[i][1] + x, m.pattern[i][2] + y));
                i++;
              }
              // polyline.rotate(-this->infill_angle, m.pattern_center);
              polylines.push_back(polyline);
              addpoint = true;
            }
          }
        }
    }

    if (true) {
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
