#ifndef YPLANE_H
#define YPLANE_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class YPlane : public Surface {
 public:
  YPlane(double x, BoundaryType bound, uint32_t i_id, std::string i_name);
  ~YPlane() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double y0;

};  // YPlane

//===========================================================================
// Non-Member Functions
std::shared_ptr<YPlane> make_yplane(YAML::Node surface_node);

#endif
