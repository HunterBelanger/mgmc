#ifndef XPLANE_H
#define XPLANE_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class XPlane : public Surface {
 public:
  XPlane(double x, BoundaryType bound, uint32_t i_id, std::string i_name);
  ~XPlane() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double x0;

};  // XPlane

//===========================================================================
// Non-Member Functions
std::shared_ptr<XPlane> make_xplane(YAML::Node surface_node);

#endif
