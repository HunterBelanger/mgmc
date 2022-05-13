#ifndef ZPLANE_H
#define ZPLANE_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class ZPlane : public Surface {
 public:
  ZPlane(double z, BoundaryType bound, uint32_t i_id, std::string i_name);
  ~ZPlane() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double z0;

};  // ZPlane

//===========================================================================
// Non-Member Functions
std::shared_ptr<ZPlane> make_zplane(YAML::Node surface_node);

#endif
