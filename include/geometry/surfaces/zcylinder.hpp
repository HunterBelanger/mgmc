#ifndef ZCYLINDER_H
#define ZCYLINDER_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class ZCylinder : public Surface {
 public:
  ZCylinder(double x_, double y_, double r_, BoundaryType bound, uint32_t i_id,
            std::string i_name);
  ~ZCylinder() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double x0, y0, R;

};  // ZCylinder

//===========================================================================
// Non-Member Functions
std::shared_ptr<ZCylinder> make_zcylinder(YAML::Node surface_node);

#endif
