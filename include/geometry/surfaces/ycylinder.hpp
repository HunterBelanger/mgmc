#ifndef YCYLINDER_H
#define YCYLINDER_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class YCylinder : public Surface {
 public:
  YCylinder(double x_, double z_, double r_, BoundaryType bound, uint32_t i_id,
            std::string i_name);
  ~YCylinder() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double x0, z0, R;

};  // YCylinder

//===========================================================================
// Non-Member Functions
std::shared_ptr<YCylinder> make_ycylinder(YAML::Node surface_node);

#endif
