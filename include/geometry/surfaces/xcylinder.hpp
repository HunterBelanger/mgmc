#ifndef XCYLINDER_H
#define XCYLINDER_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class XCylinder : public Surface {
 public:
  XCylinder(double y_, double z_, double r_, BoundaryType bound, uint32_t i_id,
            std::string i_name);
  ~XCylinder() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double y0, z0, R;

};  // XCylinder

//===========================================================================
// Non-Member Functions
std::shared_ptr<XCylinder> make_xcylinder(YAML::Node surface_node);

#endif
