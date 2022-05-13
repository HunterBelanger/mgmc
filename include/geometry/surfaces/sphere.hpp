#ifndef SPHERE_H
#define SPHERE_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class Sphere : public Surface {
 public:
  Sphere(double x_, double y_, double z_, double r_, BoundaryType bound,
         uint32_t i_id, std::string i_name);
  ~Sphere() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double x0, y0, z0, R;

};  // Sphere

//===========================================================================
// Non-Member Functions
std::shared_ptr<Sphere> make_sphere(YAML::Node surface_node);

#endif
