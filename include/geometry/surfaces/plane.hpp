#ifndef PLANE_H
#define PLANE_H

#include <yaml-cpp/yaml.h>

#include <geometry/surfaces/surface.hpp>

class Plane : public Surface {
 public:
  Plane(double A_, double B_, double C_, double D_, BoundaryType bound,
        uint32_t i_id, std::string i_name);
  ~Plane() = default;

  int sign(const Position &r, const Direction &u) const override;

  double distance(const Position &r, const Direction &u,
                  bool on_surf) const override;

  Direction norm(const Position &r) const override;

 private:
  double A, B, C, D;

};  // Plane

//===========================================================================
// Non-Member Functions
std::shared_ptr<Plane> make_plane(YAML::Node surface_node);

#endif
