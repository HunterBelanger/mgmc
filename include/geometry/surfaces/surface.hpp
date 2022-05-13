#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include <utils/direction.hpp>
#include <utils/position.hpp>

enum BoundaryType { Vacuum, Reflective, Normal };

class Surface {
 public:
  Surface(BoundaryType bound, uint32_t i_id, std::string i_name);
  virtual ~Surface() = default;

  virtual int sign(const Position &r, const Direction &u) const = 0;

  virtual double distance(const Position &r, const Direction &u,
                          bool on_surf) const = 0;

  virtual Direction norm(const Position &r) const = 0;

  BoundaryType boundary() const;
  uint32_t id() const;
  std::string name() const;

 protected:
  BoundaryType boundary_;
  uint32_t id_;
  std::string name_;

};  // Surface

#endif
