#ifndef MG_GEO_LILY_PAD_H
#define MG_GEO_LILY_PAD_H

#include <array>
#include <cstdint>
#include <utils/position.hpp>

struct GeoLilyPad {
  enum class PadType { Universe, Lattice, Cell };

  GeoLilyPad() {}
  GeoLilyPad(PadType t, uint32_t i, Position r, std::array<int32_t, 3> ti,
             bool ou)
      : type(t), id(i), r_local(r), tile(ti), in_lattice_outside_universe(ou) {}

  PadType type = PadType::Universe;
  uint32_t id = 0;

  // r_local is the position within the lattice,
  // NOT the position within the tile !!
  Position r_local = Position();
  std::array<int32_t, 3> tile = {0, 0, 0};
  bool in_lattice_outside_universe = false;
};

#endif