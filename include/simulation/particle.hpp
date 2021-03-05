/*=============================================================================*
 * Copyright (C) 2021, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est un programme informatique servant à faire des comparaisons
 * entre les méthodes de transport qui sont capable de traiter les milieux
 * continus avec la méthode Monte Carlo. Il résoud l'équation de Boltzmann
 * pour les particules neutres, à une vitesse et dans une dimension.
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#ifndef MG_PARTICLE_H
#define MG_PARTICLE_H

#include <utils/direction.hpp>
#include <utils/position.hpp>

class Particle {
 public:
  Particle(Position r, Direction u, int engy, double wgt)
      : position{r}, direction{u}, energy{engy}, weight{wgt} {};

  Position& r() { return position; }
  const Position& r() const { return position; }
  Direction& u() { return direction; }
  double wgt() const { return weight; }
  bool is_alive() const { return alive; }
  bool is_reflected() const { return reflected; }
  int E() const { return energy; }
  Position previous_r() const { return previous_position; }
  int previous_E() const { return previous_energy; }

  void set_position(Position r) {
    previous_position = position;
    position = r;
  }
  void set_direction(Direction u) { direction = u; }
  void set_weight(double w) { weight = w; }
  void set_energy(int E) {
    previous_energy = energy;
    energy = E;
  }
  void set_previous_r(Position r) { previous_position = r; }
  void set_reflected(bool ref) { reflected = ref; }

  void move(double dist) {
    if (reflected == false) {
      previous_position = position;
      position = position + dist * direction;
    } else {
      position = position + dist * direction;
      reflected = false;
    }
  }

  void kill() { alive = false; }

  Position parents_previous_position = Position();
  int parents_previous_energy = 0;
  double Esmp_parent = 1.;

 private:
  Position position;
  Direction direction;

  Position previous_position = Position();
  int previous_energy = 0;

  int energy;
  double weight;
  bool alive = true;
  bool reflected = false;

};  // Particle

#endif  // MG_PARTICLE_H
