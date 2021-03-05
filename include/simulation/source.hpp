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
#ifndef MG_SOURCE_H
#define MG_SOURCE_H

#include <fstream>
#include <memory>
#include <simulation/particle.hpp>
#include <utils/constants.hpp>
#include <utils/rng.hpp>

class Source {
 public:
  Source(double w) : weight{w} {};
  virtual ~Source() = default;
  virtual Particle generate_particle(std::shared_ptr<RNG> rng) = 0;
  virtual double wgt() const { return weight; }

 protected:
  double weight;

};  // Source

class PointSource : public Source {
 public:
  PointSource(Position i_r, std::vector<double> i_x, double w)
      : Source{w}, r{i_r}, chi{i_x} {};
  ~PointSource() = default;

  Particle generate_particle(std::shared_ptr<RNG> rng) override {
    double mu = 2. * rng->rand() - 1.;
    double phi = 2. * PI * rng->rand();
    double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
    double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
    double uz = mu;

    int E = rng->discrete(chi);

    return Particle(r, Direction(ux, uy, uz), E, 1.0);
  }

  double wgt() const override { return weight; }

 private:
  Position r;
  std::vector<double> chi;

};  // PointSource

class BoxSource : public Source {
 public:
  BoxSource(Position low, Position hi, std::vector<double> i_x, double w)
      : Source{w}, low_left{low}, up_right{hi}, chi{i_x} {}
  ~BoxSource() = default;

  Particle generate_particle(std::shared_ptr<RNG> rng) override {
    // Sample random positon in box
    double x = (up_right.x() - low_left.x()) * rng->rand() + low_left.x();
    double y = (up_right.y() - low_left.y()) * rng->rand() + low_left.y();
    double z = (up_right.z() - low_left.z()) * rng->rand() + low_left.z();
    Position r(x, y, z);

    // Sample random direction
    double mu = 2. * rng->rand() - 1.;
    double phi = 2. * PI * rng->rand();
    double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
    double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
    double uz = mu;
    Direction u(ux, uy, uz);

    // Sample random group
    int E = rng->discrete(chi);

    return Particle(r, u, E, 1.0);
  }

  double wgt() const override { return weight; }

 private:
  Position low_left;
  Position up_right;
  std::vector<double> chi;
};

#endif  // MG_SOURCE_H
