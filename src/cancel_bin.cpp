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
#include <cmath>
#include <simulation/cancel_bin.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

CancelBin::CancelBin() : bin{}, uniform_weight{0.}, material{nullptr} {}

CancelBin::CancelBin(std::shared_ptr<Material> mat)
    : bin{}, uniform_weight{0.}, material{mat} {}

CancelBin::CancelBin(const CancelBin& other)
    : bin{}, uniform_weight{0.}, material{other.material} {}

CancelBin& CancelBin::operator=(const CancelBin& other) {
  uniform_weight = 0.;
  bin = {};
  material = other.material;

  return *this;
}

void CancelBin::add_particle(Particle& p, double Xl, double Xh, double Yl,
                             double Yh, double Zl, double Zh) {
  // Make sure parent was not inside bin. If it was,
  // then dont add it to the bin to be cancelled.
  if (!is_inside(p.parents_previous_position, Xl, Xh, Yl, Yh, Zl, Zh)) {
    bin.push_back(&p);
  }
}

void CancelBin::perform_cancelation(double Xl, double Xh, double Yl, double Yh,
                                    double Zl, double Zh) {
  if (!material) {
    double x = (Xh + Xl) * 0.5;
    double y = (Yh + Yl) * 0.5;
    double z = (Zh + Zl) * 0.5;
    std::string mssg = "Cancellation bin for coordinates ";
    mssg += "x = " + std::to_string(x) + ", y = " + std::to_string(y);
    mssg += ", z = " + std::to_string(z);
    mssg += " has no material.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  uniform_weight = 0.;
  for (auto& p : bin) {
    Position r_parent = p->parents_previous_position;
    int E = p->parents_previous_energy;
    double Esmp = p->Esmp_parent;

    // Technically, this should be determined based on the corner points
    // c1, c2, etc. or by p->r(), but in this constant XS MG case, it doesn't
    // matter.
    double Ef = material->Ef(r_parent, E);

    // Determine Beta value and f value for particle
    double B = get_beta(r_parent, Esmp, Ef, Xl, Xh, Yl, Yh, Zl, Zh);
    double f = get_f(p->r(), r_parent, Esmp, Ef);

    // Tally the uniform portion of weight
    uniform_weight += p->wgt() * (B / f);

    // Adjust particle weight
    p->set_weight(p->wgt() * (f - B) / f);
  }

  // All weights have been changed. We can now reset bin
  bin.clear();
}

std::vector<Particle> CancelBin::get_uniform_particles(std::shared_ptr<RNG> rng,
                                                       double Xl, double Xh,
                                                       double Yl, double Yh,
                                                       double Zl, double Zh) {
  std::vector<Particle> prts_to_return;

  // Determine number of new particles to add
  uint32_t N = std::ceil(std::abs(uniform_weight));
  double w = uniform_weight / N;

  prts_to_return.reserve(N);

  for (size_t i = 0; i < N; i++) {
    // Sample position along z
    double x = (Xh - Xl) * rng->rand() + Xl;
    double y = (Yh - Yl) * rng->rand() + Yl;
    double z = (Zh - Zl) * rng->rand() + Zl;

    // Sample random direction
    double mu = 2. * rng->rand() - 1.;
    double phi = 2. * PI * rng->rand();
    double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
    double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
    double uz = mu;

    // Sample particle energy
    int E = rng->discrete(material->chi({x, y, z}, 0));

    // Save sampled particle
    prts_to_return.push_back(Particle({x, y, z}, {ux, uy, uz}, E, w));
  }

  uniform_weight = 0.;

  return prts_to_return;
}

bool CancelBin::is_inside(const Position& r, double Xl, double Xh, double Yl,
                          double Yh, double Zl, double Zh) {
  if (r.x() < Xl || r.x() > Xh) return false;

  if (r.y() < Yl || r.y() > Yh) return false;

  if (r.z() < Zl || r.z() > Zh) return false;

  return true;
}

double CancelBin::get_beta(Position r_parent, double Esmp, double Ef, double Xl,
                           double Xh, double Yl, double Yh, double Zl,
                           double Zh) {
  // We now need the sampling XS and total XS at the location of d.
  // This is generally position dependent, but in our multi-group,
  // constant density/temperature case, the general XSs for the
  // region will sufice.

  if (Esmp <= 0. || Ef <= 0.) {
    std::string mssg = "Encountered negative cross section in Beta sampling.";
    if (Esmp <= 0.) {
      mssg += " Esmp = " + std::to_string(Esmp) + ".";
    } else {
      mssg += " Ef = " + std::to_string(Ef) + ".";
    }
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Create 8 positions for cube corners
  Position c1(Xh, Yh, Zh);
  Position c2(Xh, Yh, Zl);
  Position c3(Xh, Yl, Zh);
  Position c4(Xh, Yl, Zl);
  Position c5(Xl, Yh, Zh);
  Position c6(Xl, Yh, Zl);
  Position c7(Xl, Yl, Zh);
  Position c8(Xl, Yl, Zl);

  double f1 = get_f(c1, r_parent, Esmp, Ef);
  double f2 = get_f(c2, r_parent, Esmp, Ef);
  double f3 = get_f(c3, r_parent, Esmp, Ef);
  double f4 = get_f(c4, r_parent, Esmp, Ef);
  double f5 = get_f(c5, r_parent, Esmp, Ef);
  double f6 = get_f(c6, r_parent, Esmp, Ef);
  double f7 = get_f(c7, r_parent, Esmp, Ef);
  double f8 = get_f(c8, r_parent, Esmp, Ef);

  double beta = std::min(std::min(std::min(f1, f2), std::min(f3, f4)),
                         std::min(std::min(f5, f6), std::min(f7, f8)));

  return beta;
}

double CancelBin::get_f(Position r, Position r_parent, double Esmp, double Ef) {
  double d = std::sqrt(std::pow(r.x() - r_parent.x(), 2.) +
                       std::pow(r.y() - r_parent.y(), 2.) +
                       std::pow(r.z() - r_parent.z(), 2.));

  // Can now calculate point fission reaction rate
  double f = (1. / (4. * PI * d * d)) * std::exp(-Esmp * d) * Ef;

  return f;
}
