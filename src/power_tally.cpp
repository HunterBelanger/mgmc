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
#include <simulation/power_tally.hpp>
#include <utils/output.hpp>

PowerTally::PowerTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                       uint64_t nz)
    : r_low{low},
      r_hi{hi},
      Nx{nx},
      Ny{ny},
      Nz{nz},
      dx(),
      dy(),
      dz(),
      g(),
      power_gen(),
      power_avg(),
      power_var() {
  dx = (r_hi.x() - r_low.x()) / static_cast<double>(Nx);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(Ny);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(Nz);

  // Allocate flux_gen and flux arrays
  power_gen.reallocate({Nx, Ny, Nz});
  power_avg.reallocate({Nx, Ny, Nz});
  power_var.reallocate({Nx, Ny, Nz});

  // Make sure both are zero
  power_gen.fill(0.);
  power_avg.fill(0.);
  power_var.fill(0.);
}

void PowerTally::score_power(const Particle& p,
                             const std::shared_ptr<Material>& mat) {
  double Et = mat->Et(p.r(), p.E());
  double Ef = mat->Ef(p.r(), p.E());

  if (Ef != 0.) {
    double scr = p.wgt() * Ef / Et;

    int i = std::floor((p.r().x() - r_low.x()) / dx);
    int j = std::floor((p.r().y() - r_low.y()) / dy);
    int k = std::floor((p.r().z() - r_low.z()) / dz);

    if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
        j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
      uint64_t ui = static_cast<uint64_t>(i);
      uint64_t uj = static_cast<uint64_t>(j);
      uint64_t uk = static_cast<uint64_t>(k);

#pragma omp atomic
      power_gen(ui, uj, uk) += scr;
    }
  }
}

void PowerTally::record_generation(double gen) {
  g = gen;
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < power_gen.size(); i++) {
    // Get new average
    double old_avg = power_avg[i];
    double val = power_gen[i];
    double avg = old_avg + (val - old_avg) / gen;
    power_avg[i] = avg;

    // Get new variance
    double var = power_var[i];
    var = var + (((val - old_avg) * (val - avg) - (var)) / gen);
    power_var[i] = var;
  }
}

void PowerTally::clear_generation() { power_gen.fill(0.); }

void PowerTally::write_power(std::string fname) {
  Output::instance()->write(" Writing power distribution file...\n");

  std::ofstream flxf(fname + "_mesh.txt");

  // First write coordinates and number of groups
  flxf << " X:";
  for (int i = 0; i < static_cast<int>(Nx) - 1; i++) {
    double x = i * dx + (0.5 * dx) + r_low.x();
    flxf << x << ",";
  }
  flxf << (Nx - 1) * dx + (0.5 * dx) + r_low.x() << "\n";

  flxf << " Y:";
  for (int i = 0; i < static_cast<int>(Ny) - 1; i++) {
    double y = i * dy + (0.5 * dy) + r_low.y();
    flxf << y << ",";
  }
  flxf << (Ny - 1) * dy + (0.5 * dy) + r_low.y() << "\n";

  flxf << " Z:";
  for (int i = 0; i < static_cast<int>(Nz) - 1; i++) {
    double z = i * dz + (0.5 * dz) + r_low.z();
    flxf << z << ",";
  }
  flxf << (Nz - 1) * dz + (0.5 * dz) + r_low.z() << "\n";

  flxf.close();

  // Convert power_var to error on the mean for saving
  for (size_t l = 0; l < power_var.size(); l++)
    power_var[l] = std::sqrt(power_var[l] / g);

  power_avg.save(fname + "_avg.npy");
  power_var.save(fname + "_err.npy");
}