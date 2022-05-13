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
#ifndef SETTINGS_H
#define SETTINGS_H

#include <pcg_random.hpp>
#include <string>
#include <utils/timer.hpp>

namespace settings {
enum class SimulationMode { K_EIGENVALUE, FIXED_SOURCE, NOISE };
enum class TrackingMode { SURFACE_TRACKING, DELTA_TRACKING, CARTER_TRACKING };
enum class EnergyMode { CE, MG };

extern int nparticles;
extern int ngenerations;
extern int nignored;
extern int nskip;
extern uint32_t ngroups;
extern int n_cancel_noise_gens;

extern Timer alpha_omega_timer;
extern double max_time;

extern double min_energy;
extern double max_energy;
extern double target_at_rest_threshold;

extern SimulationMode mode;
extern TrackingMode tracking;
extern EnergyMode energy_mode;

extern uint64_t rng_seed;
extern uint64_t rng_stride;
extern pcg32 rng;

extern double wgt_cutoff;
extern double wgt_survival;
extern double wgt_split;

extern double w_noise;
extern double eta;
extern double keff;

extern bool converged;

extern bool regional_cancellation;
extern bool regional_cancellation_noise;
extern bool inner_generations;
extern bool branchless_noise;
extern bool normalize_noise_source;
extern bool rng_stride_warnings;
extern bool save_source;
extern bool load_source_file;

// Energy bounds for multi-group mode
extern std::vector<double> energy_bounds;
extern std::vector<double> mg_speeds;
extern std::vector<double> sample_xs_ratio;

extern std::string output_file_name;
extern std::string source_file_name;
extern std::string in_source_file_name;

void initialize_global_rng();
}  // namespace settings

#endif  // MG_SETTINGS_H
