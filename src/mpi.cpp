/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
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
#ifdef MGMC_USE_MPI
#include <mpi.h>
#endif

#include <simulation/particle.hpp>
#include <type_traits>
#include <utils/error.hpp>
#include <utils/mpi.hpp>

namespace mpi {

Timer timer;

#ifdef MGMC_USE_MPI
const Com com = MPI_COMM_WORLD;
const DType Bool = MPI_C_BOOL;
const DType Int = MPI_INT;
const DType Double = MPI_DOUBLE;
const DType UInt64 = MPI_UINT64_T;
DType BParticle;

const OpType Sum = MPI_SUM;
const OpType And = MPI_LAND;
const OpType Or = MPI_LOR;
#endif

std::vector<uint64_t> node_nparticles;
std::vector<uint64_t> node_nparticles_noise;

int size = 1;
int rank = 0;

void register_banked_particle_type() {
#ifdef MGMC_USE_MPI
  // Ensure that BankedParticle is standard layout ! This is
  // required by MPI.
  static_assert(std::is_standard_layout<BankedParticle>::value);

  BankedParticle p;
  int err = 0;
  int sizes[13]{3, 3, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1};
  DType dtypes[13]{Double, Double, Double, Double, Double, UInt64, UInt64,
                   Bool,   Double, Double, Double, Double, Double};
  MPI_Aint disps[13];
  MPI_Get_address(&p.r, &disps[0]);
  MPI_Get_address(&p.u, &disps[1]);
  MPI_Get_address(&p.E, &disps[2]);
  MPI_Get_address(&p.wgt, &disps[3]);
  MPI_Get_address(&p.wgt2, &disps[4]);
  MPI_Get_address(&p.parent_history_id, &disps[5]);
  MPI_Get_address(&p.parent_daughter_id, &disps[6]);
  MPI_Get_address(&p.parents_previous_was_virtual, &disps[7]);
  MPI_Get_address(&p.parents_previous_position, &disps[8]);
  MPI_Get_address(&p.parents_previous_direction, &disps[9]);
  MPI_Get_address(&p.parents_previous_previous_energy, &disps[10]);
  MPI_Get_address(&p.parents_previous_energy, &disps[11]);
  MPI_Get_address(&p.Esmp_parent, &disps[12]);
  for (int i = 12; i >= 0; --i) {
    disps[i] -= disps[0];
  }

  DType tmp_BParticle;
  err = MPI_Type_create_struct(13, sizes, disps, dtypes, &tmp_BParticle);
  check_error(err, __FILE__, __LINE__);

  MPI_Aint lb, extnt;
  err = MPI_Type_get_extent(tmp_BParticle, &lb, &extnt);
  check_error(err, __FILE__, __LINE__);

  err = MPI_Type_create_resized(tmp_BParticle, lb, extnt, &BParticle);
  check_error(err, __FILE__, __LINE__);

  err = MPI_Type_commit(&BParticle);
  check_error(err, __FILE__, __LINE__);
#endif
}

void initialize_mpi(int* argc, char*** argv) {
  timer.reset();

#ifdef MGMC_USE_MPI
  int err = 0;

  // Initialize MPI
  err = MPI_Init(argc, argv);
  check_error(err, __FILE__, __LINE__);

  // Get worls size and our rank
  err = MPI_Comm_size(com, &size);
  check_error(err, __FILE__, __LINE__);

  err = MPI_Comm_rank(com, &rank);
  check_error(err, __FILE__, __LINE__);

  // Register baked particle type
  register_banked_particle_type();
#else
  (void)argc;
  (void)argv;
#endif
}

void finalize_mpi() {
#ifdef MGMC_USE_MPI
  MPI_Finalize();
#endif
}

void abort_mpi() {
#ifdef MGMC_USE_MPI
  MPI_Abort(com, 1);
#endif
}

void synchronize() {
#ifdef MGMC_USE_MPI
  if (size > 1) {
    timer.start();
    MPI_Barrier(com);
    timer.stop();
  }
#endif
}

void check_error(int err, const char* file, int line) {
#ifdef MGMC_USE_MPI
  if (err != MPI_SUCCESS) {
    // First, we should get the error string
    char err_str[MPI_MAX_ERROR_STRING];
    int str_len = 0;
    MPI_Error_string(err, err_str, &str_len);
    std::string mssg(err_str, str_len);
    fatal_error(mssg, file, line);
  }
#else
  (void)err;
  (void)file;
  (void)line;
#endif
}
}  // namespace mpi
