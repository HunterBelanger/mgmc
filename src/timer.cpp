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
#include <utils/timer.hpp>

//============================================================================
// Initialization of static members of Timer singleton
std::shared_ptr<Timer> Timer::timer_instance = nullptr;
std::mutex Timer::instance_mutex;

//============================================================================
// Timer Singleton Methods
Timer::Timer()
    : start_time{},
      start_parse_input{},
      end_parse_input{},
      start_transport{},
      end_transport{} {
  start_time = std::chrono::high_resolution_clock::now();
}

std::shared_ptr<Timer> Timer::instance() {
  if (timer_instance == nullptr) {
    instance_mutex.lock();
    if (timer_instance == nullptr) {
      timer_instance = std::shared_ptr<Timer>(new Timer());
    }
    instance_mutex.unlock();
  }
  return timer_instance;
}

double Timer::get_seconds(
    std::chrono::high_resolution_clock::time_point t0,
    std::chrono::high_resolution_clock::time_point t1) const {
  std::chrono::duration<double> diff = t1 - t0;
  diff = std::chrono::duration_cast<std::chrono::seconds>(diff);
  return diff.count();
}

void Timer::start_parsing_timer() {
  if (began_parsing == false) {
    start_parse_input = std::chrono::high_resolution_clock::now();
    began_parsing = true;
  }
}

void Timer::end_parsing_timer() {
  if (began_parsing && finished_parsing == false) {
    end_parse_input = std::chrono::high_resolution_clock::now();
    finished_parsing = true;
  }
}

double Timer::get_parse_time() const {
  if (began_parsing && finished_parsing) {
    return get_seconds(start_parse_input, end_parse_input);
  } else if (began_parsing) {
    std::chrono::high_resolution_clock::time_point now;
    now = std::chrono::high_resolution_clock::now();
    return get_seconds(start_parse_input, now);
  } else {
    return -1.;
  }
}

void Timer::start_transport_timer() {
  if (began_transport == false) {
    start_transport = std::chrono::high_resolution_clock::now();
    began_transport = true;
  }
}

void Timer::end_transport_timer() {
  if (began_transport && finished_transport == false) {
    end_transport = std::chrono::high_resolution_clock::now();
    finished_transport = true;
  }
}

double Timer::get_transport_time() const {
  if (began_transport && finished_transport) {
    return get_seconds(start_transport, end_transport);
  } else if (began_transport) {
    std::chrono::high_resolution_clock::time_point now;
    now = std::chrono::high_resolution_clock::now();
    return get_seconds(start_transport, now);
  } else {
    return -1.;
  }
}

double Timer::get_total_time() const {
  std::chrono::high_resolution_clock::time_point now;
  now = std::chrono::high_resolution_clock::now();
  return get_seconds(start_time, now);
}