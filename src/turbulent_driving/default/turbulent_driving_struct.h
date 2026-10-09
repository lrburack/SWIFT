/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT contributors
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_TURBULENT_DRIVING_STRUCT_DEFAULT_H
#define SWIFT_TURBULENT_DRIVING_STRUCT_DEFAULT_H

/**
 * @brief Extra turbulent driving fields carried by each hydro particle.
 *
 * All fields are transient, single-timestep scratch data: they are (re-)set
 * by the trigger task, filled in by the density/apply neighbour loops, and
 * consumed by the ghost/apply_ghost tasks within the same step. None of this
 * needs to survive a restart.
 */
struct turbulent_driving_xpart_data {

  /*! Did this particle trigger an injection event this step? */
  char triggered;

  /*! Number of gas neighbours found within the kernel of a triggered
   *  particle (filled by the density loop, valid only if triggered). */
  int ngb_count;

  /*! Energy (internal units) that each of the ngb_count + 1 recipients
   *  (the particle itself and its neighbours) receives, valid only if
   *  triggered. Computed by the ghost task once ngb_count is known. */
  float share;

  /*! Energy (internal units) accumulated from neighbouring triggered
   *  particles, to be turned into a specific internal energy increase and
   *  applied by the apply_ghost task. */
  float delta_u_accum;
};

/**
 * @brief Global turbulent driving properties, read from the parameter file.
 */
struct turbulent_driving_props {

  /*! Dimensionless efficiency entering the trigger probability
   *  (prob = efficiency * dt / t_ff). */
  float efficiency;

  /*! Density threshold, in internal (physical) density units, converted
   *  from the Hydrogen-number-density-in-cm^-3 parameter the same way
   *  GEARStarFormation does (i.e. density_threshold = n_H_threshold * m_p,
   *  no explicit correction for the Hydrogen mass fraction). */
  float density_threshold;

  /*! Minimum temperature, in internal units, above which a particle is
   *  eligible for injection. */
  float temperature_threshold;

  /*! Radius, in internal length units, within which a particle is
   *  eligible for injection. */
  float radius;

  /*! Energy injected per triggered event, in internal energy units. */
  float energy_per_event;

  /*! Is the centre specified relative to the box centre (1) or as an
   *  absolute position in the box (0)? */
  int use_relative_centre;

  /*! Centre of the region eligible for injection, in internal length
   *  units. Either an absolute position, or an offset from the box
   *  centre, depending on use_relative_centre. */
  double centre[3];
};

#endif /* SWIFT_TURBULENT_DRIVING_STRUCT_DEFAULT_H */
