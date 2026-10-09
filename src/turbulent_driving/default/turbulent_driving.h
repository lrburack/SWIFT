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
#ifndef SWIFT_TURBULENT_DRIVING_DEFAULT_H
#define SWIFT_TURBULENT_DRIVING_DEFAULT_H

/**
 * @file src/turbulent_driving/default/turbulent_driving.h
 * @brief Stochastic thermal energy injection ("SN-mimicking") turbulent
 * driving model.
 *
 * For every gas particle satisfying a density and temperature threshold
 * within a given region, each active step draws a stochastic trigger with
 * probability `efficiency * dt / t_ff` (t_ff the local free-fall time, in
 * the same spirit as the GEAR stochastic star formation model). On trigger,
 * a fixed amount of energy is deposited, split equally between the particle
 * and its SPH neighbours (found via the particle's usual smoothing length).
 */

/* Local includes */
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "pressure_floor.h"
#include "random.h"
#include "space.h"
#include "turbulent_driving_struct.h"
#include "units.h"

/**
 * @brief Compute the position of the centre of the eligible region.
 *
 * @param props The #turbulent_driving_props.
 * @param s The #space (for the box dimensions).
 * @param centre (return) The centre, in internal length units.
 */
__attribute__((always_inline)) INLINE static void
turbulent_driving_get_centre(const struct turbulent_driving_props *props,
                             const struct space *s, double centre[3]) {

  if (props->use_relative_centre) {
    centre[0] = 0.5 * s->dim[0] + props->centre[0];
    centre[1] = 0.5 * s->dim[1] + props->centre[1];
    centre[2] = 0.5 * s->dim[2] + props->centre[2];
  } else {
    centre[0] = props->centre[0];
    centre[1] = props->centre[1];
    centre[2] = props->centre[2];
  }
}

/**
 * @brief Is this particle eligible for a turbulent driving injection event?
 *
 * Checks the density, temperature and radial distance from the centre of
 * the eligible region.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param props The #turbulent_driving_props.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal unit system.
 * @param cooling The cooling function properties.
 * @param s The #space (for the box dimensions).
 */
__attribute__((always_inline)) INLINE static int turbulent_driving_is_eligible(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct turbulent_driving_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const struct unit_system *us,
    const struct cooling_function_data *cooling, const struct space *s) {

  /* Density criterion */
  const float density = hydro_get_physical_density(p, cosmo);
  if (density < props->density_threshold) return 0;

  /* Temperature criterion (note: unlike star formation, we require the gas
   * to be ABOVE the threshold, mimicking gas already heated by a previous
   * feedback-like event). */
  const float temperature = cooling_get_temperature(phys_const, hydro_props,
                                                     us, cosmo, cooling, p, xp);
  if (temperature < props->temperature_threshold) return 0;

  /* Spatial criterion: distance from the centre of the eligible region */
  double centre[3];
  turbulent_driving_get_centre(props, s, centre);

  const double dx = p->x[0] - centre[0];
  const double dy = p->x[1] - centre[1];
  const double dz = p->x[2] - centre[2];
  const double r2 = dx * dx + dy * dy + dz * dz;

  return r2 < (double)(props->radius * props->radius);
}

/**
 * @brief Decide whether an eligible particle triggers an injection event
 * this step.
 *
 * Draws a stochastic trigger with probability `efficiency * dt / t_ff`,
 * with `t_ff = sqrt(3*pi / (32*G*rho))` the local free-fall time, mirroring
 * the stochastic star formation probability used elsewhere in SWIFT.
 *
 * @param p The #part.
 * @param props The #turbulent_driving_props.
 * @param e The #engine (for physical constants, cosmology and random
 * numbers).
 * @param dt The time-step of this particle, in internal units.
 */
__attribute__((always_inline)) INLINE static int
turbulent_driving_should_trigger(const struct part *restrict p,
                                 const struct turbulent_driving_props *props,
                                 const struct engine *e, const double dt) {

  if (dt <= 0.) return 0;

  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;

  const float G = phys_const->const_newton_G;
  const float density = hydro_get_physical_density(p, cosmo);

  /* inv_t_ff = sqrt(32 * G * rho / (3 * pi)) */
  const float inv_t_ff =
      sqrtf(density * 32.f * G * 0.33333333f * (float)M_1_PI);

  const float prob = props->efficiency * inv_t_ff * (float)dt;

  const float random_number = random_unit_interval(
      p->id, e->ti_current, random_number_turbulent_driving);

  return random_number < prob;
}

/**
 * @brief Apply any turbulent driving energy accumulated onto a particle by
 * its triggering neighbours, and reset the accumulator.
 *
 * Called on every particle whenever it is drifted (see cell_drift.c),
 * regardless of whether it is currently active. This mirrors the deferred
 * feedback-application pattern used by e.g. GEAR thermal feedback, and
 * guarantees energy deposited on an inactive neighbour is neither lost nor
 * applied more than once.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
turbulent_driving_update_part(struct part *p, struct xpart *xp,
                              const struct engine *e) {

  if (xp->turbulent_driving_data.delta_u_accum == 0.f) return;

  const struct cosmology *cosmo = e->cosmology;
  const float mass = hydro_get_mass(p);
  const float du = xp->turbulent_driving_data.delta_u_accum / mass;

  const float u_old = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float u_new = u_old + du;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo,
                                             e->pressure_floor_props, u_new);

  xp->turbulent_driving_data.delta_u_accum = 0.f;
}

/**
 * @brief Prints the properties of the turbulent driving model to stdout.
 *
 * @param props The #turbulent_driving_props.
 */
INLINE static void turbulent_driving_print_backend(
    const struct turbulent_driving_props *props) {
  message("Turbulent driving model is 'default' (stochastic SN-mimicking)");
}

/**
 * @brief Initialises the turbulent driving properties from the parameter
 * file.
 *
 * @param parameter_file The parsed parameter file.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 * @param props The #turbulent_driving_props to initialise.
 */
INLINE static void turbulent_driving_init_backend(
    struct swift_params *parameter_file, const struct phys_const *phys_const,
    const struct unit_system *us, struct turbulent_driving_props *props) {


  message("TURBULENT_DRIVING: initializing");

  /* Stochastic efficiency (dimensionless, paper value: 0.02) */
  props->efficiency = parser_get_param_float(
      parameter_file, "TurbulentDriving:efficiency");

  /* Density threshold, in Hydrogen atoms per cm^3. Converted below using the
   * same convention as GEARStarFormation: density_threshold = n_H * m_p, with
   * no explicit correction for the Hydrogen mass fraction. */
  props->density_threshold = parser_get_param_float(
      parameter_file, "TurbulentDriving:density_threshold_Hpcm3");

  /* Minimum temperature, in Kelvin */
  props->temperature_threshold = parser_get_param_float(
      parameter_file, "TurbulentDriving:temperature_threshold_K");

  /* Radius of the eligible region, in kpc (internal length units are assumed
   * to already be kpc for this model; the value is used as-is). */
  props->radius =
      parser_get_param_float(parameter_file, "TurbulentDriving:radius");

  /* Energy injected per event, in erg. Read as a double: 1e51 erg overflows
   * a 32-bit float before it has a chance to be converted down to the
   * (much smaller) internal unit system below. */
  const double energy_per_event_erg = parser_get_opt_param_double(
      parameter_file, "TurbulentDriving:energy_per_event_erg", 1e51);

  /* Centre of the eligible region */
  props->use_relative_centre = parser_get_opt_param_int(
      parameter_file, "TurbulentDriving:use_relative_centre", 1);

  double centre[3] = {0., 0., 0.};
  parser_get_opt_param_double_array(parameter_file, "TurbulentDriving:centre",
                                    3, centre);
  props->centre[0] = centre[0];
  props->centre[1] = centre[1];
  props->centre[2] = centre[2];

  /* Convert the density threshold from Hydrogen atoms per cm^3 to internal
   * (physical) density units. */
  const double m_p_cgs = phys_const->const_proton_mass *
                         units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  props->density_threshold *=
      m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  /* Convert the temperature threshold to internal units. */
  props->temperature_threshold /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Convert the energy per event from erg (cgs) to internal units. */
  props->energy_per_event =
      energy_per_event_erg / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  if (engine_rank == 0) {
    message("efficiency                = %g", props->efficiency);
    message("density_threshold         = %g", props->density_threshold);
    message("temperature_threshold     = %g", props->temperature_threshold);
    message("radius                    = %g", props->radius);
    message("energy_per_event          = %g", props->energy_per_event);
    message("use_relative_centre       = %d", props->use_relative_centre);
    message("centre                    = (%g, %g, %g)", props->centre[0],
            props->centre[1], props->centre[2]);
  }
}

#endif /* SWIFT_TURBULENT_DRIVING_DEFAULT_H */
