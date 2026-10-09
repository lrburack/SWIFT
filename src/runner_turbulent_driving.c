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

/**
 * @file src/runner_turbulent_driving.c
 * @brief The three per-cell self-loop tasks of the turbulent driving
 * pipeline: trigger -> (density/apply neighbour loops, elsewhere) -> ghost
 * -> (apply neighbour loop, elsewhere) -> apply_ghost.
 *
 * The actual energy application (converting the accumulated
 * #turbulent_driving_xpart_data.delta_u_accum into a specific internal
 * energy increase) happens lazily in cell_drift.c, the same way GEAR
 * feedback defers its own energy application to the next drift. This
 * guarantees the energy is applied and the accumulator is cleared exactly
 * once, even for particles that received energy while inactive.
 */

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "black_holes.h"
#include "cell.h"
#include "cosmology.h"
#include "engine.h"
#include "mhd.h"
#include "part.h"
#include "sink.h"
#include "timeline.h"
#include "timers.h"
#include "timestep.h"
#include "turbulent_driving.h"

/**
 * @brief Decide which gas particles trigger a turbulent driving injection
 * event this step.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_turbulent_driving_trigger(struct runner *r, struct cell *c,
                                         int timer) {

  struct engine *e = r->e;

  TIMER_TIC;

  if (c->hydro.count == 0 || !cell_is_active_hydro(c, e)) {
    if (timer) TIMER_TOC(timer_do_turbulent_driving_trigger);
    return;
  }

  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_turbulent_driving_trigger(r, c->progeny[k], 0);
    return;
  }

  const struct turbulent_driving_props *props = e->turbulent_driving_props;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const struct space *s = e->s;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  for (int i = 0; i < count; i++) {

    struct part *restrict p = &parts[i];
    struct xpart *restrict xp = &xparts[i];

    if (!part_is_active(p, e)) continue;

    if (!turbulent_driving_is_eligible(p, xp, props, phys_const, cosmo,
                                       hydro_props, us, cooling, s)) {
      xp->turbulent_driving_data.triggered = 0;
      continue;
    }

    /* Time-step size for this particle, in internal units. */
    double dt;
    if (with_cosmology) {
      const integertime_t ti_step = get_integer_timestep(p->time_bin);
      const integertime_t ti_begin =
          get_integer_time_begin(ti_current - 1, p->time_bin);
      dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
    } else {
      dt = get_timestep(p->time_bin, time_base);
    }
    

    if (turbulent_driving_should_trigger(p, props, e, dt)) {
      xp->turbulent_driving_data.triggered = 1;
      xp->turbulent_driving_data.ngb_count = 0;
    } else {
      xp->turbulent_driving_data.triggered = 0;
    }
  }

  if (timer) TIMER_TOC(timer_do_turbulent_driving_trigger);
}

/**
 * @brief Finalise the per-particle energy share of triggered particles once
 * the neighbour count from the density loop is known, and apply the
 * self-term.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_turbulent_driving_ghost(struct runner *r, struct cell *c,
                                       int timer) {

  struct engine *e = r->e;

  TIMER_TIC;

  if (c->hydro.count == 0 || !cell_is_active_hydro(c, e)) {
    if (timer) TIMER_TOC(timer_do_turbulent_driving_ghost);
    return;
  }

  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_turbulent_driving_ghost(r, c->progeny[k], 0);
    return;
  }

  const struct turbulent_driving_props *props = e->turbulent_driving_props;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  for (int i = 0; i < count; i++) {

    struct part *restrict p = &parts[i];
    struct xpart *restrict xp = &xparts[i];

    if (!part_is_active(p, e)) continue;
    if (!xp->turbulent_driving_data.triggered) continue;

    const int n_total = xp->turbulent_driving_data.ngb_count + 1;
    const float share = props->energy_per_event / (float)n_total;
    xp->turbulent_driving_data.share = share;

    /* Self term: the triggering particle receives its own share too. */
    message("TURBULENT_DRIVING: adding energy");
    atomic_add_f(&xp->turbulent_driving_data.delta_u_accum, share);
  }

  if (timer) TIMER_TOC(timer_do_turbulent_driving_ghost);
}

/**
 * @brief Reset the "triggered" flag of particles that triggered an event
 * this step, now that the apply loop has consumed it.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_turbulent_driving_apply_ghost(struct runner *r, struct cell *c,
                                             int timer) {

  struct engine *e = r->e;

  TIMER_TIC;

  if (c->hydro.count == 0 || !cell_is_active_hydro(c, e)) {
    if (timer) TIMER_TOC(timer_do_turbulent_driving_apply_ghost);
    return;
  }

  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_turbulent_driving_apply_ghost(r, c->progeny[k], 0);
    return;
  }

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  for (int i = 0; i < count; i++) {

    struct part *restrict p = &parts[i];
    struct xpart *restrict xp = &xparts[i];

    if (!part_is_active(p, e)) continue;
    if (!xp->turbulent_driving_data.triggered) continue;

    xp->turbulent_driving_data.triggered = 0;
  }

  if (timer) TIMER_TOC(timer_do_turbulent_driving_apply_ghost);
}
