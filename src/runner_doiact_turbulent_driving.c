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
 * @file src/runner_doiact_turbulent_driving.c
 * @brief Gas-gas neighbour loops for the turbulent driving module.
 *
 * Unlike SWIFT's usual hydro/stars pair-interaction machinery, these loops
 * are deliberately simple, unvectorised, brute-force neighbour searches:
 * they trade some performance for a much smaller and easier to review
 * implementation, which is an acceptable trade-off for this single-node,
 * no-MPI, correctness-first module. They also assume non-periodic boundary
 * conditions (no nearest-image correction), matching the target simulation
 * this module was built for.
 */

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "part.h"
#include "timers.h"
#include "turbulent_driving_iact.h"

/**
 * @brief Brute-force gas-gas neighbour-count loop within a single cell.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_doself_turbulent_driving_density(struct runner *r, struct cell *c,
                                             int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  if (c->hydro.count == 0) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_doself_turbulent_driving_density(r, c->progeny[k], 0);
    if (timer) TIMER_TOC(timer_doself_turbulent_driving_density);
    return;
  }

  if (!cell_is_active_hydro(c, e)) {
    if (timer) TIMER_TOC(timer_doself_turbulent_driving_density);
    return;
  }

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  for (int i = 0; i < count; i++) {

    struct part *restrict pi = &parts[i];
    struct xpart *restrict xpi = &xparts[i];
    if (!xpi->turbulent_driving_data.triggered) continue;

    const float hi = pi->h;

    for (int j = 0; j < count; j++) {
      if (i == j) continue;

      const struct part *restrict pj = &parts[j];
      const float dx[3] = {pi->x[0] - pj->x[0], pi->x[1] - pj->x[1],
                           pi->x[2] - pj->x[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      runner_iact_nonsym_turbulent_driving_density(r2, hi, pi, xpi, pj);
    }
  }

  if (timer) TIMER_TOC(timer_doself_turbulent_driving_density);
}

/**
 * @brief Brute-force gas-gas neighbour-count loop between two cells.
 *
 * Recurses into progeny when either cell is split, and otherwise falls back
 * to a full double loop over the two (leaf) cells' particles.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param timer Are we timing this?
 */
void runner_dopair_turbulent_driving_density(struct runner *r,
                                             struct cell *ci, struct cell *cj,
                                             int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Nothing to do if neither side is active: no particle here could have
   * been marked as triggered this step. */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) {
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_density);
    return;
  }

  /* Recurse? */
  if (ci->split && cj->split) {
    for (int i = 0; i < 8; i++) {
      if (ci->progeny[i] == NULL) continue;
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] == NULL) continue;
        runner_dopair_turbulent_driving_density(r, ci->progeny[i],
                                                cj->progeny[j], 0);
      }
    }
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_density);
    return;
  } else if (ci->split) {
    for (int i = 0; i < 8; i++)
      if (ci->progeny[i] != NULL)
        runner_dopair_turbulent_driving_density(r, ci->progeny[i], cj, 0);
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_density);
    return;
  } else if (cj->split) {
    for (int j = 0; j < 8; j++)
      if (cj->progeny[j] != NULL)
        runner_dopair_turbulent_driving_density(r, ci, cj->progeny[j], 0);
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_density);
    return;
  }

  /* Both are leaves: brute-force double loop, checked in both directions
   * since either particle could be the triggering one. */
  struct part *restrict parts_i = ci->hydro.parts;
  struct xpart *restrict xparts_i = ci->hydro.xparts;
  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  for (int i = 0; i < count_i; i++) {

    struct part *restrict pi = &parts_i[i];
    struct xpart *restrict xpi = &xparts_i[i];
    const int i_triggered = xpi->turbulent_driving_data.triggered;
    const float hi = pi->h;

    for (int j = 0; j < count_j; j++) {

      struct part *restrict pj = &parts_j[j];
      struct xpart *restrict xpj = &xparts_j[j];

      const float dx[3] = {pi->x[0] - pj->x[0], pi->x[1] - pj->x[1],
                           pi->x[2] - pj->x[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (i_triggered)
        runner_iact_nonsym_turbulent_driving_density(r2, hi, pi, xpi, pj);

      if (xpj->turbulent_driving_data.triggered)
        runner_iact_nonsym_turbulent_driving_density(r2, pj->h, pj, xpj, pi);
    }
  }

  if (timer) TIMER_TOC(timer_dopair_turbulent_driving_density);
}

/**
 * @brief Brute-force gas-gas energy-apply loop within a single cell.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_doself_turbulent_driving_apply(struct runner *r, struct cell *c,
                                           int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  if (c->hydro.count == 0) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_doself_turbulent_driving_apply(r, c->progeny[k], 0);
    if (timer) TIMER_TOC(timer_doself_turbulent_driving_apply);
    return;
  }

  if (!cell_is_active_hydro(c, e)) {
    if (timer) TIMER_TOC(timer_doself_turbulent_driving_apply);
    return;
  }

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  for (int i = 0; i < count; i++) {

    struct part *restrict pi = &parts[i];
    struct xpart *restrict xpi = &xparts[i];
    if (!xpi->turbulent_driving_data.triggered) continue;

    const float hi = pi->h;

    for (int j = 0; j < count; j++) {
      if (i == j) continue;

      struct part *restrict pj = &parts[j];
      struct xpart *restrict xpj = &xparts[j];
      const float dx[3] = {pi->x[0] - pj->x[0], pi->x[1] - pj->x[1],
                           pi->x[2] - pj->x[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      runner_iact_nonsym_turbulent_driving_apply(r2, hi, pi, xpi, pj, xpj);
    }
  }

  if (timer) TIMER_TOC(timer_doself_turbulent_driving_apply);
}

/**
 * @brief Brute-force gas-gas energy-apply loop between two cells.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param timer Are we timing this?
 */
void runner_dopair_turbulent_driving_apply(struct runner *r, struct cell *ci,
                                           struct cell *cj, int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) {
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_apply);
    return;
  }

  /* Recurse? */
  if (ci->split && cj->split) {
    for (int i = 0; i < 8; i++) {
      if (ci->progeny[i] == NULL) continue;
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] == NULL) continue;
        runner_dopair_turbulent_driving_apply(r, ci->progeny[i],
                                              cj->progeny[j], 0);
      }
    }
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_apply);
    return;
  } else if (ci->split) {
    for (int i = 0; i < 8; i++)
      if (ci->progeny[i] != NULL)
        runner_dopair_turbulent_driving_apply(r, ci->progeny[i], cj, 0);
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_apply);
    return;
  } else if (cj->split) {
    for (int j = 0; j < 8; j++)
      if (cj->progeny[j] != NULL)
        runner_dopair_turbulent_driving_apply(r, ci, cj->progeny[j], 0);
    if (timer) TIMER_TOC(timer_dopair_turbulent_driving_apply);
    return;
  }

  struct part *restrict parts_i = ci->hydro.parts;
  struct xpart *restrict xparts_i = ci->hydro.xparts;
  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  for (int i = 0; i < count_i; i++) {

    struct part *restrict pi = &parts_i[i];
    struct xpart *restrict xpi = &xparts_i[i];
    const int i_triggered = xpi->turbulent_driving_data.triggered;
    const float hi = pi->h;

    for (int j = 0; j < count_j; j++) {

      struct part *restrict pj = &parts_j[j];
      struct xpart *restrict xpj = &xparts_j[j];

      const float dx[3] = {pi->x[0] - pj->x[0], pi->x[1] - pj->x[1],
                           pi->x[2] - pj->x[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (i_triggered)
        runner_iact_nonsym_turbulent_driving_apply(r2, hi, pi, xpi, pj, xpj);

      if (xpj->turbulent_driving_data.triggered)
        runner_iact_nonsym_turbulent_driving_apply(r2, pj->h, pj, xpj, pi,
                                                   xpi);
    }
  }

  if (timer) TIMER_TOC(timer_dopair_turbulent_driving_apply);
}
