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
#ifndef SWIFT_TURBULENT_DRIVING_IACT_DEFAULT_H
#define SWIFT_TURBULENT_DRIVING_IACT_DEFAULT_H

/* Local includes */
#include "atomic.h"
#include "kernel_hydro.h"
#include "part.h"
#include "turbulent_driving_struct.h"

/**
 * @brief Gas-gas interaction for the turbulent driving neighbour-count loop.
 *
 * If particle @p pi triggered an injection event this step, and @p pj lies
 * within pi's kernel, increment pi's neighbour counter. This is the
 * "nonsym" half of the interaction: callers must invoke this once with
 * (i, j) and once with (j, i) to cover both particles potentially being
 * triggering particles.
 *
 * @param r2 Square of the distance between the two particles.
 * @param hi Smoothing length of pi.
 * @param pi First particle (mutated: its neighbour counter is incremented).
 * @param xpi Extended data of pi.
 * @param pj Second particle (read-only).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_turbulent_driving_density(
    const float r2, const float hi, const struct part *restrict pi,
    struct xpart *restrict xpi, const struct part *restrict pj) {

  /* Nothing to do if pi is not a triggering particle this step. */
  if (!xpi->turbulent_driving_data.triggered) return;

  /* Ignore self-interactions. */
  if (pi->id == pj->id) return;

  /* Is pj within pi's kernel? */
  const float hig2 = hi * hi * kernel_gamma * kernel_gamma;
  if (r2 >= hig2) return;

  atomic_inc(&xpi->turbulent_driving_data.ngb_count);
}

/**
 * @brief Gas-gas interaction for the turbulent driving energy-apply loop.
 *
 * If particle @p pi triggered an injection event this step, and @p pj lies
 * within pi's kernel, add pi's per-neighbour energy share to pj's energy
 * accumulator. As above, this is the "nonsym" half of the interaction:
 * callers must invoke this once with (i, j) and once with (j, i).
 *
 * @param r2 Square of the distance between the two particles.
 * @param hi Smoothing length of pi.
 * @param pi First particle (read-only: its share is read, not modified).
 * @param xpi Extended data of pi.
 * @param pj Second particle (read-only, for the id check).
 * @param xpj Extended data of pj (mutated: its energy accumulator grows).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_turbulent_driving_apply(
    const float r2, const float hi, const struct part *restrict pi,
    const struct xpart *restrict xpi, const struct part *restrict pj,
    struct xpart *restrict xpj) {

  if (!xpi->turbulent_driving_data.triggered) return;

  if (pi->id == pj->id) return;

  const float hig2 = hi * hi * kernel_gamma * kernel_gamma;
  if (r2 >= hig2) return;

  atomic_add_f(&xpj->turbulent_driving_data.delta_u_accum,
              xpi->turbulent_driving_data.share);
}

#endif /* SWIFT_TURBULENT_DRIVING_IACT_DEFAULT_H */
