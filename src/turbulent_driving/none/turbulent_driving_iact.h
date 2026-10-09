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
#ifndef SWIFT_TURBULENT_DRIVING_IACT_NONE_H
#define SWIFT_TURBULENT_DRIVING_IACT_NONE_H

/* Local includes */
#include "part.h"
#include "turbulent_driving_struct.h"

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_turbulent_driving_density(
    const float r2, const float hi, const struct part *restrict pi,
    struct xpart *restrict xpi, const struct part *restrict pj) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_turbulent_driving_apply(
    const float r2, const float hi, const struct part *restrict pi,
    const struct xpart *restrict xpi, const struct part *restrict pj,
    struct xpart *restrict xpj) {}

#endif /* SWIFT_TURBULENT_DRIVING_IACT_NONE_H */
