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
#ifndef SWIFT_TURBULENT_DRIVING_NONE_H
#define SWIFT_TURBULENT_DRIVING_NONE_H

/**
 * @file src/turbulent_driving/none/turbulent_driving.h
 * @brief No-op turbulent driving model.
 *
 * Kept only so that the code compiles when no turbulent driving model is
 * selected. Enabling turbulent driving at runtime (--turbulent-driving) while
 * compiled with this backend is an error, just like STAR_FORMATION_NONE.
 */

/* Local includes */
#include "engine.h"
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "turbulent_driving_struct.h"
#include "units.h"

__attribute__((always_inline)) INLINE static int turbulent_driving_is_eligible(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct turbulent_driving_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const struct unit_system *us,
    const struct cooling_function_data *cooling, const struct space *s) {
  return 0;
}

__attribute__((always_inline)) INLINE static int
turbulent_driving_should_trigger(const struct part *restrict p,
                                 const struct turbulent_driving_props *props,
                                 const struct engine *e, const double dt) {
  return 0;
}

__attribute__((always_inline)) INLINE static void
turbulent_driving_update_part(struct part *p, struct xpart *xp,
                              const struct engine *e) {}

INLINE static void turbulent_driving_print_backend(
    const struct turbulent_driving_props *props) {
  message("Turbulent driving model is 'none'");
}

INLINE static void turbulent_driving_init_backend(
    struct swift_params *parameter_file, const struct phys_const *phys_const,
    const struct unit_system *us, struct turbulent_driving_props *props) {}

#endif /* SWIFT_TURBULENT_DRIVING_NONE_H */
