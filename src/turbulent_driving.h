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
#ifndef SWIFT_TURBULENT_DRIVING_H
#define SWIFT_TURBULENT_DRIVING_H

/**
 * @file src/turbulent_driving.h
 * @brief Branches between the different turbulent driving models.
 *
 * A turbulent driving model mimics the effect of supernova feedback without
 * running star formation: gas particles satisfying a density/temperature
 * threshold within a given region stochastically trigger a thermal energy
 * injection event, shared between the particle and its SPH neighbours.
 */

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* Import the right turbulent driving definition */
#if defined(TURBULENT_DRIVING_NONE)
#include "./turbulent_driving/none/turbulent_driving.h"
#elif defined(TURBULENT_DRIVING_DEFAULT)
#include "./turbulent_driving/default/turbulent_driving.h"
#else
#error "Invalid choice of turbulent driving function."
#endif

/* General functions defined in the source file */
void turbulent_driving_init(struct swift_params *parameter_file,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            struct turbulent_driving_props *props);

void turbulent_driving_print(const struct turbulent_driving_props *props);

#endif /* SWIFT_TURBULENT_DRIVING_H */
