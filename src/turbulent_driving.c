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

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "turbulent_driving.h"

/**
 * @brief Initialises the turbulent driving properties in the internal unit
 * system.
 *
 * @param parameter_file The parsed parameter file.
 * @param phys_const Physical constants in internal units.
 * @param us The current internal system of units.
 * @param props The properties of the turbulent driving model to initialise.
 */
void turbulent_driving_init(struct swift_params *parameter_file,
                            const struct phys_const *phys_const,
                            const struct unit_system *us,
                            struct turbulent_driving_props *props) {

  turbulent_driving_init_backend(parameter_file, phys_const, us, props);
}

/**
 * @brief Print the properties of the turbulent driving model.
 *
 * @param props The properties of the turbulent driving model.
 */
void turbulent_driving_print(const struct turbulent_driving_props *props) {

  turbulent_driving_print_backend(props);
}
