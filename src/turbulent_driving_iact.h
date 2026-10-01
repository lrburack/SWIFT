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
#ifndef SWIFT_TURBULENT_DRIVING_IACT_H
#define SWIFT_TURBULENT_DRIVING_IACT_H

/**
 * @file src/turbulent_driving_iact.h
 * @brief Branches between the different turbulent driving interactions.
 */

/* Config parameters. */
#include <config.h>

/* Import the right turbulent driving definition */
#if defined(TURBULENT_DRIVING_NONE)
#include "./turbulent_driving/none/turbulent_driving_iact.h"
#elif defined(TURBULENT_DRIVING_DEFAULT)
#include "./turbulent_driving/default/turbulent_driving_iact.h"
#else
#error "Invalid choice of turbulent driving function."
#endif

#endif /* SWIFT_TURBULENT_DRIVING_IACT_H */
