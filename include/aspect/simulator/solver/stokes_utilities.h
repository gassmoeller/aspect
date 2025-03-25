/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_simulator_solver_stokes_utilities_h
#define _aspect_simulator_solver_stokes_utilities_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator/solver/interface.h>


namespace aspect
{
  namespace StokesSolverUtilities
  {
    template<int dim>
    bool
    is_stokes_matrix_free(const StokesSolver::Interface<dim> &solver);
  }
}

#endif
