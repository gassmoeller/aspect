/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_revision_h
#define _aspect_revision_h

// This configuration file will be copied to
// ${CMAKE_BINARY_DIR}/include/aspect/revision.h
// and filled with the current values by cmake upon configuration.

/**
 * Full version number of the ASPECT version.
 */
#define ASPECT_PACKAGE_VERSION "2.4.0-pre"

/**
 * Name of the local git branch of the source directory.
 */
#define ASPECT_GIT_BRANCH "AV_3D_simpler"

/**
 * Full sha1 revision of the current git HEAD.
 */
#define ASPECT_GIT_REVISION "c8834e073a0db7500d1dc9a2715f50a2bc3f759e"

/**
 * Short sha1 revision of the current git HEAD.
 */
#define ASPECT_GIT_SHORTREV "c8834e073"

#endif
