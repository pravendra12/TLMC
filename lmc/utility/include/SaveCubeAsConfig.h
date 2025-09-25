#pragma once

#include "Config.h"
#include "TiledSupercell.h"
#include "Cube.h"

void SaveCubeAsConfig(
    const string &filename,
    const Cube &cubeObj,
    const TiledSupercell &tiledSupercell);
