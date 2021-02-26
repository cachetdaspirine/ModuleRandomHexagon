#ifndef HEADER
#define HEADER

#ifdef DEBUG
#define DEBUG_IF(X) if(X)
#else
#define DEBUG_IF(X) if(false)
#endif

#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <tuple>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <map>
#include <random>
#include <stdio.h>
#include <iomanip>
#include <set>

#include "nr3.h"
#include "bracket.h"
#include "dbrent.h"
#include "df1dim.h"
#include "dlinmin.h"
#include "frprmn.h"

#include "Site.h"
#include "Node.h"
#include "Spring.h"
#include "Spring3.h"

#include "Ham.h"
#include "CG.h"
#include "System.h"

// Include the structural part of the simulation
#include "SiteAdjacency.h"
#include "NodeAdjacency.h"
#include "SpringAdjacency.h"

//Require the structure to know how the nodes are arrange



#endif
