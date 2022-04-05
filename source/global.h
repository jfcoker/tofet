///////////////////////////////////////////////////////////////////////
//  This file is part of ToFeT.
//  
//  ToFeT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  ToFeT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public License
//  along with ToFeT.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////
#ifndef _GLOBAL_H
#define	_GLOBAL_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include <string>
#include <cstring>
#include <cmath>
#include "RandomB.h"
#ifndef RandomB
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
extern gsl_rng * gslRand;
#endif

// Clean exit on timeout or program kill
#include <mutex>
#include <chrono>
#include <thread>
#include <csignal>

template<typename T> void safe_delete(T& );
template<typename T> void clearListList(T&obj);
const double hbar_js = 1.054571628e-34;
const double hbar_eVs= 6.58211899e-16;
const double k_eVK   = 8.617343e-5;
const double pi      = 3.14159265;
const double e       = 1.60217646e-19;
extern bool VERBOSITY_HIGH;
extern int WARNINGS;
//#define printTotalOccupation

// Print an error message and exit program.
void ERROR(int code, std::string msg);

// Try to output results on receiving terminate signal.
extern bool RECEIVED_TERM_SIGNAL;
void signal_handler(int s);

#endif	/* _GLOBAL_H */
