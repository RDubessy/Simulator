/* This file is a part of Simulator. {{{
 * Copyright (C) 2010 Romain Dubessy
 *
 * findMinimum is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *
 * findMinimum is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with findMinimum.  If not, see <http://www.gnu.org/licenses/>.
 *
 * }}} */
/* Documentation {{{ */
/*! \mainpage  Simulator documentation.
 *
 * \section Installation
 * Download the latest version of this code and extract the archive.
 * In the main folder simply issue the command:
 * \code
 * make all
 * \endcode
 * If you want this program to be installed in your /usr/local/bin directory, 
 * use (as root):
 * \code
 * make install
 * \endcode
 * (of course this directory must exists).
 *
 * See the src/Makefile file for documentation on the various possible
 * optimizations.
 *
 * \section Usage
 * Usage information may be obtained issuing:
 * \code
 * findMinimum --usage
 * \endcode
 * The generic usage of findMinimum would be:
 * \code
 * findMinimum options [--key=value] file
 * \endcode
 * Where options is one of the arguments detailled hereafter and file is a 
 * configuration file.
 * The arguments may be passed inside the configuration file or as command line
 * options and in this case they override the contents of the file.
 * See the examples to find out how the configuration file is formatted.
 * The arguments may be any of those:
 * - foo
 * - blah
 *
 * \subsection Warnings
 * The warnings are displayed on the program standard error output.
 * Small ones are labelled by a <code>[W]</code> and more serious ones are
 * labelled by a <code>[E]</code>.
 */
/* }}} */
#include <iostream>
#include "common.h"
#include "integrator.h"
int main(int argc, char *argv[]) {
    ConfigMap config;
    if(!parseOptions(argc,argv,config)) {        //Parse cmd line options
        std::cerr << "==> Try '" << argv[0] << " --usage'" << std::endl;
        return -1;
    }
    if(config["usage"].size()>0||config["help"].size()>0) {
        printUsage(argv[0]);
        return -1;
    }
    if(!parseConfig(config)) {
        std::cerr << "==> Try 'man " << argv[0] << "'" << std::endl;
        return -1;
    }
    Integrator *integrator=initIntegrator(config);
    return integrator->evolve();
}
/* main.cpp */
