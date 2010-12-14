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
#include <ctime>        //For time.
#include <stdlib.h>     //For rand.
#include <iostream>     //For standard i/o: cerr, cout, cin, endl...
#include <cstring>      //For memset, memcpy.
#include "atoms.h"
#include "potential.h"
#include "constants.h"
#include "common.h"
#include "integrator.h"
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
/* Class Integrator implementation {{{ */
Integrator::Integrator(const int n) {
    _t=1;
    _dt=1e-5;
    _dtOut=1e-3;
    _atoms=new Atoms(n);
    _potential=new Quadrupole(6.7e3);
}
Integrator::~Integrator(void) {
    if(_atoms!=0)
        delete _atoms;
    if(_potential!=0)
        delete _potential;
}
int Integrator::evolve(void) {
    srand(time(0));
    double t=0.;
    _potential->ePot(_atoms);
    cout << t << " " << *_atoms << " " << _atoms->eKin() << " "
        << _atoms->ePot() << "\n";
    double tOut=_dtOut;
    while(t<_t) {
        doSteps();
        t+=_dt;
        if(t>=tOut) {
            _potential->ePot(_atoms);
            cout << t << " " << *_atoms << " " << _atoms->eKin() << " "
                << _atoms->ePot() << "\n";
            tOut+=_dtOut;
        }
    }
    cout << "nan\n";
    return 0;
}
/* }}} */
/* Class RK2 implementation {{{ */
void RK2::doSteps(void) {
    int n=_atoms->n();
    double *acc=new double[3*n];
    double *oldpos=new double[3*n];
    double *oldvel=new double[3*n];
    double *pos=_atoms->pos();
    double *vel=_atoms->vel();
    _potential->forces(_atoms,acc);
    //First step: middle point evaluation.
    memcpy(oldpos,pos,3*n*sizeof(double));
    memcpy(oldvel,vel,3*n*sizeof(double));
    double dt=_dt*0.5;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]+=dt*vel[ii+d];
            vel[ii+d]+=dt*acc[ii+d];
        }
    }
    _potential->forces(_atoms,acc);
    //Second step: update atoms positions and velocities.
    dt=_dt;
    double v2=0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]=oldpos[ii+d]+dt*vel[ii+d];
            double v=oldvel[ii+d]+dt*acc[ii+d];
            vel[ii+d]=v;
            v2+=v*v;
        }
    }
    _atoms->eKin()=_atoms->m()*(0.5*mp/h)*v2/(double)n;
    delete[] oldvel;
    delete[] oldpos;
    delete[] acc;
};
/* }}} */
Integrator *initIntegrator(ConfigMap &config) {
    return new RK2(1);
}
/* integrator.cpp */
