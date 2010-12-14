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
#include <cmath>                //For sqrt...
#include "constants.h"
#include "atoms.h"
#include "potential.h"
/* Quadrupole class implementation {{{ */
void Quadrupole::forces(Atoms *atoms, double *acc) {
    int n=atoms->n();
    double *pos=atoms->pos();
    double coeff=(-1.*h/mp)*_bp*atoms->chi()/atoms->m();
    for(int i=0;i<n;i++) {
        int ii=3*i;
        double x=pos[ii];
        double y=pos[ii+1];
        double z=pos[ii+2];
        double r2=x*x+y*y+4*z*z;
        double r=sqrt(r2);
        double invr=coeff/r;
        acc[ii]=x*invr;
        acc[ii+1]=y*invr;
        acc[ii+2]=4*z*invr-_g;
    }
    return;
}
void Quadrupole::ePot(Atoms *atoms) {
    int n=atoms->n();
    double *pos=atoms->pos();
    double epot=0;
    double epotg=0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        double x=pos[ii];
        double y=pos[ii+1];
        double z=pos[ii+2];
        double r2=x*x+y*y+4*z*z;
        double r=sqrt(r2);
        epot+=r;
        epotg+=z;
    }
    epot*=_bp*atoms->chi();
    epotg*=_g*atoms->m()*(mp/h);
    atoms->ePot()=(epot+epotg)/n;
    return;
}
/* }}} */
/* Harmonic class implementation {{{ */
Harmonic::Harmonic(const double ox, const double oy, const double oz) : 
    Potential() {
    _ox=2*pi*ox;
    _ox*=_ox;
    _oy=2*pi*oy;
    _oy*=_oy;
    _oz=2*pi*oz;
    _oz*=_oz;
}
void Harmonic::forces(Atoms *atoms, double *acc) {
    int n=atoms->n();
    double *pos=atoms->pos();
    for(int i=0;i<n;i++) {
        int ii=3*i;
        acc[ii]=-_ox*pos[ii];
        acc[ii+1]=-_ox*pos[ii+1];
        acc[ii+2]=-_ox*pos[ii+2]-_g;
    }
    return;
}
void Harmonic::ePot(Atoms *atoms) {
    int n=atoms->n();
    double *pos=atoms->pos();
    double epot=0;
    double epotg=0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        double x=pos[ii];
        double y=pos[ii+1];
        double z=pos[ii+2];
        epot+=_ox*x*x+_oy*y*y+_oz*z*z;
        epotg+=z;
    }
    epot*=0.5;
    epotg*=_g;
    atoms->ePot()=(epot+epotg)*(mp/h)*atoms->m()/(double)n;
    return;
}
/* }}} */
/* potential.cpp */
