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
#include <cstring>              //For memset, memcpy
#include <cmath>                //For sqrt...
#include <ostream>              //For ostream.
#include <stdlib.h>             //For rand.
#include "constants.h"
#include "atoms.h"
/* Default constructor {{{ */
Atoms::Atoms(const int n, const int m, const double chi) {
    _n=n;
    _m=m;
    _chi=chi;
    _pos=_vel=0;
    _ePot=_eKin=0;
    if(_n>0) {
        _pos=new double[3*_n];
        _vel=new double[3*_n];
        memset(_pos,0,3*_n*sizeof(double));
        memset(_vel,0,3*_n*sizeof(double));
        initCloud(5e-4,5e-4);
        double v2=0;
        for(int i=0;i<3*_n;i++) {
            double v=_vel[i];
            v2+=v*v;
        }
        _eKin=(0.5*mp/h)*_m*v2/(double)_n;
    }
}
/* }}} */
/* Destructor {{{ */
Atoms::~Atoms(void) {
    if(_n>0) {
        delete[] _pos;
        delete[] _vel;
    }
    _n=0;
    _pos=_vel=0;
}
/* }}} */
/* Initialization method {{{ */
void Atoms::initCloud(double T, double r) {
    double v=sqrt(kB*T/(_m*mp));
    for(int i=0;i<3*_n;i+=3) {
        for(int d=0;d<3;d++) {
            //Box-Muller method.
            double r1=sqrt(-2*log((double)rand()/RAND_MAX));
            double r2=2*pi*((double)rand()/RAND_MAX);
            _pos[i+d]=r*r1*cos(r2);
            _vel[i+d]=v*r1*sin(r2);
        }
    }
}
/* }}} */
/* Conversion to ostream {{{ */
std::ostream &operator<<(std::ostream &os, const Atoms &atoms) {
    double x=0;
    double y=0;
    double z=0;
    for(int i=0;i<atoms._n;i++) {
        int ii=3*i;
        x+=atoms._pos[ii];
        y+=atoms._pos[ii+1];
        z+=atoms._pos[ii+2];
    }
    double norm=1./atoms._n;
    x*=norm;
    y*=norm;
    z*=norm;
    os << x << " " << y << " " << z;
    return os;
}
/* }}} */
/* atoms.cpp */
