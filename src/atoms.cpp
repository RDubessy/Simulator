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
#include "coltree.h"
#include "atoms.h"
/* Atoms: {{{ */
Atoms::Atoms(const int n, const int m, const double chi) {
    _n=n;
    _m=m;
    _chi=chi;
    _Gvac=1.;
    _n0=0.;
    _sigma=0.;
    _nc=0;
    _pos=_vel=0;
    _ePot=_eKin=0;
    if(_n>0)
        initCloud(5e-4,5e-4);
}
Atoms::Atoms(ConfigMap &config) {
    _n=getConfig(config,"Atoms::n",1);
    _m=getConfig(config,"Atoms::m",83.);
    _chi=getConfig(config,"Atoms::chi",0.7e6);
    _Gvac=1.0/getConfig(config,"Atoms::lifetime",120.);
    _sigma=getConfig(config,"Atoms::sigma",7e-16);
    _n0=0.;
    _nc=0;
    _pos=_vel=0;
    _ePot=_eKin=0;
    double size=getConfig(config,"Atoms::size",5e-4);
    double T=getConfig(config,"Atoms::T",5e-4);
    if(_n>0)
        initCloud(T,size);
    double dt=getConfig(config,"Integrator::dt",1e-5);
    collisions(dt);
}
/* }}} */
/* ~Atoms: {{{ */
Atoms::~Atoms(void) {
    if(_n>0) {
        delete[] _pos;
        delete[] _vel;
    }
    _n=0;
    _pos=_vel=0;
}
/* }}} */
/* initCloud: {{{ */
void Atoms::initCloud(double T, double r) {
    _pos=new double[3*_n];
    _vel=new double[3*_n];
    memset(_pos,0,3*_n*sizeof(double));
    memset(_vel,0,3*_n*sizeof(double));
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
    double v2=0;
    for(int i=0;i<3*_n;i++) {
        double v=_vel[i];
        v2+=v*v;
    }
    _eKin=(0.5*mp/h)*_m*v2/(double)_n;
}
/* }}} */
/* lifetime: {{{ */
void Atoms::lifetime(double dt) {
    int crit=(int)(dt*_Gvac*RAND_MAX);
    int n=_n;
    for(int i=0;i<_n;i++) {
        if(rand()<crit)
            n--;
    }
    _n=n;
}
/* }}} */
/* collisions: {{{ */
void Atoms::collisions(double dt) {
    CollisionTree tree; 
    _n0=tree.init(this);
    _nc+=tree.compute(this,dt);
}
/* }}} */
/* operator<<: {{{ */
ostream &operator<<(ostream &os, const Atoms &atoms) {
    double x=0;
    double y=0;
    double z=0;
    double x2=0;
    double y2=0;
    double z2=0;
    for(int i=0;i<atoms._n;i++) {
        int ii=3*i;
        x+=atoms._pos[ii];
        x2+=atoms._pos[ii]*atoms._pos[ii];
        y+=atoms._pos[ii+1];
        y2+=atoms._pos[ii+1]*atoms._pos[ii+1];
        z+=atoms._pos[ii+2];
        z2+=atoms._pos[ii+2]*atoms._pos[ii+2];
    }
    double norm=1./atoms._n;
    x*=norm;
    y*=norm;
    z*=norm;
    x2*=norm;
    y2*=norm;
    z2*=norm;
    x2-=x*x;
    y2-=y*y;
    z2-=z*z;
    os << x << " " << y << " " << z << " " << x2 << " " << y2 << " " << z2;
    return os;
}
/* }}} */
/* atoms.cpp */
