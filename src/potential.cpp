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
/* Potential class implementation {{{ */
Potential::Potential(ConfigMap &config) {
    _g=getConfig(config,"Potential::gravity",9.81);
}
/* }}} */
/* Quadrupole class implementation {{{ */
/* Quadrupole: {{{ */
Quadrupole::Quadrupole(double bp, double U) : Potential() {
    _bp=bp;
    _U=U;
}
Quadrupole::Quadrupole(ConfigMap &config) : Potential(config) {
    _bp=getConfig(config,"Potential::gradB",6.7e3);
    _U=getConfig(config,"Potential::depth",1e7);
}
/* }}} */
/* forces: {{{ */
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
}
/* }}} */
/* ePot: {{{ */
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
}
/* }}} */
/* losses: {{{ */
void Quadrupole::losses(Atoms *atoms) {
    int n=atoms->n();
    double *pos=atoms->pos();
    double *vel=atoms->vel();
    double crit=_U/(atoms->chi()*_bp);  //RF evaporation criteria.
    crit*=crit;
    double majorana=atoms->chi()*_bp;   //Majorana losses.
    int i=0;
    while(i<n) {
        int ii=3*i;
        double x=pos[ii];
        double y=pos[ii+1];
        double z=pos[ii+2];
        double r2=x*x+y*y+4*z*z;
        if(r2>=crit) {
            n--;
            int nn=3*n;
            pos[ii]=pos[nn];
            pos[ii+1]=pos[nn+1];
            pos[ii+2]=pos[nn+2];
            vel[ii]=vel[nn];
            vel[ii+1]=vel[nn+1];
            vel[ii+2]=vel[nn+2];
        } else {
            double vx=vel[ii];
            double vy=vel[ii+1];
            double vz=vel[ii+2];
            double v2=vx*vx+vy*vy+vz*vz;
            double vinvr2=sqrt(v2)/r2;
            if(vinvr2>majorana) {
                n--;
                int nn=3*n;
                pos[ii]=pos[nn];
                pos[ii+1]=pos[nn+1];
                pos[ii+2]=pos[nn+2];
                vel[ii]=vel[nn];
                vel[ii+1]=vel[nn+1];
                vel[ii+2]=vel[nn+2];
            } else
                i++;
        }
    }
    atoms->n()=n;
}
/* }}} */
/* }}} */
/* Harmonic class implementation {{{ */
/* Harmonic: {{{ */
Harmonic::Harmonic(double ox, double oy, double oz) : Potential() {
    _ox=2*pi*ox;
    _ox*=_ox;
    _oy=2*pi*oy;
    _oy*=_oy;
    _oz=2*pi*oz;
    _oz*=_oz;
}
Harmonic::Harmonic(ConfigMap &config) : Potential(config) {
    _ox=getConfig(config,"Potential::nu_x",100);
    _ox*=2*pi*_ox;
    _oy=getConfig(config,"Potential::nu_y",100);
    _oy*=2*pi*_oy;
    _oz=getConfig(config,"Potential::nu_z",100);
    _oz*=2*pi*_oz;
}
/* }}} */
/* forces: {{{ */
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
/* }}} */
/* ePot: {{{ */
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
/* }}} */
/* }}} */
/* potential.cpp */
