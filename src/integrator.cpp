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
#include "integrator.h"
using std::cout;
using std::cerr;
using std::endl;
/* Class Integrator implementation {{{ */
/* Integrator: {{{ */
Integrator::Integrator(void) {
    _t=_dt=_dtOut=_dtOut=0;
    _seed=(int)time(0);
    srand(_seed);
    _atoms=0;
    _potential=0;
    _run=true;
}
Integrator::Integrator(ConfigMap &config) {
    _t=getConfig(config,"Integrator::t",1.);
    _dt=getConfig(config,"Integrator::dt",_t*1e-3);
    _dtEvent=10*_dt;
    _dtOut=getConfig(config,"Integrator::dtOut",_dt*10.);
    _seed=getConfig(config,"Integrator::seed",(int)time(0));
    srand(_seed);
    _atoms=new Atoms(config);
    _potential=0;
    string type=getConfig(config,"Potential::type","Quadrupole");
    if(type=="Quadrupole")
        _potential=new Quadrupole(config);
    else if(type=="Harmonic")
        _potential=new Harmonic(config);
    _run=true;
}
/* }}} */
/* ~Integrator: {{{ */
Integrator::~Integrator(void) {
    if(_atoms!=0)
        delete _atoms;
    if(_potential!=0)
        delete _potential;
}
/* }}} */
/* evolve: {{{ */
int Integrator::evolve(void) {
    double t=0.;
    double tOut=0.;
    double tEvent=0.;
    cout << "t <x> <y> <z> <x2> <y2> <z2> <Ekin> <Epot> n n0 Gc\n";
    while(_run) {
        if(t>=tEvent) {
            events();
            tEvent+=_dtEvent;
        }
        if(t>=tOut) {
            _potential->ePot(_atoms);
            measure(t);
            tOut+=_dtOut;
        }
        if(t>=_t) {
            _run=false;
            break;
        }
        doSteps();
        t+=_dt;
    }
    return 0;
}
/* }}} */
/* measure: {{{ */
void Integrator::measure(double t) {
    cout << t << " " << *_atoms << " " << _atoms->eKin() << " "
        << _atoms->ePot() << " " << _atoms->n() << " ";
    double Gc=_dtOut*(_atoms->n());
    Gc=1./Gc*(_atoms->nc());
    _atoms->nc()=0;
    cout << _atoms->n0() << " " << Gc;
    cout << "\n";
}
/* }}} */
/* events: {{{ */
void Integrator::events(void) {
    //_atoms->lifetime(_dtEvent);
    _atoms->collisions(_dtEvent);
    double nc=((double)_atoms->nc()/_atoms->n());
    if(nc>0.5)
        _dtEvent/=10;
    else if(nc<0.01)
        _dtEvent*=10;
    if(_dtEvent<_dt)
        _dtEvent=_dt;
    else if(_dtEvent>_dtOut)
        _dtEvent=_dtOut;
    //_potential->losses(_atoms);
}
/* }}} */
/* }}} */
/* Class RK2 implementation {{{ */
/* RK2: {{{ */
RK2::RK2(void) : Integrator() {
    _acc=0;
    _oldpos=0;
    _oldvel=0;
    _n=_atoms->n();
    init();
}
RK2::RK2(ConfigMap &config) : Integrator(config) {
    _acc=0;
    _oldpos=0;
    _oldvel=0;
    _n=_atoms->n();
    init();
}
/* }}} */
/* ~RK2: {{{ */
RK2::~RK2(void) {
    if(_acc!=0)
        delete[] _acc;
    if(_oldpos!=0)
        delete[] _oldpos;
    if(_oldvel!=0)
        delete[] _oldvel;
}
/* }}} */
/* doSteps: {{{ */
void RK2::doSteps(void) {
    int n=_atoms->n();
    if(n==0) {
        _run=false;
        return;
    }
    if(_n!=n) {
        _n=n;
        //init();
    }
    double *pos=_atoms->pos();
    double *vel=_atoms->vel();
    memcpy(_oldpos,pos,3*n*sizeof(double));
    memcpy(_oldvel,vel,3*n*sizeof(double));
    double dt=_dt*0.5;
    //First step.
    _potential->forces(_atoms,_acc);
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]+=dt*vel[ii+d];
            vel[ii+d]+=dt*_acc[ii+d];
        }
    }
    dt=_dt;
    //Second step.
    _potential->forces(_atoms,_acc);
    double v2=0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]=_oldpos[ii+d]+dt*vel[ii+d];
            double v=_oldvel[ii+d]+dt*_acc[ii+d];
            vel[ii+d]=v;
            v2+=v*v;
        }
    }
    _atoms->eKin()=_atoms->m()*(0.5*mp/h)*v2/(double)n;
    return;
}
/* }}} */
/* init: {{{ */
void RK2::init(void) {
    if(_acc!=0)
        delete[] _acc;
    if(_oldpos!=0)
        delete[] _oldpos;
    if(_oldvel!=0)
        delete[] _oldvel;
    _acc=new double[3*_n];
    _oldpos=new double[3*_n];
    _oldvel=new double[3*_n];
    memset(_acc,0,3*_n*sizeof(double));
    memset(_oldpos,0,3*_n*sizeof(double));
    memset(_oldvel,0,3*_n*sizeof(double));
    return;
}
/* }}} */
/* }}} */
/* Class RK4 implementation {{{ */
/* RK4: {{{ */
RK4::RK4(void) : Integrator() {
    _oldpos=0;
    _oldvel=0;
    _pos=0;
    _vel=0;
    _acc=0;
    _n=_atoms->n();
    init();
}
RK4::RK4(ConfigMap &config) : Integrator(config) {
    _oldpos=0;
    _oldvel=0;
    _pos=0;
    _vel=0;
    _acc=0;
    _n=_atoms->n();
    init();
}
/* }}} */
/* ~RK4: {{{ */
RK4::~RK4(void) {
    if(_oldpos!=0)
        delete[] _oldpos;
    if(_oldvel!=0)
        delete[] _oldvel;
    if(_pos!=0)
        delete[] _pos;
    if(_vel!=0)
        delete[] _vel;
    if(_acc!=0)
        delete[] _acc;
}
/* }}} */
/* doSteps: {{{ */
void RK4::doSteps(void) {
    int n=_atoms->n();
    if(n==0) {
        _run=false;
        return;
    }
    if(_n!=n) {
        _n=n;
        //init();
    }
    double *pos=_atoms->pos();
    double *vel=_atoms->vel();
    memcpy(_oldpos,pos,3*n*sizeof(double));
    memcpy(_pos,pos,3*n*sizeof(double));
    memcpy(_oldvel,vel,3*n*sizeof(double));
    memcpy(_vel,vel,3*n*sizeof(double));
    double dt=_dt*0.5;
    double dt0=_dt/6.;
    //First step
    _potential->forces(_atoms,_acc);
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]+=dt*vel[ii+d];
            _pos[ii+d]+=dt0*vel[ii+d];
            vel[ii+d]+=dt*_acc[ii+d];
            _vel[ii+d]+=dt0*_acc[ii+d];
        }
    }
    //Second step
    dt0=_dt/3.;
    _potential->forces(_atoms,_acc);
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]=_oldpos[ii+d]+dt*vel[ii+d];
            _pos[ii+d]+=dt0*vel[ii+d];
            vel[ii+d]=_oldvel[ii+d]+dt*_acc[ii+d];
            _vel[ii+d]+=dt0*_acc[ii+d];
        }
    }
    dt=_dt;
    //Third step
    _potential->forces(_atoms,_acc);
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]=_oldpos[ii+d]+dt*vel[ii+d];
            _pos[ii+d]+=dt0*vel[ii+d];
            vel[ii+d]=_oldvel[ii+d]+dt*_acc[ii+d];
            _vel[ii+d]+=dt0*_acc[ii+d];
        }
    }
    dt0=_dt/6.;
    //Fourth step
    _potential->forces(_atoms,_acc);
    double v2=0;
    for(int i=0;i<n;i++) {
        int ii=3*i;
        for(int d=0;d<3;d++) {
            pos[ii+d]=_pos[ii+d]+dt0*vel[ii+d];
            double v=_vel[ii+d]+dt0*_acc[ii+d];
            vel[ii+d]=v;
            v2+=v*v;
        }
    }
    _atoms->eKin()=_atoms->m()*(0.5*mp/h)*v2/(double)n;
}
/* }}} */
/* init: {{{ */
void RK4::init(void) {
    if(_oldpos!=0)
        delete[] _oldpos;
    if(_oldvel!=0)
        delete[] _oldvel;
    if(_pos!=0)
        delete[] _pos;
    if(_vel!=0)
        delete[] _vel;
    if(_acc!=0)
        delete[] _acc;
    _oldpos=new double[3*_n];
    _oldvel=new double[3*_n];
    _pos=new double[3*_n];
    _vel=new double[3*_n];
    _acc=new double[3*_n];
    memset(_oldpos,0,3*_n*sizeof(double));
    memset(_oldvel,0,3*_n*sizeof(double));
    memset(_pos,0,3*_n*sizeof(double));
    memset(_vel,0,3*_n*sizeof(double));
    memset(_acc,0,3*_n*sizeof(double));
    return;
}
/* }}} */
/* }}} */
/* initIntegrator: {{{ */
Integrator *initIntegrator(ConfigMap &config) {
    string type=getConfig(config,"Integrator::type","RungeKutta2");
    if(type=="RungeKutta2")
        return new RK2(config);
    else if(type=="RungeKutta4")
        return new RK4(config);
    return 0;
}
/* }}} */
/* integrator.cpp */
