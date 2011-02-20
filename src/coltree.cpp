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
#include <stdlib.h>             //For rand.
#include <iostream>
#include "atoms.h"
#include "coltree.h"
/* CollisionTree: {{{ */
CollisionTree::CollisionTree(void) {
    _center[0]=_center[1]=_center[2]=0;
    _size=0.1;
    _n=0;
    _i=_j=-1;
    _child=_next=_skip=0;
}
/* }}} */
/* ~CollisionTree: {{{ */
CollisionTree::~CollisionTree(void) {
    if(_child!=0)
        delete[] _child;
}
/* }}} */
/* init: {{{ */
double CollisionTree::init(Atoms *atoms) {
    int n=atoms->n();
    double *pos=atoms->pos();
    double res=_size;
    for(int i=0;i<n;i++) {
        CollisionTree* tmp=this;
        while(tmp->_size!=0) {
            if(tmp->_n==0) {     //First case: empty node, insert the atom.
                tmp->_n++;
                tmp->_i=i;
                if(tmp->_size<res)
                    res=tmp->_size;
                break;
            } else {            //Second case: occupied node.
                if(tmp->_n==1) {
                    tmp->_j=i;
                    /* Create the children {{{ */
                    tmp->_child=new CollisionTree[8];
                    double halfsize=tmp->_size/2;
                    for(int j=0;j<8;j++)
                        tmp->_child[j]._size=halfsize;
                    /* }}} */
                    /* Update the center coordinates {{{ */
                    halfsize/=2;
                    tmp->_child[0]._center[0]=tmp->_center[0]-halfsize;
                    tmp->_child[0]._center[1]=tmp->_center[1]-halfsize;
                    tmp->_child[0]._center[2]=tmp->_center[2]-halfsize;
                    tmp->_child[1]._center[0]=tmp->_center[0]-halfsize;
                    tmp->_child[1]._center[1]=tmp->_center[1]-halfsize;
                    tmp->_child[1]._center[2]=tmp->_center[2]+halfsize;
                    tmp->_child[2]._center[0]=tmp->_center[0]-halfsize;
                    tmp->_child[2]._center[1]=tmp->_center[1]+halfsize;
                    tmp->_child[2]._center[2]=tmp->_center[2]-halfsize;
                    tmp->_child[3]._center[0]=tmp->_center[0]-halfsize;
                    tmp->_child[3]._center[1]=tmp->_center[1]+halfsize;
                    tmp->_child[3]._center[2]=tmp->_center[2]+halfsize;
                    tmp->_child[4]._center[0]=tmp->_center[0]+halfsize;
                    tmp->_child[4]._center[1]=tmp->_center[1]-halfsize;
                    tmp->_child[4]._center[2]=tmp->_center[2]-halfsize;
                    tmp->_child[5]._center[0]=tmp->_center[0]+halfsize;
                    tmp->_child[5]._center[1]=tmp->_center[1]-halfsize;
                    tmp->_child[5]._center[2]=tmp->_center[2]+halfsize;
                    tmp->_child[6]._center[0]=tmp->_center[0]+halfsize;
                    tmp->_child[6]._center[1]=tmp->_center[1]+halfsize;
                    tmp->_child[6]._center[2]=tmp->_center[2]-halfsize;
                    tmp->_child[7]._center[0]=tmp->_center[0]+halfsize;
                    tmp->_child[7]._center[1]=tmp->_center[1]+halfsize;
                    tmp->_child[7]._center[2]=tmp->_center[2]+halfsize;
                    /* }}} */
                    /* Move the existing atom {{{ */
                    int ii=3*(tmp->_i);
                    CollisionTree *aux=tmp;
                    if(pos[ii]<tmp->_center[0]) {
                        if(pos[ii+1]<tmp->_center[1]) {
                            if(pos[ii+2]<tmp->_center[2])
                                aux=&(tmp->_child[0]);
                            else
                                aux=&(tmp->_child[1]);
                        } else {
                            if(pos[ii+2]<tmp->_center[2])
                                aux=&(tmp->_child[2]);
                            else
                                aux=&(tmp->_child[3]);
                        }
                    } else {
                        if(pos[ii+1]<tmp->_center[1]) {
                            if(pos[ii+2]<tmp->_center[2])
                                aux=&(tmp->_child[4]);
                            else
                                aux=&(tmp->_child[5]);
                        } else {
                            if(pos[ii+2]<tmp->_center[2])
                                aux=&(tmp->_child[6]);
                            else
                                aux=&(tmp->_child[7]);
                        }
                    }
                    aux->_n++;
                    aux->_i=tmp->_i;
                    /* }}} */
                }
                tmp->_n++;       //Increase the node weight.
                /* Insert the atom into its cell {{{ */
                int jj=3*i;
                CollisionTree *aux=tmp;
                if(pos[jj]<tmp->_center[0]) {
                    if(pos[jj+1]<tmp->_center[1]) {
                        if(pos[jj+2]<tmp->_center[2])
                            aux=&(tmp->_child[0]);
                        else
                            aux=&(tmp->_child[1]);
                    } else {
                        if(pos[jj+2]<tmp->_center[2])
                            aux=&(tmp->_child[2]);
                        else
                            aux=&(tmp->_child[3]);
                    }
                } else {
                    if(pos[jj+1]<tmp->_center[1]) {
                        if(pos[jj+2]<tmp->_center[2])
                            aux=&(tmp->_child[4]);
                        else
                            aux=&(tmp->_child[5]);
                    } else {
                        if(pos[jj+2]<tmp->_center[2])
                            aux=&(tmp->_child[6]);
                        else
                            aux=&(tmp->_child[7]);
                    }
                }
                tmp=aux;
                /* }}} */
            }
        }
    }
    updatePointers();
    return 1.0/(res*res*res);
}
/* }}} */
/* updatePointers: {{{ */
void CollisionTree::updatePointers(void) {
    if(_n<2)
        return;
    CollisionTree *tmp=this;
    do {
        if(tmp->_n==1) {
            tmp->_next=tmp->_skip;
        } else {
            int i=0;
            while(tmp->_child[i]._n==0)
                i++;
            tmp->_next=&(tmp->_child[i]);
            for(int j=i+1;j<8;j++) {
                if(tmp->_child[j]._n!=0) {
                    tmp->_child[i]._skip=&(tmp->_child[j]);
                    i=j;
                }
            }
            tmp->_child[i]._skip=tmp->_skip;
        }
        tmp=tmp->_next;
    } while(tmp!=0);
    return;
}
/* }}} */
/* print: {{{ */
void CollisionTree::print(void) {
    for(CollisionTree *tmp=this;tmp!=0;tmp=tmp->_next)
        std::cerr << tmp << "->";
    std::cerr << "0" << std::endl;
}
/* }}} */
/* compute: {{{ */
int CollisionTree::compute(Atoms* atoms, double dt) {
    double *vel=atoms->vel();
    double crit=2*dt*(atoms->sigma());
    int res=0;
    for(CollisionTree *tmp=this;tmp!=0;) {
        if(tmp->_n==2) {        //Collision ?
            //Get the indexes.
            int ii=3*(tmp->_i);
            int jj=3*(tmp->_j);
            //Compute the relative velocity norm.
            double vx=vel[ii]-vel[jj];
            double vy=vel[ii+1]-vel[jj+1];
            double vz=vel[ii+2]-vel[jj+2];
            double v=sqrt(vx*vx+vy*vy+vz*vz);
            //Compute the (local) density.
            double invrn0=tmp->_size;
            invrn0=invrn0*invrn0*invrn0;
            //p<sigma*n0*v*dt ?
            double test=invrn0*((double)rand()/RAND_MAX);
            if(test<crit*v) {
                v/=2;
                //Randomize the velocities.
                vx=(vel[ii]+vel[jj])/2;
                vy=(vel[ii+1]+vel[jj+1])/2;
                vz=(vel[ii+2]+vel[jj+2])/2;
                double ctheta=2*((double)rand()/RAND_MAX)-1;
                double stheta=sqrt(1-ctheta*ctheta);
                double cphi=2*((double)rand()/RAND_MAX)-1;
                double sphi=sqrt(1-cphi*cphi);
                vel[ii]=vx+v*stheta*cphi;
                vel[ii+1]=vy+v*stheta*sphi;
                vel[ii+2]=vz+v*ctheta;
                vel[jj]=vx-v*stheta*cphi;
                vel[jj+1]=vy-v*stheta*sphi;
                vel[jj+2]=vz-v*ctheta;
                res++;
            }
            tmp=tmp->_skip;
        } else
            tmp=tmp->_next;
    }
    return res;
}
/* }}} */
/* coltree.cpp */
