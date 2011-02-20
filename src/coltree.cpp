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
#include "coltree.h"
/* CollisionTree: {{{ */
CollisionTree::CollisionTree(void) {
    _center[0]=_center[1]=_center[2]=0;
    _size=1.0;
    _n=0;
    _i=-1;
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
void CollisionTree::init(double *pos, int n) {
    for(int i=0;i<n;i++) {
        CollisionTree* tmp=this;
        while(tmp._size!=0) {
            if(tmp._n==0) {     //First case: empty node, insert the atom.
                tmp._n++;
                tmp._i=i;
                break;
            } else {            //Second case: occupied node.
                if(tmp._n==1) {
                    /* Create the children {{{ */
                    tmp._child=new CollisionTree[8];
                    double halfsize=_size/2;
                    for(int j=0;j<8;j++)
                        tmp._child[j]._size=halfsize;
                    /* Update the center coordinates {{{ */
                    halfsize/=2;
                    tmp._child[0]._center[0]=tmp._center[0]-halfsize;
                    tmp._child[0]._center[1]=tmp._center[1]-halfsize;
                    tmp._child[0]._center[2]=tmp._center[2]-halfsize;
                    tmp._child[1]._center[0]=tmp._center[0]-halfsize;
                    tmp._child[1]._center[1]=tmp._center[1]-halfsize;
                    tmp._child[1]._center[2]=tmp._center[2]+halfsize;
                    tmp._child[2]._center[0]=tmp._center[0]-halfsize;
                    tmp._child[2]._center[1]=tmp._center[1]+halfsize;
                    tmp._child[2]._center[2]=tmp._center[2]-halfsize;
                    tmp._child[3]._center[0]=tmp._center[0]-halfsize;
                    tmp._child[3]._center[1]=tmp._center[1]+halfsize;
                    tmp._child[3]._center[2]=tmp._center[2]+halfsize;
                    tmp._child[4]._center[0]=tmp._center[0]+halfsize;
                    tmp._child[4]._center[1]=tmp._center[1]-halfsize;
                    tmp._child[4]._center[2]=tmp._center[2]-halfsize;
                    tmp._child[5]._center[0]=tmp._center[0]+halfsize;
                    tmp._child[5]._center[1]=tmp._center[1]+halfsize;
                    tmp._child[5]._center[2]=tmp._center[2]-halfsize;
                    tmp._child[6]._center[0]=tmp._center[0]+halfsize;
                    tmp._child[6]._center[1]=tmp._center[1]+halfsize;
                    tmp._child[6]._center[2]=tmp._center[2]-halfsize;
                    tmp._child[7]._center[0]=tmp._center[0]+halfsize;
                    tmp._child[7]._center[1]=tmp._center[1]+halfsize;
                    tmp._child[7]._center[2]=tmp._center[2]+halfsize;
                    /* }}} */
                    /* Update the pointers {{{ */
                    tmp._next=&(_child[0]);
                    for(int j=0;j<7;j++)
                        _child[j]._next=_child[j]._skip=&(_child[j+1]);
                    _child[7]._next=_child[7]._skip=tmp._skip;
                    /* }}} */
                    /* Move the existing atom {{{ */
                    int ii=tmp._i;
                    CollisionTree *aux=tmp;
                    if(pos[ii]<tmp._center[0]) {
                        if(pos[ii+1]<tmp._center[1]) {
                            if(pos[ii+2]<tmp._center[2])
                                aux=&(tmp._child[0]);
                            else
                                aux=&(tmp._child[1]);
                        } else {
                            if(pos[ii+2]<tmp._center[2])
                                aux=&(tmp._child[2]);
                            else
                                aux=&(tmp._child[3]);
                        }
                    } else {
                        if(pos[ii+1]<tmp._center[1]) {
                            if(pos[ii+2]<tmp._center[2])
                                aux=&(tmp._child[4]);
                            else
                                aux=&(tmp._child[5]);
                        } else {
                            if(pos[ii+2]<tmp._center[2])
                                aux=&(tmp._child[6]);
                            else
                                aux=&(tmp._child[7]);
                        }
                    }
                    aux._n++;
                    aux._i=ii;
                    tmp._i=-1;
                    /* }}} */
                }
                tmp._n++;       //Increase the node weight.
                /* Insert the atom into its cell {{{ */
                if(pos[i]<tmp._center[0]) {
                    if(pos[i+1]<tmp._center[1]) {
                        if(pos[i+2]<tmp._center[2])
                            tmp=&(tmp._child[0]);
                        else
                            tmp=&(tmp._child[1]);
                    } else {
                        if(pos[i+2]<tmp._center[2])
                            tmp=&(tmp._child[2]);
                        else
                            tmp=&(tmp._child[3]);
                    }
                } else {
                    if(pos[i+1]<tmp._center[1]) {
                        if(pos[i+2]<tmp._center[2])
                            tmp=&(tmp._child[4]);
                        else
                            tmp=&(tmp._child[5]);
                    } else {
                        if(pos[i+2]<tmp._center[2])
                            tmp=&(tmp._child[6]);
                        else
                            tmp=&(tmp._child[7]);
                    }
                }
                /* }}} */
                /* }}} */
            }
        }
    }
}
/* }}} */
/* coltree.cpp */
