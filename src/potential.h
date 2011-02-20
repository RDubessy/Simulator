/* Copyright (C) 2010 Romain Dubessy */
#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "common.h"
/*!\brief Abstract class that represents an external potential. */
class Potential {
    public:
        /*!\brief Default constructor */
        Potential(double g=9.81) { _g=g; };
        /*!\brief Constructor. */
        Potential(ConfigMap &);
        /*!\brief Computes the forces on the atoms, stored in the acc array. */
        virtual void forces(Atoms *, double *) =0;
        /*!\brief Computes the potential energy of the atoms. */
        virtual void ePot(Atoms *) =0;
        /*!\brief Potential induced losses on atoms. */
        virtual void losses(Atoms *) =0;
    protected:
        double _g;              //!<\brief Gravity [m/s^2].
};
/*!\brief Represents a quadrupole potential of finite depth. */
class Quadrupole : public Potential {
    public:
        /*!\brief Default constructor. */
        Quadrupole(double=2.4525e4, double=1e20);
        /*!\brief Constructor. */
        Quadrupole(ConfigMap &);
        void forces(Atoms *, double *);
        void ePot(Atoms *);
        void losses(Atoms *);
    private:
        double _bp;             //!<\brief Quadrupole gradient [Gauss/m].
        double _U;              //!<\brief Trap depth (given by RF) [Hz].
};
/*!\brief Represents a harmonic potential. */
class Harmonic : public Potential {
    public:
        /*!\brief Default constructor. */
        Harmonic(double=1, double=1, double=1);
        /*!\brief Constructor. */
        Harmonic(ConfigMap &);
        void forces(Atoms *, double *);
        void ePot(Atoms *);
        void losses(Atoms *) {};
    private:
        double _ox;             //!<\brief Pulsation squared [Rad^2/s^2].
        double _oy;             //!<\brief Pulsation squared [Rad^2/s^2].
        double _oz;             //!<\brief Pulsation squared [Rad^2/s^2].
};
#endif //POTENTIAL_H
/* potential.h */
