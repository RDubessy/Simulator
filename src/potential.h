/* Copyright (C) 2010 Romain Dubessy */
#ifndef POTENTIAL_H
#define POTENTIAL_H
/*!\brief Abstract class that represents an external potential. */
class Potential {
    public:
        /*!\brief Constructor. */
        Potential(void) { _g=9.81; };
        /*!\brief Computes the forces on the atoms, stored in the acc array. */
        virtual void forces(Atoms *atoms, double *acc) =0;
        /*!\brief Computes the potential energy of the atoms. */
        virtual void ePot(Atoms *atoms) =0;
    protected:
        double _g;              //!<\brief Gravity.
};
/*!\brief Represents a quadrupole potential. */
class Quadrupole : public Potential {
    public:
        /*!\brief Constructor. */
        Quadrupole(const double bp=6.7e3) : Potential() { _bp=bp; };
        void forces(Atoms *, double *);
        void ePot(Atoms *);
    private:
        double _bp;             //!<\brief Quadrupole gradient [Gauss/m].
};
/*!\brief Represents a harmonic potential. */
class Harmonic : public Potential {
    public:
        /*!\brief Constructor. */
        Harmonic(const double=1, const double=1, const double=1);
        void forces(Atoms *, double *);
        void ePot(Atoms *);
    private:
        double _ox;
        double _oy;
        double _oz;
};
#endif
/* potential.h */
