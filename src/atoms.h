/* Copyright (C) 2010 Romain Dubessy */
#ifndef ATOMS_H
#define ATOMS_H
#include <iosfwd>               //For ostream forward declaration.
/*!\brief Represents a cloud of atoms. */
class Atoms {
    public:
        /*!\brief Constructor. */
        Atoms(const int=0, const int=87, const double=1.4e6);
        /*!\brief Destructor. */
        ~Atoms(void);
        /*!\brief Initialization method. */
        void initCloud(double, double);
        /* Access member methods {{{ */
        /*!\brief Return the number of atoms. */
        int n(void) const { return _n; };
        /*!\brief Return a pointer to the array of atom's positions. */
        double *pos(void) { return _pos; };
        /*!\brief Return a pointer to the array of atom's velocities. */
        double *vel(void) { return _vel; };
        /*!\brief Return the atom mass (atomic units). */
        double m(void) const { return _m; };
        /*!\brief Return the atom susceptibility (Hz/Gauss). */
        double chi(void) const { return _chi; };
        double &ePot(void) { return _ePot; };
        double ePot(void) const { return _ePot; };
        double &eKin(void) { return _eKin; };
        double eKin(void) const { return _eKin; };
        /* }}} */
        /*!\brief Conversion to ostream operator. */
        friend std::ostream &operator<<(std::ostream &, const Atoms &);
    private:
        double _eKin;           //!<\brief Mean kinetic energy [J].
        double _ePot;           //!<\brief Mean potential energy [J].
        double _m;              //!<\brief Atom's mass [mp].
        double _chi;            //!<\brief Magnetic susceptibility [Hz/Gauss]
        double *_pos;           //!<\brief Positions (double[3*_n]) [m].
        double *_vel;           //!<\brief Velocities (double[3*_n]) [m/s].
        int _n;                 //!<\brief Atom's number.
};
#endif
/* atoms.h */
