/* Copyright (C) 2010 Romain Dubessy */
#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "common.h"             //For ConfigMap.
class Atoms;
class Potential;
/*!\brief Abstract base class that represents an integrator. */
class Integrator {
    public:
        /*!\brief Default Constructor. */
        Integrator(void);
        /*!\brief Constructor */
        Integrator(ConfigMap &);
        /*!\brief Destructor. */
        ~Integrator(void);
        /*!\brief Evolution method. */
        int evolve(void);
        /*!\brief Integrator steps method. */
        virtual void doSteps(void) =0;
        /*!\brief Initialization method. */
        virtual void init(void) =0;
        /*!\brief Events method. */
        void events(void);
        /*!\brief Measuze method. */
        void measure(double);
    protected:
        double _t;              //!<\brief Total time of the simulation.
        double _dt;             //!<\brief Temporal step size.
        double _dtOut;          //!<\brief Measurement step size.
        double _dtEvent;        //!<\brief Event step size.
        Atoms *_atoms;          //!<\brief Atoms.
        Potential *_potential;  //!<\brief Potential.
        int _seed;              //!<\brief Random number generator seed.
        bool _run; 
};
/*!\brief 2nd order Runge-Kutta integrator implementation. */
class RK2 : public Integrator {
    public:
        /*!\brief Default constructor. */
        RK2(void);
        /*!\brief Constructor. */
        RK2(ConfigMap &);
        /*!\brief Destructor. */
        ~RK2(void);
        void doSteps(void);
        void init(void);
    private:
        double *_acc;
        double *_oldpos;
        double *_oldvel;
        int _n;
};
/*!\brief 4th order Runge-Kutta integrator implementation. */
class RK4 : public Integrator {
    public:
        /*!\brief Default constructor. */
        RK4(void);
        /*!\brief Constructor. */
        RK4(ConfigMap &);
        /*!\brief Destructor. */
        ~RK4(void);
        void doSteps(void);
        void init(void);
    private:
        double *_oldpos;
        double *_oldvel;
        double *_pos;
        double *_vel;
        double *_acc;
        int _n;
};
/*!\brief Integrator initialization method. */
Integrator *initIntegrator(ConfigMap &config);
#endif //INTEGRATOR_H
/* integrator.h */
