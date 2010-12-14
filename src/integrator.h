/* Copyright (C) 2010 Romain Dubessy */
#ifndef INTEGRATOR_H
#define INTEGRATOR_H
class Atoms;
class Potential;
/*!\brief Abstract base class that represents an integrator. */
class Integrator {
    public:
        /*!\brief Constructor. */
        Integrator(const int=1);
        /*!\brief Destructor. */
        ~Integrator(void);
        /*!\brief Evolution method. */
        int evolve(void);
        /*!\brief Integrator steps method. */
        virtual void doSteps(void) =0;
    protected:
        double _t;              //!<\brief Total time of the simulation.
        double _dt;             //!<\brief Temporal step size.
        double _dtOut;          //!<\brief Measurement step size.
        Atoms *_atoms;          //!<\brief Atoms.
        Potential *_potential;  //!<\brief Potential.
};
/*!\brief 2nd order Runge-Kutta integrator implementation. */
class RK2 : public Integrator {
    public:
        /*!\brief Constructor. */
        RK2(const int n=1) : Integrator(n) {};
        void doSteps(void);
};
Integrator *initIntegrator(ConfigMap &config);
#endif
/* integrator.h */
