/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  ParallelGravityKernel.h : the in-plane gravity term (the flow down a
*                            tilted substrate)
*
*  The skeleton lives in BodyForceKernel<T>. All that is specific to the
*  in-plane gravity is that its driving field is a constant direction that
*  does not depend on position.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _PARALLEL_GRAVITY_KERNEL_INCLUDED_
#define _PARALLEL_GRAVITY_KERNEL_INCLUDED_

#include "BodyForceKernel.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  In-plane gravity. It takes effect when the substrate is not level.
//  The driving field is constant, so only its direction (its sign) is used
//  to choose the upwind side; its magnitude goes into the per-axis
//  coefficients Pgx/3 and Pgy/3.
//--------------------------------------------------------------------------
template<class T>
class ParallelGravityKernel : public BodyForceKernel<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::Boundary  Boundary;
	typedef typename T::VelField  VelField;
	typedef PhysicsFace::tEdge    tEdge;

	/** @brief Start from a level substrate until Build() reads the input. */
	ParallelGravityKernel() : BodyForceKernel<T>()
	{
		for(int a = 0; a < T::nAxis; ++a) { _g[a] = 0; _Mag[a] = 0; }
	}
	virtual ~ParallelGravityKernel(){}

	/**
	 * @brief Read the per-axis coefficient Pg/3 and its direction from the input.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         the grid
	 * @param nonDim    the set of non-dimensional numbers
	 * @param bc        the boundary condition
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.HorizontalGravity.Enable;

		for(int a = 0; a < T::nAxis; ++a)
		{
			const double pg = (a == 0) ? this->_NonDim->Pgx() : this->_NonDim->Pgy();
			_Mag[a] = pg / 3.0;
			_g[a]   = Math::Sign(pg);
		}
		this->_Magnitude = _Mag[0];

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		this->SetLinearCoef();
	}

protected:
	double _g[T::nAxis];		// the sign of the direction of gravity (used only to choose the upwind side)
	double _Mag[T::nAxis];		// the per-axis coefficient Pg/3

	/**
	 * @brief The driving field used to choose the upwind side. It does not
	 *        depend on position, so only the direction is returned.
	 */
	virtual double Drive(Site &s, tEdge te) const
	{
		(void)s;
		const int d = T::DirIndex(te);
		return d < 0 ? 0.0 : _g[d/2];
	}

	/**
	 * @brief The version measured at the face coordinate. The field is
	 *        constant, so this is the same as Drive.
	 */
	virtual double DriveFace(Site &s, tEdge te) const { return Drive(s, te); }

	//  Weight is left at the base class default (1). The magnitude is held
	//  on the Mag side instead.

	/** @brief The per-axis coefficient Pg/3. */
	virtual double Mag(int axis) const { return _Mag[axis]; }

	/** @brief Where the in-plane gravity part of the velocity is written. */
	virtual VelField& VelocityField(Field &field) const
	{ return field.UU().U_HorizontalGravity(); }

	/**
	 * @brief The sign used when writing the velocity. CalVelocity forms the
	 *        flux with CV = -(interpolation coefficients), so the sum comes
	 *        out with the sign reversed; -1 puts it back, and the velocity
	 *        points along Acceleration.
	 */
	virtual double VelocitySign() const { return -1.0; }

private:
	ParallelGravityKernel(const ParallelGravityKernel&);
	void operator=(const ParallelGravityKernel&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_PARALLEL_GRAVITY_KERNEL_INCLUDED_
