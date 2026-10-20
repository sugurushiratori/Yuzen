/**************************************************************************
*  Yuzen - thin liquid film flow solver (lubrication approximation)
*
*  CentrifugalKernel.h : the centrifugal force term (spin coating)
*
*  The skeleton lives in BodyForceKernel<T>. All that is specific to the
*  centrifugal force is that its driving field is the distance from the
*  axis of rotation.
*
*  Copyright 2017-2026 Suguru Shiratori, Tokyo City University
*  SPDX-License-Identifier: Apache-2.0
**************************************************************************/

#ifndef _CENTRIFUGAL_KERNEL_INCLUDED_
#define _CENTRIFUGAL_KERNEL_INCLUDED_

#include "BodyForceKernel.h"

#ifdef _DEF_YUZEN_
namespace yuzen{
#endif

//--------------------------------------------------------------------------
//  The centrifugal force: the effect that flings the liquid film outwards
//  when the substrate rotates, and the main driving force of spin coating.
//  The driving field is the distance r from the axis of rotation, and the
//  coefficient is Sp/3.
//--------------------------------------------------------------------------
template<class T>
class CentrifugalKernel : public BodyForceKernel<T>
{
public:
	typedef typename T::Grid      Grid;
	typedef typename T::Site      Site;
	typedef typename T::Field     Field;
	typedef typename T::Boundary  Boundary;
	typedef typename T::VelField  VelField;
	typedef PhysicsFace::tEdge    tEdge;

	/** @brief Put the rotation center at the origin until Build() reads it. */
	CentrifugalKernel() : BodyForceKernel<T>()
	{
		for(int a = 0; a < T::nAxis; ++a) _Cnt[a] = 0;
	}
	virtual ~CentrifugalKernel(){}

	/**
	 * @brief Read the rotation center and the coefficient Sp/3 from the input.
	 * @param inp_phys  the Physics block of the input file
	 * @param g         the grid
	 * @param nonDim    the set of non-dimensional numbers
	 * @param bc        the boundary condition
	 *
	 * The rotation center is given as a real length in the input, so it is
	 * non-dimensionalized the same way as the grid (divided by ExtentScale)
	 * before being stored.
	 */
	virtual void Build(Input::I_Physics &inp_phys, Grid *g, NonDimensionalNumbers *nonDim, Boundary *bc)
	{
		this->_G = g;
		this->_NonDim = nonDim;
		this->_boundary = bc;

		this->_enabled = inp_phys.Centrifugal.Enable;

		for(int a = 0; a < T::nAxis; ++a)
		{
			_Cnt[a] = inp_phys.Centrifugal.RotationAxis[a];
			_Cnt[a] /= this->_NonDim->ExtentScale();
		}

		this->_Magnitude = this->_NonDim->Sp() / 3.0;

		T::SetExtent(this->_G, this->_Nx, this->_Ny);
		this->SetLinearCoef();
	}

protected:
	double _Cnt[T::nAxis];		// the center of rotation (non-dimensional)

	/**
	 * @brief The distance of a face from the axis of rotation, as used to
	 *        choose the upwind side.
	 *
	 * It is measured as the mean of the center coordinates of the two cells
	 * that straddle the face.
	 */
	virtual double Drive(Site &s, tEdge te) const
	{
		const int d = T::DirIndex(te);
		if(d < 0) return 0.0;
		return T::CellAvgCoord(this->_G, s, te) - _Cnt[d/2];
	}

	/**
	 * @brief The distance of a face from the axis of rotation, as used by
	 *        the matrix assembly and the velocity computation.
	 *
	 * This one is measured at the face coordinate itself. On a uniform grid
	 * it is the same value as Drive, but the two do not agree down to the
	 * rounding, so they are kept separate.
	 */
	virtual double DriveFace(Site &s, tEdge te) const
	{
		const int d = T::DirIndex(te);
		if(d < 0) return 0.0;
		return T::FaceCoord(this->_G, s, te) - _Cnt[d/2];
	}

	/**
	 * @brief The weight applied to the stencil; for the centrifugal force
	 *        it is the distance itself.
	 */
	virtual double Weight(Site &s, tEdge te) const { return DriveFace(s, te); }

	/**
	 * @brief The coefficient of the given axis; for the centrifugal force
	 *        it is the same in X and in Y.
	 */
	virtual double Mag(int axis) const { (void)axis; return this->_Magnitude; }

	/**
	 * @brief Whether the outflow in a direction that is not solved for is
	 *        accounted for on the diagonal.
	 *
	 * The 1-D case does not solve for the Y direction, so the liquid flung
	 * out that way is added to the diagonal. In 2-D both directions are
	 * solved for, so this is not needed.
	 */
	virtual bool SideOutflow() const { return T::nAxis == 1; }

	/** @brief Where the centrifugal part of the velocity is written. */
	virtual VelField& VelocityField(Field &field) const
	{ return field.UU().U_Centrifugal(); }

	/**
	 * @brief The sign used when writing the velocity. CalVelocity forms the
	 *        flux with CV = -(interpolation coefficients), so the sum comes
	 *        out with the sign reversed, exactly as for the in-plane gravity;
	 *        -1 makes the velocity point outwards, the way the film goes.
	 *        (At +1 the velocity written, and carried into the concentration
	 *        equation, pointed at the axis while the film spun off.)
	 */
	virtual double VelocitySign() const { return -1.0; }

private:
	CentrifugalKernel(const CentrifugalKernel&);
	void operator=(const CentrifugalKernel&);
};

#ifdef _DEF_YUZEN_
}	// END: namespace yuzen
#endif

#endif //_CENTRIFUGAL_KERNEL_INCLUDED_
