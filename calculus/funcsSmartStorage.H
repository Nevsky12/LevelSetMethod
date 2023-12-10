#ifndef funcsSmartStorage_H
#define funcsSmartStorage_H

#include "funcs/charFuncs.H"
#include "funcs/heaviside.H"
#include "funcs/delta.H"

namespace Foam
{
namespace levelSet::calculus
{

class funcsSmartStorage
{

	volScalarField charFuncPos_;

	volScalarField charFuncNeg_;

	volScalarField heaviside_;

	volScalarField delta_;

public:

       funcsSmartStorage() = delete;
      ~funcsSmartStorage() = default;

       void updateAndCorrect(volScalarField const &psi, scalar const epsilon) noexcept
       {
	    charFuncPos_ = charFuncPos(psi         )();
	    charFuncNeg_ = charFuncNeg(psi         )();
	    heaviside_   = heaviside  (psi, epsilon)();
	    delta_       = delta      (psi, epsilon)();

	   
	    charFuncPos_.correctBoundaryConditions();
	    charFuncNeg_.correctBoundaryConditions();
	    heaviside_.  correctBoundaryConditions();
	    delta_.      correctBoundaryConditions();
       }

       funcsSmartStorage(volScalarField const &psi, scalar const epsilon) noexcept
       :
       	    charFuncPos_                                          
	    (
	    	IOobject
		(
		   "CharPosfunc",
		   psi.time().timeName(),
		   psi.mesh()
		),
		psi.mesh(),
		dimensionedScalar("Char[+]func", dimless, 0.),
		psi.boundaryField().types()
	    ),
            charFuncNeg_                                           
            (
            	IOobject
                (
                   "CharNegfunc",
                   psi.time().timeName(),
                   psi.mesh()
                ),
                psi.mesh(),
                dimensionedScalar("Char[-]func", dimless, 0.),
		psi.boundaryField().types()
            ),
            heaviside_                                                     
            (
            	IOobject
                (
                   "Heaviside",
                   psi.time().timeName(),
                   psi.mesh()
                ),
                psi.mesh(),
                dimensionedScalar("Heaviside", dimless, 0.),
		psi.boundaryField().types()
            ),
            delta_                                                
            (
            	IOobject
                (
                   "DeltaFunc",
                   psi.time().timeName(),
                   psi.mesh()
                ),
                psi.mesh(),
                dimensionedScalar("DeltaFunc", dimless, 0.),
		psi.boundaryField().types()
            )
       {
	    updateAndCorrect(psi, epsilon);
       }

       volScalarField const & gCharFuncPos() const noexcept { return charFuncPos_; }
       volScalarField const & gcharFuncNeg() const noexcept { return charFuncNeg_; }
       volScalarField const & gHeaviside  () const noexcept { return heaviside_;   }
       volScalarField const & gDelta      () const noexcept { return delta_;       }

};

} // End namespace levelSet::calculus
} // End namespace Foam

#endif
