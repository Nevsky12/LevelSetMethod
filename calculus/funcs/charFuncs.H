#ifndef charFuncs_H
#define charFuncs_H

#include "volFields.H"

namespace Foam
{ 
namespace levelSet::calculus
{

inline tmp<volScalarField> charFuncPos(volScalarField const &psi) noexcept 
{
    tmp<volScalarField> tCharFuncPos
    (
    	new volScalarField
    	(
        	"CharPosFunc",
		psi * 0.
    	)
    );
    volScalarField &self = const_cast<volScalarField &>(tCharFuncPos());

    forAll(psi.mesh().cells(), cellI)
    {
        if (psi[cellI] <= 0)
            self[cellI] = 0;
	else
	    self[cellI] = 1;
    }

    return tCharFuncPos;
}

inline tmp<volScalarField> charFuncNeg(volScalarField const &psi) noexcept 
{
    tmp<volScalarField> tCharFuncNeg
    (
        new volScalarField
    	(
        	"CharNegFunc",
		psi * 0.
    	)
    );
    volScalarField &self = const_cast<volScalarField &>(tCharFuncNeg());

    forAll(psi.mesh().cells(), cellI)
    {
        if (psi[cellI] <= 0)
            self[cellI] = 1;
	else
	    self[cellI] = 0;
    }

    return tCharFuncNeg;
}

} // End namespace levelSet::calculus
} // End namespace Foam

#endif
