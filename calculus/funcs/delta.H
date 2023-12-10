#ifndef delta_H
#define delta_H

#include "volFields.H"

namespace Foam
{ 
namespace levelSet::calculus
{

inline tmp<volScalarField> delta(volScalarField const &psi, scalar const epsilon) noexcept 
{
    tmp<volScalarField> tDelta
    (
    	new volScalarField
    	(	
        	"DeltaFunc",
		psi * 0.
    	)
    );
    volScalarField &self = const_cast<volScalarField &>(tDelta());

    forAll(psi.mesh().cells(), cellI)
    {
        if (mag(psi[cellI]) > epsilon)
               self[cellI] = 0;
	else
               self[cellI] = 1 / 2 / epsilon * (1 + Foam::sin(M_PI * psi[cellI] / epsilon));
    }

    return tDelta;
}

} // End namespace levelSet::calculus
} // End namespace Foam

#endif
