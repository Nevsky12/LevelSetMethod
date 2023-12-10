#ifndef heaviside_H
#define heaviside_H

#include "volFields.H"

namespace Foam
{ 
namespace levelSet::calculus
{

inline tmp<volScalarField> heaviside(volScalarField const &psi, scalar const epsilon) noexcept 
{
    tmp<volScalarField> tHeaviside
    (
    	new volScalarField
    	(
        	"Heaviside",
		psi * 0.
    	)
    );
    volScalarField &self = const_cast<volScalarField &>(tHeaviside());

    forAll(psi.mesh().cells(), cellI)
    {
             if (psi[cellI] < -epsilon)
                self[cellI] = 0;
        else if (psi[cellI] >  epsilon)
                self[cellI] = 1;
        else
	        self[cellI] = 1 / 2 
		              * ( 1 
				  +                           psi[cellI] / epsilon 
				  + M_1_PI * Foam::cos(M_PI * psi[cellI] / epsilon)
				);
    }

    self == min(max(self, 0.), 1.);

    return tHeaviside;
}

} // End namespace levelSet::calculus
} // End namespace Foam
  
#endif
