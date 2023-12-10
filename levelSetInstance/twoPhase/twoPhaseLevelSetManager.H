#ifndef twoPhaseLevetSetManager_H
#define twoPhaseLevelSetManager_H

#include "levelSetManager.H"
namespace Foam
{
namespace levelSet
{

class twoPhaseLevelSetManager final: public levelSetManager
{

	scalar gamma_;

public:

	void psi0FromVolumeFraction(volScalarField const &alpha) noexcept;

	twoPhaseLevelSetManager
	( 
		volScalarField const &psi,
		scalar const epsilon  = 1.5,
		scalar const gamma    = 0.75,
		scalar const deltaTau = 0.1
	) noexcept;

        tmp<volScalarField> reconstructTwoPhaseField
	(
		dimensionedScalar const prop1,
		dimensionedScalar const prop2
	) const noexcept;

	void updateState(volScalarField const &alpha, bool const usePhi0) noexcept;

	volScalarField const & psi () const noexcept {return psi_ ;}
	volScalarField const & psi0() const noexcept {return psi0_;}

};

} // End namespace levelSet
} // End namespace Foam

#endif
