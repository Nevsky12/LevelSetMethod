#ifndef shockwaveLevetSetManager_H
#define shockwaveLevelSetManager_H

#include "levelSetManager.H"
namespace Foam
{
namespace levelSet
{

class shockwaveLevelSetManager final: private levelSetManager
{

	void transportPsi0(surfaceScalarField const &phi) noexcept;

public:

	shockwaveLevelSetManager
	( 
		volScalarField const &psi,
		scalar const epsilon  = 1.5,
		scalar const deltaTau = 0.1
	) noexcept;

	void updateState(surfaceScalarField const &phi, bool const usePhi0) noexcept;

	volScalarField const & psi () const noexcept {return psi_ ;}
	volScalarField const & psi0() const noexcept {return psi0_;}

};

} // End namespace levelSet
} // End namespace Foam

#endif
