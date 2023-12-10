#ifndef levelSetManager_H
#define levelSetManager_H

#include "funcsSmartStorage.H"
namespace Foam
{
namespace levelSet
{

class levelSetManager
{

protected:

	volScalarField psi_;

	volScalarField psi0_;

        scalar epsilon_;

	scalar deltaTau_;

	autoPtr<calculus::funcsSmartStorage> funcsManager_;

public:

	levelSetManager
	( 
		volScalarField const &psi,
		scalar const epsilon  = 1.5,
		scalar const deltaTau = 0.1
	) noexcept;

	void reinitialize(bool const usePhi0) noexcept;
	tmp<volScalarField> reinitialize(volScalarField const &psi0, bool const usePhi0) noexcept;

	void transport(surfaceScalarField const &W) noexcept;

        void updatePsi(DimensionedField<scalar, volMesh> const &newLS) noexcept;
	volScalarField const &psi() const noexcept {return psi_;}
	volScalarField const &psi0() const noexcept {return psi0_;}

};

} // End namespace levelSet
} // End namespace Foam 

#endif
