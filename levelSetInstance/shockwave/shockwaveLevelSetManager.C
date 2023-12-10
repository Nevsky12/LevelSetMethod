#include "shockwaveLevelSetManager.H"
#include "minCellSize.H"
#include "fvmDdt.H"
#include "fvmDiv.H"

Foam::levelSet::shockwaveLevelSetManager::shockwaveLevelSetManager
(
	volScalarField const &psi,
	scalar const epsilon,
	scalar const deltaTau	
) noexcept
:
	Foam::levelSet::levelSetManager
	(
		psi,
		epsilon,
		deltaTau
	)
{}


void Foam::levelSet::shockwaveLevelSetManager::transportPsi0(surfaceScalarField const &phi) noexcept
{	
	solve
	(
		fvm::ddt(psi0_)
	    +   fvm::div(phi, psi0_)
	);
}


void Foam::levelSet::shockwaveLevelSetManager::updateState(surfaceScalarField const &phi, bool const usePhi0) noexcept
{
	reinitialize(usePhi0);
	transportPsi0(phi);
	funcsManager_->updateAndCorrect(psi_, epsilon_);
	reinitialize(usePhi0);
}

