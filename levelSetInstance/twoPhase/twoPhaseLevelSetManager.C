#include "twoPhaseLevelSetManager.H"
#include "minCellSize.H"

Foam::levelSet::twoPhaseLevelSetManager::twoPhaseLevelSetManager
(
	volScalarField const &psi,
	scalar const epsilon,
        scalar const gamma,
	scalar const deltaTau	
) noexcept
:
	Foam::levelSet::levelSetManager
	(
		psi,
		epsilon,
		deltaTau
	),
	gamma_(utils::minCellSize(psi.mesh()) * gamma)
{}


void Foam::levelSet::twoPhaseLevelSetManager::psi0FromVolumeFraction(volScalarField const &alpha) noexcept
{
	psi0_ == (2. * alpha - 1.) * gamma_; 
}


Foam::tmp<Foam::volScalarField> Foam::levelSet::twoPhaseLevelSetManager::reconstructTwoPhaseField
( 
	dimensionedScalar const prop1, 
	dimensionedScalar const prop2
) const noexcept                                                          
{
   tmp<volScalarField> tNewField
   (
      new volScalarField
      (
         "tmp",
	 psi_ * 0. * prop1
      )
   );
   volScalarField &self = const_cast<volScalarField &>(tNewField());

   self = prop1 *       funcsManager_->gHeaviside() 
	+ prop2 * (1. - funcsManager_->gHeaviside());

   return tNewField;
}


void Foam::levelSet::twoPhaseLevelSetManager::updateState(volScalarField const &alpha, bool const usePhi0) noexcept
{
	psi0FromVolumeFraction(alpha);
	funcsManager_->updateAndCorrect(psi_, epsilon_);
	reinitialize(usePhi0);
}

