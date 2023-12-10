#include "levelSetManager.H"
#include "minCellSize.H"
#include "fvcDiv.H"
#include "fvmDdt.H"
#include "fvcGrad.H"


void Foam::levelSet::levelSetManager::reinitialize(bool const usePhi0) noexcept
{
   psi0_.correctBoundaryConditions();
   psi_ == psi0_;

   volScalarField sgnPsi
   (
       "sgnPsi",
       psi_ * 0.
   );


   if (usePhi0)
       sgnPsi = psi0_ / mag(psi0_ + 1e-12);
   else
       sgnPsi = psi_  / mag(psi_  + 1e-12);

   dimensionedScalar const addLength("pickLength", dimLength, 1.);
   
   for (unsigned int i = 0u; i < static_cast<unsigned int>(epsilon_ / deltaTau_); ++i)
   {
	psi_ = psi_ + sgnPsi * (1. - mag(fvc::grad(psi_ * addLength))) * deltaTau_;
   	psi_.correctBoundaryConditions();
   }
}


Foam::tmp<Foam::volScalarField> Foam::levelSet::levelSetManager::reinitialize(volScalarField const &psi0, bool const usePhi0) noexcept
{
   tmp<volScalarField> psiT
   (
       new volScalarField(psi0)
   );
   volScalarField &psi = psiT.ref();

   volScalarField sgnPsi
   (
       "sgnPsi",
       psi * 0.
   );


   if (usePhi0)
       sgnPsi = psi0 / mag(psi0 + 1e-12);
   else
       sgnPsi = psi  / mag(psi  + 1e-12);

   dimensionedScalar const addLength("pickLength", dimLength, 1.);
   
   for (unsigned int i = 0u; i < static_cast<unsigned int>(epsilon_ / deltaTau_); ++i)
   {
	psi = psi + sgnPsi * (1. - mag(fvc::grad(psi * addLength))) * deltaTau_;
   	psi.correctBoundaryConditions();
   }

   return psiT;
}

void Foam::levelSet::levelSetManager::transport(surfaceScalarField const &W) noexcept
{
	solve(fvm::ddt(psi_) + fvc::div(W, psi_));
}

void Foam::levelSet::levelSetManager::updatePsi(DimensionedField<scalar, volMesh> const &newLS) noexcept
{
    psi_.internalFieldRef() = newLS;
}

Foam::levelSet::levelSetManager::levelSetManager
( 
   volScalarField const &psi, 
   scalar const epsilon,  
   scalar const deltaTau
) noexcept
:
   psi_(psi),
   psi0_
   (
      IOobject
      (
         "psi0",
	 psi.time().timeName(),
	 psi.mesh()
      ),
      psi.mesh(),
      dimensionedScalar("psi0", dimless, 0.),
      psi.boundaryField().types()
   ),
   epsilon_ (Foam::levelSet::utils::minCellSize(psi.mesh()) *  epsilon),
   deltaTau_(Foam::levelSet::utils::minCellSize(psi.mesh()) * deltaTau),
   funcsManager_(new Foam::levelSet::calculus::funcsSmartStorage(psi_, epsilon_))
{
	psi0_ == psi;
}
