#pragma once                                                                       
namespace Foam
{ 
namespace levelSet::geometry
{

inline tmp<volScalarField> gradient(volScalarField const &psi) noexcept
{
   tmp<volScalarField> tGradPsi(fvc::grad(psi));

   volScalarField& psiRef = tGradPsi.ref();
   psiRef.rename(psi.name() + " gradient");

   return tGradPsi;
}

} // End namespace levelSet::geometry
} // End namespace Foam
