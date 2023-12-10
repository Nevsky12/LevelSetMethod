#pragma once                                                                       
namespace Foam
{ 
namespace levelSet::geometry
{

inline tmp<volScalarField> normal(volScalarField const &psi) noexcept
{
   tmp<volScalarField> tNormal
   (
             fvc::grad(psi) 
      / 
        (mag(fvc::grad(psi)) + sAdd)  
   );

   volScalarField& normalRef = tNormal.ref();
   normalRef.rename(psi.name() + " normal");

   return tNormal;
}

} // End namespace levelSet::geometry
} // End namespace Foam
