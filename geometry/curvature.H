#pragma once                                                                       
namespace Foam
{ 
namespace levelSet::geometry
{

inline tmp<volScalarField> curvature(volScalarField const &psi) noexcept
{
   tmp<volScalarField> tCurvature
   (
      fvc::div
      ( 
               fvc::grad(psi) 
	/ (mag(fvc::grad(psi)) + sAdd) 
      )
   );

   volScalarField& curvatureRef = tCurvature.ref();
   curvatureRef.rename(psi.name() + " curvature");

   return tCurvature;
}

} // End namespace levelSet::geometry
} // End namespace Foam
