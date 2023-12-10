#pragma once
#include "volFields.H"
namespace Foam
{
namespace levelSet
{

class Manifold
{
private:

    fvMesh const &mesh_;

    volScalarField const &p_;

    volScalarField const &rho_;

    volVectorField const &U_;

    volScalarField const &c_;

    volVectorField gradP_;
    
    volScalarField magGradP_;

    scalar maxMagGradP_;

    scalar aveMagGradP_;

    scalar maxP_;

    volVectorField gradRho_;

    volScalarField magGradRho_;

    scalar maxMagGradRho_;

    scalar aveMagGradRho_;

    scalar maxRho_;

    volScalarField divU_;
 
    dictionary pDict_;

    dictionary rhoDict_;
    
    dictionary UDict_;

    scalar pGradMaxLimit_;
    scalar pGradSpecific_;
    scalar pGradSpecificMax_;

    scalar rhoGradMaxLimit_;
    scalar rhoGradSpecific_;
    scalar rhoGradSpecificMax_;
    
    scalar divUSpecific_;


    bool checkPressureConds(label const cellI) const noexcept;

    bool  checkDensityConds(label const cellI) const noexcept;

    bool checkVelocityConds(label const cellI) const noexcept;

public:


    static List<word> const scalarFieldConditions_;
    
    static List<word> const vectorFieldConditions_;


    Manifold( volScalarField const &p
	    , volScalarField const &rho
	    , volVectorField const &U
	    , volScalarField const &c
	    , IOdictionary   const &dict
	    ) noexcept;

    bool isSatisfied(label const cellI) const noexcept;
};

} // namespace levelSet
} // namespace Foam
