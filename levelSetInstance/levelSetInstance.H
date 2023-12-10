#pragma once
#include "volFields.H"
#include "surfaceFields.H"
namespace Foam
{
namespace levelSet
{


class levelSetInstance
{
    
private:

    fvMesh const &mesh_;

    volScalarField levelSetFunc_;

    dimensionedScalar epsilon0_;

    volScalarField epsilon_;

    volScalarField H_; 

    surfaceScalarField nSfd_;

    volScalarField K_;   


public:

    TypeName("levelSetInstance");

        
    levelSetInstance( volScalarField const &otherLS
//		    , dictionary     const &dict
		    ) noexcept;

    levelSetInstance( labelList  const &targetFaces
//		    , dictionary const &dict
		    , fvMesh     const &mesh
                    ) noexcept;

    levelSetInstance(fvMesh const &mesh) noexcept;

    virtual ~levelSetInstance();


    volScalarField     const &     levelSetFunc() const noexcept
    {
        return levelSetFunc_;
    }

    volScalarField           &     levelSetFunc()       noexcept
    {
        return levelSetFunc_;
    }

    volScalarField     const &          epsilon() const noexcept
    {
        return epsilon_;
    }

    volScalarField     const &        heaviside() const noexcept
    {
        return H_;
    }

    surfaceScalarField const & surfaceNormalDir() const noexcept
    {
        return nSfd_;
    }

    volScalarField     const &        curvature() const noexcept
    {
        return K_;
    }

    labelList                         interface     () const noexcept;

    tmp<surfaceScalarField>           interfaceField() const noexcept;

    void createLSF(labelList const &isoSurf     ) noexcept;

    void redistance(scalar   const isoValue = 0.) noexcept;

    void transportLSF(surfaceScalarField const &shockFlow) noexcept;

    tmp<volScalarField> interfaceIndicator() const noexcept;

    tmp<volVectorField>    interfaceNormal() const noexcept;

};

} // namespace levelSet
} // namespace Foam
