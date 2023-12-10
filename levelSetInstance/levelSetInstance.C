#include "levelSetInstance.H"
#include "isoSurfaceCell.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "distributedTriSurfaceMesh.H"
#include "volPointInterpolation.H"
#include "fvc.H"
#include "fvm.H"
#include "VDT.H"
#include "signing.H"

namespace Foam::levelSet
{
    defineTypeNameAndDebug(levelSetInstance, 0);
} // namespace Foam

Foam::levelSet::levelSetInstance::levelSetInstance( volScalarField const &otherLS
	//	                                  , dictionary     const &dict
						  ) noexcept
:
    mesh_(otherLS.mesh()),
    levelSetFunc_(otherLS),
    //epsilon0_("epsilon", dimless, dict),
    epsilon_
    (
        IOobject
	(
	    IOobject::groupName("levelSet::thickness", otherLS.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimLength, 0),
	extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    H_
    (
        IOobject
	(
	    IOobject::groupName("H", otherLS.group()),
	    mesh_.time().timeName(),
	    mesh_,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	mesh_,
	dimensionedScalar(dimless, 0),
	zeroGradientFvPatchScalarField::typeName
    ),
    nSfd_
    (
        IOobject
	(
	    IOobject::groupName("surfaceNormalDirection", otherLS.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimArea, 0)
    ),
    K_
    (
        IOobject
	(
	    IOobject::groupName("curvature", otherLS.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimless/dimLength, 0)
    )
{}


Foam::levelSet::levelSetInstance::levelSetInstance( labelList  const &targetFaces
//		                                  , dictionary const &dict
						  , fvMesh     const &mesh
						  ) noexcept
:
    mesh_(mesh),
    levelSetFunc_
    (
        IOobject
	(
	    IOobject::groupName("levelSetFunction", "levelSet"),
	    mesh_.time().timeName(),
	    mesh_
	),
        levelSet::construction::signUDF
	(
	    mag(levelSet::construction::VDT(mesh_, targetFaces)),
	    targetFaces
	)
    ),
  //  epsilon0_("epsilon", dimless, dict),
    epsilon_
    (
        IOobject
	(
	    IOobject::groupName("levelSet::thickness", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimLength, 0),
	extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    H_
    (
        IOobject
	(
	    IOobject::groupName("H", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	mesh_,
	dimensionedScalar(dimless, 0),
	zeroGradientFvPatchScalarField::typeName
    ),
    nSfd_
    (
        IOobject
	(
	    IOobject::groupName("surfaceNormalDirection", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimArea, 0)
    ),
    K_
    (
        IOobject
	(
	    IOobject::groupName("curvature", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimless/dimLength, 0)
    )
{
    this->redistance(0.);
}


Foam::levelSet::levelSetInstance::levelSetInstance(fvMesh const &mesh) noexcept
:
    mesh_(mesh),
    levelSetFunc_
    (
        IOobject
	(
	    IOobject::groupName("levelSetFunction", "levelSet"),
	    mesh_.time().timeName(),
	    mesh_
	),
        mesh_,
	dimensionedScalar(dimLength, 0),
	zeroGradientFvPatchScalarField::typeName
    ),
  //  epsilon0_("epsilon", dimless, dict),
    epsilon_
    (
        IOobject
	(
	    IOobject::groupName("levelSet::thickness", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimLength, 0),
	extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    H_
    (
        IOobject
	(
	    IOobject::groupName("H", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	mesh_,
	dimensionedScalar(dimless, 0),
	zeroGradientFvPatchScalarField::typeName
    ),
    nSfd_
    (
        IOobject
	(
	    IOobject::groupName("surfaceNormalDirection", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimArea, 0)
    ),
    K_
    (
        IOobject
	(
	    IOobject::groupName("curvature", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
	mesh_,
	dimensionedScalar(dimless/dimLength, 0)
    )
{}

Foam::levelSet::levelSetInstance::~levelSetInstance()
{}

void Foam::levelSet::levelSetInstance::createLSF(labelList const &isoSurf) noexcept
{
    levelSetFunc_ = levelSet::construction::signUDF
    (
        mag(levelSet::construction::VDT(mesh_, isoSurf)),
        isoSurf
    );
    this->redistance(0.);
}

void Foam::levelSet::levelSetInstance::redistance(scalar const isoValue) noexcept
{
    volScalarField const &isoField = levelSetFunc_;
    pointScalarField pointIsoField
    (
        volPointInterpolation::New(mesh_).interpolate(isoField)
    );

    isoSurfaceCell contour
    (
        mesh_,
        isoField,
        pointIsoField,
        isoValue
    );

    contour.triangulate();

    triFaceList triFaces(contour.size());
    forAll(contour, facei)
    {
        triFaces[facei][0] = contour[facei][0];
        triFaces[facei][1] = contour[facei][1];
        triFaces[facei][2] = contour[facei][2];
    }

    autoPtr<triSurfaceMesh> triMeshPtr;
    triSurface tri(triFaces, contour.points());
    
    bool const useDistributed_ = true;

    dictionary triMeshDict_;
    
    if (Pstream::parRun() && useDistributed_)
    {
        triMeshDict_.set("bounds", List<boundBox>(1, mesh_.bounds()));
        triMeshPtr.reset
        (
            new distributedTriSurfaceMesh
            (
                IOobject
                (
                    "contour_" + isoField.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                tri,
                triMeshDict_
            )
        );
    }
    else
    {
        triMeshPtr.reset
        (
            new triSurfaceMesh
            (
                IOobject
                (
                    "contour_" + isoField.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                tri
            )
        );
    }
    triSurfaceMesh& triMesh = triMeshPtr();

    if (debug)
    {
        triMesh.triSurface::write("contour_" + isoField.name() + ".stl");
    }

    volScalarField ls
    (
        IOobject
	(
	    IOobject::groupName("levelSet", levelSetFunc_.group()),
	    mesh_.time().timeName(),
	    mesh_
	),
        mesh_,
        dimensionedScalar(dimLength, -GREAT),
        zeroGradientFvPatchScalarField::typeName
    );

    pointField samples(mesh_.C());
    scalarField nearestDistSqr(samples.size(), magSqr(mesh_.bounds().span()));
    List<pointIndexHit> hitPoints(samples.size());

    triMesh.findNearest(samples, nearestDistSqr, hitPoints);

    forAll(mesh_.C(), celli)
    {
        ls[celli] =
            mag(mesh_.C()[celli] - hitPoints[celli].rawPoint())
           *(levelSetFunc_[celli] > 0. ? 1.0 : -1.0);
    }
    ls.correctBoundaryConditions();

    levelSetFunc_ = ls;
}


Foam::labelList Foam::levelSet::levelSetInstance::interface() const noexcept
{
    labelList const &ow  = mesh_.faceOwner();
    labelList const &nei = mesh_.faceNeighbour();
    labelList res; 

    forAll(mesh_.faces(), fI)
    {
        if(mesh_.faceNeighbour().size() >= fI && levelSetFunc_[ow[fI]] * levelSetFunc_[nei[fI]] < 0.)
	   res.append(fI);
    }

    return res;
}

Foam::tmp<Foam::surfaceScalarField> Foam::levelSet::levelSetInstance::interfaceField() const noexcept
{
    tmp<surfaceScalarField> const tsf
    (
        surfaceScalarField::New
	(
	    IOobject::groupName("levelSet", levelSetFunc_.group()),
	    mesh_,
	    dimensionedScalar(dimLength, 0)
	)
    );
    surfaceScalarField &self = tsf.ref();

    surfaceScalarField const lsff = fvc::interpolate(levelSetFunc_);
    labelList const isoSurf = this->interface();
    for(label const fI: isoSurf)
	self[fI] = lsff[fI];

    return tsf;
}

void Foam::levelSet::levelSetInstance::transportLSF(surfaceScalarField const &shockFlow) noexcept
{
    //solve(fvm::ddt(levelSetFunc_) + fvm::div(shockFlow, levelSetFunc_));
    //this->redistance(0.);
}

Foam::tmp<Foam::volScalarField> Foam::levelSet::levelSetInstance::interfaceIndicator() const noexcept
{
    return volScalarField::New
    (
        IOobject::groupName("interfaceIndicator", levelSetFunc_.group()),
	pos0(mag(levelSetFunc_) + 2. * epsilon_)
    );
}

Foam::tmp<Foam::volVectorField> Foam::levelSet::levelSetInstance::   interfaceNormal() const noexcept
{
    volVectorField const gradLS = fvc::grad(levelSetFunc_);
    return volVectorField
    (
         "interfaceNormal",
	 (
	       gradLS
	     /
	       max
	       (
	           mag(gradLS),
		   dimensionedScalar(gradLS.dimensions(), 1e-6)
	       )
	 )
    );
}
