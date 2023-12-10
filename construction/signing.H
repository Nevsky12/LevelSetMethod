#pragma once
#include "polyMesh.H"
#include "meshTools.H"
#include <cassert>
#include <queue>
#include <set>
namespace Foam
{
namespace levelSet::construction
{

tmp<volScalarField> signUDF( volScalarField const &UDF
		           , labelList      const &surfaceFaces
		           ) noexcept
{
    fvMesh const &mesh = UDF.mesh();
    tmp<volScalarField> SDF = tmp<volScalarField>::New(UDF);
    volScalarField &self = SDF.ref();

    std::stack<label> forColor;
    point centerMass = point::zero;
    DynamicList<label> band;
    for(label const sf: surfaceFaces)
    {
	 band.append(mesh.faceOwner()[sf]);
         centerMass += mesh.cellCentres()[mesh.faceOwner()[sf]];
    }
 
    centerMass /= surfaceFaces.size();
    centerMass += point(ROOTSMALL, ROOTSMALL, ROOTSMALL);
    assert(mesh.findCell(centerMass) != -1);
    forColor.push(mesh.findCell(centerMass));

    DynamicList<label> colored;
    while(!forColor.empty())
    {
        label const cI = forColor.top();
        forColor.pop();

        if(!band.found(cI))
        {
            if(!colored.found(cI))
            {
                self[cI] *= -1.0;
                colored.append(cI);
            }
            for(label const cJ: mesh.cellCells()[cI])
            {
                if(!colored.found(cJ))
                   forColor.push(cJ);
            }
        }
     }

    return SDF;
}	

tmp<volScalarField> signUDF( tmp<volScalarField> const &tUDF
		           , labelList           const &surfaces
		           ) noexcept
{
    return signUDF(tUDF(), surfaces);
}

} // end namespace levelSet::construction
} // end namespace Foam
