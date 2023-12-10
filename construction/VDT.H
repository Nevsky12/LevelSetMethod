#pragma once
#include "HashSet.H"
#include "fvMesh.H"
namespace Foam
{
namespace levelSet::construction
{

vectorField calcDist( vectorField   const &cellCentres
		    , labelListList const &nei
		    , HashSet<label, Hash<label>> const &start
		    , vectorField const &initDist
		    ) noexcept
{
    vectorField res = initDist;

    auto const min = [](vector const &v1, vector const &v2) noexcept
	             -> vector
    {
        return mag(v1) < mag(v2)
	     ?     v1  :     v2;
    };

    std::vector<bool> fixedDist(cellCentres.size(), false);
    HashSet<label, Hash<label>> curr = start;
    HashSet<label, Hash<label>> next;

    while(true)
    {
	forAllConstIters(curr, it)
	{
           fixedDist[*it] = true;
	}

        forAllConstIters(curr, it)
	{
	    label const i = *it;
	    
	    forAll(nei[i], j)
	    {
	        label const k = nei[i][j];

		if(!fixedDist[k])
		{
		   res[k] = min(cellCentres[k] - cellCentres[i] + res[i], res[k]);
		   next.insert(k);
		}
	    }
	}

	if(next.empty())
           return res;
	else
	{
	    curr = next;
            next = HashSet<label, Hash<label>>(0);
	}
    }
}

template<typename Domain>
auto VDT(fvMesh const &mesh, Domain const &pp) noexcept
{
    HashSet<label, Hash<label>> start;
    labelListList nei(mesh.cells().size());

    vectorField  dist(mesh.cells().size(), vector::rootMax);

    auto const begin = pp.begin();
    for(unsigned i = 0u; i < unsigned(pp.size()); ++i)
    {
        label const fi = *(begin + i);
	label const owner = mesh.faceOwner()[fi];

	dist[owner] = mesh.cellCentres()[owner] - mesh.faces()[fi].centre(mesh.points());
	start.insert(owner);
    }

    forAll(mesh.cells(), ci)
    {
        cell const c = mesh.cells()[ci];

	forAll(c, i)
	{
            label const fi = c[i];
	    if(fi > mesh.faceNeighbour().size())
	       continue;
            nei[ci].append( mesh.faceOwner    ()[fi] == ci 
			  ? mesh.faceNeighbour()[fi]
			  : mesh.faceOwner    ()[fi]      );
	}
    }

    return tmp<volVectorField>
    (
        new volVectorField
	(
	    IOobject
	    (
	        "UDF",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
            mesh,
	    dimensionSet(0, 1, 0, 0, 0, 0, 0),
	    calcDist(mesh.cellCentres(), nei, start, dist)
	)
    );
}

} // end namespace levelSet::construction
} // end namespace Foam
