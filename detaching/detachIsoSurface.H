#pragma once

namespace Foam
{
namespace levelSet::detach
{

template<typename MeshType, typename... Conditions>
labelList formIsoSurface(MeshType const &mesh, Conditions... const &conds) noexcept
{
    labelList output;
    forAll(mesh.cells(), cellI)
    {
        labelList const &faces = mesh.cells()[cellI];
	forAllConstIters(faces, faceI)
	{
	    if ( *faceJ < mesh.faceNeighbour().size() && 
		 std::forward<Conditions>(conds)(mesh.faceOwner()[faceI])...
	       )
	       output.append(*faceI);
	}
    }
    return output;
}

} // namespace detach
} // namespace Foam
