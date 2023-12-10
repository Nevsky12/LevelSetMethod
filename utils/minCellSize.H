#pragma once

namespace Foam
{
namespace levelSet::utils
{

inline scalar minCellSize(fvMesh const &mesh) noexcept
{
   edgeList const &edges = mesh.edges();
   pointField const &points = mesh.points();
   
   scalarField edgeLengths(edges.size(), 0.);

   forAll(edges, edgeI)
         edgeLengths[edgeI] = edges[edgeI].mag(points);

   return Foam::min(edgeLengths);
}

} // End namespace levelSet::utils
} // End namespace Foam
