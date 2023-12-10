#include "manifold.H"
#include "fvc.H"

Foam::List<Foam::word> const Foam::levelSet::Manifold::scalarFieldConditions_ =
{
    "positiveDerX",
    "positiveDerY",
    "positiveDerZ",

    "magGradMoreThanBoundedMax",
    "magGradMoreThanAverage",

    "magGradMoreThanSpecific",
    "magGradMoreThanSpecificMax",
};

Foam::List<Foam::word> const Foam::levelSet::Manifold::vectorFieldConditions_ =
{
    "positiveCompression",
    "compressionMoreThanSpecific",
    "compressionMoreThanSpecificMax",
};

Foam::levelSet::Manifold::Manifold( volScalarField const &p
		                  , volScalarField const &rho
				  , volVectorField const &U
				  , volScalarField const &c
				  , IOdictionary   const &dict
		                  ) noexcept
:
    mesh_(p.mesh()),

    p_  (p  ),
    rho_(rho),
    U_  (U  ),
    c_  (c  ),

    gradP_      (fvc::grad(p  )),
    magGradP_   (mag(   gradP_)),
    maxMagGradP_(max(magGradP_).value()),
    aveMagGradP_(magGradP_.weightedAverage(mesh_.V()).value()),

    gradRho_      (fvc::grad(rho  )),
    magGradRho_   (mag(   gradRho_)),
    maxMagGradRho_(max(magGradRho_).value()),
    aveMagGradRho_(magGradRho_.weightedAverage(mesh_.V()).value()),

    divU_   (fvc::div(U))
{ 
   pDict_   = dict.subDict("pFilters"  ); 
   rhoDict_ = dict.subDict("rhoFilters");
   UDict_   = dict.subDict("UFilters"  );

   auto const extractParam = [](dictionary const &fDict, word const &pName) noexcept
	                     -> scalar
   { 
       return fDict.found(pName) ? fDict.get<scalar>(pName) : -1.;
   };

   pGradMaxLimit_    = extractParam(pDict_, "gradMaxLimit"   ); 
   pGradSpecific_    = extractParam(pDict_, "gradSpecific"   );
   pGradSpecificMax_ = extractParam(pDict_, "gradSpecificMax");
   
   rhoGradMaxLimit_    = extractParam(rhoDict_, "gradMaxLimit"   ); 
   rhoGradSpecific_    = extractParam(rhoDict_, "gradSpecific"   );
   rhoGradSpecificMax_ = extractParam(rhoDict_, "gradSpecificMax");

   divUSpecific_ = extractParam(UDict_, "divSpecific");
}

bool Foam::levelSet::Manifold::checkPressureConds(label const cellI) const noexcept
{
    if(pDict_.get<bool>(scalarFieldConditions_[0]))
    if(gradP_[cellI].component(vector::X) <= 0.)
       return false;

    if(pDict_.get<bool>(scalarFieldConditions_[1]))
    if(gradP_[cellI].component(vector::Y) <= 0)
       return false;
       
    if(pDict_.get<bool>(scalarFieldConditions_[2]))
    if(gradP_[cellI].component(vector::Z) <= 0)
       return false;

    if(pDict_.get<bool>(scalarFieldConditions_[3]))
    if(magGradP_[cellI] < pGradMaxLimit_ * maxMagGradP_)
       return false;

    if(pDict_.get<bool>(scalarFieldConditions_[4]))
    if(magGradP_[cellI] < aveMagGradP_)
       return false;
 
    if(pDict_.get<bool>(scalarFieldConditions_[6]) || pDict_.get<bool>(scalarFieldConditions_[5]))
    {
        scalar d = SMALL;
        for(point const &pI: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))        
        for(point const &pJ: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))
            d = max(d, mag(pI - pJ));

         scalar const h = d / std::sqrt(2.);

	 if(pDict_.get<bool>(scalarFieldConditions_[6]))
         if(magGradP_[cellI] < pGradSpecific_ * p_[cellI] / h)
            return false;

	 if(pDict_.get<bool>(scalarFieldConditions_[5]))
         if(magGradP_[cellI] < pGradSpecificMax_ * maxP_ / h)
            return false;
    }

    return true;
}


bool Foam::levelSet::Manifold::checkDensityConds(label const cellI) const noexcept
{
    if(rhoDict_.get<bool>(scalarFieldConditions_[0]))
    if(gradRho_[cellI].component(vector::X) <= 0)
       return false;

    if(rhoDict_.get<bool>(scalarFieldConditions_[1]))
    if(gradRho_[cellI].component(vector::Y) <= 0)
       return false;
       
    if(rhoDict_.get<bool>(scalarFieldConditions_[2]))
    if(gradRho_[cellI].component(vector::Z) <= 0)
       return false;

    if(rhoDict_.get<bool>(scalarFieldConditions_[3]))
    if(magGradRho_[cellI] < rhoGradMaxLimit_ * maxMagGradRho_)
       return false;
 
    if(rhoDict_.get<bool>(scalarFieldConditions_[4]))
    if(magGradRho_[cellI] < aveMagGradRho_) 
       return false;

    if(rhoDict_.get<bool>(scalarFieldConditions_[6]) || rhoDict_.get<bool>(scalarFieldConditions_[5]))
    {
        scalar d = SMALL;
        for(point const &pI: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))        
        for(point const &pJ: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))
            d = max(d, mag(pI - pJ));

         scalar const h = d / std::sqrt(2.);

	 if(rhoDict_.get<bool>(scalarFieldConditions_[6]))
         if(magGradRho_[cellI] < rhoGradSpecific_ * rho_[cellI] / h)
            return false;

	 if(rhoDict_.get<bool>(scalarFieldConditions_[5]))
         if(magGradRho_[cellI] < rhoGradSpecificMax_ * maxRho_ / h)
            return false;
    }

    return true;
}

bool Foam::levelSet::Manifold::checkVelocityConds(label const cellI) const noexcept
{
    if(UDict_.get<bool>(vectorFieldConditions_[0]))
    if(divU_[cellI] >= 0)
       return false;

    if(UDict_.get<bool>(vectorFieldConditions_[1]))
    {
        scalar d = SMALL;
        for(point const &pI: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))        
        for(point const &pJ: mesh_.cells()[cellI].points(mesh_.faces(), mesh_.points()))
            d = max(d, mag(pI - pJ));

        scalar const h = d / std::sqrt(2.);
        if(divU_[cellI] >= -divUSpecific_ * c_[cellI] / h)
           return false;
    }

    return true;
}

bool Foam::levelSet::Manifold::isSatisfied(label const cellI) const noexcept
{
    return checkPressureConds(cellI)
        &&  checkDensityConds(cellI)
        && checkVelocityConds(cellI);
}
