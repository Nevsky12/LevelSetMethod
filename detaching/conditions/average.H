#pragma once
#include "Vector.H"
namespace Foam
{
namespace levelSet::condition
{


auto isPositive(volScalarField      const &f) noexcept
{
    return [val = f  ](label const cellI) noexcept -> bool
    {
        return val[cellI] > scalar(0.);
    };
}

auto isPositive(tmp<volScalarField> const &f) noexcept
{
    return [val = f()](label const cellI) noexcept -> bool
    {
        return val[cellI] > scalar(0.);
    };
}



auto isNegative(volScalarField      const &f) noexcept
{
    return [val = f  ](label const cellI) noexcept -> bool
    {
        return val[cellI] < scalar(0.);
    };
}

auto isNegative(tmp<volScalarField> const &f) noexcept
{
    return [val = f()](label const cellI) noexcept -> bool
    {
        return val[cellI] < scalar(0.);
    };
}



template<Vector::components Cmp>
auto isPositiveDer(volScalarField      const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) > scalar(0.);
    };
}

template<Vector::components Cmp>
auto isPositiveDer(tmp<volScalarField> const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) > scalar(0.);
    };
}



template<Vector::components Cmp>
auto isNegativeDer(volScalarField       const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) < scalar(0.);
    };
}

template<Vector::components Cmp>
auto isNegativeDer(tmp<volScalarField> const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) < scalar(0.);
    };
}



} // namespace levelSet::condition
} // namespace Foam
