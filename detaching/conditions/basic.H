#pragma once
#include "Vector.H"
namespace Foam
{
namespace levelSet::condition
{


inline auto isPositive(volScalarField      const &f) noexcept
{
    return [val = f  ](label const cellI) noexcept -> bool
    {
        return val[cellI] > scalar(0.);
    };
}

inline auto isPositive(tmp<volScalarField> const &f) noexcept
{
    return [val = f()](label const cellI) noexcept -> bool
    {
        return val[cellI] > scalar(0.);
    };
}



inline auto isNegative(volScalarField      const &f) noexcept
{
    return [val = f  ](label const cellI) noexcept -> bool
    {
        return val[cellI] < scalar(0.);
    };
}

inline auto isNegative(tmp<volScalarField> const &f) noexcept
{
    return [val = f()](label const cellI) noexcept -> bool
    {
        return val[cellI] < scalar(0.);
    };
}



template<Vector::components Cmp>
inline auto isPositiveDer(volScalarField      const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) > scalar(0.);
    };
}

template<Vector::components Cmp>
inline auto isPositiveDer(tmp<volScalarField> const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) > scalar(0.);
    };
}



template<Vector::components Cmp>
inline auto isNegativeDer(volScalarField       const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) < scalar(0.);
    };
}

template<Vector::components Cmp>
inline auto isNegativeDer(tmp<volScalarField> const &f) noexcept
{
    return [der = fvc::grad(f)()](label const cellI) noexcept -> bool
    {
        return der[cellI].component(Cmp) < scalar(0.);
    };
}

} // namespace levelSet::condition
} // namespace Foam
