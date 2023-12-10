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



template<typename Type>
auto moreThanThresholdedMax(volScalarField      const &f, scalar const threshold) noexcept
{
    assert(threshold < scalar(1.));
    return [val = f,   maxfValue = max(f)](label const cellI) noexcept -> bool
    {
        return val[cellI] > threshold * maxfValue;
    };
}

template<typename Type>
auto moreThanThresholdedMax(tmp<volScalarField> const &f, scalar const threshold) noexcept
{
    assert(threshold < scalar(1.));
    return [val = f(), maxfValue = max(f)](label const cellI) noexcept -> bool
    {
        return val[cellI] > threshold * maxfValue;
    };
}



auto lessThanMinGained(volScalarField      const &f, scalar const gain) noexcept
{
    assert(gain > scalar(1.));
    return [val = f,   minfValue = max(f)](label const cellI) noexcept -> bool
    {
        return val[cellI] < gain * minfValue;
    };
}

auto lessThanMinGained(tmp<volScalarField> const &f, scalar const gain) noexcept
{
    assert(gain > scalar(1.));
    return [val = f(), minfValue = min(f)](label const cellI) noexcept -> bool
    {
        return val[cellI] < gain * minfValue;
    };
}



auto moreThanVolumeWeightedAverage(volScalarField      const &f) noexcept
{
    return [val = f, avefValue = f.    weightedAverage(f.mesh().V())](label const cellI) noexcept -> bool
    {
        return val[cellI] > avefValue.value();
    };
}

auto moreThanVolumeWeightedAverage(tmp<volScalarField> const &f) noexcept
{
    return [val = f(), avefValue = f().weightedAverage(f.mesh().V())](label const cellI) noexcept -> bool
    {
        return val[cellI] > avefValue.value();
    };
}



auto lessThanVolumeWeightedAverage(volScalarField      const &f) noexcept
{
    return [val = f,   avefValue = f.  weightedAverage(f.mesh().V())](label const cellI) noexcept -> bool
    {
        return val[cellI] < avefValue.value();
    };
}

auto lessThanVolumeWeightedAverage(tmp<volScalarField> const &f) noexcept
{
    return [val = f(), avefValue = f().weightedAverage(f.mesh().V())](label const cellI) noexcept -> bool
    {
        return val[cellI] < avefValue.value();
    };
}

} // namespace levelSet::condition
} // namespace Foam
