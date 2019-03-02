/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

inline const vector Foam::DiscreteVelocityPoint::xi() const
{
    return xi_;
}

inline scalar Foam::DiscreteVelocityPoint::weight() const
{
    return weight_;
}

inline const volScalarField& Foam::DiscreteVelocityPoint::hVol() const
{
    return hVol_;
}

inline const volScalarField& Foam::DiscreteVelocityPoint::bVol() const
{
    return bVol_;
}

inline const volVectorField& Foam::DiscreteVelocityPoint::hGradVol() const
{
    return hGradVol_;
}

inline const volVectorField& Foam::DiscreteVelocityPoint::bGradVol() const
{
    return bGradVol_;
}

inline volScalarField& Foam::DiscreteVelocityPoint::hVol_par()
{
    return hVol_;
}

inline volScalarField& Foam::DiscreteVelocityPoint::bVol_par()
{
    return bVol_;
}

inline volVectorField& Foam::DiscreteVelocityPoint::hGradVol_par()
{
    return hGradVol_;
}

inline volVectorField& Foam::DiscreteVelocityPoint::bGradVol_par()
{
    return bGradVol_;
}

inline const surfaceScalarField& Foam::DiscreteVelocityPoint::hFluxSurf() const
{
    return hFluxSurf_;
}

inline const surfaceScalarField& Foam::DiscreteVelocityPoint::bFluxSurf() const
{
    return bFluxSurf_;
}
// ************************************************************************* //
