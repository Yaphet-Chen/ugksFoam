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

#include "fvDVM.H"
#include "DiscreteVelocityPoint.H"
#include <map>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DiscreteVelocityPoint::DiscreteVelocityPoint(
    fvDVM &dvm,
    const fvMesh &mesh,
    const Time &time,
    const scalar weight,
    const vector xi,
    const label DVid,
    const label symXtargetDVid,
    const label symYtargetDVid,
    const label symZtargetDVid)
    : dvm_(dvm),
      mesh_(mesh),
      time_(time),
      weight_(weight),
      xi_(xi),
      myDVid_(DVid),
      symXtargetDVid_(symXtargetDVid),
      symYtargetDVid_(symYtargetDVid),
      symZtargetDVid_(symZtargetDVid),
      hVol_(
          IOobject(
              "hVol_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0.0)),
      bVol_(
          IOobject(
              "bVol_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0.0)),
      hGradVol_(
          IOobject(
              "hGradVol_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", hVol_.dimensions() / dimLength, vector(0, 0, 0)),
          "zeroGradient"),
      bGradVol_(
          IOobject(
              "bGradVol_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", bVol_.dimensions() / dimLength, vector(0, 0, 0)),
          "zeroGradient"),
      hFluxSurf_(
          IOobject(
              "hFluxSurf_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0.0)),
      bFluxSurf_(
          IOobject(
              "bFluxSurf_" + name(DVid),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0.0))
{
    initDFtoEq();
    setBCtype();
    initBoundaryField();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DiscreteVelocityPoint::~DiscreteVelocityPoint()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::DiscreteVelocityPoint::initDFtoEq()
{
    DiscreteMaxwell(
        hVol_,
        bVol_,
        dvm_.rhoVol(),
        dvm_.Uvol(),
        dvm_.lambdaVol());
}

void Foam::DiscreteVelocityPoint::setBCtype()
{
    const volScalarField::Boundary &rhoBCs = dvm_.rhoVol().boundaryFieldRef();

    // Bondary condition type map from rho BC to h/bVol_ BC
    std::map<word, word> bcMap;
    bcMap["fixedValue"] = "fixedValue";
    bcMap["zeroGradient"] = "zeroGradient";
    bcMap["calculated"] = "calculated"; // Maxwell wall

    forAll(rhoBCs, patchi)
    {
        word rhoBCtype = rhoBCs[patchi].type();
        if (bcMap.find(rhoBCtype) != bcMap.end()) // found
        {
            hVol_.boundaryFieldRef().set(
                patchi,
                fvPatchField<scalar>::New(
                    bcMap[rhoBCtype], mesh_.boundary()[patchi], hVol_));
            bVol_.boundaryFieldRef().set(
                patchi,
                fvPatchField<scalar>::New(
                    bcMap[rhoBCtype], mesh_.boundary()[patchi], bVol_));
        }
    }
}

void Foam::DiscreteVelocityPoint::initBoundaryField()
{
    volScalarField::Boundary &
        hBCs = hVol_.boundaryFieldRef();
    volScalarField::Boundary &
        bBCs = bVol_.boundaryFieldRef();

    forAll(hBCs, patchi)
    {
        if (hBCs[patchi].type() == "fixedValue")
        {
            DiscreteMaxwell(
                hBCs[patchi],
                bBCs[patchi],
                dvm_.rhoVol().boundaryFieldRef()[patchi],
                dvm_.Uvol().boundaryFieldRef()[patchi],
                dvm_.lambdaVol().boundaryFieldRef()[patchi]);
        }
    }
}

void Foam::DiscreteVelocityPoint::Reconstruction()
{
    const volVectorField &C = mesh_.C();
    // Boundary faces
    forAll(hVol_.boundaryField(), patchi)
    {
        word type = hVol_.boundaryField()[patchi].type();
        fvPatchField<scalar> &hVolPatch = hVol_.boundaryFieldRef()[patchi];
        fvPatchField<scalar> &bVolPatch = bVol_.boundaryFieldRef()[patchi];
        const fvsPatchField<vector> &CfPatch = mesh_.Cf().boundaryField()[patchi];
        const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();

        if (type == "zeroGradient")
        {
            hVolPatch = hVolPatch.patchInternalField();
            bVolPatch = bVolPatch.patchInternalField();
        }
        else if (type == "calculated")
        {
            forAll(hVolPatch, pfacei)
            {
                scalar own = pOwner[pfacei];
                hVolPatch[pfacei] = hVol_[own] + (hGradVol_[own] & (CfPatch[pfacei] - C[own]));
                bVolPatch[pfacei] = bVol_[own] + (bGradVol_[own] & (CfPatch[pfacei] - C[own]));
            }
        }
    }
    hGradVol_ = fvc::grad(hVol_);
    bGradVol_ = fvc::grad(bVol_);
}

void Foam::DiscreteVelocityPoint::Update(scalar h, scalar b, label celli)
{
    hVol_[celli] = h;
    bVol_[celli] = b;
}

surfaceScalarField &Foam::DiscreteVelocityPoint::UpdatehFlux()
{
    return hFluxSurf_;
}

surfaceScalarField &Foam::DiscreteVelocityPoint::UpdatebFlux()
{
    return bFluxSurf_;
}

template <template <class> class PatchType, class GeoMesh>
void Foam::DiscreteVelocityPoint::DiscreteMaxwell(
    GeometricField<scalar, PatchType, GeoMesh> &hEq,
    GeometricField<scalar, PatchType, GeoMesh> &bEq,
    const GeometricField<scalar, PatchType, GeoMesh> &rho,
    const GeometricField<vector, PatchType, GeoMesh> &U,
    const GeometricField<scalar, PatchType, GeoMesh> &lambda)
{
    label D = mesh_.nSolutionD();
    label K = dvm_.KInner();

    hEq = rho * pow(sqrt(lambda / pi), D) * exp(-lambda * magSqr(U - xi_));
    bEq = (K + 3.0 - D) * hEq / (2.0 * lambda);
}

void Foam::DiscreteVelocityPoint::DiscreteMaxwell(
    Foam::fvPatchScalarField &hEq,
    Foam::fvPatchScalarField &bEq,
    const Foam::fvPatchScalarField &rho,
    const Foam::fvPatchVectorField &U,
    const Foam::fvPatchScalarField &lambda)
{
    label D = mesh_.nSolutionD();
    label K = dvm_.KInner();

    hEq = rho * pow(sqrt(lambda / pi), D) * exp(-lambda * magSqr(U - xi_));
    bEq = (K + 3.0 - D) * hEq / (2.0 * lambda);
}
// ************************************************************************* //
