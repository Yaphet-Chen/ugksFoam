/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::fvDVM

Description
    Define a distribution function of Finite Volume Discrete Velocity Method.

SourceFiles
    fvDVM.C

\*---------------------------------------------------------------------------*/

#ifndef fvDVM_H
#define fvDVM_H

#include "DiscreteVelocityPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class DiscreteVelocityPoint;

/*---------------------------------------------------------------------------*\
                           Class fvDVM Declaration
\*---------------------------------------------------------------------------*/

class fvDVM
{
  private:
    // Private data

    //- Reference to the mesh database
    const fvMesh &mesh_;

    //- Reference to the time database
    const Time &time_;

    //- Reference to the density field
    volScalarField &rhoVol_;

    //- Reference to the macro velocity field
    volVectorField &Uvol_;

    //- Reference to the lambda field
    volScalarField &lambdaVol_;

    //- Number of discrete velocity points in each direction
    const label nXiPerDim_;

    //- Total number of discrete velocities
    label nXi_;

    label nXiX_;

    label nXiY_;

    label nXiZ_;

    //- Maximum discrete velocity
    const scalar xiMax_;

    //- Minimum discrete velocity
    const scalar xiMin_;

    //- grid space in velocity space
    const scalar dXi_;

    //- Gas: omega
    const scalar omega_;

    //- Gas: muRef
    const scalar muRef_;

    //- Gas: Pr
    const scalar Pr_;

    //- Gas molecular innter degree of freedom
    const label KInner_;

    //- List of pointers to discrete velocities
    PtrList<DiscreteVelocityPoint> DV_;

    //- surface marco rho flux
    surfaceScalarField rhoFluxSurf_;

    //- surface marco momentum flux
    surfaceVectorField rhoUfluxSurf_;

    //- surface marco energy flux
    surfaceScalarField rhoEfluxSurf_;

    //- volume heat flux vector
    volVectorField qVol_;

    //- surface velocity, used for getCoNum
    surfaceVectorField Usurf_;

    //- Local coordinate at each face
    surfaceTensorField LocalFrameSurf_;

    // Private Member Functions

    //- Initialise
    //  Create the discrete velocity(DV) set and Initialize
    //  each of the DV the with their discrete velocity value and weight
    void initialiseDV(scalarList &weights1D, scalarList &Xis1D);

    // Set up local coordinate
    void SetLocalFrame();
    tensor GramSchmidtProcess(const vector a);

    // Reconstruct the gradient of initial distribution function
    void Reconstruction();

    // Macro and micro flux calculation
    void CalcFluxSurf();

    // Update macro cell averaged values and distribution function
    void Update();

    void PrimToConserved(scalar &rhoP, vector &U, scalar &lambda, scalar &rho, vector &rhoU, scalar &rhoE);
    void ConservedToPrim(scalar &rhoC, vector &rhoU, scalar &rhoE, scalar &rho, vector &U, scalar &lambda);
    void MicroSlope(scalar &srho, vector &srhoU, scalar &srhoE, scalar &rho, vector &U, scalar &lambda, scalar &grho, vector &gU, scalar &glambda);
    void CalcMoment(const scalar &rho, const vector &U, const scalar &lambda,
                    scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                    scalarField &MuL, scalarField &MuR);
    scalarField Moment_uv1v2xi(scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                               label alpha, label beta, label gamma, label delta);
    scalarField Moment_auv1v2xi(scalar &srho, vector &sU, scalar &sLambda,
                                scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                                label alpha, label beta, label gamma);
    scalar GetTau(const scalar &rho, const scalar &lambda);

    void equilibriumShakhov(
        scalar &hEq,
        scalar &bEq,
        const scalar &rho,
        const vector &U,
        const scalar &lambda,
        const vector &q,
        const vector &xi);

    void ShakhovPart(
        const scalar &H,
        const scalar &B,
        const scalar &rho,
        const vector &U,
        const scalar &lambda,
        const vector &q,
        const vector &xi,
        scalar &Hplus,
        scalar &Bplus);

    void DiscreteMaxwell(
        scalar &hEq,
        scalar &bEq,
        const scalar &rho,
        const vector &U,
        const scalar &lambda,
        const vector &xi);

    //- Disallow default bitwise copy construct
    fvDVM(const fvDVM &);

    //- Disallow default bitwise assignment
    void operator=(const fvDVM &);

  public:
    //- Runtime type information
    TypeName("fvDVM");

    // Constructors

    //- Construct from components
    fvDVM(
        volScalarField &rho,
        volVectorField &U,
        volScalarField &lambda,
        dictionary &DVMProperties);

    //- Destructor
    virtual ~fvDVM();

    // Member functions

    // Edit

    //- Solve DVB equation(s), call the updateXxx functions.
    void evolution();

    // Access

    //- Co num, max and mean
    void getCoNum(scalar &maxCoNum, scalar &meanCoNum);

    //- DiscreteVelocityPoint for I DiscreteVelocityPoint
    inline const DiscreteVelocityPoint &DVi(const label I) const;

    //- DiscreteVelocityPoint for DiscreteVelocityPoint Ix, Iy, Iz
    inline const DiscreteVelocityPoint &DVi(
        const label ix,
        const label iy,
        const label iz) const;

    //- Density field at cell center
    inline volScalarField &rhoVol();

    //- Lambda field at cell center
    inline volScalarField &lambdaVol() const;

    //- Macro velocity field at cell center
    inline volVectorField &Uvol() const;

    //- Heat flux field at cell center
    inline const volVectorField &qVol() const;

    //- Local coordinate system at face center
    inline const surfaceTensorField &LocalFrameSurf() const;

    //- Total number of discrete velocities
    inline label nXi() const;

    //- Number of discrete velocities in each direction
    inline label nXiPerDim() const;

    inline scalar xiMax() const;

    inline scalar xiMin() const;

    inline scalar omega() const;

    inline scalar muRef() const;

    inline scalar Pr() const;

    inline label KInner() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDVMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
