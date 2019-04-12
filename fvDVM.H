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

    //- Ratio of specific heat
    const scalar Gamma_;

    //- indicator of Global Order, second order by default
    label orderGlobal;

    //- Use first order to calculate the initial field during the firstOrderSteps
    label firstOrderSteps;

    //- List of pointers to discrete velocities
    PtrList<DiscreteVelocityPoint> DV_;

    //- Surface marco rho flux
    surfaceScalarField rhoFluxSurf_;

    //- Surface marco momentum flux
    surfaceVectorField rhoUfluxSurf_;

    //- Surface marco energy flux
    surfaceScalarField rhoEfluxSurf_;

    //- Gradient of conserved variable rho
    volVectorField rhoGradVol_;

    //- Gradient of conserved variable momentum
    volTensorField rhoUgradVol_;

    //- Gradient of conserved variable energy
    volVectorField rhoEgradVol_;

    //- Volume heat flux vector
    volVectorField qVol_;

    //- Volume pressure
    volScalarField pVol_;

    //- Projections of the control volume on the y-z, x-z, x-y plane
    volVectorField VolPro_;

    //- Local coordinate at each face
    surfaceTensorField LocalFrameSurf_;

    //- rho left and rho right at surface
    surfaceScalarField L_rhoSurf_, R_rhoSurf_;

    //- U left and rho right at surface
    surfaceVectorField L_Usurf_, R_Usurf_;

    //- lambda left and lambda right at surface
    surfaceScalarField L_lambdaSurf_, R_lambdaSurf_;

    //- Cell order of reconstruction
    volScalarField cellOrderVol_;

    // Private Member Functions

    //- Initialise
    //  Create the discrete velocity(DV) set and Initialize
    //  each of the DV the with their discrete velocity value and weight
    void initialiseDV(scalarList &weights1D, scalarList &Xis1D);

    // Set up local coordinate
    void SetLocalFrame();
    tensor GramSchmidtProcess(const vector a);

    // Set up the projections of the control volume on the y-z, x-z, x-y plane
    void SetVolPro();

    // Reconstruct the gradient of initial distribution function
    void Reconstruction();

    // Inter-processor communication of hVol_, bVol_
    void hbVolCommunicate();

    // Check if the reconstruction is reliable and assign the primL and primR at interface
    // If the reconstructed density or temperature is negative, set the cell order to first order
    void CheckReconstruction();

    // Macro and micro flux calculation
    void CalcFluxSurf();

    // Update macro cell averaged values and distribution function
    void Update();

    void PrimToConserved(scalar &rhoP, vector &U, scalar &lambda, scalar &rho, vector &rhoU, scalar &rhoE);
    void ConservedToPrim(scalar &rhoC, vector &rhoU, scalar &rhoE, scalar &rho, vector &U, scalar &lambda);
    scalarField ConservedToPrim(scalar rhoC, vector rhoU, scalar rhoE);
    label ConservedToPrim(scalar rhoC, vector rhoU, scalar rhoE, scalarField &prim);
    scalarField VariablesToField(const scalar &head, const vector &V, const scalar &tail);
    void FieldToVariables(const scalarField &s, scalar &head, vector &V, scalar &tail);
    void MicroSlope(const scalar &srho, const vector &srhoU, const scalar &srhoE, const scalar &rho, const vector &U, const scalar &lambda, scalar &grho, vector &gU, scalar &glambda);
    scalarField MicroSlope(const scalarField &slope, const scalarField &prim);
    void CalcMoment(const scalarField &prim, scalarField &Mu, scalarField &Mv1, scalarField &Mv2,
                    scalarField &Mxi, scalarField &MuL, scalarField &MuR);
    scalarField Moment_uv1v2xi(const scalarField &Mu, const scalarField &Mv1, const scalarField &Mv2, const scalarField &Mxi,
                               label alpha, label beta, label gamma, label delta);
    scalarField Moment_auv1v2xi(const scalarField &s, const scalarField &Mu, const scalarField &Mv1, const scalarField &Mv2,
                                const scalarField &Mxi, label alpha, label beta, label gamma);
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

    tmp<scalarField> getSoundSpeed();

    label CheckInterfacePrim(const scalar &rho, const vector &rhoU, const scalar &rhoE, const label &celli, const vector d,
                             scalar &rho_face, vector &U_face, scalar &lambda_face);
    scalarField GetTimeIntegration(const scalarField &prim, const scalar &dt, const scalarField &primL, const scalarField &primR);
    scalarField GetTimeIntegration(const scalar &tau, const scalar &dt);
    scalar artificialViscosity(const scalarField &primL, const scalarField &primR);
    void moment_accumulator_first(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, scalar &b, vector &u, const scalar &hl, const scalar &bl, const tensor &frame, const vector &xi, const scalar &weight, const scalar &hr, const scalar &br);
    void moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                   const scalar &hl, const vector &dhl, const scalar &bl, const vector &dbl, const vector &dl,
                                   const tensor &frame, const vector &xi, const scalar &weight,
                                   const scalar &hr, const vector &dhr, const scalar &br, const vector &dbr, const vector &dr);
    void moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                   const scalar &hl, const vector &dhl, const scalar &bl, const vector &dbl, const vector &dl,
                                   const tensor &frame, const vector &xi, const scalar &weight,
                                   const scalar &hr, const scalar &br);
    void moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                   const scalar &hl, const scalar &bl,
                                   const tensor &frame, const vector &xi, const scalar &weight,
                                   const scalar &hr, const vector &dhr, const scalar &br, const vector &dbr, const vector &dr);
    void moment_accumulator(scalar &rho, vector &rhoU, scalar &rhoE, const scalar &h, const scalar &b, const vector &u, const scalar &weight);
    void calc_g0_slope(const scalarField &primL, const vector &rhoGrad_L, const tensor &rhoUgrad_L, const vector &rhoEgrad_L,
                       const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                       const scalarField &primR, const vector &rhoGrad_R, const tensor &rhoUgrad_R, const vector &rhoEgrad_R);
    void calc_g0_slope(const scalarField &primL, const vector &rhoGrad_L, const tensor &rhoUgrad_L, const vector &rhoEgrad_L,
                       const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT);
    void calc_g0_slope(const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                       const scalarField &primR, const vector &rhoGrad_R, const tensor &rhoUgrad_R, const vector &rhoEgrad_R);
    void get_abcT_slope(const scalarField &prim, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                        const scalarField &N, const scalarField &T1, const scalarField &T2);
    vector GetHeatFlux(const vector &U, const PtrList<DiscreteVelocityPoint> &DV_, const scalarField &h, const scalarField &b, const vectorField &u);
    void flux_calculator_second(const scalar &hi, const vector &dhi, const scalar &bi, const vector &dbi, const vector &xii,
                                const vector &qf, const scalar &weight, const scalar &magSf, const scalarField &Mt,
                                const scalarField &prim, const scalarField &aBar, const scalarField &bBar, const scalarField &cBar, const scalarField &aT,
                                scalar &hFlux, scalar &bFlux, scalar &rhoFlux, vector &rhoUflux, scalar &rhoEflux);
    void flux_calculator_first(const scalar &hi, const scalar &bi, const vector &xii,
                               const vector &qf, const scalar &weight, const scalar &magSf, const scalarField &Mt, const scalarField &prim,
                               scalar &hFlux, scalar &bFlux, scalar &rhoFlux, vector &rhoUflux, scalar &rhoEflux);

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

    //- Projections of the control volume on the y-z, x-z, x-y plane
    inline const volVectorField &VolPro() const;

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

    inline scalar Gamma() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDVMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
