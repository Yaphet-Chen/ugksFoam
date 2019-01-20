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

Class
    Foam::DiscreteVelocityPoint

Description
    Discrete Velocity Point

SourceFiles
    DiscreteVelocityPoint.C

\*---------------------------------------------------------------------------*/

#ifndef DiscreteVelocityPoint_H
#define DiscreteVelocityPoint_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvDVM;

/*---------------------------------------------------------------------------*\
                    Class DiscreteVelocityPoint Declaration
\*---------------------------------------------------------------------------*/

class DiscreteVelocityPoint
{
public:

private:

    // Private data

        //- Refence to the owner fvDVM object
        fvDVM& dvm_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the time
        const Time& time_;

        //- weight
        const scalar weight_;

        //- discrete velocity 
        const vector xi_;

        //- My DV Id
        const label myDVid_;

        //- 
        const label symXtargetDVid_;
        const label symYtargetDVid_;
        const label symZtargetDVid_;

        //- volume data
        volScalarField hVol_;

        volScalarField bVol_;

        //- gradient of g/h bar B
        volVectorField hGradVol_;

        volVectorField bGradVol_;

        //- surface micro h flux
        surfaceScalarField hFluxSurf_;

        //- surface micro b flux
        surfaceScalarField bFluxSurf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DiscreteVelocityPoint(const DiscreteVelocityPoint&);

        //- Disallow default bitwise assignment
        void operator=(const DiscreteVelocityPoint&);

        //- Initialise DF to equilibrium
        void initDFtoEq();

        //- set BC type for DF
        void setBCtype();

        void initBoundaryField();

        //- Obtain discretized Maxwellian distribution
        template <template<class> class PatchType, class GeoMesh> 
        void DiscreteMaxwell
        (
            GeometricField<scalar, PatchType, GeoMesh>& hEq,
            GeometricField<scalar, PatchType, GeoMesh>& bEq,
            const GeometricField<scalar, PatchType, GeoMesh>& rho,
            const GeometricField<vector, PatchType, GeoMesh>& U,
            const GeometricField<scalar, PatchType, GeoMesh>& lambda
        );

        //- Obtain discretized Maxwellian distribution at boundary patches
        void DiscreteMaxwell
        (
            fvPatchScalarField& hEq,
            fvPatchScalarField& bEq,
            const fvPatchScalarField&  rho,
            const fvPatchVectorField&  U,
            const fvPatchScalarField&  lambda
        );
public:

    // Constructors

        //- Construct form components
        DiscreteVelocityPoint
        (
            fvDVM& dvm,
            const fvMesh& mesh,
            const Time& time,
            const scalar weight,
            const vector xi,
            const label DVid,
            const label symXtargetDVid,
            const label symYtargetDVid,
            const label symZtargetDVid
        );

    //- Destructor
    ~DiscreteVelocityPoint();


    // Member functions

        // Edit

            //- Reconstruction the gradient of Discrete Velocity Point
            void Reconstruction();

            //- Update distribution function of Discrete Velocity Point
            void Update(scalar h, scalar b, label celli);

            surfaceScalarField& UpdatehFlux();
            surfaceScalarField& UpdatebFlux();

        // Access

            //- Return xi
            inline const vector xi() const;

            //- Return weight
            inline scalar weight() const;

            //- Return h
            inline const volScalarField& hVol() const;

            //- Return b
            inline const volScalarField& bVol() const;

            //- Return grad h
            inline const volVectorField& hGradVol() const;

            //- Return grad b
            inline const volVectorField& bGradVol() const;

            //- Return hFlux
            inline const surfaceScalarField& hFluxSurf() const;

            //- Return bFlux
            inline const surfaceScalarField& bFluxSurf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DiscreteVelocityPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
