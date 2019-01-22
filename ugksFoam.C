/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Application
    ugksFoam

Description
    Unified Gas Kinetic Scheme solver.
    Author      : Yipei Chen
    Email       : ychendh@connect.ust.hk
    Date        : 2018-12-07

\*---------------------------------------------------------------------------*/

#include "fvDVM.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run() && ((Uchange > convergeTol) || (lambdaChange > convergeTol) || (rhoChange > convergeTol)))
    {
        #include "CourantNo.H"
        #include "setDeltaT.H"
        runTime++;
        dvm.evolution();

        if (runTime.timeIndex() % convergeCheckSteps == 0 && runTime.timeIndex() >= convergeCheckSteps)
        {
            runTime.write();
            Info << "Step =" << runTime.timeIndex() << nl
                 << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << nl
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh>>
                deltaLam = mag(lambda - lambdaOld);
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh>>
                deltaRho = mag(rho - rhoOld);
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh>>
                deltaU = mag(U - Uold);

            lambdaChange = gSum(deltaLam()) / gSum(lambda);
            rhoChange = gSum(deltaRho()) / gSum(rho);
            Uchange = gSum(deltaU()) / gSum(mag(U)());

            Info << "Lambda changes = " << lambdaChange << endl;
            Info << "Density     changes = " << rhoChange << endl;
            Info << "Velocity    changes = " << Uchange << nl << endl;
            lambdaOld = lambda;
            rhoOld = rho;
            Uold = U;
            T = 1.0 / lambda;
        }
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
