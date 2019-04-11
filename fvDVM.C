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

\*---------------------------------------------------------------------------*/

#include "fvDVM.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(fvDVM, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::fvDVM::initialiseDV(scalarList &weights1D, scalarList &Xis1D)
{
    // Check scalarList consistency
    if (Xis1D.size() != nXiPerDim_ || weights1D.size() != nXiPerDim_)
    {
        WarningIn("initialiseDV")
            << "Num of discreteVelocityPoint not consistent" << endl;
        std::exit(-1);
    }
    else if (Xis1D[0] != xiMin_ || Xis1D[nXiPerDim_ - 1] != xiMax_)
    {
        WarningIn("initialiseDV")
            << "xi not consistant" << endl;
        std::exit(-1);
    }

    scalarField weightsGlobal;
    vectorField XisGlobal;
    labelField symmXtgID;
    labelField symmYtgID;
    labelField symmZtgID;

    if (mesh_.nSolutionD() == 3) //3D(X & Y & Z)
    {
        nXiX_ = nXiY_ = nXiZ_ = nXiPerDim_;
        nXi_ = nXiX_ * nXiY_ * nXiZ_;

        weightsGlobal.setSize(nXi_);
        XisGlobal.setSize(nXi_);
        symmXtgID.setSize(nXi_);
        symmYtgID.setSize(nXi_);
        symmZtgID.setSize(nXi_);

        label i = 0;
        for (label iz = 0; iz < nXiZ_; iz++)
        {
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiZ_; ix++)
                {
                    weightsGlobal[i] = weights1D[iz] * weights1D[iy] * weights1D[ix];
                    vector xi(Xis1D[ix], Xis1D[iy], Xis1D[iz]);
                    XisGlobal[i] = xi;
                    symmXtgID[i] = iz * nXiY_ * nXiX_ + iy * nXiX_ + (nXiX_ - ix - 1);
                    symmYtgID[i] = iz * nXiY_ * nXiX_ + (nXiY_ - iy - 1) * nXiX_ + ix;
                    symmZtgID[i] = (nXiZ_ - iz - 1) * nXiY_ * nXiX_ + iy * nXiX_ + ix;
                    i++;
                }
            }
        }
    }
    else if (mesh_.nSolutionD() == 2) //2D (X & Y)
    {
        nXiX_ = nXiY_ = nXiPerDim_;
        nXiZ_ = 1;
        nXi_ = nXiX_ * nXiY_ * nXiZ_;
        weightsGlobal.setSize(nXi_);
        XisGlobal.setSize(nXi_);
        symmXtgID.setSize(nXi_);
        symmYtgID.setSize(nXi_);
        symmZtgID.setSize(nXi_);
        label i = 0;
        for (label iy = 0; iy < nXiY_; iy++)
        {
            for (label ix = 0; ix < nXiX_; ix++)
            {
                weightsGlobal[i] = weights1D[iy] * weights1D[ix] * 1.0;
                vector xi(Xis1D[ix], Xis1D[iy], 0.0);
                XisGlobal[i] = xi;
                symmXtgID[i] = iy * nXiX_ + (nXiX_ - ix - 1);
                symmYtgID[i] = (nXiY_ - iy - 1) * nXiX_ + ix;
                symmZtgID[i] = 0;
                i++;
            }
        }
    }
    else //1D (X)
    {
        nXiX_ = nXiPerDim_;
        nXiY_ = nXiZ_ = 1;
        nXi_ = nXiX_ * nXiY_ * nXiZ_;
        weightsGlobal.setSize(nXi_);
        XisGlobal.setSize(nXi_);
        symmXtgID.setSize(nXi_);
        symmYtgID.setSize(nXi_);
        symmZtgID.setSize(nXi_);
        label i = 0;
        for (label ix = 0; ix < nXiX_; ix++)
        {
            weightsGlobal[i] = weights1D[ix] * 1.0 * 1.0;
            vector xi(Xis1D[ix], 0.0, 0.0);
            XisGlobal[i] = xi;
            symmXtgID[i] = (nXiX_ - ix - 1);
            symmYtgID[i] = 0;
            symmZtgID[i] = 0;
            i++;
        }
    }

    DV_.setSize(nXi_);
    forAll(DV_, i)
    {
        DV_.set(
            i,
            new DiscreteVelocityPoint(
                *this,
                mesh_,
                time_,
                weightsGlobal[i],
                XisGlobal[i],
                i,
                symmXtgID[i],
                symmYtgID[i],
                symmZtgID[i]));
    }
}

void Foam::fvDVM::SetLocalFrame()
{
    const surfaceVectorField nf(mesh_.Sf() / mesh_.magSf());
    // Internal faces
    forAll(nf, facei)
    {
        LocalFrameSurf_[facei] = GramSchmidtProcess(nf[facei]);
    }
    // Boundary faces
    forAll(LocalFrameSurf_.boundaryField(), patchi)
    {
        fvsPatchField<tensor> &LocalFrameSurfPatch =
            LocalFrameSurf_.boundaryFieldRef()[patchi];
        const fvsPatchField<vector> &nfPatch =
            nf.boundaryField()[patchi];

        forAll(LocalFrameSurfPatch, pFacei)
        {
            LocalFrameSurfPatch[pFacei] = GramSchmidtProcess(nfPatch[pFacei]);
        }
    }
}

void Foam::fvDVM::SetVolPro()
{
    const surfaceVectorField &ssf = mesh_.Sf();
    const labelUList &owner = mesh_.owner();
    const labelUList &neighbour = mesh_.neighbour();

    forAll(owner, facei)
    {
        VolPro_[owner[facei]] += cmptMag(ssf[facei]);
        VolPro_[neighbour[facei]] += cmptMag(ssf[facei]);
    }

    forAll(mesh_.boundary(), patchi)
    {
        const labelUList &pFaceCells =
            mesh_.boundary()[patchi].faceCells();

        const fvsPatchField<vector> &pssf = ssf.boundaryField()[patchi];

        forAll(mesh_.boundary()[patchi], facei)
        {
            VolPro_[pFaceCells[facei]] += cmptMag(pssf[facei]);
        }
    }
    VolPro_ = VolPro_ * 0.5;
}

tensor Foam::fvDVM::GramSchmidtProcess(const vector a)
{
    vector t;
    if ((std::abs(a.x()) < std::abs(a.y())) && (std::abs(a.x()) < std::abs(a.z())))
    {
        t = vector(1, 0, 0);
        t = t - (t & a) * a;
    }
    else if (std::abs(a.y()) < std::abs(a.z()))
    {
        t = vector(0, 1, 0);
        t = t - (t & a) * a;
    }
    else
    {
        t = vector(0, 0, 1);
        t = t - (t & a) * a;
    }

    if (mag(t) < SMALL)
    {
        FatalErrorInFunction
            << "axis1, axis2 appear to be co-linear: "
            << a << ", " << t << endl
            << abort(FatalError);
    }

    const vector b = t / mag(t);
    const vector c = a ^ b;
    return tensor(a, b, c);
}

void Foam::fvDVM::Reconstruction()
{
    // Prepare date for processor boundary
    rhoVol_.correctBoundaryConditions();
    Uvol_.correctBoundaryConditions();
    lambdaVol_.correctBoundaryConditions();

    // Boundary faces
    forAll(rhoVol_.boundaryField(), patchi)
    {
        word type = rhoVol_.boundaryField()[patchi].type();
        fvPatchField<scalar> &rhoVolPatch = rhoVol_.boundaryFieldRef()[patchi];
        fvPatchField<scalar> &lambdaVolPatch = lambdaVol_.boundaryFieldRef()[patchi];
        const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();

        if (type == "calculated")
        {
            forAll(rhoVolPatch, pfacei)
            {
                const label own = pOwner[pfacei];
                rhoVolPatch[pfacei] = rhoVol()[own] / lambdaVol()[own] * lambdaVolPatch[pfacei]; // Isobaric in boundary cell
            }
        }
    }

    // Calculate gradient
    forAll(DV_, DVid)
    {
        DV_[DVid].Reconstruction();
    }
    rhoGradVol_ = fvc::grad(rhoVol_);
    rhoUgradVol_ = fvc::grad(rhoVol_ * Uvol_);
    rhoEgradVol_ = fvc::grad(0.5 * rhoVol_ * (magSqr(Uvol_) + (KInner() + 3) / (2.0 * lambdaVol_)));

    // Prepare date for processor boundary
    rhoGradVol_.correctBoundaryConditions();
    rhoUgradVol_.correctBoundaryConditions();
    rhoEgradVol_.correctBoundaryConditions();
}

/*---------------------------------------------------------------------------*\
 * @Brief       : Reconstruct primL primR at interface in global frame
 *                according to global order with blow-up guard
\*---------------------------------------------------------------------------*/
void Foam::fvDVM::CheckReconstruction()
{
    // Global geometric information
    const cellList &cells = mesh_.cells();
    const labelUList &owner = mesh_.faceOwner();
    const surfaceVectorField &Cf = mesh_.Cf();
    const volVectorField &C = mesh_.C();

    forAll(cells, celli)
    {
        // Auxiliary variables
        scalar rho;
        vector rhoU;
        scalar rhoE;
        label myOrder = 1;

        if (time_.timeIndex() > firstOrderSteps)
            myOrder = orderGlobal;

        // Cell-faces information
        const cell myface = cells[celli];

        // Obtain conVars
        PrimToConserved(rhoVol_[celli], Uvol_[celli], lambdaVol_[celli], rho, rhoU, rhoE);

        forAll(myface, index)
        {
            if (myOrder == 2)
            {
                const label facei = myface[index];
                const label own = owner[facei];

                // This gives you the patch ID for a face. If the face is internal, you get -1
                const label patchID = mesh_.boundaryMesh().whichPatch(facei);
                if (patchID == -1)
                { // Internal
                    const vector d = Cf[facei] - C[celli];
                    if (own == celli)
                    {
                        myOrder = CheckInterfacePrim(rho, rhoU, rhoE, celli, d,
                                                     L_rhoSurf_[facei], L_Usurf_[facei], L_lambdaSurf_[facei]);
                    }
                    else
                    {
                        myOrder = CheckInterfacePrim(rho, rhoU, rhoE, celli, d,
                                                     R_rhoSurf_[facei], R_Usurf_[facei], R_lambdaSurf_[facei]);
                    }
                }
                else if (mesh_.boundaryMesh().types()[patchID] != "empty") // Boundary: exclude empty cases
                {
                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(facei);
                    const vector d = Cf.boundaryField()[patchID][faceID] - C[celli];
                    myOrder = CheckInterfacePrim(rho, rhoU, rhoE, celli, d,
                                                 L_rhoSurf_.boundaryFieldRef()[patchID][faceID],
                                                 L_Usurf_.boundaryFieldRef()[patchID][faceID],
                                                 L_lambdaSurf_.boundaryFieldRef()[patchID][faceID]);
                }
            }
        }

        if (myOrder != 2)
        {
            forAll(myface, index)
            {
                const label facei = myface[index];
                const label own = owner[facei];

                // This gives you the patch ID for a face. If the face is internal, you get -1
                const label patchID = mesh_.boundaryMesh().whichPatch(facei);
                if (patchID == -1)
                { // Internal
                    if (own == celli)
                    {
                        L_rhoSurf_[facei] = rhoVol_[celli];
                        L_Usurf_[facei] = Uvol_[celli];
                        L_lambdaSurf_[facei] = lambdaVol_[celli];
                    }
                    else
                    {
                        R_rhoSurf_[facei] = rhoVol_[celli];
                        R_Usurf_[facei] = Uvol_[celli];
                        R_lambdaSurf_[facei] = lambdaVol_[celli];
                    }
                }
                else if (mesh_.boundaryMesh().types()[patchID] != "empty")
                { // Boundary: exclude empty cases
                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(facei);
                    L_rhoSurf_.boundaryFieldRef()[patchID][faceID] = rhoVol_[celli];
                    L_Usurf_.boundaryFieldRef()[patchID][faceID] = Uvol_[celli];
                    L_lambdaSurf_.boundaryFieldRef()[patchID][faceID] = lambdaVol_[celli];
                }
            }
        }
        cellOrderVol_[celli] = static_cast<scalar>(myOrder);
    }
    cellOrderVol_.correctBoundaryConditions();
}

void Foam::fvDVM::CalcFluxSurf()
{
    // Init Flux to zero
    rhoFluxSurf_ = dimensionedScalar("0", rhoFluxSurf_.dimensions(), 0);
    rhoUfluxSurf_ = dimensionedVector("0", rhoUfluxSurf_.dimensions(), vector(0, 0, 0));
    rhoEfluxSurf_ = dimensionedScalar("0", rhoEfluxSurf_.dimensions(), 0);

    // Global geometric information
    const labelUList &owner = mesh_.owner();
    const labelUList &neighbour = mesh_.neighbour();
    const surfaceScalarField &magSf = mesh_.magSf();
    const surfaceVectorField &Cf = mesh_.Cf();
    const volVectorField &C = mesh_.C();
    const surfaceTensorField &localFrame = LocalFrameSurf();
    const scalar dt = time_.deltaTValue();

    // Internal faces first
    forAll(owner, facei)
    {
        // Local geometric information
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const tensor frame = localFrame[facei];
        const label ownOrder = round(cellOrderVol_[own]);
        const label neiOrder = round(cellOrderVol_[nei]);
        label faceOrder = orderGlobal;

        // Face variable
        scalarField h(nXi()), b(nXi());
        vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
        scalar lambda, rho, rhoE;
        vector U, rhoU;

        // Macro slope of g in local frame and time coefficients
        scalarField aBar(5), bBar(5), cBar(5), aT(5), Mt(5);

        // Prim^L and Prim^R in local frame at cell interface
        scalarField primL = VariablesToField(L_rhoSurf_[facei], frame & L_Usurf_[facei], L_lambdaSurf_[facei]);
        scalarField primR = VariablesToField(R_rhoSurf_[facei], frame & R_Usurf_[facei], R_lambdaSurf_[facei]);

        // Obtain the conserved variables in local frame at cell interface by collision
        scalarField prim(5);
        if (faceOrder == 2)
        {
            rho = 0.0;
            rhoU = vector(0, 0, 0);
            rhoE = 0.0;
            if (ownOrder == 2 && neiOrder == 2)
            {
                vector dl = Cf[facei] - C[own];
                vector dr = Cf[facei] - C[nei];
                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                              dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                              frame, dv.xi(), dv.weight(),
                                              dv.hVol()[nei], dv.hGradVol()[nei], dv.bVol()[nei], dv.bGradVol()[nei], dr);
                }
                faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                if (faceOrder == 2)
                {
                    calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                  prim, frame, aBar, bBar, cBar, aT,
                                  primR, rhoGradVol_[nei], rhoUgradVol_[nei], rhoEgradVol_[nei]);
                }
            }
            else if (ownOrder == 2 && neiOrder != 2)
            {
                vector dl = Cf[facei] - C[own];
                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                              dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                              frame, dv.xi(), dv.weight(),
                                              dv.hVol()[nei], dv.bVol()[nei]);
                }
                faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                if (faceOrder == 2)
                {
                    calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                  prim, frame, aBar, bBar, cBar, aT);
                }
            }
            else if (ownOrder != 2 && neiOrder == 2)
            {
                vector dr = Cf[facei] - C[nei];
                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                              dv.hVol()[own], dv.bVol()[own],
                                              frame, dv.xi(), dv.weight(),
                                              dv.hVol()[nei], dv.hGradVol()[nei], dv.bVol()[nei], dv.bGradVol()[nei], dr);
                }
                faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                if (faceOrder == 2)
                {
                    calc_g0_slope(prim, frame, aBar, bBar, cBar, aT,
                                  primR, rhoGradVol_[nei], rhoUgradVol_[nei], rhoEgradVol_[nei]);
                }
            }
            else
            {
                faceOrder = 1;
            }
        }
        if (faceOrder != 2)
        {
            rho = 0.0;
            rhoU = vector(0, 0, 0);
            rhoE = 0.0;
            forAll(DV_, dvi)
            {
                DiscreteVelocityPoint &dv = DV_[dvi];
                moment_accumulator_first(rho, rhoU, rhoE, h[dvi], b[dvi], ui[dvi],
                                         dv.hVol()[own], dv.bVol()[own],
                                         frame, dv.xi(), dv.weight(),
                                         dv.hVol()[nei], dv.bVol()[nei]);
            }
            prim = ConservedToPrim(rho, rhoU, rhoE);
            if (prim[0] < SMALL || prim[4] < SMALL)
                Pout << "Error: moment_accumulator_first failed at internal face!" << tab
                     << prim[0] << tab << prim[4] << tab << Cf[facei] << nl;
        }

        // Calculate collision time and some time integration terms
        Mt = GetTimeIntegration(prim, dt, primL, primR);
        FieldToVariables(prim, rho, U, lambda);

        // Get heatflux at interface
        vector qf = GetHeatFlux(U, DV_, h, b, ui);

        if (faceOrder == 2)
        {
            forAll(DV_, dvi)
            {
                DiscreteVelocityPoint &dv = DV_[dvi];
                flux_calculator_second(h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                       qf, DV_[dvi].weight(), magSf[facei], Mt,
                                       prim, aBar, bBar, cBar, aT,
                                       dv.UpdatehFlux()[facei], dv.UpdatebFlux()[facei],
                                       rhoFluxSurf_[facei], rhoUfluxSurf_[facei], rhoEfluxSurf_[facei]);
            }
        }
        else
        {
            forAll(DV_, dvi)
            {
                DiscreteVelocityPoint &dv = DV_[dvi];
                flux_calculator_first(h[dvi], b[dvi], ui[dvi],
                                      qf, DV_[dvi].weight(), magSf[facei], Mt, prim,
                                      dv.UpdatehFlux()[facei], dv.UpdatebFlux()[facei],
                                      rhoFluxSurf_[facei], rhoUfluxSurf_[facei], rhoEfluxSurf_[facei]);
            }
        }
        rhoUfluxSurf_[facei] = frame.T() & rhoUfluxSurf_[facei]; // Momentum flux in global frame
        rhoEfluxSurf_[facei] = 0.5 * rhoEfluxSurf_[facei];       // Add the 0.5 back in energy flux
    }

    // boundary faces
    forAll(rhoVol_.boundaryField(), patchi)
    {
        word type = rhoVol_.boundaryField()[patchi].type();
        fvsPatchField<scalar> &rhoFluxPatch = rhoFluxSurf_.boundaryFieldRef()[patchi];
        fvsPatchField<vector> &rhoUfluxPatch = rhoUfluxSurf_.boundaryFieldRef()[patchi];
        fvsPatchField<scalar> &rhoEfluxPatch = rhoEfluxSurf_.boundaryFieldRef()[patchi];

        // Prim^L at boundary interface
        const fvsPatchField<scalar> &L_rhoPatch = L_rhoSurf_.boundaryField()[patchi];
        const fvsPatchField<vector> &L_Upatch = L_Usurf_.boundaryField()[patchi];
        const fvsPatchField<scalar> &L_lambdaPatch = L_lambdaSurf_.boundaryField()[patchi];

        // Obtain fixed value at boundary
        const fvPatchField<scalar> &rhoPatch = rhoVol_.boundaryField()[patchi];
        const fvPatchField<vector> &Upatch = Uvol_.boundaryField()[patchi];
        const fvPatchField<scalar> &lambdaPatch = lambdaVol_.boundaryField()[patchi];

        // Patch geometric information
        const fvsPatchField<scalar> &magSfPatch = magSf.boundaryField()[patchi];
        const fvsPatchField<vector> &CfPatch = Cf.boundaryField()[patchi];
        const fvsPatchField<tensor> &FramePatch = localFrame.boundaryField()[patchi];
        const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();

        if (type == "calculated")
        {
            rhoFluxPatch = 0.0;
            rhoUfluxPatch = vector(0, 0, 0);
            rhoEfluxPatch = 0.0;
            const label &D = mesh_.nSolutionD();

            // Wall data output initilization
            fileName outputDir;
            if (time_.outputTime() && pOwner.size() != 0)
            {
                if (Pstream::parRun())
                {
                    outputDir = time_.path() / ".." / "wallData/";
                }
                else
                {
                    outputDir = time_.path() / "wallData/";
                }
                if (!isDir(outputDir))
                {
                    mkDir(outputDir);
                }
                outputDir = outputDir + rhoPatch.patch().name() + "_" + time_.timeName();
                if (Pstream::parRun())
                {
                    outputDir = outputDir + "_" + std::to_string(Pstream::myProcNo());
                }
                outputDir = outputDir + ".dat";
            }
            OFstream wallFile(outputDir);

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                const tensor frame = FramePatch[pFacei];
                const label ownOrder = round(cellOrderVol_[own]);

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField dh(nXi(), vector(0, 0, 0)), db(nXi(), vector(0, 0, 0)); // Micro gradient in local frame
                scalarField Hg(nXi()), Bg(nXi());
                vectorField ui(nXi());

                // Wall primary variables at interface
                scalar rho_w;
                vector U_w = frame & Upatch[pFacei]; // Macro velocity in local frame
                scalar lambda_w = lambdaPatch[pFacei];

                // Prim^L at interface
                scalar rho_g = L_rhoPatch[pFacei];
                vector U_g = frame & L_Upatch[pFacei]; // Macro velocity in local frame
                scalar lambda_g = L_lambdaPatch[pFacei];

                // Calculate collision time and some time integration terms
                scalar tau = GetTau(rho_g, lambda_g);
                scalarField Mt = GetTimeIntegration(tau, dt);

                // Get the incidence and reflection density flux
                scalar incidence = 0.0, reflection = 0.0;
                if (ownOrder == 2)
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        const scalar &weight = dv.weight();
                        const vector xii = frame & dv.xi(); // Micro velocity in local frame
                        const scalar &vn = xii.x();
                        ui[dvi] = xii;

                        if (vn >= VSMALL) // Comming from own
                        {
                            h[dvi] = dv.hVol()[own] + (dv.hGradVol()[own] & (CfPatch[pFacei] - C[own]));
                            b[dvi] = dv.bVol()[own] + (dv.bGradVol()[own] & (CfPatch[pFacei] - C[own]));
                            dh[dvi] = frame & dv.hGradVol()[own];
                            db[dvi] = frame & dv.bGradVol()[own];
                            DiscreteMaxwell(Hg[dvi], Bg[dvi], rho_g, U_g, lambda_g, xii);
                            incidence += weight * vn * (Mt[0] * Hg[dvi] + Mt[3] * h[dvi] - Mt[4] * (xii & dh[dvi]));
                        }
                        else // Comming form nei
                        {
                            reflection += dt * weight * vn * pow(sqrt(lambda_w / pi), D) * exp(-lambda_w * magSqr(U_w - xii));
                        }
                    }
                }
                else
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        const scalar &weight = dv.weight();
                        const vector xii = frame & dv.xi(); // Micro velocity in local frame
                        const scalar &vn = xii.x();
                        ui[dvi] = xii;

                        if (vn >= VSMALL) // Comming from own
                        {
                            h[dvi] = dv.hVol()[own];
                            b[dvi] = dv.bVol()[own];
                            DiscreteMaxwell(Hg[dvi], Bg[dvi], rho_g, U_g, lambda_g, xii);
                            incidence += weight * vn * (Mt[0] * Hg[dvi] + Mt[3] * h[dvi]);
                        }
                        else // Comming form nei
                        {
                            reflection += dt * weight * vn * pow(sqrt(lambda_w / pi), D) * exp(-lambda_w * magSqr(U_w - xii));
                        }
                    }
                }
                rho_w = -incidence / reflection;
                if (rho_w < 0.0)
                    Pout << "Error: negative wall density!" << tab << rho_w << tab << CfPatch[pFacei] << nl;

                // Update the micro and macro flux, combined first and second order flux
                scalar my_magSf = magSfPatch[pFacei];
                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const scalar &weight = dv.weight();
                    const vector &xii = ui[dvi];
                    const scalar &vn = xii.x();
                    scalar hFlux, bFlux;

                    if (vn >= VSMALL) // Comming from own
                    {
                        hFlux = vn * (Mt[0] * Hg[dvi] + Mt[3] * h[dvi] - Mt[4] * (xii & dh[dvi])) * my_magSf;
                        bFlux = vn * (Mt[0] * Bg[dvi] + Mt[3] * b[dvi] - Mt[4] * (xii & db[dvi])) * my_magSf;

                        dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = hFlux;
                        dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = bFlux;
                    }
                    else // Comming form nei
                    {
                        scalar H_w, B_w;
                        DiscreteMaxwell(H_w, B_w, rho_w, U_w, lambda_w, xii);
                        scalar temp = dt * vn * my_magSf;
                        hFlux = temp * H_w;
                        bFlux = temp * B_w;

                        dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = hFlux;
                        dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = bFlux;
                    }
                    scalar tmpwh = weight * hFlux;
                    scalar tmpwb = weight * bFlux;

                    rhoFluxPatch[pFacei] += tmpwh;
                    rhoUfluxPatch[pFacei] += xii * tmpwh;
                    rhoEfluxPatch[pFacei] += magSqr(xii) * tmpwh + tmpwb; // lack of 0.5, will locate after the loop
                }
                rhoEfluxPatch[pFacei] = 0.5 * rhoEfluxPatch[pFacei]; // Add the 0.5 back in energy flux

                if (time_.outputTime())
                {
                    //- Wall heat flux and stress, only defined at wall boundary faces
                    //- meanningless at other kind of faces
                    scalar theta_ = -atan2(CfPatch[pFacei].y(), CfPatch[pFacei].x()) * 180.0 / pi + 180;
                    scalar normalizedFactor = dt * my_magSf;
                    scalar qWall_ = rhoEfluxPatch[pFacei] / normalizedFactor;
                    scalar normalStressWall_ = rhoUfluxPatch[pFacei].x() / normalizedFactor;
                    scalar shearStressWall_ = (rhoUfluxPatch[pFacei].y() + rhoUfluxPatch[pFacei].z()) / normalizedFactor;
                    wallFile << theta_ << tab << qWall_ << tab << normalStressWall_ << tab << shearStressWall_ << endl;
                    wallFile.flush();
                }
                rhoUfluxPatch[pFacei] = frame.T() & rhoUfluxPatch[pFacei]; // Momentum flux in global frame
            }
        }
        else if (type == "fixedValue")
        {
            rhoFluxPatch = 0.0;
            rhoUfluxPatch = vector(0, 0, 0);
            rhoEfluxPatch = 0.0;

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                const tensor frame = FramePatch[pFacei];
                label faceOrder = round(cellOrderVol_[own]);

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
                scalar lambda, rho, rhoE;
                vector U, rhoU;

                // Macro slope of g in local frame and time coefficients
                scalarField aBar(5), bBar(5), cBar(5), aT(5), Mt(5);

                // Prim^L and Prim^R in local frame at cell interface
                scalarField primL = VariablesToField(L_rhoPatch[pFacei], frame & L_Upatch[pFacei], L_lambdaPatch[pFacei]);
                scalarField primR = VariablesToField(rhoPatch[pFacei], frame & Upatch[pFacei], lambdaPatch[pFacei]);

                // Obtain the conserved variables in local frame at cell interface by collision
                scalarField prim(5);
                if (faceOrder == 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;

                    vector dl = CfPatch[pFacei] - C[own];
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                                  dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                                  frame, dv.xi(), dv.weight(),
                                                  dv.hVol().boundaryField()[patchi][pFacei], dv.bVol().boundaryField()[patchi][pFacei]);
                    }
                    faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                    if (faceOrder == 2)
                    {
                        calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                      prim, frame, aBar, bBar, cBar, aT);
                    }
                }
                if (faceOrder != 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;

                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        moment_accumulator_first(rho, rhoU, rhoE, h[dvi], b[dvi], ui[dvi],
                                                 dv.hVol()[own], dv.bVol()[own],
                                                 frame, dv.xi(), dv.weight(),
                                                 dv.hVol().boundaryField()[patchi][pFacei], dv.bVol().boundaryField()[patchi][pFacei]);
                    }
                    prim = ConservedToPrim(rho, rhoU, rhoE);
                    if (prim[0] < SMALL || prim[4] < SMALL)
                        Pout << "Error: moment_accumulator_first failed at fixedValue patch face!" << tab
                             << prim[0] << tab << prim[4] << tab << CfPatch[pFacei] << nl;
                }

                // Calculate collision time and some time integration terms
                Mt = GetTimeIntegration(prim, dt, primL, primR);
                FieldToVariables(prim, rho, U, lambda);

                // Get heatflux at interface
                vector qf = GetHeatFlux(U, DV_, h, b, ui);

                if (faceOrder == 2)
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_second(h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                               qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt,
                                               prim, aBar, bBar, cBar, aT,
                                               dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                               rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                else
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_first(h[dvi], b[dvi], ui[dvi],
                                              qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt, prim,
                                              dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                              rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                rhoUfluxPatch[pFacei] = frame.T() & rhoUfluxPatch[pFacei]; // Momentum flux in global frame
                rhoEfluxPatch[pFacei] = 0.5 * rhoEfluxPatch[pFacei];       // Add the 0.5 back in energy flux
            }
        }
        else if (type == "zeroGradient")
        {
            rhoFluxPatch = 0.0;
            rhoUfluxPatch = vector(0, 0, 0);
            rhoEfluxPatch = 0.0;

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                const tensor frame = FramePatch[pFacei];
                label faceOrder = round(cellOrderVol_[own]);

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
                scalar lambda, rho, rhoE;
                vector U, rhoU;

                // Macro slope of g in local frame and time coefficients
                scalarField aBar(5), bBar(5), cBar(5), aT(5), Mt(5);

                // Prim^L and Prim^R in local frame at cell interface
                scalarField primL = VariablesToField(L_rhoPatch[pFacei], frame & L_Upatch[pFacei], L_lambdaPatch[pFacei]);
                scalarField primR = VariablesToField(rhoVol()[own], frame & Uvol()[own], lambdaVol()[own]);

                // Obtain the conserved variables in local frame at cell interface by collision
                scalarField prim(5);
                if (faceOrder == 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;

                    vector dl = CfPatch[pFacei] - C[own];
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                                  dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                                  frame, dv.xi(), dv.weight(),
                                                  dv.hVol()[own], dv.bVol()[own]);
                    }
                    faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                    if (faceOrder == 2)
                    {
                        calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                      prim, frame, aBar, bBar, cBar, aT);
                    }
                }
                if (faceOrder != 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;

                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        moment_accumulator_first(rho, rhoU, rhoE, h[dvi], b[dvi], ui[dvi],
                                                 dv.hVol()[own], dv.bVol()[own],
                                                 frame, dv.xi(), dv.weight(),
                                                 dv.hVol()[own], dv.bVol()[own]);
                    }
                    prim = ConservedToPrim(rho, rhoU, rhoE);
                    if (prim[0] < SMALL || prim[4] < SMALL)
                        Pout << "Error: moment_accumulator_first failed at zeroGradient patch face!" << tab
                             << prim[0] << tab << prim[4] << tab << CfPatch[pFacei] << nl;
                }

                // Calculate collision time and some time integration terms
                Mt = GetTimeIntegration(prim, dt, primL, primR);
                FieldToVariables(prim, rho, U, lambda);

                // Get heatflux at interface
                vector qf = GetHeatFlux(U, DV_, h, b, ui);

                if (faceOrder == 2)
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_second(h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                               qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt,
                                               prim, aBar, bBar, cBar, aT,
                                               dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                               rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                else
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_first(h[dvi], b[dvi], ui[dvi],
                                              qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt, prim,
                                              dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                              rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                rhoUfluxPatch[pFacei] = frame.T() & rhoUfluxPatch[pFacei]; // Momentum flux in global frame
                rhoEfluxPatch[pFacei] = 0.5 * rhoEfluxPatch[pFacei];       // Add the 0.5 back in energy flux
            }
        }
        else if (rhoVol_.boundaryField()[patchi].coupled())
        {
            Field<scalar> rhoVolNei(
                rhoVol_.boundaryField()[patchi].patchNeighbourField());
            Field<vector> UvolNei(
                Uvol_.boundaryField()[patchi].patchNeighbourField());
            Field<scalar> lambdaVolNei(
                lambdaVol_.boundaryField()[patchi].patchNeighbourField());
            const Field<vector> Cnei(
                C.boundaryField()[patchi].patchNeighbourField());
            const Field<vector> rhoGradVolNei(
                rhoGradVol_.boundaryField()[patchi].patchNeighbourField());
            const Field<tensor> rhoUgradVolNei(
                rhoUgradVol_.boundaryField()[patchi].patchNeighbourField());
            const Field<vector> rhoEgradVolNei(
                rhoEgradVol_.boundaryField()[patchi].patchNeighbourField());
            const Field<scalar> cellOrderVolNei(
                cellOrderVol_.boundaryField()[patchi].patchNeighbourField());

            PtrList<scalarField> hVolNei(nXi()), bVolNei(nXi());
            PtrList<vectorField> hGradVolNei(nXi()), bGradVolNei(nXi());
            forAll(DV_, dvi)
            {
                DiscreteVelocityPoint &dv = DV_[dvi];
                hVolNei.set(dvi, dv.hVol().boundaryField()[patchi].patchNeighbourField());
                bVolNei.set(dvi, dv.bVol().boundaryField()[patchi].patchNeighbourField());
                hGradVolNei.set(dvi, dv.hGradVol().boundaryField()[patchi].patchNeighbourField());
                bGradVolNei.set(dvi, dv.bGradVol().boundaryField()[patchi].patchNeighbourField());
            }

            rhoFluxPatch = 0.0;
            rhoUfluxPatch = vector(0, 0, 0);
            rhoEfluxPatch = 0.0;

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                const tensor frame = FramePatch[pFacei];
                const label ownOrder = round(cellOrderVol_[own]);
                const label neiOrder = round(cellOrderVolNei[pFacei]);
                label faceOrder = orderGlobal;

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
                scalar lambda, rho, rhoE;
                vector U, rhoU;

                // Macro slope of g in local frame and time coefficients
                scalarField aBar(5), bBar(5), cBar(5), aT(5), Mt(5);

                // Prim^L and Prim^R in local frame at cell interface
                scalarField primL = VariablesToField(L_rhoPatch[pFacei], frame & L_Upatch[pFacei], L_lambdaPatch[pFacei]);
                scalarField primR;
                if (faceOrder == 2 && neiOrder == 2)
                {
                    vector dr = CfPatch[pFacei] - Cnei[pFacei];
                    PrimToConserved(rhoVolNei[pFacei], UvolNei[pFacei], lambdaVolNei[pFacei], rho, rhoU, rhoE);
                    rho = rho + (rhoGradVolNei[pFacei] & dr);
                    rhoU = frame & (rhoU + (rhoUgradVolNei[pFacei].T() & dr));
                    rhoE = rhoE + (rhoEgradVolNei[pFacei] & dr);
                    primR = ConservedToPrim(rho, rhoU, rhoE);
                }
                else
                {
                    primR = VariablesToField(rhoVolNei[pFacei], frame & UvolNei[pFacei], lambdaVolNei[pFacei]);
                }

                // Obtain the conserved variables in local frame at cell interface by collision
                scalarField prim(5);
                if (faceOrder == 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;
                    if (ownOrder == 2 && neiOrder == 2)
                    {
                        vector dl = CfPatch[pFacei] - C[own];
                        vector dr = CfPatch[pFacei] - Cnei[pFacei];
                        forAll(DV_, dvi)
                        {
                            DiscreteVelocityPoint &dv = DV_[dvi];
                            moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                                      dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                                      frame, dv.xi(), dv.weight(),
                                                      hVolNei[dvi][pFacei], hGradVolNei[dvi][pFacei], bVolNei[dvi][pFacei], bGradVolNei[dvi][pFacei], dr);
                        }
                        faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                        if (faceOrder == 2)
                        {
                            calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                          prim, frame, aBar, bBar, cBar, aT,
                                          primR, rhoGradVolNei[pFacei], rhoUgradVolNei[pFacei], rhoEgradVolNei[pFacei]);
                        }
                    }
                    else if (ownOrder == 2 && neiOrder != 2)
                    {
                        vector dl = CfPatch[pFacei] - C[own];
                        forAll(DV_, dvi)
                        {
                            DiscreteVelocityPoint &dv = DV_[dvi];
                            moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                                      dv.hVol()[own], dv.hGradVol()[own], dv.bVol()[own], dv.bGradVol()[own], dl,
                                                      frame, dv.xi(), dv.weight(),
                                                      hVolNei[dvi][pFacei], bVolNei[dvi][pFacei]);
                        }
                        faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                        if (faceOrder == 2)
                        {
                            calc_g0_slope(primL, rhoGradVol_[own], rhoUgradVol_[own], rhoEgradVol_[own],
                                          prim, frame, aBar, bBar, cBar, aT);
                        }
                    }
                    else if (ownOrder != 2 && neiOrder == 2)
                    {
                        vector dr = CfPatch[pFacei] - Cnei[pFacei];
                        forAll(DV_, dvi)
                        {
                            DiscreteVelocityPoint &dv = DV_[dvi];
                            moment_accumulator_second(rho, rhoU, rhoE, h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                                      dv.hVol()[own], dv.bVol()[own],
                                                      frame, dv.xi(), dv.weight(),
                                                      hVolNei[dvi][pFacei], hGradVolNei[dvi][pFacei], bVolNei[dvi][pFacei], bGradVolNei[dvi][pFacei], dr);
                        }
                        faceOrder = ConservedToPrim(rho, rhoU, rhoE, prim);
                        if (faceOrder == 2)
                        {
                            calc_g0_slope(prim, frame, aBar, bBar, cBar, aT,
                                          primR, rhoGradVolNei[pFacei], rhoUgradVolNei[pFacei], rhoEgradVolNei[pFacei]);
                        }
                    }
                    else
                    {
                        faceOrder = 1;
                    }
                }
                if (faceOrder != 2)
                {
                    rho = 0.0;
                    rhoU = vector(0, 0, 0);
                    rhoE = 0.0;
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        moment_accumulator_first(rho, rhoU, rhoE, h[dvi], b[dvi], ui[dvi],
                                                 dv.hVol()[own], dv.bVol()[own],
                                                 frame, dv.xi(), dv.weight(),
                                                 hVolNei[dvi][pFacei], bVolNei[dvi][pFacei]);
                    }
                    prim = ConservedToPrim(rho, rhoU, rhoE);
                    if (prim[0] < SMALL || prim[4] < SMALL)
                        Pout << "Error: moment_accumulator_first failed at coupled patch face!" << tab
                             << prim[0] << tab << prim[4] << tab << CfPatch[pFacei] << nl;
                }

                // Calculate collision time and some time integration terms
                Mt = GetTimeIntegration(prim, dt, primL, primR);
                FieldToVariables(prim, rho, U, lambda);

                // Get heatflux at interface
                vector qf = GetHeatFlux(U, DV_, h, b, ui);

                if (faceOrder == 2)
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_second(h[dvi], dh[dvi], b[dvi], db[dvi], ui[dvi],
                                               qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt,
                                               prim, aBar, bBar, cBar, aT,
                                               dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                               rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                else
                {
                    forAll(DV_, dvi)
                    {
                        DiscreteVelocityPoint &dv = DV_[dvi];
                        flux_calculator_first(h[dvi], b[dvi], ui[dvi],
                                              qf, DV_[dvi].weight(), magSfPatch[pFacei], Mt, prim,
                                              dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei], dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei],
                                              rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);
                    }
                }
                rhoUfluxPatch[pFacei] = frame.T() & rhoUfluxPatch[pFacei]; // Momentum flux in global frame
                rhoEfluxPatch[pFacei] = 0.5 * rhoEfluxPatch[pFacei];       // Add the 0.5 back in energy flux
            }
        }
    }
}

void Foam::fvDVM::Update()
{
    // Global geometric information
    const cellList &cells = mesh_.cells();
    const labelUList &owner = mesh_.faceOwner();
    const scalarField &Vol = mesh_.cellVolumes();
    const scalar dt = time_.deltaTValue();
    const volVectorField &C = mesh_.C();

    forAll(cells, celli)
    {
        // Store prim^n
        scalar rhoOld = rhoVol_[celli];
        vector Uold = Uvol_[celli];
        scalar lambdaOld = lambdaVol_[celli];
        scalar tauOld = GetTau(rhoOld, lambdaOld);

        // Auxiliary variables
        scalar rho;
        vector rhoU;
        scalar rhoE;

        // Obtain conVars^n
        PrimToConserved(rhoOld, Uold, lambdaOld, rho, rhoU, rhoE);

        // Cell-faces information
        const cell myface = cells[celli];
        const scalar V = Vol[celli];
        forAll(myface, index)
        {
            const label facei = myface[index];
            const label own = owner[facei];

            // This gives you the patch ID for a face. If the face is internal, you get -1
            const label patchID = mesh_.boundaryMesh().whichPatch(facei);
            if (patchID == -1)
            { // Internal
                if (own == celli)
                {
                    rho -= rhoFluxSurf_[facei] / V;
                    rhoU -= rhoUfluxSurf_[facei] / V;
                    rhoE -= rhoEfluxSurf_[facei] / V;
                }
                else
                {
                    rho += rhoFluxSurf_[facei] / V;
                    rhoU += rhoUfluxSurf_[facei] / V;
                    rhoE += rhoEfluxSurf_[facei] / V;
                }
            }
            else if (mesh_.boundaryMesh().types()[patchID] != "empty") // Boundary: exclude empty cases
            {
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(facei);
                rho -= rhoFluxSurf_.boundaryField()[patchID][faceID] / V;
                rhoU -= rhoUfluxSurf_.boundaryField()[patchID][faceID] / V;
                rhoE -= rhoEfluxSurf_.boundaryField()[patchID][faceID] / V;
            }
        }

        // Update prim^{n+1}
        ConservedToPrim(rho, rhoU, rhoE, rhoVol_[celli], Uvol_[celli], lambdaVol_[celli]);
        if (rhoVol_[celli] <= 0.0 || lambdaVol_[celli] <= 0.0)
            Pout << "Error: negative density or lambda after update!"
                 << tab << rhoVol_[celli] << tab << lambdaVol_[celli]
                 << tab << C[celli] << nl;

        // Calculate heat flux at t=t^n
        vector qf = vector(0, 0, 0);
        scalar pressure = 0.0;
        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
            vector c = dv.xi() - Uold;
            scalar temp = (magSqr(c) * dv.hVol()[celli] + dv.bVol()[celli]);
            qf += 0.5 * dv.weight() * c * temp;
            pressure += (1.0 / 3.0) * dv.weight() * temp;
        }
        qVol_[celli] = qf;
        pVol_[celli] = pressure;

        scalar tau = GetTau(rhoVol_[celli], lambdaVol_[celli]);

        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
            scalar Hold, Bold;
            scalar H, B;
            scalar h = dv.hVol()[celli], b = dv.bVol()[celli];
            scalar timeFactor = 1.0 + 0.5 * dt / tau;
            equilibriumShakhov(Hold, Bold, rhoOld, Uold, lambdaOld, qf, dv.xi());
            equilibriumShakhov(H, B, rhoVol_[celli], Uvol_[celli], lambdaVol_[celli], qf, dv.xi());

            h = 0.5 * dt * (H / tau + (Hold - h) / tauOld);
            b = 0.5 * dt * (B / tau + (Bold - b) / tauOld);

            forAll(myface, index)
            {
                const label facei = myface[index];
                const label own = owner[facei];

                // This gives you the patch ID for a face. If the face is internal, you get -1
                const label patchID = mesh_.boundaryMesh().whichPatch(facei);
                if (patchID == -1)
                { // Internal
                    if (own == celli)
                    {
                        h -= dv.hFluxSurf()[facei] / V;
                        b -= dv.bFluxSurf()[facei] / V;
                    }
                    else
                    {
                        h += dv.hFluxSurf()[facei] / V;
                        b += dv.bFluxSurf()[facei] / V;
                    }
                }
                else if (mesh_.boundaryMesh().types()[patchID] != "empty")
                { // Boundary: exclude empty cases
                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(facei);
                    h -= dv.hFluxSurf().boundaryField()[patchID][faceID] / V;
                    b -= dv.bFluxSurf().boundaryField()[patchID][faceID] / V;
                }
            }
            dv.Update((dv.hVol()[celli] + h) / timeFactor, (dv.bVol()[celli] + b) / timeFactor, celli);
        }
    }
}

void Foam::fvDVM::PrimToConserved(scalar &rhoP, vector &U, scalar &lambda, scalar &rho, vector &rhoU, scalar &rhoE)
{
    const label &K = KInner_;
    rho = rhoP;
    rhoU = rho * U;
    rhoE = 0.5 * rho * (magSqr(U) + (K + 3) / (2.0 * lambda));
}

void Foam::fvDVM::ConservedToPrim(scalar &rhoC, vector &rhoU, scalar &rhoE, scalar &rho, vector &U, scalar &lambda)
{
    const label &K = KInner_;
    rho = rhoC;
    U = rhoU / rho;
    lambda = rho * (K + 3) / (4.0 * (rhoE - 0.5 * (rhoU & U)));
}

scalarField Foam::fvDVM::ConservedToPrim(scalar rhoC, vector rhoU, scalar rhoE)
{
    scalar rho, lambda;
    vector U;
    ConservedToPrim(rhoC, rhoU, rhoE, rho, U, lambda);
    return VariablesToField(rho, U, lambda);
}

label Foam::fvDVM::ConservedToPrim(scalar rhoC, vector rhoU, scalar rhoE, scalarField &prim)
{
    scalar rho, lambda;
    vector U;
    ConservedToPrim(rhoC, rhoU, rhoE, rho, U, lambda);
    prim = VariablesToField(rho, U, lambda);
    if (prim[0] < SMALL || prim[4] < SMALL)
        return 1;
    return 2;
}

scalarField Foam::fvDVM::VariablesToField(const scalar &head, const vector &V, const scalar &tail)
{
    scalarField s(5);
    s[0] = head;
    s[1] = V.x();
    s[2] = V.y();
    s[3] = V.z();
    s[4] = tail;
    return s;
}

void Foam::fvDVM::FieldToVariables(const scalarField &s, scalar &head, vector &V, scalar &tail)
{
    head = s[0];
    V.x() = s[1];
    V.y() = s[2];
    V.z() = s[3];
    tail = s[4];
}

void Foam::fvDVM::MicroSlope(const scalar &srho, const vector &srhoU, const scalar &srhoE, const scalar &rho, const vector &U, const scalar &lambda, scalar &grho, vector &gU, scalar &glambda)
{
    const label &K = KInner_;
    const vector dU = (srhoU - U * srho) / rho;
    const scalar dE = (srhoE - 0.5 * (magSqr(U) + (K + 3) / (2.0 * lambda)) * srho) / rho;

    glambda = 4.0 * lambda * lambda / (K + 3) * (2.0 * dE - 2.0 * (U & dU));
    gU = 2.0 * lambda * dU - U * glambda;
    grho = srho / rho - (U & gU) - 0.5 * (magSqr(U) + (K + 3) / (2.0 * lambda)) * glambda;
}

scalarField Foam::fvDVM::MicroSlope(const scalarField &slope, const scalarField &prim)
{
    scalar srho = slope[0];
    vector srhoU(slope[1], slope[2], slope[3]);
    scalar srhoE = slope[4];
    scalar rho = prim[0];
    vector U(prim[1], prim[2], prim[3]);
    scalar lambda = prim[4];
    scalar grho;
    vector gU;
    scalar glambda;
    MicroSlope(srho, srhoU, srhoE, rho, U, lambda, grho, gU, glambda);
    return VariablesToField(grho, gU, glambda);
}

void Foam::fvDVM::CalcMoment(const scalarField &prim, scalarField &Mu, scalarField &Mv1, scalarField &Mv2,
                             scalarField &Mxi, scalarField &MuL, scalarField &MuR)
{
    const label &K = KInner_;
    const label &D = mesh_.nSolutionD();
    label innerK = K + 3 - D;

    scalar rho, lambda;
    vector U;
    FieldToVariables(prim, rho, U, lambda);

    // Moments of normal velocity
    MuL[0] = 0.5 * erfc(-sqrt(lambda) * U.x());
    MuL[1] = U.x() * MuL[0] + 0.5 * exp(-lambda * U.x() * U.x()) / sqrt(pi * lambda);
    MuR[0] = 0.5 * erfc(sqrt(lambda) * U.x());
    MuR[1] = U.x() * MuR[0] - 0.5 * exp(-lambda * U.x() * U.x()) / sqrt(pi * lambda);
    for (label i = 2; i < (MuL).size(); ++i)
    {
        MuL[i] = U.x() * MuL[i - 1] + 0.5 * (i - 1) * MuL[i - 2] / lambda;
        MuR[i] = U.x() * MuR[i - 1] + 0.5 * (i - 1) * MuR[i - 2] / lambda;
    }
    Mu = MuL + MuR;

    // Moments of tangential velocity
    Mv1[0] = 1.0;
    Mv1[1] = U.y();
    Mv2[0] = 1.0;
    Mv2[1] = U.z();
    for (label i = 2; i < (Mv1).size(); ++i)
    {
        Mv1[i] = U.y() * Mv1[i - 1] + 0.5 * (i - 1) * Mv1[i - 2] / lambda;
        Mv2[i] = U.z() * Mv2[i - 1] + 0.5 * (i - 1) * Mv2[i - 2] / lambda;
    }

    // Moments of \xi
    Mxi[0] = 1.0;                                               // <\xi^0>
    Mxi[1] = 0.5 * innerK / lambda;                             // <\xi^2>
    Mxi[2] = innerK * (innerK + 2.0) / (4.0 * lambda * lambda); // <\xi^4>
}

scalarField Foam::fvDVM::Moment_uv1v2xi(const scalarField &Mu, const scalarField &Mv1, const scalarField &Mv2, const scalarField &Mxi,
                                        label alpha, label beta, label gamma, label delta)
{
    scalarField uv1v2xi(5);
    uv1v2xi[0] = Mu[alpha] * Mv1[beta] * Mv2[gamma] * Mxi[delta / 2];
    uv1v2xi[1] = Mu[alpha + 1] * Mv1[beta] * Mv2[gamma] * Mxi[delta / 2];
    uv1v2xi[2] = Mu[alpha] * Mv1[beta + 1] * Mv2[gamma] * Mxi[delta / 2];
    uv1v2xi[3] = Mu[alpha] * Mv1[beta] * Mv2[gamma + 1] * Mxi[delta / 2];
    uv1v2xi[4] = 0.5 * (Mu[alpha + 2] * Mv1[beta] * Mv2[gamma] * Mxi[delta / 2] +
                        Mu[alpha] * Mv1[beta + 2] * Mv2[gamma] * Mxi[delta / 2] +
                        Mu[alpha] * Mv1[beta] * Mv2[gamma + 2] * Mxi[delta / 2] +
                        Mu[alpha] * Mv1[beta] * Mv2[gamma] * Mxi[(delta + 2) / 2]);
    return uv1v2xi;
}

scalarField Foam::fvDVM::Moment_auv1v2xi(const scalarField &s, const scalarField &Mu, const scalarField &Mv1, const scalarField &Mv2,
                                         const scalarField &Mxi, label alpha, label beta, label gamma)
{
    scalarField auv1v2xi(5);
    auv1v2xi = s[0] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 0, 0) +
               s[1] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 1, beta + 0, gamma + 0, 0) +
               s[2] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 1, gamma + 0, 0) +
               s[3] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 1, 0) +
               0.5 * s[4] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 2, beta + 0, gamma + 0, 0) +
               0.5 * s[4] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 2, gamma + 0, 0) +
               0.5 * s[4] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 2, 0) +
               0.5 * s[4] * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 0, 2);
    return auv1v2xi;
}

scalar Foam::fvDVM::GetTau(const scalar &rho, const scalar &lambda)
{
    return muRef_ * 2.0 * exp((1 - omega_) * log(lambda)) / rho;
}

void Foam::fvDVM::equilibriumShakhov(
    scalar &hEq,
    scalar &bEq,
    const scalar &rho,
    const vector &U,
    const scalar &lambda,
    const vector &q,
    const vector &xi)
{
    const scalar &Pr = Pr_;
    const label &D = mesh_.nSolutionD();
    const label &K = KInner_;

    vector c = xi - U;
    scalar cSqrByRT = magSqr(c) * 2.0 * lambda;
    scalar cqBy5pRT = 0.8 * (1.0 - Pr) * pow(lambda, 2) * (c & q) / rho;
    scalar gEqBGK = rho * pow(sqrt(lambda / pi), D) * exp(-cSqrByRT / 2.0);

    hEq = (1.0 + cqBy5pRT * (cSqrByRT + K - D - 2.0)) * gEqBGK;
    bEq = (1.0 + cqBy5pRT * (cSqrByRT + K - D)) * gEqBGK * (K + 3.0 - D) / (2.0 * lambda);
}

void Foam::fvDVM::ShakhovPart(
    const scalar &H,
    const scalar &B,
    const scalar &rho,
    const vector &U,
    const scalar &lambda,
    const vector &q,
    const vector &xi,
    scalar &Hplus,
    scalar &Bplus)
{
    const scalar &Pr = Pr_;
    const label &D = mesh_.nSolutionD();
    const label &K = KInner_;

    vector c = xi - U;
    scalar cSqrByRT = magSqr(c) * 2.0 * lambda;
    scalar cqBy5pRT = 0.8 * (1.0 - Pr) * pow(lambda, 2) * (c & q) / rho;

    Hplus = cqBy5pRT * (cSqrByRT + K - D - 2.0) * H;
    Bplus = cqBy5pRT * (cSqrByRT + K - D) * B;
}

void Foam::fvDVM::DiscreteMaxwell(
    scalar &hEq,
    scalar &bEq,
    const scalar &rho,
    const vector &U,
    const scalar &lambda,
    const vector &xi)
{
    const label &D = mesh_.nSolutionD();
    const label &K = KInner_;

    hEq = rho * pow(sqrt(lambda / pi), D) * exp(-lambda * magSqr(xi - U));
    bEq = (K + 3.0 - D) * hEq / (2.0 * lambda);
}

tmp<scalarField> Foam::fvDVM::getSoundSpeed()
{
    tmp<scalarField> tvf(sqrt(0.5 * Gamma() / lambdaVol()));
    return tvf;
}

label Foam::fvDVM::CheckInterfacePrim(const scalar &rho, const vector &rhoU, const scalar &rhoE, const label &celli, const vector d,
                                      scalar &rho_face, vector &U_face, scalar &lambda_face)
{
    scalar rho_temp = rho + (rhoGradVol_[celli] & d);
    vector rhoU_temp = rhoU + (rhoUgradVol_[celli].T() & d);
    scalar rhoE_temp = rhoE + (rhoEgradVol_[celli] & d);
    ConservedToPrim(rho_temp, rhoU_temp, rhoE_temp, rho_face, U_face, lambda_face);
    if (rho_face < SMALL || lambda_face < SMALL)
        return 1;
    return orderGlobal;
}

scalarField Foam::fvDVM::GetTimeIntegration(const scalarField &prim, const scalar &dt, const scalarField &primL, const scalarField &primR)
{
    scalar tau = GetTau(prim[0], prim[4]) + artificialViscosity(primL, primR) * dt; // Add additional dissipation in case of strong shock
    return GetTimeIntegration(tau, dt);
}

scalarField Foam::fvDVM::GetTimeIntegration(const scalar &tau, const scalar &dt)
{
    scalarField Mt(5);

    scalar rt = tau / dt;
    if (rt > 1.0e3)
    {
        scalar x[7];
        x[0] = 1.0;
        x[1] = dt / tau;
        for (label xi = 2; xi <= 6; xi++)
            x[xi] = x[xi - 1] * x[1];

        Mt[0] = x[1] / 2.0 - x[2] / 6.0 + x[3] / 24.0 - x[4] / 120.0 + x[5] / 720.0 - x[6] / 5040.0;
        Mt[1] = -x[1] / 6.0 + x[2] / 12.0 - x[3] / 40.0 + x[4] / 180.0 - x[5] / 1008.0 + x[6] / 6720.0;
        Mt[2] = +x[1] / 6.0 - x[2] / 24.0 + x[3] / 120.0 - x[4] / 720.0 + x[5] / 5040.0 - x[6] / 40320.0;
        Mt[3] = 1.0 - Mt[0];
        Mt[4] = -0.5 + Mt[2] - Mt[1];

        Mt[0] = Mt[0] * dt;
        Mt[1] = Mt[1] * dt * dt;
        Mt[2] = Mt[2] * dt * dt;
        Mt[3] = Mt[3] * dt;
        Mt[4] = Mt[4] * dt * dt;
    }
    else
    {
        scalar expp = exp(-dt / tau);
        Mt[3] = tau * (1.0 - expp);
        Mt[4] = -tau * dt * expp + tau * Mt[3];
        Mt[0] = dt - Mt[3];
        Mt[1] = -tau * Mt[0] + Mt[4];
        Mt[2] = 0.5 * dt * dt - tau * Mt[0];
    }
    return Mt;
}

scalar Foam::fvDVM::artificialViscosity(const scalarField &primL, const scalarField &primR)
{
    return std::abs(primL[0] / primL[4] - primR[0] / primR[4]) / std::abs(primL[0] / primL[4] + primR[0] / primR[4]);
}

void Foam::fvDVM::moment_accumulator_first(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, scalar &b, vector &u,
                                           const scalar &hl, const scalar &bl,
                                           const tensor &frame, const vector &xi, const scalar &weight,
                                           const scalar &hr, const scalar &br)
{
    u = frame & xi;      // Micro velocity in local frame
    if (u.x() >= VSMALL) // Comming from own
    {
        h = hl;
        b = bl;
    }
    else if (u.x() < -VSMALL) // Comming form nei
    {
        h = hr;
        b = br;
    }
    else
    {
        h = 0.5 * (hl + hr);
        b = 0.5 * (bl + br);
    }

    // Obtain the conserved variables in local frame by collision
    moment_accumulator(rho, rhoU, rhoE, h, b, u, weight);
}

void Foam::fvDVM::moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                            const scalar &hl, const vector &dhl, const scalar &bl, const vector &dbl, const vector &dl,
                                            const tensor &frame, const vector &xi, const scalar &weight,
                                            const scalar &hr, const vector &dhr, const scalar &br, const vector &dbr, const vector &dr)
{
    u = frame & xi; // Micro velocity in local frame

    if (u.x() >= VSMALL) // Comming from own
    {
        h = hl + (dhl & dl);
        b = bl + (dbl & dl);
        dh = frame & dhl;
        db = frame & dbl;
    }
    else if (u.x() < -VSMALL) // Comming form nei
    {
        h = hr + (dhr & dr);
        b = br + (dbr & dr);
        dh = frame & dhr;
        db = frame & dbr;
    }
    else
    {
        h = 0.5 * (hl + (dhl & dl) + hr + (dhr & dr));
        b = 0.5 * (bl + (dbl & dl) + br + (dbr & dr));
        dh = 0.5 * frame & (dhl + dhr);
        db = 0.5 * frame & (dbl + dbr);
    }
    // Obtain the conserved variables in local frame by collision
    moment_accumulator(rho, rhoU, rhoE, h, b, u, weight);
}

void Foam::fvDVM::moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                            const scalar &hl, const vector &dhl, const scalar &bl, const vector &dbl, const vector &dl,
                                            const tensor &frame, const vector &xi, const scalar &weight,
                                            const scalar &hr, const scalar &br)
{
    u = frame & xi; // Micro velocity in local frame

    if (u.x() >= VSMALL) // Comming from own
    {
        h = hl + (dhl & dl);
        b = bl + (dbl & dl);
        dh = frame & dhl;
        db = frame & dbl;
    }
    else if (u.x() < -VSMALL) // Comming form nei
    {
        h = hr;
        b = br;
        dh = Zero;
        db = Zero;
    }
    else
    {
        h = 0.5 * (hl + (dhl & dl) + hr);
        b = 0.5 * (bl + (dbl & dl) + br);
        dh = 0.5 * frame & dhl;
        db = 0.5 * frame & dbl;
    }
    // Obtain the conserved variables in local frame by collision
    moment_accumulator(rho, rhoU, rhoE, h, b, u, weight);
}

void Foam::fvDVM::moment_accumulator_second(scalar &rho, vector &rhoU, scalar &rhoE, scalar &h, vector &dh, scalar &b, vector &db, vector &u,
                                            const scalar &hl, const scalar &bl,
                                            const tensor &frame, const vector &xi, const scalar &weight,
                                            const scalar &hr, const vector &dhr, const scalar &br, const vector &dbr, const vector &dr)
{
    u = frame & xi; // Micro velocity in local frame

    if (u.x() >= VSMALL) // Comming from own
    {
        h = hl;
        b = bl;
        dh = Zero;
        db = Zero;
    }
    else if (u.x() < -VSMALL) // Comming form nei
    {
        h = hr + (dhr & dr);
        b = br + (dbr & dr);
        dh = frame & dhr;
        db = frame & dbr;
    }
    else
    {
        h = 0.5 * (hl + hr + (dhr & dr));
        b = 0.5 * (bl + br + (dbr & dr));
        dh = 0.5 * frame & dhr;
        db = 0.5 * frame & dbr;
    }
    // Obtain the conserved variables in local frame by collision
    moment_accumulator(rho, rhoU, rhoE, h, b, u, weight);
}

void Foam::fvDVM::moment_accumulator(scalar &rho, vector &rhoU, scalar &rhoE, const scalar &h, const scalar &b, const vector &u, const scalar &weight)
{
    rho += weight * h;
    rhoU += weight * u * h; // local frame
    rhoE += 0.5 * weight * (magSqr(u) * h + b);
}

void Foam::fvDVM::calc_g0_slope(const scalarField &primL, const vector &rhoGrad_L, const tensor &rhoUgrad_L, const vector &rhoEgrad_L,
                                const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                                const scalarField &primR, const vector &rhoGrad_R, const tensor &rhoUgrad_R, const vector &rhoEgrad_R)
{
    // Moment
    scalarField Mu_L(7), MuL_L(7), MuR_L(7), Mv1_L(6), Mv2_L(6), Mxi_L(3);
    scalarField Mu_R(7), MuL_R(7), MuR_R(7), Mv1_R(6), Mv2_R(6), Mxi_R(3);
    CalcMoment(primL, Mu_L, Mv1_L, Mv2_L, Mxi_L, MuL_L, MuR_L);
    CalcMoment(primR, Mu_R, Mv1_R, Mv2_R, Mxi_R, MuL_R, MuR_R);

    // Obtain the normal and tangential gradient of conserved variables in local frame by collision
    vector tempVector1 = frame & rhoGrad_L;
    tensor tempTensor = frame & rhoUgrad_L & frame.T();
    vector tempVector2 = frame & rhoEgrad_L;
    scalarField N = VariablesToField(tempVector1.x(), tempTensor.x(), tempVector2.x());
    scalarField T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
    scalarField T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
    scalarField aL = MicroSlope(N, primL);
    scalarField bL = MicroSlope(T1, primL);
    scalarField cL = MicroSlope(T2, primL);

    tempVector1 = frame & rhoGrad_R;
    tempTensor = frame & rhoUgrad_R & frame.T();
    tempVector2 = frame & rhoEgrad_R;
    N = VariablesToField(tempVector1.x(), tempTensor.x(), tempVector2.x());
    T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
    T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
    scalarField aR = MicroSlope(N, primR);
    scalarField bR = MicroSlope(T1, primR);
    scalarField cR = MicroSlope(T2, primR);

    N = primL[0] * Moment_auv1v2xi(aL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(aR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
    T1 = primL[0] * Moment_auv1v2xi(bL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(bR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
    T2 = primL[0] * Moment_auv1v2xi(cL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(cR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);

    get_abcT_slope(prim, aBar, bBar, cBar, aT, N, T1, T2);
}

void Foam::fvDVM::calc_g0_slope(const scalarField &primL, const vector &rhoGrad_L, const tensor &rhoUgrad_L, const vector &rhoEgrad_L,
                                const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT)
{
    // Moment
    scalarField Mu_L(7), MuL_L(7), MuR_L(7), Mv1_L(6), Mv2_L(6), Mxi_L(3);
    CalcMoment(primL, Mu_L, Mv1_L, Mv2_L, Mxi_L, MuL_L, MuR_L);

    // Obtain the normal and tangential gradient of conserved variables in local frame by collision
    vector tempVector1 = frame & rhoGrad_L;
    tensor tempTensor = frame & rhoUgrad_L & frame.T();
    vector tempVector2 = frame & rhoEgrad_L;
    scalarField N = VariablesToField(tempVector1.x(), tempTensor.x(), tempVector2.x());
    scalarField T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
    scalarField T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
    scalarField aL = MicroSlope(N, primL);
    scalarField bL = MicroSlope(T1, primL);
    scalarField cL = MicroSlope(T2, primL);

    N = primL[0] * Moment_auv1v2xi(aL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0);
    T1 = primL[0] * Moment_auv1v2xi(bL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0);
    T2 = primL[0] * Moment_auv1v2xi(cL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0);

    get_abcT_slope(prim, aBar, bBar, cBar, aT, N, T1, T2);
}

void Foam::fvDVM::calc_g0_slope(const scalarField &prim, const tensor &frame, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                                const scalarField &primR, const vector &rhoGrad_R, const tensor &rhoUgrad_R, const vector &rhoEgrad_R)
{
    // Moment
    scalarField Mu_R(7), MuL_R(7), MuR_R(7), Mv1_R(6), Mv2_R(6), Mxi_R(3);
    CalcMoment(primR, Mu_R, Mv1_R, Mv2_R, Mxi_R, MuL_R, MuR_R);

    // Obtain the normal and tangential gradient of conserved variables in local frame by collision
    vector tempVector1 = frame & rhoGrad_R;
    tensor tempTensor = frame & rhoUgrad_R & frame.T();
    vector tempVector2 = frame & rhoEgrad_R;
    scalarField N = VariablesToField(tempVector1.x(), tempTensor.x(), tempVector2.x());
    scalarField T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
    scalarField T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
    scalarField aR = MicroSlope(N, primR);
    scalarField bR = MicroSlope(T1, primR);
    scalarField cR = MicroSlope(T2, primR);

    N = primR[0] * Moment_auv1v2xi(aR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
    T1 = primR[0] * Moment_auv1v2xi(bR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
    T2 = primR[0] * Moment_auv1v2xi(cR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);

    get_abcT_slope(prim, aBar, bBar, cBar, aT, N, T1, T2);
}

void Foam::fvDVM::get_abcT_slope(const scalarField &prim, scalarField &aBar, scalarField &bBar, scalarField &cBar, scalarField &aT,
                                 const scalarField &N, const scalarField &T1, const scalarField &T2)
{
    // Moment
    scalarField Mu(7), MuL(7), MuR(7), Mv1(6), Mv2(6), Mxi(3);
    CalcMoment(prim, Mu, Mv1, Mv2, Mxi, MuL, MuR);

    // Calculate aBar, bBar, cBar
    aBar = MicroSlope(N, prim);
    bBar = MicroSlope(T1, prim);
    cBar = MicroSlope(T2, prim);

    // Calculate aT
    scalarField Mau = Moment_auv1v2xi(aBar, Mu, Mv1, Mv2, Mxi, 1, 0, 0);
    scalarField Mbv1 = Moment_auv1v2xi(bBar, Mu, Mv1, Mv2, Mxi, 0, 1, 0);
    scalarField Mcv2 = Moment_auv1v2xi(cBar, Mu, Mv1, Mv2, Mxi, 0, 0, 1);
    scalarField sw(-prim[0] * (Mau + Mbv1 + Mcv2));
    aT = MicroSlope(sw, prim);
}

vector Foam::fvDVM::GetHeatFlux(const vector &U, const PtrList<DiscreteVelocityPoint> &DV_, const scalarField &h, const scalarField &b, const vectorField &u)
{
    vector qf = vector(0, 0, 0);
    forAll(DV_, dvi)
    {
        vector c = u[dvi] - U;
        qf += 0.5 * DV_[dvi].weight() * c * (magSqr(c) * h[dvi] + b[dvi]);
    }
    return qf;
}

void Foam::fvDVM::flux_calculator_second(const scalar &hi, const vector &dhi, const scalar &bi, const vector &dbi, const vector &xii,
                                         const vector &qf, const scalar &weight, const scalar &magSf, const scalarField &Mt,
                                         const scalarField &prim, const scalarField &aBar, const scalarField &bBar, const scalarField &cBar, const scalarField &aT,
                                         scalar &hFlux, scalar &bFlux, scalar &rhoFlux, vector &rhoUflux, scalar &rhoEflux)
{
    // Auxiliary variables of slope
    scalar rho, aRho, bRho, cRho, aTrho;
    vector U, aU, bU, cU, aTU;
    scalar lambda, aLambda, bLambda, cLambda, aTlambda;
    FieldToVariables(prim, rho, U, lambda);
    FieldToVariables(aBar, aRho, aU, aLambda);
    FieldToVariables(bBar, bRho, bU, bLambda);
    FieldToVariables(cBar, cRho, cU, cLambda);
    FieldToVariables(aT, aTrho, aTU, aTlambda);

    // <\xi^4>
    label innerK = KInner_ + 3 - mesh_.nSolutionD();
    scalar E4 = innerK * (innerK + 2.0) / (4.0 * lambda * lambda);

    const scalar &vn = xii.x();
    scalar H0, B0, Hplus, Bplus;
    DiscreteMaxwell(H0, B0, rho, U, lambda, xii);
    ShakhovPart(H0, B0, rho, U, lambda, qf, xii, Hplus, Bplus);

    scalar huuvvwwxixi = 0.5 * (magSqr(xii) * H0 + B0);
    hFlux = vn * (Mt[0] * (H0 + Hplus) + Mt[1] * (xii.x() * (aRho * H0 + (aU & xii) * H0 + aLambda * huuvvwwxixi) + xii.y() * (bRho * H0 + (bU & xii) * H0 + bLambda * huuvvwwxixi) + xii.z() * (cRho * H0 + (cU & xii) * H0 + cLambda * huuvvwwxixi)) + Mt[2] * (aTrho * H0 + (aTU & xii) * H0 + aTlambda * huuvvwwxixi) + Mt[3] * hi - Mt[4] * (xii & dhi)) * magSf;

    scalar buuvvwwxixi = 0.5 * (magSqr(xii) * B0 + E4 * H0);
    bFlux = vn * (Mt[0] * (B0 + Bplus) + Mt[1] * (xii.x() * (aRho * B0 + (aU & xii) * B0 + aLambda * buuvvwwxixi) + xii.y() * (bRho * B0 + (bU & xii) * B0 + bLambda * buuvvwwxixi) + xii.z() * (cRho * B0 + (cU & xii) * B0 + cLambda * buuvvwwxixi)) + Mt[2] * (aTrho * B0 + (aTU & xii) * B0 + aTlambda * buuvvwwxixi) + Mt[3] * bi - Mt[4] * (xii & dbi)) * magSf;

    scalar tmpwh = weight * hFlux;
    scalar tmpwb = weight * bFlux;

    rhoFlux += tmpwh;
    rhoUflux += tmpwh * xii;
    rhoEflux += (magSqr(xii) * tmpwh + tmpwb); // lack of 0.5, will locate after the function call
}

void Foam::fvDVM::flux_calculator_first(const scalar &hi, const scalar &bi, const vector &xii,
                                        const vector &qf, const scalar &weight, const scalar &magSf, const scalarField &Mt, const scalarField &prim,
                                        scalar &hFlux, scalar &bFlux, scalar &rhoFlux, vector &rhoUflux, scalar &rhoEflux)
{
    // Auxiliary variables of slope
    scalar rho, lambda;
    vector U;
    FieldToVariables(prim, rho, U, lambda);

    const scalar &vn = xii.x();
    scalar H0, B0, Hplus, Bplus;
    DiscreteMaxwell(H0, B0, rho, U, lambda, xii);
    ShakhovPart(H0, B0, rho, U, lambda, qf, xii, Hplus, Bplus);

    hFlux = vn * (Mt[0] * (H0 + Hplus) + Mt[3] * hi) * magSf;
    bFlux = vn * (Mt[0] * (B0 + Bplus) + Mt[3] * bi) * magSf;

    scalar tmpwh = weight * hFlux;
    scalar tmpwb = weight * bFlux;

    rhoFlux += tmpwh;
    rhoUflux += tmpwh * xii;
    rhoEflux += (magSqr(xii) * tmpwh + tmpwb); // lack of 0.5, will locate after the function call
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvDVM::fvDVM(
    volScalarField &rho,
    volVectorField &U,
    volScalarField &lambda,
    dictionary &DVMProperties)
    : mesh_(rho.mesh()),
      time_(rho.time()),
      rhoVol_(rho),
      Uvol_(U),
      lambdaVol_(lambda),
      nXiPerDim_(readLabel(DVMProperties.subDict("fvDVMparas").lookup("nDV"))),
      xiMax_(readScalar(DVMProperties.subDict("fvDVMparas").lookup("xiMax"))),
      xiMin_(readScalar(DVMProperties.subDict("fvDVMparas").lookup("xiMin"))),
      dXi_((xiMax_ - xiMin_) / (nXiPerDim_ - 1)),
      omega_(readScalar(DVMProperties.subDict("gasProperties").lookup("omega"))),
      muRef_(readScalar(DVMProperties.subDict("gasProperties").lookup("muRef"))),
      Pr_(readScalar(DVMProperties.subDict("gasProperties").lookup("Pr"))),
      KInner_(DVMProperties.subDict("gasProperties").lookupOrDefault<label>("KInner", 0)),
      Gamma_(scalar(KInner_ + 5) / scalar(KInner_ + 3)),
      orderGlobal(DVMProperties.lookupOrDefault<label>("orderGlobal", 2)),
      firstOrderSteps(DVMProperties.lookupOrDefault<label>("firstOrderSteps", 0)),
      DV_(0),
      rhoFluxSurf_(
          IOobject(
              "rhoFluxSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0)),
      rhoUfluxSurf_(
          IOobject(
              "rhoUfluxSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      rhoEfluxSurf_(
          IOobject(
              "rhoEfluxSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0)),
      rhoGradVol_(
          IOobject(
              "rhoGradVol",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", rhoVol_.dimensions() / dimLength, vector(0, 0, 0)),
          "zeroGradient"),
      rhoUgradVol_(
          IOobject(
              "rhoUgradVol",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedTensor("0", rhoVol_.dimensions() * Uvol_.dimensions() / dimLength, Zero),
          "zeroGradient"),
      rhoEgradVol_(
          IOobject(
              "rhoEgradVol",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", rhoVol_.dimensions() * Uvol_.dimensions() * Uvol_.dimensions() / dimLength, vector(0, 0, 0)),
          "zeroGradient"),
      qVol_(
          IOobject(
              "q",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      pVol_(
          IOobject(
              "P",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 0)),
      VolPro_(
          IOobject(
              "VolPro",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      LocalFrameSurf_(
          IOobject(
              "LocalFrameSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedTensor("0", dimless, I)),
      L_rhoSurf_(
          IOobject(
              "L_rhoSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 1.0)),
      R_rhoSurf_(
          IOobject(
              "R_rhoSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 1.0)),
      L_Usurf_(
          IOobject(
              "L_Usurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      R_Usurf_(
          IOobject(
              "R_Usurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      L_lambdaSurf_(
          IOobject(
              "L_lambdaSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 1.0)),
      R_lambdaSurf_(
          IOobject(
              "R_lambdaSurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 1.0)),
      cellOrderVol_(
          IOobject(
              "cellOrderVol",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedScalar("0", dimless, 2.0))
{
    scalarList weights = DVMProperties.lookup("weights");
    scalarList Xis = DVMProperties.lookup("Xis");
    initialiseDV(weights, Xis);
    SetLocalFrame();
    SetVolPro();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::getCoNum(scalar &maxCoNum, scalar &meanCoNum)
{
    scalar dt = time_.deltaTValue();
    scalarField UbyDxMacro(((cmptMag(Uvol()) + getSoundSpeed() * vector(1.0, 1.0, 1.0)) & VolPro()) / mesh_.V());
    scalarField UbyDxMicro = mesh_.deltaCoeffs() * sqrt(scalar(mesh_.nSolutionD())) * xiMax();

    label UbyDxMacro_size = UbyDxMacro.size();
    reduce<label>(UbyDxMacro_size, sumOp<label>(), Pstream::msgType(), UPstream::worldComm);
    label UbyDxMicro_size = UbyDxMicro.size();
    reduce<label>(UbyDxMicro_size, sumOp<label>(), Pstream::msgType(), UPstream::worldComm);

    meanCoNum = max(gSum(UbyDxMacro) / UbyDxMacro_size, gSum(UbyDxMicro) / UbyDxMicro_size) * dt;
    maxCoNum = max(gMax(UbyDxMacro), gMax(UbyDxMicro)) * dt;
}

void Foam::fvDVM::evolution()
{
    if (orderGlobal == 2 && time_.timeIndex() > firstOrderSteps)
        Reconstruction();
    CheckReconstruction();
    CalcFluxSurf();
    Update();
}

// ************************************************************************* //