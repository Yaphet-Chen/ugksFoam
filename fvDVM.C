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
    forAll(DV_, DVid)
        DV_[DVid].Reconstruction();

    // Boundary faces
    forAll(rhoVol_.boundaryField(), patchi)
    {
        word type = rhoVol_.boundaryField()[patchi].type();
        fvPatchField<scalar> &rhoVolPatch = rhoVol_.boundaryFieldRef()[patchi];
        fvPatchField<vector> &UvolPatch = Uvol_.boundaryFieldRef()[patchi];
        fvPatchField<scalar> &lambdaVolPatch = lambdaVol_.boundaryFieldRef()[patchi];
        const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();

        if (type == "zeroGradient")
        {
            rhoVolPatch = rhoVolPatch.patchInternalField();
            UvolPatch = UvolPatch.patchInternalField();
            lambdaVolPatch = lambdaVolPatch.patchInternalField();
        }
        else if (type == "calculated")
        {
            forAll(rhoVolPatch, pfacei)
            {
                scalar own = pOwner[pfacei];
                rhoVolPatch[pfacei] = rhoVol()[own] / lambdaVol()[own] * lambdaVolPatch[pfacei]; // Isobaric in boundary cell
            }
        }
    }
    rhoGradVol_ = fvc::grad(rhoVol_);
    rhoUgradVol_ = fvc::grad(rhoVol_ * Uvol_);
    rhoEgradVol_ = fvc::grad(0.5 * rhoVol_ * (magSqr(Uvol_) + (KInner() + 3) / (2.0 * lambdaVol_)));
}

void Foam::fvDVM::CalcFluxSurf()
{   
    // Prepare date for processor boundary
    rhoVol_.correctBoundaryConditions();
    Uvol_.correctBoundaryConditions();
    lambdaVol_.correctBoundaryConditions();
    rhoGradVol_.correctBoundaryConditions();
    rhoUgradVol_.correctBoundaryConditions();
    rhoEgradVol_.correctBoundaryConditions();
    forAll(DV_, dvi)
    {
        DiscreteVelocityPoint &dv = DV_[dvi];
        dv.hVol_par().correctBoundaryConditions();
        dv.bVol_par().correctBoundaryConditions();
        dv.hGradVol_par().correctBoundaryConditions();
        dv.bGradVol_par().correctBoundaryConditions();
    }

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
        label own = owner[facei];
        label nei = neighbour[facei];
        tensor frame = localFrame[facei];

        // Face variable
        scalarField h(nXi()), b(nXi());
        vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
        scalar lambda, rho, rhoL, rhoR, rhoE, rhoEL, rhoER;
        vector U, rhoU, rhoUL, rhoUR;

        // Gradient and micro slope
        scalar aRho, bRho, cRho, aTrho;
        vector aU, bU, cU, aTU;
        scalar aLambda, bLambda, cLambda, aTlambda;

        // Moment
        scalarField Mu(7), MuL(7), MuR(7), Mv1(6), Mv2(6), Mxi(3);
        scalarField Mu_L(7), MuL_L(7), MuR_L(7), Mv1_L(6), Mv2_L(6), Mxi_L(3);
        scalarField Mu_R(7), MuL_R(7), MuR_R(7), Mv1_R(6), Mv2_R(6), Mxi_R(3);
        scalarField M0u(5), Mau(5), Mbv1(5), Mcv2(5), MaTu(5), Mt(5);

        PrimToConserved(rhoVol_[own], Uvol_[own], lambdaVol_[own], rhoL, rhoUL, rhoEL); // Left cell-centered conservered value
        PrimToConserved(rhoVol_[nei], Uvol_[nei], lambdaVol_[nei], rhoR, rhoUR, rhoER); // Left cell-centered conservered value

        // Obtain the conserved variables in local frame at cell interface by central difference
        scalar dL = frame.x() & (Cf[facei] - C[own]);
        scalar dR = frame.x() & (C[nei] - Cf[facei]);
        scalar d = dL + dR;
        rho = (dR * rhoL + dL * rhoR) / d;
        rhoU = frame & (dR * rhoUL + dL * rhoUR) / d;
        rhoE = (dR * rhoEL + dL * rhoER) / d;

        // Normal derivative 
        scalarField N(VariablesToField(rhoR - rhoL, frame & (rhoUR - rhoUL), rhoER - rhoEL) / d);

        // W^L in local frame at cell interface
        rhoL = rhoL + (rhoGradVol_[own] & (Cf[facei] - C[own]));
        rhoUL = frame & (rhoUL + (rhoUgradVol_[own].T() & (Cf[facei] - C[own])));
        rhoEL = rhoEL + (rhoEgradVol_[own] & (Cf[facei] - C[own]));
        scalarField primL = ConservedToPrim(rhoL, rhoUL, rhoEL);
        CalcMoment(primL, Mu_L, Mv1_L, Mv2_L, Mxi_L, MuL_L, MuR_L);
        
        // W^R in local frame at cell interface
        rhoR = rhoR + (rhoGradVol_[nei] & (Cf[facei] - C[nei]));
        rhoUR = frame & (rhoUR + (rhoUgradVol_[nei].T() & (Cf[facei] - C[nei])));
        rhoER = rhoER + (rhoEgradVol_[nei] & (Cf[facei] - C[nei]));
        scalarField primR = ConservedToPrim(rhoR, rhoUR, rhoER);
        CalcMoment(primR, Mu_R, Mv1_R, Mv2_R, Mxi_R, MuL_R, MuR_R);

        // Obtain the tangential gradient of conserved variables in local frame by collision
        vector tempVector1 = frame & rhoGradVol_[own];
        tensor tempTensor = frame & rhoUgradVol_[own] & frame.T();
        vector tempVector2 = frame & rhoEgradVol_[own];
        scalarField T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
        scalarField T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
        scalarField bL = MicroSlope(T1, primL);
        scalarField cL = MicroSlope(T2, primL);

        tempVector1 = frame & rhoGradVol_[nei];
        tempTensor = frame & rhoUgradVol_[nei] & frame.T();
        tempVector2 = frame & rhoEgradVol_[nei];
        T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
        T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
        scalarField bR = MicroSlope(T1, primR);
        scalarField cR = MicroSlope(T2, primR);

        T1 = primL[0] * Moment_auv1v2xi(bL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(bR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
        T2 = primL[0] * Moment_auv1v2xi(cL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(cR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);

        // Calculate aBar, bBar, cBar and related moment
        scalarField prim = ConservedToPrim(rho, rhoU, rhoE);
        scalarField aBar = MicroSlope(N, prim);
        scalarField bBar = MicroSlope(T1, prim);
        scalarField cBar = MicroSlope(T2, prim);
        FieldToVariables(aBar, aRho, aU, aLambda);
        FieldToVariables(bBar, bRho, bU, bLambda);
        FieldToVariables(cBar, cRho, cU, cLambda);

        CalcMoment(prim, Mu, Mv1, Mv2, Mxi, MuL, MuR);
        Mau = Moment_auv1v2xi(aBar, Mu, Mv1, Mv2, Mxi, 1, 0, 0);
        Mbv1 = Moment_auv1v2xi(bBar, Mu, Mv1, Mv2, Mxi, 0, 1, 0);
        Mcv2 = Moment_auv1v2xi(cBar, Mu, Mv1, Mv2, Mxi, 0, 0, 1);
        scalarField sw(-rho * (Mau + Mbv1 + Mcv2));
        scalarField aT = MicroSlope(sw, prim);
        FieldToVariables(aT, aTrho, aTU, aTlambda);

        // Calculate collision time and some time integration terms
        FieldToVariables(prim, rho, U, lambda);
        scalar tau = GetTau(rho, lambda) + std::abs(primL[0] / primL[4] - primR[0] / primR[4]) / std::abs(primL[0] / primL[4] + primR[0] / primR[4]) * dt; // Add additional dissipation in case of strong shock

        Mt[3] = tau * (1.0 - exp(-dt / tau));
        Mt[4] = -tau * dt * exp(-dt / tau) + tau * Mt[3];
        Mt[0] = dt - Mt[3];
        Mt[1] = -tau * Mt[0] + Mt[4];
        Mt[2] = 0.5 * dt * dt - tau * Mt[0];

        M0u = Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, 1, 0, 0, 0);
        Mau = Moment_auv1v2xi(aBar, Mu, Mv1, Mv2, Mxi, 2, 0, 0);
        Mbv1 = Moment_auv1v2xi(bBar, Mu, Mv1, Mv2, Mxi, 1, 1, 0);
        Mcv2 = Moment_auv1v2xi(cBar, Mu, Mv1, Mv2, Mxi, 1, 0, 1);
        MaTu = Moment_auv1v2xi(aT, Mu, Mv1, Mv2, Mxi, 1, 0, 0);

        scalarField flux(rho * (Mt[0] * M0u + Mt[1] * (Mau + Mbv1 + Mcv2) + Mt[2] * MaTu));
        FieldToVariables(flux, rhoFluxSurf_[facei], rhoUfluxSurf_[facei], rhoEfluxSurf_[facei]);

        vector qf = vector(0, 0, 0);
        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
            const vector xii = frame & dv.xi(); // Micro velocity in local frame
            ui[dvi] = xii;

            if (xii.x() >= VSMALL) // Comming from own
            {
                h[dvi] = dv.hVol()[own] + (dv.hGradVol()[own] & (Cf[facei] - C[own]));
                b[dvi] = dv.bVol()[own] + (dv.bGradVol()[own] & (Cf[facei] - C[own]));
                dh[dvi] = frame & dv.hGradVol()[own];
                db[dvi] = frame & dv.bGradVol()[own];
            }
            else if (xii.x() < -VSMALL) // Comming form nei
            {
                h[dvi] = dv.hVol()[nei] + (dv.hGradVol()[nei] & (Cf[facei] - C[nei]));
                b[dvi] = dv.bVol()[nei] + (dv.bGradVol()[nei] & (Cf[facei] - C[nei]));
                dh[dvi] = frame & dv.hGradVol()[nei];
                db[dvi] = frame & dv.bGradVol()[nei];
            }
            else
            {
                h[dvi] = 0.5 * (dv.hVol()[own] + (dv.hGradVol()[own] & (Cf[facei] - C[own])) + dv.hVol()[nei] + (dv.hGradVol()[nei] & (Cf[facei] - C[nei])));
                b[dvi] = 0.5 * (dv.bVol()[own] + (dv.bGradVol()[own] & (Cf[facei] - C[own])) + dv.bVol()[nei] + (dv.bGradVol()[nei] & (Cf[facei] - C[nei])));
                dh[dvi] = 0.5 * frame & (dv.hGradVol()[own] + dv.hGradVol()[nei]);
                db[dvi] = 0.5 * frame & (dv.bGradVol()[own] + dv.bGradVol()[nei]);
            }
            vector c = ui[dvi] - U;
            qf += 0.5 * dv.weight() * c * (magSqr(c) * h[dvi] + b[dvi]);
        }

        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
            const vector &xii = ui[dvi];
            const scalar &vn = xii.x();
            const scalar &weight = dv.weight();
            const scalar &hi = h[dvi];
            const scalar &bi = b[dvi];
            const vector &dhi = dh[dvi];
            const vector &dbi = db[dvi];

            scalar H0, B0, Hplus, Bplus;
            DiscreteMaxwell(H0, B0, rho, U, lambda, xii);
            ShakhovPart(H0, B0, rho, U, lambda, qf, xii, Hplus, Bplus);

            rhoFluxSurf_[facei] += weight * vn * (Mt[0] * Hplus + Mt[3] * hi - Mt[4] * (xii & dhi));
            rhoUfluxSurf_[facei] += weight * vn * xii * (Mt[0] * Hplus + Mt[3] * hi - Mt[4] * (xii & dhi));
            rhoEfluxSurf_[facei] += 0.5 * weight * vn * (Mt[0] * (magSqr(xii) * Hplus + Bplus) + Mt[3] * (magSqr(xii) * hi + bi) - Mt[4] * (magSqr(xii) * (xii & dhi) + (xii & dbi)));

            dv.UpdatehFlux()[facei] = vn * (Mt[0] * (H0 + Hplus) + Mt[1] * (xii.x() * (aRho * H0 + (aU & xii) * H0 + 0.5 * aLambda * (magSqr(xii) * H0 + B0)) + xii.y() * (bRho * H0 + (bU & xii) * H0 + 0.5 * bLambda * (magSqr(xii) * H0 + B0)) + xii.z() * (cRho * H0 + (cU & xii) * H0 + 0.5 * cLambda * (magSqr(xii) * H0 + B0))) + Mt[2] * (aTrho * H0 + (aTU & xii) * H0 + 0.5 * aTlambda * (magSqr(xii) * H0 + B0)) + Mt[3] * hi - Mt[4] * (xii & dhi)) * magSf[facei];

            dv.UpdatebFlux()[facei] = vn * (Mt[0] * (B0 + Bplus) + Mt[1] * (xii.x() * (aRho * B0 + (aU & xii) * B0 + 0.5 * aLambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + xii.y() * (bRho * B0 + (bU & xii) * B0 + 0.5 * bLambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + xii.z() * (cRho * B0 + (cU & xii) * B0 + 0.5 * cLambda * (magSqr(xii) * B0 + Mxi[2] * H0))) + Mt[2] * (aTrho * B0 + (aTU & xii) * B0 + 0.5 * aTlambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + Mt[3] * bi - Mt[4] * (xii & dbi)) * magSf[facei];
        }
        rhoFluxSurf_[facei] = rhoFluxSurf_[facei] * magSf[facei];
        rhoUfluxSurf_[facei] = (frame.T() & rhoUfluxSurf_[facei]) * magSf[facei];
        rhoEfluxSurf_[facei] = rhoEfluxSurf_[facei] * magSf[facei];
    }

    // boundary faces
    forAll(rhoVol_.boundaryField(), patchi)
    {
        word type = rhoVol_.boundaryField()[patchi].type();
        fvsPatchField<scalar> &rhoFluxPatch = rhoFluxSurf_.boundaryFieldRef()[patchi];
        fvsPatchField<vector> &rhoUfluxPatch = rhoUfluxSurf_.boundaryFieldRef()[patchi];
        fvsPatchField<scalar> &rhoEfluxPatch = rhoEfluxSurf_.boundaryFieldRef()[patchi];

        // Obtain fixed value at boundary
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

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                tensor frame = FramePatch[pFacei];

                // Face variable
                scalarField h(nXi(), 0.0), b(nXi(), 0.0);
                vectorField dh(nXi(), vector(0, 0, 0)), db(nXi(), vector(0, 0, 0)); // Micro gradient in local frame
                scalarField H_g(nXi(), 0.0), B_g(nXi(), 0.0);
                vectorField ui(nXi());

                scalar rho_w, incidence = 0.0, reflection = 0.0;
                vector U_w = frame & Upatch[pFacei]; // Macro velocity in local frame
                scalar lambda_w = lambdaPatch[pFacei];
                scalar rho_g = rhoVol()[own] / lambdaVol()[own] * lambda_w;

                scalar tau = GetTau(rho_g, lambda_w);
                scalar T3 = tau * (1.0 - exp(-dt / tau));
                scalar T4 = -tau * dt * exp(-dt / tau) + tau * T3;
                scalar T0 = dt - T3;

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
                        DiscreteMaxwell(H_g[dvi], B_g[dvi], rho_g, U_w, lambda_w, xii);
                        incidence += weight * vn * (T0 * H_g[dvi] + T3 * h[dvi] - T4 * (xii & dh[dvi]));
                    }
                    else // Comming form nei
                    {
                        reflection += dt * weight * vn * pow(sqrt(lambda_w / pi), D) * exp(-lambda_w * magSqr(U_w - xii));
                    }
                }
                rho_w = -incidence / reflection;

                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const scalar &weight = dv.weight();
                    const vector &xii = ui[dvi];
                    const scalar &vn = xii.x();

                    if (vn >= VSMALL) // Comming from own
                    {
                        rhoFluxPatch[pFacei] += weight * vn * (T0 * H_g[dvi] + T3 * h[dvi] - T4 * (xii & dh[dvi]));
                        rhoUfluxPatch[pFacei] += weight * vn * xii * (T0 * H_g[dvi] + T3 * h[dvi] - T4 * (xii & dh[dvi]));
                        rhoEfluxPatch[pFacei] += 0.5 * weight * vn * (T0 * (magSqr(xii) * H_g[dvi] + B_g[dvi]) + T3 * (magSqr(xii) * h[dvi] + b[dvi]) - T4 * (magSqr(xii) * (xii & dh[dvi]) + (xii & db[dvi])));
                        dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = vn * (T0 * H_g[dvi] + T3 * h[dvi] - T4 * (xii & dh[dvi])) * magSfPatch[pFacei];
                        dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = vn * (T0 * B_g[dvi] + T3 * b[dvi] - T4 * (xii & db[dvi])) * magSfPatch[pFacei];
                    }
                    else // Comming form nei
                    {
                        scalar H_w, B_w;
                        DiscreteMaxwell(H_w, B_w, rho_w, U_w, lambda_w, xii);
                        scalar temp = dt * weight * vn;
                        rhoFluxPatch[pFacei] += temp * H_w;
                        rhoUfluxPatch[pFacei] += temp * xii * H_w;
                        rhoEfluxPatch[pFacei] += temp * 0.5 * (magSqr(xii) * H_w + B_w);
                        dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = dt * vn * H_w * magSfPatch[pFacei];
                        dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = dt * vn * B_w * magSfPatch[pFacei];
                    }
                }
                rhoFluxPatch[pFacei] = rhoFluxPatch[pFacei] * magSfPatch[pFacei];
                rhoUfluxPatch[pFacei] = (frame.T() & rhoUfluxPatch[pFacei]) * magSfPatch[pFacei];
                rhoEfluxPatch[pFacei] = rhoEfluxPatch[pFacei] * magSfPatch[pFacei];
            }
        }
        else if (type == "fixedValue")
        {
        }
        else if (type == "zeroGradient")
        {
        }
        else if (rhoVol_.boundaryField()[patchi].coupled())
        {
            Field<scalar> rhoVolNei
            (
                rhoVol_.boundaryField()[patchi].patchNeighbourField()
            );
            Field<vector> UvolNei
            (
                Uvol_.boundaryField()[patchi].patchNeighbourField()
            );
            Field<scalar> lambdaVolNei
            (
                lambdaVol_.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<vector> Cnei
            (
                C.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<vector> rhoGradVolNei
            (
                rhoGradVol_.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<tensor> rhoUgradVolNei
            (
                rhoUgradVol_.boundaryField()[patchi].patchNeighbourField()
            );
            const Field<vector> rhoEgradVolNei
            (
                rhoEgradVol_.boundaryField()[patchi].patchNeighbourField()
            );

            PtrList<scalarField> hVolNei(nXi()), bVolNei(nXi());
            PtrList<vectorField> hGradVolNei(nXi()), bGradVolNei(nXi());
            forAll(DV_, dvi)
            {
                DiscreteVelocityPoint &dv = DV_[dvi];
                hVolNei.set(dvi, dv.hVol_par().boundaryField()[patchi].patchNeighbourField());
                bVolNei.set(dvi, dv.bVol_par().boundaryField()[patchi].patchNeighbourField());
                hGradVolNei.set(dvi, dv.hGradVol_par().boundaryField()[patchi].patchNeighbourField());
                bGradVolNei.set(dvi, dv.bGradVol_par().boundaryField()[patchi].patchNeighbourField());
            }

            rhoFluxPatch = 0.0;
            rhoUfluxPatch = vector(0, 0, 0);
            rhoEfluxPatch = 0.0;

            forAll(pOwner, pFacei)
            {
                // Local geometric information
                const label own = pOwner[pFacei];
                tensor frame = FramePatch[pFacei];

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
                scalar lambda, rho, rhoL, rhoR, rhoE, rhoEL, rhoER;
                vector U, rhoU, rhoUL, rhoUR;

                // Gradient and micro slope
                scalar aRho, bRho, cRho, aTrho;
                vector aU, bU, cU, aTU;
                scalar aLambda, bLambda, cLambda, aTlambda;

                // Moment
                scalarField Mu(7), MuL(7), MuR(7), Mv1(6), Mv2(6), Mxi(3);
                scalarField Mu_L(7), MuL_L(7), MuR_L(7), Mv1_L(6), Mv2_L(6), Mxi_L(3);
                scalarField Mu_R(7), MuL_R(7), MuR_R(7), Mv1_R(6), Mv2_R(6), Mxi_R(3);
                scalarField M0u(5), Mau(5), Mbv1(5), Mcv2(5), MaTu(5), Mt(5);

                PrimToConserved(rhoVol_[own], Uvol_[own], lambdaVol_[own], rhoL, rhoUL, rhoEL); // Left cell-centered conservered value
                PrimToConserved(rhoVolNei[pFacei], UvolNei[pFacei], lambdaVolNei[pFacei], rhoR, rhoUR, rhoER); // Left cell-centered conservered value

                // Obtain the conserved variables in local frame at cell interface by central difference
                scalar dL = frame.x() & (CfPatch[pFacei] - C[own]);
                scalar dR = frame.x() & (Cnei[pFacei] - CfPatch[pFacei]);
                scalar d = dL + dR;
                rho = (dR * rhoL + dL * rhoR) / d;
                rhoU = frame & (dR * rhoUL + dL * rhoUR) / d;
                rhoE = (dR * rhoEL + dL * rhoER) / d;

                // Normal derivative 
                scalarField N(VariablesToField(rhoR - rhoL, frame & (rhoUR - rhoUL), rhoER - rhoEL) / d);

                // W^L in local frame at cell interface
                rhoL = rhoL + (rhoGradVol_[own] & (CfPatch[pFacei] - C[own]));
                rhoUL = frame & (rhoUL + (rhoUgradVol_[own].T() & (CfPatch[pFacei] - C[own])));
                rhoEL = rhoEL + (rhoEgradVol_[own] & (CfPatch[pFacei] - C[own]));
                scalarField primL = ConservedToPrim(rhoL, rhoUL, rhoEL);
                CalcMoment(primL, Mu_L, Mv1_L, Mv2_L, Mxi_L, MuL_L, MuR_L);
            
                // W^R in local frame at cell interface
                rhoR = rhoR + (rhoGradVolNei[pFacei] & (CfPatch[pFacei] - Cnei[pFacei]));
                rhoUR = frame & (rhoUR + (rhoUgradVolNei[pFacei].T() & (CfPatch[pFacei] - Cnei[pFacei])));
                rhoER = rhoER + (rhoEgradVolNei[pFacei] & (CfPatch[pFacei] - Cnei[pFacei]));
                scalarField primR = ConservedToPrim(rhoR, rhoUR, rhoER);
                CalcMoment(primR, Mu_R, Mv1_R, Mv2_R, Mxi_R, MuL_R, MuR_R);

                // Obtain the tangential gradient of conserved variables in local frame by collision
                vector tempVector1 = frame & rhoGradVol_[own];
                tensor tempTensor = frame & rhoUgradVol_[own] & frame.T();
                vector tempVector2 = frame & rhoEgradVol_[own];
                scalarField T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
                scalarField T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
                scalarField bL = MicroSlope(T1, primL);
                scalarField cL = MicroSlope(T2, primL);

                tempVector1 = frame & rhoGradVolNei[pFacei];
                tempTensor = frame & rhoUgradVolNei[pFacei] & frame.T();
                tempVector2 = frame & rhoEgradVolNei[pFacei];
                T1 = VariablesToField(tempVector1.y(), tempTensor.y(), tempVector2.y());
                T2 = VariablesToField(tempVector1.z(), tempTensor.z(), tempVector2.z());
                scalarField bR = MicroSlope(T1, primR);
                scalarField cR = MicroSlope(T2, primR);

                T1 = primL[0] * Moment_auv1v2xi(bL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(bR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);
                T2 = primL[0] * Moment_auv1v2xi(cL, MuL_L, Mv1_L, Mv2_L, Mxi_L, 0, 0, 0) + primR[0] * Moment_auv1v2xi(cR, MuR_R, Mv1_R, Mv2_R, Mxi_R, 0, 0, 0);

                // Calculate aBar, bBar, cBar and related moment
                scalarField prim = ConservedToPrim(rho, rhoU, rhoE);
                scalarField aBar = MicroSlope(N, prim);
                scalarField bBar = MicroSlope(T1, prim);
                scalarField cBar = MicroSlope(T2, prim);
                FieldToVariables(aBar, aRho, aU, aLambda);
                FieldToVariables(bBar, bRho, bU, bLambda);
                FieldToVariables(cBar, cRho, cU, cLambda);

                CalcMoment(prim, Mu, Mv1, Mv2, Mxi, MuL, MuR);
                Mau = Moment_auv1v2xi(aBar, Mu, Mv1, Mv2, Mxi, 1, 0, 0);
                Mbv1 = Moment_auv1v2xi(bBar, Mu, Mv1, Mv2, Mxi, 0, 1, 0);
                Mcv2 = Moment_auv1v2xi(cBar, Mu, Mv1, Mv2, Mxi, 0, 0, 1);
                scalarField sw(-rho * (Mau + Mbv1 + Mcv2));
                scalarField aT = MicroSlope(sw, prim);
                FieldToVariables(aT, aTrho, aTU, aTlambda);

                // Calculate collision time and some time integration terms
                FieldToVariables(prim, rho, U, lambda);
                scalar tau = GetTau(rho, lambda) + std::abs(primL[0] / primL[4] - primR[0] / primR[4]) / std::abs(primL[0] / primL[4] + primR[0] / primR[4]) * dt; // Add additional dissipation in case of strong shock

                Mt[3] = tau * (1.0 - exp(-dt / tau));
                Mt[4] = -tau * dt * exp(-dt / tau) + tau * Mt[3];
                Mt[0] = dt - Mt[3];
                Mt[1] = -tau * Mt[0] + Mt[4];
                Mt[2] = 0.5 * dt * dt - tau * Mt[0];

                M0u = Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, 1, 0, 0, 0);
                Mau = Moment_auv1v2xi(aBar, Mu, Mv1, Mv2, Mxi, 2, 0, 0);
                Mbv1 = Moment_auv1v2xi(bBar, Mu, Mv1, Mv2, Mxi, 1, 1, 0);
                Mcv2 = Moment_auv1v2xi(cBar, Mu, Mv1, Mv2, Mxi, 1, 0, 1);
                MaTu = Moment_auv1v2xi(aT, Mu, Mv1, Mv2, Mxi, 1, 0, 0);

                scalarField flux(rho * (Mt[0] * M0u + Mt[1] * (Mau + Mbv1 + Mcv2) + Mt[2] * MaTu));
                FieldToVariables(flux, rhoFluxPatch[pFacei], rhoUfluxPatch[pFacei], rhoEfluxPatch[pFacei]);

                vector qf = vector(0, 0, 0);
                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const vector xii = frame & dv.xi(); // Micro velocity in local frame
                    ui[dvi] = xii;

                    if (xii.x() >= VSMALL) // Comming from own
                    {
                        h[dvi] = dv.hVol()[own] + (dv.hGradVol()[own] & (CfPatch[pFacei] - C[own]));
                        b[dvi] = dv.bVol()[own] + (dv.bGradVol()[own] & (CfPatch[pFacei] - C[own]));
                        dh[dvi] = frame & dv.hGradVol()[own];
                        db[dvi] = frame & dv.bGradVol()[own];
                    }
                    else if (xii.x() < -VSMALL) // Comming form nei
                    {
                        h[dvi] = hVolNei[dvi][pFacei] + (hGradVolNei[dvi][pFacei] & (CfPatch[pFacei] - Cnei[pFacei]));
                        b[dvi] = bVolNei[dvi][pFacei] + (bGradVolNei[dvi][pFacei] & (CfPatch[pFacei] - Cnei[pFacei]));
                        dh[dvi] = frame & hGradVolNei[dvi][pFacei];
                        db[dvi] = frame & bGradVolNei[dvi][pFacei];
                    }
                    else
                    {
                        h[dvi] = 0.5 * (dv.hVol()[own] + (dv.hGradVol()[own] & (CfPatch[pFacei] - C[own])) + hVolNei[dvi][pFacei] + (hGradVolNei[dvi][pFacei] & (CfPatch[pFacei] - Cnei[pFacei])));
                        b[dvi] = 0.5 * (dv.bVol()[own] + (dv.bGradVol()[own] & (CfPatch[pFacei] - C[own])) + bVolNei[dvi][pFacei] + (bGradVolNei[dvi][pFacei] & (CfPatch[pFacei] - Cnei[pFacei])));
                        dh[dvi] = 0.5 * frame & (dv.hGradVol()[own] + hGradVolNei[dvi][pFacei]);
                        db[dvi] = 0.5 * frame & (dv.bGradVol()[own] + bGradVolNei[dvi][pFacei]);
                    }
                    vector c = ui[dvi] - U;
                    qf += 0.5 * dv.weight() * c * (magSqr(c) * h[dvi] + b[dvi]);
                }

                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const vector &xii = ui[dvi];
                    const scalar &vn = xii.x();
                    const scalar &weight = dv.weight();
                    const scalar &hi = h[dvi];
                    const scalar &bi = b[dvi];
                    const vector &dhi = dh[dvi];
                    const vector &dbi = db[dvi];

                    scalar H0, B0, Hplus, Bplus;
                    DiscreteMaxwell(H0, B0, rho, U, lambda, xii);
                    ShakhovPart(H0, B0, rho, U, lambda, qf, xii, Hplus, Bplus);

                    rhoFluxPatch[pFacei] += weight * vn * (Mt[0] * Hplus + Mt[3] * hi - Mt[4] * (xii & dhi));
                    rhoUfluxPatch[pFacei] += weight * vn * xii * (Mt[0] * Hplus + Mt[3] * hi - Mt[4] * (xii & dhi));
                    rhoEfluxPatch[pFacei] += 0.5 * weight * vn * (Mt[0] * (magSqr(xii) * Hplus + Bplus) + Mt[3] * (magSqr(xii) * hi + bi) - Mt[4] * (magSqr(xii) * (xii & dhi) + (xii & dbi)));

                    dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = vn * (Mt[0] * (H0 + Hplus) + Mt[1] * (xii.x() * (aRho * H0 + (aU & xii) * H0 + 0.5 * aLambda * (magSqr(xii) * H0 + B0)) + xii.y() * (bRho * H0 + (bU & xii) * H0 + 0.5 * bLambda * (magSqr(xii) * H0 + B0)) + xii.z() * (cRho * H0 + (cU & xii) * H0 + 0.5 * cLambda * (magSqr(xii) * H0 + B0))) + Mt[2] * (aTrho * H0 + (aTU & xii) * H0 + 0.5 * aTlambda * (magSqr(xii) * H0 + B0)) + Mt[3] * hi - Mt[4] * (xii & dhi)) * magSfPatch[pFacei];

                    dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = vn * (Mt[0] * (B0 + Bplus) + Mt[1] * (xii.x() * (aRho * B0 + (aU & xii) * B0 + 0.5 * aLambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + xii.y() * (bRho * B0 + (bU & xii) * B0 + 0.5 * bLambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + xii.z() * (cRho * B0 + (cU & xii) * B0 + 0.5 * cLambda * (magSqr(xii) * B0 + Mxi[2] * H0))) + Mt[2] * (aTrho * B0 + (aTU & xii) * B0 + 0.5 * aTlambda * (magSqr(xii) * B0 + Mxi[2] * H0)) + Mt[3] * bi - Mt[4] * (xii & dbi)) * magSfPatch[pFacei];
                }
                rhoFluxPatch[pFacei] = rhoFluxPatch[pFacei] * magSfPatch[pFacei];
                rhoUfluxPatch[pFacei] = (frame.T() & rhoUfluxPatch[pFacei]) * magSfPatch[pFacei];
                rhoEfluxPatch[pFacei] = rhoEfluxPatch[pFacei] * magSfPatch[pFacei];
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
        const cell mycell = cells[celli];
        const scalar V = Vol[celli];
        forAll(mycell, index)
        {
            const label facei = mycell[index];
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

        // Calculate heat flux at t=t^n
        vector qf = vector(0, 0, 0);
        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
            vector c = dv.xi() - Uold;
            qf += 0.5 * dv.weight() * c * (magSqr(c) * dv.hVol()[celli] + dv.bVol()[celli]);
        }
        qVol_[celli] = qf;

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

            forAll(mycell, index)
            {
                const label facei = mycell[index];
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

scalarField Foam::fvDVM::VariablesToField(scalar head, vector V, scalar tail)
{
    scalarField s(5);
    s[0] = head;
    s[1] = V.x();
    s[2] = V.y();
    s[3] = V.z();
    s[4] = tail;
    return s;
}

void Foam::fvDVM::FieldToVariables(scalarField s, scalar &head, vector &V, scalar &tail)
{
    head = s[0];
    V.x() = s[1];
    V.y() = s[2];
    V.z() = s[3];
    tail = s[4];
}

void Foam::fvDVM::MicroSlope(scalar &srho, vector &srhoU, scalar &srhoE, scalar &rho, vector &U, scalar &lambda, scalar &grho, vector &gU, scalar &glambda)
{
    const label &K = KInner_;
    const vector dU = (srhoU - U * srho) / rho;
    const scalar dE = (srhoE - 0.5 * (magSqr(U) + (K + 3) / (2.0 * lambda)) * srho) / rho;

    glambda = 4.0 * lambda * lambda / (K + 3) * (2.0 * dE - 2.0 * (U & dU));
    gU = 2.0 * lambda * dU - U * glambda;
    grho = srho / rho - (U & gU) - 0.5 * (magSqr(U) + (K + 3) / (2.0 * lambda)) * glambda;
}

scalarField Foam::fvDVM::MicroSlope(scalarField &slope, scalarField &prim)
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

scalarField Foam::fvDVM::Moment_uv1v2xi(scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
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

scalarField Foam::fvDVM::Moment_auv1v2xi(scalarField &s, scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                                         label alpha, label beta, label gamma)
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
          dimensionedTensor("0", dimless, I))
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

tmp<scalarField> Foam::fvDVM::getSoundSpeed()
{
    tmp<scalarField> tvf(sqrt(0.5 * Gamma() / lambdaVol()));
    return tvf;
}

void Foam::fvDVM::evolution()
{
    Reconstruction();
    CalcFluxSurf();
    Update();
}

// ************************************************************************* //
