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
    const surfaceVectorField nf = mesh_.Sf() / mesh_.magSf();
    // Internal faces
    forAll(nf, facei)
    {
        const vector a = nf[facei];
        vector t;
        if ((abs(a.x()) < abs(a.y())) && (abs(a.x()) < abs(a.z())))
        {
            t = vector(1, 0, 0);
            t = t - (t & a) * a;
        }
        else if (abs(a.y()) < abs(a.z()))
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
        LocalFrameSurf_[facei] = tensor(a, b, c);
    }
    // Boundary faces
    forAll(LocalFrameSurf_.boundaryField(), patchi)
    {
        fvsPatchField<tensor> &LocalFrameSurfPatch =
            LocalFrameSurf_.boundaryFieldRef()[patchi];
        const fvsPatchField<vector> &nfPatch =
            nf.boundaryField()[patchi];

        // const labelUList &pOwner = mesh_.boundary()[patchi].faceCells();
        forAll(LocalFrameSurfPatch, pFacei)
        {
            const vector a = nfPatch[pFacei];
            vector t;
            if ((abs(a.x()) < abs(a.y())) && (abs(a.x()) < abs(a.z())))
            {
                t = vector(1, 0, 0);
                t = t - (t & a) * a;
            }
            else if (abs(a.y()) < abs(a.z()))
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
            LocalFrameSurfPatch[pFacei] = tensor(a, b, c);
        }
    }
}

void Foam::fvDVM::Reconstruction()
{
    forAll(DV_, DVid)
        DV_[DVid]
            .Reconstruction();
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

    // Info << "internal" << endl;
    // Internal faces first
    forAll(owner, facei)
    {
        // Info << facei << endl;
        // Local geometric information
        label own = owner[facei];
        label nei = neighbour[facei];
        tensor frame = localFrame[facei];

        // Face variable
        scalarField h(nXi()), b(nXi());
        vectorField dh(nXi()), db(nXi()), ui(nXi()); // Micro gradient and velocity in local frame
        scalar rho, lambda, rhoE;
        vector U, rhoU;

        // Gradient and micro slope
        scalar Nrho, T1rho, T2rho;    // Normal derivate of conserved variables
        vector NrhoU, T1rhoU, T2rhoU; // First tangential derivate of conserved variables
        scalar NrhoE, T1rhoE, T2rhoE; // Second tangetial derivate of conserved variables
        scalar aRho, bRho, cRho, aTrho;
        vector aU, bU, cU, aTU;
        scalar aLambda, bLambda, cLambda, aTlambda;

        // Moment
        scalarField Mu(7), MuL(7), MuR(7), Mv1(6), Mv2(6), Mxi(3);
        scalarField Mau0(5), Mau(5), Mbv1(5), Mcv2(5), MauT(5), Mt(5);

        // Initilization
        rho = 0.0;
        U = rhoU = vector(0, 0, 0);
        lambda = rhoE = 0.0;
        Nrho = T1rho = T2rho = 0.0;
        NrhoU = T1rhoU = T2rhoU = vector(0, 0, 0);
        NrhoE = T1rhoE = T2rhoE = 0.0;

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
            // Calculate conVars at cell interface in local frame
            const scalar &weight = dv.weight();
            const scalar &hi = h[dvi];
            const scalar &bi = b[dvi];
            const vector &dhi = dh[dvi]; // Micro gradient dh already in local frame
            const vector &dbi = db[dvi]; // Micro gradient db already in local frame

            // Obtain the conserved variables in local frame by collision
            rho += weight * hi;
            rhoU += weight * xii * hi; // local frame
            rhoE += 0.5 * weight * (magSqr(xii) * hi + bi);

            // Obtain the gradient of conserved variables by collision
            Nrho += weight * dhi.x();
            NrhoU += weight * xii * dhi.x();
            NrhoE += 0.5 * weight * (magSqr(xii) * dhi.x() + dbi.x());

            T1rho += weight * dhi.y();
            T1rhoU += weight * xii * dhi.y();
            T1rhoE += 0.5 * weight * (magSqr(xii) * dhi.y() + dbi.y());

            T2rho += weight * dhi.z();
            T2rhoU += weight * xii * dhi.z();
            T2rhoE += 0.5 * weight * (magSqr(xii) * dhi.z() + dbi.z());
        }
        ConservedToPrim(rho, rhoU, rhoE, rho, U, lambda);
        // Info << rho << tab << rhoU << tab << rhoE << endl;
        // Info << rho << tab << U << tab << lambda << endl;
        // Usurf_[facei] = frame.T() & U; // Used for coNum evaluation
        MicroSlope(Nrho, NrhoU, NrhoE, rho, U, lambda, aRho, aU, aLambda);
        MicroSlope(T1rho, T1rhoU, T1rhoE, rho, U, lambda, bRho, bU, bLambda);
        MicroSlope(T2rho, T2rhoU, T2rhoE, rho, U, lambda, cRho, cU, cLambda);

        CalcMoment(rho, U, lambda, Mu, Mv1, Mv2, Mxi, MuL, MuR);
        Mau = Moment_auv1v2xi(aRho, aU, aLambda, Mu, Mv1, Mv2, Mxi, 1, 0, 0);
        Mbv1 = Moment_auv1v2xi(bRho, bU, bLambda, Mu, Mv1, Mv2, Mxi, 0, 1, 0);
        Mcv2 = Moment_auv1v2xi(cRho, cU, cLambda, Mu, Mv1, Mv2, Mxi, 0, 0, 1);
        scalar swrho = -rho * (Mau[0] + Mbv1[0] + Mcv2[0]);
        vector swrhoU = vector(-rho * (Mau[1] + Mbv1[1] + Mcv2[1]), -rho * (Mau[2] + Mbv1[2] + Mcv2[2]), -rho * (Mau[3] + Mbv1[3] + Mcv2[3]));
        scalar swrhoE = -rho * (Mau[4] + Mbv1[4] + Mcv2[4]);
        MicroSlope(swrho, swrhoU, swrhoE, rho, U, lambda, aTrho, aTU, aTlambda);

        // Calculate collision time and some time integration terms
        scalar tau = GetTau(rho, lambda);

        Mt[3] = tau * (1.0 - exp(-dt / tau));
        Mt[4] = -tau * dt * exp(-dt / tau) + tau * Mt[3];
        Mt[0] = dt - Mt[3];
        Mt[1] = -tau * Mt[0] + Mt[4];
        Mt[2] = 0.5 * dt * dt - tau * Mt[0];

        Mau0 = Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, 1, 0, 0, 0);
        Mau = Moment_auv1v2xi(aRho, aU, aLambda, Mu, Mv1, Mv2, Mxi, 2, 0, 0);
        Mbv1 = Moment_auv1v2xi(bRho, bU, bLambda, Mu, Mv1, Mv2, Mxi, 1, 1, 0);
        Mcv2 = Moment_auv1v2xi(cRho, cU, cLambda, Mu, Mv1, Mv2, Mxi, 1, 0, 1);
        MauT = Moment_auv1v2xi(aTrho, aTU, aTlambda, Mu, Mv1, Mv2, Mxi, 1, 0, 0);

        rhoFluxSurf_[facei] = rho * (Mt[0] * Mau0[0] + Mt[1] * (Mau[0] + Mbv1[0] + Mcv2[0]) + Mt[2] * MauT[0]);
        rhoUfluxSurf_[facei].x() = rho * (Mt[0] * Mau0[1] + Mt[1] * (Mau[1] + Mbv1[1] + Mcv2[1]) + Mt[2] * MauT[1]);
        rhoUfluxSurf_[facei].y() = rho * (Mt[0] * Mau0[2] + Mt[1] * (Mau[2] + Mbv1[2] + Mcv2[2]) + Mt[2] * MauT[2]);
        rhoUfluxSurf_[facei].z() = rho * (Mt[0] * Mau0[3] + Mt[1] * (Mau[3] + Mbv1[3] + Mcv2[3]) + Mt[2] * MauT[3]);
        rhoEfluxSurf_[facei] = rho * (Mt[0] * Mau0[4] + Mt[1] * (Mau[4] + Mbv1[4] + Mcv2[4]) + Mt[2] * MauT[4]);

        vector qf = vector(0, 0, 0);
        forAll(DV_, dvi)
        {
            DiscreteVelocityPoint &dv = DV_[dvi];
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

    // Info << "boundary" << endl;
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
                Usurf_.boundaryFieldRef()[patchi][pFacei] = Upatch[pFacei];
                // Local geometric information
                const label own = pOwner[pFacei];
                tensor frame = FramePatch[pFacei];

                // Face variable
                scalarField h(nXi()), b(nXi());
                vectorField ui(nXi());
                scalar rho, incidence = 0.0, reflection = 0.0;
                vector U = frame & Upatch[pFacei]; // Macro velocity in local frame
                scalar lambda = lambdaPatch[pFacei];

                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const vector xii = frame & dv.xi(); // Micro velocity in local frame
                    ui[dvi] = xii;
                    const scalar &weight = dv.weight();
                    const scalar &vn = xii.x();

                    h[dvi] = dv.hVol()[own] + (dv.hGradVol()[own] & (CfPatch[pFacei] - C[own]));
                    b[dvi] = dv.bVol()[own] + (dv.bGradVol()[own] & (CfPatch[pFacei] - C[own]));

                    if (vn >= VSMALL) // Comming from own
                    {
                        incidence += weight * vn * h[dvi];
                    }
                    else // Comming form nei
                    {
                        reflection += pow(sqrt(lambda / pi), D) * weight * vn * exp(-lambda * magSqr(U - xii));
                    }
                }
                rho = -incidence / reflection;

                forAll(DV_, dvi)
                {
                    DiscreteVelocityPoint &dv = DV_[dvi];
                    const scalar &weight = dv.weight();
                    const vector &xii = ui[dvi];
                    const scalar &vn = xii.x();
                    scalar hi;
                    scalar bi;

                    if (vn >= VSMALL) // Comming from own
                    {
                        hi = h[dvi];
                        bi = b[dvi];
                    }
                    else // Comming form nei
                    {
                        DiscreteMaxwell(hi, bi, rho, U, lambda, xii);
                    }
                    rhoFluxPatch[pFacei] += weight * vn * hi;
                    rhoUfluxPatch[pFacei] += weight * vn * xii * hi;
                    rhoEfluxPatch[pFacei] += 0.5 * weight * vn * (magSqr(xii) * hi + bi);

                    dv.UpdatehFlux().boundaryFieldRef()[patchi][pFacei] = dt * vn * hi * magSfPatch[pFacei];
                    dv.UpdatebFlux().boundaryFieldRef()[patchi][pFacei] = dt * vn * bi * magSfPatch[pFacei];
                }
                rhoFluxPatch[pFacei] = dt * rhoFluxPatch[pFacei] * magSfPatch[pFacei];
                rhoUfluxPatch[pFacei] = dt * (frame.T() & rhoUfluxPatch[pFacei]) * magSfPatch[pFacei];
                rhoEfluxPatch[pFacei] = dt * rhoEfluxPatch[pFacei] * magSfPatch[pFacei];
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
    const label &D = mesh_.nSolutionD();
    rho = rhoP;
    rhoU = rho * U;
    rhoE = 0.5 * rho * (magSqr(U) + (K + D) / (2.0 * lambda));
}

void Foam::fvDVM::ConservedToPrim(scalar &rhoC, vector &rhoU, scalar &rhoE, scalar &rho, vector &U, scalar &lambda)
{
    const label &K = KInner_;
    const label &D = mesh_.nSolutionD();
    rho = rhoC;
    U = rhoU / rho;
    lambda = rho * (K + D) / (4.0 * (rhoE - 0.5 * (rhoU & U)));
}

void Foam::fvDVM::MicroSlope(scalar &srho, vector &srhoU, scalar &srhoE, scalar &rho, vector &U, scalar &lambda, scalar &grho, vector &gU, scalar &glambda)
{
    const label &K = KInner_;
    const label &D = mesh_.nSolutionD();
    const vector dU = (srhoU - U * srho) / rho;
    const scalar dE = (srhoE - 0.5 * (magSqr(U) + (K + D) / (2.0 * lambda)) * srho) / rho;

    glambda = 4.0 * lambda * lambda / (K + D) * (2.0 * dE - 2.0 * (U & dU));
    gU = 2.0 * lambda * dU - U * glambda;
    grho = srho / rho - (U & gU) - 0.5 * (magSqr(U) + (K + D) / (2.0 * lambda)) * glambda;
}

void Foam::fvDVM::CalcMoment(const scalar &rho, const vector &U, const scalar &lambda,
                             scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                             scalarField &MuL, scalarField &MuR)
{
    const label &K = KInner_;
    const label &D = mesh_.nSolutionD();
    label innerK = K + 3 - D;

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
scalarField Foam::fvDVM::Moment_auv1v2xi(scalar &srho, vector &sU, scalar &sLambda,
                                         scalarField &Mu, scalarField &Mv1, scalarField &Mv2, scalarField &Mxi,
                                         label alpha, label beta, label gamma)
{
    scalarField auv1v2xi(5);
    auv1v2xi = srho * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 0, 0) +
               sU.x() * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 1, beta + 0, gamma + 0, 0) +
               sU.y() * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 1, gamma + 0, 0) +
               sU.z() * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 1, 0) +
               0.5 * sLambda * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 2, beta + 0, gamma + 0, 0) +
               0.5 * sLambda * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 2, gamma + 0, 0) +
               0.5 * sLambda * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 2, 0) +
               0.5 * sLambda * Moment_uv1v2xi(Mu, Mv1, Mv2, Mxi, alpha + 0, beta + 0, gamma + 0, 2);
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
      qVol_(
          IOobject(
              "q",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE),
          mesh_,
          dimensionedVector("0", dimless, vector(0, 0, 0))),
      Usurf_(
          IOobject(
              "Usurf",
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh_,
          dimensionedVector("0", U.dimensions(), vector(0, 0, 0))),
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
    Usurf_ = fvc::interpolate(Uvol_, "linear"); // for first time Dt calculation.
    SetLocalFrame();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::getCoNum(scalar &maxCoNum, scalar &meanCoNum)
{
    scalar dt = time_.deltaTValue();
    Usurf_ = fvc::interpolate(Uvol_, "linear");
    scalarField UbyDx =
        mesh_.deltaCoeffs() * (mag(Usurf_) + sqrt(scalar(mesh_.nSolutionD())) * xiMax_);
    maxCoNum = gMax(UbyDx) * dt;
    meanCoNum = gSum(UbyDx) / UbyDx.size() * dt;
}

void Foam::fvDVM::evolution()
{
    // Info << "Begin evolution" << endl;
    Reconstruction();
    // Info << "Done reconstruction of nonequilibrium distribution" << endl;
    CalcFluxSurf();
    // Info << "Done Calculation of macro and micro flux" << endl;
    Update();
    // Info << "Done update macro and micro variables" << endl;
}

// ************************************************************************* //
