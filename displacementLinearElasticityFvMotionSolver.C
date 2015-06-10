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

#include "displacementLinearElasticityFvMotionSolver.H"
#include "motionDiffusivity.H"
//#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"

// additional headers for elasticity solve
#include "fvcReconstruct.H"
#include "fvCFD.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLinearElasticityFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLinearElasticityFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLinearElasticityFvMotionSolver::displacementLinearElasticityFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolverCore(mesh),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellDisplacement",
            pointDisplacement_.dimensions(),
            vector::zero
        ),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(NULL),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID(coeffDict().lookup("frozenPointsZone"))
      : -1
    )//,
    ///////////////////////////////////////////////////////////////////////////
    //
    // from solidEquilibriumDisplacementFoam
    //
//    mu_(
//        IOobject
//        (
//            "mu",
//            mesh.time().timeName(),
//            mesh
//        ),
//        diffusivityPtr_.operator()() // <----- THIS DOESNT WORK
//    ),
//    Dcorr_
//    (
////        IOobject
////        (
////            "Dcorr",
////            mesh.time().timeName(),
////            mesh
////        ),
////        cellDisplacement_
//        IOobject
//        (
//            "Dcorr",
//            mesh.time().timeName(),
//            mesh
//        ),
//        fvMesh_,
//        dimensionedVector
//        (
//            "Dcorr",
//            pointDisplacement_.dimensions(),
//            vector::zero
//        ),
//        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
//    )//,
//    sigmaD_
//    (
//        IOobject
//        (
//            "sigmaD",
//            mesh.time().timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        //mu*twoSymm(fvc::grad(D)) + (lambda*I)*tr(fvc::grad(D))
//        mu_*twoSymm(fvc::grad(cellDisplacement_))
//    ),
//    sigmaExp_
//    (
//        IOobject
//        (
//            "sigmaExp",
//            mesh.time().timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//      //  (lambda - mu)*fvc::grad(Dcorr) + mu*fvc::grad(Dcorr)().T()
//      //+ (lambda*I)*tr(fvc::grad(Dcorr))
//        mu_*( -fvc::grad(Dcorr_) + fvc::grad(Dcorr_)().T() )
//    )
    ///////////////////////////////////////////////////////////////////////////
{

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

//    if (debug)
//    {
        Info<< "displacementLinearElasticityFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
//    }


    if (io.headerOk())
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "displacementLinearElasticityFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }

    Info<< "Constructed displacementLinearElasticityFvMotionSolver (in fvMotionSolver/fvMotionSolvers)\n" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLinearElasticityFvMotionSolver::
~displacementLinearElasticityFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::displacementLinearElasticityFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_.valid())
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }
    return diffusivityPtr_();
}


Foam::tmp<Foam::pointField>
Foam::displacementLinearElasticityFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "displacementLinearElasticityFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().internalField() =
            points0()
          + pointDisplacement_.internalField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().internalField());

        return tmp<pointField>(pointLocation_().internalField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            points0() + pointDisplacement_.internalField()
        );

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                tcurPoints()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(tcurPoints());

        return tcurPoints;
    }
}


void Foam::displacementLinearElasticityFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivity().correct();
    pointDisplacement_.boundaryField().updateCoeffs();

////////////////////////////////////////////////////////////////////////////////
//
// Laplacian solver:
// solves div( Gamma * grad(U) ) == 0
//
// - Definition of laplacian is in $FOAM_SRC/finiteVolume/finiteVolume/fvm/fvmLaplacian.C
// - diffusivity() returns pointer to diffusivity field, Foam::motionDiffusivity
// - diffusivity().operator() returns a temporary surfaceScalarField, a type of GeometricField

//    Foam::solve
//    (
//        fvm::laplacian
//        (
//            diffusivity().operator()(),
//            cellDisplacement_,
//            "laplacian(diffusivity,cellDisplacement)"
//        )
//    );

////////////////////////////////////////////////////////////////////////////////
//
// from solidEquilibriumDisplacementFoam
//
//    {
//        Info<< "Reconstructing mu" << endl;
//        volScalarField mu
//        //volVectorField mu // reconstruct returns a volVectorField(?)
//        (
//            IOobject
//            (
//                "mu",
//                fvMesh_.time().timeName(),
//                fvMesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            //fvc::reconstruct( diffusivity().operator()() )
//            mag(fvc::reconstruct( diffusivity().operator()() )) // TODO: DOES THIS WORK???
//        );
//
////        Info<< "Reading/setting displacement correction field Dcorr" << endl;
////        volVectorField Dcorr
////        (
////            IOobject
////            (
////                "Dcorr",
////                fvMesh_.time().timeName(),
////                fvMesh_
////            ),
////            cellDisplacement_ // this is a volVectorField
////        );
////        Info<< Dcorr << endl;
////        Dcorr *= 0.0;
//
//        Info<< "Calculating stress field sigmaD" << endl;
//        volSymmTensorField sigmaD
//        (
//            IOobject
//            (
//                "sigmaD",
//                fvMesh_.time().timeName(),
//                fvMesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            //mu*twoSymm(fvc::grad(D)) + (lambda*I)*tr(fvc::grad(D))
//            mu*twoSymm(fvc::grad(cellDisplacement_))
//        );
//
//        Info<< "Calculating stress field sigmaExp" << endl;
//        volTensorField sigmaExp
//        (
//            IOobject
//            (
//                "sigmaExp",
//                fvMesh_.time().timeName(),
//                fvMesh_,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//        //  (lambda - mu)*fvc::grad(Dcorr) + mu*fvc::grad(Dcorr)().T()
//        //+ (lambda*I)*tr(fvc::grad(Dcorr))
//            mu*( -fvc::grad(Dcorr_) + fvc::grad(Dcorr_)().T() )
//        );
//
//        // DEBUG
//        //Info<< "cellDisplacement_\n" << cellDisplacement_ << endl;
//        //Info<< "Dcorr_\n" << Dcorr_ << endl;
//
//        volVectorField test
//        (
//            IOobject
//            (
//                "test",
//                fvMesh_.time().timeName(),
//                fvMesh_
//            ),
//            cellDisplacement_
//        );
//        test *= 0.0;
//
//
//        //-- main linear elasticity solver loop
//        Info<< "Begin linear elasticity solve" << endl;
//        solverPerformance solverPerf;
//        scalar nIter=0;
//        {
//            ++nIter;
//
//            Info<< "--executing solver" << endl;
//            solverPerf = Foam::solve
//            (
//                //fvm::laplacian(2*mu + lambda, Dcorr, "laplacian(DD,Dcorr)")
////                fvm::laplacian
////                (
////                    2*mu, // =Gamma, a volScalarField         (typedef GeometricField< scalar, fvPatchField, volMesh >)
////                    Dcorr // needs to be a pointVectorField?? (typedef GeometricField< vector, pointPatchField, pointMesh >)
////                          //      instead of a volVectorField (typedef GeometricField< vector, fvPatchField, volMesh >)
////                )
////                fvm::laplacian(Dcorr)
////              + fvc::div(sigmaExp + sigmaD)
//                
////                fvm::laplacian
////                (
////                    2*mu,
////                    //cellDisplacement_,
////                    Dcorr_,
////                    "laplacian(two_mu,Dcorr)"
////                )
//
//                fvm::laplacian(test)
//            );
//
//            Info<< "--updating cell displacement" << endl;
//            //cellDisplacement_ += accFac*Dcorr;
//            cellDisplacement_ += Dcorr_;
//
//            Info<< "--updating stress fields" << endl;
//            {
//                volTensorField gradDcorr(fvc::grad(Dcorr_));
//
//                //sigmaExp =
//                //    (lambda - mu)*gradDcorr + mu*gradDcorr.T()
//                //  + (lambda*I)*tr(gradDcorr);
//                sigmaExp =
//                    mu*(-gradDcorr + gradDcorr.T());
//
//                //sigmaD += accFac*(mu*twoSymm(gradDcorr) + (lambda*I)*tr(gradDcorr));
//                sigmaD += mu*twoSymm(gradDcorr);
//            }
//
//            //#include "calculateStress.H" // for output
//            //#include "kineticEnergyLimiter.H"
//        }
//        Info<< "Linear elasticity solver performed " << nIter << " iterations\n" << endl;
//        
//    }//end of linear elasticity analogy solver

////////////////////////////////////////////////////////////////////////////////
//
// same approach as displacementSBRStressFvMotionSolver
//

    surfaceScalarField mu(diffusivityPtr_->operator()());

    volTensorField DU(fvc::grad(cellDisplacement_));

    Foam::solve
    (
        fvm::laplacian
        (
            2*mu,
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )

      + fvc::div
        (
            mu
           *(
                (
                    cellDisplacement_.mesh().Sf()
                  & fvc::interpolate(DU.T() - DU)
                )

                // Solid-body rotation "lambda" term
//              - cellDisplacement_.mesh().Sf()*fvc::interpolate(tr(gradCd))
            )
        )
    );

}


void Foam::displacementLinearElasticityFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.clear();
}


// ************************************************************************* //
