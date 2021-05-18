/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2010-08-02 Eelco van Vliet: 1st public version of wallHeatFluxIncompressible:
  http://www.cfd-online.com/Forums/openfoam-solving/66705-wallheatflux-bc-not-constant-after-restart.html#post269812

2012-05-21 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post362191
  «modified the standard wallHeatflux utility which comes default with OF into
  a version for incompressible flows. Also removed a bug out of the code.»

2012-06-26 Eelco van Vliet:
  Quoting: http://www.cfd-online.com/Forums/openfoam-post-processing/101972-wallheatflux-utility-incompressible-case.html#post368330
  «p is now not required anymore.»

2014-06-22: Bruno Santos: Adapted to OpenFOAM 2.2.x.

2018-06-15: Bruno Santos @ FSD blueCAPE Lda: Adapted to OpenFOAM 5.x.

-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wallHeatFluxIncompressible

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.
    Based on wallHeatFlux with changes to allow it on incompressible flows
    Also removed a bug at the typeid checkline
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

#       include "createFields.H"
#       include "readTransportProperties.H"

         // update the turbulence fields
       turbulence->read();
}
      
{
       
    IOobject UPrime2MeanHeader
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // IOobject kMeanHeader
   // (
   //     "kMean",
   //     runTime.timeName(),
   //     mesh,
   //     IOobject::MUST_READ
   // );

    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
          
    IOobject turbDiffMeanHeader
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );         
    
    IOobject SGSDiffMeanHeader
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );     
   // if
   // (
  //      UPrime2MeanHeader.headerOk()
  //   && kMeanHeader.headerOk()
 //    && UMeanHeader.headerOk()
  //   && turbDiffMeanHeader.headerOk()
 //    && SGSDiffMeanHeader.headerOk()     
 //   )
    
    {
        Info<< "    Reading average field UPrime2Mean" << endl;
        const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);
         Info<< "    Reading average field nut" << endl;
     

       // Info<< "    Reading average field kMean" << endl;
      //  const volScalarField kMean(kMeanHeader, mesh);

        Info<< "    Reading average field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);                          
        
        Info<< "    Reading average field turbDiffMean" << endl;
        const volVectorField turbDiffMean(turbDiffMeanHeader, mesh);                                  
          
        Info<< "    Reading average field SGSDiffMean" << endl;
        const volVectorField SGSDiffMean(SGSDiffMeanHeader, mesh); 
        
        volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);                 


///////////////////////////////////////////////////////////////////////        
 volScalarField TKEMean
        (
            IOobject
            (
                "TKEMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("TKEMean",UPrime2Mean.dimensions(),0)
        );
       


volTensorField GU
(
    IOobject
    (
        "GU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedTensor
      (
           "GU", 
           dimVelocity/(dimLength),
           tensor::zero
      )
);
 volScalarField nutMean
    (
        IOobject
        (
            "nutMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


volScalarField produu
        (
            IOobject
            (
                "produu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("produu",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField prodvv
        (
            IOobject
            (
                "prodvv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("prodvv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField prodww
        (
            IOobject
            (
                "prodww",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("prodww",UPrime2Mean.dimensions()/dimTime,0)
        );

volScalarField produv
        (
            IOobject
            (
                "produv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("produv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField produw
        (
            IOobject
            (
                "produw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("produw",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField prodvw
        (
            IOobject
            (
                "prodvw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("prodvw",UPrime2Mean.dimensions()/dimTime,0)
        );
  #include "createPhi.H"  
        //~ #include "createFields.H" 
        
          singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu());                             

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );
     
        
        Info<< "    Calculating ProdMean" << endl;
	GU=fvc::grad(UMean);
        
    produu = -2*(UPrime2Mean.component(0)*GU.component(0)+UPrime2Mean.component(1)*GU.component(3)+UPrime2Mean.component(2)*GU.component(6)); 
prodvv = -2*(UPrime2Mean.component(1)*GU.component(1)+UPrime2Mean.component(3)*GU.component(4)+UPrime2Mean.component(4)*GU.component(7)); 
prodww = -2*(UPrime2Mean.component(2)*GU.component(2)+UPrime2Mean.component(4)*GU.component(5)+UPrime2Mean.component(5)*GU.component(8));

produv = -((UPrime2Mean.component(1)*GU.component(0)+UPrime2Mean.component(0)*GU.component(1)+UPrime2Mean.component(3)*GU.component(3))+(UPrime2Mean.component(1)*GU.component(4)+UPrime2Mean.component(4)*GU.component(6)+UPrime2Mean.component(2)*GU.component(7))); 

produw = -((UPrime2Mean.component(2)*GU.component(0)+UPrime2Mean.component(0)*GU.component(2)+UPrime2Mean.component(4)*GU.component(3))+(UPrime2Mean.component(1)*GU.component(5)+UPrime2Mean.component(5)*GU.component(6)+UPrime2Mean.component(2)*GU.component(8))); 

prodvw = -((UPrime2Mean.component(2)*GU.component(1)+UPrime2Mean.component(1)*GU.component(2)+UPrime2Mean.component(4)*GU.component(4))+(UPrime2Mean.component(3)*GU.component(5)+UPrime2Mean.component(5)*GU.component(7)+UPrime2Mean.component(4)*GU.component(8))); 
                
        
       
    
      produu.write();
	prodvv.write();
	prodww.write();
 produv.write();
	produw.write();
	prodvw.write();
      
   
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
  
    return 0;
}

// ************************************************************************* //
