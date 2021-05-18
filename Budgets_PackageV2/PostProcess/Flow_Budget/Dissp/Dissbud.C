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



    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
          
    

    {
        Info<< "    Reading average field UPrime2Mean" << endl;
        const volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);
         Info<< "    Reading average field nut" << endl;
     

       // Info<< "    Reading average field kMean" << endl;
      //  const volScalarField kMean(kMeanHeader, mesh);

        Info<< "    Reading average field UMean" << endl;
        const volVectorField UMean(UMeanHeader, mesh);                          
        
       
        
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

volScalarField dissuu
        (
            IOobject
            (
                "dissuu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissuu",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissvv
        (
            IOobject
            (
                "dissvv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissvv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissww
        (
            IOobject
            (
                "dissww",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissww",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissuv
        (
            IOobject
            (
                "dissuv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissuv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissuw
        (
            IOobject
            (
                "dissuw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissuw",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissvw
        (
            IOobject
            (
                "dissvw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dissvw",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField dissUVMean
(
    IOobject
    (
        "dissUVMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
         
);
volScalarField dissUWMean
(
    IOobject
    (
        "dissUWMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
volScalarField dissVWMean
(
    IOobject
    (
        "dissVWMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
volTensorField tkeDissMean
(
    IOobject
    (
        "tkeDissMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
  #include "createPhi.H"  
        //~ #include "createFields.H" 
        
          singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu());                             

      
       
          
        
        Info<< "    Calculating DissipationMean" << endl;
        
        dissuu = 2*(nuLam)*(tkeDissMean.component(0)+tkeDissMean.component(3)+tkeDissMean.component(6))/(nuLam); 
        dissvv = 2*(nuLam)*(tkeDissMean.component(1)+tkeDissMean.component(4)+tkeDissMean.component(7))/(nuLam);
        dissww = 2*(nuLam)*(tkeDissMean.component(2)+tkeDissMean.component(5)+tkeDissMean.component(8))/(nuLam);
        dissuv = 2*(nuLam)*(dissUVMean)/(nuLam); 
        dissuw = 2*(nuLam)*(dissUWMean)/(nuLam); 
        dissvw = 2*(nuLam)*(dissVWMean)/(nuLam);         
        
    
      dissuu.write();
      dissvv.write();
      dissww.write();
      dissuv.write();
      dissuw.write();
      dissvw.write();
   
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
  
    return 0;
}

// ************************************************************************* //
