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
     //   Info<< "    Reading average field nut" << endl;
      //  const volScalarField nutMean(nutMean, mesh);


   
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
        
       
       
      
  

volScalarField UUUMean
(
    IOobject
    (
        "UUUMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
volScalarField UUVMean
(
    IOobject
    (
        "UUVMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField UUWMean
(
    IOobject
    (
        "UUWMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
         
);
volScalarField VVUMean
(
    IOobject
    (
        "VVUMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
         
);
volScalarField VVVMean
(
    IOobject
    (
        "VVVMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
);
volScalarField VVWMean
(
    IOobject
    (
        "VVWMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
volScalarField WWUMean
(
    IOobject
    (
        "WWUMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
          
);
volScalarField WWVMean
(
    IOobject
    (
        "WWVMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
           
);
volScalarField WWWMean
(
    IOobject
    (
        "WWWMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    
           
);
 Info<< "Reading sdfield TUmean\n" << endl;
volScalarField UVWMean
(
    IOobject
    (
        "UVWMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
    
           
);
     Info<< "Reading field TUmean\n" << endl;
    volVectorField TUMean
    (
        IOobject
        (
            "TUMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
      Info<< "Reading field TTUmean\n" << endl;
    volVectorField TTUMean
    (
        IOobject
        (
            "TTUMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
  Info<< "Reading field TTDissMean\n" << endl;
    volVectorField TTDissMean
    (
        IOobject
        (
            "TTDissMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
  
  
     

      


volScalarField Tdiffuu
        (
            IOobject
            (
                "Tdiffuu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffuu",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField Tdiffvv
        (
            IOobject
            (
                "Tdiffvv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffvv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField Tdiffww
        (
            IOobject
            (
                "Tdiffww",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffww",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField Tdiffuv
        (
            IOobject
            (
                "Tdiffuv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffuv",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField Tdiffuw
        (
            IOobject
            (
                "Tdiffuw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffuw",UPrime2Mean.dimensions()/dimTime,0)
        );
volScalarField Tdiffvw
        (
            IOobject
            (
                "Tdiffvw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdiffvw",UPrime2Mean.dimensions()/dimTime,0)
        );

volScalarField Tdifftke
        (
            IOobject
            (
                "Tdifftke",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tdifftke",UPrime2Mean.dimensions()/dimTime,0)
        );

volVectorField Tdiffuuug
        (
            IOobject
            (
                "Tdiffuuug",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffuuug",UPrime2Mean.dimensions()/dimTime,vector (0,0,0) )
      
); 
volVectorField Tdiffuuvg
        (
            IOobject
            (
                "Tdiffuuvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffuuvg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffuuwg
        (
            IOobject
            (
                "Tdiffuuwg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffuuwg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffvvug
        (
            IOobject
            (
                "Tdiffvvug",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffvvug", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffvvvg
        (
            IOobject
            (
                "Tdiffvvvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffvvvg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffvvwg
        (
            IOobject
            (
                "Tdiffvvwg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffvvwg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 

volVectorField Tdiffwwug
        (
            IOobject
            (
                "Tdiffwwug",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffwwug", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffwwvg
        (
            IOobject
            (
                "Tdiffwwvg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffwwvg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdiffwwwg
        (
            IOobject
            (
                "Tdiffwwwg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffwwwg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 

volVectorField Tdiffuvwg
        (
            IOobject
            (
                "Tdiffuvwg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdiffuvwg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdifftkexg
        (
            IOobject
            (
                "Tdifftkexg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdifftkexg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdifftkeyg
        (
            IOobject
            (
                "Tdifftkeyg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdifftkeyg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 
volVectorField Tdifftkezg
        (
            IOobject
            (
                "Tdifftkezg",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
dimensionedVector
      (
           "Tdifftkezg", 
          UPrime2Mean.dimensions()/dimTime,
           vector (0,0,0)
      )
); 


  #include "createPhi.H"  
        //~ #include "createFields.H" 
        
          singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::turbulenceModel> turbulence
        (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );
        
        const volScalarField nuLam(turbulence->nu()); 
       
                
    const volSymmTensorField G(turbulence->R());
    
    const volScalarField nu_SGS(turbulence->nut());                            

     

     
        
        Info<< "    Calculating viscDiffusionMean" << endl;
         //add nut to all
        Tdiffuuug = fvc::grad(UUUMean);
        Tdiffuuvg = fvc::grad(UUVMean);
 	Tdiffuuwg = fvc::grad(UUWMean);
	Tdiffvvug = fvc::grad(VVUMean);
	Tdiffvvvg = fvc::grad(VVVMean);
	Tdiffvvwg = fvc::grad(VVWMean);
	Tdiffwwug = fvc::grad(WWUMean);
	Tdiffwwvg = fvc::grad(WWVMean);
	Tdiffwwwg = fvc::grad(WWWMean);
	Tdiffuvwg = fvc::grad(UVWMean);
        Tdifftkexg = fvc::grad((UUUMean+VVUMean+WWUMean)*0.5);
        Tdifftkeyg = fvc::grad((UUVMean+VVVMean+WWVMean)*0.5);
	Tdifftkezg = fvc::grad((UUWMean+VVWMean+WWWMean)*0.5);
        Tdiffuu = -(Tdiffuuug.component(vector::X)+Tdiffuuvg.component(vector::Y)+Tdiffuuwg.component(vector::Z));
        Tdiffvv = -(Tdiffvvug.component(vector::X)+Tdiffvvvg.component(vector::Y)+Tdiffvvwg.component(vector::Z));
 	Tdiffww = -(Tdiffwwug.component(vector::X)+Tdiffwwvg.component(vector::Y)+Tdiffwwwg.component(vector::Z));
	Tdiffuv = -(Tdiffuuvg.component(vector::X)+Tdiffvvug.component(vector::Y)+Tdiffuvwg.component(vector::Z));
	Tdiffuw = -(Tdiffuuwg.component(vector::X)+Tdiffuvwg.component(vector::Y)+Tdiffwwug.component(vector::Z));
	Tdiffvw = -(Tdiffuvwg.component(vector::X)+Tdiffvvwg.component(vector::Y)+Tdiffwwvg.component(vector::Z));
	Tdifftke = -(Tdifftkexg.component(vector::X)+Tdifftkeyg.component(vector::Y)+Tdifftkezg.component(vector::Z));
	
      
        
      
        Tdiffuu.write();
        Tdiffvv.write();
 	Tdiffww.write();
	Tdiffuv.write();
	Tdiffuw.write();
	Tdiffvw.write();
	Tdifftke.write();
   
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
