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
            dimensionedScalar("dimkMean",UPrime2Mean.dimensions(),0)
        );

        volScalarField resLES
        (
            IOobject
            (
                "resLES",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimless",dimless,0)
        );

        volScalarField TKEMeanProd
        (
            IOobject
            (
                "TKEMeanProd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("dimTKEMeanProd",UPrime2Mean.dimensions()/dimTime,0)
        );
        
        volScalarField turbDiffusionMean
        (
            IOobject
            (
                "turbDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("turbDiffusionMean",UPrime2Mean.dimensions()/dimTime,0)
        );        
        
        volScalarField SGSDiffusionMean
        (
            IOobject
            (
                "SGSDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("SGSDiffusionMean",UPrime2Mean.dimensions()/dimTime,0)
        );
        
        volVectorField Gt
(
    IOobject
    (
        "Gt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Gt",dimTemperature/dimLength,vector (0,0,0))
);      


  volScalarField Ctt
        (
            IOobject
            (
                "Ctt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ctt",dimTemperature*dimTemperature/dimTime,0)
        );
     
         volScalarField Ptt
        (
            IOobject
            (
                "Ptt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ptt",dimTemperature*dimTemperature/dimTime,0)
        );
          volScalarField Dtt
        (
            IOobject
            (
                "Dtt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dtt",dimTemperature*dimTemperature/(dimLength*dimLength),0)
        );
           volScalarField Ttt
        (
            IOobject
            (
                "Ttt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Ttt",dimTemperature*dimTemperature/(dimTime),0)
        );
           volScalarField Mtt
        (
            IOobject
            (
                "Mtt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Mtt",dimTemperature*dimTemperature/(dimLength*dimLength),0)
        );
        volScalarField Cut
        (
            IOobject
            (
                "Cut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cut",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
         volScalarField Cvt
        (
            IOobject
            (
                "Cvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
         volScalarField Cwt
        (
            IOobject
            (
                "Cwt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Cwt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
        volScalarField Pvt
        (
            IOobject
            (
                "Pvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Pvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
          volScalarField Dvt
        (
            IOobject
            (
                "Dvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dvt",dimTemperature/(dimLength*dimTime),0)
        );
           volScalarField Dut
        (
            IOobject
            (
                "Dut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Dut",dimTemperature/(dimLength*dimTime),0)
        );
          volScalarField Tvt
        (
            IOobject
            (
                "Tvt",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("Tvt",dimTemperature*dimLength/(dimTime*dimTime),0)
        );
        volScalarField viscDiffusionMean
        (
            IOobject
            (
                "viscDiffusionMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("viscDiffusionMean",UPrime2Mean.dimensions()/dimTime,0)
        ); 
          Info<< "Reading field Tmean\n" << endl;
    volScalarField TMean
    (
        IOobject
        (
            "TMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
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
    volVectorField GradTUPMean
    (
        IOobject
        (
            "GradTUPMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
       volVectorField GradTP
    (
        IOobject
        (
            "GradTP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
  Info<< "Reading field TPrime2Mean\n" << endl;
    volScalarField TPrime2Mean
    (
        IOobject
        (
            "TPrime2Mean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
     Info<< "Reading field UUTMean\n" << endl;
    volTensorField UUTMean
    (
        IOobject
        (
            "UUTMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volTensorField GradUP
    (
        IOobject
        (
            "GradUP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
      
  volScalarField test
(
    IOobject
    (
        "test",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
           dimensionedScalar
      (
           "test", 
           dimTemperature/(dimLength*dimTime),
           0
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

        Info<< "    Calculating mean turbulent kinetic energy TKEMean" << endl;
        TKEMean = 0.5 * (
                            UPrime2Mean.component(0)
                          + UPrime2Mean.component(3)
                          + UPrime2Mean.component(5)
                        );

       // Info<< "    Calculating Mean LES resolutness resLES" << endl;
      //  forAll(resLES, i)
      //  {
      //      resLES[i] = TKEMean[i] / (TKEMean[i] + kMean[i]);
     //   }

        Info<< "    Calculating Mean TKE production TKEMeanProd" << endl;
        TKEMeanProd = - UPrime2Mean && fvc::grad(UMean);
        
        Info<< "    Calculating turbDiffusionMean" << endl;
        turbDiffusionMean = fvc::div(turbDiffMean);

        Info<< "    Calculating SGSDiffusionMean" << endl;
        SGSDiffusionMean = fvc::div(SGSDiffMean);        
        
        Info<< "    Calculating viscDiffusionMean" << endl;
        viscDiffusionMean = nuLam*fvc::laplacian(TKEMean);         
        
        // calculating temperature variance TT budgets//
        
        
        Info<< "    Calculating convection" << endl;
        Gt= fvc::grad(TMean);
        Ctt =  UMean & fvc::grad(TPrime2Mean);
        Info<< "    Calculating Production" << endl;
       Ptt = -2* (TUMean) & fvc::grad(TMean); 
        Info<< "Calculating Dissipation" << endl;
       Dtt = -2*(TTDissMean.component(0)+TTDissMean.component(1)+TTDissMean.component(2));//2.18591*10^-5 multiply it (value of alpha)
        Info<< "Calculating Turbulent Diffusion" << endl;
       Ttt = - fvc::div(TTUMean);
        Info<< "    Calculating Molecular DiffusionMean" << endl;
        Mtt = fvc::laplacian(TPrime2Mean);//2.18591*10^-5 multiply it (value of alpha)
        
        
        // calculating variance ut budgets//
        Info<< "    Calculating convection ut" << endl;
        Cut= UMean & fvc::grad(TUMean.component(0));
        Cvt= UMean & fvc::grad(TUMean.component(1));
        Cwt= UMean & fvc::grad(TUMean.component(2));
        Info<< "    Calculating Production ut" << endl;
        //Pvt= -((TUMean & fvc::grad(UMean.component(0)))+(UPrime2Mean.component(0)*Gt.component(0)+UPrime2Mean.component(1)*Gt.component(1)+UPrime2Mean.component(2)*Gt.component(2)));
        Pvt= -((TUMean & fvc::grad(UMean.component(1)))+(UPrime2Mean.component(1)*Gt.component(0)+UPrime2Mean.component(3)*Gt.component(1)+UPrime2Mean.component(4)*Gt.component(2)));
        //Pwt= -((TUMean & fvc::grad(UMean.component(2)))+(UPrime2Mean.component(2)*Gt.component(0)+UPrime2Mean.component(4)*Gt.component(1)+UPrime2Mean.component(5)*Gt.component(2)));
        Info<< "    Calculating Dissipation ut" << endl;
        Dut= -GradTUPMean.component(0); //(nu+alpha) multiplied
        Dvt= -GradTUPMean.component(1); //(nu+alpha) multiplied
        Info<< "    Calculating T diffusion ut" << endl;
    //    Tvt= -fvc::grad(UUTMean.component(3)+UUTMean.component(4)+UUTMean.component(5)); //(nu+alpha) multiplied
        //Tut= -fvc::grad(UUTMean); //(nu+alpha) multiplied
        Info<< "    Calculating Temperature pressure diffusion ut" << endl;
        //component zero for the first ut
        Info<< "    Calculating molecular diffusion ut" << endl;
       //test= GradTUMean.component(tensor::XX)*GradTMean.component(vector::X)+GradUPMean.component(tensor::XY)*GradTMean.component(vector::Y)+GradUPMean.component(tensor::XZ)*GradTMean.component(vector::Z);
      
        resLES.write();
        TKEMean.write();
        TKEMeanProd.write(); 
        turbDiffusionMean.write();        
        SGSDiffusionMean.write();
        viscDiffusionMean.write();  
        test.write();
        Gt.write();
       Ctt.write();  
       Ptt.write();
       Dtt.write(); 
       Ttt.write();
       Mtt.write();              
       Cut.write(); 
       Cvt.write();
       Cwt.write();
       Pvt.write(); 
       Dvt.write();    
       Dut.write();  
    }
   
    
    Info<< "End" << endl;
}


// ************************************************************************* //
    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
