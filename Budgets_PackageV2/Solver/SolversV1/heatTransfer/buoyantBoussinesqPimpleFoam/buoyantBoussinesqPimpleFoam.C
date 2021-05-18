/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    buoyantBoussinesqPimpleFoam

Description
    Transient solver for buoyant, turbulent flow of incompressible fluids.

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        \frac{beta(T - T_{ref})}{rho_{ref}} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "pimpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"
    #include "products.H"
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
         
           const volScalarField nuLam(turbulence->nu());        
                
    const volSymmTensorField G(turbulence->R());
    
    const volScalarField nu_SGS(turbulence->nut());
    
   // ..........................................Flow budgets....................................................// 
            UPrime = (U-UMean);         //Resolved Velocity Fluctuation Vector
            pPrime = (p-pMean); 
          //  symmTensor Ei(0,0,0,0,0,0);
          //  D=(symm(fvc::grad(U)));//Sij
            //Tkk= cmptMultiply(D,Ei);//trace //will be zero for incompressible flows??
          // Tij= 2.0*nu_SGS*(D); //not nu_sgs but nu?
          //tensor E(1,0,0,0,1,0,0,0,1);
         //  SGSdiffuT= (UPrime.component(vector::X)*Tij);
         //  SGSdiffvT= (UPrime.component(vector::Y)*Tij);
          // SGSdiffwT= (UPrime.component(vector::Z)*Tij);   //Components for SGS diffusion due to modelling
           GradUP= (fvc::grad(UPrime));
           //SGSdissuxT= (GradUP.component(tensor::XX)*Tij);
           //SGSdissuyT= (GradUP.component(tensor::YX)*Tij);
           //SGSdissuzT= (GradUP.component(tensor::ZX)*Tij);
           //SGSdissvxT= (GradUP.component(tensor::XY)*Tij);
           //SGSdissvyT= (GradUP.component(tensor::YY)*Tij);
           //SGSdissvzT= (GradUP.component(tensor::ZY)*Tij);
           //SGSdisswxT= (GradUP.component(tensor::XZ)*Tij);
           //SGSdisswyT= (GradUP.component(tensor::YZ)*Tij);
           //SGSdisswzT= (GradUP.component(tensor::ZZ)*Tij);  //Components for SGS dissipation due to modelling
           
           DP=(symm(fvc::grad(UPrime)));//Sij
            //Tkk= cmptMultiply(D,Ei);//trace //will be zero for incompressible flows??
           TijP= 2.0*nu_SGS*(DP); //not nu_sgs but nu?
           SGSdiffuTP= (UPrime.component(vector::X)*TijP);
           SGSdiffvTP= (UPrime.component(vector::Y)*TijP);
           SGSdiffwTP= (UPrime.component(vector::Z)*TijP);   //Components for SGS diffusion due to modelling
         //  GradUP= (fvc::grad(UPrime));
           SGSdissuxTP= (GradUP.component(tensor::XX)*TijP);
           SGSdissuyTP= (GradUP.component(tensor::YX)*TijP);
           SGSdissuzTP= (GradUP.component(tensor::ZX)*TijP);
           SGSdissvxTP= (GradUP.component(tensor::XY)*TijP);
           SGSdissvyTP= (GradUP.component(tensor::YY)*TijP);
           SGSdissvzTP= (GradUP.component(tensor::ZY)*TijP);
           SGSdisswxTP= (GradUP.component(tensor::XZ)*TijP);
           SGSdisswyTP= (GradUP.component(tensor::YZ)*TijP);
           SGSdisswzTP= (GradUP.component(tensor::ZZ)*TijP);  //Components for SGS dissipation due to modelling
           
          
          
           // Ti= (D.DiagTensor());           //Resolved Pressure Fluctuation
            B = -2.0*nu_SGS*(symm(fvc::grad(U)));
        //    X=magSqr(UPrime);
            turbDiff = -0.5*(UPrime*magSqr(UPrime));
          //  turbDiffTest = -0.5*(UPrime*UPrime);     //Turbulent Diffusion Term--divergence operator will be applied afterwards (uuu, uvv, uWw)
            pressDiff = -UPrime & fvc::grad(pPrime);           //Pressure Diffusion Term (make this universal), can be done in postLES
            SGSstrainTensor = symm(fvc::grad(UPrime));       //Tensor of Strain Rate of Resolved Fluctuations
          viscDiss = -2*nuLam*(SGSstrainTensor && SGSstrainTensor); //Viscous Dissipation of Resolved Fluctuations
            SGSDiff = -UPrime & G;                        //SGS Diffusion Term--divergence operator will be applied afterwards
             SGSDiss = B && SGSstrainTensor;   
             PU=  pPrime*UPrime;              //SGS Dissipation Term
             PD = -2*symm(fvc::grad(pPrime*UPrime)); //Pressure Diffusion (2 because symm func has half multiplied)
             PS = 2*pPrime*symm(fvc::grad(UPrime));  // Pressure-strain correlation... PD+PS= vel-pressure grad corre
           //  UUU= UPrime2Mean.component(0)*UMean.component(0);
             UUU = (UPrime.component(0)*UPrime.component(0)*UPrime.component(0));
             UUV = (UPrime.component(0)*UPrime.component(0)*UPrime.component(1));
             UUW = (UPrime.component(0)*UPrime.component(0)*UPrime.component(2));
             VVU = (UPrime.component(1)*UPrime.component(1)*UPrime.component(0));
             VVV = (UPrime.component(1)*UPrime.component(1)*UPrime.component(1));
             VVW = (UPrime.component(1)*UPrime.component(1)*UPrime.component(2));
             WWU = (UPrime.component(2)*UPrime.component(2)*UPrime.component(0));
             WWV = (UPrime.component(2)*UPrime.component(2)*UPrime.component(1));
             WWW = (UPrime.component(2)*UPrime.component(2)*UPrime.component(2));
             UVW = (UPrime.component(0)*UPrime.component(1)*UPrime.component(2));
           //  GradUP= (fvc::grad(UPrime));
             //Info<< "nulam = " << nuLam << " s";
             //tkeDiss= -2*nuLam*(fvc::grad(UPrime));
             tkeDiss=  -nuLam*cmptMultiply(GradUP,GradUP);//check this (should get uu,vv,and ww terms for dissipation at 0,4 and 8) 
            // X=nuLam;
            
             dissUV = -nuLam*(GradUP.component(tensor::XX)*GradUP.component(tensor::XY)+GradUP.component(tensor::YX)*GradUP.component(tensor::YY)+GradUP.component(tensor::ZX)*GradUP.component(tensor::ZY));
            
             dissUW = -nuLam*(GradUP.component(tensor::XX)*GradUP.component(tensor::XZ)+GradUP.component(tensor::YX)*GradUP.component(tensor::YZ)+GradUP.component(tensor::ZX)*GradUP.component(tensor::ZZ));
             dissVW = -nuLam*(GradUP.component(tensor::XY)*GradUP.component(tensor::XZ)+GradUP.component(tensor::YY)*GradUP.component(tensor::YZ)+GradUP.component(tensor::ZY)*GradUP.component(tensor::ZZ));
        // ..........................................Flow budgets_end....................................................// 
       
       // ..........................................TT budgets....................................................//
       // UPrime = (U-UMean);         //Resolved Velocity Fluctuation Vector
        TPrime = (T-TMean);
        
      //  pPrime = (p-pMean);
        TU= UPrime*TPrime; //production
        GradTP= (fvc::grad(TPrime));
        TTDiss=  cmptMultiply(GradTP,GradTP);
        TTU= UPrime*TPrime*TPrime;
	TTDissL= cmptMultiply((fvc::grad(T)),GradTP);  //Extra disspation term (LES)
	TTDiffL= TPrime*(fvc::grad(T));		 //Extra diffussion term (LES)


        //GradTTU= (fvc::grad(TTU));
        //GradUP= (fvc::grad(UPrime));
       
         // ..........................................TT budgets_end....................................................//
        
        // ..........................................heat flux budgets....................................................//
         GradTUP=( GradTP & GradUP ); // for dissipation
         GradTUP1=( GradUP & GradTP); //not correct
        GradTU= GradUP.component(tensor::XX)*GradTP.component(vector::X)+GradUP.component(tensor::YX)*GradTP.component(vector::Y)+GradUP.component(tensor::ZX)*GradTP.component(vector::Z);
       GradTV= GradUP.component(tensor::XY)*GradTP.component(vector::X)+GradUP.component(tensor::YY)*GradTP.component(vector::Y)+GradUP.component(tensor::ZY)*GradTP.component(vector::Z);
       GradTW=GradUP.component(tensor::XZ)*GradTP.component(vector::X)+GradUP.component(tensor::YZ)*GradTP.component(vector::Y)+GradUP.component(tensor::ZZ)*GradTP.component(vector::Z);
         UUT= UPrime*UPrime*TPrime; // for turbulent diffusion
         TGradp=-( TPrime*fvc::grad(pPrime)); // for Temp-Pressure diffusion
      //  Gradp= (fvc::grad(pPrime));
        MD1= TPrime*GradUP; //Molecular diffusion part 1
        MD2= UPrime*GradTP; 
        
        MD1L= TPrime*(fvc::grad(U)); //MD 1st term LES
        MD2L=UPrime*(fvc::grad(T));   //MD 2nd termLES
        GradTUPL= ( GradTP & (fvc::grad(U)));  //LES dissipation term
	GradTUPL2= ( (fvc::grad(T)) & GradUP);  //LES dissipation term
        
        
        
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
