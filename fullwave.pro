
Include "defines.pro";

Group{
    SkinFeed = Region[ SKINFEED ];
    SkinConductor = Region[ SKINCONDUCTOR ];
    SigmaInf = Region[ SIGMAINF ];
    Air = Region[ AIR ];
    Pml = Region[ PML ];
    Pec = Region[ PEC ];
    SurBC = Region[ {SkinFeed} ];
    DomainTot = Region[ {SkinFeed, Air, Pml, SigmaInf, Pec} ];
    DomainC = Region[ {} ];
    Domain = Region[ {Air, Pml} ];
    DefineGroup[DomainS];
    DefineGroup[SurS];
    TrGr = ElementsOf[ Domain , OnOneSideOf SkinFeed];
}

Function{
    mu0 = 1.25663706212e-06 ;
    nu0 = 795774.7150262763 ;
    ep0 = 8.8541878128e-12 ;
    gap = GetNumber["gap"] ;
    pml_delta = GetNumber["pml_delta"] ;
    air_boundary = GetNumber["air_boundary"] ;
    zl = 50.0 ;
    I[] = Complex[0.0, 1.0] ;
    epsilon[#{Air, SkinFeed, SigmaInf}] = ep0 ;
    nu[#{Air, SkinFeed, SigmaInf}] = nu0 ;
    sigma[] = 6.0e7 ;
    DefineFunction[js0];
    DefineFunction[ks0];
    DefineFunction[nxh];
    r[] = Sqrt[X[]^2 + Y[]^2 + Z[]^2] ;
    dumping_r[] = (r[] >= air_boundary) ? 1.0 / (pml_delta - (r[] - air_boundary)) : 0.0 ;
    cx[] = Complex[1.0, -dumping_r[] / k0] ;
    cy[] = Complex[1.0, -dumping_r[] / k0] ;
    cz[] = Complex[1.0, -dumping_r[] / k0] ;
    tens[] = TensorDiag[cy[] * cz[] / cx[],cx[] * cz[] / cy[],cx[] * cy[] / cz[]] ;
    epsilon[Pml] = ep0 * tens[] ;
    nu[Pml] = nu0 / tens[] ;
    x_feed = GetNumber["x_feed"];
    y_feed = GetNumber["y_feed"];
    r_xy[] = Sqrt[(X[] - x_feed)^2 + (Y[] - y_feed)^2] * gap;
    BC_Fct_e[] = Vector[(X[] - x_feed) / r_xy[], (Y[] - y_feed) / r_xy[], 0.0];
    dr[SkinFeed] = Vector[-(Y[] - y_feed) / r_xy[], (X[] - x_feed) / r_xy[], 0.0] ;
    // BC_Fct_e[] = Vector[1.0 / gap, 0.0, 0.0] ;
    // dr[SkinFeed] = Vector[0.0, 1.0, 0.0] ;
}

Constraint{
    { Name ElectricField; 
         Case  {  
           { Region SkinFeed; Type AssignFromResolution; NameOfResolution Microwave_e_BC; } 
           { Region SkinConductor; Type Assign; Value 0.0; } 
           { Region SigmaInf; Type Assign; Value 0.0; } 
           { Region Pec; Type Assign; Value 0.0; } 
         }
    }
}

FunctionSpace{
       { Name Hcurl_e; Type Form1; 
           BasisFunction {  
            {  Name se; NameOfCoef ee; Function BF_Edge; Support DomainTot; Entity EdgesOf[All]; } 
            }
           Constraint {  
            {  NameOfCoef ee; EntityType EdgesOf; NameOfConstraint ElectricField; } 
            }
       }
       { Name Hcurl_h; Type Form1; 
           BasisFunction {  
            {  Name sh; NameOfCoef he; Function BF_Edge; Support DomainTot; Entity EdgesOf[All]; } 
            }
       }
}

Jacobian{
    { Name JVol; 
         Case  {  
           { Region All; Jacobian Vol; } 
         }
    }
    { Name JSur; 
         Case  {  
           { Region All; Jacobian Sur; } 
         }
    }
}

Integration{
    { Name I1; 
         Case {  
         {  Type Gauss ; 
         Case  { 
            {  GeoElement Point; NumberOfPoints 1; } 
            {  GeoElement Line; NumberOfPoints 3; } 
            {  GeoElement Triangle; NumberOfPoints 4; } 
            {  GeoElement Quadrangle; NumberOfPoints 4; } 
            {  GeoElement Tetrahedron; NumberOfPoints 4; } 
            {  GeoElement Hexahedron; NumberOfPoints 6; } 
            {  GeoElement Prism; NumberOfPoints 6; }  
          }
    }
    }
    }
    { Name I2; 
         Case {  
         {  Type Gauss ; 
         Case  { 
            {  GeoElement Point; NumberOfPoints 1; } 
            {  GeoElement Line; NumberOfPoints 4; } 
            {  GeoElement Triangle; NumberOfPoints 7; } 
            {  GeoElement Quadrangle; NumberOfPoints 7; } 
            {  GeoElement Tetrahedron; NumberOfPoints 15; } 
            {  GeoElement Hexahedron; NumberOfPoints 34; } 
            {  GeoElement Prism; NumberOfPoints 21; }  
          }
    }
    }
    }
}

Formulation{
    { Name Microwave_e_BC; Type FemEquation; 
         Quantity {  
            {  Name e; Type Local; NameOfSpace Hcurl_e; } 
            }
         Equation {  
           Galerkin { [ Dof{e} , {e} ]; In SurBC; Integration I2; Jacobian JSur; } 
           Galerkin { [ -BC_Fct_e[] , {e} ]; In SurBC; Integration I2; Jacobian JSur; } 
            }
    }
    { Name Microwave_e; Type FemEquation; 
        Quantity {  
            {  Name e; Type Local; NameOfSpace Hcurl_e; } 
            {  Name h; Type Local; NameOfSpace Hcurl_h; } 
        }
        Equation {  
           Galerkin { [ nu[] * Dof{d e} , {d e} ]; In Domain; Integration I1; Jacobian JVol; } 
           Galerkin { DtDof[ sigma[] * Dof{e}, {e} ]; In DomainC; Integration I1; Jacobian JVol; } 
           Galerkin { DtDtDof[ epsilon[] * Dof{e} , {e} ]; In Domain; Integration I1; Jacobian JVol; } 
           Galerkin { [ Dof{h} , {h} ]; In TrGr; Integration I1; Jacobian JVol; } 
           Galerkin { [ -I[] * nu[] * Dof{d e} / (2.0 * Pi * freq), {h} ]; In TrGr; Integration I1; Jacobian JVol; } 
        }
    }
}

Resolution{
    { Name Microwave_e_BC; Hidden 1; 
         System {  
            { Name B; NameOfFormulation Microwave_e_BC; DestinationSystem A; } 
            }
         Operation {  
           Generate[B]; Solve[B]; TransferSolution[B];  }
    }
    { Name Analysis; 
         System {  
            { Name A; NameOfFormulation Microwave_e; Type Complex; Frequency freq; } 
            }
         Operation {  
           CreateDirectory['build']; Generate[A]; Solve[A]; SaveSolution[A];  }
    }
}

PostProcessing{
    { Name Microwave_e; NameOfFormulation Microwave_e; 
         Quantity {  
            {  Name e; Value  {  Local {  [ {e} ]; In DomainTot; Jacobian JVol; } } }
            {  Name h; Value  {  Local {  [ I[] * nu[] * {d e} / (2.0 * Pi * freq) ]; In Domain; Jacobian JVol; } } }
            {  Name y; Value  {  Integral {  [ {h} * dr[] ]; In SkinFeed; Jacobian JSur; Integration I2; } } }
            {  Name exh ; Value{ Local{ [CrossProduct[ {e}, Conj[I[] * nu[] * {d e} / (2.0 * Pi * freq)]]]; In Domain;  Jacobian JVol;} } }
            // {  Name e_norm; Value  {  Local { [Norm[{e}]]; In Domain; Jacobian JVol; } } }
            // {  Name h_norm; Value  {  Local { [Norm[{h}]]; In Domain; Jacobian JVol; } } }
            {  Name s11; Value  { Local { Type Global; [ 20.0 * Log10[Norm[(1.0 - zl * $y) / (1.0 + zl * $y)]] ]; In SkinFeed; } } }
          }
    }
}

PostOperation{
    { Name Microwave_e; NameOfPostProcessing Microwave_e; 
         Operation {  
            // Print [ e, OnElementsOf Region[{Domain}], File "./build/e.pos" ]; 
            // Print [ h, OnElementsOf Region[{Domain}], File "./build/h.pos" ]; 
            // Print [ exh, OnElementsOf Region[{Domain}], File "./build/exh.pos" ]; 
            // Print [ e_norm, OnSection {{0.0, 0.0, 0.0} {1.0, 0.0, 0.0} {0.0, 1.0, 0.0}}, File "./build/e_norm.pos" ]; 
            // Print [ e_norm, OnSection {{0.0, 0.0, 0.01} {1.0, 0.0, 0.01} {0.0, 1.0, 0.01}}, File "./build/e_norm1.pos" ]; 
            // Print [ e_norm, OnSection {{0.0, 0.0, -0.01} {1.0, 0.0, -0.01} {0.0, 1.0, -0.01}}, File "./build/e_norm2.pos" ]; 
            // Print [ e_norm, OnGrid SkinConductor, File "./build/e_norm3.pos" ]; 
            Print [ y[SkinFeed], OnGlobal , Format FrequencyTable, StoreInVariable $y, File "./build/y.txt" ]; 
            Print [ s11, OnRegion SkinFeed, Format FrequencyTable, StoreInVariable $s11, SendToServer "s11", File "./build/s11.txt" ]; 
            }
    }
}
    // # poi0.add(
    //   #     'e', OnSection='{{0.0, 0.0, 0.0} {1.0, 0.0, 0.0} {0.0, 1.0, 0.0}}', File='./build/e_norm.pos')
      