/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | Copyright (C) 2018 Ehsan Madadi-Kandjani        |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'    ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([{,}])
define(calc, [{esyscmd(perl -e 'use Math::Trig; print ($1)')}])
define(VCOUNT, 0)
define(vlabel, [{[{// }]Vertex $1 = VCOUNT define($1, VCOUNT)define([{VCOUNT}], incr(VCOUNT))}])

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(btQuad, ($1b $2b $2t $1t))
define(topQuad, ($1t $4t $3t $2t))
define(bottomQuad, ($1b $2b $3b $4b))

// Returns a zero array
define(
         zeroArray, 
         [{
             esyscmd(perl -e 'use Math::Trig;
             use 5.010;
             for $i (0 .. $1)
             {
                @array[$i] = 0
             };
             print "(@array)"; ')
         }]
      )

// Returns a zero array with commas in the middle of it
define(
         zeroArrayCommas, 
         [{ 
            esyscmd(perl -e 'use Math::Trig;
            use 5.010;
            $\="\n";
            for $i (0 .. $1)
            {
               @array[$i] = $2
            };
            {
               local $" = [{", "}];
               print "(@array)"
            }')}]
      )

// Returns a linear  array with commas (initial point, last point, number of segments, c)
define(
         XComma, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $3)
             {
                 @arrayBeta[$i] = $1 + (($2-$1)/$3)*($i);
                 @arrayX[$i] = ($4)*0.5*(1.0 - cos(@arrayBeta[$i]))
             };
             {
               local $" = [{", "}];
               print "(@arrayX)"
            } ')
         }]
      )
  
// Returns a y_t with commas (arrayX, number of segments, t, c)
define(
         yt, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $2)
             {
                 @arrayYt[$i] = (5.0*$3*$4)*(0.2969*($1[$i]/$4)**0.5 - 
0.1260*($1[$i]/$4) - 0.3516*($1[$i]/$4)**2.0 + 0.2843*($1[$i]/$4)**3.0 - 
0.1015*($1[$i]/$4)**4.0)
             };
             {
               local $" = [{", "}];
               print "(@arrayYt)"
            } ')
         }]
      )

// Returns a y_c with commas (arrayX, number of segments, m, p, c)
define(
         yc, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $2)
             {
                 if ($1[$i] < $4*$5)
                 {
                    @arrayYc[$i] = ($3*$5/$4**2)*(2*$4*($1[$i]/$5) - ($1[$i]/$5)**2)
                 }
                 else
                 {
                    @arrayYc[$i] = ($3*$5/((1.0 - $4)**2))*((1.0 - 2.0*$4) + 
2.0*$4*($1[$i]/$5) - ($1[$i]/$5)**2)
                 };
             };
             {
               local $" = [{", "}];
               print "(@arrayYc)"
            } ')
         }]
      )

// Returns a theta with commas (arrayX, number of segments, m, p, c)
define(
         theta, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $2)
             {
                 if ($1[$i] < $4*$5)
                 {
                    @arraydYcdX[$i] = (2.0*$3/$4**2)*($4 - ($1[$i]/$5))
                 }
                 else
                 {
                    @arraydYcdX[$i] = (2.0*$3/(1.0 - $4)**2)*($4 - ($1[$i]/$5))
                 };
                 @arrayTheta[$i] = atan(@arraydYcdX[$i])
             };
             {
               local $" = [{", "}];
               print "(@arrayTheta)"
            } ')
         }]
      )
      
// Returns a xu with commas (arrayX, arrayYt, arrayTheta, number of segments)
define(
         xu, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayXu[$i] = $1[$i] - $2[$i]*sin($3[$i])
             };
             {
               local $" = [{", "}];
               print "(@arrayXu)"
            } ')
         }]
      )

// Returns a yu with commas (arrayYc, arrayYt, arrayTheta, number of segments)
define(
         yu, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayYu[$i] = $1[$i] + $2[$i]*cos($3[$i])
             };
             {
               local $" = [{", "}];
               print "(@arrayYu)"
            } ')
         }]
      )
      
// Returns a xl with commas (arrayX, arrayYt, arrayTheta, number of segments)
define(
         xl, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayXl[$i] = $1[$i] + $2[$i]*sin($3[$i])
             };
             {
               local $" = [{", "}];
               print "(@arrayXl)"
            } ')
         }]
      )

// Returns a yl with commas (arrayYc, arrayYt, arrayTheta, number of segments)
define(
         yl, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayYl[$i] = $1[$i] - $2[$i]*cos($3[$i])
             };
             {
               local $" = [{", "}];
               print "(@arrayYl)"
            } ')
         }]
      )
      
// Include angle of attack (arrayXu, arrayYu, alpha, number of segments)
define(
         xuAlpha, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i (0 .. $4)
             {
                 @arrayXuAlpha[$i] = $1[$i]*cos($3) + $2[$i]*sin($3)
             };
             print "(@arrayXuAlpha)"; ')
         }]
      )
      
define(
         yuAlpha, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i (0 .. $4)
             {
                 @arrayYuAlpha[$i] = -1.0*$1[$i]*sin($3) + $2[$i]*cos($3)
             };
             print "(@arrayYuAlpha)"; ')
         }]
      )
      
// Include angle of attack (arrayXl, arrayYl, alpha, number of segments)
define(
         xlAlpha, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i (0 .. $4)
             {
                 @arrayXlAlpha[$i] = $1[$i]*cos($3) + $2[$i]*sin($3)
             };
             print "(@arrayXlAlpha)"; ')
         }]
      )
      
define(
         ylAlpha, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i (0 .. $4)
             {
                 @arrayYlAlpha[$i] = -1.0*$1[$i]*sin($3) + $2[$i]*cos($3)
             };
             print "(@arrayYlAlpha)"; ')
         }]
      )
      
// Xu with comma to obtain Cmax (arrayXu, arrayYu, alpha, number of segments)
define(
         xuAlphaC, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayXuAlphaC[$i] = $1[$i]*cos($3) + $2[$i]*sin($3)
             };
             {
               local $" = [{", "}];
               print "(@arrayXuAlphaC)"
             } ')
         }]
      )

define(
         yuAlphaC, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayYuAlphaC[$i] = -1.0*$1[$i]*sin($3) + $2[$i]*cos($3)
             };
             {
               local $" = [{", "}];
               print "(@arrayYuAlphaC)"
             } ')
         }]
      )
      
// Include angle of attack (arrayXl, arrayYl, alpha, number of segments)
define(
         xlAlphaC, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayXlAlphaC[$i] = $1[$i]*cos($3) + $2[$i]*sin($3)
             };
             {
               local $" = [{", "}];
               print "(@arrayXlAlphaC)"
             } ')
         }]
      )
      
define(
         ylAlphaC, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             $\="\n";
             for $i (0 .. $4)
             {
                 @arrayYlAlphaC[$i] = -1.0*$1[$i]*sin($3) + $2[$i]*cos($3)
             };
             {
               local $" = [{", "}];
               print "(@arrayYlAlphaC)"
             } ')
         }]
      )
      
// Returns the index of maximum array (Yc, Yt, p, number of segments)
define(
         maxIndex, 
         [{
            esyscmd(perl -e 'use Math::Trig;
            use 5.010;
            $index = 0;
            if ($3 > 0)
            {
                $maxval = $1[ $index ];
                for $i (0 .. $4)
                {
                    if ( $maxval < $1[$i] )
                    {
                        $index = $i;
                        $maxval = $1[$i];
                    }
                }
            }
            else
            {
                $maxval = $2[ $index ];
                for $i (0 .. $4)
                {
                    if ( $maxval < $2[$i] )
                    {
                        $index = $i;
                        $maxval = $2[$i];
                    }
                }
            }
            print "$index";')
         }]
      )
      
define(
         spl, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i ($2 .. $3)
             {
                 @arrayXlAlpha[$i] = $1[$i]
             };
             print "(@arrayXlAlpha)"; ')
         }]
      )
      
define(
         splMerge, 
         [{
             esyscmd(perl -e 'use Math::Trig; 
             use 5.010;
             for $i ($4 .. $5)
             {
                 push(@merged, ("(", $1[$i], $2[$i], $3[$i], ")"))
             };
             print "(@merged)"; ')
         }]
      )

      
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1.0;

//NACA airfoil geometry

//Chord length
define(c, 1.0)

//Angle of attack (in radian)
define(alpha, 0.3)

//NACA digits
define(NACA1, 4)
define(NACA2, 4)
define(NACA3, 1)
define(NACA4, 2)

//Mesh

//Height of mesh in y direction
define(L1, 8.0)

//Length of downstream
define(L2, 16.0)

//Foil depth in z direction perpendicular to x-y surface
define(L3, 0.5)

//Scale factor (the mesh is normalized to have chord length equal to one)
define(sc, 1.0)

//Number of cells and points at each direction and element

//Number of cells in y direction
define(Nl1, 400)

//Number of cells in downstream
define(Nl2, 1500)

//Number of cells in z direction
define(Nl3, 1)

//Number of meshes on the front part of airfoil edges p8-p9 and p8-p10b
define(Nl4, 100)

//Number of meshes on the back part of airfoil edges p9-p11 and p10-p11
define(Nl5, 70)

//Number of interpolation points along the airfoil for defining the splines
define(Naf, 99)

//Cell expansion ratios

//Expansion ratio in y direction
define(E1, 100)

//Expansion ratio in downstream side
define(E2, 100)

//Expansion ratio in inlet
define(E3, 10)

// Base z
define(Zb, 0)

// Depth of airfoil
define(z, 0.01)

// Front z
define(Zt, calc(Zb + z))

// m,p, and t values for the airfoil
define(m, calc(NACA1/100))
define(p, calc(NACA2/10))
define(t, calc((NACA3*10+NACA4)/100))

// Calculate x array
define(ax, XComma(0, pi, Naf, c))

// Calculate yt array
define(ayt, yt(ax, Naf, t, c))

// Calculate camber
define(ayc, yc(ax, Naf, m, p, c))

// Calculate theta
define(atheta, theta(ax, Naf, m, p, c))

// Calculate upper surface x components
define(axu, xu(ax, ayt, atheta, Naf))

// Calculate upper surface y components
define(ayu, yu(ayc, ayt, atheta, Naf))

// Calculate lower surface x components
define(axl, xl(ax, ayt, atheta, Naf))

// Calculate lower surface y components
define(ayl, yl(ayc, ayt, atheta, Naf))

define(axuAlphaC, xuAlphaC(axu, ayu, alpha, Naf))
define(ayuAlphaC, yuAlphaC(axu, ayu, alpha, Naf))
define(axlAlphaC, xlAlphaC(axl, ayl, alpha, Naf))
define(aylAlphaC, ylAlphaC(axl, ayl, alpha, Naf))

// define Cmax
define(Cmax, maxIndex(ayc, ayt, p, Naf))

// Move the airfoil
define(noseX, calc((-1*L1 + axuAlphaC[Cmax])*cos(alpha)))
define(noseY, calc((L1 - axuAlphaC[Cmax])*sin(alpha)))

define(x00, noseX)
define(x01, calc(axlAlphaC[Cmax]))
define(x02, calc(axuAlphaC[Naf]))
define(x03, L2)
define(x04, L2)
define(x05, L2)
define(x06, calc(axuAlphaC[Naf]))
define(x07, calc(axuAlphaC[Cmax]))
define(x08, 0)
define(x09, calc(axuAlphaC[Cmax]))
define(x10, calc(axlAlphaC[Cmax]))
define(x11, calc(axuAlphaC[Naf]))

define(y00, noseY)
define(y01, calc(-1*L1))
define(y02, calc(-1*L1))
define(y03, calc(-1*L1))
define(y04, calc(aylAlphaC[Naf]))
define(y05, L1)
define(y06, L1)
define(y07, L1)
define(y08, 0)
define(y09, calc(ayuAlphaC[Cmax]))
define(y10, calc(aylAlphaC[Cmax]))
define(y11, calc(aylAlphaC[Naf]))

//Define the points for the inlet arc allocation
define(p12, 45)
define(p13, 45)

define(cp12, calc(cos((pi/180)*p12)))
define(cp13, calc(cos((pi/180)*p13)))

define(sp12, calc(sin((pi/180)*p12)))
define(sp13, calc(sin((pi/180)*p13)))

define(x12, calc(-1*L1*cp12 + axuAlphaC[Cmax]))
define(x13, calc(-1*L1*cp13 + axuAlphaC[Cmax]))

define(y12, calc(L1*sp12))
define(y13, calc(-1*L1*sp13))

define(azbAlphaC, zeroArrayCommas(Naf, Zb))
define(azAlphaC, zeroArrayCommas(Naf, z))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

define(vert, (x$1$2 y$1$2 $3))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(p00b)
    vert(0, 1, Zb) vlabel(p01b)
    vert(0, 2, Zb) vlabel(p02b)
    vert(0, 3, Zb) vlabel(p03b)
    vert(0, 4, Zb) vlabel(p04b)
    vert(0, 5, Zb) vlabel(p05b)
    vert(0, 6, Zb) vlabel(p06b)
    vert(0, 7, Zb) vlabel(p07b)
    vert(0, 8, Zb) vlabel(p08b)
    vert(0, 9, Zb) vlabel(p09b)
    vert(1, 0, Zb) vlabel(p10b)
    vert(1, 1, Zb) vlabel(p11b)
    
    vert(0, 0, Zt) vlabel(p00t)
    vert(0, 1, Zt) vlabel(p01t)
    vert(0, 2, Zt) vlabel(p02t)
    vert(0, 3, Zt) vlabel(p03t)
    vert(0, 4, Zt) vlabel(p04t)
    vert(0, 5, Zt) vlabel(p05t)
    vert(0, 6, Zt) vlabel(p06t)
    vert(0, 7, Zt) vlabel(p07t)
    vert(0, 8, Zt) vlabel(p08t)
    vert(0, 9, Zt) vlabel(p09t)
    vert(1, 0, Zt) vlabel(p10t)
    vert(1, 1, Zt) vlabel(p11t)
);

blocks
(
    //B0
    hex2D(p00, p01, p10, p08)
    square
    (Nl4 Nl1 Nl3)
    simpleGrading (1 1 1)
    
    //B1
    hex2D(p01, p02, p11, p10)
    square
    (Nl5 Nl1 Nl3)
    simpleGrading (1 1 1)
    
    //B2
    hex2D(p02, p03, p04, p11)
    square
    (Nl2 Nl1 Nl3)
    simpleGrading (1 1 1)
    
    //B3
    hex2D(p11, p04, p05, p06)
    square
    (Nl2 Nl1 Nl3)
    simpleGrading (1 1 1)
    
    //B4
    hex2D(p09, p11, p06, p07)
    square
    (Nl5 Nl1 Nl3)
    simpleGrading (1 1 1)
    
    //B5
    hex2D(p08, p09, p07, p00)
    square
    (Nl4 Nl1 Nl3)
    simpleGrading (1 1 1)
);

edges
(
    //Inlet arc
    arc p00b p07b vert(1, 2, Zb)
    arc p00b p01b vert(1, 3, Zb)

    arc p00t p07t vert(1, 2, Zt)
    arc p00t p01t vert(1, 3, Zt)
    
    //Airfoil surface spline
    spline p08b p10b splMerge(axlAlphaC, aylAlphaC, azbAlphaC, 1, Cmax - 1)
    spline p10b p11b splMerge(axlAlphaC, aylAlphaC, azbAlphaC, Cmax + 1, Naf - 1)
    spline p08b p09b splMerge(axuAlphaC, ayuAlphaC, azbAlphaC, 1, Cmax - 1)
    spline p09b p11b splMerge(axuAlphaC, ayuAlphaC, azbAlphaC, Cmax + 1, Naf - 1)
    
    spline p08t p10t splMerge(axlAlphaC, aylAlphaC, azAlphaC, 1, Cmax - 1)
    spline p10t p11t splMerge(axlAlphaC, aylAlphaC, azAlphaC, Cmax + 1, Naf - 1)
    spline p08t p09t splMerge(axuAlphaC, ayuAlphaC, azAlphaC, 1, Cmax - 1)
    spline p09t p11t splMerge(axuAlphaC, ayuAlphaC, azAlphaC, Cmax + 1, Naf - 1)
    
);

patches
(
    wall airfoil
    (
        btQuad(p08, p10)
        btQuad(p10, p11)
        btQuad(p08, p09)
        btQuad(p09, p11)
    )
    
    patch farBoundaries
    (
        btQuad(p06, p07)
        btQuad(p05, p06)
        btQuad(p02, p03)
        btQuad(p01, p02)
    )
    
    patch inlet
    (
        btQuad(p07, p00)
        btQuad(p00, p01)
    )
    
    patch outlet
    (
        btQuad(p04, p05)
        btQuad(p03, p04)
    )
);

mergePatchPairs
(
);
