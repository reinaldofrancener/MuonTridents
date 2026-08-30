#include "functions.h"
#include "globais.h"

#include <iostream>
#include <cmath>



//*****************************************************************
// squared matrix elements for lepton tridents (Bethe-Heitler + bremss with interference)
//*****************************************************************

double SquaredMatrixElementPLP_BHbremss_interf(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    double m32,m12;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;
    m32 = m3*m3;
    m12 = m1*m1;

      result = (4*M2*m32*m32*p1p2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*m32*p1p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*m32*p1p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1P*p2P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*m32*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*p1p2*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*p1p3*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p4*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p2p3*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*m32*p1p2*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p4*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1P*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p2*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p4*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*p1p3*p2p4*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1P*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p2*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*m32*p1p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p3*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p4*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p2*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p3*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p4*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p4*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*p1p2*p3p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2P*p3p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*m32*p1p2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*m32*p1p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*m32*p1p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*m32*p1P*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p2*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*m32*p1p3*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p4*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*p1p4*p2p3*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p2*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p1p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*m32*p1p4*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p2*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p3*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p4*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2P*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*p1p3*p2p4*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1P*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1P*p2p3*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p4*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*m32*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p2*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*p1p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*p1p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p2*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p3*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2P*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p2*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*p1p4*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2P*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p3P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p3P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p3p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*m32*p1p2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*m32*p1p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*m32*p1P*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*m32*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p2*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p4*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p2p3*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*M2*m32*m32*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*p1p2*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*p1p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*p1p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 32*m32*p1P*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p2*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p4*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*M2*m32*p2p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*p1p3*p2p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 32*m32*p1P*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*p1P*p2p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*M2*m32*m32*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p2*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p1p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*m32*p2p3*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*M2*m32*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p2*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*M2*p1p4*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*M2*m32*p3p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2p4*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1P*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p3*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1P*p2p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p3p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p4*p3p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 32*m32*p1P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*p1P*p2p4*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 m32*m32*p1p2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*m32*p1p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*m32*p1p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*m32*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*p1p2*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*p1p3*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p4*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p2p3*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*m32*p1p2*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p4*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p2*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p4*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 p1p3*p2p4*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p2*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*m32*p1p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p3*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p4*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p2*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p3*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p4*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 p1p2*p3p4*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*m32*p1p2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*m32*p1p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*m32*p1p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p2*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*m32*p1p3*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p4*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 p1p4*p2p3*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p2*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p1p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*m32*p1p4*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p2*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p3*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p4*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 p1p3*p2p4*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*m32*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p2*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*p1p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*p1p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p2*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p3*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p2*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 p1p4*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p3p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*m32*p1p2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*m32*p1p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*m32*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p2*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p4*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p2p3*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 4*m32*m32*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*p1p2*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*p1p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*p1p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p2*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p4*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 m32*p2p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*p1p3*p2p4*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*m32*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p2*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p1p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*m32*p2p3*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 2*m32*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p2*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 2*p1p4*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 m32*p3p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 4*m32*m32*p1p2*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*p1p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*m32*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1p2*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1p3*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p4*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p3*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1p2*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p4*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p2*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p4*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1p3*p2p4*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p2*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1p3*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p3*p2p3*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p4*p2p3*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p2*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p3*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p4*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1p2*p3p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1p2*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*p1p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*p1p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p2*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1p3*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p4*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1p4*p2p3*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p2*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1p4*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p2*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p3*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p4*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1p3*p2p4*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*m32*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p2*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1p3*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1p4*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p2*p2p3*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p3*p2p3*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p2*p2p4*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1p4*p2p4*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p3p4*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1p2*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*m32*p1p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*m32*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p2*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p4*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p3*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*m32*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1p2*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*p1p4*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p2*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p4*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p4*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1p3*p2p4*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 16*m32*m32*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p2*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p3*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p2*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1p4*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p3p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p2P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2P*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2P*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2P*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p3*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p4*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2P*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p4*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2P*p3p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p2P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2P*p2p3*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2P*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2P*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p3*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p4*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2P*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2P*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p2P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2P*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p2P*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2P*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p2p4*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2P*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2p4*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p3*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p2p4*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p3*p2p4*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p4*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2p3*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p2p4*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32)/((m32 + p3p4)*(m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 2*m32*p1p2)/((m32 + p3p4)*(m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 2*p1p4*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 2*p1p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 2*m32*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 64*m32*m32*p2P*p2P)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*m32*p1p2*p2P*p2P)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p4*p2P*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p3*p2P*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*p2P*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*M2*m32*m32)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*M2*m32*p1p2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*m32*p2P)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*p1P*p2P)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*m32*p1p2*p2P)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*M2*p1p4*p2p3)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p4*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*M2*p1p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*p1p4*p2P*p3P)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*M2*m32*p3p4)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*p1p3*p2P*p4P)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 64*m32*m32*p2P*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*p1p2*p2P*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*m32*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*m32*m32*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*m32*p1p2*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p4*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p3*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*p2P*p2P*qp1)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*m32*p2P*qp1)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*m32*p2P*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*m32*p1P*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*m32*p2P*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1p3*p2P*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1p4*p2P*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 32*m32*p1P*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p4*p2P*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 32*m32*p1P*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p3*p2P*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 32*p1P*p2P*p2p3*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 32*m32*m32*p2P*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p1p2*p2P*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p4*p2P*p2p3*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p2*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p4*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1P*p2P*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1P*p2P*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p3P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*m32*p1P*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*m32*p2P*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 24*m32*p1p3*p2P*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1p4*p2P*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1P*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 40*m32*p1P*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p3*p2P*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1P*p2P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 24*m32*m32*p2P*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 24*m32*p1p2*p2P*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1p4*p2P*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p2p3*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 24*m32*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p2*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p4*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1P*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1P*p2P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 32*m32*m32*p2P*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p1p2*p2P*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p2P*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p2*p2P*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p3*p2P*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p3*p2P*p2p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p3p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*m32*p1p3*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*m32*p1p4*p2P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p2p3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 12*m32*m32*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p4*p2p3*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2p4*p3P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*p1p2*p2P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p3*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p4*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p3P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*p1p2*p3P*p3p4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p1p2*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 12*m32*p1p3*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p4*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*m32*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p1p3*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*p1p2*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 12*m32*m32*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p4*p2p3*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*p1p2*p3p4*p4P)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*m32*m32*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*m32*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*m32*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 16*m32*p1p2*p2P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p1p3*p2P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p1p4*p2P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 24*m32*m32*p3P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p3P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2p3*p3P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p4*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p2*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1P*p2p3*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1P*p2p4*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1p2*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 24*m32*p1p3*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1p3*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p4*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1P*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 24*m32*m32*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p1p3*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*m32*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*m32*p2P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p2P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p2p3*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*M2*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p2P*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2p4*p3P*p3P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p2P*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*p2P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2P*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 20*m32*p2P*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p2P*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2P*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p3P*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p2P*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2p3*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2p4*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2P*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2p3*p4P*p4P*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*m32*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 5*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 8*m32*m32*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2p3*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*p2p3*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 12*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p2p3*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*p1P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1P*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1P*p3p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p1P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p4*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p1P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*m32*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*p2p3*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p4*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*p3p4*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p3*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p3*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p4*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p3*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p3*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p4*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 64*m32*m32*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*m32*p1P*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*p1p2*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*p1p4*p2P*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*p1p3*p2P*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p4*p2P*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*m32*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*m32*m32*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*p1P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*m32*p1p2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*p1p4*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*p1p3*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 32*m32*m32*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*m32*p1p2*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*M2*m32*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*m32*p2P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*m32*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*m32*p1P*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*p1p2*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*p1p4*p2p3*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*p1p3*p2p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*p1p4*p3P*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*m32*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 16*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*M2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*m32*p1P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*M2*m32*p1p2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*M2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*M2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*m32*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*p1P*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*p1p3*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*p1p4*p2P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*m32*p1P*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*M2*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*p1p4*p2P*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*m32*p1P*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*M2*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*p1p3*p2P*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*p1P*p2p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*m32*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 24*m32*p1P*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p1p2*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p4*p2p3*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p2*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p4*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p4*p3P*p3P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*M2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p1P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*M2*p1p2*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*p1P*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1P*p2p3*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1P*p2p4*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*m32*p1P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*M2*m32*p1p2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 12*M2*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*M2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*m32*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p1P*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 24*m32*p1p3*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p1p4*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1P*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*M2*p1p4*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 20*m32*p1P*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 16*p1p3*p2P*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1P*p2p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 12*m32*m32*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p1P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 12*m32*p1p2*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1p4*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p4*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p2p3*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 12*m32*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p2*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p4*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1P*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*M2*p1p2*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1P*p2p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*m32*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 24*m32*p1P*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p1p2*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p2*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p3*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p3*p2p4*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p3*p3P*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p4*p3P*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p3p4*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p3p4*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1p2*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p2p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1p4*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p3*p4P*p4P*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) + (
 2*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) + (
 2*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 m32*m32*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 3*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) + (
 m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) + (
 2*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) - (
 m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)) + (
 8*m32*m32*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p1p2*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1p2*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1p2*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p4*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*m32*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p2P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p2p3*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*m32*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p2p3*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p2P*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p2p3*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p2p4*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*m32*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p2P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p2p3*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 10*m32*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p2p3*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p2P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p2p4*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p2p3*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 12*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p2p3*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*m32*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p1P*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1p3*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p4*p2p3*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*p1P*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*m32*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 6*m32*p1p3*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 2*m32*p1p4*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2p4*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p4*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1P*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 16*m32*m32)/((q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4)/((q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3)/((q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4)/((q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p3p4)/((q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*m32*p2P*p2P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p1p4*p2P*p2P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p2P*p2P*p2p3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*p2P*p2p4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p2*p2P*p2P*p3p4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*m32*m32)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*m32*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p4*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1P*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p2*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p3*p2P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 256*m32*m32*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p1p4*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p2P*p2p3*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p2*p2P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p2P*p2P*p3p4*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2P*p3p4*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p2P*p3p4*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 256*m32*m32*p2P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p1p4*p2P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p2P*p2p3*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p2P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1P*p2P*p3p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p2*p2P*p3p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*m32*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p2*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p2P*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2p4*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p3P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1P*p3p4*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p3p4*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p3*p4P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p3p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*m32*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p1p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p2p3*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p2*p3P*p3P*p3p4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*m32*m32)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*m32*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2p3*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1P*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p2*p3P*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p2*p3P*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 256*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p1p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p2*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p2p4*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p2p4*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*qp2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p3p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (512*m32*m32*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*p1p4*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (256*m32*p2P*p2p3*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*p1p3*p2P*p2p4*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*p1p2*p2P*p3P*p3p4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*M2*m32*p1p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p1p4*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*M2*m32*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p2P*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p1p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p4*p2P*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p2p3*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1P*p2p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1P*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p1P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1P*p2p3*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*m32*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p1p4*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2P*p2p3*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p1p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p2*p2P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p2*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*m32*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p1p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p2p3*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1p2*p3p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p4*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p3*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p2*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p2P*p2p4*p3P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p2P*p3P*p3p4*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p2P*p3p4*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p2p3*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*m32*m32*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p1p4*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2p3*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p2p4*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1P*p3P*p3p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p2*p3P*p3p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p3*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*m32*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*m32*p1p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*m32*p2p3*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p3*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*p1P*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p2*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p1p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p2p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p4*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1P*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p2p4*p3P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*M2*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p3P*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p3P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p2)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p3)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p3)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p3p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*m32*p2P*p2P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p2*p2P*p2P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p3*p2P*p2P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p4*p2P*p2P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2P*p2p3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p2P*p2p4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*p2P*p2p4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2P*p3p4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*m32*m32)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p1p2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p1p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p1p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p3*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p3*p2P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*m32*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p2*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p3*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p4*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p2p3*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*p1P*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*p2P*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p2P*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1P*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1P*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2P*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2P*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1P*p2P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p2P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p1P*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p2P*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p2P*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p2P*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 80*m32*p1P*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2P*p2p3*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*m32*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*p1p4*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2P*p2p4*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p2*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*p1p3*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*p1P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*m32*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p2*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p4*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p1P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p4*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p2P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2P*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*p2P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 40*m32*p2P*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p4P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 10*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*m32*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1P*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p2*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p3*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p4*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p2p3*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2p4*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2p4*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p3p4*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p4P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3p4*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p4P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*m32*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p4*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2p4*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3p4*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p4P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p1P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p3*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p3*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1P*p2p4*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p1P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 80*m32*p1P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p3*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*m32*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*p1p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p4*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*p1p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*p1P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*m32*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2p4*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*p1P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p3*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*m32*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p1p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p1P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p1P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 48*m32*m32*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 48*m32*p1p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p4*p3P*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2P*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p2p4*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*M2*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 40*m32*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p3*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p4*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p3*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p4P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 10*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p3*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p3*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p3*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p2)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p3)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p3)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p3p4)/((m32 + p3p4)*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*m32*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p2*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p3*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1p4*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p2p3*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2p4*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2p4*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p3P*p3p4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p1p2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p1p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p1p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*p1p4*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2P*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p2*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p3*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p3*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*m32*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2p4*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p2*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p3*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p4*p2p3*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p3P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p4P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p4P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2p4*p3P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p4P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2p4*p3P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p2P*p3p4*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2p3*p4P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*m32*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p4*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p2P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p4*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p3P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p3p4*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p4P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p4P*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p4P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p4P*qp2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p3P*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p4P*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p4P*qP*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*p3p4*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp1*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p3P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p4P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp2*qp2)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*m32*p1P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1P*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1P*p2p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p4*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p1P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p2p3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1P*p2p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1P*p3p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p3p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*p1P*p2p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*m32*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*p2p3*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*p2p4*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p4*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*m32*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p3*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*p2p4*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p2p4*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 32*p1p4*p2P*p2P*qp3)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p4*p2P*qp3)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p4*p2P*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*M2*m32*m32*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*m32*p2P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p2P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2P*p2p3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p3*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p4*p2P*p3P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1P*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p2*p3p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p3p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*M2*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*m32*p1p3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*M2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 12*m32*m32*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 12*m32*p1p2*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p1p4*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p2p3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 12*m32*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*p3P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p3P*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3P*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p4*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p4P*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p4P*p4P*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p1p4*p2p3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 3*m32*m32*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 3*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 3*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*m32*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p4*p2p3*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p3*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p3p4*qP*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1p2*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p3*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 12*m32*p2P*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p2p4*p3P*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2P*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p3P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2p3*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*M2*p1p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p4*p2P*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*m32*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1P*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1p2*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p4*p2p3*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p2*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p4*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p4*p3P*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1P*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 6*m32*m32*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 8*m32*p1P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 6*m32*p1p2*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 2*m32*p1p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p4*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p2p3*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 6*m32*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1P*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p2*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p4*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p4P*qp2*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (128*m32*p2P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*M2*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 256*m32*m32*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p1p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*m32*p2p3*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1P*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p2*p3P*p3p4*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p2*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*m32*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2p3*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p4P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p1p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p2*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*p2p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p2p4*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*M2*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 256*m32*m32*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p1p4*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*p2P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2P*p2p3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1P*p2P*p2p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p2P*p2p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p2*p2P*p3p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p2*p2P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*m32*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*m32*p1p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*m32*p2p3*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p2*p3p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (32*p1p2*p4P*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p1p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1p4*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1P*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p2*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*M2*p2p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p2P*p3p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p2P*p4P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*m32*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*M2*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p2P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*m32*p2p3*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1P*p2p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*m32*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1P*p3p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*p1P*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1p2*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (64*p1p3*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p2*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p1p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*m32*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p1p4*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p4*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p4P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p2*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 6*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 6*m32*p1p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 6*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p2P*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p2*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1p3*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p2p3*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*p3P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*p3p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p4P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*p4P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p1P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p4*p2P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p4P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p4P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*m32*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p1p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2p3*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p4*p2P*qp3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1P*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p1P*qp3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p4P*qp3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p2*p4P*qp3*qp3)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (64*m32*m32*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2p3*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1P*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p3*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p3p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p1p2*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p4P*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*qp3*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*qp3*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p4P*qp3*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*m32*p1P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*m32*m32*p2P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p1p4*p2P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1P*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p4*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p2P*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p1P*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p3*p2P*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1P*p2p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*m32*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p4*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2p3*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p2*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p4*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p1P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*p1p2*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*p1p3*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*p1p4*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1P*p2p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1P*p2p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*p1P*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p2*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*M2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*p1P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p1p2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*M2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*M2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*m32*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1P*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p3*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p4*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1P*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p2P*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 40*m32*p1P*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p2P*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p3*p2P*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1P*p2p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 24*m32*m32*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p2*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 24*m32*p1p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p4*p2P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 24*m32*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p2*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p4*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*M2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*p1p2*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1P*p2p4*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p2*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2p3*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p4*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 48*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p2*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p2*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 24*m32*p1p3*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p3*p2p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p4*p3P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p3p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p4P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 2*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 2*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p2*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2p3*p2p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2p3*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p2p4*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p3P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*m32*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p2p3*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 20*m32*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2p3*p2p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2P*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2p4*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p3*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p3p4*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p3*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*m32*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) +
  (8*p1p3*p2p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*m32*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p1p3*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p1p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*m32*m32*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*p1p3*p2p4*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*p1p3*p3P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p2*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p3*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*p3P*p4P*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*m32*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 12*m32*m32*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p1p2*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 12*m32*p1p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 12*m32*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p4P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*p1P*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*m32*p1P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*m32*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p4*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1P*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*p1p4*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2P*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p3*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2p3*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*m32*m32*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p4*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p2p3*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1P*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p2*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p4*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p3*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p4*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2p3*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*p1P*p2p4*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*p1p2*p3P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 24*M2*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*p1P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*M2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*M2*m32*p1p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 24*M2*m32*p1p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p3*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p4*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2P*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*M2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 40*m32*p1P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2p3*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 48*m32*m32*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1P*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*p1p2*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 48*m32*p1p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 48*m32*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p2*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1p4*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*M2*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1P*p2p4*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*m32*p3P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1p2*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p4*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 48*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p2*p3P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p2*p3p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 32*m32*p1P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p2*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 24*m32*p1p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 32*p1P*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p1p3*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*p1p2*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 6*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 6*m32*p1p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 6*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p3*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p2P*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p3*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p4*p3P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2P*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p3*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*p2p4*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*m32*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p2P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p2p3*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 20*m32*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p2P*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p3*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p3P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p2p4*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p2p3*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p2p4*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*p2p4*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*m32*m32*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 64*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*p1p3*p2P*p2p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p2*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p3*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*m32*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*m32*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p1p3*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p1p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) -
  (8*m32*p2P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 16*m32*m32*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 12*m32*m32*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p1p2*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 12*m32*p1p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 12*m32*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p3p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p4P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (4*m32*m32)/((m32 + p2p3)*(m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 2*m32*p1p4)/((m32 + p2p3)*(m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 2*m32*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 2*p1p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 2*p1p2*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 64*m32*m32*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*m32*p1p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*p2p3*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*M2*m32*m32)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*M2*m32*p1p4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*M2*m32*p2p3)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*M2*p1p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*M2*p1p2*p3p4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*m32*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*m32*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*p1p2*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 64*m32*m32*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*m32*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*m32*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*m32*p1p4*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*m32*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p2*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*m32*p4P*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p3*p4P*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p2*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p2*p4P*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*m32*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 32*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1p4*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2P*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 32*m32*p1P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1P*p2p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 32*m32*m32*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p2*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p4*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 32*m32*p1P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p3*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p4*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1P*p2p3*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 32*p1P*p2p4*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p2*p3P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1P*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2P*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 40*m32*p1P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1P*p2p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 24*m32*m32*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1p2*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 24*m32*p1p4*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 24*m32*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p2*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p4*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1P*p2p4*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p3P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*m32*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1p2*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1p3*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2p3*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p2*p3p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*m32*p4P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1p2*p4P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 24*m32*p1p3*p4P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2p3*p4P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*p4P*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 12*m32*m32*p2P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*p1p4*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*p1P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 12*m32*m32*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*p1p4*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*p1p2*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p2*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*p1P*p2p3*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p2*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*m32*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p1p2*p2P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*p1p4*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*p1P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*m32*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p1p3*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*p1p2*p2p4*p3P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p2*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p3P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 8*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 12*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 4*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 24*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1P*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 24*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p2*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p3*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1P*p2p3*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p3P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p2*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1P*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1p3*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p2*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3P*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 24*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2P*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2p4*p3P*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2P*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*M2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2P*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p2P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2P*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p3P*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2P*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2P*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p3*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p3*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p3*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p2p4*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2P*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 20*m32*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p3*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p4*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p2p3*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*m32*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 5*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 8*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 12*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 64*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*p1p2*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*p1p4*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p2P*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 64*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*M2*m32*m32)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*M2*m32*p1p2)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*M2*m32*p1p3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*M2*m32*p1p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*m32*p2P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*p2P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p2*p2P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p3*p2P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*M2*m32*p2p3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2P*p2p3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*M2*m32*p2p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*M2*p1p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*p3P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*M2*m32*p3p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2P*p3p4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*m32*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2p3*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p3p4*p4P)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1p2*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p3*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*m32*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*m32*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p2*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p3*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2p3*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*p1p3*p2p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p3p4*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*m32*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p2*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p4*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2p3*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p3p4*qP*qP)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p2P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2P*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p3P*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2P*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p4P*qP*qp1)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 8*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p4*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p2*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1P*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1P*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 12*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*m32*p1p2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1P*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p4*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 12*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p2*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1P*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p3*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p3*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p3*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 24*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p3*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p4*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*m32*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*m32*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p1p2*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p1p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p1p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p2P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p2P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 12*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2p3*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2p3*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p4*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1P*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p4*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1P*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 2*m32*p1p2)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 2*m32*p1p3)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 2*m32*p1p4)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 2*m32*p2p3)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 2*m32*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 4*p1p3*p2p4)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 2*m32*p3p4)/((m32 + p2p3)*(m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*p2P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p2*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2p3*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*m32*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 4*m32*p1p2*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 4*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 4*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 4*m32*p2p3*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 4*m32*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p3P*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 4*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*m32*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p2*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p3P*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 32*p1p3*p4P*qP*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*M2*m32*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*m32)/((m32 + p2p3)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 256*m32*m32*p2P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*p1p4*p2P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*m32*p2P*p2p3*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*p1p2*p2P*p3p4*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*m32*m32)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*m32*p1p4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*m32*p2P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1p4*p2P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*m32*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*p1p4*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2P*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*p2p4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p4*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1P*p2p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p2P*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*m32*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1P*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*m32*p2P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p1p4*p2P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p2P*p2p3*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*m32*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p1p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p2p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*m32*qP)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 96*m32*m32*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p1p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p4*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*p2P*p3p4*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*p4P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p2p4*p3P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p2P*p3p4*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p2p3*p4P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*m32*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*m32*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*p1p4*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*p2p3*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*p1p3*p2p4*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*p3P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*p1P*p3p4*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*p1p2*p3p4*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2p3*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p4*p3P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1P*p3p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p4P*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*m32)/((m32 + p2p3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 256*m32*m32*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*m32*p1p4*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*p2p3*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*p1p3*p2p4*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*p1p2*p3P*p3p4*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*m32*m32)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*m32*p1p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*m32*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*p1p4*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*m32*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p4*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2p3*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1P*p2P*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3P*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*m32*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1P*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*m32*m32*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p1p4*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2p3*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p3*p2p4*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p2*p3P*p3p4*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*m32*m32*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p1p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*m32*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 96*m32*m32*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p1p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p4*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p3P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*p2p4*p3P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*p4P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p2p4*p3P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p2P*p3p4*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p2p3*p4P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p3P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p3*p3P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*p3P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*p1p3*p4P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p4*p3P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1P*p3p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*p3p4*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p3P*p4P*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 8*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p4*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p3*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1P*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*M2*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*m32*p1p3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*M2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*p3P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p1p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 16*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p2*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p3P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p2*p3p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 12*m32*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p1p2*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 12*m32*p1p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 12*m32*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p2*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p4*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p3p4*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 2*p1p2*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 3*m32*m32*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 3*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 3*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) + (
 2*m32*p3p4*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)) - (
 4*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p3*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p2*p3p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p2*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p1p3*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*p1p4*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p3p4*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p2p4*p3P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2P*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 8*m32*p2P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p3P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 12*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p2p3*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 16*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 32*m32*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p1P*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p2P*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p4P*qP*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) + (
 16*M2*m32*qp1*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)) - (
 4*m32*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1P*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1P*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p4*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1P*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p3*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 8*m32*p1P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*p2P*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*p1p2*p4P*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p4*p2P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1P*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p4P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*p2p4*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 64*p2P*p4P*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*m32*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*p1p4*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 128*m32*m32*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*m32*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*p1p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*m32*p2P*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*p2p3*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 16*p1P*p2p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 16*p1p2*p3p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p4*p2P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2p3*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1P*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p4P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p4P*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*m32*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*p1p4*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 4*m32*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*p1p4*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1P*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p1P*qp3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) + (
 4*p1p2*p4P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 4*m32*p4P*qp3*qp3)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)) - (
 32*p1p2*p4P*qp3*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 8*m32*p1p4)/((m32 + p2p3)*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 8*m32*p2p3)/((m32 + p2p3)*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 8*p1p3*p2p4)/((m32 + p2p3)*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 8*p1p2*p3p4)/((m32 + p2p3)*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p1p4*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p2p3*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p3P*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p2*p3p4*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*m32*qP*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p3P*qP*qp2)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1P*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 8*m32*p1p4)/((m32 + p2p3)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 8*m32*p2p3)/((m32 + p2p3)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 8*p1p3*p2p4)/((m32 + p2p3)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 8*p1p2*p3p4)/((m32 + p2p3)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p1p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2p3*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p3*p2p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p2*p3P*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1p2*p3p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*m32*qP*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*p2P*qP*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1P*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 16*m32*m32)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2p3)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*m32*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p2*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p4*p2p3*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3P*p3P*p3p4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*m32)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p2p4*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p2p4*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p3P*qp2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p4*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p4*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p2p3*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1P*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3P*p3p4*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p2p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2p4*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32)/((q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2)/((q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2p3)/((q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4)/((q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4)/((q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*m32*p4P*p4P)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p2*p4P*p4P)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p4*p2p3*p4P*p4P)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p4P*p4P)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3p4*p4P*p4P)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*m32)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p3*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*p4P)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p4P*qP)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p4P*qP)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p2p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p2p3*p4P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p3*p4P*qp1)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p2p3*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p4P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p4P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp3)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32)/((q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2)/((q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p3)/((q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4)/((q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4)/((q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (512*m32*m32*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*p1p2*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p1p4*p2p3*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p1p3*p2p4*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (256*m32*p3P*p3p4*p4P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*m32*p1p2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1P*p2P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p2*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2P*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2p3*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1P*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*m32*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1P*p2P*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3P*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p2*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2P*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1P*p2p3*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2p3*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p2*p3P*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3p4*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p2*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*m32*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p2*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p2p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*m32*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1p2*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1p4*p2p3*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p4*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p3*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p2*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p2p3*p3P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p2p4*p3P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p2p3*p4P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p2P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p2P*p3p4*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p4*p3P*qp2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1P*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*M2*m32*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*M2*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p3P*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*m32*m32*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p2*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p4*p2P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p2p3*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1P*p2p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p4P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*m32*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p2*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*p1p4*p2P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p4*p2p3*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p3p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p2*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1P*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p2*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*p2p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p2P*p4P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p2p3*p4P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1P*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p1P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*p1p3*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p2P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 80*m32*p1P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*m32*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*p1p2*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p4*p3P*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*p3P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p1P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p2P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p3*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p4*p3P*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p1P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*m32*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p1p2*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p1P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2p3*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*m32*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p1P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*p1p2*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p1P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*m32*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p2P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 40*m32*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p3*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p2P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p3*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p3*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 10*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p2p3*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p2P*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p2P*p2p3*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p2p4*p3P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2P*p3P*p3p4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*m32)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p1p2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p2*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p3P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 96*m32*m32*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p4*p3P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2P*p3p4*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p3*p4P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p1p2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p1p3*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p1p2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p3P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p2*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p3P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p3P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 6*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p2P*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p2p3*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p3P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1P*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2p3)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2p3*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p3P*p3P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p3p4*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p3P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p3*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p4P*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p3P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p3p4*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p3P*p4P*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p3P*qp2*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4)/((q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*m32*m32*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p2*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1p3*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p2p3*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p2p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (256*p1p3*p2p4*p3P*p3P)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p3P*p3p4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*m32*m32)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*M2*m32*p1p3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*M2*m32*p2p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*m32*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1P*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1p2*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p3*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p2P*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p2p3*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1P*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3P*p3p4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3P*p4P)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (256*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1p2*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*m32*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p3*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p2p4*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (256*p2p4*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p3P*qp2)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p3)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4)/((q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p3*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p4*p2P*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2P*p2p3*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*p2p4*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*p1p3*p2P*p2p4*p3P)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2P*p3P*p3p4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*p2P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p2p3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2P*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3p4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p3*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p2*p3P*p4P)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2p3*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*m32*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p2*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p2P*p2p4*p3P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp1)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2p4*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2P*p3p4*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2p3*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2p4*p3P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p3P*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p2P*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p4P*qP*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp1*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p3P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1p2*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p3*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p1p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p2p3*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (256*p1P*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p1p3*p2p4*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p3P*p3p4*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*m32*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1p2*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p3*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p2P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p2p3*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p3p4*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p4P*qp3)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1p2*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*m32*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*M2*p2p4*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p2p4*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*M2*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*m32*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2p3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1P*p2P*p2p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p2p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3p4*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp3)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p2*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1P*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*qp3*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p1P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p2*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1p3*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p2p3*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1P*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1p3*p2p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p3p4*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p1P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p2)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p3)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p4)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p2p3)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p2p4)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p3p4)/((m32 + p2p3)*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*p4P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*m32)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p1p2)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p1p3)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p1p4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p2p4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p3*p2p4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p1P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*p1p3*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p2P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 80*m32*p1P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 48*m32*m32*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*p1p2*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p4*p3P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*p1P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*p2P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p4*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p4*p3P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p2p3*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p3*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*p2P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2P*p2p3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2P*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2p3*p2p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2p4*p3P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p2P*p3p4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*p1P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p2P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2P*p2p3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2P*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2p3*p2p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2p4*p3P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p2P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3P*p3p4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2p3*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p3p4*p4P)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*m32*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2p3*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p3p4*p4P)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p2P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2P*p2p3*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p2p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p4*p3P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2P*p3p4*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p1P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p2P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p2p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p4*p3P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*p3p4*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p3p4*p4P*qP)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 48*m32*m32*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p3p4*p4P*qP)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p2p4*p3P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p2P*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*p2P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p2p4*p3P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p2P*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p2p3*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 40*m32*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p2p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3p4*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2P*p2p3*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p3*p2p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2P*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3p4*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p4P*p4P*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p4P*p4P*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*m32*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*p2p3*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 10*m32*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p2p3*p2p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p2p4*p3p4*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*m32*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p2p3*p2p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p2p3*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p2p4*p3p4*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*p2p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3p4*qP*qp1)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*p2p3*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p3*p2p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p3*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3p4*qP*qp1)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*m32*m32*p2P*p4P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p2*p2P*p4P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p4*p2P*p2p3*p4P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p2p4*p4P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*p3p4*p4P)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*m32)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p2)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p2P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*p2p3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*p2p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*p3P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p1p2*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2P*p3p4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3P*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p2P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*p2p3*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*p2p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3p4*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p4P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p4P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2p3*p4P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*p4P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*p4P*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 96*m32*m32*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3p4*qP*qP)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p4P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2P*p2p3*p4P*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p3*p4P*qp1)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p4*p3P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2P*p3p4*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p3*p4P*qP*qp1)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p4*p2p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*p3P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*M2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p2*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2P*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p1p2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p2P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2P*p3P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p2*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2P*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*p1p3*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p3P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p3p4*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p3*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p3P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p3p4*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p4P*p4P*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p4P*p4P*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*m32*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 6*m32*p1p3*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p1p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p3*p2p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 2*m32*p3p4*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*m32*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p4*p2p3*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p3*p2p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p3p4*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2p3*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p3p4*qP*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p3p4*qP*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p2p4*p3P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p2P*p3p4*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p2p4*p3P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p2P*p3p4*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p2p3*p4P*qp1*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*p2p3*p4P*qp1*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p3P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1P*p3p4*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p3P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1P*p3p4*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p4P*qp2*qp2)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p4*p2p3)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p3*p2p4)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4)/((m32 + p3p4)*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p3P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p4P*p4P*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p3*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p3P*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp2)/((m32 + p3p4)*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p3P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p3p4*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p4P*qP*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p3p4*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p3*p4P*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p3P*p4P*qp1*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2*qp2)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p2)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p3)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p4)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p2p3)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p2p4)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p3p4)/((m32 + p2p3)*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p3P*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p3p4*p4P)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p1p2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p1p3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p1p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p2P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p4*p2p3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p2p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*p3P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p1p2*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2P*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*p3p4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3P*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*p4P)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p3P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p3p4*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qP)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3p4*qP*qP)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2p4*p3P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*p4P*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p4*p3P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp1)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p4*p3P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2P*p3p4*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p3*p4P*qP*qp1)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p3P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p3P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*p4P*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p4P*qp2)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p3P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p3p4*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qP*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p3p4*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p3P*p4P*qp1*qp2)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p4)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p4)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4)/((q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*m32*m32*p2P*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p2*p2P*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p3*p2P*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p4*p2P*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*p2p3*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2P*p2p4*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*p1p3*p2P*p2p4*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*p3p4*p4P)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*m32)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p3)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p2P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p2P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p2P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p2P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p2p3)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p2p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p3*p2p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*p2p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p3P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3p4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p2P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2p3*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p2p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*m32*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp1)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp1)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp1)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p4P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp2)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p4P*qp2*qp2)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p3)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p3)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p4)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4)/((q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p1p3*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p4*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2p3*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2p4*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*p1p3*p2p4*p3P*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3P*p3p4*p4P)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p2)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p1p3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p1p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*p2p3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p1p4*p2p3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p2p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*p3P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*p3p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p2*p3p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2P*p3p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3p4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p3*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p2*p3P*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*p4P)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p3P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*p4P*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*m32*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2p3*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p2*p3p4*qP*qP)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p2p4*p3P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp1)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p2p4*p3P*qp1)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p2P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2p4*p3P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2P*p3p4*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2p3*p4P*qP*qp1)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p3P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp2)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp2)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p4P*qp2)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p3P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p3p4*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p4P*qP*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p3p4*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p3P*p4P*qp1*qp2)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p4*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 24*m32*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p1p2*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 24*m32*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p4*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p2*p4P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 6*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 6*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 6*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p4*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p1p2*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p2p3*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p4*p2P*p4P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2P*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp3)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p4P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p2p4*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2P*p4P*qp1*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1P*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p3P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p4*p2P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p4P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p2*p4P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp2*qp3)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p3P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*p4P*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*qp3)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p4P*qP*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp1*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p1p4*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p2P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp3)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp1*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1P*p2p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p4P*p4P*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp3)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p2*p4P*qP*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp1*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p4*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1P*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p4P*qp2*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p4P*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*qp3*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp3*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p4P*qp3*qp3)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp3*qp3)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p4*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*m32*p1P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p1P*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p1P*p2p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*m32*p1P*p2p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p1P*p3p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*p1P*p2p4*p3p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 16*m32*p1P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p1P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 8*p1P*p2p4*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p3*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p4*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*m32*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*m32*p2p3*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*p2p3*p2p3*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*m32*p2p4*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p3*p2p4*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 2*m32*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 2*p2p4*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p2p3*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 8*m32*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p3*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 4*m32*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 4*p2p4*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 32*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*p1p3*p2P*qp4)/((m32 + p3p4)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 32*p1p3*p2P*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) - (
 8*M2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 16*m32*m32*p2P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p2P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1p2*p2P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p4*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*p1P*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p3*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p2p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p3*p2P*p3P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1P*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p3P*p3P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p3P*p3P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p2*p3p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p3p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 12*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*M2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 24*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 12*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1P*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*M2*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p4P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1P*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p3P*p4P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3P*p4P*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*m32*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*m32*p3p4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*m32*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*m32*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 2*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) + (
 2*m32*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)) - (
 4*m32*m32*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1p2*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 8*p1p2*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p4*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2p4*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p3p4*qP*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p1p2*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p1p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 12*m32*p2P*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2p4*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2P*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p2p3*p4P*qp1*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 16*M2*p1p3*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 16*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)) + (
 8*m32*m32*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p1P*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1p2*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 8*p1P*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p2*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p3*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2p4*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p4*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1P*p3p4*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p1p2*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 2*m32*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1P*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*m32*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (
 4*p1p4*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p3*p4P*qp2*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*m32*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) + (
 4*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)) - (128*p1p3*p2P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p2*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p2*p3P*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p2*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p2*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1p3*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (32*p1p2*p3P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*m32*p1P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (128*p1P*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*M2*m32*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*M2*p2p3*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (256*p2P*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*M2*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1P*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p2*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p2P*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 128*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*m32*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 128*p1p3*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p2*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p2*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p1p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*m32*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*m32*p1p3*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*M2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*p1p4*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p2P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p2*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1P*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 6*m32*p1p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p2*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p2*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p3p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 24*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p2p4*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p2P*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p2p3*p4P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*m32*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*p1p3*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p2*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*M2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p1p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*m32*p1p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*M2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*M2*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p2*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*m32*p1p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1P*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 6*m32*p1p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 4*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p2*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p1p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 24*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*p2p4*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*m32*p4P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p2P*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p1P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p3*p2P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*m32*p3P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p1P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p2P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*p1p3*p2P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*m32*p1P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1P*p2p3*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*m32*p3P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*p1p2*p3P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*qp1*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*p2p3*qp1*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*p2P*p3P*qp1*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p4*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p2P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 16*m32*p3P*qp2*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 4*p1p3*p2P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p4*p2P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1P*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1P*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p3P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p3P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p4*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p4P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p4P*qp3*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*m32*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (64*M2*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 64*p1p2*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 64*M2*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (128*p1P*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (64*p1p3*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p4*p2P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1P*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p3P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p2*p4P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 32*M2*m32*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p1P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*M2*p1p2*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 32*p1p3*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*p3P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p2*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*m32*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p1p2*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 12*m32*p1p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p4*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p3p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p4P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*m32*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*m32*p1P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p2*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q\
2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p1p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p3*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p4*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*m32*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1P*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p2*p3p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*m32*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 16*m32*p1P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 4*m32*p1p2*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 12*m32*p1p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 4*m32*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) - (
 8*p1p3*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 8*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*p1P*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)) + (
 64*m32*m32*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*p1p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*m32*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 32*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*p1p2*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 32*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*M2*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*M2*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*M2*p1p2*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 16*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p2*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 32*M2*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p1p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 8*m32*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 16*m32*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 8*m32*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) + (
 16*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)) - (
 32*m32*p1P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*M2*p1p3*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1p4*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2p3*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p3P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3p4*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p1P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*M2*p1p2*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p1p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*m32*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p3*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*m32*p2p3*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p3*p2p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p3P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 32*p1p2*p3p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) + (
 32*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)) - (
 64*m32*m32*qp2*qp4)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*p1P*p3P*qp2*qp4)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp2)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 64*m32*m32*qp3*qp4)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) - (
 64*p1P*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p2p3 - q2 + 2*qp2 + 2*qp3)*(q2 - 2*qp4)) + (
 128*m32*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 256*m32*m32*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p1p2*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1P*p2p3*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p4*p2p3*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*p3P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3p4*p4P*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qp4)/(
 q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*qP*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*qP*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p4*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p2p3*qp1*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p2p3*p4P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p4P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p4P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 256*m32*m32*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p1p2*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p3*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1P*p2p3*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p4*p2p3*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2p4*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p3p4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*m32*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1p2*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*p1P*p2p3*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p4*p2p3*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*p1p3*p2p4*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (32*m32*p3p4*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p2*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1p3*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1P*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*p1p2*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*p2p3*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p2P*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p2p4*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*M2*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*p1P*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*m32*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (128*M2*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (256*p1P*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1p3*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1p4*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1P*p2p3*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p4*p2p3*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*p1P*p2p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*p1p3*p2p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (64*m32*p3p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p4P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*M2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*M2*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1P*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*M2*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p2*p3P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 4*m32*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*m32*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p1p2*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p2p4*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1P*p2p3*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p2p3*qp1*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2P*p3P*qp1*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1P*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*p1p2*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*m32*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p3P*qp4)/(
 q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (128*m32*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p1P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp4)/(
 q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p2*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p2p3*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p2P*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p2*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (64*m32*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (32*m32*p1P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1P*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p3*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p4*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2p4*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3p4*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*m32*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p1p3*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1P*p2p3*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p4*p2p3*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2p4*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp4)/((m32 + p3p4)*q2*(q2 - 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p2P*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p2*p3P*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp1*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p1p3*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2p3*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p3P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p3P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p4P*qp2*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p3p4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*m32*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*p1p3*p2P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*p1p3*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p2*p3P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p2p3*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*m32*m32*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p2*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p3*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p4*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2P*p2P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p2p3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p2p4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p2p4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p3p4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp4)/(
 q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*m32*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*M2*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p4*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p3*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p4*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2p4*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p4P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*m32*m32*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p2*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p1p3*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1p4*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1p3*p2P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p2p3*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2p4*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p1p3*p2p4*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*p3P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p3p4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*m32*p3P*p4P*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp4)/(
 q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p1P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p3*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p2*p3P*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qP*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*p2p3*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 128*p2P*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p2p4*p3P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp1*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp2*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*m32*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*M2*p1p2*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p4*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*m32*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*M2*p1p2*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p2*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p3*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p3*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3p4*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p4P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*M2*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p4P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*m32*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*M2*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p2*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p3*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 128*p1P*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2p3*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1P*p2p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3p4*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p4P*qp3*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 4*p1p3*p2P*qp4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1P*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p3P*qp4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p2*p3P*qp4*qp4)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 4*m32*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) - (
 4*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)) + (
 8*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1P*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) + (
 8*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 8*m32*p1P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)) - (
 16*m32*m32*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*m32*p1p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 8*m32*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*p1p3*p2p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 16*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) + (
 8*p1p2*p3p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)) - (
 64*m32*m32*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p1p2*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1P*p2p3*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p4*p2p3*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3p4*qp4*qp4)/(
 q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*m32*p3P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 8*m32*p3P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*m32*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p2*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p1p3*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1p4*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 64*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p2p3*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2p4*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*p1p3*p2p4*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 16*m32*p3p4*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp4)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*p1p3*p2P*qp4*qp4)/((m32 + p3p4)*q2*q2*(-q2 + 2*qp2)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p1P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p2P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 16*m32*p3P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p2P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 64*p1p3*p2P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p3P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p1P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 32*m32*p2P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) + (
 32*m32*p3P*qp4*qp4)/(
 q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(q2 - 2*qp4)*(2*m32 + 2*p3p4 + q2 - 2*qp3 - 2*qp4)) - (
 8*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*M2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*M2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 24*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1p4*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p1P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*m32*p1P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*M2*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1P*p2p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 24*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p3*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p2*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p4*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p2*p3P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*m32*p1P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*M2*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p3*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p4*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1P*p2p3*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*p1P*p2p4*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p2*p3P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*M2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 12*M2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1p4*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 20*m32*p1P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1P*p2p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 12*m32*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1p2*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 12*m32*p1p4*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 12*m32*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p2*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p4*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p3P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*M2*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1P*p2p4*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p3P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*m32*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*p1p2*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*p1p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*m32*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 16*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*p1p3*p2p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*p1p2*p3p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*m32*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1p2*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 24*m32*p1p3*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 16*p1p3*p2p4*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) + (
 2*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) + (
 2*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 m32*m32*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) + (
 m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 3*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) - (
 m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) + (
 2*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)) + (
 8*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p2p3*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p2p3*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p2p4*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 12*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p2p3*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 10*m32*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p2p3*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p2p4*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 32*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p1p2*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p1p4*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 32*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 32*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 32*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*m32*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 4*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 4*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 4*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 4*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 4*m32*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 4*m32*p3p4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 16*m32*p1P*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p1p2*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p1p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 32*p1p3*p2P*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*m32*p3P*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*M2*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1P*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p3*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p4*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p1p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1P*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p4*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 128*m32*m32*p2P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*m32*p1P*p2P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*m32*p1p4*p2P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*m32*p2P*p2p3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p1p2*p2P*p3p4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*m32*m32*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*m32*p1P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*m32*p1p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*m32*p2p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*p1p3*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*p1p2*p3P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*p1p2*p3p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*m32*m32*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*m32*p1P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*m32*p1p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1p3*p2P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*m32*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1P*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1p2*p3P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*M2*m32*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*M2*p2p3*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p2P*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 128*m32*m32*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*m32*p1P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*m32*p1p4*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p1p3*p2P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*m32*p2p3*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p1p3*p2p4*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 64*p1p2*p3P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p1p2*p3P*p3p4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*m32*m32*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*m32*p1P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*m32*p1p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*p1p3*p2P*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 16*m32*p2p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*p1p3*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 16*p1p2*p3p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*m32*m32*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*m32*p1P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*m32*p1p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p1p3*p2P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*m32*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p1P*p2p3*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*p1p2*p3P*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*M2*m32*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*M2*p2p3*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 64*p2P*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) - (
 32*m32*p3P*qp2*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 8*m32*m32*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p1p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1P*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p2*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p4*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p2*p3P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p2*p3p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 6*m32*m32*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 2*m32*p1p2*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 6*m32*p1p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p4*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 6*m32*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p2*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p4*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*m32*p3P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p3p4*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 8*m32*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) + (
 32*m32*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 4*m32*m32*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1p2*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 8*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p2*p3p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*m32*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 2*m32*p1p2*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 6*m32*p1p3*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 2*m32*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) + (
 4*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p3*p2p4*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 4*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)) - (
 8*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 16*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 8*m32*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)) - (
 32*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp2)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 32*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(2*m32 + 2*p2p3 + q2 - 2*qp2 - 2*qp3)*(-q2 + 2*qp4)) + (
 24*M2*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*p1P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 24*M2*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p1p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*M2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1P*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 24*m32*p1p3*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p4*p2P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q\
2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*M2*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 40*m32*p1P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1P*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p3*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 48*m32*m32*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1P*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 48*m32*p1p2*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p2p3*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 48*m32*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p2*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p4*p2p4*p3P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2P*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p4*p3p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*m32*p1P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*m32*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1p3*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1p4*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p4*p2P*p2P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1P*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p4*p2P*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p3*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 48*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*p1P*p2p3*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*m32*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1p2*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p4*p2P*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*p1p4*p2p3*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1P*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p2*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p4*p2p4*p3P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*p1p2*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p2P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p3*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p4*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p3P*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p2*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p3p4*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*m32*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p2*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p2*p2P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p2*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p3*p2p3*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1P*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p2*p3P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p3p4*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p4P*p4P*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p4P*p4P*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 6*m32*m32*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 6*m32*p1p2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 2*m32*p1p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p2p3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 6*m32*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p2*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p4*p2p4*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p1p2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*p1p4*p2p3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p2*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p4*p2p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p3p4*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p3*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p4*qP*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p2p3*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3p4*qP*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*m32*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p2p3*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 20*m32*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p2P*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p3*p2p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p4*p3p4*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p2p3*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p2P*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p2p3*p2p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p4*p3P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2P*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p3*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p4*p3p4*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p2p4*p4P*qp1*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p3*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p2p4*p4P*qp1*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*p1p2*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1p4*p2P*p2P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p4*p2P*p2p3*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1P*p2P*p2p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p3*p2P*p2p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p4*p2P*p3P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*p2P*p3p4*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*p2P*p4P*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*m32*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*p4P*qP*qp3)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*m32*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1P*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 12*m32*p1p3*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p1p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p3p4*qp2*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p3*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p2p3*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p3P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p3p4*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p4P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p4P*qp2*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 12*m32*m32*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 12*m32*p1p2*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p1p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p3*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 12*m32*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p4*qp3*qp3)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*m32*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p2*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2P*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p2p3*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3p4*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p4P*qp3*qp3)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32)/((m32 + p2p3)*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4)/((m32 + p2p3)*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*m32*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1P*p2p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p3*p2p4*p4P*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p2*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p3*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p4*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p4*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp3)/((m32 + p2p3)*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1P*p4P*qp2*qp3)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*m32*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p1p2*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p1p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*m32*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p2*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p4P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1P*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp3)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1P*p2P*qp3*qp4)/((m32 + p2p3)*q2*q2*(-q2 + 2*qp3)*(q2 - 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*M2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*p1P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*M2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*M2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p1p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p1P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 24*m32*p1p3*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p4*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p2P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*M2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*p1p4*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2P*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 40*m32*p1P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p3*p2P*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1P*p2p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 24*m32*m32*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 24*m32*p1p2*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p2*p2P*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p3*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 24*m32*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p2*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p4*p2p4*p3P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p2*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2P*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1P*p2p4*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*m32*p1P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*m32*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p3*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p4*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p1P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*p1p4*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 48*p1P*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p4*p2P*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p1P*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p3*p2P*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1P*p2p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*m32*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p1p2*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p2*p2P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*p1P*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p4*p2p3*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p2*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p4*p2p4*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*p3P*p3P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*M2*m32*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1P*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*M2*p1p2*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2P*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1P*p2p3*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1P*p2p4*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p3P*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*m32*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p2*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p1p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p3*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p2p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p3*p2p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p2*p3P*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p3p4*p4P*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*m32*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p1p2*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1p3*p2P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*m32*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p1P*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*p1p2*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*p1p3*p2p3*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 32*p1p3*p2p4*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p1p2*p3P*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 32*m32*p3p4*p4P*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 2*m32*m32*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 2*m32*p1p2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 2*m32*p1p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 2*m32*p2p3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p3*p2p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p3p4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*m32*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p1p2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*p1p2*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*p1p3*p2p3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*p1p3*p2p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p3p4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p3*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p3p4*qP*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p3*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p3p4*qP*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*m32*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p2p3*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 20*m32*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p3*p2p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p4*p3P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2P*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p4*p3p4*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p3*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p2P*p2p3*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*m32*p2p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 16*p2p3*p2p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*p2p3*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2p4*p3P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p2P*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p3*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p4*p3p4*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p3*p4P*qp1*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p2p3*p4P*qp1*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*p1p2*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p3*p2P*p2P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1P*p2P*p2p3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1p4*p2P*p2p3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1p3*p2P*p2p4*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*p2P*p3P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*p2P*p3p4*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1p3*p2P*p4P*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*m32*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*m32*p3P*qP*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*m32*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 12*m32*p1p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p1p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p2p3*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p3P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p3p4*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*m32*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p3*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p3P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p3p4*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p4P*qp2*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p4P*qp2*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*m32*m32*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1p3*p2P*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p1p3*p2p4*p3P*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p2*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp4)/((m32 + p2p3)*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 64*p2P*p3P*qp1*qp4)/((m32 + p2p3)*q2*q2*(q2 - 2*qp3)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 12*m32*m32*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 12*m32*p1p2*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 4*m32*p1p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p2p3*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 12*m32*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p2p4*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*m32*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p2*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p4*p2p3*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1P*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p4*p2p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p3P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p4P*qp3*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p4P*qp3*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 64*p1P*p2P*qp3*qp4)/((m32 + p3p4)*q2*q2*(q2 - 2*qp2)*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*m32*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p1P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p1p2*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p1p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 4*m32*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p2p3*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3P*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp4*qp4)/((m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 16*m32*m32*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p1p2*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*m32*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1P*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p2*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*p1p3*p2p3*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p3*p2p4*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) + (
 8*p1p2*p3P*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 8*m32*p3p4*qp4*qp4)/((m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2*(-q2 + 2*qp4)*(-2*m32 - 2*p3p4 - q2 + 2*qp3 + 2*qp4)) - (
 12*M2*m32*m32*m32)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 12*M2*m32*m32*m32)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 12*M2*m32*m32*m32)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) - (
 3*m32*m32*m32)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 3*m32*m32*m32)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) - (
 3*m32*m32*m32)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2) + (
 12*qP*m32*m32*m32)/((m32 + p2p3)*(m32 + p2p3)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 12*qP*m32*m32*m32)/((m32 + p3p4)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2) + (
 12*qP*m32*m32*m32)/((m32 + p2p3)*(m32 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*(m32 + p2p3 + p2p4 + p3p4)*q2*q2);  
         
      if(result < 0.0){result = 0.0;}       

    return result;
}


