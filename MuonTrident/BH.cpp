#include "functions.h"
#include "globais.h"

#include <iostream>
#include <cmath>



//*****************************************************************
// squared matrix elements for lepton tridents (Bethe-Heitler)
//*****************************************************************

double SquaredMatrixElementPLP_BH(double GFV, double GFA){

    double GFA2, GFV2, GFAGFV, result;
    double m32,m12;
    
    GFA2 = GFA*GFA;
    GFV2 = GFV*GFV;
    GFAGFV = GFA*GFV;
    m32 = m3*m3;
    m12 = m1*m1;
    
      result = (256*(2*M2*(qp4)*m12*(qp3)*(qp3)*(qp3)-
                8*(p1P)*(p2P)*(qp4)*(qp4)*(qp3)*(qp3)+
                4*M2*(p1p2)*(qp4)*(qp4)*(qp3)*(qp3)-
                M2*(p3p4)*q2*m12*(qp3)*(qp3)-
                2*M2*(qp4)*q2*m12*(qp3)*(qp3)+
                4*(p2P)*(Pq)*(p1p4)*(qp4)*(qp3)*(qp3)-
                2*M2*(qp1)*(p2p3)*(qp4)*(qp3)*(qp3)+
                4*(p1P)*(Pq)*(p2p4)*(qp4)*(qp3)*(qp3)-
                2*M2*(qp1)*(p2p4)*(qp4)*(qp3)*(qp3)-
                2*M2*(p1p3)*(qp2)*(qp4)*(qp3)*(qp3)-
                2*M2*(p1p4)*(qp2)*(qp4)*(qp3)*(qp3)-
                2*(p2P)*(Pq)*(p1p4)*q2*(qp3)*(qp3)+
                M2*(p1p4)*(p2p3)*q2*(qp3)*(qp3)-
                2*(p1P)*(Pq)*(p2p4)*q2*(qp3)*(qp3)+
                M2*(p1p3)*(p2p4)*q2*(qp3)*(qp3)+
                M2*(qp1)*(p2p4)*q2*(qp3)*(qp3)+
                M2*(p1p4)*(qp2)*q2*(qp3)*(qp3)+
                4*(p1P)*(p2P)*(qp4)*q2*(qp3)*(qp3)-
                2*M2*(p1p2)*(qp4)*q2*(qp3)*(qp3)+
                4*(p2P)*(Pq)*(p1p3)*(qp4)*(qp4)*(qp3)+
                4*(p1P)*(Pq)*(p2p3)*(qp4)*(qp4)*(qp3)-
                2*M2*(qp1)*(p2p3)*(qp4)*(qp4)*(qp3)-
                2*M2*(qp1)*(p2p4)*(qp4)*(qp4)*(qp3)-
                2*M2*(p1p3)*(qp2)*(qp4)*(qp4)*(qp3)-
                2*M2*(p1p4)*(qp2)*(qp4)*(qp4)*(qp3)+
                (p2P)*(Pq)*(p1p4)*q2*q2*(qp3)-
                2*M2*(p1p4)*(p2p3)*q2*q2*(qp3)+
                (p1P)*(Pq)*(p2p4)*q2*q2*(qp3)-
                2*M2*(p1p3)*(p2p4)*q2*q2*(qp3)-
                M2*(qp1)*(p2p4)*q2*q2*(qp3)-
                M2*(p1p4)*(qp2)*q2*q2*(qp3)-
                2*(p1P)*(p2P)*(p3p4)*q2*q2*(qp3)+
                2*M2*(p1p2)*(p3p4)*q2*q2*(qp3)-
                2*(p1P)*(p2P)*(qp4)*q2*q2*(qp3)+
                M2*(p1p2)*(qp4)*q2*q2*(qp3)+
                2*M2*(qp4)*(qp4)*(qp4)*m12*(qp3)+
                M2*(qp4)*q2*q2*m12*(qp3)-
                2*M2*(qp4)*(qp4)*q2*m12*(qp3)+
                2*M2*(p3p4)*(qp4)*q2*m12*(qp3)-
                4*(qP)*(qP)*(p1p4)*(p2p3)*(qp4)*(qp3)-
                4*(qP)*(qP)*(p1p3)*(p2p4)*(qp4)*(qp3)+
                4*(qP)*(qP)*(p1p2)*(p3p4)*(qp4)*(qp3)-
                4*(p2P)*(Pq)*(qp1)*(p3p4)*(qp4)*(qp3)-
                4*(p1P)*(Pq)*(qp2)*(p3p4)*(qp4)*(qp3)+
                4*M2*(qp1)*(qp2)*(p3p4)*(qp4)*(qp3)+
                4*(p1P)*(p2P)*(qp4)*(qp4)*q2*(qp3)-
                2*M2*(p1p2)*(qp4)*(qp4)*q2*(qp3)+
                2*(qP)*(qP)*(p1p4)*(p2p3)*q2*(qp3)+
                2*(qP)*(qP)*(p1p3)*(p2p4)*q2*(qp3)-
                2*(qP)*(qP)*(p1p2)*(p3p4)*q2*(qp3)+
                2*(p2P)*(Pq)*(qp1)*(p3p4)*q2*(qp3)+
                2*(p1P)*(Pq)*(qp2)*(p3p4)*q2*(qp3)-
                2*M2*(qp1)*(qp2)*(p3p4)*q2*(qp3)-
                2*(p2P)*(Pq)*(p1p3)*(qp4)*q2*(qp3)-
                2*(p2P)*(Pq)*(p1p4)*(qp4)*q2*(qp3)-
                2*(p1P)*(Pq)*(p2p3)*(qp4)*q2*(qp3)+
                2*M2*(p1p4)*(p2p3)*(qp4)*q2*(qp3)+
                3*M2*(qp1)*(p2p3)*(qp4)*q2*(qp3)-
                2*(p1P)*(Pq)*(p2p4)*(qp4)*q2*(qp3)+
                2*M2*(p1p3)*(p2p4)*(qp4)*q2*(qp3)+
                3*M2*(qp1)*(p2p4)*(qp4)*q2*(qp3)+
                3*M2*(p1p3)*(qp2)*(qp4)*q2*(qp3)+
                3*M2*(p1p4)*(qp2)*(qp4)*q2*(qp3)+
                4*(p1P)*(p2P)*(p3p4)*(qp4)*q2*(qp3)-
                4*M2*(p1p2)*(p3p4)*(qp4)*q2*(qp3)+
                M2*(p1p4)*(p2p3)*q2*q2*q2+
                M2*(p1p3)*(p2p4)*q2*q2*q2+
                (p1P)*(p2P)*(p3p4)*q2*q2*q2-
                M2*(p1p2)*(p3p4)*q2*q2*q2-
                (qP)*(qP)*(p1p4)*(p2p3)*q2*q2-
                (qP)*(qP)*(p1p3)*(p2p4)*q2*q2+
                (qP)*(qP)*(p1p2)*(p3p4)*q2*q2-
                (p2P)*(Pq)*(qp1)*(p3p4)*q2*q2-
                (p1P)*(Pq)*(qp2)*(p3p4)*q2*q2+
                M2*(qp1)*(qp2)*(p3p4)*q2*q2+
                (p2P)*(Pq)*(p1p3)*(qp4)*q2*q2+
                (p1P)*(Pq)*(p2p3)*(qp4)*q2*q2-
                2*M2*(p1p4)*(p2p3)*(qp4)*q2*q2-
                M2*(qp1)*(p2p3)*(qp4)*q2*q2-
                2*M2*(p1p3)*(p2p4)*(qp4)*q2*q2-
                M2*(p1p3)*(qp2)*(qp4)*q2*q2-
                2*(p1P)*(p2P)*(p3p4)*(qp4)*q2*q2+
                2*M2*(p1p2)*(p3p4)*(qp4)*q2*q2-
                M2*(p3p4)*(qp4)*(qp4)*q2*m12+
                (((p1p2)*(qp3-qp4)*(qp3-qp4)*q2-
                (qp1)*(qp2)*(2*(qp3)-q2)*(q2-2*(qp4)))*M2+
                (-q2*(q2*q2-2*(qp3+qp4)*q2+2*((qp3)*(qp3)+
                (qp4)*(qp4)))*M2-
                (qP)*(qP)*(2*(qp3)-q2)*(q2-2*(qp4)))*m12+
                (p1P)*(Pq)*(qp2)*(2*(qp3)-q2)*(q2-2*(qp4))-
                (p2P)*(2*(qp3)-q2)*(q2-2*(qp4))*((p1P)*q2-(Pq)*(qp1)))*m32-
                2*(p2P)*(Pq)*(p1p3)*(qp4)*(qp4)*q2-
                2*(p1P)*(Pq)*(p2p3)*(qp4)*(qp4)*q2+
                M2*(p1p4)*(p2p3)*(qp4)*(qp4)*q2+
                M2*(qp1)*(p2p3)*(qp4)*(qp4)*q2+
                M2*(p1p3)*(p2p4)*(qp4)*(qp4)*q2+
                M2*(p1p3)*(qp2)*(qp4)*(qp4)*q2+
                2*(qP)*(qP)*(p1p4)*(p2p3)*(qp4)*q2+
                2*(qP)*(qP)*(p1p3)*(p2p4)*(qp4)*q2-
                2*(qP)*(qP)*(p1p2)*(p3p4)*(qp4)*q2+
                2*(p2P)*(Pq)*(qp1)*(p3p4)*(qp4)*q2+
                2*(p1P)*(Pq)*(qp2)*(p3p4)*(qp4)*q2-
                2*M2*(qp1)*(qp2)*(p3p4)*(qp4)*q2+
                (p4P)*(p4P)*(2*(qp3)-q2)*(2*(qp3)*(qp3)*m12-
                2*(p3p4)*(qp3)*m12+2*(qp3)*(qp4)*m12+
                (p3p4)*q2*m12-2*(qp3)*q2*m12+
                (2*(qp3)-q2)*(p1p2-2*m12)*m32+
                2*(p1p3)*(p2p4)*(qp3)-2*(p1p3)*(qp2)*(qp3)+
                (p1p4)*(p2p3)*(2*(qp3)-q2)-(p1p3)*(p2p4)*q2+
                (p1p3)*(qp2)*q2+(qp1)*(p2p3)*(q2-2*(qp3)))+
                (p3P)*(p3P)*(2*(qp4)-q2)*(2*(qp4)*(qp4)*m12-
                2*(p3p4)*(qp4)*m12+2*(qp3)*(qp4)*m12+(p3p4)*q2*m12-
                2*(qp4)*q2*m12+(2*(qp4)-q2)*(p1p2-2*m12)*m32-2*(qp1)*(p2p4)*(qp4)+
                (p1p3)*(p2p4)*(2*(qp4)-q2)+(p1p4)*(p2p3-qp2)*(2*(qp4)-q2)+
                (qp1)*(p2p4)*q2)+(p3P)*(2*(qp4)-q2)*((Pq)*(-2*(qp4)*(qp4)*m12-
                2*(p3p4)*(qp3)*m12+2*(p3p4)*(qp4)*m12+(qp4)*q2*m12+
                2*(qp3-qp4)*(p1p2-2*m12)*m32+2*(p1p3)*(p2p4)*(qp3-qp4)+
                2*(qp1)*(p2p4)*(qp4)-2*(p1p2)*(qp3)*(qp4)+
                (p1p4)*(2*(p2p3)*(qp3-qp4)+(qp2)*(2*(qp4)-q2))-
                (qp1)*(p2p4)*q2+(p1p2)*(qp4)*q2)-
                (2*(qp3)-q2)*((p2P)*((p1p3-p1p4-qp1)*(qp4)+
                (p1p4)*q2)+(p1P)*((p2p3-p2p4-qp2)*(qp4)+
                (p2p4)*q2)))-(p4P)*(2*(qp3)-q2)*((Pq)*(2*(qp3)*(qp3)*m12-
                2*(p3p4)*(qp3)*m12+2*(p3p4)*(qp4)*m12-(qp3)*q2*m12+2*(qp3-
                qp4)*(p1p2-2*m12)*m32+2*(p1p3)*(p2p4)*(qp3)-
                2*(p1p3)*(qp2)*(qp3)+2*(p1p4)*(p2p3)*(qp3-qp4)-
                2*(p1p3)*(p2p4)*(qp4)+2*(p1p2)*(qp3)*(qp4)+
                (p1p3)*(qp2)*q2-(p1p2)*(qp3)*q2+
                (qp1)*(p2p3)*(q2-2*(qp3)))+(2*(qp4)-q2)*(-(((p2P)*(p1p3-p1p4+qp1)+
                (p1P)*(p2p3-p2p4+qp2))*(qp3))+
                ((p2P)*(p1p3)+(p1P)*(p2p3))*q2+(p3P)*(2*(-(p3p4)+qp3+qp4)*m12+
                2*(p1p2-2*m12)*m32-(qp1)*(p2p3+p2p4-2*(qp2))+(p1p4)*(2*(p2p3)-qp2)+
                (p1p3)*(2*(p2p4)-qp2)-q2*(m12 + p1p2))))))/(q2*q2*(q2-2*(qp3))*(q2-2*(qp3))*(q2-2*(qp4))*(q2-2*(qp4))*(2*m32+2*(p3p4)-
                2*(qp3)-2*(qp4)+q2)*(2*m32+2*(p3p4)-2*(qp3)-2*(qp4)+q2));
      
      /* (32*(4*M2*m12*m32*m32*m32*m32 - 8*M2*m12*m32*m32*m32*m32 - 4*M2*m32*m32*m32*m32*p1p2 + 
     12*M2*m32*m32*m32*m32*p1p2 - 8*m32*m32*m32*m32*p1P*p2P + 32*m32*m32*m32*m32*p1P*p2P - 
     48*m32*m32*m32*m32*p1P*p2P - 4*M2*m32*m32*m32*p1p4*p2p3 + 
     12*M2*m32*m32*m32*p1p4*p2p3 - 12*M2*m32*m32*m32*p1p4*p2p3 - 
     4*M2*m32*m32*m32*p1p3*p2p4 + 12*M2*m32*m32*m32*p1p3*p2p4 - 
     12*M2*m32*m32*m32*p1p3*p2p4 + 4*M2*m12*m32*m32*m32*p3p4 - 
     12*M2*m12*m32*m32*m32*p3p4 + 12*M2*m12*m32*m32*m32*p3p4 + 
     8*M2*m12*m32*m32*m32*q2 + m12*m32*m32*m32*m32*q2 - 16*M2*m12*m32*m32*m32*q2 - 
     2*m12*m32*m32*m32*m32*q2 + 8*M2*m12*m32*m32*m32*q2 - 
     8*M2*m32*m32*m32*p1p2*q2 - m32*m32*m32*m32*p1p2*q2 + 
     16*M2*m32*m32*m32*p1p2*q2 + 3*m32*m32*m32*m32*p1p2*q2 - 
     16*m32*m32*m32*p1P*p2P*q2 + 56*m32*m32*m32*p1P*p2P*q2 - 
     64*m32*m32*m32*p1P*p2P*q2 + 2*m32*m32*m32*p1p3*p2P*q2 - 
     6*m32*m32*m32*p1p3*p2P*q2 + 6*m32*m32*m32*p1p3*p2P*q2 + 
     2*m32*m32*m32*p1p4*p2P*q2 - 6*m32*m32*m32*p1p4*p2P*q2 + 
     6*m32*m32*m32*p1p4*p2P*q2 + 2*m32*m32*m32*p1P*p2p3*q2 - 
     6*m32*m32*m32*p1P*p2p3*q2 + 6*m32*m32*m32*p1P*p2p3*q2 - 
     m32*m32*m32*p1p4*p2p3*q2 + 3*m32*m32*m32*p1p4*p2p3*q2 - 
     3*m32*m32*m32*p1p4*p2p3*q2 + 2*m32*m32*m32*p1P*p2p4*q2 - 
     6*m32*m32*m32*p1P*p2p4*q2 + 6*m32*m32*m32*p1P*p2p4*q2 - 
     m32*m32*m32*p1p3*p2p4*q2 + 3*m32*m32*m32*p1p3*p2p4*q2 - 
     3*m32*m32*m32*p1p3*p2p4*q2 + 8*M2*m12*m32*m32*p3p4*q2 + 
     m12*m32*m32*m32*p3p4*q2 - 16*M2*m12*m32*m32*p3p4*q2 - 
     3*m12*m32*m32*m32*p3p4*q2 + 8*M2*m12*m32*m32*p3p4*q2 + 
     3*m12*m32*m32*m32*p3p4*q2 - 8*M2*m32*m32*p1p2*p3p4*q2 + 
     16*M2*m32*m32*p1p2*p3p4*q2 - 8*M2*m32*m32*p1p2*p3p4*q2 + 
     8*m32*m32*p1P*p2P*p3p4*q2 - 16*m32*m32*p1P*p2P*p3p4*q2 + 
     8*m32*m32*p1P*p2P*p3p4*q2 + 4*M2*m12*m32*m32*q2*q2 + 
     2*m12*m32*m32*m32*q2*q2 - 16*M2*m12*m32*m32*q2*q2 - 
     2*m12*m32*m32*m32*q2*q2 + 12*M2*m12*m32*m32*q2*q2 - 
     2*m12*m32*m32*m32*q2*q2 - 4*M2*m32*m32*p1p2*q2*q2 - 
     2*m32*m32*m32*p1p2*q2*q2 + 4*M2*m32*m32*p1p2*q2*q2 + 
     4*m32*m32*m32*p1p2*q2*q2 - 8*m32*m32*p1P*p2P*q2*q2 + 
     32*m32*m32*p1P*p2P*q2*q2 - 24*m32*m32*p1P*p2P*q2*q2 + 
     4*m32*m32*p1p3*p2P*q2*q2 - 8*m32*m32*p1p3*p2P*q2*q2 + 
     4*m32*m32*p1p3*p2P*q2*q2 + 4*m32*m32*p1p4*p2P*q2*q2 - 
     8*m32*m32*p1p4*p2P*q2*q2 + 4*m32*m32*p1p4*p2P*q2*q2 + 
     4*m32*m32*p1P*p2p3*q2*q2 - 8*m32*m32*p1P*p2p3*q2*q2 + 
     4*m32*m32*p1P*p2p3*q2*q2 + 12*M2*m32*p1p4*p2p3*q2*q2 - 
     2*m32*m32*p1p4*p2p3*q2*q2 - 12*M2*m32*p1p4*p2p3*q2*q2 + 
     4*m32*m32*p1p4*p2p3*q2*q2 - 2*m32*m32*p1p4*p2p3*q2*q2 + 
     4*m32*m32*p1P*p2p4*q2*q2 - 8*m32*m32*p1P*p2p4*q2*q2 + 
     4*m32*m32*p1P*p2p4*q2*q2 + 12*M2*m32*p1p3*p2p4*q2*q2 - 
     2*m32*m32*p1p3*p2p4*q2*q2 - 12*M2*m32*p1p3*p2p4*q2*q2 + 
     4*m32*m32*p1p3*p2p4*q2*q2 - 2*m32*m32*p1p3*p2p4*q2*q2 + 
     4*M2*m12*m32*p3p4*q2*q2 + 2*m12*m32*m32*p3p4*q2*q2 - 
     4*M2*m12*m32*p3p4*q2*q2 - 4*m12*m32*m32*p3p4*q2*q2 + 
     2*m12*m32*m32*p3p4*q2*q2 - 16*M2*m32*p1p2*p3p4*q2*q2 + 
     16*M2*m32*p1p2*p3p4*q2*q2 + 16*m32*p1P*p2P*p3p4*q2*q2 - 
     16*m32*p1P*p2P*p3p4*q2*q2 + m12*m32*m32*q2*q2*q2 - 8*M2*m12*m32*q2*q2*q2 - 
     m12*m32*m32*q2*q2*q2 - m32*m32*p1p2*q2*q2*q2 + m32*m32*p1p2*q2*q2*q2 + 
     8*m32*p1P*p2P*q2*q2*q2 + 2*m32*p1p3*p2P*q2*q2*q2 - 
     2*m32*p1p3*p2P*q2*q2*q2 + 2*m32*p1p4*p2P*q2*q2*q2 - 
     2*m32*p1p4*p2P*q2*q2*q2 + 2*m32*p1P*p2p3*q2*q2*q2 - 
     2*m32*p1P*p2p3*q2*q2*q2 + 8*M2*p1p4*p2p3*q2*q2*q2 - 
     m32*p1p4*p2p3*q2*q2*q2 + m32*p1p4*p2p3*q2*q2*q2 + 
     2*m32*p1P*p2p4*q2*q2*q2 - 2*m32*p1P*p2p4*q2*q2*q2 + 
     8*M2*p1p3*p2p4*q2*q2*q2 - m32*p1p3*p2p4*q2*q2*q2 + m32*p1p3*p2p4*q2*q2*q2 +
      m12*m32*p3p4*q2*q2*q2 - m12*m32*p3p4*q2*q2*q2 - 8*M2*p1p2*p3p4*q2*q2*q2 + 
     8*p1P*p2P*p3p4*q2*q2*q2 - 4*m12*m32*m32*m32*m32*qP + 8*m12*m32*m32*m32*m32*qP + 
     4*m32*m32*m32*m32*p1p2*qP - 12*m32*m32*m32*m32*p1p2*qP - 
     4*m32*m32*m32*p1p3*p2P*qP + 12*m32*m32*m32*p1p3*p2P*qP - 
     12*m32*m32*m32*p1p3*p2P*qP - 4*m32*m32*m32*p1p4*p2P*qP + 
     12*m32*m32*m32*p1p4*p2P*qP - 12*m32*m32*m32*p1p4*p2P*qP - 
     4*m32*m32*m32*p1P*p2p3*qP + 12*m32*m32*m32*p1P*p2p3*qP - 
     12*m32*m32*m32*p1P*p2p3*qP + 4*m32*m32*m32*p1p4*p2p3*qP - 
     12*m32*m32*m32*p1p4*p2p3*qP + 12*m32*m32*m32*p1p4*p2p3*qP - 
     4*m32*m32*m32*p1P*p2p4*qP + 12*m32*m32*m32*p1P*p2p4*qP - 
     12*m32*m32*m32*p1P*p2p4*qP + 4*m32*m32*m32*p1p3*p2p4*qP - 
     12*m32*m32*m32*p1p3*p2p4*qP + 12*m32*m32*m32*p1p3*p2p4*qP - 
     4*m12*m32*m32*m32*p3p4*qP + 12*m12*m32*m32*m32*p3p4*qP - 
     12*m12*m32*m32*m32*p3p4*qP - 8*m12*m32*m32*m32*q2*qP + 
     8*m12*m32*m32*m32*q2*qP + 8*m12*m32*m32*m32*q2*qP + 
     8*m32*m32*m32*p1p2*q2*qP - 16*m32*m32*m32*p1p2*q2*qP - 
     8*m32*m32*p1p3*p2P*q2*qP + 16*m32*m32*p1p3*p2P*q2*qP - 
     8*m32*m32*p1p3*p2P*q2*qP - 8*m32*m32*p1p4*p2P*q2*qP + 
     16*m32*m32*p1p4*p2P*q2*qP - 8*m32*m32*p1p4*p2P*q2*qP - 
     8*m32*m32*p1P*p2p3*q2*qP + 16*m32*m32*p1P*p2p3*q2*qP - 
     8*m32*m32*p1P*p2p3*q2*qP + 8*m32*m32*p1p4*p2p3*q2*qP - 
     16*m32*m32*p1p4*p2p3*q2*qP + 8*m32*m32*p1p4*p2p3*q2*qP - 
     8*m32*m32*p1P*p2p4*q2*qP + 16*m32*m32*p1P*p2p4*q2*qP - 
     8*m32*m32*p1P*p2p4*q2*qP + 8*m32*m32*p1p3*p2p4*q2*qP - 
     16*m32*m32*p1p3*p2p4*q2*qP + 8*m32*m32*p1p3*p2p4*q2*qP - 
     8*m12*m32*m32*p3p4*q2*qP + 16*m12*m32*m32*p3p4*q2*qP - 
     8*m12*m32*m32*p3p4*q2*qP - 4*m12*m32*m32*q2*q2*qP + 
     4*m12*m32*m32*q2*q2*qP + 4*m32*m32*p1p2*q2*q2*qP - 
     4*m32*m32*p1p2*q2*q2*qP - 4*m32*p1p3*p2P*q2*q2*qP + 
     4*m32*p1p3*p2P*q2*q2*qP - 4*m32*p1p4*p2P*q2*q2*qP + 
     4*m32*p1p4*p2P*q2*q2*qP - 4*m32*p1P*p2p3*q2*q2*qP + 
     4*m32*p1P*p2p3*q2*q2*qP + 4*m32*p1p4*p2p3*q2*q2*qP - 
     4*m32*p1p4*p2p3*q2*q2*qP - 4*m32*p1P*p2p4*q2*q2*qP + 
     4*m32*p1P*p2p4*q2*q2*qP + 4*m32*p1p3*p2p4*q2*q2*qP - 
     4*m32*p1p3*p2p4*q2*q2*qP - 4*m12*m32*p3p4*q2*q2*qP + 
     4*m12*m32*p3p4*q2*q2*qP + 8*m12*m32*m32*m32*qP*qP - 
     16*m12*m32*m32*m32*qP*qP - 8*m32*m32*p1p4*p2p3*qP*qP + 
     16*m32*m32*p1p4*p2p3*qP*qP - 8*m32*m32*p1p4*p2p3*qP*qP - 
     8*m32*m32*p1p3*p2p4*qP*qP + 16*m32*m32*p1p3*p2p4*qP*qP - 
     8*m32*m32*p1p3*p2p4*qP*qP + 8*m32*m32*p1p2*p3p4*qP*qP - 
     16*m32*m32*p1p2*p3p4*qP*qP + 8*m32*m32*p1p2*p3p4*qP*qP + 
     16*m12*m32*m32*q2*qP*qP - 16*m12*m32*m32*q2*qP*qP - 
     16*m32*p1p4*p2p3*q2*qP*qP + 16*m32*p1p4*p2p3*q2*qP*qP - 
     16*m32*p1p3*p2p4*q2*qP*qP + 16*m32*p1p3*p2p4*q2*qP*qP + 
     16*m32*p1p2*p3p4*q2*qP*qP - 16*m32*p1p2*p3p4*q2*qP*qP + 
     8*m12*m32*q2*q2*qP*qP - 8*p1p4*p2p3*q2*q2*qP*qP - 
     8*p1p3*p2p4*q2*q2*qP*qP + 8*p1p2*p3p4*q2*q2*qP*qP + 4*m32*m32*m32*m32*p2P*qp1 - 
     16*m32*m32*m32*m32*p2P*qp1 + 24*m32*m32*m32*m32*p2P*qp1 + 
     4*M2*m32*m32*m32*p2p3*qp1 - 12*M2*m32*m32*m32*p2p3*qp1 + 
     12*M2*m32*m32*m32*p2p3*qp1 + 4*M2*m32*m32*m32*p2p4*qp1 - 
     12*M2*m32*m32*m32*p2p4*qp1 + 12*M2*m32*m32*m32*p2p4*qp1 + 
     8*m32*m32*m32*p2P*q2*qp1 - 24*m32*m32*m32*p2P*q2*qp1 + 
     24*m32*m32*m32*p2P*q2*qp1 + 8*M2*m32*m32*p2p3*q2*qp1 - 
     16*M2*m32*m32*p2p3*q2*qp1 + 8*M2*m32*m32*p2p3*q2*qp1 + 
     8*M2*m32*m32*p2p4*q2*qp1 - 16*M2*m32*m32*p2p4*q2*qp1 + 
     8*M2*m32*m32*p2p4*q2*qp1 + 4*m32*m32*p2P*q2*q2*qp1 - 
     8*m32*m32*p2P*q2*q2*qp1 + 4*m32*m32*p2P*q2*q2*qp1 + 
     4*M2*m32*p2p3*q2*q2*qp1 - 4*M2*m32*p2p3*q2*q2*qp1 + 
     4*M2*m32*p2p4*q2*q2*qp1 - 4*M2*m32*p2p4*q2*q2*qp1 - 
     8*m32*m32*m32*p2P*qP*qp1 + 16*m32*m32*m32*p2P*qP*qp1 - 
     2*m32*m32*m32*p2p3*qP*qp1 + 6*m32*m32*m32*p2p3*qP*qp1 - 
     6*m32*m32*m32*p2p3*qP*qp1 - 2*m32*m32*m32*p2p4*qP*qp1 + 
     6*m32*m32*m32*p2p4*qP*qp1 - 6*m32*m32*m32*p2p4*qP*qp1 - 
     8*m32*m32*p2P*p3p4*qP*qp1 + 16*m32*m32*p2P*p3p4*qP*qp1 - 
     8*m32*m32*p2P*p3p4*qP*qp1 - 16*m32*m32*p2P*q2*qP*qp1 + 
     16*m32*m32*p2P*q2*qP*qp1 - 4*m32*m32*p2p3*q2*qP*qp1 + 
     8*m32*m32*p2p3*q2*qP*qp1 - 4*m32*m32*p2p3*q2*qP*qp1 - 
     4*m32*m32*p2p4*q2*qP*qp1 + 8*m32*m32*p2p4*q2*qP*qp1 - 
     4*m32*m32*p2p4*q2*qP*qp1 - 16*m32*p2P*p3p4*q2*qP*qp1 + 
     16*m32*p2P*p3p4*q2*qP*qp1 - 8*m32*p2P*q2*q2*qP*qp1 - 
     2*m32*p2p3*q2*q2*qP*qp1 + 2*m32*p2p3*q2*q2*qP*qp1 - 
     2*m32*p2p4*q2*q2*qP*qp1 + 2*m32*p2p4*q2*q2*qP*qp1 - 
     8*p2P*p3p4*q2*q2*qP*qp1 + 4*m32*m32*m32*m32*p1P*qp2 - 16*m32*m32*m32*m32*p1P*qp2 + 
     24*m32*m32*m32*m32*p1P*qp2 + 4*M2*m32*m32*m32*p1p3*qp2 - 
     12*M2*m32*m32*m32*p1p3*qp2 + 12*M2*m32*m32*m32*p1p3*qp2 + 
     4*M2*m32*m32*m32*p1p4*qp2 - 12*M2*m32*m32*m32*p1p4*qp2 + 
     12*M2*m32*m32*m32*p1p4*qp2 + 8*m32*m32*m32*p1P*q2*qp2 - 
     24*m32*m32*m32*p1P*q2*qp2 + 24*m32*m32*m32*p1P*q2*qp2 + 
     8*M2*m32*m32*p1p3*q2*qp2 - 16*M2*m32*m32*p1p3*q2*qp2 + 
     8*M2*m32*m32*p1p3*q2*qp2 + 8*M2*m32*m32*p1p4*q2*qp2 - 
     16*M2*m32*m32*p1p4*q2*qp2 + 8*M2*m32*m32*p1p4*q2*qp2 + 
     4*m32*m32*p1P*q2*q2*qp2 - 8*m32*m32*p1P*q2*q2*qp2 + 
     4*m32*m32*p1P*q2*q2*qp2 + 4*M2*m32*p1p3*q2*q2*qp2 - 
     4*M2*m32*p1p3*q2*q2*qp2 + 4*M2*m32*p1p4*q2*q2*qp2 - 
     4*M2*m32*p1p4*q2*q2*qp2 - 8*m32*m32*m32*p1P*qP*qp2 + 
     16*m32*m32*m32*p1P*qP*qp2 - 2*m32*m32*m32*p1p3*qP*qp2 + 
     6*m32*m32*m32*p1p3*qP*qp2 - 6*m32*m32*m32*p1p3*qP*qp2 - 
     2*m32*m32*m32*p1p4*qP*qp2 + 6*m32*m32*m32*p1p4*qP*qp2 - 
     6*m32*m32*m32*p1p4*qP*qp2 - 8*m32*m32*p1P*p3p4*qP*qp2 + 
     16*m32*m32*p1P*p3p4*qP*qp2 - 8*m32*m32*p1P*p3p4*qP*qp2 - 
     16*m32*m32*p1P*q2*qP*qp2 + 16*m32*m32*p1P*q2*qP*qp2 - 
     4*m32*m32*p1p3*q2*qP*qp2 + 8*m32*m32*p1p3*q2*qP*qp2 - 
     4*m32*m32*p1p3*q2*qP*qp2 - 4*m32*m32*p1p4*q2*qP*qp2 + 
     8*m32*m32*p1p4*q2*qP*qp2 - 4*m32*m32*p1p4*q2*qP*qp2 - 
     16*m32*p1P*p3p4*q2*qP*qp2 + 16*m32*p1P*p3p4*q2*qP*qp2 - 
     8*m32*p1P*q2*q2*qP*qp2 - 2*m32*p1p3*q2*q2*qP*qp2 + 
     2*m32*p1p3*q2*q2*qP*qp2 - 2*m32*p1p4*q2*q2*qP*qp2 + 
     2*m32*p1p4*q2*q2*qP*qp2 - 8*p1P*p3p4*q2*q2*qP*qp2 - 
     2*m32*m32*m32*m32*qp1*qp2 + 8*M2*m32*m32*m32*qp1*qp2 + 8*m32*m32*m32*m32*qp1*qp2 - 
     16*M2*m32*m32*m32*qp1*qp2 - 12*m32*m32*m32*m32*qp1*qp2 + 
     8*M2*m32*m32*p3p4*qp1*qp2 - 16*M2*m32*m32*p3p4*qp1*qp2 + 
     8*M2*m32*m32*p3p4*qp1*qp2 - 4*m32*m32*m32*q2*qp1*qp2 + 
     16*M2*m32*m32*q2*qp1*qp2 + 12*m32*m32*m32*q2*qp1*qp2 - 
     16*M2*m32*m32*q2*qp1*qp2 - 12*m32*m32*m32*q2*qp1*qp2 + 
     16*M2*m32*p3p4*q2*qp1*qp2 - 16*M2*m32*p3p4*q2*qp1*qp2 - 
     2*m32*m32*q2*q2*qp1*qp2 + 8*M2*m32*q2*q2*qp1*qp2 + 
     4*m32*m32*q2*q2*qp1*qp2 - 2*m32*m32*q2*q2*qp1*qp2 + 
     8*M2*p3p4*q2*q2*qp1*qp2 - 16*M2*m12*m32*m32*m32*qp3 + 
     32*M2*m12*m32*m32*m32*qp3 - 16*M2*m12*m32*m32*m32*qp3 + 
     8*M2*m32*m32*m32*p1p2*qp3 - 16*M2*m32*m32*m32*p1p2*qp3 + 
     24*m32*m32*m32*p1P*p2P*qp3 - 72*m32*m32*m32*p1P*p2P*qp3 + 
     72*m32*m32*m32*p1P*p2P*qp3 - 4*m32*m32*m32*p1p3*p2P*qp3 + 
     12*m32*m32*m32*p1p3*p2P*qp3 - 12*m32*m32*m32*p1p3*p2P*qp3 + 
     2*m32*m32*m32*p1p4*p2P*qp3 - 6*m32*m32*m32*p1p4*p2P*qp3 + 
     6*m32*m32*m32*p1p4*p2P*qp3 - 4*m32*m32*m32*p1P*p2p3*qp3 + 
     12*m32*m32*m32*p1P*p2p3*qp3 - 12*m32*m32*m32*p1P*p2p3*qp3 + 
     8*M2*m32*m32*p1p4*p2p3*qp3 - 16*M2*m32*m32*p1p4*p2p3*qp3 + 
     8*M2*m32*m32*p1p4*p2p3*qp3 + 2*m32*m32*m32*p1P*p2p4*qp3 - 
     6*m32*m32*m32*p1P*p2p4*qp3 + 6*m32*m32*m32*p1P*p2p4*qp3 + 
     8*M2*m32*m32*p1p3*p2p4*qp3 - 16*M2*m32*m32*p1p3*p2p4*qp3 + 
     8*M2*m32*m32*p1p3*p2p4*qp3 - 8*M2*m12*m32*m32*p3p4*qp3 + 
     16*M2*m12*m32*m32*p3p4*qp3 - 8*M2*m12*m32*m32*p3p4*qp3 - 
     24*M2*m12*m32*m32*q2*qp3 - 2*m12*m32*m32*m32*q2*qp3 + 
     48*M2*m12*m32*m32*q2*qp3 + 2*m12*m32*m32*m32*q2*qp3 - 
     24*M2*m12*m32*m32*q2*qp3 + 2*m12*m32*m32*m32*q2*qp3 + 
     8*M2*m32*m32*p1p2*q2*qp3 + 2*m32*m32*m32*p1p2*q2*qp3 - 
     8*M2*m32*m32*p1p2*q2*qp3 - 4*m32*m32*m32*p1p2*q2*qp3 + 
     32*m32*m32*p1P*p2P*q2*qp3 - 80*m32*m32*p1P*p2P*q2*qp3 + 
     48*m32*m32*p1P*p2P*q2*qp3 - 12*m32*m32*p1p3*p2P*q2*qp3 + 
     24*m32*m32*p1p3*p2P*q2*qp3 - 12*m32*m32*p1p3*p2P*q2*qp3 - 
     12*m32*m32*p1P*p2p3*q2*qp3 + 24*m32*m32*p1P*p2p3*q2*qp3 - 
     12*m32*m32*p1P*p2p3*q2*qp3 - 8*M2*m32*p1p4*p2p3*q2*qp3 + 
     2*m32*m32*p1p4*p2p3*q2*qp3 + 8*M2*m32*p1p4*p2p3*q2*qp3 - 
     4*m32*m32*p1p4*p2p3*q2*qp3 + 2*m32*m32*p1p4*p2p3*q2*qp3 - 
     8*M2*m32*p1p3*p2p4*q2*qp3 + 2*m32*m32*p1p3*p2p4*q2*qp3 + 
     8*M2*m32*p1p3*p2p4*q2*qp3 - 4*m32*m32*p1p3*p2p4*q2*qp3 + 
     2*m32*m32*p1p3*p2p4*q2*qp3 - 8*M2*m12*m32*p3p4*q2*qp3 - 
     2*m12*m32*m32*p3p4*q2*qp3 + 8*M2*m12*m32*p3p4*q2*qp3 + 
     4*m12*m32*m32*p3p4*q2*qp3 - 2*m12*m32*m32*p3p4*q2*qp3 + 
     16*M2*m32*p1p2*p3p4*q2*qp3 - 16*M2*m32*p1p2*p3p4*q2*qp3 - 
     16*m32*p1P*p2P*p3p4*q2*qp3 + 16*m32*p1P*p2P*p3p4*q2*qp3 - 
     8*M2*m12*m32*q2*q2*qp3 - 2*m12*m32*m32*q2*q2*qp3 + 
     24*M2*m12*m32*q2*q2*qp3 + 2*m12*m32*m32*q2*q2*qp3 + 
     2*m32*m32*p1p2*q2*q2*qp3 - 2*m32*m32*p1p2*q2*q2*qp3 + 
     8*m32*p1P*p2P*q2*q2*qp3 - 24*m32*p1P*p2P*q2*q2*qp3 - 
     8*m32*p1p3*p2P*q2*q2*qp3 + 8*m32*p1p3*p2P*q2*q2*qp3 - 
     2*m32*p1p4*p2P*q2*q2*qp3 + 2*m32*p1p4*p2P*q2*q2*qp3 - 
     8*m32*p1P*p2p3*q2*q2*qp3 + 8*m32*p1P*p2p3*q2*q2*qp3 - 
     16*M2*p1p4*p2p3*q2*q2*qp3 + 2*m32*p1p4*p2p3*q2*q2*qp3 - 
     2*m32*p1p4*p2p3*q2*q2*qp3 - 2*m32*p1P*p2p4*q2*q2*qp3 + 
     2*m32*p1P*p2p4*q2*q2*qp3 - 16*M2*p1p3*p2p4*q2*q2*qp3 + 
     2*m32*p1p3*p2p4*q2*q2*qp3 - 2*m32*p1p3*p2p4*q2*q2*qp3 - 
     2*m12*m32*p3p4*q2*q2*qp3 + 2*m12*m32*p3p4*q2*q2*qp3 + 
     16*M2*p1p2*p3p4*q2*q2*qp3 - 16*p1P*p2P*p3p4*q2*q2*qp3 + 
     12*m12*m32*m32*m32*qP*qp3 - 20*m12*m32*m32*m32*qP*qp3 + 
     4*m12*m32*m32*m32*qP*qp3 - 8*m32*m32*m32*p1p2*qP*qp3 + 
     16*m32*m32*m32*p1p2*qP*qp3 + 8*m32*m32*p1p3*p2P*qP*qp3 - 
     16*m32*m32*p1p3*p2P*qP*qp3 + 8*m32*m32*p1p3*p2P*qP*qp3 + 
     16*m32*m32*p1p4*p2P*qP*qp3 - 32*m32*m32*p1p4*p2P*qP*qp3 + 
     16*m32*m32*p1p4*p2P*qP*qp3 + 8*m32*m32*p1P*p2p3*qP*qp3 - 
     16*m32*m32*p1P*p2p3*qP*qp3 + 8*m32*m32*p1P*p2p3*qP*qp3 - 
     8*m32*m32*p1p4*p2p3*qP*qp3 + 16*m32*m32*p1p4*p2p3*qP*qp3 - 
     8*m32*m32*p1p4*p2p3*qP*qp3 + 16*m32*m32*p1P*p2p4*qP*qp3 - 
     32*m32*m32*p1P*p2p4*qP*qp3 + 16*m32*m32*p1P*p2p4*qP*qp3 - 
     8*m32*m32*p1p3*p2p4*qP*qp3 + 16*m32*m32*p1p3*p2p4*qP*qp3 - 
     8*m32*m32*p1p3*p2p4*qP*qp3 + 8*m12*m32*m32*p3p4*qP*qp3 - 
     16*m12*m32*m32*p3p4*qP*qp3 + 8*m12*m32*m32*p3p4*qP*qp3 + 
     16*m12*m32*m32*q2*qP*qp3 - 16*m12*m32*m32*q2*qP*qp3 - 
     8*m32*m32*p1p2*q2*qP*qp3 + 8*m32*m32*p1p2*q2*qP*qp3 + 
     8*m32*p1p3*p2P*q2*qP*qp3 - 8*m32*p1p3*p2P*q2*qP*qp3 + 
     24*m32*p1p4*p2P*q2*qP*qp3 - 24*m32*p1p4*p2P*q2*qP*qp3 + 
     8*m32*p1P*p2p3*q2*qP*qp3 - 8*m32*p1P*p2p3*q2*qP*qp3 - 
     8*m32*p1p4*p2p3*q2*qP*qp3 + 8*m32*p1p4*p2p3*q2*qP*qp3 + 
     24*m32*p1P*p2p4*q2*qP*qp3 - 24*m32*p1P*p2p4*q2*qP*qp3 - 
     8*m32*p1p3*p2p4*q2*qP*qp3 + 8*m32*p1p3*p2p4*q2*qP*qp3 + 
     8*m12*m32*p3p4*q2*qP*qp3 - 8*m12*m32*p3p4*q2*qP*qp3 + 
     4*m12*m32*q2*q2*qP*qp3 - 4*m12*m32*q2*q2*qP*qp3 + 
     8*p1p4*p2P*q2*q2*qP*qp3 + 8*p1P*p2p4*q2*q2*qP*qp3 - 
     16*m12*m32*m32*qP*qP*qp3 + 16*m12*m32*m32*qP*qP*qp3 + 
     16*m32*p1p4*p2p3*qP*qP*qp3 - 16*m32*p1p4*p2p3*qP*qP*qp3 + 
     16*m32*p1p3*p2p4*qP*qP*qp3 - 16*m32*p1p3*p2p4*qP*qP*qp3 - 
     16*m32*p1p2*p3p4*qP*qP*qp3 + 16*m32*p1p2*p3p4*qP*qP*qp3 - 
     16*m12*m32*q2*qP*qP*qp3 + 16*p1p4*p2p3*q2*qP*qP*qp3 + 
     16*p1p3*p2p4*q2*qP*qP*qp3 - 16*p1p2*p3p4*q2*qP*qP*qp3 - 
     10*m32*m32*m32*p2P*qp1*qp3 + 30*m32*m32*m32*p2P*qp1*qp3 - 
     30*m32*m32*m32*p2P*qp1*qp3 - 16*M2*m32*m32*p2p3*qp1*qp3 + 
     2*m32*m32*m32*p2p3*qp1*qp3 + 32*M2*m32*m32*p2p3*qp1*qp3 - 
     6*m32*m32*m32*p2p3*qp1*qp3 - 16*M2*m32*m32*p2p3*qp1*qp3 + 
     6*m32*m32*m32*p2p3*qp1*qp3 - 8*M2*m32*m32*p2p4*qp1*qp3 - 
     m32*m32*m32*p2p4*qp1*qp3 + 16*M2*m32*m32*p2p4*qp1*qp3 + 
     3*m32*m32*m32*p2p4*qp1*qp3 - 8*M2*m32*m32*p2p4*qp1*qp3 - 
     3*m32*m32*m32*p2p4*qp1*qp3 - 12*m32*m32*p2P*q2*qp1*qp3 + 
     24*m32*m32*p2P*q2*qp1*qp3 - 12*m32*m32*p2P*q2*qp1*qp3 - 
     16*M2*m32*p2p3*q2*qp1*qp3 + 2*m32*m32*p2p3*q2*qp1*qp3 + 
     16*M2*m32*p2p3*q2*qp1*qp3 - 4*m32*m32*p2p3*q2*qp1*qp3 + 
     2*m32*m32*p2p3*q2*qp1*qp3 - 16*M2*m32*p2p4*q2*qp1*qp3 + 
     16*M2*m32*p2p4*q2*qp1*qp3 - 2*m32*p2P*q2*q2*qp1*qp3 + 
     2*m32*p2P*q2*q2*qp1*qp3 - 8*M2*p2p4*q2*q2*qp1*qp3 + 
     m32*p2p4*q2*q2*qp1*qp3 - m32*p2p4*q2*q2*qp1*qp3 + 
     16*m32*m32*p2P*qP*qp1*qp3 - 16*m32*m32*p2P*qP*qp1*qp3 + 
     12*m32*m32*p2p3*qP*qp1*qp3 - 24*m32*m32*p2p3*qP*qp1*qp3 + 
     12*m32*m32*p2p3*qP*qp1*qp3 - 4*m32*m32*p2p4*qP*qp1*qp3 + 
     8*m32*m32*p2p4*qP*qp1*qp3 - 4*m32*m32*p2p4*qP*qp1*qp3 + 
     16*m32*p2P*p3p4*qP*qp1*qp3 - 16*m32*p2P*p3p4*qP*qp1*qp3 + 
     16*m32*p2P*q2*qP*qp1*qp3 + 12*m32*p2p3*q2*qP*qp1*qp3 - 
     12*m32*p2p3*q2*qP*qp1*qp3 - 4*m32*p2p4*q2*qP*qp1*qp3 + 
     4*m32*p2p4*q2*qP*qp1*qp3 + 16*p2P*p3p4*q2*qP*qp1*qp3 - 
     10*m32*m32*m32*p1P*qp2*qp3 + 30*m32*m32*m32*p1P*qp2*qp3 - 
     30*m32*m32*m32*p1P*qp2*qp3 - 16*M2*m32*m32*p1p3*qp2*qp3 + 
     2*m32*m32*m32*p1p3*qp2*qp3 + 32*M2*m32*m32*p1p3*qp2*qp3 - 
     6*m32*m32*m32*p1p3*qp2*qp3 - 16*M2*m32*m32*p1p3*qp2*qp3 + 
     6*m32*m32*m32*p1p3*qp2*qp3 - 8*M2*m32*m32*p1p4*qp2*qp3 - 
     m32*m32*m32*p1p4*qp2*qp3 + 16*M2*m32*m32*p1p4*qp2*qp3 + 
     3*m32*m32*m32*p1p4*qp2*qp3 - 8*M2*m32*m32*p1p4*qp2*qp3 - 
     3*m32*m32*m32*p1p4*qp2*qp3 - 12*m32*m32*p1P*q2*qp2*qp3 + 
     24*m32*m32*p1P*q2*qp2*qp3 - 12*m32*m32*p1P*q2*qp2*qp3 - 
     16*M2*m32*p1p3*q2*qp2*qp3 + 2*m32*m32*p1p3*q2*qp2*qp3 + 
     16*M2*m32*p1p3*q2*qp2*qp3 - 4*m32*m32*p1p3*q2*qp2*qp3 + 
     2*m32*m32*p1p3*q2*qp2*qp3 - 16*M2*m32*p1p4*q2*qp2*qp3 + 
     16*M2*m32*p1p4*q2*qp2*qp3 - 2*m32*p1P*q2*q2*qp2*qp3 + 
     2*m32*p1P*q2*q2*qp2*qp3 - 8*M2*p1p4*q2*q2*qp2*qp3 + 
     m32*p1p4*q2*q2*qp2*qp3 - m32*p1p4*q2*q2*qp2*qp3 + 
     16*m32*m32*p1P*qP*qp2*qp3 - 16*m32*m32*p1P*qP*qp2*qp3 + 
     12*m32*m32*p1p3*qP*qp2*qp3 - 24*m32*m32*p1p3*qP*qp2*qp3 + 
     12*m32*m32*p1p3*qP*qp2*qp3 - 4*m32*m32*p1p4*qP*qp2*qp3 + 
     8*m32*m32*p1p4*qP*qp2*qp3 - 4*m32*m32*p1p4*qP*qp2*qp3 + 
     16*m32*p1P*p3p4*qP*qp2*qp3 - 16*m32*p1P*p3p4*qP*qp2*qp3 + 
     16*m32*p1P*q2*qP*qp2*qp3 + 12*m32*p1p3*q2*qP*qp2*qp3 - 
     12*m32*p1p3*q2*qP*qp2*qp3 - 4*m32*p1p4*q2*qP*qp2*qp3 + 
     4*m32*p1p4*q2*qP*qp2*qp3 + 16*p1P*p3p4*q2*qP*qp2*qp3 + 
     4*m32*m32*m32*qp1*qp2*qp3 - 16*M2*m32*m32*qp1*qp2*qp3 - 
     12*m32*m32*m32*qp1*qp2*qp3 + 16*M2*m32*m32*qp1*qp2*qp3 + 
     12*m32*m32*m32*qp1*qp2*qp3 - 16*M2*m32*p3p4*qp1*qp2*qp3 + 
     16*M2*m32*p3p4*qp1*qp2*qp3 + 4*m32*m32*q2*qp1*qp2*qp3 - 
     16*M2*m32*q2*qp1*qp2*qp3 - 8*m32*m32*q2*qp1*qp2*qp3 + 
     4*m32*m32*q2*qp1*qp2*qp3 - 16*M2*p3p4*q2*qp1*qp2*qp3 + 
     24*M2*m12*m32*m32*qp3*qp3 - 2*m12*m32*m32*m32*qp3*qp3 - 
     32*M2*m12*m32*m32*qp3*qp3 + 2*m12*m32*m32*m32*qp3*qp3 + 
     8*M2*m12*m32*m32*qp3*qp3 + 2*m12*m32*m32*m32*qp3*qp3 - 
     8*M2*m32*m32*p1p2*qp3*qp3 + 2*m32*m32*m32*p1p2*qp3*qp3 + 
     8*M2*m32*m32*p1p2*qp3*qp3 - 4*m32*m32*m32*p1p2*qp3*qp3 - 
     16*m32*m32*p1P*p2P*qp3*qp3 + 32*m32*m32*p1P*p2P*qp3*qp3 - 
     16*m32*m32*p1P*p2P*qp3*qp3 + 8*m32*m32*p1p3*p2P*qp3*qp3 - 
     16*m32*m32*p1p3*p2P*qp3*qp3 + 8*m32*m32*p1p3*p2P*qp3*qp3 - 
     4*m32*m32*p1p4*p2P*qp3*qp3 + 8*m32*m32*p1p4*p2P*qp3*qp3 - 
     4*m32*m32*p1p4*p2P*qp3*qp3 + 8*m32*m32*p1P*p2p3*qp3*qp3 - 
     16*m32*m32*p1P*p2p3*qp3*qp3 + 8*m32*m32*p1P*p2p3*qp3*qp3 - 
     8*M2*m32*p1p4*p2p3*qp3*qp3 + 2*m32*m32*p1p4*p2p3*qp3*qp3 + 
     8*M2*m32*p1p4*p2p3*qp3*qp3 - 4*m32*m32*p1p4*p2p3*qp3*qp3 + 
     2*m32*m32*p1p4*p2p3*qp3*qp3 - 4*m32*m32*p1P*p2p4*qp3*qp3 + 
     8*m32*m32*p1P*p2p4*qp3*qp3 - 4*m32*m32*p1P*p2p4*qp3*qp3 - 
     8*M2*m32*p1p3*p2p4*qp3*qp3 + 2*m32*m32*p1p3*p2p4*qp3*qp3 + 
     8*M2*m32*p1p3*p2p4*qp3*qp3 - 4*m32*m32*p1p3*p2p4*qp3*qp3 + 
     2*m32*m32*p1p3*p2p4*qp3*qp3 + 8*M2*m12*m32*p3p4*qp3*qp3 - 
     2*m12*m32*m32*p3p4*qp3*qp3 - 8*M2*m12*m32*p3p4*qp3*qp3 + 
     4*m12*m32*m32*p3p4*qp3*qp3 - 2*m12*m32*m32*p3p4*qp3*qp3 + 
     24*M2*m12*m32*q2*qp3*qp3 - 2*m12*m32*m32*q2*qp3*qp3 - 
     40*M2*m12*m32*q2*qp3*qp3 + 8*m12*m32*m32*q2*qp3*qp3 - 
     6*m12*m32*m32*q2*qp3*qp3 + 8*M2*m32*p1p2*q2*qp3*qp3 - 
     2*m32*m32*p1p2*q2*qp3*qp3 + 2*m32*m32*p1p2*q2*qp3*qp3 - 
     16*m32*p1P*p2P*q2*qp3*qp3 + 16*m32*p1P*p2P*q2*qp3*qp3 + 
     8*m32*p1p3*p2P*q2*qp3*qp3 - 8*m32*p1p3*p2P*q2*qp3*qp3 - 
     4*m32*p1p4*p2P*q2*qp3*qp3 + 4*m32*p1p4*p2P*q2*qp3*qp3 + 
     8*m32*p1P*p2p3*q2*qp3*qp3 - 8*m32*p1P*p2p3*q2*qp3*qp3 + 
     8*M2*p1p4*p2p3*q2*qp3*qp3 - 2*m32*p1p4*p2p3*q2*qp3*qp3 + 
     2*m32*p1p4*p2p3*q2*qp3*qp3 - 4*m32*p1P*p2p4*q2*qp3*qp3 + 
     4*m32*p1P*p2p4*q2*qp3*qp3 + 8*M2*p1p3*p2p4*q2*qp3*qp3 - 
     2*m32*p1p3*p2p4*q2*qp3*qp3 + 2*m32*p1p3*p2p4*q2*qp3*qp3 - 
     8*M2*m12*p3p4*q2*qp3*qp3 + 2*m12*m32*p3p4*q2*qp3*qp3 - 
     2*m12*m32*p3p4*q2*qp3*qp3 - 16*m12*m32*m32*qP*qp3*qp3 + 
     16*m12*m32*m32*qP*qp3*qp3 + 16*m32*m32*p1p2*qP*qp3*qp3 - 
     16*m32*m32*p1p2*qP*qp3*qp3 - 16*m32*p1p4*p2P*qP*qp3*qp3 + 
     16*m32*p1p4*p2P*qP*qp3*qp3 + 16*m32*p1p4*p2p3*qP*qp3*qp3 - 
     16*m32*p1p4*p2p3*qP*qp3*qp3 - 16*m32*p1P*p2p4*qP*qp3*qp3 + 
     16*m32*p1P*p2p4*qP*qp3*qp3 + 16*m32*p1p3*p2p4*qP*qp3*qp3 - 
     16*m32*p1p3*p2p4*qP*qp3*qp3 - 16*m12*m32*p3p4*qP*qp3*qp3 + 
     16*m12*m32*p3p4*qP*qp3*qp3 - 16*m12*m32*q2*qP*qp3*qp3 + 
     16*m12*m32*q2*qP*qp3*qp3 - 16*p1p4*p2P*q2*qP*qp3*qp3 - 
     16*p1P*p2p4*q2*qP*qp3*qp3 + 4*m32*m32*p2P*qp1*qp3*qp3 - 
     8*m32*m32*p2P*qp1*qp3*qp3 + 4*m32*m32*p2P*qp1*qp3*qp3 + 
     16*M2*m32*p2p3*qp1*qp3*qp3 - 4*m32*m32*p2p3*qp1*qp3*qp3 - 
     16*M2*m32*p2p3*qp1*qp3*qp3 + 8*m32*m32*p2p3*qp1*qp3*qp3 - 
     4*m32*m32*p2p3*qp1*qp3*qp3 + 8*M2*m32*p2p4*qp1*qp3*qp3 - 
     8*M2*m32*p2p4*qp1*qp3*qp3 + 4*m32*p2P*q2*qp1*qp3*qp3 - 
     4*m32*p2P*q2*qp1*qp3*qp3 + 8*M2*p2p4*q2*qp1*qp3*qp3 - 
     16*m32*p2p3*qP*qp1*qp3*qp3 + 16*m32*p2p3*qP*qp1*qp3*qp3 + 
     4*m32*m32*p1P*qp2*qp3*qp3 - 8*m32*m32*p1P*qp2*qp3*qp3 + 
     4*m32*m32*p1P*qp2*qp3*qp3 + 16*M2*m32*p1p3*qp2*qp3*qp3 - 
     4*m32*m32*p1p3*qp2*qp3*qp3 - 16*M2*m32*p1p3*qp2*qp3*qp3 + 
     8*m32*m32*p1p3*qp2*qp3*qp3 - 4*m32*m32*p1p3*qp2*qp3*qp3 + 
     8*M2*m32*p1p4*qp2*qp3*qp3 - 8*M2*m32*p1p4*qp2*qp3*qp3 + 
     4*m32*p1P*q2*qp2*qp3*qp3 - 4*m32*p1P*q2*qp2*qp3*qp3 + 
     8*M2*p1p4*q2*qp2*qp3*qp3 - 16*m32*p1p3*qP*qp2*qp3*qp3 + 
     16*m32*p1p3*qP*qp2*qp3*qp3 - 16*M2*m12*m32*qp3*qp3*qp3 + 
     4*m12*m32*m32*qp3*qp3*qp3 + 16*M2*m12*m32*qp3*qp3*qp3 - 8*m12*m32*m32*qp3*qp3*qp3 + 
     4*m12*m32*m32*qp3*qp3*qp3 + 16*m12*m32*qP*qp3*qp3*qp3 - 16*m12*m32*qP*qp3*qp3*qp3 - 
     16*M2*m12*m32*m32*m32*qp4 + 32*M2*m12*m32*m32*m32*qp4 - 
     16*M2*m12*m32*m32*m32*qp4 + 8*M2*m32*m32*m32*p1p2*qp4 - 
     16*M2*m32*m32*m32*p1p2*qp4 + 24*m32*m32*m32*p1P*p2P*qp4 - 
     72*m32*m32*m32*p1P*p2P*qp4 + 72*m32*m32*m32*p1P*p2P*qp4 + 
     2*m32*m32*m32*p1p3*p2P*qp4 - 6*m32*m32*m32*p1p3*p2P*qp4 + 
     6*m32*m32*m32*p1p3*p2P*qp4 - 4*m32*m32*m32*p1p4*p2P*qp4 + 
     12*m32*m32*m32*p1p4*p2P*qp4 - 12*m32*m32*m32*p1p4*p2P*qp4 + 
     2*m32*m32*m32*p1P*p2p3*qp4 - 6*m32*m32*m32*p1P*p2p3*qp4 + 
     6*m32*m32*m32*p1P*p2p3*qp4 + 8*M2*m32*m32*p1p4*p2p3*qp4 - 
     16*M2*m32*m32*p1p4*p2p3*qp4 + 8*M2*m32*m32*p1p4*p2p3*qp4 - 
     4*m32*m32*m32*p1P*p2p4*qp4 + 12*m32*m32*m32*p1P*p2p4*qp4 - 
     12*m32*m32*m32*p1P*p2p4*qp4 + 8*M2*m32*m32*p1p3*p2p4*qp4 - 
     16*M2*m32*m32*p1p3*p2p4*qp4 + 8*M2*m32*m32*p1p3*p2p4*qp4 - 
     8*M2*m12*m32*m32*p3p4*qp4 + 16*M2*m12*m32*m32*p3p4*qp4 - 
     8*M2*m12*m32*m32*p3p4*qp4 - 24*M2*m12*m32*m32*q2*qp4 - 
     2*m12*m32*m32*m32*q2*qp4 + 48*M2*m12*m32*m32*q2*qp4 + 
     2*m12*m32*m32*m32*q2*qp4 - 24*M2*m12*m32*m32*q2*qp4 + 
     2*m12*m32*m32*m32*q2*qp4 + 8*M2*m32*m32*p1p2*q2*qp4 + 
     2*m32*m32*m32*p1p2*q2*qp4 - 8*M2*m32*m32*p1p2*q2*qp4 - 
     4*m32*m32*m32*p1p2*q2*qp4 + 32*m32*m32*p1P*p2P*q2*qp4 - 
     80*m32*m32*p1P*p2P*q2*qp4 + 48*m32*m32*p1P*p2P*q2*qp4 - 
     12*m32*m32*p1p4*p2P*q2*qp4 + 24*m32*m32*p1p4*p2P*q2*qp4 - 
     12*m32*m32*p1p4*p2P*q2*qp4 - 8*M2*m32*p1p4*p2p3*q2*qp4 + 
     2*m32*m32*p1p4*p2p3*q2*qp4 + 8*M2*m32*p1p4*p2p3*q2*qp4 - 
     4*m32*m32*p1p4*p2p3*q2*qp4 + 2*m32*m32*p1p4*p2p3*q2*qp4 - 
     12*m32*m32*p1P*p2p4*q2*qp4 + 24*m32*m32*p1P*p2p4*q2*qp4 - 
     12*m32*m32*p1P*p2p4*q2*qp4 - 8*M2*m32*p1p3*p2p4*q2*qp4 + 
     2*m32*m32*p1p3*p2p4*q2*qp4 + 8*M2*m32*p1p3*p2p4*q2*qp4 - 
     4*m32*m32*p1p3*p2p4*q2*qp4 + 2*m32*m32*p1p3*p2p4*q2*qp4 - 
     8*M2*m12*m32*p3p4*q2*qp4 - 2*m12*m32*m32*p3p4*q2*qp4 + 
     8*M2*m12*m32*p3p4*q2*qp4 + 4*m12*m32*m32*p3p4*q2*qp4 - 
     2*m12*m32*m32*p3p4*q2*qp4 + 16*M2*m32*p1p2*p3p4*q2*qp4 - 
     16*M2*m32*p1p2*p3p4*q2*qp4 - 16*m32*p1P*p2P*p3p4*q2*qp4 + 
     16*m32*p1P*p2P*p3p4*q2*qp4 - 8*M2*m12*m32*q2*q2*qp4 - 
     2*m12*m32*m32*q2*q2*qp4 + 24*M2*m12*m32*q2*q2*qp4 + 
     2*m12*m32*m32*q2*q2*qp4 + 2*m32*m32*p1p2*q2*q2*qp4 - 
     2*m32*m32*p1p2*q2*q2*qp4 + 8*m32*p1P*p2P*q2*q2*qp4 - 
     24*m32*p1P*p2P*q2*q2*qp4 - 2*m32*p1p3*p2P*q2*q2*qp4 + 
     2*m32*p1p3*p2P*q2*q2*qp4 - 8*m32*p1p4*p2P*q2*q2*qp4 + 
     8*m32*p1p4*p2P*q2*q2*qp4 - 2*m32*p1P*p2p3*q2*q2*qp4 + 
     2*m32*p1P*p2p3*q2*q2*qp4 - 16*M2*p1p4*p2p3*q2*q2*qp4 + 
     2*m32*p1p4*p2p3*q2*q2*qp4 - 2*m32*p1p4*p2p3*q2*q2*qp4 - 
     8*m32*p1P*p2p4*q2*q2*qp4 + 8*m32*p1P*p2p4*q2*q2*qp4 - 
     16*M2*p1p3*p2p4*q2*q2*qp4 + 2*m32*p1p3*p2p4*q2*q2*qp4 - 
     2*m32*p1p3*p2p4*q2*q2*qp4 - 2*m12*m32*p3p4*q2*q2*qp4 + 
     2*m12*m32*p3p4*q2*q2*qp4 + 16*M2*p1p2*p3p4*q2*q2*qp4 - 
     16*p1P*p2P*p3p4*q2*q2*qp4 + 12*m12*m32*m32*m32*qP*qp4 - 
     20*m12*m32*m32*m32*qP*qp4 + 4*m12*m32*m32*m32*qP*qp4 - 
     8*m32*m32*m32*p1p2*qP*qp4 + 16*m32*m32*m32*p1p2*qP*qp4 + 
     16*m32*m32*p1p3*p2P*qP*qp4 - 32*m32*m32*p1p3*p2P*qP*qp4 + 
     16*m32*m32*p1p3*p2P*qP*qp4 + 8*m32*m32*p1p4*p2P*qP*qp4 - 
     16*m32*m32*p1p4*p2P*qP*qp4 + 8*m32*m32*p1p4*p2P*qP*qp4 + 
     16*m32*m32*p1P*p2p3*qP*qp4 - 32*m32*m32*p1P*p2p3*qP*qp4 + 
     16*m32*m32*p1P*p2p3*qP*qp4 - 8*m32*m32*p1p4*p2p3*qP*qp4 + 
     16*m32*m32*p1p4*p2p3*qP*qp4 - 8*m32*m32*p1p4*p2p3*qP*qp4 + 
     8*m32*m32*p1P*p2p4*qP*qp4 - 16*m32*m32*p1P*p2p4*qP*qp4 + 
     8*m32*m32*p1P*p2p4*qP*qp4 - 8*m32*m32*p1p3*p2p4*qP*qp4 + 
     16*m32*m32*p1p3*p2p4*qP*qp4 - 8*m32*m32*p1p3*p2p4*qP*qp4 + 
     8*m12*m32*m32*p3p4*qP*qp4 - 16*m12*m32*m32*p3p4*qP*qp4 + 
     8*m12*m32*m32*p3p4*qP*qp4 + 16*m12*m32*m32*q2*qP*qp4 - 
     16*m12*m32*m32*q2*qP*qp4 - 8*m32*m32*p1p2*q2*qP*qp4 + 
     8*m32*m32*p1p2*q2*qP*qp4 + 24*m32*p1p3*p2P*q2*qP*qp4 - 
     24*m32*p1p3*p2P*q2*qP*qp4 + 8*m32*p1p4*p2P*q2*qP*qp4 - 
     8*m32*p1p4*p2P*q2*qP*qp4 + 24*m32*p1P*p2p3*q2*qP*qp4 - 
     24*m32*p1P*p2p3*q2*qP*qp4 - 8*m32*p1p4*p2p3*q2*qP*qp4 + 
     8*m32*p1p4*p2p3*q2*qP*qp4 + 8*m32*p1P*p2p4*q2*qP*qp4 - 
     8*m32*p1P*p2p4*q2*qP*qp4 - 8*m32*p1p3*p2p4*q2*qP*qp4 + 
     8*m32*p1p3*p2p4*q2*qP*qp4 + 8*m12*m32*p3p4*q2*qP*qp4 - 
     8*m12*m32*p3p4*q2*qP*qp4 + 4*m12*m32*q2*q2*qP*qp4 - 
     4*m12*m32*q2*q2*qP*qp4 + 8*p1p3*p2P*q2*q2*qP*qp4 + 
     8*p1P*p2p3*q2*q2*qP*qp4 - 16*m12*m32*m32*qP*qP*qp4 + 
     16*m12*m32*m32*qP*qP*qp4 + 16*m32*p1p4*p2p3*qP*qP*qp4 - 
     16*m32*p1p4*p2p3*qP*qP*qp4 + 16*m32*p1p3*p2p4*qP*qP*qp4 - 
     16*m32*p1p3*p2p4*qP*qP*qp4 - 16*m32*p1p2*p3p4*qP*qP*qp4 + 
     16*m32*p1p2*p3p4*qP*qP*qp4 - 16*m12*m32*q2*qP*qP*qp4 + 
     16*p1p4*p2p3*q2*qP*qP*qp4 + 16*p1p3*p2p4*q2*qP*qP*qp4 - 
     16*p1p2*p3p4*q2*qP*qP*qp4 - 10*m32*m32*m32*p2P*qp1*qp4 + 
     30*m32*m32*m32*p2P*qp1*qp4 - 30*m32*m32*m32*p2P*qp1*qp4 - 
     8*M2*m32*m32*p2p3*qp1*qp4 - m32*m32*m32*p2p3*qp1*qp4 + 
     16*M2*m32*m32*p2p3*qp1*qp4 + 3*m32*m32*m32*p2p3*qp1*qp4 - 
     8*M2*m32*m32*p2p3*qp1*qp4 - 3*m32*m32*m32*p2p3*qp1*qp4 - 
     16*M2*m32*m32*p2p4*qp1*qp4 + 2*m32*m32*m32*p2p4*qp1*qp4 + 
     32*M2*m32*m32*p2p4*qp1*qp4 - 6*m32*m32*m32*p2p4*qp1*qp4 - 
     16*M2*m32*m32*p2p4*qp1*qp4 + 6*m32*m32*m32*p2p4*qp1*qp4 - 
     12*m32*m32*p2P*q2*qp1*qp4 + 24*m32*m32*p2P*q2*qp1*qp4 - 
     12*m32*m32*p2P*q2*qp1*qp4 - 16*M2*m32*p2p3*q2*qp1*qp4 + 
     16*M2*m32*p2p3*q2*qp1*qp4 - 16*M2*m32*p2p4*q2*qp1*qp4 + 
     2*m32*m32*p2p4*q2*qp1*qp4 + 16*M2*m32*p2p4*q2*qp1*qp4 - 
     4*m32*m32*p2p4*q2*qp1*qp4 + 2*m32*m32*p2p4*q2*qp1*qp4 - 
     2*m32*p2P*q2*q2*qp1*qp4 + 2*m32*p2P*q2*q2*qp1*qp4 - 
     8*M2*p2p3*q2*q2*qp1*qp4 + m32*p2p3*q2*q2*qp1*qp4 - 
     m32*p2p3*q2*q2*qp1*qp4 + 16*m32*m32*p2P*qP*qp1*qp4 - 
     16*m32*m32*p2P*qP*qp1*qp4 - 4*m32*m32*p2p3*qP*qp1*qp4 + 
     8*m32*m32*p2p3*qP*qp1*qp4 - 4*m32*m32*p2p3*qP*qp1*qp4 + 
     12*m32*m32*p2p4*qP*qp1*qp4 - 24*m32*m32*p2p4*qP*qp1*qp4 + 
     12*m32*m32*p2p4*qP*qp1*qp4 + 16*m32*p2P*p3p4*qP*qp1*qp4 - 
     16*m32*p2P*p3p4*qP*qp1*qp4 + 16*m32*p2P*q2*qP*qp1*qp4 - 
     4*m32*p2p3*q2*qP*qp1*qp4 + 4*m32*p2p3*q2*qP*qp1*qp4 + 
     12*m32*p2p4*q2*qP*qp1*qp4 - 12*m32*p2p4*q2*qP*qp1*qp4 + 
     16*p2P*p3p4*q2*qP*qp1*qp4 - 10*m32*m32*m32*p1P*qp2*qp4 + 
     30*m32*m32*m32*p1P*qp2*qp4 - 30*m32*m32*m32*p1P*qp2*qp4 - 
     8*M2*m32*m32*p1p3*qp2*qp4 - m32*m32*m32*p1p3*qp2*qp4 + 
     16*M2*m32*m32*p1p3*qp2*qp4 + 3*m32*m32*m32*p1p3*qp2*qp4 - 
     8*M2*m32*m32*p1p3*qp2*qp4 - 3*m32*m32*m32*p1p3*qp2*qp4 - 
     16*M2*m32*m32*p1p4*qp2*qp4 + 2*m32*m32*m32*p1p4*qp2*qp4 + 
     32*M2*m32*m32*p1p4*qp2*qp4 - 6*m32*m32*m32*p1p4*qp2*qp4 - 
     16*M2*m32*m32*p1p4*qp2*qp4 + 6*m32*m32*m32*p1p4*qp2*qp4 - 
     12*m32*m32*p1P*q2*qp2*qp4 + 24*m32*m32*p1P*q2*qp2*qp4 - 
     12*m32*m32*p1P*q2*qp2*qp4 - 16*M2*m32*p1p3*q2*qp2*qp4 + 
     16*M2*m32*p1p3*q2*qp2*qp4 - 16*M2*m32*p1p4*q2*qp2*qp4 + 
     2*m32*m32*p1p4*q2*qp2*qp4 + 16*M2*m32*p1p4*q2*qp2*qp4 - 
     4*m32*m32*p1p4*q2*qp2*qp4 + 2*m32*m32*p1p4*q2*qp2*qp4 - 
     2*m32*p1P*q2*q2*qp2*qp4 + 2*m32*p1P*q2*q2*qp2*qp4 - 
     8*M2*p1p3*q2*q2*qp2*qp4 + m32*p1p3*q2*q2*qp2*qp4 - 
     m32*p1p3*q2*q2*qp2*qp4 + 16*m32*m32*p1P*qP*qp2*qp4 - 
     16*m32*m32*p1P*qP*qp2*qp4 - 4*m32*m32*p1p3*qP*qp2*qp4 + 
     8*m32*m32*p1p3*qP*qp2*qp4 - 4*m32*m32*p1p3*qP*qp2*qp4 + 
     12*m32*m32*p1p4*qP*qp2*qp4 - 24*m32*m32*p1p4*qP*qp2*qp4 + 
     12*m32*m32*p1p4*qP*qp2*qp4 + 16*m32*p1P*p3p4*qP*qp2*qp4 - 
     16*m32*p1P*p3p4*qP*qp2*qp4 + 16*m32*p1P*q2*qP*qp2*qp4 - 
     4*m32*p1p3*q2*qP*qp2*qp4 + 4*m32*p1p3*q2*qP*qp2*qp4 + 
     12*m32*p1p4*q2*qP*qp2*qp4 - 12*m32*p1p4*q2*qP*qp2*qp4 + 
     16*p1P*p3p4*q2*qP*qp2*qp4 + 4*m32*m32*m32*qp1*qp2*qp4 - 
     16*M2*m32*m32*qp1*qp2*qp4 - 12*m32*m32*m32*qp1*qp2*qp4 + 
     16*M2*m32*m32*qp1*qp2*qp4 + 12*m32*m32*m32*qp1*qp2*qp4 - 
     16*M2*m32*p3p4*qp1*qp2*qp4 + 16*M2*m32*p3p4*qp1*qp2*qp4 + 
     4*m32*m32*q2*qp1*qp2*qp4 - 16*M2*m32*q2*qp1*qp2*qp4 - 
     8*m32*m32*q2*qp1*qp2*qp4 + 4*m32*m32*q2*qp1*qp2*qp4 - 
     16*M2*p3p4*q2*qp1*qp2*qp4 + 40*M2*m12*m32*m32*qp3*qp4 + 
     2*m12*m32*m32*m32*qp3*qp4 - 80*M2*m12*m32*m32*qp3*qp4 + 
     2*m12*m32*m32*m32*qp3*qp4 + 40*M2*m12*m32*m32*qp3*qp4 - 
     10*m12*m32*m32*m32*qp3*qp4 + 8*M2*m32*m32*p1p2*qp3*qp4 - 
     16*M2*m32*m32*p1p2*qp3*qp4 - 4*m32*m32*m32*p1p2*qp3*qp4 + 
     8*M2*m32*m32*p1p2*qp3*qp4 + 8*m32*m32*m32*p1p2*qp3*qp4 - 
     80*m32*m32*p1P*p2P*qp3*qp4 + 160*m32*m32*p1P*p2P*qp3*qp4 - 
     80*m32*m32*p1P*p2P*qp3*qp4 + 4*m32*m32*p1p3*p2P*qp3*qp4 - 
     8*m32*m32*p1p3*p2P*qp3*qp4 + 4*m32*m32*p1p3*p2P*qp3*qp4 + 
     4*m32*m32*p1p4*p2P*qp3*qp4 - 8*m32*m32*p1p4*p2P*qp3*qp4 + 
     4*m32*m32*p1p4*p2P*qp3*qp4 + 4*m32*m32*p1P*p2p3*qp3*qp4 - 
     8*m32*m32*p1P*p2p3*qp3*qp4 + 4*m32*m32*p1P*p2p3*qp3*qp4 - 
     4*m32*m32*p1p4*p2p3*qp3*qp4 + 8*m32*m32*p1p4*p2p3*qp3*qp4 - 
     4*m32*m32*p1p4*p2p3*qp3*qp4 + 4*m32*m32*p1P*p2p4*qp3*qp4 - 
     8*m32*m32*p1P*p2p4*qp3*qp4 + 4*m32*m32*p1P*p2p4*qp3*qp4 - 
     4*m32*m32*p1p3*p2p4*qp3*qp4 + 8*m32*m32*p1p3*p2p4*qp3*qp4 - 
     4*m32*m32*p1p3*p2p4*qp3*qp4 + 4*m12*m32*m32*p3p4*qp3*qp4 - 
     8*m12*m32*m32*p3p4*qp3*qp4 + 4*m12*m32*m32*p3p4*qp3*qp4 + 
     32*M2*m12*m32*q2*qp3*qp4 + 4*m12*m32*m32*q2*qp3*qp4 - 
     32*M2*m12*m32*q2*qp3*qp4 - 8*m12*m32*m32*q2*qp3*qp4 + 
     4*m12*m32*m32*q2*qp3*qp4 + 16*M2*m32*p1p2*q2*qp3*qp4 - 
     32*M2*m32*p1p2*q2*qp3*qp4 - 64*m32*p1P*p2P*q2*qp3*qp4 + 
     96*m32*p1P*p2P*q2*qp3*qp4 + 12*m32*p1p3*p2P*q2*qp3*qp4 - 
     12*m32*p1p3*p2P*q2*qp3*qp4 + 12*m32*p1p4*p2P*q2*qp3*qp4 - 
     12*m32*p1p4*p2P*q2*qp3*qp4 + 12*m32*p1P*p2p3*q2*qp3*qp4 - 
     12*m32*p1P*p2p3*q2*qp3*qp4 + 16*M2*p1p4*p2p3*q2*qp3*qp4 + 
     12*m32*p1P*p2p4*q2*qp3*qp4 - 12*m32*p1P*p2p4*q2*qp3*qp4 + 
     16*M2*p1p3*p2p4*q2*qp3*qp4 + 16*M2*m12*p3p4*q2*qp3*qp4 - 
     32*M2*p1p2*p3p4*q2*qp3*qp4 + 32*p1P*p2P*p3p4*q2*qp3*qp4 + 
     8*M2*m12*q2*q2*qp3*qp4 - 2*m12*m32*q2*q2*qp3*qp4 + 
     2*m12*m32*q2*q2*qp3*qp4 + 8*M2*p1p2*q2*q2*qp3*qp4 - 
     16*p1P*p2P*q2*q2*qp3*qp4 - 16*m12*m32*m32*qP*qp3*qp4 + 
     64*m12*m32*m32*qP*qp3*qp4 - 48*m12*m32*m32*qP*qp3*qp4 - 
     16*m32*m32*p1p2*qP*qp3*qp4 + 16*m32*m32*p1p2*qP*qp3*qp4 - 
     32*m32*p1p3*p2P*qP*qp3*qp4 + 32*m32*p1p3*p2P*qP*qp3*qp4 - 
     32*m32*p1p4*p2P*qP*qp3*qp4 + 32*m32*p1p4*p2P*qP*qp3*qp4 - 
     32*m32*p1P*p2p3*qP*qp3*qp4 + 32*m32*p1P*p2p3*qP*qp3*qp4 - 
     16*m32*p1p4*p2p3*qP*qp3*qp4 + 16*m32*p1p4*p2p3*qP*qp3*qp4 - 
     32*m32*p1P*p2p4*qP*qp3*qp4 + 32*m32*p1P*p2p4*qP*qp3*qp4 - 
     16*m32*p1p3*p2p4*qP*qp3*qp4 + 16*m32*p1p3*p2p4*qP*qp3*qp4 + 
     16*m12*m32*p3p4*qP*qp3*qp4 - 16*m12*m32*p3p4*qP*qp3*qp4 - 
     16*p1p3*p2P*q2*qP*qp3*qp4 - 16*p1p4*p2P*q2*qP*qp3*qp4 - 
     16*p1P*p2p3*q2*qP*qp3*qp4 - 16*p1P*p2p4*q2*qP*qp3*qp4 + 
     32*m12*m32*qP*qP*qp3*qp4 - 32*p1p4*p2p3*qP*qP*qp3*qp4 - 
     32*p1p3*p2p4*qP*qP*qp3*qp4 + 32*p1p2*p3p4*qP*qP*qp3*qp4 + 
     24*m32*m32*p2P*qp1*qp3*qp4 - 48*m32*m32*p2P*qp1*qp3*qp4 + 
     24*m32*m32*p2P*qp1*qp3*qp4 + 24*M2*m32*p2p3*qp1*qp3*qp4 - 
     24*M2*m32*p2p3*qp1*qp3*qp4 + 24*M2*m32*p2p4*qp1*qp3*qp4 - 
     24*M2*m32*p2p4*qp1*qp3*qp4 + 8*m32*p2P*q2*qp1*qp3*qp4 - 
     8*m32*p2P*q2*qp1*qp3*qp4 + 24*M2*p2p3*q2*qp1*qp3*qp4 - 
     4*m32*p2p3*q2*qp1*qp3*qp4 + 4*m32*p2p3*q2*qp1*qp3*qp4 + 
     24*M2*p2p4*q2*qp1*qp3*qp4 - 4*m32*p2p4*q2*qp1*qp3*qp4 + 
     4*m32*p2p4*q2*qp1*qp3*qp4 - 32*m32*p2P*qP*qp1*qp3*qp4 + 
     8*m32*p2p3*qP*qp1*qp3*qp4 - 8*m32*p2p3*qP*qp1*qp3*qp4 + 
     8*m32*p2p4*qP*qp1*qp3*qp4 - 8*m32*p2p4*qP*qp1*qp3*qp4 - 
     32*p2P*p3p4*qP*qp1*qp3*qp4 + 24*m32*m32*p1P*qp2*qp3*qp4 - 
     48*m32*m32*p1P*qp2*qp3*qp4 + 24*m32*m32*p1P*qp2*qp3*qp4 + 
     24*M2*m32*p1p3*qp2*qp3*qp4 - 24*M2*m32*p1p3*qp2*qp3*qp4 + 
     24*M2*m32*p1p4*qp2*qp3*qp4 - 24*M2*m32*p1p4*qp2*qp3*qp4 + 
     8*m32*p1P*q2*qp2*qp3*qp4 - 8*m32*p1P*q2*qp2*qp3*qp4 + 
     24*M2*p1p3*q2*qp2*qp3*qp4 - 4*m32*p1p3*q2*qp2*qp3*qp4 + 
     4*m32*p1p3*q2*qp2*qp3*qp4 + 24*M2*p1p4*q2*qp2*qp3*qp4 - 
     4*m32*p1p4*q2*qp2*qp3*qp4 + 4*m32*p1p4*q2*qp2*qp3*qp4 - 
     32*m32*p1P*qP*qp2*qp3*qp4 + 8*m32*p1p3*qP*qp2*qp3*qp4 - 
     8*m32*p1p3*qP*qp2*qp3*qp4 + 8*m32*p1p4*qP*qp2*qp3*qp4 - 
     8*m32*p1p4*qP*qp2*qp3*qp4 - 32*p1P*p3p4*qP*qp2*qp3*qp4 - 
     8*m32*m32*qp1*qp2*qp3*qp4 + 32*M2*m32*qp1*qp2*qp3*qp4 + 
     16*m32*m32*qp1*qp2*qp3*qp4 - 8*m32*m32*qp1*qp2*qp3*qp4 + 
     32*M2*p3p4*qp1*qp2*qp3*qp4 - 32*M2*m12*m32*qp3*qp3*qp4 + 
     32*M2*m12*m32*qp3*qp3*qp4 - 16*M2*m32*p1p2*qp3*qp3*qp4 + 
     16*M2*m32*p1p2*qp3*qp3*qp4 + 64*m32*p1P*p2P*qp3*qp3*qp4 - 
     64*m32*p1P*p2P*qp3*qp3*qp4 - 16*m32*p1p3*p2P*qp3*qp3*qp4 + 
     16*m32*p1p3*p2P*qp3*qp3*qp4 + 8*m32*p1p4*p2P*qp3*qp3*qp4 - 
     8*m32*p1p4*p2P*qp3*qp3*qp4 - 16*m32*p1P*p2p3*qp3*qp3*qp4 + 
     16*m32*p1P*p2p3*qp3*qp3*qp4 + 8*m32*p1P*p2p4*qp3*qp3*qp4 - 
     8*m32*p1P*p2p4*qp3*qp3*qp4 - 16*M2*m12*q2*qp3*qp3*qp4 + 
     4*m12*m32*q2*qp3*qp3*qp4 - 4*m12*m32*q2*qp3*qp3*qp4 - 
     16*M2*p1p2*q2*qp3*qp3*qp4 + 32*p1P*p2P*q2*qp3*qp3*qp4 + 
     32*p1p4*p2P*qP*qp3*qp3*qp4 + 32*p1P*p2p4*qP*qp3*qp3*qp4 - 
     8*m32*p2P*qp1*qp3*qp3*qp4 + 8*m32*p2P*qp1*qp3*qp3*qp4 - 
     16*M2*p2p3*qp1*qp3*qp3*qp4 + 4*m32*p2p3*qp1*qp3*qp3*qp4 - 
     4*m32*p2p3*qp1*qp3*qp3*qp4 - 16*M2*p2p4*qp1*qp3*qp3*qp4 - 
     8*m32*p1P*qp2*qp3*qp3*qp4 + 8*m32*p1P*qp2*qp3*qp3*qp4 - 
     16*M2*p1p3*qp2*qp3*qp3*qp4 + 4*m32*p1p3*qp2*qp3*qp3*qp4 - 
     4*m32*p1p3*qp2*qp3*qp3*qp4 - 16*M2*p1p4*qp2*qp3*qp3*qp4 + 
     16*M2*m12*qp3*qp3*qp3*qp4 - 4*m12*m32*qp3*qp3*qp3*qp4 + 
     4*m12*m32*qp3*qp3*qp3*qp4 + 24*M2*m12*m32*m32*qp4*qp4 - 2*m12*m32*m32*m32*qp4*qp4 - 
     32*M2*m12*m32*m32*qp4*qp4 + 2*m12*m32*m32*m32*qp4*qp4 + 
     8*M2*m12*m32*m32*qp4*qp4 + 2*m12*m32*m32*m32*qp4*qp4 - 
     8*M2*m32*m32*p1p2*qp4*qp4 + 2*m32*m32*m32*p1p2*qp4*qp4 + 
     8*M2*m32*m32*p1p2*qp4*qp4 - 4*m32*m32*m32*p1p2*qp4*qp4 - 
     16*m32*m32*p1P*p2P*qp4*qp4 + 32*m32*m32*p1P*p2P*qp4*qp4 - 
     16*m32*m32*p1P*p2P*qp4*qp4 - 4*m32*m32*p1p3*p2P*qp4*qp4 + 
     8*m32*m32*p1p3*p2P*qp4*qp4 - 4*m32*m32*p1p3*p2P*qp4*qp4 + 
     8*m32*m32*p1p4*p2P*qp4*qp4 - 16*m32*m32*p1p4*p2P*qp4*qp4 + 
     8*m32*m32*p1p4*p2P*qp4*qp4 - 4*m32*m32*p1P*p2p3*qp4*qp4 + 
     8*m32*m32*p1P*p2p3*qp4*qp4 - 4*m32*m32*p1P*p2p3*qp4*qp4 - 
     8*M2*m32*p1p4*p2p3*qp4*qp4 + 2*m32*m32*p1p4*p2p3*qp4*qp4 + 
     8*M2*m32*p1p4*p2p3*qp4*qp4 - 4*m32*m32*p1p4*p2p3*qp4*qp4 + 
     2*m32*m32*p1p4*p2p3*qp4*qp4 + 8*m32*m32*p1P*p2p4*qp4*qp4 - 
     16*m32*m32*p1P*p2p4*qp4*qp4 + 8*m32*m32*p1P*p2p4*qp4*qp4 - 
     8*M2*m32*p1p3*p2p4*qp4*qp4 + 2*m32*m32*p1p3*p2p4*qp4*qp4 + 
     8*M2*m32*p1p3*p2p4*qp4*qp4 - 4*m32*m32*p1p3*p2p4*qp4*qp4 + 
     2*m32*m32*p1p3*p2p4*qp4*qp4 + 8*M2*m12*m32*p3p4*qp4*qp4 - 
     2*m12*m32*m32*p3p4*qp4*qp4 - 8*M2*m12*m32*p3p4*qp4*qp4 + 
     4*m12*m32*m32*p3p4*qp4*qp4 - 2*m12*m32*m32*p3p4*qp4*qp4 + 
     24*M2*m12*m32*q2*qp4*qp4 - 2*m12*m32*m32*q2*qp4*qp4 - 
     40*M2*m12*m32*q2*qp4*qp4 + 8*m12*m32*m32*q2*qp4*qp4 - 
     6*m12*m32*m32*q2*qp4*qp4 + 8*M2*m32*p1p2*q2*qp4*qp4 - 
     2*m32*m32*p1p2*q2*qp4*qp4 + 2*m32*m32*p1p2*q2*qp4*qp4 - 
     16*m32*p1P*p2P*q2*qp4*qp4 + 16*m32*p1P*p2P*q2*qp4*qp4 - 
     4*m32*p1p3*p2P*q2*qp4*qp4 + 4*m32*p1p3*p2P*q2*qp4*qp4 + 
     8*m32*p1p4*p2P*q2*qp4*qp4 - 8*m32*p1p4*p2P*q2*qp4*qp4 - 
     4*m32*p1P*p2p3*q2*qp4*qp4 + 4*m32*p1P*p2p3*q2*qp4*qp4 + 
     8*M2*p1p4*p2p3*q2*qp4*qp4 - 2*m32*p1p4*p2p3*q2*qp4*qp4 + 
     2*m32*p1p4*p2p3*q2*qp4*qp4 + 8*m32*p1P*p2p4*q2*qp4*qp4 - 
     8*m32*p1P*p2p4*q2*qp4*qp4 + 8*M2*p1p3*p2p4*q2*qp4*qp4 - 
     2*m32*p1p3*p2p4*q2*qp4*qp4 + 2*m32*p1p3*p2p4*q2*qp4*qp4 - 
     8*M2*m12*p3p4*q2*qp4*qp4 + 2*m12*m32*p3p4*q2*qp4*qp4 - 
     2*m12*m32*p3p4*q2*qp4*qp4 - 16*m12*m32*m32*qP*qp4*qp4 + 
     16*m12*m32*m32*qP*qp4*qp4 + 16*m32*m32*p1p2*qP*qp4*qp4 - 
     16*m32*m32*p1p2*qP*qp4*qp4 - 16*m32*p1p3*p2P*qP*qp4*qp4 + 
     16*m32*p1p3*p2P*qP*qp4*qp4 - 16*m32*p1P*p2p3*qP*qp4*qp4 + 
     16*m32*p1P*p2p3*qP*qp4*qp4 + 16*m32*p1p4*p2p3*qP*qp4*qp4 - 
     16*m32*p1p4*p2p3*qP*qp4*qp4 + 16*m32*p1p3*p2p4*qP*qp4*qp4 - 
     16*m32*p1p3*p2p4*qP*qp4*qp4 - 16*m12*m32*p3p4*qP*qp4*qp4 + 
     16*m12*m32*p3p4*qP*qp4*qp4 - 16*m12*m32*q2*qP*qp4*qp4 + 
     16*m12*m32*q2*qP*qp4*qp4 - 16*p1p3*p2P*q2*qP*qp4*qp4 - 
     16*p1P*p2p3*q2*qP*qp4*qp4 + 4*m32*m32*p2P*qp1*qp4*qp4 - 
     8*m32*m32*p2P*qp1*qp4*qp4 + 4*m32*m32*p2P*qp1*qp4*qp4 + 
     8*M2*m32*p2p3*qp1*qp4*qp4 - 8*M2*m32*p2p3*qp1*qp4*qp4 + 
     16*M2*m32*p2p4*qp1*qp4*qp4 - 4*m32*m32*p2p4*qp1*qp4*qp4 - 
     16*M2*m32*p2p4*qp1*qp4*qp4 + 8*m32*m32*p2p4*qp1*qp4*qp4 - 
     4*m32*m32*p2p4*qp1*qp4*qp4 + 4*m32*p2P*q2*qp1*qp4*qp4 - 
     4*m32*p2P*q2*qp1*qp4*qp4 + 8*M2*p2p3*q2*qp1*qp4*qp4 - 
     16*m32*p2p4*qP*qp1*qp4*qp4 + 16*m32*p2p4*qP*qp1*qp4*qp4 + 
     4*m32*m32*p1P*qp2*qp4*qp4 - 8*m32*m32*p1P*qp2*qp4*qp4 + 
     4*m32*m32*p1P*qp2*qp4*qp4 + 8*M2*m32*p1p3*qp2*qp4*qp4 - 
     8*M2*m32*p1p3*qp2*qp4*qp4 + 16*M2*m32*p1p4*qp2*qp4*qp4 - 
     4*m32*m32*p1p4*qp2*qp4*qp4 - 16*M2*m32*p1p4*qp2*qp4*qp4 + 
     8*m32*m32*p1p4*qp2*qp4*qp4 - 4*m32*m32*p1p4*qp2*qp4*qp4 + 
     4*m32*p1P*q2*qp2*qp4*qp4 - 4*m32*p1P*q2*qp2*qp4*qp4 + 
     8*M2*p1p3*q2*qp2*qp4*qp4 - 16*m32*p1p4*qP*qp2*qp4*qp4 + 
     16*m32*p1p4*qP*qp2*qp4*qp4 - 32*M2*m12*m32*qp3*qp4*qp4 + 
     32*M2*m12*m32*qp3*qp4*qp4 - 16*M2*m32*p1p2*qp3*qp4*qp4 + 
     16*M2*m32*p1p2*qp3*qp4*qp4 + 64*m32*p1P*p2P*qp3*qp4*qp4 - 
     64*m32*p1P*p2P*qp3*qp4*qp4 + 8*m32*p1p3*p2P*qp3*qp4*qp4 - 
     8*m32*p1p3*p2P*qp3*qp4*qp4 - 16*m32*p1p4*p2P*qp3*qp4*qp4 + 
     16*m32*p1p4*p2P*qp3*qp4*qp4 + 8*m32*p1P*p2p3*qp3*qp4*qp4 - 
     8*m32*p1P*p2p3*qp3*qp4*qp4 - 16*m32*p1P*p2p4*qp3*qp4*qp4 + 
     16*m32*p1P*p2p4*qp3*qp4*qp4 - 16*M2*m12*q2*qp3*qp4*qp4 + 
     4*m12*m32*q2*qp3*qp4*qp4 - 4*m12*m32*q2*qp3*qp4*qp4 - 
     16*M2*p1p2*q2*qp3*qp4*qp4 + 32*p1P*p2P*q2*qp3*qp4*qp4 + 
     32*p1p3*p2P*qP*qp3*qp4*qp4 + 32*p1P*p2p3*qP*qp3*qp4*qp4 - 
     8*m32*p2P*qp1*qp3*qp4*qp4 + 8*m32*p2P*qp1*qp3*qp4*qp4 - 
     16*M2*p2p3*qp1*qp3*qp4*qp4 - 16*M2*p2p4*qp1*qp3*qp4*qp4 + 
     4*m32*p2p4*qp1*qp3*qp4*qp4 - 4*m32*p2p4*qp1*qp3*qp4*qp4 - 
     8*m32*p1P*qp2*qp3*qp4*qp4 + 8*m32*p1P*qp2*qp3*qp4*qp4 - 
     16*M2*p1p3*qp2*qp3*qp4*qp4 - 16*M2*p1p4*qp2*qp3*qp4*qp4 + 
     4*m32*p1p4*qp2*qp3*qp4*qp4 - 4*m32*p1p4*qp2*qp3*qp4*qp4 + 
     32*M2*p1p2*qp3*qp3*qp4*qp4 - 64*p1P*p2P*qp3*qp3*qp4*qp4 - 
     16*M2*m12*m32*qp4*qp4*qp4 + 4*m12*m32*m32*qp4*qp4*qp4 + 
     16*M2*m12*m32*qp4*qp4*qp4 - 8*m12*m32*m32*qp4*qp4*qp4 + 
     4*m12*m32*m32*qp4*qp4*qp4 + 16*m12*m32*qP*qp4*qp4*qp4 - 16*m12*m32*qP*qp4*qp4*qp4 + 
     16*M2*m12*qp3*qp4*qp4*qp4 - 4*m12*m32*qp3*qp4*qp4*qp4 + 
     4*m12*m32*qp3*qp4*qp4*qp4 + 
     2*p3P*(m32 - m32 + q2 - 2*qp4)*(-2*m12*m32*m32*q2 + 
        4*m12*m32*m32*q2 - 2*m12*m32*m32*q2 - 2*m12*m32*q2*q2 + 
        2*m12*m32*q2*q2 + 4*m12*m32*m32*qP - 8*m12*m32*m32*qP + 
        4*m12*m32*m32*qP + 4*m12*m32*q2*qP - 4*m12*m32*q2*qP - 
        2*m32*m32*p2p3*qp1 + 4*m32*m32*p2p3*qp1 - 2*m32*m32*p2p3*qp1 + 
        m32*m32*p2p4*qp1 - 2*m32*m32*p2p4*qp1 + m32*m32*p2p4*qp1 - 
        2*m32*p2p3*q2*qp1 + 2*m32*p2p3*q2*qp1 + m32*p2p4*q2*qp1 - 
        m32*p2p4*q2*qp1 + 4*m32*p2p4*qP*qp1 - 4*m32*p2p4*qP*qp1 + 
        4*p2p4*q2*qP*qp1 - 2*m32*m32*p1p3*qp2 + 4*m32*m32*p1p3*qp2 - 
        2*m32*m32*p1p3*qp2 + m32*m32*p1p4*qp2 - 2*m32*m32*p1p4*qp2 + 
        m32*m32*p1p4*qp2 - 2*m32*p1p3*q2*qp2 + 2*m32*p1p3*q2*qp2 + 
        m32*p1p4*q2*qp2 - m32*p1p4*q2*qp2 + 4*m32*p1p4*qP*qp2 - 
        4*m32*p1p4*qP*qp2 + 4*p1p4*q2*qP*qp2 + 2*m32*m32*qp1*qp2 - 
        4*m32*m32*qp1*qp2 + 2*m32*m32*qp1*qp2 + 2*m32*q2*qp1*qp2 - 
        2*m32*q2*qp1*qp2 + 4*m12*m32*m32*qp3 - 4*m12*m32*m32*qp3 - 
        4*m32*m32*p1p2*qp3 + 4*m32*m32*p1p2*qp3 - 
        4*m32*p1p4*p2p3*qp3 + 4*m32*p1p4*p2p3*qp3 - 
        4*m32*p1p3*p2p4*qp3 + 4*m32*p1p3*p2p4*qp3 + 
        4*m12*m32*p3p4*qp3 - 4*m12*m32*p3p4*qp3 + 
        8*m12*m32*q2*qp3 - 8*m12*m32*q2*qp3 - 8*m12*m32*qP*qp3 + 
        24*m12*m32*qP*qp3 - 8*m32*p1p2*qP*qp3 - 8*p1p4*p2p3*qP*qp3 - 
        8*p1p3*p2p4*qP*qp3 + 8*m12*p3p4*qP*qp3 + 
        4*m32*p2p3*qp1*qp3 - 4*m32*p2p3*qp1*qp3 + 
        2*m32*p2p4*qp1*qp3 - 2*m32*p2p4*qp1*qp3 + 
        4*m32*p1p3*qp2*qp3 - 4*m32*p1p3*qp2*qp3 + 
        2*m32*p1p4*qp2*qp3 - 2*m32*p1p4*qp2*qp3 - 
        4*m32*qp1*qp2*qp3 + 4*m32*qp1*qp2*qp3 - 8*m12*m32*qp3*qp3 + 
        8*m12*m32*qp3*qp3 - 2*m12*m32*m32*qp4 - 4*m12*m32*m32*qp4 + 
        6*m12*m32*m32*qp4 + 4*m32*m32*p1p2*qp4 - 4*m32*m32*p1p2*qp4 + 
        4*m32*p1p4*p2p3*qp4 - 4*m32*p1p4*p2p3*qp4 + 
        4*m32*p1p3*p2p4*qp4 - 4*m32*p1p3*p2p4*qp4 - 
        4*m12*m32*p3p4*qp4 + 4*m12*m32*p3p4*qp4 - 
        2*m12*m32*q2*qp4 + 2*m12*m32*q2*qp4 - 4*m12*m32*qP*qp4 - 
        12*m12*m32*qP*qp4 - 4*m32*p1p2*qP*qp4 + 12*m32*p1p2*qP*qp4 + 
        8*p1p4*p2p3*qP*qp4 + 8*p1p3*p2p4*qP*qp4 - 8*m12*p3p4*qP*qp4 - 
        4*m12*q2*qP*qp4 - 4*p1p2*q2*qP*qp4 - 4*m32*p2p4*qp1*qp4 + 
        4*m32*p2p4*qp1*qp4 - 8*p2p4*qP*qp1*qp4 - 
        4*m32*p1p4*qp2*qp4 + 4*m32*p1p4*qp2*qp4 - 
        8*p1p4*qP*qp2*qp4 + 8*p1p2*qP*qp3*qp4 + 4*m12*m32*qp4*qp4 - 
        4*m12*m32*qp4*qp4 + 8*m12*qP*qp4*qp4 + 
        2*p2P*(m32 - m32 + q2 - 2*qp3)*(-2*p1p4*q2 + 
           m32*(2*p1p3 - p1p4 - qp1) + m32*(-2*p1p3 + p1p4 + qp1) - 
           2*p1p3*qp4 + 2*p1p4*qp4 + 2*qp1*qp4) + 
        2*p1P*(m32 - m32 + q2 - 2*qp3)*(-2*p2p4*q2 + 
           m32*(2*p2p3 - p2p4 - qp2) + m32*(-2*p2p3 + p2p4 + qp2) - 
           2*p2p3*qp4 + 2*p2p4*qp4 + 2*qp2*qp4)) + 
     8*p3P*p3P*(m32 - m32 + q2 - 
        2*qp4)*((m32*p1p2 + p2p4*(p1p3 - qp1) + 
           p1p4*(p2p3 - qp2))*(m32 - m32 + q2 - 2*qp4) - 
        m12*(m32*m32 - m32*m32 + p3p4*q2 + 
           m32*(-p3p4 + q2 + 2*qp3 - 2*qp4) - 2*p3p4*qp4 - 2*q2*qp4 + 
           2*qp3*qp4 + 2*qp4*qp4 + m32*(p3p4 + q2 - 2*(qp3 + qp4)))) + 
     8*p4P*p4P*(m32 - m32 + q2 - 
        2*qp3)*((m32*p1p2 + p1p4*p2p3 + p1p3*p2p4 - p2p3*qp1 - 
           p1p3*qp2)*(m32 - m32 + q2 - 2*qp3) - 
        m12*(m32*m32 - m32*m32 + p3p4*q2 - 2*p3p4*qp3 - 2*q2*qp3 + 
           2*qp3*qp3 + 2*qp3*qp4 + m32*(-p3p4 + q2 - 2*qp3 + 2*qp4) + 
           m32*(p3p4 + q2 - 2*(qp3 + qp4)))) - 
     2*p4P*(m32 - m32 + q2 - 
        2*qp3)*(m32*m32*(2*p1p3*p2P - 4*p1p4*p2P + 2*p1P*p2p3 - 
           4*p1P*p2p4 - 4*m12*p3P + 2*m12*q2 - 4*m12*qP + 2*p2P*qp1 - 
           p2p3*qp1 + 2*p2p4*qp1 + 2*p1P*qp2 - p1p3*qp2 + 
           2*p1p4*qp2 - 2*qp1*qp2 + 2*m12*qp3 - 4*m12*qp4) + 
        m32*m32*(2*p1p3*p2P - 4*p1p4*p2P + 2*p1P*p2p3 - 4*p1P*p2p4 + 
           12*m12*p3P - 8*p1p2*p3P + 2*m12*q2 - 4*m12*qP + 2*p2P*qp1 -
            p2p3*qp1 + 2*p2p4*qp1 + 2*p1P*qp2 - p1p3*qp2 + 
           2*p1p4*qp2 - 2*qp1*qp2 - 6*m12*qp3 + 4*p1p2*qp3 + 
           4*m12*qp4 - 4*p1p2*qp4) - 
        4*(2*m12*p3P*p3p4*q2 - p1P*p2p3*q2*q2 + m12*p3P*q2*q2 + 
           p1p2*p3P*q2*q2 + p2p3*p3P*q2*qp1 + p2p4*p3P*q2*qp1 + 
           p2p3*q2*qP*qp1 - 2*p3P*q2*qp1*qp2 + p1P*p2p3*q2*qp3 - 
           p1P*p2p4*q2*qp3 - 2*m12*p3P*q2*qp3 - 2*m12*p3p4*qP*qp3 - 
           m12*q2*qP*qp3 - p1p2*q2*qP*qp3 + p2P*q2*qp1*qp3 - 
           2*p2p3*qP*qp1*qp3 + p1P*q2*qp2*qp3 + 2*m12*qP*qp3*qp3 + 
           p1p3*qp2*(qP*(q2 - 2*qp3) + p3P*(q2 - 2*qp4)) - 
           p1p3*p2P*(q2 - qp3)*(q2 - 2*qp4) + 
           p1p4*(p3P*qp2 - p2P*qp3)*(q2 - 2*qp4) - 
           4*m12*p3P*p3p4*qp4 + 2*p1P*p2p3*q2*qp4 - 
           4*m12*p3P*q2*qp4 - 2*p1p2*p3P*q2*qp4 + 2*m12*p3p4*qP*qp4 - 
           2*p2p3*p3P*qp1*qp4 - 2*p2p4*p3P*qp1*qp4 + 
           4*p3P*qp1*qp2*qp4 - 2*p1P*p2p3*qp3*qp4 + 
           2*p1P*p2p4*qp3*qp4 + 4*m12*p3P*qp3*qp4 + 
           2*p1p2*qP*qp3*qp4 - 2*p2P*qp1*qp3*qp4 - 
           2*p1P*qp2*qp3*qp4 + 4*m12*p3P*qp4*qp4 - 
           2*p1p4*p2p3*(p3P*(q2 - 2*qp4) + qP*(-qp3 + qp4)) - 
           2*p1p3*p2p4*(p3P*(q2 - 2*qp4) + qP*(-qp3 + qp4))) - 
        m32*(-8*p1p3*p2p4*p3P + 8*m12*p3P*p3p4 - 6*p1p3*p2P*q2 - 
           6*p1P*p2p3*q2 + 4*p1P*p2p4*q2 + 8*m12*p3P*q2 + 
           4*p1p2*p3P*q2 - 2*m12*q2*q2 + 4*m12*q2*qP + 
           4*p2p3*p3P*qp1 + 4*p2p4*p3P*qp1 - 2*p2P*q2*qp1 + 
           p2p3*q2*qp1 - 2*p2p4*q2*qp1 + 4*p2p3*qP*qp1 + 
           4*p1p3*p3P*qp2 - 2*p1P*q2*qp2 + p1p3*q2*qp2 + 
           4*p1p3*qP*qp2 - 8*p3P*qp1*qp2 + 2*q2*qp1*qp2 + 
           4*p1p3*p2P*qp3 + 4*p1P*p2p3*qp3 - 4*p1P*p2p4*qp3 + 
           4*p1p3*p2p4*qp3 - 8*m12*p3P*qp3 - 4*m12*p3p4*qp3 - 
           2*m12*q2*qp3 - 4*m12*qP*qp3 - 4*p1p2*qP*qp3 + 
           4*p2P*qp1*qp3 - 4*p2p3*qp1*qp3 + 4*p1P*qp2*qp3 - 
           4*p1p3*qp2*qp3 + 4*m12*qp3*qp3 + 4*p1p3*p2P*qp4 + 
           4*p1P*p2p3*qp4 - 8*p1P*p2p4*qp4 - 4*p1p3*p2p4*qp4 - 
           16*m12*p3P*qp4 + 4*m12*p3p4*qp4 + 8*m12*q2*qp4 - 
           8*m12*qP*qp4 + 4*p2P*qp1*qp4 + 2*p2p3*qp1*qp4 + 
           4*p2p4*qp1*qp4 + 4*p1P*qp2*qp4 + 2*p1p3*qp2*qp4 - 
           4*qp1*qp2*qp4 - 8*m12*qp4*qp4 + 
           m32*(4*p1p3*p2P - 8*p1p4*p2P + 4*p1P*p2p3 - 8*p1P*p2p4 + 
              8*m12*p3P - 8*p1p2*p3P + 4*m12*q2 - 8*m12*qP + 
              4*p2P*qp1 - 2*p2p3*qp1 + 4*p2p4*qp1 + 4*p1P*qp2 - 
              2*p1p3*qp2 + 4*p1p4*qp2 - 4*qp1*qp2 - 4*m12*qp3 + 
              4*p1p2*qp3 - 4*p1p2*qp4) - 
           2*p1p4*(qp2*(-2*p3P + q2 - 2*qp4) + 
              2*p2p3*(2*p3P - qp3 + qp4) + 
              p2P*(-2*q2 + 2*qp3 + 4*qp4))) + 
        m32*(8*m12*p3P*p3p4 - 6*p1P*p2p3*q2 + 4*p1P*p2p4*q2 - 
           8*m12*p3P*q2 + 12*p1p2*p3P*q2 - 2*m12*q2*q2 + 4*m12*q2*qP + 
           4*p2p3*p3P*qp1 + 4*p2p4*p3P*qp1 - 2*p2P*q2*qp1 + 
           p2p3*q2*qp1 - 2*p2p4*q2*qp1 + 4*p2p3*qP*qp1 - 
           2*p1P*q2*qp2 - 8*p3P*qp1*qp2 + 2*q2*qp1*qp2 + 
           4*p1P*p2p3*qp3 - 4*p1P*p2p4*qp3 - 8*m12*p3P*qp3 - 
           4*m12*p3p4*qp3 - 2*m12*q2*qp3 + 12*m12*qP*qp3 - 
           12*p1p2*qP*qp3 + 4*p2P*qp1*qp3 - 4*p2p3*qp1*qp3 + 
           4*p1P*qp2*qp3 + 4*m12*qp3*qp3 + 4*p1P*p2p3*qp4 - 
           8*p1P*p2p4*qp4 + 16*m12*p3P*qp4 - 16*p1p2*p3P*qp4 + 
           4*m12*p3p4*qp4 + 8*m12*q2*qp4 - 24*m12*qP*qp4 + 
           8*p1p2*qP*qp4 + 4*p2P*qp1*qp4 + 2*p2p3*qp1*qp4 + 
           4*p2p4*qp1*qp4 + 4*p1P*qp2*qp4 - 4*qp1*qp2*qp4 - 
           8*m12*qp4*qp4 + 
           p1p3*(-6*p2P*q2 - 4*p2p4*(2*p3P - qp3 + qp4) + 
              4*p2P*(qp3 + qp4) + 
              qp2*(4*p3P + q2 + 4*qP - 4*qp3 + 2*qp4)) - 
           2*p1p4*(qp2*(-2*p3P + q2 - 2*qp4) + 
              2*p2p3*(2*p3P - qp3 + qp4) + 
              p2P*(-2*q2 + 2*qp3 + 4*qp4)))) + 
     8*M2*m12*m32*m32*m32*m32 - 12*M2*m32*p1p2*m32*m32*m32 + 
     32*m32*p1P*p2P*m32*m32*m32 + 4*M2*p1p4*p2p3*m32*m32*m32 + 
     4*M2*p1p3*p2p4*m32*m32*m32 - 4*M2*m12*p3p4*m32*m32*m32 + 
     2*m12*m32*q2*m32*m32*m32 - 8*M2*p1p2*q2*m32*m32*m32 - 
     3*m32*p1p2*q2*m32*m32*m32 + 24*p1P*p2P*q2*m32*m32*m32 - 
     2*p1p3*p2P*q2*m32*m32*m32 - 2*p1p4*p2P*q2*m32*m32*m32 - 
     2*p1P*p2p3*q2*m32*m32*m32 + p1p4*p2p3*q2*m32*m32*m32 - 
     2*p1P*p2p4*q2*m32*m32*m32 + p1p3*p2p4*q2*m32*m32*m32 - 
     m12*p3p4*q2*m32*m32*m32 + 2*m12*q2*q2*m32*m32*m32 - 
     2*p1p2*q2*q2*m32*m32*m32 - 8*m12*m32*qP*m32*m32*m32 + 
     12*m32*p1p2*qP*m32*m32*m32 + 4*p1p3*p2P*qP*m32*m32*m32 + 
     4*p1p4*p2P*qP*m32*m32*m32 + 4*p1P*p2p3*qP*m32*m32*m32 - 
     4*p1p4*p2p3*qP*m32*m32*m32 + 4*p1P*p2p4*qP*m32*m32*m32 - 
     4*p1p3*p2p4*qP*m32*m32*m32 + 4*m12*p3p4*qP*m32*m32*m32 - 
     8*m12*q2*qP*m32*m32*m32 + 8*p1p2*q2*qP*m32*m32*m32 + 
     8*m12*qP*qP*m32*m32*m32 - 16*m32*p2P*qp1*m32*m32*m32 - 
     4*M2*p2p3*qp1*m32*m32*m32 - 4*M2*p2p4*qp1*m32*m32*m32 - 
     8*p2P*q2*qp1*m32*m32*m32 - 8*p2P*qP*qp1*m32*m32*m32 + 
     2*p2p3*qP*qp1*m32*m32*m32 + 2*p2p4*qP*qp1*m32*m32*m32 - 
     16*m32*p1P*qp2*m32*m32*m32 - 4*M2*p1p3*qp2*m32*m32*m32 - 
     4*M2*p1p4*qp2*m32*m32*m32 - 8*p1P*q2*qp2*m32*m32*m32 - 
     8*p1P*qP*qp2*m32*m32*m32 + 2*p1p3*qP*qp2*m32*m32*m32 + 
     2*p1p4*qP*qp2*m32*m32*m32 + 8*M2*qp1*qp2*m32*m32*m32 + 
     8*m32*qp1*qp2*m32*m32*m32 + 4*q2*qp1*qp2*m32*m32*m32 + 
     8*M2*p1p2*qp3*m32*m32*m32 - 24*p1P*p2P*qp3*m32*m32*m32 + 
     4*p1p3*p2P*qp3*m32*m32*m32 - 2*p1p4*p2P*qp3*m32*m32*m32 + 
     4*p1P*p2p3*qp3*m32*m32*m32 - 2*p1P*p2p4*qp3*m32*m32*m32 - 
     2*m12*q2*qp3*m32*m32*m32 + 2*p1p2*q2*qp3*m32*m32*m32 + 
     4*m12*qP*qp3*m32*m32*m32 - 8*p1p2*qP*qp3*m32*m32*m32 + 
     10*p2P*qp1*qp3*m32*m32*m32 - 2*p2p3*qp1*qp3*m32*m32*m32 + 
     p2p4*qp1*qp3*m32*m32*m32 + 10*p1P*qp2*qp3*m32*m32*m32 - 
     2*p1p3*qp2*qp3*m32*m32*m32 + p1p4*qp2*qp3*m32*m32*m32 - 
     4*qp1*qp2*qp3*m32*m32*m32 - 2*m12*qp3*qp3*m32*m32*m32 + 
     2*p1p2*qp3*qp3*m32*m32*m32 + 8*M2*p1p2*qp4*m32*m32*m32 - 
     24*p1P*p2P*qp4*m32*m32*m32 - 2*p1p3*p2P*qp4*m32*m32*m32 + 
     4*p1p4*p2P*qp4*m32*m32*m32 - 2*p1P*p2p3*qp4*m32*m32*m32 + 
     4*p1P*p2p4*qp4*m32*m32*m32 - 2*m12*q2*qp4*m32*m32*m32 + 
     2*p1p2*q2*qp4*m32*m32*m32 + 4*m12*qP*qp4*m32*m32*m32 - 
     8*p1p2*qP*qp4*m32*m32*m32 + 10*p2P*qp1*qp4*m32*m32*m32 + 
     p2p3*qp1*qp4*m32*m32*m32 - 2*p2p4*qp1*qp4*m32*m32*m32 + 
     10*p1P*qp2*qp4*m32*m32*m32 + p1p3*qp2*qp4*m32*m32*m32 - 
     2*p1p4*qp2*qp4*m32*m32*m32 - 4*qp1*qp2*qp4*m32*m32*m32 + 
     6*m12*qp3*qp4*m32*m32*m32 - 4*p1p2*qp3*qp4*m32*m32*m32 - 
     2*m12*qp4*qp4*m32*m32*m32 + 2*p1p2*qp4*qp4*m32*m32*m32 - 
     4*M2*m12*m32*m32*m32*m32 + 4*M2*p1p2*m32*m32*m32*m32 - 
     8*p1P*p2P*m32*m32*m32*m32 - m12*q2*m32*m32*m32*m32 + 
     p1p2*q2*m32*m32*m32*m32 + 4*m12*qP*m32*m32*m32*m32 - 
     4*p1p2*qP*m32*m32*m32*m32 + 4*p2P*qp1*m32*m32*m32*m32 + 
     4*p1P*qp2*m32*m32*m32*m32 - 2*qp1*qp2*m32*m32*m32*m32))/(q2*q2*(m32 -
      m32 + q2 - 2*qp3)*(m32 -
      m32 + q2 - 2*qp3)*(m32 - m32 + q2 - 2*qp4)*(m32 - m32 + q2 - 2*qp4)*(2*m32 + 
     2*p3p4 + q2 - 2*qp3 - 2*qp4)*(2*m32 + 
     2*p3p4 + q2 - 2*qp3 - 2*qp4));  */
         
      if(result < 0.0){result = 0.0;}       

    return result;
}


