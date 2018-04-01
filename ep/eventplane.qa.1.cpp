hQred1->Fill( 0.0, q2bb[0].Reduced() );
hQred1->Fill( 1.0, q2bb[1].Reduced() );
hQred1->Fill( 2.0, q2bb[2].Reduced() );
hQred1->Fill( 3.0, q2fv[0].Reduced() );
hQred1->Fill( 4.0, q2ex[0].Reduced() );
hQred1->Fill( 5.0, q2ex[1].Reduced() );
hQred1->Fill( 6.0, q2ex[2].Reduced() );
hQred1->Fill( 7.0, q2ex[3].Reduced() );
hQred1->Fill( 8.0, q2ex[4].Reduced() );
hQred1->Fill( 9.0, q2ex[5].Reduced() );
hQred1->Fill( 10., q2ex[6].Reduced() );
hQred1->Fill( 11., q2ex[7].Reduced() );
hQred1->Fill( 12., q2ex[8].Reduced() );
hQred1->Fill( 13., q2ca[0].Reduced() );

hPSI1->Fill( 0.0, q2bb[0].Psi2Pi() );
hPSI1->Fill( 1.0, q2bb[1].Psi2Pi() );
hPSI1->Fill( 2.0, q2bb[2].Psi2Pi() );
hPSI1->Fill( 3.0, q2fv[0].Psi2Pi() );
hPSI1->Fill( 4.0, q2ex[0].Psi2Pi() );
hPSI1->Fill( 5.0, q2ex[1].Psi2Pi() );
hPSI1->Fill( 6.0, q2ex[2].Psi2Pi() );
hPSI1->Fill( 7.0, q2ex[3].Psi2Pi() );
hPSI1->Fill( 8.0, q2ex[4].Psi2Pi() );
hPSI1->Fill( 9.0, q2ex[5].Psi2Pi() );
hPSI1->Fill( 10., q2ex[6].Psi2Pi() );
hPSI1->Fill( 11., q2ex[7].Psi2Pi() );
hPSI1->Fill( 12., q2ex[8].Psi2Pi() );
hPSI1->Fill( 13., q2ca[0].Psi2Pi() );

// BBAC BBCBC BBAB
// CAAC CABC CAAB
// FVBB FVCA
// EX0BB EX0CA EX0FV (x9)
// = 32
hEVC[1][ 0]->Fill( q2bb[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][ 1]->Fill( q2bb[1].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][ 2]->Fill( q2bb[0].Psi2Pi(), q2bb[1].Psi2Pi() );
//-
hEVC[1][ 3]->Fill( q2fv[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][ 4]->Fill( q2fv[0].Psi2Pi(), q2ca[0].Psi2Pi() );
//-
hEVC[1][ 5]->Fill( q2ex[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][ 6]->Fill( q2ex[0].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][ 7]->Fill( q2ex[0].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][ 8]->Fill( q2ex[1].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][ 9]->Fill( q2ex[1].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][10]->Fill( q2ex[1].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][11]->Fill( q2ex[2].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][12]->Fill( q2ex[2].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][13]->Fill( q2ex[2].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][14]->Fill( q2ex[3].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][15]->Fill( q2ex[3].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][16]->Fill( q2ex[3].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][17]->Fill( q2ex[4].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][18]->Fill( q2ex[4].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][19]->Fill( q2ex[4].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][20]->Fill( q2ex[5].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][21]->Fill( q2ex[5].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][22]->Fill( q2ex[5].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][23]->Fill( q2ex[6].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][24]->Fill( q2ex[6].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][25]->Fill( q2ex[6].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[1][26]->Fill( q2ex[7].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[1][27]->Fill( q2ex[7].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[1][28]->Fill( q2ex[7].Psi2Pi(), q2fv[0].Psi2Pi() );

//hEVC[1][29]->Fill( q2ex[8].Psi2Pi(), q2bb[2].Psi2Pi() );
//hEVC[1][30]->Fill( q2ex[8].Psi2Pi(), q2ca[0].Psi2Pi() );
//hEVC[1][31]->Fill( q2ex[8].Psi2Pi(), q2fv[0].Psi2Pi() );
