hBBCX->Fill( q2bb[0].X(), q2bb[1].X() );
hBBCY->Fill( q2bb[0].Y(), q2bb[1].Y() );
hBBCL->Fill( TMath::Sqrt( q2bb[0].X()*q2bb[0].X()+q2bb[0].Y()*q2bb[0].Y() ),
	     TMath::Sqrt( q2bb[1].X()*q2bb[1].X()+q2bb[1].Y()*q2bb[1].Y() ) );
hBBCH->Fill( TMath::Sqrt( (q2bb[0].X()*q2bb[0].X()+q2bb[0].Y()*q2bb[0].Y())/q2bb[0].M() ),
	     TMath::Sqrt( (q2bb[1].X()*q2bb[1].X()+q2bb[1].Y()*q2bb[1].Y())/q2bb[1].M() ) );
hBBCS->Fill( q2bb[0].Y()/(q2bb[0].X()+1e-8), q2bb[1].Y()/(q2bb[0].X()+1e-8) );
hBBCM->Fill( q2bb[0].M(), q2bb[1].M() );
hBBCNP->Fill( q2bb[0].NP(), q2bb[1].NP() );
hBBCatan01->Fill( TMath::ATan2( q2bb[0].Y(), q2bb[0].X() )/2.0,
		  TMath::ATan2( q2bb[1].Y(), q2bb[1].X() )/2.0 );
hBBCatan02->Fill( TMath::ATan2( q2bb[0].Y(), q2bb[0].X() )/2.0,
		  TMath::ATan2( q2bb[2].Y(), q2bb[2].X() )/2.0 );
hBBCatan12->Fill( TMath::ATan2( q2bb[1].Y(), q2bb[1].X() )/2.0,
		  TMath::ATan2( q2bb[2].Y(), q2bb[2].X() )/2.0 );

hQM->Fill(0.0,q2bb[0].M());
hQM->Fill(1.0,q2bb[1].M());
hQM->Fill(2.0,q2bb[2].M());
hQM->Fill(3.0,q2fv[0].M());
hQM->Fill(4.0,q2ex[0].M());
hQM->Fill(5.0,q2ex[1].M());
hQM->Fill(6.0,q2ex[2].M());
hQM->Fill(7.0,q2ex[3].M());
hQM->Fill(8.0,q2ex[4].M());
hQM->Fill(9.0,q2ex[5].M());
hQM->Fill(10.,q2ex[6].M());
hQM->Fill(11.,q2ex[7].M());
//hQM->Fill(12.,q2ex[8].M());
hQM->Fill(13.,q2ca[0].M());

hQred0->Fill( 0.0, q2bb[0].Reduced() );
hQred0->Fill( 1.0, q2bb[1].Reduced() );
hQred0->Fill( 2.0, q2bb[2].Reduced() );
hQred0->Fill( 3.0, q2fv[0].Reduced() );
hQred0->Fill( 4.0, q2ex[0].Reduced() );
hQred0->Fill( 5.0, q2ex[1].Reduced() );
hQred0->Fill( 6.0, q2ex[2].Reduced() );
hQred0->Fill( 7.0, q2ex[3].Reduced() );
hQred0->Fill( 8.0, q2ex[4].Reduced() );
hQred0->Fill( 9.0, q2ex[5].Reduced() );
hQred0->Fill( 10., q2ex[6].Reduced() );
hQred0->Fill( 11., q2ex[7].Reduced() );
// Qred 12 not stored yet
hQred0->Fill( 13., q2ca[0].Reduced() );

hPSI0->Fill( 0.0, q2bb[0].Psi2Pi() );
hPSI0->Fill( 1.0, q2bb[1].Psi2Pi() );
hPSI0->Fill( 2.0, q2bb[2].Psi2Pi() );
hPSI0->Fill( 3.0, q2fv[0].Psi2Pi() );
hPSI0->Fill( 4.0, q2ex[0].Psi2Pi() );
hPSI0->Fill( 5.0, q2ex[1].Psi2Pi() );
hPSI0->Fill( 6.0, q2ex[2].Psi2Pi() );
hPSI0->Fill( 7.0, q2ex[3].Psi2Pi() );
hPSI0->Fill( 8.0, q2ex[4].Psi2Pi() );
hPSI0->Fill( 9.0, q2ex[5].Psi2Pi() );
hPSI0->Fill( 10., q2ex[6].Psi2Pi() );
hPSI0->Fill( 11., q2ex[7].Psi2Pi() );
// PSI 12 not stored yet
hPSI0->Fill( 13., q2ca[0].Psi2Pi() );

// BBAC BBCBC BBAB                                                                                                                                                            
// FVBB FVCA                                                                                                                                                                  
// EX0BB EX0CA EX0FV (x9)                                                                                                                                                     
// = 32
hEVC[0][ 0]->Fill( q2bb[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][ 1]->Fill( q2bb[1].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][ 2]->Fill( q2bb[0].Psi2Pi(), q2bb[1].Psi2Pi() );
//-
hEVC[0][ 3]->Fill( q2fv[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][ 4]->Fill( q2fv[0].Psi2Pi(), q2ca[0].Psi2Pi() );
//-
hEVC[0][ 5]->Fill( q2ex[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][ 6]->Fill( q2ex[1].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][ 7]->Fill( q2ex[2].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][ 8]->Fill( q2ex[3].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][ 9]->Fill( q2ex[4].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][10]->Fill( q2ex[5].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][11]->Fill( q2ex[6].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][12]->Fill( q2ex[7].Psi2Pi(), q2fv[0].Psi2Pi() );
//hEVC[0][13]->Fill( q2ex[8].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][14]->Fill( q2ex[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][15]->Fill( q2ex[1].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][16]->Fill( q2ex[2].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][17]->Fill( q2ex[3].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][18]->Fill( q2ex[4].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][19]->Fill( q2ex[5].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][20]->Fill( q2ex[6].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][21]->Fill( q2ex[7].Psi2Pi(), q2fv[0].Psi2Pi() );
//hEVC[0][22]->Fill( q2ex[8].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][23]->Fill( q2ex[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][24]->Fill( q2ex[1].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][25]->Fill( q2ex[2].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][26]->Fill( q2ex[3].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][27]->Fill( q2ex[4].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[0][28]->Fill( q2ex[5].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[0][29]->Fill( q2ex[6].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[0][30]->Fill( q2ex[7].Psi2Pi(), q2fv[0].Psi2Pi() );
//hEVC[0][31]->Fill( q2ex[8].Psi2Pi(), q2ca[0].Psi2Pi() );


