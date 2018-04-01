hPSI2->Fill( 0.0, q2bb[0].Psi2Pi() );
hPSI2->Fill( 1.0, q2bb[1].Psi2Pi() );
hPSI2->Fill( 2.0, q2bb[2].Psi2Pi() );
hPSI2->Fill( 3.0, q2fv[0].Psi2Pi() );
hPSI2->Fill( 4.0, q2ex[0].Psi2Pi() );
hPSI2->Fill( 5.0, q2ex[1].Psi2Pi() );
hPSI2->Fill( 6.0, q2ex[2].Psi2Pi() );
hPSI2->Fill( 7.0, q2ex[3].Psi2Pi() );
hPSI2->Fill( 8.0, q2ex[4].Psi2Pi() );
hPSI2->Fill( 9.0, q2ex[5].Psi2Pi() );
hPSI2->Fill( 10., q2ex[6].Psi2Pi() );
hPSI2->Fill( 11., q2ex[7].Psi2Pi() );
hPSI2->Fill( 12., q2ex[8].Psi2Pi() );
hPSI2->Fill( 13., q2ca[0].Psi2Pi() );

// BBAC BBCBC BBAB
// CAAC CABC CAAB
// FVBB FVCA
// EX0BB EX0CA EX0FV (x9)
// = 32
hEVC[2][ 0]->Fill( q2bb[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][ 1]->Fill( q2bb[1].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][ 2]->Fill( q2bb[0].Psi2Pi(), q2bb[1].Psi2Pi() );
//-
hEVC[2][ 3]->Fill( q2fv[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][ 4]->Fill( q2fv[0].Psi2Pi(), q2ca[0].Psi2Pi() );
//-
hEVC[2][ 5]->Fill( q2ex[0].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][ 6]->Fill( q2ex[0].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][ 7]->Fill( q2ex[0].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][ 8]->Fill( q2ex[1].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][ 9]->Fill( q2ex[1].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][10]->Fill( q2ex[1].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][11]->Fill( q2ex[2].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][12]->Fill( q2ex[2].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][13]->Fill( q2ex[2].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][14]->Fill( q2ex[3].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][15]->Fill( q2ex[3].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][16]->Fill( q2ex[3].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][17]->Fill( q2ex[4].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][18]->Fill( q2ex[4].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][19]->Fill( q2ex[4].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][20]->Fill( q2ex[5].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][21]->Fill( q2ex[5].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][22]->Fill( q2ex[5].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][23]->Fill( q2ex[6].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][24]->Fill( q2ex[6].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][25]->Fill( q2ex[6].Psi2Pi(), q2fv[0].Psi2Pi() );
hEVC[2][26]->Fill( q2ex[7].Psi2Pi(), q2bb[2].Psi2Pi() );
hEVC[2][27]->Fill( q2ex[7].Psi2Pi(), q2ca[0].Psi2Pi() );
hEVC[2][28]->Fill( q2ex[7].Psi2Pi(), q2fv[0].Psi2Pi() );
//hEVC[2][29]->Fill( q2ex[8].Psi2Pi(), q2bb[2].Psi2Pi() );
//hEVC[2][30]->Fill( q2ex[8].Psi2Pi(), q2ca[0].Psi2Pi() );
//hEVC[2][31]->Fill( q2ex[8].Psi2Pi(), q2fv[0].Psi2Pi() );

//==                                                                                                                                                                        
// 32 +
// EX04 EX15 EX26 EX37
// = 36
double two;
two = TMath::Cos( 2.0*(bbEP[0]-bbEP[2]) ); hRes->Fill( 0.0, two); hResD->Fill( 0.0, two );//0  BBA BBC
two = TMath::Cos( 2.0*(bbEP[1]-bbEP[2]) ); hRes->Fill( 1.0, two); hResD->Fill( 1.0, two );//1  BBB BBC
two = TMath::Cos( 2.0*(bbEP[0]-bbEP[1]) ); hRes->Fill( 2.0, two); hResD->Fill( 2.0, two );//2  BBA BBB

two = TMath::Cos( 2.0*(fvEP[0]-bbEP[2]) ); hRes->Fill( 3.0, two); hResD->Fill( 3.0, two );//3  FVT BB
two = TMath::Cos( 2.0*(fvEP[0]-caEP[2]) ); hRes->Fill( 4.0, two); hResD->Fill( 4.0, two );//4  FVT CA

two = TMath::Cos( 2.0*(exEP[0]-bbEP[2]) ); hRes->Fill( 5.0, two); hResD->Fill( 5.0, two );//5  EX0 BBC
two = TMath::Cos( 2.0*(exEP[1]-bbEP[2]) ); hRes->Fill( 6.0, two); hResD->Fill( 6.0, two );//6  EX1 BBC
two = TMath::Cos( 2.0*(exEP[2]-bbEP[2]) ); hRes->Fill( 7.0, two); hResD->Fill( 7.0, two );//7  EX2 BBC
two = TMath::Cos( 2.0*(exEP[3]-bbEP[2]) ); hRes->Fill( 8.0, two); hResD->Fill( 8.0, two );//8  EX3 BBC
two = TMath::Cos( 2.0*(exEP[4]-bbEP[2]) ); hRes->Fill( 9.0, two); hResD->Fill( 9.0, two );//9  EX4 BBC
two = TMath::Cos( 2.0*(exEP[5]-bbEP[2]) ); hRes->Fill(10.0, two); hResD->Fill(10.0, two );//10 EX5 BBC
two = TMath::Cos( 2.0*(exEP[6]-bbEP[2]) ); hRes->Fill(11.0, two); hResD->Fill(11.0, two );//11 EX6 BBC
two = TMath::Cos( 2.0*(exEP[7]-bbEP[2]) ); hRes->Fill(12.0, two); hResD->Fill(12.0, two );//12 EX7 BBC

two = TMath::Cos( 2.0*(exEP[0]-fvEP[0]) ); hRes->Fill(13.0, two); hResD->Fill(13.0, two );//13 EX0 FVT
two = TMath::Cos( 2.0*(exEP[1]-fvEP[0]) ); hRes->Fill(14.0, two); hResD->Fill(14.0, two );//14 EX1 FVT
two = TMath::Cos( 2.0*(exEP[2]-fvEP[0]) ); hRes->Fill(15.0, two); hResD->Fill(15.0, two );//15 EX2 FVT
two = TMath::Cos( 2.0*(exEP[3]-fvEP[0]) ); hRes->Fill(16.0, two); hResD->Fill(16.0, two );//16 EX3 FVT
two = TMath::Cos( 2.0*(exEP[4]-fvEP[0]) ); hRes->Fill(17.0, two); hResD->Fill(17.0, two );//17 EX4 FVT
two = TMath::Cos( 2.0*(exEP[5]-fvEP[0]) ); hRes->Fill(18.0, two); hResD->Fill(18.0, two );//18 EX5 FVT
two = TMath::Cos( 2.0*(exEP[6]-fvEP[0]) ); hRes->Fill(19.0, two); hResD->Fill(19.0, two );//19 EX6 FVT
two = TMath::Cos( 2.0*(exEP[7]-fvEP[0]) ); hRes->Fill(20.0, two); hResD->Fill(20.0, two );//20 EX7 FVT

two = TMath::Cos( 2.0*(exEP[0]-caEP[2]) ); hRes->Fill(21.0, two); hResD->Fill(21.0, two );//21 EX0 CAC
two = TMath::Cos( 2.0*(exEP[1]-caEP[2]) ); hRes->Fill(22.0, two); hResD->Fill(22.0, two );//22 EX1 CAC
two = TMath::Cos( 2.0*(exEP[2]-caEP[2]) ); hRes->Fill(23.0, two); hResD->Fill(23.0, two );//23 EX2 CAC
two = TMath::Cos( 2.0*(exEP[3]-caEP[2]) ); hRes->Fill(24.0, two); hResD->Fill(24.0, two );//24 EX3 CAC
two = TMath::Cos( 2.0*(exEP[4]-caEP[2]) ); hRes->Fill(25.0, two); hResD->Fill(25.0, two );//25 EX4 CAC
two = TMath::Cos( 2.0*(exEP[5]-caEP[2]) ); hRes->Fill(26.0, two); hResD->Fill(26.0, two );//26 EX5 CAC
two = TMath::Cos( 2.0*(exEP[6]-caEP[2]) ); hRes->Fill(27.0, two); hResD->Fill(27.0, two );//27 EX6 CAC
two = TMath::Cos( 2.0*(exEP[7]-caEP[2]) ); hRes->Fill(28.0, two); hResD->Fill(28.0, two );//28 EX7 CAC

two = TMath::Cos( 2.0*(exEP[8]-bbEP[2]) ); hRes->Fill(29.0, two); hResD->Fill(29.0, two );//29 EX8 BBC
two = TMath::Cos( 2.0*(exEP[8]-fvEP[2]) ); hRes->Fill(30.0, two); hResD->Fill(30.0, two );//30 EX8 FVT
two = TMath::Cos( 2.0*(exEP[8]-caEP[2]) ); hRes->Fill(31.0, two); hResD->Fill(31.0, two );//31 EX8 CAC

two = TMath::Cos( 2.0*(exEP[0]-exEP[4]) ); hRes->Fill(32.0, two); hResD->Fill(32.0, two );//32 EX0 EX4
two = TMath::Cos( 2.0*(exEP[1]-exEP[5]) ); hRes->Fill(33.0, two); hResD->Fill(33.0, two );//33 EX1 EX5
two = TMath::Cos( 2.0*(exEP[2]-exEP[6]) ); hRes->Fill(34.0, two); hResD->Fill(34.0, two );//34 EX2 EX6
two = TMath::Cos( 2.0*(exEP[3]-exEP[7]) ); hRes->Fill(35.0, two); hResD->Fill(35.0, two );//35 EX3 EX7
