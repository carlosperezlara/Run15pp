int findtwidrange() {
  ifstream fin("twid.dat");
  int min[8];
  int max[8];
  int id, sec, iy, iz, secant;
  secant = -1;
  for(;;) {
    fin >> id;
    if(!fin.good()) break;
    fin >> sec >> iy >> iz;
    if(secant<sec) {
      secant = sec;
      min[sec] = id;
      cout << " moving to sec " << sec;
    }
    max[sec] = id;
  }
  for(int i=0; i!=8; ++i)
    cout << "SEC " << i << "  MIN " << min[i] << " " << max[i] << endl;
  return 0;
}
