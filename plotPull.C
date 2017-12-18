{
  //MPW, vg=c/n, n vs Pull along +X
  //simulation: 6 MeV e- at center and travelling along +x
  double vg[] = {1.34,1.35,1.36,1.37,1.38,1.38486, 1.39,1.40,1.41,1.42,1.43};
  double pull_vg[] = {66.3323, 65.9211, 65.3811, 65.1651, 64.3964, 64.2264, 64.5379, 63.9455, 63.3104, 63.0291, 62.7428};
  double pull_rat = 68.13;//6MeV e-, along +x, just one value
  

  //Radius vs pull
  //simulation: 6 MeV e- at thinshell: 1 to 6 m, travelling along +x
  double radius = {1,2,3,4,5,5.5};
  double pull_r_rat[]={66.0887, 71.7865, 63.2652, 61.6347, 83.8165};
  double pull_r_rat_cor[]={11.4227, 17.5684, 8.76436, 7.37302, 29.602};//after drive cor.
  double pull_r_mpw[]={};
  double pull_r_mpw_cor[]={};//after drive cor.

  TGraph *gr = new TGraph(11,vg,pull_vg);
  gr->GetXaxis()->SetTitle("n");
  gr->GetYaxis()->SetTitle("pull along direction [mm]");
  gr->Draw("a*");


} 
