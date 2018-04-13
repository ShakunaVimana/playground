{
  TFile *ff1 = new TFile("NMRdata1.root");
  TFile *ff2 = new TFile("NMRdata2.root");
  TFile *ff3 = new TFile("NMRdata3.root");
  TFile *ff4 = new TFile("NMRdata4.root");
  TFile *ff5 = new TFile("NMRdata5.root");
  TFile *ff6 = new TFile("NMRdata6.root");
  TFile *ff7 = new TFile("NMRdata7.root");
  TFile *ff8 = new TFile("NMRdata8.root");

  TH1F *h1 = (TH1F*)ff1->Get("hspec_1");
  TH1F *h2 = (TH1F*)ff2->Get("hspec_1");
  TH1F *h3 = (TH1F*)ff3->Get("hspec_1");
  TH1F *h4 = (TH1F*)ff4->Get("hspec_1");
  TH1F *h5 = (TH1F*)ff5->Get("hspec_1");
  TH1F *h6 = (TH1F*)ff6->Get("hspec_1");
  TH1F *h7 = (TH1F*)ff7->Get("hspec_1");
  TH1F *h8 = (TH1F*)ff8->Get("hspec_1");

  h1->GetXaxis()->SetTitle("[ppm]*-1");h1->GetYaxis()->SetTitle("Intensity");
  h2->GetXaxis()->SetTitle("[ppm]*-1");h2->GetYaxis()->SetTitle("Intensity");
  h3->GetXaxis()->SetTitle("[ppm]*-1");h3->GetYaxis()->SetTitle("Intensity");
  h4->GetXaxis()->SetTitle("[ppm]*-1");h4->GetYaxis()->SetTitle("Intensity");
  h5->GetXaxis()->SetTitle("[ppm]*-1");h5->GetYaxis()->SetTitle("Intensity");
  h6->GetXaxis()->SetTitle("[ppm]*-1");h6->GetYaxis()->SetTitle("Intensity");
  h7->GetXaxis()->SetTitle("[ppm]*-1");h7->GetYaxis()->SetTitle("Intensity");
  h8->GetXaxis()->SetTitle("[ppm]*-1");h8->GetYaxis()->SetTitle("Intensity");
  
  h1->SetTitle("#1"); 
  h2->SetTitle("#2"); 
  h3->SetTitle("#3"); 
  h4->SetTitle("#4"); 
  h5->SetTitle("#5"); 
  h6->SetTitle("#6"); 
  h7->SetTitle("#7"); 
  h8->SetTitle("#8"); 

  TCanvas c("c","samples",800,600);
  c.Divide(4,2);
  c.cd(1);h1->Draw();
  c.cd(2);h2->Draw();
  c.cd(3);h3->Draw();
  c.cd(4);h4->Draw();
  c.cd(5);h5->Draw();
  c.cd(6);h6->Draw();
  c.cd(7);h7->Draw();
  c.cd(8);h8->Draw();

}
