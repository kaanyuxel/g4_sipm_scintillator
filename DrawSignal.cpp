void DrawSignal(){
   TFile *f    = new TFile("time.root");
   TTree *tree = (TTree *) f->Get("test");
   TH2F *hist_signal_1  = new TH2F("hist_signal_1","Original Distribution;Time[ns];Voltage[mV]",250, 0, 500, 350, 0, 900);
   TH2F *hist_signal_2  = new TH2F("hist_signal_2","Original Distribution;Time[ns];Voltage[mV]",250, 0, 500, 350, 0, 900);
   TH1F *hist_signal_1_new  = new TH1F("hist_signal_1_new",";Time[ns];Voltage[mV]",250, 0, 500);
   TH1F *hist_signal_2_new  = new TH1F("hist_signal_2_new",";Time[ns];Voltage[mV]",250, 0, 500);

   tree->Draw("(voltage*1e9):time>>hist_signal_1","plane == 1 && eventID == 1");
   tree->Draw("(voltage*1e9):time>>hist_signal_2","plane == 2 && eventID == 1");


   for(int xbin=0; xbin<500; xbin++){
      for(int ybin=0; ybin<350; ybin++){
         double nBinCont = hist_signal_1->GetBinContent(xbin, ybin);
         if(nBinCont > 0)
           hist_signal_1_new->SetBinContent(xbin, hist_signal_1->GetYaxis()->GetBinLowEdge(ybin));
         nBinCont = hist_signal_2->GetBinContent(xbin, ybin);
         if(nBinCont > 0)
           hist_signal_2_new->SetBinContent(xbin, hist_signal_2->GetYaxis()->GetBinLowEdge(ybin));        
      }
   }

  hist_signal_1->SetMarkerColor(kBlue);
  hist_signal_2->SetMarkerColor(kRed);
  hist_signal_1->SetMarkerStyle(8);
  hist_signal_2->SetMarkerStyle(8); 
  hist_signal_1->SetMarkerSize(0.35);
  hist_signal_2->SetMarkerSize(0.35);

  hist_signal_1_new->SetLineColor(kBlue);
  hist_signal_2_new->SetLineColor(kRed);  
  hist_signal_1_new->SetLineWidth(2);
  hist_signal_2_new->SetLineWidth(2);  

  gStyle->SetOptStat(0000); 
  TCanvas *c1 = new TCanvas("c1","Summary",600,600);
  c1->cd(1);
  c1->cd(1)->SetLeftMargin(0.1320195);
  auto leg1 = new TLegend(0.635,0.76,0.86,0.88);
  leg1->SetNColumns(1);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hist_signal_1_new, "1^{st} Plane", "l");  
  leg1->AddEntry(hist_signal_2_new, "2^{nd} Plane", "l");  
  THStack *ths1 = new THStack("ths1","Histogram Distribution;Time[ns];Voltage[mV]"); 
  ths1->Add(hist_signal_1_new,"HIST sames");
  ths1->Add(hist_signal_2_new,"HIST sames"); 
  ths1->Draw("nostack");
  leg1->Draw("same");

  TCanvas *c2 = new TCanvas("c2","Summary",600,600);
  c2->cd(1);
  c2->cd(1)->SetLeftMargin(0.1320195);
  auto leg2 = new TLegend(0.635,0.76,0.86,0.88);
  leg2->SetNColumns(1);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hist_signal_1, "1^{st} Plane", "p");  
  leg2->AddEntry(hist_signal_2, "2^{nd} Plane", "p");  
  hist_signal_1->Draw("same");
  hist_signal_2->Draw("same");
  leg2->Draw("same");  
}