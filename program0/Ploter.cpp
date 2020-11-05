using namespace std;

void hist(vector<double> var,
          const char* c_no,  const char *hist_no,
          const char *title, const char *name, const char *units,
          double low_b,      double high_b){
  // plot the results into histograms
  TGaxis::SetMaxDigits(4);
  TCanvas *c1 = new TCanvas(c_no, "c1", 150, 10, 990, 660);
  TH1F *hist  = new TH1F(hist_no, " ", 200, low_b, high_b);
  for (uint i=0;i<var.size();i++) {hist->Fill(var[i]);}

  hist->GetXaxis()->SetTitle(title);
  hist->GetYaxis()->SetTitle(Form("Candidates/(%.2f %s)",hist->GetBinWidth(1),units));
  hist->Draw("E");
  c1->SaveAs(name);
}
