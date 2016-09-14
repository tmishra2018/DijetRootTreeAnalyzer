double RSS(TH1D*, double, double, double, double, double, double, double,
           double);
double getN(TH1D*, double, double);

void Fisher_test(int dataset)
{
    switch (dataset) {
    case 0:
        cout << "Fisher test for PF RECO" << endl;
        break;
    case 1:
        cout << "Fisher test for Calo Scouting" << endl;
        break;
    default:
        cout << "Error: Not a valid value for dataset" << endl;
        return;
    }

    double lumi = 12910.0; // Integrated luminosity [/pb]
    // Fit range [GeV]
    double minX = 1.0;
    double maxX = 14000.0;
    switch (dataset) {
    case 0:
        minX = 1058.0;
        maxX = 8152.0;
        break;
    case 1:
        minX = 453.0;
        maxX = 2037.0;
        break;
    }

    TFile *f_in;
    TH1D *h_data_1GeV;
    switch (dataset) {
    case 0:
        f_in = new TFile("rawhistV10_RECO_2016BCD_rounds45678_Spring16_25nsV6_ICHEP_12.root");
        f_in->GetObject("mjj_gev", h_data_1GeV);
        break;
    case 1:
        f_in = new TFile("data_CaloScoutingHT_Run2016BCD_NewBiasCorrectedFlat_Golden12910pb_CaloDijet2016.root");
        f_in->GetObject("h_mjj_HLTpass_HT250_1GeVbin", h_data_1GeV);
        break;
    }

    const int nMassBins = 103;
    double massBoundaries[nMassBins+1]
        = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176,
           197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606,
           649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246,
           1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132,
           2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416,
           3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253,
           5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866,
           8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179,
           11571, 11977, 12395, 12827, 13272, 13732, 14000};

    TH1D *h_data = static_cast<TH1D*>(h_data_1GeV->Rebin(nMassBins, "",
                                                         massBoundaries));

    double N_data = getN(h_data, minX, maxX);

    double parameters[2][4][5] = {
        // p0, p1, p2, p3, p4
        { // PF RECO
            {3.16766240835e-07, 0.0, 6.27895318762, 0.0, 0.0}, // f2
            {1.51335086414e-05, 9.05904632837, 5.02571873249, 0.0, 0.0}, // f3
            {3.46397990934e-06, 7.34961109006, 5.96523987147, 0.163614755034,
             0.0}, // f4
            {2.38002753645e-05, 8.87638444646, 4.18327834544, -0.457885218595,
             -0.0784871491669} // f5
        },
        { // Calo Scouting
            {3.3110576299e-06, 0.0, 5.40474642534, 0.0, 0.0}, // f2
            {8.58951893126e-05, 14.3144967628, 4.57510842619, 0.0, 0.0}, // f3
            {9.49994522633e-07, 5.84090837867, 6.83225395106, 0.299718539449,
             0.0}, // f4
            {9.52739808228e-07, 5.84779556315, 6.83125819802, 0.299773743681,
             3.51085801569e-05} // f5
        }
    };

    double RSS2 = RSS(h_data, parameters[dataset][0][0],
                      parameters[dataset][0][1], parameters[dataset][0][2],
                      parameters[dataset][0][3], parameters[dataset][0][4],
                      minX, maxX, lumi);
    double RSS3 = RSS(h_data, parameters[dataset][1][0],
                      parameters[dataset][1][1], parameters[dataset][1][2],
                      parameters[dataset][1][3], parameters[dataset][1][4],
                      minX, maxX, lumi);
    double RSS4 = RSS(h_data, parameters[dataset][2][0],
                      parameters[dataset][2][1], parameters[dataset][2][2],
                      parameters[dataset][2][3], parameters[dataset][2][4],
                      minX, maxX, lumi);
    double RSS5 = RSS(h_data, parameters[dataset][3][0],
                      parameters[dataset][3][1], parameters[dataset][3][2],
                      parameters[dataset][3][3], parameters[dataset][3][4],
                      minX, maxX, lumi);

    cout << endl;
    cout << "RSS2: " << RSS2 << endl;
    cout << "RSS3: " << RSS3 << endl;
    cout << "RSS4: " << RSS4 << endl;
    cout << "RSS5: " << RSS5 << endl;

    double F32 = (RSS2 - RSS3)/(RSS3/(N_data - 3.0));
    double F43 = (RSS3 - RSS4)/(RSS4/(N_data - 4.0));
    double F54 = (RSS4 - RSS5)/(RSS5/(N_data - 5.0));

    cout << endl;
    cout << "F32: " << F32 << endl;
    cout << "F43: " << F43 << endl;
    cout << "F54: " << F54 << endl;

    double CL32 = 1.0 - TMath::FDistI(F32, 1.0, N_data - 3.0);
    double CL43 = 1.0 - TMath::FDistI(F43, 1.0, N_data - 4.0);
    double CL54 = 1.0 - TMath::FDistI(F54, 1.0, N_data - 5.0);

    cout << endl;
    cout << "CL32: " << CL32 << endl;
    cout << "CL43: " << CL43 << endl;
    cout << "CL54: " << CL54 << endl;

    return;
}

double RSS(TH1D *hist, double p0, double p1, double p2, double p3, double p4,
           double minX, double maxX, double lumi)
{
    TF1 *func = new TF1("func",
                        "[0]*TMath::Power(1.0 - x/13000.0, [1])/TMath::Power(x/13000.0, [2] + [3]*TMath::Log(x/13000) + [4]*TMath::Log(x/13000)*TMath::Log(x/13000))",
                        minX, maxX);
    func->SetParameter(0, p0);
    func->SetParameter(1, p1);
    func->SetParameter(2, p2);
    func->SetParameter(3, p3);
    func->SetParameter(4, p4);

    double rss = 0.0;
    for (int i=1; i<=hist->GetXaxis()->GetNbins(); ++i) {
        double value = hist->GetBinContent(i);
        double lowX = hist->GetXaxis()->GetBinLowEdge(i);
        double upX = hist->GetXaxis()->GetBinUpEdge(i);
        if (lowX < minX || upX > maxX)
            continue;
        if (value == 0.0)
            continue;
        double fit = func->Integral(lowX, upX)*lumi;
        rss += (value - fit)*(value - fit);
    }

    return rss;
}

double getN(TH1D *hist, double minX, double maxX)
{
    int N = 0;
    int zeros = 0;
    for (int i=1; i<=hist->GetXaxis()->GetNbins(); ++i) {
        if (hist->GetXaxis()->GetBinLowEdge(i) < minX
            || hist->GetXaxis()->GetBinUpEdge(i) > maxX)
            continue;

        if (hist->GetBinContent(i) > 0.0)
            ++N;
        else
            ++zeros;
    }

    cout << endl;
    cout << N << " data points between " << minX << " and " << maxX
         << " GeV. " << zeros << " bins with zero in that range." << endl;

    return static_cast<double>(N);
}
