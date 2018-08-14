#include "rootheader.h"

int main(int argc,char **argv)
{

    int l_border = 0;
    int r_border = 0;

    if (argc!=5)
    {
        cout<<"Syntax:"<<argv[0]<<" [run num] [1 or 0] [ch1] [ch2] "<<endl;
        exit(1);
    }


    TString runno = TString(argv[1]);
    int ch_anode   = atoi(argv[3]);
    int ch_dy9    = atoi(argv[4]);

//LEDDR==1 is led test,LEDDR==0 is dr test;
    int LEDDR = atoi(argv[2]);

    if(LEDDR==1)
    {
        l_border = 161;
        r_border = 215;
    }
    else if(LEDDR==0)
    {
        l_border = 86;
        r_border = 98;
    }
    else
    {
        cout<<"integration border is not set correctly!!"<<endl;
        exit(1);
    }

    TString dirname = Form("/data1/pmt_testing/data");
    TString keychar = Form("adc_run%s",runno.Data());

    int c=0;
    TString file1;

    void *dirp = gSystem->OpenDirectory(dirname);
    char *direntry;
    while ((direntry=(char*)gSystem->GetDirEntry(dirp)))
    {
        TString filename(direntry);
        if(filename.Contains(keychar))
        {
            file1 = Form("%s/%s",dirname.Data(),filename.Data());
            c=c+1;
        }
        if(c==1)
            break;
    }

    TFile *file = new TFile(file1,"read");
    TTree *LSDetector = (TTree*)file->Get("LSDetector");

//int Entries = LSDetector->GetEntries();
    vector<float> eventnum;//定义一个float的容器存eventnum

    int nsamps_1;
    int sum_no;
    int base_no;
    int samp_min;
    int samp_max;
    double ch1_min;
    double ch1_max;
    double ch1[1280];
    double x_num[1280];
    double *base,*base_1,*area_1,*sum_peak,*sum;
    base = new double[500000];
    base_1 = new double[500000];
    area_1 = new double[500000];
    sum_peak = new double[500000];
    sum = new double[500000];
    int countnozle =0;
    int time[1280];
    int nhits_new = 0;
//int nsamps_new = 0;
    double darkcount[16];
    double apcount[16];
    float trigtime = 0.0;
    trigtime = LSDetector->GetMaximum("TriggerTime");

    TString fout_cut;
    if(LEDDR == 1)
        fout_cut = Form("~/out_data/root/%s-area9-led.root",runno.Data());
    else if(LEDDR == 0)
        fout_cut = Form("~/out_data/root/%s-area9-dr.root",runno.Data());
    TFile *fout = new TFile(fout_cut,"recreate");

    Int_t evtNo = 0;
    Double_t area_new = 0.;
    Double_t area_old = 0.;
    Double_t base_new = 0.;
    Double_t base_old = 0.;
    Float_t rms_new = 0.0;
    Float_t overshoot_new = 0.0;
    Double_t amp = 0.;

//for containing the value of the baseline samples
    vector<int> samps;

    TCanvas *myc = new TCanvas("myc","myc",1200,600);
    myc->Divide(2,1);
    TPaveStats *ps1;
    gStyle->SetOptFit(1);

    LSDetector->GetEntry(0);
    const int cc = LSDetector->GetLeaf("nchs")->GetValue();
    int *chan = new int[cc];

    TTree *out_tree = new TTree("out_tree","multi read out");

///
    double anew[cc],bnew[cc],aold[cc],bold[cc],rnew[cc],onew[cc],Amp[cc],Ch1[cc][1280],X_num[cc][1280];
//
    TH1F *h[cc];
    for(int gg=0; gg<cc; gg++)
    {
        h[gg] = new TH1F(Form("h_ch%d",chan[gg]),"",100,0,0);
        chan[gg] = LSDetector->GetLeaf("chId")->GetValue(gg);
        out_tree->Branch("evtNo",&evtNo,"evtNo/I");
        out_tree->Branch(Form("area_new%d",chan[gg]),(&anew[gg]),Form("area_new%d/D",chan[gg]));
        out_tree->Branch(Form("area_old%d",chan[gg]),(&aold[gg]),Form("area_old%d/D",chan[gg]));
        out_tree->Branch(Form("base_new%d",chan[gg]),(&bnew[gg]),Form("base_new%d/D",chan[gg]));
        out_tree->Branch(Form("base_old%d",chan[gg]),(&bold[gg]),Form("base_old%d/D",chan[gg]));
        out_tree->Branch(Form("rms_%d",chan[gg]),(&rnew[gg]),Form("rms_%d/F",chan[gg]));
        out_tree->Branch(Form("overshoot_%d",chan[gg]),(&onew[gg]),Form("overshoot_%d/F",chan[gg]));
        out_tree->Branch(Form("amp_%d",chan[gg]),(&Amp[gg]),Form("amp%d/D",chan[gg]));
        out_tree->Branch(Form("ch%d",chan[gg]),Ch1[gg],Form("ch%d[320]/D",chan[gg]));
        out_tree->Branch(Form("x_numch%d",chan[gg]),X_num[gg],Form("x_numch%d[320]/D",chan[gg]));
    }

    int Entries = LSDetector->GetEntries();
    for (int ii=0; ii< Entries; ii++)
    {
        base[ii] =0;
        sum_peak[ii] =0;
        sum[ii] =0;

        base_new = 0;
        base_old = 0;
        area_new = 10000;
        area_old = 10000;
        evtNo = 0;

        //TH1F *h = new TH1F(Form("h_ch%d",chan[gg]),"",100,0,0);
        for(int gg=0; gg<cc; gg++)
        {
            chan[gg] = LSDetector->GetLeaf("chId")->GetValue(gg);
            if (chan[gg] != ch_anode && chan[gg] != ch_dy9) continue;

            //out_tree = new TTree(Form("out_tree%i",chan[gg]),"new_area");
            //myc->cd(chan[gg]+1);
            TString overshootcut = Form("overshoot_%d",chan[gg]);
            TString rmscut       = Form("rms_%d",chan[gg]);
            TString chtrigposcut = Form("chtrigpos_%d",chan[gg]);
            TString area0 = Form("area_%d",chan[gg]);
            //TString area0_cut = Form("area_%d[0]<0&&rms_%d<5&&overshoot_%d<15",chan[gg],chan[gg],chan[gg]);
            darkcount[chan[gg]] = 0;
            apcount[chan[gg]] = 0;
            countnozle = 0;

            //if the effective count is lower than 300, this channel is empty
            //if(LSDetector->Draw(area0,area0_cut)<300) continue;

            //////////////////////////////////////////////////
            TString nsampscut = Form("nsamps_%d",chan[gg]);
            TString basecut = Form("base_%d",chan[gg]);
            TString areacut = Form("area_%d",chan[gg]);
            TString chcut = Form("ch%d",chan[gg]);
            TString x_numbchcut = Form("x_numbch%d",chan[gg]);
            TString nhits_cut = Form("nhits_%d",chan[gg]);

            for(int xx=0; xx<1280; xx++)
            {
                time[xx]=0;
                ch1[xx]=0;
                x_num[xx]=0;
                X_num[gg][xx]=0;
            }
            LSDetector->GetEntry(ii);

            //if(LSDetector->GetLeaf(rmscut)->GetValue(0)>=5.0) continue;
            //if (LSDetector->GetLeaf(overshootcut)->GetValue(0)>=15.0) continue;

            //find peak
            nsamps_1 = LSDetector->GetLeaf(nsampscut)->GetValue(0);
            //if(nsamps_1>320||nsamps_1<=0) continue;
            base_1[ii] = LSDetector->GetLeaf(basecut)->GetValue(0);
            area_1[ii] = LSDetector->GetLeaf(areacut)->GetValue(0);
            sum_no = 0;
            base_no = 0;
            ch1_min = 10000;
            ch1_max = 8000;
            samp_min = 0;
            samp_max = 0;

            //find the peak and calculate the baseline in the defined time region

            /////////////find -pulse for anode////////////////////////////
            if (chan[gg] == ch_anode)
            {
                for (int j=0; j<nsamps_1; j++)
                {
                    ch1[j] = LSDetector->GetLeaf(chcut)->GetValue(j);
                    x_num[j] = LSDetector->GetLeaf(x_numbchcut)->GetValue(j);
                    //find the valley
                    if (x_num[j]>l_border&&x_num[j]<r_border&&ch1_min>ch1[j]&&ch1[j]>0)
                    {
                        ch1_min = ch1[j];
                        samp_min = j;
                    }
                    //find the peak
                    if (x_num[j]>l_border&&x_num[j]<r_border&&ch1_max<ch1[j]&&ch1[j]>0)
                    {
                        ch1_max = ch1[j];
                        //samp_max = j;
                    }
                }
                if (samp_min==0) continue;
                for (int k = samp_min-10; k<samp_min+11; k++)
                {
                    //calculate the baseline
                    if (k<samp_min-3 || k>samp_min+5)
                    {
                        if (k>=0)
                        {
                            base[ii] += ch1[k];
                            base_no ++;
                            //saving the baseline samples in the vector
                            samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                        }
                    }
                    //integral the pulse area
                    else if (k>=samp_min-3 && k<=samp_min+5 && k>=0 )
                    {
                        sum_peak[ii] += ch1[k];
                        sum_no ++;
                    }
                }
                base[ii] /=base_no;
                ////// add cut to recompute baseline without real signal but has a little glitch;
                //if (base[ii]-ch1_min < 6){
                if (ch1_max-ch1_min < 15)
                {
                    sum_no = 0;
                    base_no = 0;
                    base[ii] =0;
                    sum_peak[ii] =0;
                    for (int k = samp_min-10; k<samp_min+11; k++)
                    {
                        base[ii] += ch1[k];
                        base_no ++;
                        samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                    }
                    for (int m = l_border+20; m< r_border-20; m++ )
                    {
                        sum_peak[ii] += ch1[m];
                        sum_no ++;
                    }
                    base[ii] /=base_no;
                }
                ////// add cut to recompute baseline without real signal but has a little glitch;
            }
            ///////////////find +pulse for dynode 9///////////////////////
            else if (chan[gg] == ch_dy9)
            {
                for (int j=0; j<nsamps_1; j++)
                {
                    ch1[j] = LSDetector->GetLeaf(chcut)->GetValue(j);
                    x_num[j] = LSDetector->GetLeaf(x_numbchcut)->GetValue(j);
                    //find the valley
                    if (x_num[j]>l_border&&x_num[j]<r_border&&ch1_min>ch1[j]&&ch1[j]>0)
                    {
                        ch1_min = ch1[j];
                        //samp_min = j;
                    }
                    //find the peak
                    if (x_num[j]>l_border&&x_num[j]<r_border&&ch1_max<ch1[j]&&ch1[j]>0)
                    {
                        ch1_max = ch1[j];
                        samp_max = j;
                    }
                }
                if (samp_max==0) continue;
                for (int k = samp_max-10; k<samp_max+11; k++)
                {
                    //calculate the baseline
                    if (k<samp_max-3 || k>samp_max+3)
                    {
                        if (k>=0)
                        {
                            base[ii] += ch1[k];
                            base_no ++;
                            //saving the baseline samples in the vector
                            samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                        }
                    }
                    //integral the pulse area
                    else if (k>=samp_max-3 && k<=samp_max+3 && k>=0 )
                    {
                        sum_peak[ii] += ch1[k];
                        sum_no ++;
                    }
                }
                base[ii] /=base_no;
                ////// add cut to recompute baseline without real signal but has a little glitch;
                //if (base[ii]-ch1_min < 6){
                if (ch1_max-ch1_min < 6)
                {
                    sum_no = 0;
                    base_no = 0;
                    base[ii] =0;
                    sum_peak[ii] =0;
                    for (int k = samp_max-10; k<samp_max+11; k++)
                    {
                        base[ii] += ch1[k];
                        base_no ++;
                        samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                    }
                    for (int m = l_border+20; m< r_border-20; m++ )
                    {
                        sum_peak[ii] += ch1[m];
                        sum_no ++;
                    }

                    base[ii] /=base_no;
                }
                ////// add cut to recompute baseline without real signal but has a little glitch;
            }
            sum[ii] = sum_peak[ii]-sum_no*base[ii];
            rms_new =TMath::RMS(samps.size(),&samps[0]);
            rnew[gg] = rms_new; 
            overshoot_new = TMath::MaxElement(samps.size(),&samps[0]) - base[ii];
            onew[gg] = overshoot_new;
            samps.clear();
            //if(rms_new>=3) continue;
            //if(overshoot_new>=5) continue;
            amp = base[ii]-ch1_min;
            Amp[gg] = amp;
            evtNo = ii;
            base_new = base[ii];
            base_old = base_1[ii];
            area_new = sum[ii];
            area_old = area_1[ii];
            
            bnew[gg] = base_new;
            bold[gg] = base_old;
            anew[gg] = area_new;
            aold[gg] = area_old;

            nhits_new = LSDetector->GetLeaf(nhits_cut)->GetValue(0);
            for(int hit=0; hit<nhits_new; hit++)
                time[hit]=LSDetector->GetLeaf(Form("time_%d",chan[gg]))->GetValue(hit);
            if(area_new< -40)
            {
                darkcount[chan[gg]]++;
                if(nhits_new==2&&time[1]-time[0]>100&&time[1]-time[0]<1200)
                    apcount[chan[gg]]++;
            }
            if(area_new<-40)
                countnozle++;
            h[gg]->Fill(sum[ii]);

            //Ch1[gg] = ch1;
            //X_num[gg] = x_num;

        } //////end of the event loop//////
        out_tree->Fill();

//get the pointer to the current file
//writing the final tree header to the file
    }//////end of channel loop///////

    out_tree->Write();

    for(int gg=0; gg<cc; gg++)
    {
	if(chan[gg]!= ch_anode && chan[gg]!= ch_dy9)
		continue;
        if (chan[gg] == ch_anode)
            myc->cd(1);
        else if (chan[gg] == ch_dy9)
            myc->cd(2);

        h[gg]->Draw();
        gPad->Update();
        int l_fit = 0;
        int r_fit = 0;

        if(LEDDR == 1)
        {
            l_fit = -10000;
            r_fit = 500;
        }
        else if(LEDDR == 0)
        {
            l_fit = -200;
            r_fit = -60;
        }

        //h[gg]->Fit("gaus","Q","",l_fit,r_fit);
        //h[gg]->SetStats(1);
        //ps1 = (TPaveStats*)h[gg]->FindObject("stats");
        //ps1->SetX1NDC(0.75);
        //ps1->SetX2NDC(0.95);

        //double spe1;
        //double sigma1;
        //double spe1_error;
        //TF1 *f = h[gg]->GetFunction("gaus");
        //spe1 = f->GetParameter(1);
        //sigma1 = f->GetParameter(2);
        //spe1_error = f->GetParError(1);

        //h[gg]->Fit("gaus","Q","",spe1-2*sigma1,spe1+2*sigma1);
        //f = h[gg]->GetFunction("gaus");
        //spe1 = f->GetParameter(1);
        //sigma1 = f->GetParameter(2);
        //spe1_error = f->GetParError(1);


        //cout<<runno.Data()<<" "<<chan[gg]<<" "<<spe1<<" "<<"+/- "<<spe1_error<<" "<<sigma1<<" "<<(-1)*sigma1/spe1;
        //cout<<" "<<darkcount[chan[gg]]/trigtime<<" "<<"+/- "<<sqrt(darkcount[chan[gg]])/trigtime<<" "<<apcount[chan[gg]]<<" "<<apcount[chan[gg]]/darkcount[chan[gg]];
        //cout<<" pulse count (noZLE): "<<countnozle<<endl;
    }

//cout<<"base "<<base_new<<"sum "<<area_new<<endl;

    TString fname;
    if(LEDDR == 1)
    {
        fname = Form("~/out_data/pic/%s-area9-led.C",runno.Data());
        //f2name = Form("~/out_data/pic/%s-%i-dy9-led.C",runno.Data(),ch_dy9);
    }
    else if(LEDDR == 0)
    {
        fname = Form("~/out_data/pic/%s-anode-dr.C",runno.Data());
        //f2name = Form("~/out_data/pic/%s-%i-dy9-dr.C",runno.Data(),ch_dy9);
    }
    myc->Print(fname);

    fout->Close();

    return 0;
}

