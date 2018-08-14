#include "rootheader.h"

// modified by Feiyang Zhang 2018.8.14

int l_border = 0;
int r_border = 0;

////// aut to recompute baseline without real signal but has a little glitch;
void recal_baseline(int samp_max,double &base,vector<int>& samps,double* ch1,TTree *LSDetector,TString chcut){
   int sum_no = 0;
   int base_no = 0;
   base =0;
   double sum_peak =0;
   for (int k = samp_max-10; k<samp_max+11; k++)
   {
       base += ch1[k];
       base_no ++;
       samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
   }
   for (int m = l_border+20; m< r_border-20; m++ )
   {
       sum_peak += ch1[m];
       sum_no ++;
   }

   base /=base_no;

}

// calculate area, base, rms, overshoot and amp
// ii is the event id and gg is the index of channel id
void calculate(TTree *LSDetector,
		int ii,
		int gg,
		double& area_new,
		double& area_old,
		double& base_new,
		double& base_old,
		double& rms_new,
		double& overshoot_new,
		double& amp,
		double* ch1,
		double* x_num,
		int ch_anode,
		int ch_dy9)
{

	vector<int> samps;

        LSDetector->GetEntry(0);
    	int cc = LSDetector->GetLeaf("nchs")->GetValue();
	int *chan = new int[cc];

        chan[gg] = LSDetector->GetLeaf("chId")->GetValue(gg);
        TString nsampscut = Form("nsamps_%d",chan[gg]);
        TString basecut = Form("base_%d",chan[gg]);
        TString areacut = Form("area_%d",chan[gg]);
        TString chcut = Form("ch%d",chan[gg]);
        TString x_numbchcut = Form("x_numbch%d",chan[gg]);

        double base =0;
        double sum_peak =0;
        double sum =0;

	int time[1280];

        int nsamps_1;
        double base_1;
        double area_1;
        int sum_no;
        int base_no;
        double ch1_min;
        double ch1_max;
        int samp_min;
        int samp_max;

        if(chan[gg] != ch_anode && chan[gg] != ch_dy9) goto end;   // goto just like continue in a loop, it will skip to end of the function

        for(int xx=0; xx<1280; xx++)
        {
            time[xx]=0;
            ch1[xx]=0;
            x_num[xx]=0;
        }
        LSDetector->GetEntry(ii);

        //find peak
        nsamps_1 = LSDetector->GetLeaf(nsampscut)->GetValue(0);
        //if(nsamps_1>320||nsamps_1<=0) continue;
        base_1 = LSDetector->GetLeaf(basecut)->GetValue(0);
        area_1 = LSDetector->GetLeaf(areacut)->GetValue(0);
        sum_no = 0;
        base_no = 0;
        ch1_min = 10000;
        ch1_max = 8000;
        samp_min = 0;
        samp_max = 0;


//find teak and calculate the baseline in the defined time region

/////////find -pulse for anode////////////////////////////
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
            if (samp_min==0) goto end;
            for (int k = samp_min-10; k<samp_min+11; k++)
            {
                //calculate the baseline
                if (k<samp_min-3 || k>samp_min+5)
                {
                    if (k>=0)
                    {
                        base += ch1[k];
                        base_no ++;
                        //saving the baseline samples in the vector
                        samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                    }
                }
                //integral the pulse area
                else if (k>=samp_min-3 && k<=samp_min+5 && k>=0 )
                {
                    sum_peak += ch1[k];
                    sum_no ++;
                }
            }
            base /=base_no;
////// aut to recompute baseline without real signal but has a little glitch;
            //if (base-ch1_min < 6){
            if (ch1_max-ch1_min < 15){
		    recal_baseline(samp_max,base,samps,ch1,LSDetector,chcut);
            }
////// aut to recompute baseline without real signal but has a little glitch;
        }
///////////find +pulse for dynode 9///////////////////////
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
            if (samp_max==0) goto end;
            for (int k = samp_max-10; k<samp_max+11; k++)
            {
                //calculate the baseline
                if (k<samp_max-3 || k>samp_max+9)
                {
                    if (k>=0)
                    {
                        base += ch1[k];
                        base_no ++;
                        //saving the baseline samples in the vector
                        samps.push_back(LSDetector->GetLeaf(chcut)->GetValue(k));
                    }
                }
                //integral the pulse area
                else if (k>=samp_max-3 && k<=samp_max+9 && k>=0 )
                {
                    sum_peak += ch1[k];
                    sum_no ++;
                }
            }
            base /=base_no;
////// aut to recompute baseline without real signal but has a little glitch;
            //if (base-ch1_min < 6){
            if (ch1_max-ch1_min < 6){
		recal_baseline(samp_max,base,samps,ch1,LSDetector,chcut);
            }
        }
        sum = sum_peak-sum_no*base;
        rms_new =TMath::RMS(samps.size(),&samps[0]);
        overshoot_new = TMath::MaxElement(samps.size(),&samps[0]) - base;
        samps.clear();

        amp = base-ch1_min;
        base_new = base;
        base_old = base_1;
        area_new = sum;
        area_old = area_1;

end:;
     //   cout << "the end\n";
}

int main(int argc,char **argv)
{


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
        //r_border = 215;
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

    TString dirname = Form("/home/zhangfy/huangdi");
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

    TString fout_cut;
    if(LEDDR == 1)
        fout_cut = Form("/home/zhangfy/huangdi/root/%s-area9-led.root",runno.Data());
    else if(LEDDR == 0)
        fout_cut = Form("/home/zhangfy/huangdi/root/%s-area9-dr.root",runno.Data());
    TFile *fout = new TFile(fout_cut,"recreate");

    LSDetector->GetEntry(0);
    int cc = LSDetector->GetLeaf("nchs")->GetValue();   // get the number of channels
    int *chan = new int[cc];                            // array of channel id 

///// this is about definition of TTree and its branch for all channels
//
    TTree *out_tree = new TTree("out_tree","multi read out");

// define the varible of TTree branch
    int evtNo = 0;
    double* Area_new = new double[cc];
    double* Area_old = new double[cc];
    double* Base_new = new double[cc];
    double* Base_old = new double[cc];
    double* Rms_new = new double[cc];
    double* Overshoot_new = new double[cc];
    double* Amp = new double[cc];
    double** Ch1 = new double*[cc];
    double** X_num = new double*[cc];

// define the branch of TTree
    out_tree->Branch("evtNo",&evtNo,"evtNo/I");
    for(int gg=0; gg<cc; gg++){

	chan[gg] = LSDetector->GetLeaf("chId")->GetValue(gg);
        if (chan[gg] != ch_anode && chan[gg] != ch_dy9) continue;
	Ch1[gg] = new double[1280];
	X_num[gg] = new double[1280];
        out_tree->Branch(Form("area_new%d",chan[gg]),(&Area_new[gg]),Form("area_new%d/D",chan[gg]));
        out_tree->Branch(Form("area_old%d",chan[gg]),(&Area_old[gg]),Form("area_old%d/D",chan[gg]));
        out_tree->Branch(Form("base_new%d",chan[gg]),(&Base_new[gg]),Form("base_new%d/D",chan[gg]));
        out_tree->Branch(Form("base_old%d",chan[gg]),(&Base_old[gg]),Form("base_old%d/D",chan[gg]));
        out_tree->Branch(Form("rms_%d",chan[gg]),(&Rms_new[gg]),Form("rms_%d/D",chan[gg]));
        out_tree->Branch(Form("overshoot_%d",chan[gg]),(&Overshoot_new[gg]),Form("overshoot_%d/D",chan[gg]));
        out_tree->Branch(Form("amp_%d",chan[gg]),(&Amp[gg]),Form("amp%d/D",chan[gg]));
        out_tree->Branch(Form("ch%d",chan[gg]),Ch1[gg],Form("ch%d[320]/D",chan[gg]));
        out_tree->Branch(Form("x_numch%d",chan[gg]),X_num[gg],Form("x_numch%d[320]/D",chan[gg]));
    }
///////////////////////////////////////////////////////////////////
////////////////////////////


    int Entries = LSDetector->GetEntries();

    TH1F **h;
    h = new TH1F*[cc];

    for(int gg=0;gg<cc;gg++) h[gg] = new TH1F(Form("ch%i",gg),"",100,0,0);

    for(int ii=0;ii<Entries;ii++){
	evtNo = ii;
	for(int gg=0;gg<cc;gg++){
    		calculate(LSDetector,
				ii,
				gg,
				Area_new[gg],
				Area_old[gg],
				Base_new[gg],
				Base_old[gg],
				Rms_new[gg],
				Overshoot_new[gg],
				Amp[gg],Ch1[gg],
				X_num[gg],
				ch_anode,
				ch_dy9);

		h[gg]->Fill(Area_new[gg]);
	}
	out_tree->Fill();
    }

    for(int gg=0;gg<cc;gg++) h[gg]->Write();
    out_tree->Write();

    fout->Close();
    return 0;
}

