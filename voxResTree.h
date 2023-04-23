//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan  9 14:39:33 2023 by ROOT version 6.24/06
// from TTree voxResTree/Voxel results and statistics
// found on file: TPCDebugVoxRes_1661234819633_1661235419619_421704_474416.root
//////////////////////////////////////////////////////////

// voxRes_527976.root for high rate
// voxRes_520259.root for low (22f)

#ifndef voxResTree_h
#define voxResTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <iomanip>
using namespace std;
#include "functions.h"
#include "TGComboBox.h"
#include "TRootEmbeddedCanvas.h"
#include "TGSlider.h"
#include "TGNumberEntry.h"
#include "TGButton.h"
#include "TGButtonGroup.h"
#include "TGStatusBar.h"
#include "TGLabel.h"
#include "TGTextEntry.h"

#include "GPU/TPCFastTransform.h"
#include "TPCReconstruction/TPCFastTransformHelperO2.h"
#include "TPCSpaceCharge/SpaceCharge.h"
#include "SpacePoints/TrackResiduals.h"
#include "CommonUtils/TreeStreamRedirector.h"

#include "Riostream.h"

using namespace o2::tpc;
using namespace o2::gpu;

// Function to create Gaussian filter
void FilterCreation(double GKernel[][5])
{
    // initialising standard deviation to 1.0
    double sigma = 1.0;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
 
    // generating 5x5 kernel
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + 2][y + 2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            GKernel[i][j] /= sum;
}

// Function to create Gaussian filter
void vec_FilterCreation(vector<vector<vector<Double_t>>>& vec_GKernel, Int_t Delta_X, Int_t Delta_Y, Int_t Delta_Z, Double_t sigma)
{
    // initialising standard deviation to 1.0
    //double sigma = 1.0;
    double r, s = 2.0 * sigma * sigma;
 
    // sum is for normalization
    double sum = 0.0;
 
    // generating 5x5 kernel
    for(int x = -Delta_X; x <= Delta_X; x++)
    {
        for(int y = -Delta_Y; y <= Delta_Y; y++)
        {
            for(int z = -Delta_Z; z <= Delta_Z; z++)
            {
                r = sqrt(x * x + y * y + z * z);
                vec_GKernel[x + Delta_X][y + Delta_Y][z + Delta_Z] = (exp(-(r * r) / s)) / (M_PI * s);
                sum += vec_GKernel[x + Delta_X][y + Delta_Y][z + Delta_Z];
            }
        }
    }
 
    // normalising the Kernel
    for(int i = 0; i < (Delta_X*2+1); ++i)
    {
        for(int j = 0; j < (Delta_Y*2+1); ++j)
        {
            for(int k = 0; k < (Delta_Z*2+1); ++k)
            {
                vec_GKernel[i][j][k] /= sum;
                //printf("internal X/Y: {%d, %d}, value: %4.3f \n",i,j,vec_GKernel[i][j]);
            }
        }
    }
}


// Header file for the classes stored in the TTree if any.

class voxResTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //o2::tpc::TrackResiduals::VoxRes *voxRes;
   Float_t         D[4];
   Float_t         E[4];
   Float_t         DS[4];
   Float_t         DC[4];
   Float_t         EXYCorr;
   Float_t         dYSigMAD;
   Float_t         dZSigLTM;
   Float_t         stat[4];  // z/x, y/x, x, entries
   UChar_t         bvox[3];
   UChar_t         bsec;
   UChar_t         flags;

   // List of branches
   TBranch        *b_voxRes_D;   //!
   TBranch        *b_voxRes_E;   //!
   TBranch        *b_voxRes_DS;   //!
   TBranch        *b_voxRes_DC;   //!
   TBranch        *b_voxRes_EXYCorr;   //!
   TBranch        *b_voxRes_dYSigMAD;   //!
   TBranch        *b_voxRes_dZSigLTM;   //!
   TBranch        *b_voxRes_stat;   //!
   TBranch        *b_voxRes_bvox;   //!
   TBranch        *b_voxRes_bsec;   //!
   TBranch        *b_voxRes_flags;   //!

   TFile* inputfile;

   vector< vector< vector< vector<TH2D*> > > > vec_h2D_DY_X_vs_Z; // sector (0..8), phi (0..14), xyz
   vector< vector< vector< vector<TH2D*> > > > vec_h2D_DSY_X_vs_Z; // sector (0..8), phi (0..14), xyz
   vector< vector< vector< vector<TH2D*> > > > vec_h2D_stat_X_vs_Z; // sector (0..8), phi (0..14), xyz

   vector< vector< vector<TH2D*> > > vec_h2D_DY_Y_vs_X; // sector (0..8), phi (0..14), xyz
   vector< vector< vector<TH2D*> > > vec_h2D_DSY_Y_vs_X; // sector (0..8), phi (0..14), xyz
   vector< vector< vector<TH2D*> > > vec_h2D_stat_Y_vs_X; // sector (0..8), phi (0..14), xyz


   vector<TH1D*> vec_h_Distortions;
   vector<TH2D*> vec_h2D_Distortions_vs_voxZ;

   TH2D *h2D_left, *h2D_right, *h2D_right_copy;

   TH2D *h2D_DX_vs_sector, *h2D_DY_vs_DX, *h2D_DZ_vs_DX, *h2D_DY_vs_sector, *h2D_DZ_vs_sector, *h2D_DY_X_vs_Z, *h2D_Y_vs_X_TPC_sector, *h2D_DZ_vs_Z, *h2D_DZ_vs_Z_trunc;
   TH2D *h2D_DX_vs_stat, *h2D_DY_vs_stat, *h2D_DXS_vs_radius;
   TH2D *h2D_DX_vs_radius;
   vector< vector<TH2D*> > vec_h2D_DZ_vs_Z;
   TProfile *TP_DZ_vs_Z, *TP_DZ_vs_Z_trunc;
   vector< vector<TProfile*> >  vec_TP_DZ_vs_Z;
   TProfile *TP_DX_vs_R, *TP_DY_vs_R;
   vector< vector<TProfile*> > vec_TP_DZ_vs_DX_tanTheta;
   vector< vector<TProfile*> > vec_TP_DY_vs_DX_tanTheta;

   TCanvas *can_TP_DX_vs_R, *can_TP_DY_vs_R, *can_vec_TP_DZ_vs_DX_tanTheta, *can_vec_TP_DY_vs_DX_tanTheta;
   TCanvas* can_vec_h_Distortions, *can_vec_h2D_Distortions_vs_voxZ, *can_h2D_DX_vs_sector, *can_h2D_DY_vs_DX, *can_h2D_DZ_vs_DX, *can_h2D_DY_vs_sector, *can_h2D_DZ_vs_sector, *can_h2D_DZ_vs_Z;
   TCanvas *can_h2D_DX_vs_stat, *can_h2D_DY_vs_stat, *can_h2D_DXS_vs_radius;
   TCanvas *can_h2D_and_h1D_DX_vs_radius;
   TString arr_label_xyz[3] = {"#DeltaX (cm)","#DeltaY (cm)","#DeltaZ (cm)"};

   TFile* outputfile;
   TTree output_tree;
   Float_t p0negZ[4], p1negZ[4], p0posZ[4], p1posZ[4], start_TF, end_TF;
   ULong64_t start_time, end_time;


   vector<TString> vec_TS_voxfiles;

   ULong_t bcolor, ycolor, gcolor, wcolor;

   TGGroupFrame* GRF_lower;
   TGVerticalFrame*   TGV_lower;
   TGHorizontalFrame* TGH_lowerA;
   TGHorizontalFrame* TGH_lowerAa;
   TGHorizontalFrame* TGH_lowerB;
   TGComboBox *fCombo_positions[3], *fCombo_file[3];
   //TGComboBox *fCombo_positionsa, *fCombo_filea;
   TGMainFrame* Frame_Setup;
   vector<TGHorizontalFrame*> vec_TGH_general;
   vector<TGHorizontalFrame*> vec_TGH_lower_split;
   TRootEmbeddedCanvas *emb_can_h2D_DY_X_vs_Z, *emb_can_h2D_DY_Y_vs_X;
   TCanvas *can_h2D_DY_X_vs_Z, *can_h2D_DY_Y_vs_X;
   TGGroupFrame *GR_bottom_sliders_sector, *GR_bottom_sliders_phi, *GR_bottom_sliders_zbin, *GR_bottom_range, *GR_Exit, *GR_delta, *GR_bottom_GF, *GR_data_selection_A, *GR_data_selection_B, *GR_data_output, *GR_data_invert, *GR_ratio, *GR_Select;
   TGHorizontalFrame *TGH_slider_sector, *TGH_slider_phi, *TGH_slider_zbin, *TGH_DeltaX_GF, *TGH_DeltaY_GF, *TGH_DeltaZ_GF, *TGH_sigma_GF;
   TGLabel *TGL_DeltaX_GF, *TGL_DeltaY_GF, *TGL_DeltaZ_GF, *TGL_sigma_GF;
   TGVerticalFrame *TGV_range, *TGV_Exit, *TGV_delta, *TGV_Select, *TGV_GF, *TGV_ExSel_master, *TGV_DSel_master, *TGV_data_selection_A, *TGV_data_selection_B, *TGV_data_selection_master, *TGV_data_output, *TGV_data_invert, *TGV_ratio, *TGV_data_output_master, *TGV_data_invert_master;
   TGHSlider *slider_sector, *slider_phi, *slider_zbin;
   TGNumberEntry *TGNum_sector, *TGNum_phi, *TGNum_zbin, *TGNum_zmin, *TGNum_zmax, *TGNum_DeltaX_GF, *TGNum_DeltaY_GF, *TGNum_DeltaZ_GF, *TGNum_sigma_GF;
   TGTextButton  *Button_exit;
   TGTextButton  *Button_export;
   TGTextButton  *Button_applyGF;
   TGCheckButton* CheckBox_scanData;
   TGCheckButton* CheckBox_getRatio;
   TGHorizontalFrame *TGH_scale_X, *TGH_scale_Y, *TGH_scale_Z;
   TGHorizontalFrame *TGH_select_DX, *TGH_select_DY, *TGH_select_DZ;
   TGNumberEntry *TGNum_scale_X, *TGNum_scale_Y, *TGNum_scale_Z;
   TGCheckButton *CheckBox_invert_X, *CheckBox_invert_Y, *CheckBox_invert_Z;
   TGCheckButton *CheckBox_gaussfilter, *CheckBox_sectoraverage;
   TGButtonGroup *TGB_group_alignment_A, *TGB_group_alignment_B, *TGB_group_alignment_C;
   TGTextEntry* TGText_outputname;
   Int_t sector_plot = 0;
   Int_t phi_plot    = 0;
   Int_t zbin_plot    = 0;
   vector<TGRadioButton*> vec_TRB_plot_data;
   vector<TGRadioButton*> select_button;

   vector< vector<TGRadioButton*>> delta_button;

   Int_t plot_data_use = 0;
   Int_t delta_pressed[3] = {0,0,0};
   Int_t selected_file = 0;
   Int_t flag_plot_data_type = 0;
   TGStatusBar          *fStatusBar;

   TF1* func_PolyFitFunc;

   vector<TString> vec_TS_afile;

   Bool_t flag_select1 = false;
   Bool_t flag_select2 = false;
   Bool_t ratio_pressed = false;
   Bool_t draw_ratio = false;
   Bool_t just_started = true;
   Int_t global_position = 0;
      Int_t N_bins_Y_GF = 3;
   Int_t N_bins_X_GF = 3;
   Int_t N_bins_Z_GF = 1;
   Double_t sigma_GF = 1.0;

   Float_t scale_XYZ[3] = {1.0,1.0,1.0};

   Int_t save_delta[2] = {0,0};                                         // Dx:0 - dYSigMa: 8
   Int_t save_selection[2][5] = {{0,0,0,-3,3},{0,0,0,-3,3}};            // Sector, Phi, zbin, zrange_min, zrange_max 
   Int_t save_filter[2][11] = {{0,0,3,3,3,0,0,0},{0,0,3,3,3,0,0,0}};    // Gaussbox, Averagebox, DX, DY, DZ, InvBox_X, InvBox_Y, InvBox_Z
   Double_t save_invScale[2][4] = {{1.3,1.0,1.0,1.0}};                  // Sigma, InvScale: X, Y, Z

   //------------------------------------
   // output map
   TFile* outputfileMap;
   o2::tpc::TrackResiduals::VoxRes mVoxelResultsOut{};                                                                ///< the results from mVoxelResults are copied in here to be able to stream them
   o2::tpc::TrackResiduals::VoxRes* mVoxelResultsOutPtr{&mVoxelResultsOut};                                           ///< pointer to set the branch address to for the output
   std::unique_ptr<TTree> mTreeOut; ///< tree holding debug output
   vector<o2::tpc::TrackResiduals::VoxRes> vec_VoxRes;
   vector<Float_t> vec_entries_voxel;
   //------------------------------------


   voxResTree(TTree *tree=0);
   virtual ~voxResTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop(Int_t N_events, Int_t file_selected);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void     Change_slider();
   void     Change_number_entry();
   void     Update_DY_X_vs_Z();
   void     Update_data_buttons(Int_t i_button_A, Int_t i_button_B);
   void     Update_delta_buttons(Int_t button_A, Int_t button_xyz);
   void     Get_directory_list();
   void     DoExport();
   void     DoNewDataSelection(Int_t i_position);
   void     func_dummy(Int_t i_select);
   void     DoNewDataSelection();
   //void     DoNewFileSelection(Int_t i_select);
   void     DoNewFileSelection(Int_t i_file);
   void     GetRatio();
   void     InitChain();
   void     SetStatusText(const char *txt, Int_t pi);
   void     EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);

};

#endif

#ifdef voxResTree_cxx
voxResTree::voxResTree(TTree *tree) : fChain(0) 
{
    printf("constructor \n");
}

voxResTree::~voxResTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t voxResTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t voxResTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void voxResTree::SetStatusText(const char *txt, Int_t pi)
{
   // Set text in status bar.
   fStatusBar->SetText(txt,pi);
}


void voxResTree::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
    //  Writes the event status in the status bar parts

    const char *text0, *text1, *text3;
    char text2[50];
    text0 = selected->GetTitle();
    SetStatusText(text0,0);
    text1 = selected->GetName();
    SetStatusText(text1,1);
    if (event == kKeyPress)
        sprintf(text2, "%c", (char) px);
    else
        sprintf(text2, "%d,%d", px, py);
    SetStatusText(text2,2);
    text3 = selected->GetObjectInfo(px,py);
    SetStatusText(text3,3);
}

void voxResTree::InitChain()
{
    printf("InitChain \n");
    inputfile->GetObject("voxResTree",fChain);

    // Set branch addresses and branch pointers
    fCurrent = -1;
    fChain->SetMakeClass(1);

    // https://github.com/AliceO2Group/AliceO2/blob/6520fa9b5dafdf8351687e7f99a3bdedc0678df5/Detectors/TPC/calibration/SpacePoints/include/SpacePoints/TrackResiduals.h
    fChain->SetBranchAddress("D[4]", D, &b_voxRes_D); ///< values of extracted distortions
    fChain->SetBranchAddress("E[4]", E, &b_voxRes_E); ///< their errors
    fChain->SetBranchAddress("DS[4]", DS, &b_voxRes_DS); ///< smoothed residual
    fChain->SetBranchAddress("DC[4]", DC, &b_voxRes_DC); ///< Cheb parameterized residual
    fChain->SetBranchAddress("EXYCorr", &EXYCorr, &b_voxRes_EXYCorr); ///< correlation between extracted X and Y
    fChain->SetBranchAddress("dYSigMAD", &dYSigMAD, &b_voxRes_dYSigMAD); ///< MAD estimator of dY sigma (dispersion after slope removal)
    fChain->SetBranchAddress("dZSigLTM", &dZSigLTM, &b_voxRes_dZSigLTM); ///< Z sigma from unbinned LTM estimator
    fChain->SetBranchAddress("stat[4]", stat, &b_voxRes_stat); ///< statistics: averages of each voxel dimension + entries

    // VoxZ: 0..4, VoxF: 0..14, VoxX: 0..151
    fChain->SetBranchAddress("bvox[3]", bvox, &b_voxRes_bvox); ///< voxel identifier: VoxZ, VoxF, VoxX
    fChain->SetBranchAddress("bsec", &bsec, &b_voxRes_bsec); ///< sector ID (0-35)
    fChain->SetBranchAddress("flags", &flags, &b_voxRes_flags); ///< status flag
    Notify();

    printf("End initChain \n");
}

void voxResTree::Init()
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    printf("voxResTree::Init() \n");

    vec_VoxRes.resize(410400);
    vec_entries_voxel.resize(410400);
    for(Int_t ientry = 0; ientry < (Int_t)vec_VoxRes.size(); ientry++)
    {
        vec_entries_voxel[ientry] = 0.0;
    }


    SetRootGraphicStyle();

    gClient->GetColorByName("yellow", ycolor);
    gClient->GetColorByName("blue", bcolor);
    gClient->GetColorByName("green", gcolor);
    gClient->GetColorByName("white", wcolor);


    Get_directory_list();

    func_PolyFitFunc = new TF1("func_PolyFitFunc",PolyFitFunc,-300,300,6);

    vec_h_Distortions.resize(3);
    vec_h2D_Distortions_vs_voxZ.resize(3);
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        HistName = "vec_h_Distortions_";
        HistName += i_xyz;
        vec_h_Distortions[i_xyz] = new TH1D(HistName.Data(),HistName.Data(),200,-70.0,70.0);

        HistName = "vec_h2D_Distortions_vs_voxZ";
        HistName += i_xyz;
        vec_h2D_Distortions_vs_voxZ[i_xyz] = new TH2D(HistName.Data(),HistName.Data(),10,-5,5,200,-7.0,7.0);
    }



    TP_DX_vs_R              = new TProfile("TP_DX_vs_R","TP_DX_vs_R",100,-250,250);
    TP_DY_vs_R              = new TProfile("TP_DY_vs_R","TP_DY_vs_R",100,-250,250);
    h2D_DX_vs_stat          = new TH2D("h2D_DX_vs_stat","h2D_DX_vs_stat",400,0,5000,400,-10,10);
    h2D_DXS_vs_radius       = new TH2D("h2D_DXS_vs_radius","h2D_DXS_vs_radius",500,-250,250,400,-10,10);
    h2D_DX_vs_radius        = new TH2D("h2D_DX_vs_radius","h2D_DX_vs_radius",500,-250,250,400,-10,10);
    h2D_DY_vs_stat          = new TH2D("h2D_DY_vs_stat","h2D_DY_vs_stat",400,0,5000,400,-10,10);
    h2D_DX_vs_sector        = new TH2D("h2D_DX_vs_sector","h2D_DX_vs_sector",540,0,36,400,-10,10);
    h2D_DY_vs_sector        = new TH2D("h2D_DY_vs_sector","h2D_DY_vs_sector",540,0,36,400,-10,10);
    h2D_DZ_vs_sector        = new TH2D("h2D_DZ_vs_sector","h2D_DZ_vs_sector",540,0,36,400,-10,10);
    h2D_DY_vs_DX            = new TH2D("h2D_DY_vs_DX","h2D_DY_vs_DX",400,-10,10,400,-10,10);
    h2D_DZ_vs_DX            = new TH2D("h2D_DZ_vs_DX","h2D_DZ_vs_DX",400,-10,10,400,-10,10);
    h2D_DY_X_vs_Z           = new TH2D("h2D_DY_X_vs_Z","h2D_DY_X_vs_Z",500,-250,250,600,-250,250);
    h2D_Y_vs_X_TPC_sector   = new TH2D("h2D_Y_vs_X_TPC_sector","h2D_Y_vs_X_TPC_sector",200,-50,50,200,80,260);
    h2D_DZ_vs_Z             = new TH2D("h2D_DZ_vs_Z","h2D_DZ_vs_Z",200,-250,250,200,-20,20);
    TP_DZ_vs_Z              = new TProfile("TP_DZ_vs_Z","TP_DZ_vs_Z",200,-250,250);
    h2D_DZ_vs_Z_trunc       = new TH2D("h2D_DZ_vs_Z_trunc","h2D_DZ_vs_Z_trunc",200,-250,250,200,-20,20);
    //TP_DZ_vs_Z_trunc      = new TProfile("TP_DZ_vs_Z_trunc","TP_DZ_vs_Z_trunc",200,-250,250);

    vec_TP_DZ_vs_DX_tanTheta.resize(2); // z
    vec_TP_DY_vs_DX_tanTheta.resize(2); // z
    for(Int_t i_z = 0; i_z < 2; i_z++)
    {
        vec_TP_DZ_vs_DX_tanTheta[i_z].resize(5); // z/x
        vec_TP_DY_vs_DX_tanTheta[i_z].resize(5); // y/x
        for(Int_t i_zx = 0; i_zx < 5; i_zx++)
        {
            HistName = "vec_TP_DZ_vs_DX_tanTheta";
            HistName += i_z;
            HistName += "_";
            HistName += i_zx;
            vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] = new TProfile(HistName.Data(),HistName.Data(),50,-8,8);

            HistName = "vec_TP_DY_vs_DX_tanTheta";
            HistName += i_z;
            HistName += "_";
            HistName += i_zx;
            vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] = new TProfile(HistName.Data(),HistName.Data(),50,-8,8);
        }
    }

    vec_h2D_DZ_vs_Z.resize(2); // normal, truncated
    vec_TP_DZ_vs_Z.resize(2); // normal, truncated
    for(Int_t i_trunc = 0; i_trunc < 2; i_trunc++)
    {
        vec_h2D_DZ_vs_Z[i_trunc].resize(4); // all radii, inner, middle, outer
        vec_TP_DZ_vs_Z[i_trunc].resize(4); // all radii, inner, middle, outer
        for(Int_t i_radius = 0; i_radius < 4; i_radius++)
        {
            HistName = "vec_h2D_DZ_vs_Z_";
            HistName += i_trunc;
            HistName += "_";
            HistName += i_radius;
            vec_h2D_DZ_vs_Z[i_trunc][i_radius] = new TH2D(HistName.Data(),HistName.Data(),200,-250,250,200,-20,20);

            HistName = "vec_TP_DZ_vs_Z";
            HistName += i_trunc;
            HistName += "_";
            HistName += i_radius;
            vec_TP_DZ_vs_Z[i_trunc][i_radius] = new TProfile(HistName.Data(),HistName.Data(),200,-250,250);
        }
    }


    // -----------------------------------------------------------------------------------------------------
    //                          X vs. Z plot
    // -----------------------------------------------------------------------------------------------------

    vec_h2D_DY_X_vs_Z.resize(2); // File 1&2
    vec_h2D_DSY_X_vs_Z.resize(2); 
    vec_h2D_stat_X_vs_Z.resize(2);

    for(Int_t i_file = 0; i_file < 2; i_file++)
    {
        vec_h2D_DY_X_vs_Z[i_file].resize(9); // sector (0..8), phi (0..14)
        vec_h2D_DSY_X_vs_Z[i_file].resize(9); // sector (0..8), phi (0..14)
        vec_h2D_stat_X_vs_Z[i_file].resize(9); // sector (0..8), phi (0..14)

        for(Int_t i_sector = 0; i_sector < 9; i_sector++)
        {
            vec_h2D_DY_X_vs_Z[i_file][i_sector].resize(15);
            vec_h2D_DSY_X_vs_Z[i_file][i_sector].resize(15);
            vec_h2D_stat_X_vs_Z[i_file][i_sector].resize(15);
            for(Int_t i_phi = 0; i_phi < 15; i_phi++)
            {
                vec_h2D_DY_X_vs_Z[i_file][i_sector][i_phi].resize(3);
                vec_h2D_DSY_X_vs_Z[i_file][i_sector][i_phi].resize(3);
                vec_h2D_stat_X_vs_Z[i_file][i_sector][i_phi].resize(3);
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    HistName = "vec_h2D_DY_X_vs_Z_file";
                    HistName += i_file;
                    HistName += "_sec";
                    HistName += i_sector;
                    HistName += "_phi";
                    HistName += i_phi;
                    HistName += "_xyz";
                    HistName += i_xyz;
                    vec_h2D_DY_X_vs_Z[i_file][i_sector][i_phi][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,600,-250,250);

                    HistName = "vec_h2D_DSY_X_vs_Z_file";
                    HistName += i_file;
                    HistName += "_sec";
                    HistName += i_sector;
                    HistName += "_phi";
                    HistName += i_phi;
                    HistName += "_xyz";
                    HistName += i_xyz;
                    vec_h2D_DSY_X_vs_Z[i_file][i_sector][i_phi][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,600,-250,250);

                    HistName = "vec_h2D_stat_X_vs_Z_file";
                    HistName += i_file;
                    HistName += "_sec";
                    HistName += i_sector;
                    HistName += "_phi";
                    HistName += i_phi;
                    HistName += "_xyz";
                    HistName += i_xyz;
                    vec_h2D_stat_X_vs_Z[i_file][i_sector][i_phi][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,600,-250,250);
                }
            }
        }


    } 



    // -----------------------------------------------------------------------------------------------------
    //                          X vs. Y plot
    // -----------------------------------------------------------------------------------------------------

    
    vec_h2D_DY_Y_vs_X.resize(2); // File 1&2
    vec_h2D_DSY_Y_vs_X.resize(2); 
    vec_h2D_stat_Y_vs_X.resize(2); 

    for(Int_t i_file = 0; i_file < 2; i_file++)
    {
        vec_h2D_DY_Y_vs_X[i_file].resize(10); // z/x bins
        vec_h2D_DSY_Y_vs_X[i_file].resize(10); // z/x bins
        vec_h2D_stat_Y_vs_X[i_file].resize(10); // z/x bins
        for(Int_t i_z_bin = 0; i_z_bin < 10; i_z_bin++)
        {
            vec_h2D_DY_Y_vs_X[i_file][i_z_bin].resize(3);
            vec_h2D_DSY_Y_vs_X[i_file][i_z_bin].resize(3);
            vec_h2D_stat_Y_vs_X[i_file][i_z_bin].resize(3);
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                HistName = "vec_h2D_DY_Y_vs_X_file";
                HistName += i_file;
                HistName += "_sec";
                HistName += i_z_bin;
                HistName += "_xyz";
                HistName += i_xyz;
                vec_h2D_DY_Y_vs_X[i_file][i_z_bin][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,500,-250,250);

                HistName = "vec_h2D_DSY_Y_vs_X_file";
                HistName += i_file;
                HistName += "_sec";
                HistName += i_z_bin;
                HistName += "_xyz";
                HistName += i_xyz;
                vec_h2D_DSY_Y_vs_X[i_file][i_z_bin][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,500,-250,250);

                HistName = "vec_h2D_stat_Y_vs_X_file";
                HistName += i_file;
                HistName += "_sec";
                HistName += i_z_bin;
                HistName += "_xyz";
                HistName += i_xyz;
                vec_h2D_stat_Y_vs_X[i_file][i_z_bin][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),500,-250,250,500,-250,250);
            }
        }   
    }

    can_vec_h_Distortions = new TCanvas("can_vec_h_Distortions","can_vec_h_Distortions",10,10,1300,500);
    can_vec_h_Distortions ->Divide(3,1);

    can_vec_h2D_Distortions_vs_voxZ = new TCanvas("can_vec_h2D_Distortions_vs_voxZ","can_vec_h2D_Distortions_vs_voxZ",10,10,1300,500);
    can_vec_h2D_Distortions_vs_voxZ ->Divide(3,1);

    can_TP_DX_vs_R       = new TCanvas("can_TP_DX_vs_R","can_TP_DX_vs_R",10,10,800,500);
    can_TP_DY_vs_R       = new TCanvas("can_TP_DY_vs_R","can_TP_DY_vs_R",10,10,800,500);

    can_vec_TP_DZ_vs_DX_tanTheta = new TCanvas("can_vec_TP_DZ_vs_DX_tanTheta","can_vec_TP_DZ_vs_DX_tanTheta",10,10,800,500);
    can_vec_TP_DY_vs_DX_tanTheta = new TCanvas("can_vec_TP_DY_vs_DX_tanTheta","can_vec_TP_DY_vs_DX_tanTheta",10,10,800,500);

    can_h2D_DX_vs_stat   = new TCanvas("can_h2D_DX_vs_stat","can_h2D_DX_vs_stat",10,10,800,500);
    can_h2D_DXS_vs_radius  = new TCanvas("can_h2D_DXS_vs_radius","can_h2D_DXS_vs_radius",10,10,800,500);
    can_h2D_and_h1D_DX_vs_radius  = new TCanvas("can_h2D_and_h1D_DX_vs_radius","can_h2D_and_h1D_DX_vs_radius",10,10,800,500);           //NEW 1D and 2D Histogram

    can_h2D_DY_vs_stat   = new TCanvas("can_h2D_DY_vs_stat","can_h2D_DY_vs_stat",10,10,800,500);
    can_h2D_DX_vs_sector = new TCanvas("can_h2D_DX_vs_sector","can_h2D_DX_vs_sector",10,10,800,500);
    can_h2D_DY_vs_sector = new TCanvas("can_h2D_DY_vs_sector","can_h2D_DY_vs_sector",10,10,800,500);
    can_h2D_DZ_vs_sector = new TCanvas("can_h2D_DZ_vs_sector","can_h2D_DZ_vs_sector",10,10,800,500);
    can_h2D_DY_vs_DX     = new TCanvas("can_h2D_DY_vs_DX","can_h2D_DY_vs_DX",10,10,800,500);
    can_h2D_DZ_vs_DX     = new TCanvas("can_h2D_DZ_vs_DX","can_h2D_DZ_vs_DX",10,10,800,500);

    can_h2D_DZ_vs_Z = new TCanvas("can_h2D_DZ_vs_Z","can_h2D_DZ_vs_Z",10,10,1250,500);
    can_h2D_DZ_vs_Z ->Divide(2,1);
    for(Int_t i_pad = 0; i_pad < 2; i_pad++)
    {
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetFillColor(10);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetTopMargin(0.04);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetBottomMargin(0.15);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetRightMargin(0.15);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetLeftMargin(0.15);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetTicks(1,1);
        can_h2D_DZ_vs_Z->cd(i_pad+1) ->SetGrid(0,0);
    }

    //------------------------------------------------------------------------
    printf("Start defining GUI \n");
    // GUI
    // https://root.cern.ch/root/htmldoc/guides/users-guide/WritingGUI.html

    Frame_Setup = new TGMainFrame(gClient->GetRoot(), 400, 100);
    Frame_Setup ->SetWindowName("Setup");

    vec_TGH_general.resize(2); // top, down
    vec_TGH_lower_split.resize(2); // down split into 2
    vec_TGH_general[0]            = new TGHorizontalFrame(Frame_Setup);
    vec_TGH_general[1]            = new TGHorizontalFrame(Frame_Setup);
    vec_TGH_lower_split[0]        = new TGHorizontalFrame(vec_TGH_general[1]);
    vec_TGH_lower_split[1]        = new TGHorizontalFrame(vec_TGH_general[1]);
    Frame_Setup ->AddFrame(vec_TGH_general[0], new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,5, 5, 5, 5));
    Frame_Setup ->AddFrame(vec_TGH_general[1], new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,5, 5, 5, 5));
    //vec_TGH_general[1] ->AddFrame(vec_TGH_lower_split[0], new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,5, 5, 5, 5));
    //vec_TGH_general[1] ->AddFrame(vec_TGH_lower_split[1], new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,5, 5, 5, 5));


    GRF_lower   = new TGGroupFrame(vec_TGH_general[1],"Input data",kHorizontalFrame);
    TGV_lower   = new TGVerticalFrame(GRF_lower);
    GRF_lower   ->AddFrame(TGV_lower, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGH_lowerA  = new TGHorizontalFrame(TGV_lower,200,100);
    TGH_lowerAa = new TGHorizontalFrame(TGV_lower,200,100);
    TGH_lowerB  = new TGHorizontalFrame(TGV_lower,200,100);
    TGV_lower   ->AddFrame(TGH_lowerA, new TGLayoutHints(kLHintsLeft,5,5,3,4)); // kLHintsCenterX
    TGV_lower   ->AddFrame(TGH_lowerAa, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_lower   ->AddFrame(TGH_lowerB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    vec_TGH_general[1] ->AddFrame(GRF_lower, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY,5, 5, 5, 5));


    emb_can_h2D_DY_X_vs_Z = new TRootEmbeddedCanvas("emb_can_h2D_DY_X_vs_Z",vec_TGH_general[0],200,400);
    can_h2D_DY_X_vs_Z = emb_can_h2D_DY_X_vs_Z->GetCanvas();
    can_h2D_DY_X_vs_Z->cd();

    can_h2D_DY_X_vs_Z ->SetFillColor(10);
    can_h2D_DY_X_vs_Z ->SetTopMargin(0.04);
    can_h2D_DY_X_vs_Z ->SetBottomMargin(0.18);
    can_h2D_DY_X_vs_Z ->SetRightMargin(0.2);
    can_h2D_DY_X_vs_Z ->SetLeftMargin(0.2);
    can_h2D_DY_X_vs_Z ->SetTicks(1,1);
    can_h2D_DY_X_vs_Z ->SetGrid(0,0);

    emb_can_h2D_DY_Y_vs_X = new TRootEmbeddedCanvas("emb_can_h2D_DY_Y_vs_X",vec_TGH_general[0],200,400);
    //Int_t wid = emb_can_h2D_DY_Y_vs_X->GetCanvasWindowId();
    //TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
    //emb_can_h2D_DY_Y_vs_X->AdoptCanvas(myc);
    //myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,"EventInfo(Int_t,Int_t,Int_t,TObject*)");

    can_h2D_DY_Y_vs_X = emb_can_h2D_DY_Y_vs_X->GetCanvas();
    can_h2D_DY_Y_vs_X->cd();

    can_h2D_DY_Y_vs_X ->SetFillColor(10);
    can_h2D_DY_Y_vs_X ->SetTopMargin(0.05);
    can_h2D_DY_Y_vs_X ->SetBottomMargin(0.2);
    can_h2D_DY_Y_vs_X ->SetRightMargin(0.28); // 0.22
    can_h2D_DY_Y_vs_X ->SetLeftMargin(0.14); // 0.2
    can_h2D_DY_Y_vs_X ->SetTicks(1,1);
    can_h2D_DY_Y_vs_X ->SetGrid(0,0);


    //-------------------
    TGB_group_alignment_A = new TGButtonGroup(TGH_lowerB,"Deltas");
    TGB_group_alignment_B = new TGButtonGroup(TGH_lowerB,"Smoothed");
    TGB_group_alignment_C = new TGButtonGroup(TGH_lowerB,"Statistics");

    // control horizontal position of the text
    TGB_group_alignment_A ->SetTitlePos(TGGroupFrame::kCenter);
    TGB_group_alignment_B ->SetTitlePos(TGGroupFrame::kCenter);
    TGB_group_alignment_C ->SetTitlePos(TGGroupFrame::kCenter);

    vec_TRB_plot_data.resize(9);
    vec_TRB_plot_data[0] = new TGRadioButton(TGB_group_alignment_A, "Dx", kTextRight);
    vec_TRB_plot_data[1] = new TGRadioButton(TGB_group_alignment_A, "Dy", kTextRight);
    vec_TRB_plot_data[2] = new TGRadioButton(TGB_group_alignment_A, "Dz", kTextRight);
    vec_TRB_plot_data[3] = new TGRadioButton(TGB_group_alignment_B, "DSx", kTextRight);
    vec_TRB_plot_data[4] = new TGRadioButton(TGB_group_alignment_B, "DSy", kTextRight);
    vec_TRB_plot_data[5] = new TGRadioButton(TGB_group_alignment_B, "DSz", kTextRight);
    vec_TRB_plot_data[6] = new TGRadioButton(TGB_group_alignment_C, "stat.", kTextRight);
    vec_TRB_plot_data[7] = new TGRadioButton(TGB_group_alignment_C, "flag", kTextRight);
    vec_TRB_plot_data[8] = new TGRadioButton(TGB_group_alignment_C, "dYSigMAD", kTextRight);

    for(Int_t i_button = 0; i_button < (Int_t)vec_TRB_plot_data.size(); i_button++)
    {
        vec_TRB_plot_data[i_button] ->Connect("Pressed()", "voxResTree", this, Form("Update_data_buttons(Int_t=%d, %d)", i_button,0));
        printf("Button Initiated: %d\n",i_button);
    }

    TGB_group_alignment_A ->SetButton(kTextCenterX);
    TGB_group_alignment_B ->SetButton(kTextCenterX);
    TGB_group_alignment_C ->SetButton(kTextCenterX);
    TGH_lowerB ->AddFrame(TGB_group_alignment_A, new TGLayoutHints(kLHintsExpandX));
    TGH_lowerB ->AddFrame(TGB_group_alignment_B, new TGLayoutHints(kLHintsExpandX));
    TGH_lowerB ->AddFrame(TGB_group_alignment_C, new TGLayoutHints(kLHintsExpandX));

    //-------------------


    //-------------------
    vec_TGH_general[0] ->AddFrame(emb_can_h2D_DY_X_vs_Z, new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 10, 10, 4, 4));
    vec_TGH_general[0] ->AddFrame(emb_can_h2D_DY_Y_vs_X, new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 10, 10, 4, 4));

    // status bar
    //Int_t parts[] = {45, 15, 10, 30};
    //fStatusBar = new TGStatusBar(vec_TGH_general[0], 50, 10, kVerticalFrame);
    //fStatusBar->SetParts(parts, 4);
    //fStatusBar->Draw3DCorner(kFALSE);
    //vec_TGH_general[0] ->AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
    //-------------------



    // Lower split top frame
    TGV_ExSel_master          = new TGVerticalFrame(TGH_lowerA);
    TGV_data_selection_master = new TGVerticalFrame(TGH_lowerA);
    TGV_DSel_master           = new TGVerticalFrame(TGH_lowerA);
    TGV_data_output_master    = new TGVerticalFrame(TGH_lowerA);
    TGV_data_invert_master    = new TGVerticalFrame(TGH_lowerA);


    //--------------------------------------
    GR_data_selection_A  = new TGGroupFrame(TGV_data_selection_master,"Data selection A",kVerticalFrame);
    TGV_data_selection_A = new TGVerticalFrame(GR_data_selection_A);
    GR_data_selection_A  ->AddFrame(TGV_data_selection_A, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    // data run selection A
    fCombo_positions[0] = new TGComboBox(TGV_data_selection_A, 88);
    TGV_data_selection_A->AddFrame(fCombo_positions[0], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    for(Int_t i_pos_entry = 0; i_pos_entry < (Int_t)vec_TS_afile.size(); i_pos_entry++)
    {
        HistName = "./Data/";
        HistName += vec_TS_afile[i_pos_entry];
        HistName += "/";
        fCombo_positions[0]->AddEntry(HistName.Data(),i_pos_entry);
    }
    fCombo_positions[0]->Resize(520, 20);
    fCombo_positions[0]->Select(0);
    fCombo_positions[0]->SetBackgroundColor(wcolor);
    fCombo_positions[0]->Connect("Changed()", "voxResTree", this,Form("DoNewDataSelection(Int_t=%d)",0));


    // file selection
    fCombo_file[0] = new TGComboBox(TGV_data_selection_A, 88);
    TGV_data_selection_A->AddFrame(fCombo_file[0], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    fCombo_file[0]->Resize(520, 20);
    fCombo_file[0]->Select(0);
    fCombo_file[0]->SetBackgroundColor(wcolor);

    //fCombo_file[0]->Connect("ReturnPressed()", "voxResTree", this, "HandleButtonsGeneral(=14)");
    fCombo_file[0]->Connect("Selected(Int_t)", "voxResTree", this, Form("DoNewFileSelection(Int_t=%d)",0));
    fCombo_file[0]->Connect("Changed()", "voxResTree", this, Form("func_dummy(Int_t=%d)",5));
    //fCombo_file[0]->Connect("Changed()", "voxResTree", this, "DoNewFileSelection()");
    //fCombo_file[0]->Connect("Changed()", "voxResTree", this, Form("DoNewFileSelection(Int_t=%d)",0));
    
    // exit button
    // Button_exit = new TGTextButton(TGH_lowerA, "&Exit ","gApplication->Terminate(0)");
    // TGH_lowerA->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // check box scan data
    // CheckBox_scanData  = new TGCheckButton(TGH_lowerA, new TGHotString("Scan data"), -1);
    // CheckBox_scanData ->SetState(kButtonUp);
    // //CheckBox_scanData ->Connect("Clicked()", "voxResTree", this, "DoNewDataSelection()");
    // CheckBox_scanData ->Connect("Clicked()", "voxResTree", this,Form("DoNewDataSelection(Int_t=%d)",0));
    // TGH_lowerA    ->AddFrame(CheckBox_scanData, new TGLayoutHints(kLHintsCenterX,5,5,3,4));



    //--------------------------------------
    // data run selection B
    GR_data_selection_B  = new TGGroupFrame(TGV_data_selection_master,"Data selection B",kVerticalFrame);
    TGV_data_selection_B = new TGVerticalFrame(GR_data_selection_B);
    GR_data_selection_B  ->AddFrame(TGV_data_selection_B, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    fCombo_positions[1] = new TGComboBox(TGV_data_selection_B, 88);
    TGV_data_selection_B->AddFrame(fCombo_positions[1], new TGLayoutHints(kLHintsLeft,5,5,3,4));
    for(Int_t i_pos_entry = 0; i_pos_entry < (Int_t)vec_TS_afile.size(); i_pos_entry++)
    {
        HistName = "./Data/";
        HistName += vec_TS_afile[i_pos_entry];
        HistName += "/";
        fCombo_positions[1]->AddEntry(HistName.Data(),i_pos_entry);
    }
    fCombo_positions[1]->Resize(520, 20);
    fCombo_positions[1]->Select(0);
    fCombo_positions[1]->SetBackgroundColor(wcolor);
    fCombo_positions[1]->Connect("Changed()", "voxResTree", this, Form("DoNewDataSelection(Int_t=%d)",1));

    // file selection
    fCombo_file[1] = new TGComboBox(TGV_data_selection_B, 88);
    TGV_data_selection_B->AddFrame(fCombo_file[1], new TGLayoutHints(kLHintsLeft,5,5,3,4));
    fCombo_file[1]->Resize(520, 20);
    fCombo_file[1]->Select(0);
    fCombo_file[1]->SetBackgroundColor(wcolor);

    //fCombo_file[1]->Connect("ReturnPressed()", "voxResTree", this, "HandleButtonsGeneral(=14)");
    fCombo_file[1]->Connect("Selected(Int_t)", "voxResTree", this, Form("DoNewFileSelection(Int_t=%d)",1));
    fCombo_file[1]->Connect("Changed()", "voxResTree", this, Form("func_dummy(Int_t=%d)",7));
    //fCombo_file[1]->Connect("Changed()", "voxResTree", this, "DoNewFileSelection()");
    //fCombo_file[1]->Connect("Changed()", "voxResTree", this, Form("DoNewFileSelection(Int_t=%d)",0));



    //--------------------------------------
    // invert sign(s)
    GR_data_invert  = new TGGroupFrame(TGV_data_invert_master,"Invert signs and scale",kVerticalFrame);
    TGV_data_invert = new TGVerticalFrame(GR_data_invert);
    GR_data_invert  ->AddFrame(TGV_data_invert, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    GR_ratio  = new TGGroupFrame(TGV_data_invert_master,"Get Ratio",kVerticalFrame);
    TGV_ratio = new TGVerticalFrame(GR_ratio);
    GR_ratio  ->AddFrame(TGV_ratio, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    // // check box get ratio
    CheckBox_getRatio  = new TGCheckButton(TGV_ratio, new TGHotString("Ratio X/B"), -1);
    // CheckBox_getRatio ->SetState(kButtonUp);
    //CheckBox_scanData ->Connect("Clicked()", "voxResTree", this, "DoNewDataSelection()");
    CheckBox_getRatio ->Connect("Clicked()", "voxResTree", this,Form("GetRatio()"));
    TGV_ratio    ->AddFrame(CheckBox_getRatio, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    // check box invert sign
    TGH_scale_X = new TGHorizontalFrame(TGV_data_invert);
    CheckBox_invert_X  = new TGCheckButton(TGH_scale_X, new TGHotString("X"), -1);
    CheckBox_invert_X ->SetState(kButtonUp);
    TGNum_scale_X = new TGNumberEntry(TGH_scale_X, 1.0, 3,(TGNumberFormat::EStyle) 1);
    TGH_scale_X ->AddFrame(CheckBox_invert_X, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_scale_X ->AddFrame(TGNum_scale_X, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_invert   ->AddFrame(TGH_scale_X, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGH_scale_Y = new TGHorizontalFrame(TGV_data_invert);
    CheckBox_invert_Y  = new TGCheckButton(TGH_scale_Y, new TGHotString("Y"), -1);
    CheckBox_invert_Y ->SetState(kButtonUp);
    TGNum_scale_Y = new TGNumberEntry(TGH_scale_Y, 1.0, 3,(TGNumberFormat::EStyle) 1);
    TGH_scale_Y ->AddFrame(CheckBox_invert_Y, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_scale_Y ->AddFrame(TGNum_scale_Y, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_invert   ->AddFrame(TGH_scale_Y, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGH_scale_Z = new TGHorizontalFrame(TGV_data_invert);
    CheckBox_invert_Z  = new TGCheckButton(TGH_scale_Z, new TGHotString("Z"), -1);
    CheckBox_invert_Z ->SetState(kButtonUp);
    TGNum_scale_Z = new TGNumberEntry(TGH_scale_Z, 1.0, 3,(TGNumberFormat::EStyle) 1);
    TGH_scale_Z ->AddFrame(CheckBox_invert_Z, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_scale_Z ->AddFrame(TGNum_scale_Z, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_invert   ->AddFrame(TGH_scale_Z, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    //--------------------------------------




    //--------------------------------------
    // output
    GR_data_output  = new TGGroupFrame(TGV_data_output_master,"Data output - export",kVerticalFrame);
    TGV_data_output = new TGVerticalFrame(GR_data_output);
    GR_data_output  ->AddFrame(TGV_data_output, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    // data run output A
    fCombo_positions[2] = new TGComboBox(TGV_data_output, 88);
    TGV_data_output->AddFrame(fCombo_positions[2], new TGLayoutHints(kLHintsLeft,5,5,3,4));
    for(Int_t i_pos_entry = 0; i_pos_entry < (Int_t)vec_TS_afile.size(); i_pos_entry++)
    {
        HistName = "./Data/";
        HistName += vec_TS_afile[i_pos_entry];
        HistName += "/";
        fCombo_positions[2]->AddEntry(HistName.Data(),i_pos_entry);
    }
    fCombo_positions[2]->Resize(200, 20);
    fCombo_positions[2]->Select(0);
    fCombo_positions[2]->SetBackgroundColor(wcolor);

    // file output
    TGText_outputname = new TGTextEntry(TGV_data_output,"Map_export.root");
    TGText_outputname ->SetTextColor(kRed);
    TGText_outputname ->SetToolTipText("Output/export map name");
    TGV_data_output   ->AddFrame(TGText_outputname, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGText_outputname ->Resize(200, 20);


    //fCombo_file[2] = new TGComboBox(TGV_data_output, 88);
    //TGV_data_output->AddFrame(fCombo_file[2], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    //fCombo_file[2]->Resize(520, 20);
    //fCombo_file[2]->Select(0);
    //fCombo_file[2]->SetBackgroundColor(wcolor);

    //fCombo_file[0]->Connect("ReturnPressed()", "voxResTree", this, "HandleButtonsGeneral(=14)");
    //fCombo_file[2]->Connect("Selected(Int_t)", "voxResTree", this, "DoNewFileSelection()");
    //fCombo_file[2]->Connect("Changed()", "voxResTree", this, Form("func_dummy(Int_t=%d)",5));
    //fCombo_file[0]->Connect("Changed()", "voxResTree", this, "DoNewFileSelection()");
    //fCombo_file[0]->Connect("Changed()", "voxResTree", this, Form("DoNewFileSelection(Int_t=%d)",0));
    //--------------------------------------




    // Group Exit Export
    GR_Exit  = new TGGroupFrame(TGV_ExSel_master,"Exit Export",kVerticalFrame);
    TGV_Exit = new TGVerticalFrame(GR_Exit);
    GR_Exit  ->AddFrame(TGV_Exit, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    // exit button
    Button_exit = new TGTextButton(TGV_Exit, "&Exit ","gApplication->Terminate(0)");
    TGV_Exit->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // export button
    Button_export = new TGTextButton(TGV_Exit, "&Export");
    Button_export->Connect("Clicked()", "voxResTree", this, "DoExport()");
    Button_export->SetToolTipText("Export map modified via inverse and/or Gaussian filter");
    TGV_Exit->AddFrame(Button_export, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    // Show Buttons
    GR_Select  = new TGGroupFrame(TGV_ExSel_master,"Select",kVerticalFrame);
    TGV_Select = new TGVerticalFrame(GR_Select);
    GR_Select  ->AddFrame(TGV_Select, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    
    select_button.resize(2);
    select_button[0] = new TGRadioButton(TGV_Select, "Show A", kTextRight);
    select_button[1] = new TGRadioButton(TGV_Select, "Show B", kTextRight);
    TGV_Select->AddFrame(select_button[0], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGV_Select->AddFrame(select_button[1], new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    select_button[0] ->Connect("Pressed()", "voxResTree", this, Form("Update_data_buttons(Int_t=%d, %d)", 0,1));
    select_button[1] ->Connect("Pressed()", "voxResTree", this, Form("Update_data_buttons(Int_t=%d, %d)", 1,1));

    // check box scan data
    //CheckBox_scanData  = new TGCheckButton(TGH_lowerA, new TGHotString("Scan data"), -1);
    //CheckBox_scanData ->SetState(kButtonUp);
    //CheckBox_scanData ->Connect("Clicked()", "voxResTree", this,Form("DoNewDataSelection(Int_t=%d)",0));

    // Delta - Selection
    GR_delta  = new TGGroupFrame(TGV_DSel_master,"Delta Selection",kVerticalFrame);
    TGV_delta = new TGVerticalFrame(GR_delta);
    GR_delta  ->AddFrame(TGV_delta, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGH_select_DX = new TGHorizontalFrame(TGV_delta);
    TGH_select_DY = new TGHorizontalFrame(TGV_delta);
    TGH_select_DZ = new TGHorizontalFrame(TGV_delta);
    
    delta_button.resize(2);
    delta_button[0].resize(3);
    delta_button[1].resize(3);

    delta_button[0][0] = new TGRadioButton(TGH_select_DX, "DX", kTextRight);
    delta_button[0][1] = new TGRadioButton(TGH_select_DY, "DY", kTextRight);
    delta_button[0][2] = new TGRadioButton(TGH_select_DZ, "DZ", kTextRight);
    delta_button[1][0] = new TGRadioButton(TGH_select_DX, "DX", kTextRight);
    delta_button[1][1] = new TGRadioButton(TGH_select_DY, "DY", kTextRight);
    delta_button[1][2] = new TGRadioButton(TGH_select_DZ, "DZ", kTextRight);

    TGH_select_DX->AddFrame(delta_button[0][0], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGH_select_DX->AddFrame(delta_button[1][0], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGV_delta ->AddFrame(TGH_select_DX, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGH_select_DY->AddFrame(delta_button[0][1], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGH_select_DY->AddFrame(delta_button[1][1], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGV_delta ->AddFrame(TGH_select_DY, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGH_select_DZ->AddFrame(delta_button[0][2], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGH_select_DZ->AddFrame(delta_button[1][2], new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    TGV_delta ->AddFrame(TGH_select_DZ, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    delta_button[0][0] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 0,0));
    delta_button[0][1] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 0,1));
    delta_button[0][2] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 0,2));
    delta_button[1][0] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 1,0));
    delta_button[1][1] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 1,1));
    delta_button[1][2] ->Connect("Pressed()", "voxResTree", this, Form("Update_delta_buttons(Int_t=%d, %d)", 1,2));





    // Group Gaussian filter
    GR_bottom_GF  = new TGGroupFrame(TGH_lowerA,"Gauss filter",kVerticalFrame);
    TGV_GF        = new TGVerticalFrame(GR_bottom_GF);
    GR_bottom_GF  ->AddFrame(TGV_GF, new TGLayoutHints(kLHintsLeft,5,5,3,4));


    // check box Gaussian filter
    //CheckBox_gaussfilter  = new TGCheckButton(TGH_lowerA, new TGHotString("Gauss filter"), -1);
    CheckBox_gaussfilter  = new TGCheckButton(TGV_GF, new TGHotString("Gauss filter"), -1);
    CheckBox_gaussfilter ->SetState(kButtonUp);
    //CheckBox_gaussfilter ->Connect("Clicked()", "voxResTree", this,Form("DoNewDataSelection(Int_t=%d)",0));

    // check box sector average
    CheckBox_sectoraverage  = new TGCheckButton(TGV_GF, new TGHotString("sector average"), -1);
    CheckBox_sectoraverage ->SetState(kButtonUp);

    // Apply Gaussian filter button
    Button_applyGF = new TGTextButton(TGV_GF, "&Apply filter");
    //Button_applyGF->Connect("Clicked()", "voxResTree", this,Form("DoNewDataSelection(Int_t=%d)",0));
    Button_applyGF->Connect("Clicked()", "voxResTree", this,Form("Loop(Int_t=%d,%d)",-1,global_position));
    Button_applyGF->SetToolTipText("Apply Gaussian filter");
    //TGH_lowerA->AddFrame(Button_export, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    TGH_DeltaX_GF          = new TGHorizontalFrame(TGV_GF);
    TGNum_DeltaX_GF       = new TGNumberEntry(TGH_DeltaX_GF, 3, 3,(TGNumberFormat::EStyle) 0);
    TGH_DeltaX_GF ->AddFrame(TGNum_DeltaX_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGL_DeltaX_GF  = new TGLabel(TGH_DeltaX_GF,"DX");
    TGH_DeltaX_GF ->AddFrame(TGL_DeltaX_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));


    TGH_DeltaY_GF          = new TGHorizontalFrame(TGV_GF);
    TGNum_DeltaY_GF       = new TGNumberEntry(TGH_DeltaY_GF, 3, 3,(TGNumberFormat::EStyle) 0);
    TGH_DeltaY_GF ->AddFrame(TGNum_DeltaY_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGL_DeltaY_GF  = new TGLabel(TGH_DeltaY_GF,"DY");
    TGH_DeltaY_GF ->AddFrame(TGL_DeltaY_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));


    TGH_DeltaZ_GF          = new TGHorizontalFrame(TGV_GF);
    TGNum_DeltaZ_GF       = new TGNumberEntry(TGH_DeltaZ_GF, 3, 3,(TGNumberFormat::EStyle) 0);
    TGH_DeltaZ_GF ->AddFrame(TGNum_DeltaZ_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGL_DeltaZ_GF  = new TGLabel(TGH_DeltaZ_GF,"DZ");
    TGH_DeltaZ_GF ->AddFrame(TGL_DeltaZ_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));


    TGH_sigma_GF          = new TGHorizontalFrame(TGV_GF);
    TGNum_sigma_GF       = new TGNumberEntry(TGH_sigma_GF, 1.2, 3,(TGNumberFormat::EStyle) 1);
    TGH_sigma_GF ->AddFrame(TGNum_sigma_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGL_sigma_GF  = new TGLabel(TGH_sigma_GF,"sigma");
    TGH_sigma_GF ->AddFrame(TGL_sigma_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));


    //TGNum_DeltaX_GF       ->Connect("ValueSet(Long_t)", "voxResTree",this, "Update_DY_X_vs_Z()");
    //TGNum_DeltaY_GF       ->Connect("ValueSet(Long_t)", "voxResTree",this, "Update_DY_X_vs_Z()");
    //TGV_GF        ->AddFrame(CheckBox_invert,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(CheckBox_gaussfilter,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(CheckBox_sectoraverage,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(TGH_DeltaX_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(TGH_DeltaY_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(TGH_DeltaZ_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(TGH_sigma_GF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_GF        ->AddFrame(Button_applyGF,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

    TGV_data_selection_master    ->AddFrame(GR_data_selection_A, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_selection_master    ->AddFrame(GR_data_selection_B, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGV_ExSel_master->AddFrame(GR_Exit, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_ExSel_master->AddFrame(GR_Select, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGV_DSel_master->AddFrame(GR_delta, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    TGV_data_output_master       ->AddFrame(GR_data_output, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_invert_master       ->AddFrame(GR_data_invert, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGV_data_invert_master       ->AddFrame(GR_ratio, new TGLayoutHints(kLHintsLeft,5,5,3,4));

    // TGH_lowerA    ->AddFrame(GR_Exit, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    // TGH_lowerA    ->AddFrame(GR_Select, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(TGV_ExSel_master, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(TGV_data_selection_master, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(TGV_DSel_master, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    //TGH_lowerA    ->AddFrame(CheckBox_scanData, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    //TGH_lowerA    ->AddFrame(CheckBox_invert, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(GR_bottom_GF, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(TGV_data_invert_master, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGH_lowerA    ->AddFrame(TGV_data_output_master, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    //--------------------------------------



    GR_bottom_sliders_sector  = new TGGroupFrame(TGH_lowerB,"sector",kHorizontalFrame);
    TGH_slider_sector         = new TGHorizontalFrame(GR_bottom_sliders_sector);
    GR_bottom_sliders_sector ->AddFrame(TGH_slider_sector, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    slider_sector = new TGHSlider(TGH_slider_sector,180,kSlider2|kScaleDownRight,1); // (kSlider1,kSlider2), (kScaleDownRight, kScaleNo, kScaleBoth)
    slider_sector ->SetRange(0,8);
    slider_sector ->SetScale(5); // number of ticks
    slider_sector ->SetPosition(0);
    slider_sector ->Connect("PositionChanged(Int_t)", "voxResTree",this, "Change_number_entry()");
    TGNum_sector = new TGNumberEntry(TGH_slider_sector, 0, 12,(TGNumberFormat::EStyle) 0);
    TGNum_sector ->Connect("ValueSet(Long_t)", "voxResTree",this, "Change_slider()");

    GR_bottom_sliders_phi  = new TGGroupFrame(TGH_lowerB,"phi",kHorizontalFrame);
    TGH_slider_phi         = new TGHorizontalFrame(GR_bottom_sliders_phi);
    GR_bottom_sliders_phi->AddFrame(TGH_slider_phi, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    slider_phi = new TGHSlider(TGH_slider_phi,180,kSlider2|kScaleDownRight,1); // (kSlider1,kSlider2), (kScaleDownRight, kScaleNo, kScaleBoth)
    slider_phi ->SetRange(0,14);
    slider_phi ->SetScale(5); // number of ticks
    slider_phi ->SetPosition(0);
    slider_phi ->Connect("PositionChanged(Int_t)", "voxResTree",this, "Change_number_entry()");
    TGNum_phi = new TGNumberEntry(TGH_slider_phi, 0, 12,(TGNumberFormat::EStyle) 0);
    TGNum_phi ->Connect("ValueSet(Long_t)", "voxResTree",this, "Change_slider()");


    GR_bottom_sliders_zbin  = new TGGroupFrame(TGH_lowerB,"z bin",kHorizontalFrame);
    TGH_slider_zbin         = new TGHorizontalFrame(GR_bottom_sliders_zbin);
    GR_bottom_sliders_zbin->AddFrame(TGH_slider_zbin, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    slider_zbin = new TGHSlider(TGH_slider_zbin,180,kSlider2|kScaleDownRight,1); // (kSlider1,kSlider2), (kScaleDownRight, kScaleNo, kScaleBoth)
    slider_zbin ->SetRange(0,9);
    slider_zbin ->SetScale(5); // number of ticks
    slider_zbin ->SetPosition(0);
    slider_zbin ->Connect("PositionChanged(Int_t)", "voxResTree",this, "Change_number_entry()");
    TGNum_zbin = new TGNumberEntry(TGH_slider_zbin, 0, 12,(TGNumberFormat::EStyle) 0);
    TGNum_zbin ->Connect("ValueSet(Long_t)", "voxResTree",this, "Change_slider()");

    GR_bottom_range  = new TGGroupFrame(TGH_lowerB,"z axis range",kVerticalFrame);
    TGV_range        = new TGVerticalFrame(GR_bottom_range);
    GR_bottom_range  ->AddFrame(TGV_range, new TGLayoutHints(kLHintsLeft,5,5,3,4));
    TGNum_zmin       = new TGNumberEntry(TGV_range, -3.0, 12,(TGNumberFormat::EStyle) 1);
    TGNum_zmax       = new TGNumberEntry(TGV_range, +3.0, 12,(TGNumberFormat::EStyle) 1);
    TGNum_zmin       ->Connect("ValueSet(Long_t)", "voxResTree",this, "Update_DY_X_vs_Z()");
    TGNum_zmax       ->Connect("ValueSet(Long_t)", "voxResTree",this, "Update_DY_X_vs_Z()");
    TGV_range        ->AddFrame(TGNum_zmin,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGV_range        ->AddFrame(TGNum_zmax,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));



    vec_TGH_general[0] ->Resize(1600,550);
    vec_TGH_general[1] ->Resize(20,20);

    //vec_TGH_lower_split[0] ->Resize(1450,550);

    TGH_slider_phi     ->AddFrame(slider_phi,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_slider_phi     ->AddFrame(TGNum_phi,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_slider_zbin    ->AddFrame(slider_zbin,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_slider_zbin    ->AddFrame(TGNum_zbin,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_slider_sector  ->AddFrame(slider_sector,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_slider_sector  ->AddFrame(TGNum_sector,new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    TGH_lowerB ->AddFrame(GR_bottom_sliders_sector, new TGLayoutHints(kLHintsLeft,2,2,2,2));
    TGH_lowerB ->AddFrame(GR_bottom_sliders_phi, new TGLayoutHints(kLHintsLeft,2,2,2,2));
    TGH_lowerB ->AddFrame(GR_bottom_sliders_zbin, new TGLayoutHints(kLHintsLeft,2,2,2,2));
    TGH_lowerB ->AddFrame(GR_bottom_range, new TGLayoutHints(kLHintsLeft,2,2,2,2));


    Frame_Setup ->Resize(1600,1050); // size of frame
    Frame_Setup ->MapSubwindows();
    Frame_Setup ->MapWindow();
    Frame_Setup ->Move(0,5); // position of frame
    //------------------------------------------------------------------------

    DoNewDataSelection(0);
    //DoNewDataSelection();
    Update_data_buttons(0,0);
    Update_data_buttons(0,1); // init


    outputfile = new TFile("dZ_over_z_fit_params.root","RECREATE");
    //output_tree = new TTree("output_tree","fit params");
    output_tree.SetName("output_tree");
    output_tree.SetTitle("fit params");

    output_tree.Branch("p0negZ_all",&p0negZ[0]);
    output_tree.Branch("p1negZ_all",&p1negZ[0]);
    output_tree.Branch("p0posZ_all",&p0posZ[0]);
    output_tree.Branch("p1posZ_all",&p1posZ[0]);
    output_tree.Branch("p0negZ_inner",&p0negZ[1]);
    output_tree.Branch("p1negZ_inner",&p1negZ[1]);
    output_tree.Branch("p0posZ_inner",&p0posZ[1]);
    output_tree.Branch("p1posZ_inner",&p1posZ[1]);
    output_tree.Branch("p0negZ_middle",&p0negZ[2]);
    output_tree.Branch("p1negZ_middle",&p1negZ[2]);
    output_tree.Branch("p0posZ_middle",&p0posZ[2]);
    output_tree.Branch("p1posZ_middle",&p1posZ[2]);
    output_tree.Branch("p0negZ_outer",&p0negZ[3]);
    output_tree.Branch("p1negZ_outer",&p1negZ[3]);
    output_tree.Branch("p0posZ_outer",&p0posZ[3]);
    output_tree.Branch("p1posZ_outer",&p1posZ[3]);
    output_tree.Branch("start_time",&start_time,"start_time/l");
    output_tree.Branch("end_time",&end_time,"end_time/l");
    output_tree.Branch("start_TF",&start_TF);
    output_tree.Branch("end_TF",&end_TF);

}

Bool_t voxResTree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void voxResTree::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t voxResTree::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

//---------------------------------------------------------------------------------
void voxResTree::Change_slider()
{
    sector_plot = TGNum_sector->GetNumberEntry()->GetNumber();
    slider_sector ->SetPosition(sector_plot);

    phi_plot = TGNum_phi->GetNumberEntry()->GetNumber();
    slider_phi ->SetPosition(phi_plot);

    zbin_plot = TGNum_zbin->GetNumberEntry()->GetNumber();
    slider_zbin ->SetPosition(zbin_plot);

    Update_DY_X_vs_Z();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
void voxResTree::Change_number_entry()
{
    sector_plot = slider_sector ->GetPosition();
    TGNum_sector ->SetNumber(sector_plot);

    phi_plot = slider_phi ->GetPosition();
    TGNum_phi ->SetNumber(phi_plot);

    zbin_plot = slider_zbin ->GetPosition();
    TGNum_zbin ->SetNumber(zbin_plot);

    Update_DY_X_vs_Z();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
void voxResTree::Update_DY_X_vs_Z()
{
    printf("voxResTree::Update_DY_X_vs_Z() \n");
    Int_t flag_loop = 0;
    //if(CheckBox_scanData ->GetState() == kButtonDown) flag_loop = 1;


    Int_t xyz_bin = plot_data_use%3;
    Int_t data_type = (Int_t)(plot_data_use/3.0);

    printf("\n------------------------------------\n");
    printf("Button pressed: %d\n",plot_data_use);
    printf("Use File: %d\n",selected_file);
    printf("Flag 1 set to: %d; Flag 2 set to: %d\n",flag_select1,flag_select2);
    printf("\n------------------------------------\n");


    Double_t range_plot_z_min = TGNum_zmin ->GetNumberEntry()->GetNumber();
    Double_t range_plot_z_max = TGNum_zmax ->GetNumberEntry()->GetNumber();
    if(range_plot_z_max <= range_plot_z_min)
    {
        printf("WARNING: z-axis range not valid!, set to -3.0..+3.0 \n");
        range_plot_z_min = -3.0;
        range_plot_z_max = 3.0;
        TGNum_zmin ->SetNumber(range_plot_z_min);
        TGNum_zmax ->SetNumber(range_plot_z_max);
    }


    if(data_type < 2 && flag_plot_data_type >= 20)
    {
        range_plot_z_min = -3.0;
        range_plot_z_max = 3.0;
        TGNum_zmin ->SetNumber(range_plot_z_min);
        TGNum_zmax ->SetNumber(range_plot_z_max);
        flag_plot_data_type = 0;
    }


    if(data_type == 2 && xyz_bin == 0 && flag_plot_data_type != 20) // statistics
    {
        range_plot_z_min = 0.0;
        range_plot_z_max = 1000.0;
        TGNum_zmin ->SetNumber(range_plot_z_min);
        TGNum_zmax ->SetNumber(range_plot_z_max);
        flag_plot_data_type = 20;
    }

    if(data_type == 2 && xyz_bin == 1 && flag_plot_data_type != 21) // flag
    {
        range_plot_z_min = 0.0;
        range_plot_z_max = 200.0;
        TGNum_zmin ->SetNumber(range_plot_z_min);
        TGNum_zmax ->SetNumber(range_plot_z_max);
        flag_plot_data_type = 21;
    }




    if(sector_plot >= 0 && phi_plot >= 0 && sector_plot < 9 && phi_plot < 15)
    {
        //printf("sector_plot: %d, phi_plot: %d \n",sector_plot,phi_plot);

        can_h2D_DY_X_vs_Z ->cd();  //Here the histogram gets drawn in the left canvas
        if(data_type == 0)
        {
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetStats(0);
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetTitle("");
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->CenterTitle();
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->CenterTitle();
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetTitle("Z (cm)");
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitle("R (cm)");
            if(xyz_bin == 0) vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaX (cm)");
            if(xyz_bin == 1) vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaY (cm)");
            if(xyz_bin == 2) vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaZ (cm)");
            vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
            // vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->DrawCopy("colz2");

            TH1D *h2D_left = (TH1D*)vec_h2D_DY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin]->Clone("h2D_left");  
            if(ratio_pressed) h2D_left->Divide(vec_h2D_DY_X_vs_Z[1][sector_plot][phi_plot][xyz_bin]);
            h2D_left->Draw("colz2");
        }
        if(data_type == 1)
        {
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetStats(0);
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetTitle("");
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->CenterTitle();
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->CenterTitle();
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetTitle("Z (cm)");
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitle("R (cm)");
            if(xyz_bin == 0) vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSX (cm)");
            if(xyz_bin == 1) vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSY (cm)");
            if(xyz_bin == 2) vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSZ (cm)");
            vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
            // vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->DrawCopy("colz2");

            TH1D *h2D_left = (TH1D*)vec_h2D_DSY_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin]->Clone("h2D_left");  
            if(ratio_pressed) h2D_left->Divide(vec_h2D_DSY_X_vs_Z[1][sector_plot][phi_plot][xyz_bin]);
            h2D_left->Draw("colz2");
        }
        if(data_type == 2)
        {
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetStats(0);
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->SetTitle("");
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->CenterTitle();
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->CenterTitle();
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetXaxis()->SetTitle("Z (cm)");
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetYaxis()->SetTitle("R (cm)");
            if(xyz_bin == 0) vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("entries");
            if(xyz_bin == 1) vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("flag");
            if(xyz_bin == 2) vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetTitle("dYSigMAD");
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
            vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin] ->DrawCopy("colz2");

            TH1D *h2D_left = (TH1D*)vec_h2D_stat_X_vs_Z[selected_file][sector_plot][phi_plot][xyz_bin]->Clone("h2D_left");  
            if(ratio_pressed) h2D_left->Divide(vec_h2D_stat_X_vs_Z[1][sector_plot][phi_plot][xyz_bin]);
            h2D_left->Draw("colz2");
        }

        Draw_R_Z_TPC(zbin_plot,0,kWhite);
        Draw_R_Z_TPC(zbin_plot,1,kBlue);
        can_h2D_DY_X_vs_Z ->Update();
    }


    can_h2D_DY_Y_vs_X ->cd();
    if(data_type == 0)
    {
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetStats(0);
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetTitle("");
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->CenterTitle();
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->CenterTitle();
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetTitle("X (cm)");
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitle("Y (cm)");
        if(xyz_bin == 0) vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaX (cm)");
        if(xyz_bin == 1) vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaY (cm)");
        if(xyz_bin == 2) vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaZ (cm)");
        vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
        // vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->DrawCopy("colz2");

        TH1D *h2D_right = (TH1D*)vec_h2D_DY_Y_vs_X[selected_file][zbin_plot][xyz_bin]->Clone("h2D_right");  
        if(ratio_pressed) h2D_right->Divide(vec_h2D_DY_Y_vs_X[1][zbin_plot][xyz_bin]);
        h2D_right->Draw("colz2");

    }
    if(data_type == 1)
    {
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetStats(0);
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetTitle("");
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->CenterTitle();
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->CenterTitle();
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetTitle("X (cm)");
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitle("Y (cm)");
        if(xyz_bin == 0) vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSX (cm)");
        if(xyz_bin == 1) vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSY (cm)");
        if(xyz_bin == 2) vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("#DeltaSZ (cm)");
        vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
        // vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->DrawCopy("colz2");

        TH1D *h2D_right = (TH1D*)vec_h2D_DSY_Y_vs_X[selected_file][zbin_plot][xyz_bin]->Clone("h2D_right");  
        if(ratio_pressed) h2D_right->Divide(vec_h2D_DSY_Y_vs_X[1][zbin_plot][xyz_bin]);
        h2D_right->Draw("colz2");
    }
    if(data_type == 2)
    {
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetStats(0);
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->SetTitle("");
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetNdivisions(505,'N');
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetNdivisions(505,'N');
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->CenterTitle();
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->CenterTitle();
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitleOffset(1.0);
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetXaxis()->SetTitle("X (cm)");
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetYaxis()->SetTitle("Y (cm)");
        if(xyz_bin == 0) vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("entries");
        if(xyz_bin == 1) vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("flag");
        if(xyz_bin == 2) vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetTitle("dYSigMAD");
        vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->GetZaxis()->SetRangeUser(range_plot_z_min,range_plot_z_max);
        // vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin] ->DrawCopy("colz2");

        TH1D *h2D_right = (TH1D*)vec_h2D_stat_Y_vs_X[selected_file][zbin_plot][xyz_bin]->Clone("h2D_right");  
        if(ratio_pressed) h2D_right->Divide(vec_h2D_stat_Y_vs_X[1][zbin_plot][xyz_bin]);
        h2D_right->Draw("colz2");
    }

    Draw_phi_sector_TPC(sector_plot,phi_plot,kWhite);
    Draw_phi_sector_TPC(sector_plot+9,phi_plot,kBlue);
    can_h2D_DY_Y_vs_X ->Update();


    can_h2D_DZ_vs_Z ->cd(1);
    can_h2D_DZ_vs_Z ->cd(1)->SetLogz(1);
    vec_h2D_DZ_vs_Z[0][0] ->SetStats(0);
    vec_h2D_DZ_vs_Z[0][0] ->SetTitle("");
    vec_h2D_DZ_vs_Z[0][0] ->GetXaxis()->SetNdivisions(505,'N');
    vec_h2D_DZ_vs_Z[0][0] ->GetYaxis()->SetNdivisions(505,'N');
    vec_h2D_DZ_vs_Z[0][0] ->GetXaxis()->CenterTitle();
    vec_h2D_DZ_vs_Z[0][0] ->GetYaxis()->CenterTitle();
    vec_h2D_DZ_vs_Z[0][0] ->GetXaxis()->SetTitleOffset(1.1);
    vec_h2D_DZ_vs_Z[0][0] ->GetYaxis()->SetTitleOffset(1.1);
    vec_h2D_DZ_vs_Z[0][0] ->GetXaxis()->SetTitle("Z (cm)");
    vec_h2D_DZ_vs_Z[0][0] ->GetYaxis()->SetTitle("#DeltaZ (cm)");
    vec_h2D_DZ_vs_Z[0][0] ->GetZaxis()->SetTitle("entries");
    //vec_h2D_DZ_vs_Z[0][0] ->GetZaxis()->SetRangeUser(-5,5);
    vec_h2D_DZ_vs_Z[0][0] ->DrawCopy("colz");
    vec_TP_DZ_vs_Z[0][0]  ->SetLineColor(kGray+2);
    vec_TP_DZ_vs_Z[0][0]  ->SetLineWidth(3);
    vec_TP_DZ_vs_Z[0][0]  ->SetLineStyle(1);
    vec_TP_DZ_vs_Z[0][0]  ->DrawCopy("same hist");

    //--------------------------------------
    Double_t z_range_min[2] = {-200.0,20.0};
    Double_t z_range_max[2] = {-20.0,200.0};
    Int_t z_color_fit[2] = {kRed,kBlue};
    Float_t fit_params[4][2][2] = {0.0}; // radius, z, par0/par1
    for(Int_t i_radius = 0; i_radius < 4; i_radius++)
    {
        for(Int_t i_z = 0; i_z < 2; i_z++)
        {
            for(Int_t i = 0; i < 6; i++)
            {
                func_PolyFitFunc ->SetParameter(i,0.0);
                func_PolyFitFunc ->SetParError(i,0.0);
            }
            func_PolyFitFunc ->SetParameter(0,0.0);
            func_PolyFitFunc ->SetParameter(1,0.0);
            func_PolyFitFunc ->FixParameter(2,0.0);
            func_PolyFitFunc ->FixParameter(3,0.0);
            func_PolyFitFunc ->FixParameter(4,0.0);
            func_PolyFitFunc ->FixParameter(5,0.0);

            vec_TP_DZ_vs_Z[0][i_radius] ->Fit("func_PolyFitFunc","QWMN","",z_range_min[i_z],z_range_max[i_z]);

            func_PolyFitFunc ->SetLineColor(z_color_fit[i_z]);
            func_PolyFitFunc ->SetLineStyle(1);
            func_PolyFitFunc ->SetLineWidth(4);
            func_PolyFitFunc ->SetRange(z_range_min[i_z],z_range_max[i_z]);
            if(i_radius == 0)
            {
                func_PolyFitFunc ->DrawCopy("same");

                HistName = "p0 = ";
                sprintf(NoP,"%4.3f",(Double_t)func_PolyFitFunc->GetParameter(0));
                HistName += NoP;
                HistName += ", p1 = ";
                sprintf(NoP,"%4.3f",(Double_t)func_PolyFitFunc->GetParameter(1));
                HistName += NoP;
                if(i_z == 0) plotTopLegend((char*)HistName.Data(),0.18,0.3,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                if(i_z == 1) plotTopLegend((char*)HistName.Data(),0.54,0.7,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            }
        }
        //--------------------------------------


        // Truncation
        TH1D* h_dummy_proj;
        for(Int_t i_bin_x = 1; i_bin_x <= vec_h2D_DZ_vs_Z[0][i_radius]->GetNbinsX(); i_bin_x++)
        {
            Double_t x_val = vec_h2D_DZ_vs_Z[0][i_radius]->GetXaxis()->GetBinCenter(i_bin_x);
            h_dummy_proj = (TH1D*)vec_h2D_DZ_vs_Z[0][i_radius] ->ProjectionY("dummy",i_bin_x,i_bin_x);
            Double_t rms_proj  = h_dummy_proj ->GetRMS();
            Double_t mean_proj = h_dummy_proj ->GetMean();
            for(Int_t i_bin_y = 1; i_bin_y <= vec_h2D_DZ_vs_Z[0][i_radius]->GetNbinsY(); i_bin_y++)
            {
                Double_t y_val = vec_h2D_DZ_vs_Z[0][i_radius]->GetYaxis()->GetBinCenter(i_bin_y);
                if(fabs(y_val - mean_proj) < 2.0*rms_proj)
                {
                    Double_t bin_cont = vec_h2D_DZ_vs_Z[0][i_radius] ->GetBinContent(i_bin_x,i_bin_y);
                    vec_h2D_DZ_vs_Z[1][i_radius] ->SetBinContent(i_bin_x,i_bin_y,bin_cont);
                }
            }
            delete h_dummy_proj;
        }

        HistName = "vec_TP_DZ_vs_Z_trunc_";
        HistName += i_radius;
        vec_TP_DZ_vs_Z[1][i_radius] = vec_h2D_DZ_vs_Z[1][i_radius] ->ProfileX(HistName.Data());

        if(i_radius == 0)
        {
            can_h2D_DZ_vs_Z ->cd(2)->SetLogz(1);
            vec_h2D_DZ_vs_Z[1][i_radius] ->SetStats(0);
            vec_h2D_DZ_vs_Z[1][i_radius] ->SetTitle("");
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetXaxis()->SetNdivisions(505,'N');
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetYaxis()->SetNdivisions(505,'N');
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetXaxis()->CenterTitle();
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetYaxis()->CenterTitle();
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetXaxis()->SetTitleOffset(1.1);
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetYaxis()->SetTitleOffset(1.1);
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetXaxis()->SetTitle("Z (cm)");
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetYaxis()->SetTitle("#DeltaZ (cm)");
            vec_h2D_DZ_vs_Z[1][i_radius] ->GetZaxis()->SetTitle("entries");
            vec_h2D_DZ_vs_Z[1][i_radius] ->DrawCopy("colz");
            vec_TP_DZ_vs_Z[1][i_radius]  ->SetLineColor(kGray+2);
            vec_TP_DZ_vs_Z[1][i_radius]  ->SetLineWidth(3);
            vec_TP_DZ_vs_Z[1][i_radius]  ->SetLineStyle(1);
            vec_TP_DZ_vs_Z[1][i_radius]  ->DrawCopy("same hist");
        }


        //--------------------------------------
        for(Int_t i_z = 0; i_z < 2; i_z++)
        {
            for(Int_t i = 0; i < 6; i++)
            {
                func_PolyFitFunc ->SetParameter(i,0.0);
                func_PolyFitFunc ->SetParError(i,0.0);
            }
            func_PolyFitFunc ->SetParameter(0,0.0);
            func_PolyFitFunc ->SetParameter(1,0.0);
            func_PolyFitFunc ->FixParameter(2,0.0);
            func_PolyFitFunc ->FixParameter(3,0.0);
            func_PolyFitFunc ->FixParameter(4,0.0);
            func_PolyFitFunc ->FixParameter(5,0.0);

            vec_TP_DZ_vs_Z[1][i_radius] ->Fit("func_PolyFitFunc","QWMN","",z_range_min[i_z],z_range_max[i_z]);

            func_PolyFitFunc ->SetLineColor(z_color_fit[i_z]);
            func_PolyFitFunc ->SetLineStyle(1);
            func_PolyFitFunc ->SetLineWidth(4);
            func_PolyFitFunc ->SetRange(z_range_min[i_z],z_range_max[i_z]);

            if(i_radius == 0)
            {
                func_PolyFitFunc ->DrawCopy("same");

                HistName = "p0 = ";
                sprintf(NoP,"%4.3f",(Double_t)func_PolyFitFunc->GetParameter(0));
                HistName += NoP;
                HistName += ", p1 = ";
                sprintf(NoP,"%4.3f",(Double_t)func_PolyFitFunc->GetParameter(1));
                HistName += NoP;
                if(i_z == 0) plotTopLegend((char*)HistName.Data(),0.18,0.3,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                if(i_z == 1) plotTopLegend((char*)HistName.Data(),0.54,0.7,0.035,z_color_fit[i_z],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
            }

            fit_params[i_radius][i_z][0] = (Double_t)func_PolyFitFunc->GetParameter(0);
            fit_params[i_radius][i_z][1] = (Double_t)func_PolyFitFunc->GetParameter(1);
        }

        //--------------------------------------
    } // end of radius loop

    Double_t average_p0 = (-fit_params[0][0][0] + fit_params[0][1][0])/2.0;
    Double_t average_p1 = (fit_params[0][0][1] + fit_params[0][1][1])/2.0;
    HistName = "<p0> = ";
    sprintf(NoP,"%4.3f",(Double_t)average_p0);
    HistName += NoP;
    HistName += ", <p1> = ";
    sprintf(NoP,"%4.3f",(Double_t)average_p1);
    HistName += NoP;
    plotTopLegend((char*)HistName.Data(),0.32,0.85,0.035,kGreen+2,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    //--------------------------------------
    TString TS_inputpos = "";
    TS_inputpos =  fCombo_file[selected_file]->GetSelectedEntry()->GetTitle();

    for(Int_t i_radius = 0; i_radius < 4; i_radius++)
    {
        printf("Fit values for root file: %s, radius: %d, z < 0 {p0: %4.6f, p1: %4.6f}, z > 0 {p0: %4.6f, p1: %4.6f} \n",TS_inputpos.Data(),i_radius,fit_params[i_radius][0][0],fit_params[i_radius][0][1],fit_params[i_radius][1][0],fit_params[i_radius][1][1]);
    }


    can_h2D_DZ_vs_Z ->Update();


    if(flag_loop)
    {
        //---------------------
        TString st_name = TS_inputpos;
        Int_t first = st_name.First("_");
        st_name.Remove(0,first+1);

        first = st_name.First("_");
        st_name.Remove(0,first+1);

        TString st_name_A = st_name;
        first = st_name_A.First("_");
        st_name_A.Remove(first,st_name_A.Length());
        st_name.Remove(0,first+1);
        ULong64_t long_start_time = st_name_A.Atoll();

        st_name_A = st_name;
        first = st_name_A.First("_");
        st_name_A.Remove(first,st_name_A.Length());
        st_name.Remove(0,first+1);
        ULong64_t long_end_time = st_name_A.Atoll();

        st_name_A = st_name;
        first = st_name_A.First("_");
        st_name_A.Remove(first,st_name_A.Length());
        st_name.Remove(0,first+1);
        Long64_t long_start_TF = st_name_A.Atoll();

        st_name_A = st_name;
        first = st_name_A.First(".");
        st_name_A.Remove(first,st_name_A.Length());
        st_name.Remove(0,first+1);
        Long64_t long_end_TF = st_name_A.Atoll();

        printf("time values from file name: %lld, %lld, %lld, %lld \n",long_start_time,long_end_time,long_start_TF,long_end_TF);
        start_time = (ULong64_t)long_start_time;
        end_time   = (ULong64_t)long_end_time;
        start_TF   = (Float_t)long_start_TF;
        end_TF     = (Float_t)long_end_TF;
        //---------------------

        for(Int_t i_radius = 0; i_radius < 4; i_radius++)
        {
            p0negZ[i_radius] = fit_params[i_radius][0][0];
            p1negZ[i_radius] = fit_params[i_radius][0][1];
            p0posZ[i_radius] = fit_params[i_radius][1][0];
            p1posZ[i_radius] = fit_params[i_radius][1][1];

            //printf("i_radius: %d, p0negZ: %4.5f \n",i_radius,p0negZ[i_radius]);
        }

        output_tree.Fill();
 

        Int_t N_entries_list = fCombo_file[selected_file]->GetNumberOfEntries();
        Int_t selected = fCombo_file[selected_file]->GetSelected();
        if(selected < (N_entries_list -1))
        {
            fCombo_file[selected_file]->Select(selected+1);
        }
        else
        {
            outputfile ->cd();
            output_tree.Write();
        }
    }
    //--------------------------------------
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
void voxResTree::Update_data_buttons(Int_t i_button_A, Int_t i_button_B)
{
    //printf("A,B: {%d, %d} \n",i_button_A,i_button_B);
    //enum  	EButtonState { kButtonUp , kButtonDown , kButtonEngaged , kButtonDisabled }
    if(i_button_B == 0)
    {
        for(Int_t i_button = 0; i_button < (Int_t)vec_TRB_plot_data.size(); i_button++)
        {
            vec_TRB_plot_data[i_button] ->SetState(kButtonUp);
        }
        vec_TRB_plot_data[i_button_A] ->SetState(kButtonDown);

        plot_data_use = i_button_A;
    }
    else if(i_button_B == 1 && flag_select1 && flag_select2)
    {
        for(Int_t i_button = 0; i_button < (Int_t)select_button.size(); i_button++)
        {
            select_button[i_button] ->SetState(kButtonUp);
        }
        select_button[i_button_A] ->SetState(kButtonDown);

        selected_file = i_button_A;
    }
    else if(i_button_B == 1 && flag_select1 && !flag_select2)
    {
        for(Int_t i_button = 0; i_button < (Int_t)select_button.size(); i_button++)
        {
            select_button[i_button] ->SetState(kButtonUp);
        }
        select_button[i_button_A] ->SetState(kButtonDown);
    }
    else{
        selected_file = 0;
    }
    Update_DY_X_vs_Z();


}
//---------------------------------------------------------------------------------

void voxResTree::Update_delta_buttons(Int_t button_A, Int_t button_xyz)
{
    if(button_A == 0 && button_xyz == 0) delta_pressed[0] = 0;
    else if(button_A == 0 && button_xyz == 1) delta_pressed[1] = 0;
    else if(button_A == 0 && button_xyz == 2) delta_pressed[2] = 0;
    else if(button_A == 1 && button_xyz == 0) delta_pressed[0] = 1;
    else if(button_A == 1 && button_xyz == 1) delta_pressed[1] = 1;
    else if(button_A == 1 && button_xyz == 2) delta_pressed[2] = 1;

    // if(button_A != 2){
    //     delta_button[button_A][button_xyz] ->SetState(kButtonDown);
    //     delta_button[abs(button_A-1)][button_xyz] ->SetState(kButtonUp);
    // }
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++){
        delta_button[delta_pressed[i_xyz]][i_xyz] ->SetState(kButtonDown);
        delta_button[abs(delta_pressed[i_xyz]-1)][i_xyz] ->SetState(kButtonUp);
    }
    

}

//---------------------------------------------------------------------------------
void voxResTree::Get_directory_list()
{
    printf("voxResTree::Get_directory_list() \n");
    // Get a list of directories in subdirectory "Positions"
    TString dir = gSystem->pwd();
    dir += "/Data";
    void *dirp = gSystem->OpenDirectory(dir);
    Char_t *afile;
    TString TS_afile;
    while((afile = const_cast<Char_t *>(gSystem->GetDirEntry(dirp))))
    {
        TS_afile = afile;
        Long_t id;
        Long_t size;
        Long_t flags;
        Long_t modtime;

        const TString path = dir+"/"+afile;
        if(gSystem->GetPathInfo(path.Data(),&id,&size,&flags,&modtime))
        {
            continue;
        }
        Int_t isDir = (02 & flags);
        if(!isDir)
        {
            continue;
        }
        if(TS_afile == "." || TS_afile == "..") continue;
        cout << path.Data() << endl;
        vec_TS_afile.push_back(TS_afile);
    }
}
//---------------------------------------------------------------------------------


void voxResTree::func_dummy(Int_t i_select = -1)
{
    printf("voxResTree::func_dummy(%d) \n",i_select);
    if(i_select == 5) flag_select1 = true;
    else if(i_select == 7) flag_select2 = true;
}
//---------------------------------------------------------------------------------
void voxResTree::DoExport()
{
    printf("Start map export \n");

    //------------------------------------
    // output average map
    TString TS_export =  fCombo_positions[2]->GetSelectedEntry()->GetTitle();
    TS_export += TGText_outputname ->GetDisplayText();
    cout << "Export map name: " << TS_export.Data() << endl;
    outputfileMap = new TFile(TS_export.Data(),"RECREATE");
    mTreeOut = std::make_unique<TTree>("voxResTree", "Voxel results and statistics");
    mTreeOut->Branch("voxRes", &mVoxelResultsOutPtr);

    for(Int_t ientry = 0; ientry < (Int_t)vec_VoxRes.size(); ientry++)
    {
        mVoxelResultsOut = vec_VoxRes[ientry];
        mTreeOut ->Fill();
    }

    outputfileMap ->cd();
    mTreeOut ->Write();
    mTreeOut.reset();
    outputfileMap ->Close();
    //------------------------------------
    printf("Map exported \n");
}
//---------------------------------------------------------------------------------


//---------------------------------------------------------------------------------
void voxResTree::DoNewDataSelection(Int_t i_position)
//void voxResTree::DoNewDataSelection()
{
    //Int_t i_select = 0;
    // To do: check for .root at the end of path

    //printf("voxResTree::DoNewDataSelection() for i_select: %d, i_selectB: %d \n",i_select,i_selectB);
    printf("voxResTree::DoNewDataSelection(%d) \n",i_position);
    
    TString TS_inputpos = "";
    TS_inputpos =  fCombo_positions[i_position]->GetSelectedEntry()->GetTitle();
    printf("   directory: %s \n",TS_inputpos.Data());

    TString dir = TS_inputpos;
    void *dirp = gSystem->OpenDirectory(dir);
    Char_t *afile;
    TString TS_afile;
    vec_TS_voxfiles.clear();
    vec_TS_voxfiles.resize(0);
    while((afile = const_cast<Char_t *>(gSystem->GetDirEntry(dirp))))
    {
        TS_afile = afile;
        Long_t id;
        Long_t size;
        Long_t flags;
        Long_t modtime;

        const TString path = dir+"/"+afile;
        //cout << path.Data() << endl;
        if(gSystem->GetPathInfo(path.Data(),&id,&size,&flags,&modtime))
        {
            continue;
        }
        Int_t isDir = (02 & flags);
        if(isDir)
        {
            continue;
        }
        if(TS_afile == "." || TS_afile == "..") continue;
        cout << path.Data() << endl;
        vec_TS_voxfiles.push_back(path);
    }

    fCombo_file[i_position] ->RemoveAll();
    fCombo_file[i_position] ->Clear();
    fCombo_file[i_position] ->Changed();
    for(Int_t i_pos_entry = 0; i_pos_entry < (Int_t)vec_TS_voxfiles.size(); i_pos_entry++)
    {
        fCombo_file[i_position]->AddEntry(vec_TS_voxfiles[i_pos_entry].Data(),i_pos_entry);
    }

    // Load fist file from the position:
    // fCombo_file[i_position]->Select(0);
    if(just_started){
        just_started = false;
        fCombo_file[i_position]->Select(0);
    }
    global_position = i_position;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
//void voxResTree::DoNewFileSelection(Int_t i_select)
void voxResTree::DoNewFileSelection(Int_t i_file)
{
    fCombo_file[global_position]->Select(i_file);

    
    printf("voxResTree::DoNewFileSelection(%d) \n",global_position);
    //printf("voxResTree::DoNewFileSelection(%d) \n",i_select);
    TString TS_inputpos = "";
    TS_inputpos =  fCombo_file[global_position]->GetSelectedEntry()->GetTitle();
    printf("Load root file: %s \n",TS_inputpos.Data());
    //inputfile = (TFile*)gROOT->GetListOfFiles()->FindObject(vec_TS_voxfiles[0]);
    if(inputfile) delete inputfile;
    printf("inputfile deleted \n");
    //if (!inputfile || !inputfile->IsOpen())
    {
        printf("new inputfile \n");
        inputfile = new TFile(TS_inputpos);
    }
    printf("Root file loaded \n");
    inputfile->GetObject("voxResTree",fChain);

    
    // if(selected_file == 1){TS_inputpos.Data()}

    // selected_file = i_file;

    InitChain();
    Loop(-1,global_position);
}
//---------------------------------------------------------------------------------


void voxResTree::GetRatio()
{
    ratio_pressed = !ratio_pressed;
    printf("\n\nRATIOBUTTON PRESSED\n\n");
    printf("Ratiobutton: %d",ratio_pressed);

    if(!flag_select1){
        printf("\nPlease select a FILE 1\n");
    }
    if(!flag_select2){
        printf("\nPlease select a FILE 2\n");
    }
    if(ratio_pressed && flag_select1 && flag_select2){
        printf("Getting Ratios of both files ready"); 
        Update_DY_X_vs_Z();
    }
}




#endif // #ifdef voxResTree_cxx
