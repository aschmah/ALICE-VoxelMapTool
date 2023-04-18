
using namespace std;

#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGaxis.h"
//#include <iostream.h>
#include <fstream>
#include <sstream>
#include "TMath.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TVirtualPad.h"
#include "TPolyLine.h"
#include "TVector3.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLCamera.h"
#include "TGLPerspectiveCamera.h"
#include "TGFrame.h"
#include "TGLUtil.h"
#include "TGLLightSet.h"
#include "TGLCameraOverlay.h"
#include "TLorentzVector.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"
#include "TList.h"
#include "TChain.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/GSLSimAnMinimizer.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"
//#include <GSLMultiMinimizer.h>
//#include "Math/GSLMultiMinimizer.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "TRotation.h"
#include "TSVDUnfold.h"
#include "TSystemDirectory.h"
#include "TMatrixT.h"
#include "TVector2.h"
#include "TEllipse.h"
#include "TPaveLabel.h"
#include "TBox.h"
#include "TSystemDirectory.h"
#include "TAttImage.h"
#include "TImage.h"
#include "TPaletteAxis.h"
#include "THelix.h"
#include "TView.h"
#include "TSystem.h"
#include "TASImage.h"
#include "TImage.h"
#include "TArrayD.h"
//#include "TPython.h"
#include "TArrow.h"
#include "TKey.h"
#include "TSpectrum.h"
#include "TNtuple.h"
#include "TDatime.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


//------------------------------------------------------------------------------------------------------------
static const Float_t Pi = TMath::Pi();
static TRandom ran;
static TString HistName, HistNameB;
static char NoP[50];
static const Double_t mu_zero = 4.0*TMath::Pi()*1.0E-07; // N/A^{2}, magnetic constant

static const Double_t kAlmost1=1. - Double_t(FLT_EPSILON);
static const Double_t kAlmost0=Double_t(FLT_MIN);
static const Double_t kB2C=-0.299792458e-3;
static Float_t fHelix[9];
static Float_t B_field = -0.001; // -5.0 kG
static Float_t track_pos[3];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
static const Int_t N_v2_vs_pt_BW = 14;
static Double_t mt_m0_cut;
static Double_t pt_cut_BW_global;
static Int_t flag_v2_BW_use[N_v2_vs_pt_BW];
static TGraphAsymmErrors* tgae_v2_stat_BW[N_v2_vs_pt_BW];
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
static const Float_t RowX[153] = {
   85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725, 96.475,
   97.225, 97.975, 98.725, 99.475, 100.225, 100.975, 101.725, 102.475, 103.225, 103.975, 104.725, 105.475, 106.225, 106.975, 107.725,
   108.475, 109.225, 109.975, 110.725, 111.475, 112.225, 112.975, 113.725, 114.475, 115.225, 115.975, 116.725, 117.475, 118.225, 118.975,
   119.725, 120.475, 121.225, 121.975, 122.725, 123.475, 124.225, 124.975, 125.725, 126.475, 127.225, 127.975, 128.725, 129.475, 130.225,
   130.975, 131.725, 135.200, 136.200, 137.200, 138.200, 139.200, 140.200, 141.200, 142.200, 143.200, 144.200, 145.200, 146.200, 147.200,
   148.200, 149.200, 150.200, 151.200, 152.200, 153.200, 154.200, 155.200, 156.200, 157.200, 158.200, 159.200, 160.200, 161.200, 162.200,
   163.200, 164.200, 165.200, 166.200, 167.200, 168.200, 171.400, 172.600, 173.800, 175.000, 176.200, 177.400, 178.600, 179.800, 181.000,
   182.200, 183.400, 184.600, 185.800, 187.000, 188.200, 189.400, 190.600, 191.800, 193.000, 194.200, 195.400, 196.600, 197.800, 199.000,
   200.200, 201.400, 202.600, 203.800, 205.000, 206.200, 209.650, 211.150, 212.650, 214.150, 215.650, 217.150, 218.650, 220.150, 221.650,
   223.150, 224.650, 226.150, 227.650, 229.150, 230.650, 232.150, 233.650, 235.150, 236.650, 238.150, 239.650, 241.150, 242.650, 244.150,
   245.650,246.650}; // last value added

static const Float_t z_over_x_voxel_binning[6]  = {0.0,0.2,0.4,0.6,0.8,1.0};
//------------------------------------------------------------------------------------------------------------


void set_helix(Float_t x, Float_t alpha, Float_t param[5], Float_t Bz)
{
    // http://alidoc.cern.ch/AliRoot/master/_ali_helix_8cxx_source.html
    // AliHelix::AliHelix(const AliExternalTrackParam &t)

    Float_t cs,sn;
    for(Int_t i = 0; i < 5; i++)
    {
        fHelix[i] = param[i];
    }

    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...
    fHelix[4]=fHelix[4]/(-1000/0.299792458/Bz);    // C
    //  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
    cs=cosf(alpha); sn=sinf(alpha); // RS use float versions: factor 2 in CPU speed

    Float_t xc, yc, rc;
    rc  =  1/fHelix[4];
    xc  =  x-fHelix[2]*rc;
    Float_t dummy = 1-(x-xc)*(x-xc)*fHelix[4]*fHelix[4];
    if (dummy<0)
    {
        dummy = 0;
    }
    yc  =  fHelix[0]+TMath::Sqrt(dummy)/fHelix[4];

    fHelix[6] = xc*cs - yc*sn;
    fHelix[7] = xc*sn + yc*cs;
    fHelix[8] =  TMath::Abs(rc);
    //
    //
    fHelix[5]=x*cs - fHelix[0]*sn;            // x0
    fHelix[0]=x*sn + fHelix[0]*cs;            // y0
    //fHelix[1]=                               // z0
    //  fHelix[2]=TMath::ASin(fHelix[2]) + alpha; // phi0
    float sphi = TMath::Abs(fHelix[2])<kAlmost1 ? asinf(fHelix[2]) : TMath::Sign(TMath::Pi()/2, fHelix[2]);
    fHelix[2]=sphi + alpha; // phi0 // RS : use float version
    //fHelix[3]=                               // tgl
    //
    //
    fHelix[5]   = fHelix[6];
    fHelix[0]   = fHelix[7];
}


void evaluate_helix(Float_t t,Float_t r[3])
{
    float phase=fHelix[4]*t+fHelix[2];
    Float_t sn=sinf(phase), cs=cosf(phase);
    //  Float_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

    r[0] = fHelix[5] + sn/fHelix[4];
    r[1] = fHelix[0] - cs/fHelix[4];
    r[2] = fHelix[1] + fHelix[3]*t;
}


Float_t get_helix_cluster_residuum(Float_t clus_x, Float_t clus_y, Float_t clus_z, Float_t track_path_start, Float_t delta_pos[3], Float_t pos_helix[3], Float_t delta_track_path)
{
    Float_t track_pos[3] = {0.0,0.0,0.0};
    Float_t track_pos_previous[3] = {0.0,0.0,0.0};
    Float_t delta_clus_previous = 99999.0;
    Float_t track_path_minimum = 0.0;
    for(Float_t track_path = track_path_start; fabs(track_path) < 350.0; track_path += delta_track_path)
    {
        evaluate_helix(track_path,track_pos);
        Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);

        //printf("track_path: %4.3f, radius_track: %4.3f \n",track_path,radius_track);

        if(radius_track > 250.0) break;

        Double_t delta_clus = TMath::Sqrt(TMath::Power(clus_x - track_pos[0],2) + TMath::Power(clus_y - track_pos[1],2) +TMath::Power(clus_z - track_pos[2],2));
        if(delta_clus > delta_clus_previous) // minimum reached
        {
            delta_pos[0] = track_pos_previous[0] - clus_x;
            delta_pos[1] = track_pos_previous[1] - clus_y;
            delta_pos[2] = track_pos_previous[2] - clus_z;

            pos_helix[0] = track_pos_previous[0];
            pos_helix[1] = track_pos_previous[1];
            pos_helix[2] = track_pos_previous[2];

            //printf("  --> track_path: %4.3f, delta_clus: %4.3f, delta: {%4.3f, %4.3f, %4.3f} \n",track_path,delta_clus,delta_pos[0],delta_pos[1],delta_pos[2]);

            return track_path_minimum;
        }
        else
        {
            delta_clus_previous = delta_clus;
            track_pos_previous[0] = track_pos[0];
            track_pos_previous[1] = track_pos[1];
            track_pos_previous[2] = track_pos[2];
            track_path_minimum = track_path;
        }
    }

    return track_path_minimum;
}


Float_t get_helix_cluster_residuum_2D(Float_t clus_x, Float_t clus_y, Float_t track_path_start, Float_t delta_pos[3], Float_t pos_helix[3], Float_t delta_track_path, Double_t& delta_clus)
{
    Float_t track_pos[3] = {0.0,0.0,0.0};
    Float_t track_pos_previous[3] = {0.0,0.0,0.0};
    Float_t delta_clus_previous = 99999.0;
    Float_t track_path_minimum = 0.0;
    for(Float_t track_path = track_path_start; fabs(track_path) < 350.0; track_path += delta_track_path)
    {
        evaluate_helix(track_path,track_pos);
        Double_t radius_track = TMath::Sqrt(track_pos[0]*track_pos[0] + track_pos[1]*track_pos[1]);

        //printf("track_path: %4.3f, radius_track: %4.3f \n",track_path,radius_track);

        if(radius_track > 250.0) break;

        delta_clus = TMath::Sqrt(TMath::Power(clus_x - track_pos[0],2) + TMath::Power(clus_y - track_pos[1],2));
        if(delta_clus > delta_clus_previous) // minimum reached
        {
            delta_pos[0] = track_pos_previous[0] - clus_x;
            delta_pos[1] = track_pos_previous[1] - clus_y;

            pos_helix[0] = track_pos_previous[0];
            pos_helix[1] = track_pos_previous[1];
            pos_helix[2] = track_pos_previous[2];

            //printf("  --> track_path: %4.3f, delta_clus: %4.3f, delta: {%4.3f, %4.3f, %4.3f} \n",track_path,delta_clus,delta_pos[0],delta_pos[1],delta_pos[2]);

            return track_path_minimum;
        }
        else
        {
            delta_clus_previous = delta_clus;
            track_pos_previous[0] = track_pos[0];
            track_pos_previous[1] = track_pos[1];
            track_pos_previous[2] = track_pos[2];
            track_path_minimum = track_path;
        }
    }

    return track_path_minimum;
}


TLorentzVector get_TLV_helix(Float_t Bz, Double_t& pt)
{
    TLorentzVector TLV_helix;
    Double_t charge = 1.0;
    Double_t conversion = -1000/0.299792458/Bz;
    pt = charge/(conversion*fHelix[4]);

    //printf("pt: %4.3f \n",pt);
    //fHelix[4] = charge/(conversion*pt); // C

    Float_t track_posA[3];
    evaluate_helix(-85.0,track_posA);
    Float_t track_posB[3];
    evaluate_helix(-85.0+0.1,track_posB);
    TVector3 TV3_dir;
    TV3_dir.SetXYZ(track_posB[0]-track_posA[0],track_posB[1]-track_posA[1],track_posB[2]-track_posA[2]);
    //printf(" \n");
    //printf("posA: {%4.3f, %4.3f}, posB: {%4.3f, %4.3f}, dirA: {%4.3f, %4.3f} \n",track_posA[0],track_posA[1],track_posB[0],track_posB[1],TV3_dir.X(),TV3_dir.Y());
    if(TV3_dir.Mag() > 0.0) TV3_dir *= 1.0/TV3_dir.Mag();
    //printf("dirB: {%4.3f, %4.3f} \n",TV3_dir.X(),TV3_dir.Y());
    Double_t perp = TV3_dir.Perp();
    //printf("perp: %4.3f, pt: %4.3f \n",perp,pt);
    if(perp > 0.0) TV3_dir *= fabs(pt)/perp;
    //printf("dirC: {%4.3f, %4.3f} \n",TV3_dir.X(),TV3_dir.Y());
    //printf("pt: %4.3f, perp: %4.3f \n",pt,TV3_dir.Perp());

    TLV_helix.SetXYZM(TV3_dir.X(),TV3_dir.Y(),TV3_dir.Z(),0.139);
    return TLV_helix;
}



//----------------------------------------------------------------------------------------
void Calc_Intersect_Straight_Y(TVector3 vec_spaceA, TVector3 vec_spaceB, Double_t y_pos, TVector3 &vec_intersect)
{
    // Calculates the intersect vec_intersect of a straight line, defined by two space points vec_spaceA and vec_spaceB in 3D,
    // with an plane perpendicular to the y-axis at pos y_pos

    vec_intersect.SetXYZ(-999,-999,-999);

    if( (vec_spaceB.Y()-vec_spaceA.Y()) != 0 )
    {
        Double_t Slope_x = (vec_spaceB.X()-vec_spaceA.X())/(vec_spaceB.Y()-vec_spaceA.Y());
        Double_t axis_x  = vec_spaceB.X() - Slope_x*vec_spaceB.Y();
        Double_t x_hit   = Slope_x*y_pos + axis_x;

        Double_t Slope_z = (vec_spaceB.Z()-vec_spaceA.Z())/(vec_spaceB.Y()-vec_spaceA.Y());
        Double_t axis_z  = vec_spaceB.Z() - Slope_z*vec_spaceB.Y();
        Double_t z_hit   = Slope_z*y_pos + axis_z;

        Double_t y_hit   = y_pos;

        vec_intersect.SetXYZ(x_hit,y_hit,z_hit);
    }

}
//----------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
TVector3 Calc_Circle_Center(TVector3 v_pointA, TVector3 v_pointB, TVector3 v_pointC)
{
    TVector3 v_center_of_circle;

    // Calculate center of circle
    Double_t xy_circle[2] = {0.0,0.0};
    Double_t x12  = v_pointA.X() - v_pointB.X();
    Double_t x13  = v_pointA.X() - v_pointC.X();
    Double_t y12  = v_pointA.Y() - v_pointB.Y();
    Double_t y13  = v_pointA.Y() - v_pointC.Y();
    Double_t x12p = v_pointA.X() + v_pointB.X();
    Double_t x13p = v_pointA.X() + v_pointC.X();
    Double_t y12p = v_pointA.Y() + v_pointB.Y();
    Double_t y13p = v_pointA.Y() + v_pointC.Y();

    xy_circle[0] = (y12*(x13*x13p + y13*y13p) - y13*(x12*x12p + y12*y12p))/(2.0*(x13*y12 - x12*y13));
    xy_circle[1] = (x13*(x12*x12p + y12*y12p) - x12*(x13*x13p + y13*y13p))/(2.0*(x13*y12 - x12*y13));

    v_center_of_circle.SetXYZ(xy_circle[0],xy_circle[1],0.0);
    return v_center_of_circle;
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void Calc_average_circle_center(std::vector<TVector3> TV3_points, Double_t &average_radius, TVector3 &average_circle_center)
{
    Int_t n_combinations = 0;
    average_radius = 0.0;
    average_circle_center.SetXYZ(0.0,0.0,0.0);
    for(Int_t i_pointA = 0; i_pointA < TV3_points.size(); i_pointA++)
    {
        for(Int_t i_pointB = (i_pointA+1); i_pointB < TV3_points.size(); i_pointB++)
        {
            for(Int_t i_pointC = (i_pointB+1); i_pointC < TV3_points.size(); i_pointC++)
            {
                TVector3 circle_center = Calc_Circle_Center(TV3_points[i_pointA],TV3_points[i_pointB],TV3_points[i_pointC]);
                average_circle_center += circle_center;
                TVector3 vec_radius = TV3_points[i_pointA] - circle_center;
                average_radius += vec_radius.Mag();
                n_combinations++;
            }
        }
    }

    if(n_combinations > 0)
    {
        average_radius        /= ((Double_t)n_combinations);
        average_circle_center *= 1.0/((Double_t)n_combinations);
    }
}
//------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------
void Convert_NDC_to_User(TPad* iPad, Double_t x_NDC, Double_t y_NDC,
                         Double_t &x_User, Double_t &y_User,
                         Int_t flag_log_x, Int_t flag_log_y)
{
    // Convert NDC coordinates to User coordinates, if x/y axes are in log scale set flag_log to 1
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_User = x_NDC*xr + iPad->GetX1();
    y_User = y_NDC*yr + iPad->GetY1();

    if(flag_log_x == 1) x_User = TMath::Power(10,x_User);
    if(flag_log_y == 1) y_User = TMath::Power(10,y_User);
    //cout << "X1: " << iPad->GetX1() << ", X2: " << iPad->GetX2() << ", Y1: " <<
    //    iPad->GetY1() << ", Y2: " << iPad->GetY2() << ", xr: " << xr << ", yr: " << yr <<
    //    ", x_User: " << x_User << ", y_User: " << y_User << endl;
}

void Convert_User_to_NDC(TPad* iPad, Double_t x_User, Double_t y_User, Double_t &x_NDC, Double_t &y_NDC)
{
    // Convert User coordinates to NDC coordinates
    Double_t xr = iPad->GetX2()-iPad->GetX1();
    Double_t yr = iPad->GetY2()-iPad->GetY1();
    x_NDC = (x_User-iPad->GetX1())/ xr;
    y_NDC = (y_User-iPad->GetY1())/ xr;
    cout << "X1: " << iPad->GetX2() << ", X2: " << iPad->GetX2() << endl;
}
//------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,
                      Float_t size=0.06,Int_t color=1,Float_t angle=0.0,
                      Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color
    // align: 1 left aligned, 32, right aligned

    // align = 10*HorizontalAlign + VerticalAlign
    // For horizontal alignment the following convention applies:
    // 1=left adjusted, 2=centered, 3=right adjusted
    // For vertical alignment the following convention applies:
    // 1=bottom adjusted, 2=centered, 3=top adjusted

    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_PHOS_in_eta_vs_phi(Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    PlotLine(-100.0,-40.0,-0.12,-0.12,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(-100.0,-40.0,0.12,0.12,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(-100.0,-100.0,-0.12,0.12,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(-40.0,-40.0,-0.12,0.12,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_phi_sector_TPC(Int_t sector, Int_t phi_bin, Int_t color)
{
    Int_t Line_Col  = color;
    Int_t LineWidth = 1;
    Int_t LineStyle = 1;

    vector<TVector2> vec_TV2_corner_point_local;
    vector<TVector2> vec_TV2_corner_point_global;
    vec_TV2_corner_point_local.resize(4);
    vec_TV2_corner_point_global.resize(4);

    Double_t y_over_x_value_min = (0.34/15.0)*phi_bin - 0.17;
    Double_t y_over_x_value_max = (0.34/15.0)*(phi_bin+1) - 0.17;
    Double_t Y_pad_position_min_lower = y_over_x_value_min*85.225;
    Double_t Y_pad_position_max_lower = y_over_x_value_max*85.225;
    Double_t Y_pad_position_min_upper = y_over_x_value_min*246.650;
    Double_t Y_pad_position_max_upper = y_over_x_value_max*246.650;

    vec_TV2_corner_point_local[0].SetX(Y_pad_position_min_lower); // lower left
    vec_TV2_corner_point_local[0].SetY(85.225); // lower left

    vec_TV2_corner_point_local[1].SetX(Y_pad_position_min_upper); // upper left
    vec_TV2_corner_point_local[1].SetY(246.650); // upper left

    vec_TV2_corner_point_local[2].SetX(Y_pad_position_max_upper); // upper right
    vec_TV2_corner_point_local[2].SetY(246.650); // upper right

    vec_TV2_corner_point_local[3].SetX(Y_pad_position_max_lower); // lower right
    vec_TV2_corner_point_local[3].SetY(85.225); // lower right


    Double_t phi_val = (-90.0 + (sector%18)*20.0 + 10.0)*TMath::DegToRad(); // rotation from local to global TPC coordinate system
    for(Int_t i_point = 0; i_point < 4; i_point++)
    {
        vec_TV2_corner_point_local[i_point].SetX(-1.0*vec_TV2_corner_point_local[i_point].X());
        vec_TV2_corner_point_global[i_point] = vec_TV2_corner_point_local[i_point].Rotate(phi_val);
        //printf("i_point: %d \n",i_point);
        //vec_TV2_corner_point_local[i_point].Print();
        //vec_TV2_corner_point_global[i_point].Print();
    }

    //printf(" \n");
    //printf("min_lower: %4.3f, max_lower: %4.3f, min_upper: %4.3f, max_upper: %4.3f \n",Y_pad_position_min_lower,Y_pad_position_max_lower,Y_pad_position_min_upper,Y_pad_position_max_upper);
    //printf("min_lower: %4.3f, max_lower: %4.3f, min_upper: %4.3f, max_upper: %4.3f \n",vec_TV2_corner_point_global[0].X(),vec_TV2_corner_point_global[3].X(),vec_TV2_corner_point_global[1].X(),vec_TV2_corner_point_global[2].X());


    PlotLine(vec_TV2_corner_point_global[0].X(),vec_TV2_corner_point_global[1].X(),vec_TV2_corner_point_global[0].Y(),vec_TV2_corner_point_global[1].Y(),Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(vec_TV2_corner_point_global[1].X(),vec_TV2_corner_point_global[2].X(),vec_TV2_corner_point_global[1].Y(),vec_TV2_corner_point_global[2].Y(),Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(vec_TV2_corner_point_global[2].X(),vec_TV2_corner_point_global[3].X(),vec_TV2_corner_point_global[2].Y(),vec_TV2_corner_point_global[3].Y(),Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(vec_TV2_corner_point_global[3].X(),vec_TV2_corner_point_global[0].X(),vec_TV2_corner_point_global[3].Y(),vec_TV2_corner_point_global[0].Y(),Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_R_Z_TPC(Int_t z_bin, Int_t top_bottom, Int_t color)
{
    Int_t Line_Col  = color;
    Int_t LineWidth = 3;
    Int_t LineStyle = 1;

    Float_t sign_top_bottom = 1.0;
    if(top_bottom != 0) sign_top_bottom = -1.0;

    Float_t sign_z = 1.0;
    if(z_bin < 5)
    {
        sign_z = -1.0;
        z_bin = 4 - z_bin;
    }
    if(z_bin >= 5)
    {
        z_bin = z_bin - 5;
    }

    Float_t X_pad_position_min = 85.225;
    Float_t X_pad_position_max = 246.650;
    Float_t Z_pad_position_lower_min = sign_z*(z_over_x_voxel_binning[z_bin]*X_pad_position_min);
    Float_t Z_pad_position_lower_max = sign_z*(z_over_x_voxel_binning[z_bin]*X_pad_position_max);
    Float_t Z_pad_position_upper_min = sign_z*(z_over_x_voxel_binning[z_bin+1]*X_pad_position_min);
    Float_t Z_pad_position_upper_max = sign_z*(z_over_x_voxel_binning[z_bin+1]*X_pad_position_max);

    //printf("Z_pad_position_lower_min: %4.3f, Z_pad_position_lower_max: %4.3f, X_pad_position_min: %4.3f, X_pad_position_max: %4.3f \n",Z_pad_position_lower_min,Z_pad_position_lower_max,X_pad_position_min,X_pad_position_max);
    PlotLine(Z_pad_position_lower_min,Z_pad_position_lower_max,sign_top_bottom*X_pad_position_min,sign_top_bottom*X_pad_position_max,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(Z_pad_position_lower_max,Z_pad_position_upper_max,sign_top_bottom*X_pad_position_max,sign_top_bottom*X_pad_position_max,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(Z_pad_position_upper_max,Z_pad_position_upper_min,sign_top_bottom*X_pad_position_max,sign_top_bottom*X_pad_position_min,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    PlotLine(Z_pad_position_upper_min,Z_pad_position_lower_min,sign_top_bottom*X_pad_position_min,sign_top_bottom*X_pad_position_min,Line_Col,LineWidth,LineStyle); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistLine2(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Float_t x_start, Float_t x_stop, Float_t y_min, Float_t y_max)
{
    TPolyLine* Hist_line = new TPolyLine();
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y >= y_min && y < y_max && y != 0 && x >= x_start && x <= x_stop)
        {
            //cout << "x = " << x << ", y = " << y << endl;
            Hist_line->SetNextPoint(x,y);
        }
    }
    Hist_line -> Draw();
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    // Modified March 14th, which year?
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y != 0 && x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line       = new TGraph(N_total_points);
    TGraph* Hist_line_upper = new TGraph(N_points);
    TGraph* Hist_line_lower = new TGraph(N_points);
    Hist_line ->SetLineWidth(LineWidth);
    Hist_line ->SetLineStyle(LineStyle);
    Hist_line ->SetLineColor(Line_Col);
    Hist_line ->SetFillStyle(FillStyle);
    Hist_line ->SetFillColor(FillColor);

    Hist_line_upper ->SetLineWidth(LineWidth);
    Hist_line_upper ->SetLineStyle(LineStyle);
    Hist_line_upper ->SetLineColor(Line_Col);
    Hist_line_lower ->SetLineWidth(LineWidth);
    Hist_line_lower ->SetLineStyle(LineStyle);
    Hist_line_lower ->SetLineColor(Line_Col);

    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            if(x >= x_start && x <= x_stop)
            {
                //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
                Hist_line->SetPoint(N_point,x,y-y_err);
                Hist_line_upper->SetPoint(N_point,x,y+y_err);
                Hist_line_lower->SetPoint(N_point,x,y-y_err);
                Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
                if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
                N_point++;
            }
        }
    }
    Hist_line       -> Draw("f");
    Hist_line_upper -> Draw("l");
    Hist_line_lower -> Draw("l");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotGraphErrorBand(TGraphAsymmErrors* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor, Float_t x_start, Float_t x_stop)
{
    Int_t N_points = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);

        if(x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }

    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line       = new TGraph(N_total_points);
    TGraph* Hist_line_upper = new TGraph(N_points);
    TGraph* Hist_line_lower = new TGraph(N_points);
    Hist_line ->SetLineWidth(LineWidth);
    Hist_line ->SetLineStyle(LineStyle);
    Hist_line ->SetLineColor(Line_Col);
    Hist_line ->SetFillStyle(FillStyle);
    Hist_line ->SetFillColor(FillColor);

    Hist_line_upper ->SetLineWidth(LineWidth);
    Hist_line_upper ->SetLineStyle(LineStyle);
    Hist_line_upper ->SetLineColor(Line_Col);
    Hist_line_lower ->SetLineWidth(LineWidth);
    Hist_line_lower ->SetLineStyle(LineStyle);
    Hist_line_lower ->SetLineColor(Line_Col);

    Int_t N_point = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);
        Double_t Xhigh = Histo->GetErrorXhigh(i_point);
        Double_t Xlow  = Histo->GetErrorXlow(i_point);
        Double_t Yhigh = Histo->GetErrorYhigh(i_point);
        Double_t Ylow  = Histo->GetErrorYlow(i_point);

        if(x >= x_start && x <= x_stop)
        {
            //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
            Hist_line->SetPoint(N_point,x,y-Ylow);
            Hist_line_upper->SetPoint(N_point,x,y+Yhigh);
            Hist_line_lower->SetPoint(N_point,x,y-Ylow);
            Hist_line->SetPoint(N_total_points-2-N_point,x,y+Yhigh);
            if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-Ylow);
            N_point++;
        }
    }
    Hist_line       -> Draw("f");
    Hist_line_upper -> Draw("l");
    Hist_line_lower -> Draw("l");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void PlotGraphErrorBandExt(TGraphAsymmErrors* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle,
                           Int_t FillColor, Float_t x_start, Float_t x_stop, Float_t transparency)
{
    Int_t N_points = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);

        if(x >= x_start && x <= x_stop)
        {
            N_points++;
        }
    }

    Int_t N_total_points = N_points*2+1;
    cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line       = new TGraph(N_total_points);
    TGraph* Hist_line_upper = new TGraph(N_points);
    TGraph* Hist_line_lower = new TGraph(N_points);
    Hist_line ->SetLineWidth(LineWidth);
    Hist_line ->SetLineStyle(LineStyle);
    Hist_line ->SetLineColor(Line_Col);
    Hist_line ->SetFillStyle(FillStyle);
    if(transparency > 1.0)
    {
        Hist_line ->SetFillColor(FillColor);
    }
    else
    {
        Hist_line ->SetFillColorAlpha(FillColor,transparency); // Make background transparent
    }

    Hist_line_upper ->SetLineWidth(LineWidth);
    Hist_line_upper ->SetLineStyle(LineStyle);
    Hist_line_upper ->SetLineColor(Line_Col);
    Hist_line_lower ->SetLineWidth(LineWidth);
    Hist_line_lower ->SetLineStyle(LineStyle);
    Hist_line_lower ->SetLineColor(Line_Col);

    Int_t N_point = 0;
    for(Int_t i_point = 0; i_point < Histo->GetN(); i_point++)
    {
        Double_t x,y;
        Histo->GetPoint(i_point,x,y);
        Double_t Xhigh = Histo->GetErrorXhigh(i_point);
        Double_t Xlow  = Histo->GetErrorXlow(i_point);
        Double_t Yhigh = Histo->GetErrorYhigh(i_point);
        Double_t Ylow  = Histo->GetErrorYlow(i_point);

        if(x >= x_start && x <= x_stop)
        {
            //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
            Hist_line->SetPoint(N_point,x,y-Ylow);
            Hist_line_upper->SetPoint(N_point,x,y+Yhigh);
            Hist_line_lower->SetPoint(N_point,x,y-Ylow);
            Hist_line->SetPoint(N_total_points-2-N_point,x,y+Yhigh);
            if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-Ylow);
            N_point++;
        }
    }
    Hist_line       -> Draw("f");
    Hist_line_upper -> Draw("l");
    Hist_line_lower -> Draw("l");
    //delete Hist_line;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t fEPD_Eta_vs_z(Double_t* z_val, Double_t* par)
{
    Double_t z_vertex, r, EPD_z_pos, flag_dir;
    r         = par[0];
    EPD_z_pos = par[1];
    flag_dir  = par[2];
    z_vertex  = z_val[0];

    Double_t theta = TMath::ATan(r/(EPD_z_pos-flag_dir*z_vertex));
    Double_t eta   = -TMath::Log(TMath::Tan(theta/2.0));
    return eta;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t ExpFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, Amplitude, T;
    Amplitude  = par[0];
    T          = par[1];
    x          = x_val[0];
    y          = Amplitude*TMath::Exp(-x/T);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    par5  = par[5];
    x = x_val[0];
    y = par0 + par1*x + par2*x*x + par3*x*x*x + par4*x*x*x*x + par5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t ImpactFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = fabs(par0) + par1*x + fabs(par2)*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PoissonFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    y = par0*TMath::Poisson(x-par1,par2);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LandauFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = par[0]; // amplitude
    par1  = par[1]; // mean
    par2  = par[2]; // sigma
    x = x_val[0];
    y = par0*TMath::Landau(x,par1,par2);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LandauGausFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, landau_amp, landau_mpv, landau_sigma, gaus_amp, gaus_mean, gaus_sigma;
    landau_amp    = par[0]; // landau amplitude
    landau_mpv    = par[1]; // ladau mpv
    landau_sigma  = par[2]; // landau sigma
    gaus_amp      = par[3]; // gaus amplitude
    gaus_mean     = par[4]; // gaus mean
    gaus_sigma    = par[5]; // gaus sigma
    x = x_val[0];
    y = landau_amp*TMath::Landau(x,landau_mpv,landau_sigma) + gaus_amp*TMath::Gaus(x,gaus_mean,gaus_sigma);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t One_over_x_FitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    y = 0.0;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    x = x_val[0];
    //if(x != 0.0 && !((par0/x) < 0.0 && par2 < 1.0)) y = pow(par0/x,par2) + par1;
    if(x != 0.0) y = par0*TMath::Power(1.0/x,par2) + par1;
    //if(x != 0.0) y =par0/x + par1 + 0.0001*par2;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TwoGaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4, par5;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    par3  = fabs(par[3]);
    par4  = par[4];
    par5  = fabs(par[5]);
    x = x_val[0];
    y = fabs(par0)*TMath::Gaus(x,par1,par2,0) + fabs(par3)*TMath::Gaus(x,par4,par5,0);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TwoGaussPedFitFunc(Double_t* x_val, Double_t* par)
{
    // With pedestal
    Double_t x, y, par0, par1, par2, par3, par4, par5, par6;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    par3  = fabs(par[3]);
    par4  = par[4];
    par5  = fabs(par[5]);
    par6  = par[6];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par3*TMath::Gaus(x,par4,par5,0) + par6;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, par3, par4;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    par3  = par[3];
    par4  = par[4];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + par4*(1.0 + 2.0*par3*TMath::Cos(2.0*x));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussPolyFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2, pol0, pol1, pol2, pol3, pol4, pol5;
    par0  = par[0];
    par1  = par[1];
    par2  = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussAreaPolyFitFunc(Double_t* x_val, Double_t* par)
{
    // From http://www.originlab.com/doc/Origin-Help/Gaussian-Function-FitFunc
    // Is using the area instead of amplitude
    // Area: area
    // xc: center
    // w: FWHM
    Double_t x, y, Area, xc, w, pol0, pol1, pol2, pol3, pol4, pol5;
    Area  = par[0];
    xc    = par[1];
    w     = par[2];
    pol0  = par[3];
    pol1  = par[4];
    pol2  = par[5];
    pol3  = par[6];
    pol4  = par[7];
    pol5  = par[8];
    x = x_val[0];
    y = Area*TMath::Exp(-4.0*TMath::Log(2)*TMath::Power(x-xc,2)/TMath::Power(w,2))/(w*TMath::Sqrt(TMath::Pi()/(4.0*TMath::Log(2)))) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t DoubleDiffGauss_dNdetaFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, A1, A2, s1, s2;
    A1 = par[0];
    A2 = par[1];
    s1 = par[2];
    s2 = par[3];
    x  = x_val[0];

    y = A1*TMath::Exp(-x*x/(2.0*s1*s1)) - A2*TMath::Exp(-x*x/(2.0*s2*s2));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t ThetaStarFunc(Double_t* x_val, Double_t* par)
{
    Double_t Theta, y, rho, A;
    rho   = par[0];
    A     = par[1];
    Theta = x_val[0];
    y = A*0.75*((1.0-rho) + (3.0*rho-1.0)*TMath::Power(TMath::Cos(Theta),2.0));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t CosThetaStarFunc(Double_t* x_val, Double_t* par)
{
    Double_t CosTheta, y, rho, A;
    rho      = par[0];
    A        = par[1];
    CosTheta = x_val[0];
    y = A*0.75*((1.0-rho) + (3.0*rho-1.0)*TMath::Power(CosTheta,2.0));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Impact_Angle_ADC_TRD_Func(Double_t* x_val, Double_t* par)
{
    Double_t Impact_angle, ADC, ADC_offset, drift_length;
    Impact_angle = x_val[0];
    drift_length = par[0];
    ADC_offset   = par[1];

    ADC = ADC_offset + drift_length/TMath::Cos(Impact_angle); // 3.0 cm drift length

    return ADC;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t FlowFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t phi, y, v0, v1, v2, v3, v4;
    v0  = par[0];
    v1  = par[1];
    v2  = par[2];
    v3  = par[3];
    v4  = par[4];
    phi = x_val[0];
    y = v0 * (1.0 + 2.0*v1*TMath::Cos(phi) + 2.0*v2*TMath::Cos(2.0*phi)
              + 2.0*v3*TMath::Cos(3.0*phi) + 2.0*v4*TMath::Cos(4.0*phi));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
    Double_t x, y, m0, Temp, Ampl, shift;
    m0    = par[0];
    Temp  = par[1];
    Ampl  = par[2];
    shift = par[3];
    x = x_val[0];
    y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp));
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t LevyFitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for d2N/(2pi*pT dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t LevyFitFunc_pT(Double_t* x_val, Double_t* par)
{
    // One over pT term is removed -> original pT distribution
    // Fit function for d2N/(2pi dpT dy)
    // taken from here: http://sampa.if.usp.br/~suaide/blog/files/papers/PLB6372006.pdf
    Double_t pT, y, B, T, n, m0;
    B    = par[0];
    T    = par[1];
    n    = par[2];
    m0   = par[3];
    pT   = x_val[0];
    Double_t mT = TMath::Sqrt(pT*pT+m0*m0);
    y = pT*B/TMath::Power(1.0+(mT-m0)/(n*T),n);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2 vs. pT
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT, a, b, c, d, n;
    pT = x_val[0];
    n  = par[0]; // number-of-constituent quarks
    a  = par[1];
    b  = par[2];
    c  = par[3];
    d  = par[4];

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_pT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. pT/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, pT_ncq, a, b, c, d, n;
    pT_ncq = x_val[0];
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT_ncq - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_ncq_FitFunc_Poly(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_ncq, a, b, c, d, n, m0, par6, par7;
    mT_ncq = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass
    par6   = par[6];
    par7   = par[7];

    Double_t mT = mT_ncq*n + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = (a/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d) + par6*mT_ncq + par7*mT_ncq*mT_ncq;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t v2_mT_FitFunc(Double_t* x_val, Double_t* par)
{
    // Fit function for v2/ncq vs. (mT-m0)/ncq
    // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
    Double_t v2, mT_m0, a, b, c, d, n, m0;
    mT_m0  = x_val[0]; // (mT-m0)/ncq
    n      = par[0]; // number-of-constituent quarks
    a      = par[1];
    b      = par[2];
    c      = par[3];
    d      = par[4];
    m0     = par[5]; // particle mass

    Double_t mT = mT_m0 + m0;
    Double_t pT = TMath::Sqrt(mT*mT-m0*m0);

    if(c != 0.0)
    {
        v2 = a*n/(1.0 + TMath::Exp(-(pT/n - b)/c)) - d*n;
    }
    else v2 = 0.0;

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncE(Double_t* x_val, Double_t* par)
{
    // Original
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        //for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            //r = r_start + j*delta_r;
            r = 0.01;

            //r = R; //

            rho = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_radial(Double_t* x_val, Double_t* par)
{
    // Blast wave function with radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Int_t nbins_r = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start = 0.0;
    Double_t delta_r = (par[4] - r_start)/nbins_r;
    Double_t r, R;
    T = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2 = par[3];
    R = par[4];

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi = phi_start + i*delta_phi;
        for(Int_t j = 0; j < nbins_r; j++)
        //for(Int_t j = 1; j < 2; j++)
        {
            //delta_r = 1.0; //

            r = r_start + j*delta_r;

            //r = R; //

            rho  = TMath::ATanH(TMath::TanH(rho_0)*r/R) + TMath::ATanH(TMath::TanH(rho_a)*r/R)*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc_no_mass_array(Double_t* x_val, Double_t* par)
{
    // Removed mass array
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Double_t Mass = par[5];
    Double_t mt = TMath::Sqrt(pt*pt + Mass*Mass);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 +=
            delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0
                                                                                         + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 +=
            delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 +
                                                                    2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[15] =
    {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677,9.46};
    //Double_t Mass[14] =
    //{1.019460,1.32131,1.32131,1.67245,1.67245,5.46,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,9.46};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    T = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 +=
            delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0
                                                                                         + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 +=
            delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 +
                                                                    2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------



/*
//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFunc(Double_t* x_val, Double_t* par)
{
    // Original function without radial dependence
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14]  = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt        = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    Int_t nbins_phi    = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;

    T          = par[0];
    rho_0      = par[1];
    rho_a      = par[2];
    s2         = par[3];
    Double_t R = par[4]; // to make it compatible with the R-dependent function

    Inte1 = 0.0;
    Inte2 = 0.0;

    for(Int_t i = 0; i < nbins_phi + 1; i++)
    {
        phi   = phi_start + i*delta_phi;
        rho   = rho_0 + rho_a*TMath::Cos(2.0*phi);
        alpha = (pt/T)*TMath::SinH(rho);
        beta  = (mt/T)*TMath::CosH(rho);

        Inte1 += delta_phi*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        Inte2 += delta_phi*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
    }

    if(Inte2 != 0)
    {
        v2 = Inte1/Inte2;
    }
    return v2;

}
//----------------------------------------------------------------------------------------
*/



//----------------------------------------------------------------------------------------
Double_t BesselFunc(Double_t* x_val, Double_t* par)
{
    Double_t x    = x_val[0];
    Double_t A    = par[0];
    Double_t par1 = par[1];
    Double_t par2 = par[2];
    Double_t y = A*TMath::BesselJ0(x)*(par1*x+par2);
    return y;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BesselFunc2D(Double_t* x_val, Double_t* par)
{
    Double_t x    = x_val[0] - par[3];
    Double_t y    = x_val[1] - par[4];
    Double_t r    = TMath::Sqrt(x*x + y*y);
    Double_t A    = par[0];
    Double_t par1 = par[1];
    Double_t par2 = par[2];
    Double_t z = A*TMath::BesselJ0(r)*(par1*r+par2);
    return z;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncC(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 100;
    Int_t nbins_r      = 100;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}
    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t BlastWaveFitFuncD(Double_t* x_val, Double_t* par)
{
    // Blast Wave Fit for v2
    Double_t pt, v2, alpha, beta, rho, rho_0, rho_a, phi, T, s2, Inte1, Inte2;
    v2 = 0.0;
    pt = x_val[0];
    Int_t PID = (Int_t)par[5];
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};
    Double_t NCQ[14]  = {2.0,3.0,3.0,3.0,3.0,3.0,3.0,2.0,3.0,3.0,2.0,2.0,2.0,2.0};
    pt /= NCQ[PID];
    Mass[PID] /= NCQ[PID];
    //pt *= NCQ[PID];
    Double_t mt = TMath::Sqrt(pt*pt + Mass[PID]*Mass[PID]);
    //cout << "mt = " << mt << endl;
    Int_t nbins_phi    = 25;
    Int_t nbins_r      = 25;
    Double_t phi_start = 0.0;
    Double_t phi_stop  = 2.0*TMath::Pi();
    Double_t delta_phi = (phi_stop - phi_start)/nbins_phi;
    Double_t r_start   = 0.0;
    Double_t r, R;
    T     = par[0];
    rho_0 = par[1];
    rho_a = par[2];
    s2    = par[3];
    R     = par[4];

    Double_t delta_r   = (R - r_start)/nbins_r;

    for(Int_t j = 0; j < nbins_r; j++)
    //for(Int_t j = 1; j < 2; j++)
    {
        r = r_start + j*delta_r;
        //r = R;

        Inte1 = 0.0;
        Inte2 = 0.0;

        for(Int_t i = 0; i < nbins_phi + 1; i++)
        {
            phi = phi_start + i*delta_phi;

            rho   = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
            alpha = (pt/T)*TMath::SinH(rho);
            beta  = (mt/T)*TMath::CosH(rho);

            //cout << "beta = " << beta << ", mt = " << mt << ", rho = " << rho << ", CosH = " << TMath::CosH(rho) << ", T = " << T << endl;

            Inte1 += delta_r*delta_phi*r*mt*TMath::Cos(2.0*phi)*TMath::BesselI(2,alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
            Inte2 += delta_r*delta_phi*r*mt*TMath::BesselI0(alpha)*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi));
        }

        if(Inte2 != 0)
        {
            v2 += Inte1/Inte2;
            //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
        }

    }

    v2 /= nbins_r;

    //if(Inte2 != 0)
    //{
    //    v2 = Inte1/Inte2;
        //cout << "PID = " << PID << ", pt = " << pt << ", mt = " << mt << ", v2 = " << v2 << endl;
    //}

    v2 *= NCQ[PID];
    //v2 /= NCQ[PID];
    //v2 *= (1.0/(1.0 + TMath::Exp((pt-4.0/NCQ[PID])/0.6)));

    return v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void BlastWaveSimultaneous(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
    Int_t npfits      = 0;
    Double_t chi2     = 0.0;
    Double_t Mass[14] = {1.019460,1.32131,1.32131,1.67245,1.67245,1.115683,1.115683,0.497648,0.938272,0.938272,0.13957,0.13957,0.493677,0.493677};

    for(Int_t i = 0; i < 14; i++) // loop over PIDs
    {
        p[5] = i; // PID
        //Double_t pt_cut = TMath::Sqrt((Mass[i]+mt_m0_cut)*(Mass[i]+mt_m0_cut) - Mass[i]*Mass[i]); // mt-m0 cut
        Double_t pt_cut = pt_cut_BW_global;
        if(flag_v2_BW_use[i] == 1)
        {
            for(Int_t ix = 0; ix < tgae_v2_stat_BW[i]->GetN(); ix++)
            {
                Double_t x[] = {tgae_v2_stat_BW[i]->GetX()[ix]};
                if(x[0] < pt_cut)
                {
                    Double_t y      = tgae_v2_stat_BW[i]->GetY()[ix];
                    Double_t ye     = tgae_v2_stat_BW[i]->GetErrorYhigh(ix);
                    //ye = 0.01;
                    //cout << "ix = " << ix << ", x = " << x[0] << ", y = " << y << ", ye = " << ye << endl;
                    Double_t bw_val = BlastWaveFitFunc(x,p);
                    //Double_t bw_val = 0.1;
                    //ye += 0.003;
                    Double_t diff   = (y - bw_val)/ye;
                    chi2 += diff*diff;
                    npfits++;
                }
                else break;
            }
        }
    }

    fval = chi2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Hist_interpolate_and_error(TH1F* hist, Double_t x, Double_t &Int_val, Double_t &Int_err)
{
    // Linear interpolation of a histogram
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[2]   = {0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetEntries() == 0) // no data poins to interpolate
    {
        return_val = -1;
        Int_val = 0;
        Int_err = 0;
        //cout << "No entries in histogram" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t binx = 1; binx < hist->GetNbinsX(); binx++)
        {
            Double_t bin_error_val   = hist->GetBinError(binx);
            Double_t bin_x_pos       = hist->GetBinCenter(binx);
            if(bin_error_val != 0)
            {
                err_counter++;
                bin_entries[1] = hist->GetBinContent(binx);
                bin_error[1]   = hist->GetBinError(binx);
                bin_x_val[1]   = hist->GetBinCenter(binx);
                if(bin_x_pos >= x)
                {
                    binx_high = binx;
                    flag_max  = 1;
                    break;
                }
                else flag_max = 0;
            }
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val = bin_entries[1];
            Int_err = bin_error[1];
            return return_val;
        }
        for(Int_t binx_low = binx_high; binx_low > 0; binx_low--)
        {
            bin_entries[0] = hist->GetBinContent(binx_low);
            bin_error[0]   = hist->GetBinError(binx_low);
            bin_x_val[0]   = hist->GetBinCenter(binx_low);
            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        if(bin_error[0] != 0 && bin_error[1] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;
                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[1];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;
                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val = bin_entries[0];
                Int_err = bin_error[0];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t TGraphAsymmErrors_interpolate_and_error(TGraphAsymmErrors* hist, Double_t x, Double_t &Int_val, Double_t &Int_err_low,Double_t &Int_err_high)
{
    // V2: 03.12.2012 -> bug fixed to calculate error bars
    // Linear interpolation of a TGraphAsymmErrors
    // Be careful, extrapolation is not done! -> closest data point is returned
    // Error calculation included
    // x-bin error is not taken into account

    Int_t return_val        = 0;
    Double_t bin_entries[2] = {0,0};
    Double_t bin_error[4]   = {0,0,0,0};
    Double_t bin_x_val[2]   = {0,0};
    Int_t binx_high;
    Int_t flag_max          = 0;

    if(hist->GetN() == 0) // no data poins to interpolate
    {
        return_val   = -1;
        Int_val      = 0;
        Int_err_low  = 0;
        Int_err_high = 0;
        //cout << "No entries in TGraphAsymmErrors" << endl;
        return 0;
    }
    else
    {
        Int_t err_counter = 0;
        for(Int_t epoint = 0; epoint < hist->GetN(); epoint++)
        {
            hist->GetPoint(epoint,bin_x_val[1],bin_entries[1]);
            bin_error[2] = hist->GetErrorYlow(epoint);
            bin_error[3] = hist->GetErrorYhigh(epoint);

            err_counter++;
            if(bin_x_val[1] >= x)
            {
                binx_high = epoint;
                flag_max  = 1;
                break;
            }
            else flag_max = 0;
        }
        if(err_counter == 1 || flag_max == 0) // There is no lower/uppper data point -> extrapolation, return closest values
        {
            return_val = 0;
            Int_val      = bin_entries[1];
            Int_err_low  = bin_error[2];
            Int_err_high = bin_error[3];
            return return_val;
        }
        for(Int_t epoint = binx_high; epoint >= 0; epoint--)
        {
            hist->GetPoint(epoint,bin_x_val[0],bin_entries[0]);
            bin_error[0] = hist->GetErrorYlow(epoint);
            bin_error[1] = hist->GetErrorYhigh(epoint);

            if(bin_x_val[0] < x && bin_error[0] != 0)
            {
                break;
            }
        }

        //cout << "bin0 = " << bin_error[0] << ", bin2 = " << bin_error[2] << endl;

        if(bin_error[0] != 0 && bin_error[2] != 0)
        {
            return_val = 1;
            if(bin_x_val[0] != bin_x_val[1])
            {
                Double_t slope = (bin_entries[1]-bin_entries[0])/(bin_x_val[1]-bin_x_val[0]);
                Double_t t_val = bin_entries[1]-slope*bin_x_val[1];
                Int_val        = slope*x+t_val;

                Double_t x1     = bin_x_val[0];
                Double_t x2     = bin_x_val[1];
                Double_t y1     = bin_entries[0];
                Double_t y2     = bin_entries[1];
                Double_t x1_err = 0.0;
                Double_t x2_err = 0.0;

                //cout << "x1 = " << x1 << ", x2 = " << x2 << ", y1 = " << y1 << ", y2 = " << y2 << endl;

                Double_t y1_err = bin_error[0];
                Double_t y2_err = bin_error[2];

                Double_t termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                Double_t termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                Double_t termC  = -(x-x2)/(x2-x1);
                Double_t termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_low = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);

                y1_err = bin_error[1];
                y2_err = bin_error[3];

                termA  = ((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2);
                termB  = -((y2-y1)/((x2-x1)*(x2-x1)))*(x-x2) - (y2-y1)/(x2-x1);
                termC  = -(x-x2)/(x2-x1);
                termD  = (x-x2)/(x2-x1) + 1.0;

                termA *= x1_err;
                termB *= x2_err;
                termC *= y1_err;
                termD *= y2_err;
                Int_err_high = TMath::Sqrt(termA*termA+termB*termB+termC*termC+termD*termD);
                return return_val;
            }
            else
            {
                Int_val      = bin_entries[0];
                Int_err_low  = bin_error[0];
                Int_err_high = bin_error[1];
                return return_val;
            }
        }
    }

    return return_val;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Get_T_muB_from_SHM(Double_t sqrt_sNN, Double_t &T, Double_t &muB)
{
    // Calculates baryon chemical potential muB and chemical freeze-out temperature T
    // based on fits to statistical hardronization model (SHM) calculations
    // from here: http://arxiv.org/pdf/hep-ph/0511094.pdf, Phys.Rev. C73 (2006) 034905
    const Double_t a = 0.166; // +/- 0.002
    const Double_t b = 0.139; // +/- 0.016
    const Double_t c = 0.053; // +/- 0.021
    const Double_t d = 1.308; // +/- 0.028
    const Double_t e = 0.273; // +/- 0.008

    if(sqrt_sNN > 0.0)
    {
        muB = d/(1.0 + e*sqrt_sNN);
        T   = a - b*muB*muB - c*muB*muB*muB*muB;
    }
    else
    {
        muB = -999.0;
        T   = -999.0;
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Get_muB(Double_t sqrt_sNN, Int_t flag_input, Double_t &muB_err)
{
    // Returns baryon chemical potential as a function of the center-of-mass energy sqrt(sNN)
    // from http://arxiv.org/pdf/1111.2406v1.pdf
    // flag_input: 0 = parametrization, 1 = STAR fits, 2 = corrected parametrization for 0-80% (see below)

    // corrected parametrization: Centrality data taken from: http://arxiv.org/pdf/0808.2041.pdf
    Double_t centrality_vals[9] = {75,65,55,45,35,25,15,7.5,2.5};
    Double_t muB_Energy[2][2][9] = // [62.4,200][values,error][70-80%,...,0-5%]
    {
        {   // Au+Au @ 62.4 GeV
            {37.7,42.5,47.0,51.3,54.2,54.5,59.4,61.0,62.7}, // values
            {6.5,5.8,5.1,5.2,5.2,5.2,5.4,5.7,6.0}  // errors
        },
        {   // Au+Au @ 200 GeV
            {14.1,15.3,17.7,18.9,18.6,21.3,21.0,22.8,21.9}, // values
            {4.2,4.2,4.2,4.2,4.2,4.2,4.2,4.5,4.5}  // errors
        }
    };

    Double_t muB     = 0.0;
    muB_err          = 0.0;
    Double_t a       = 1.482;  // +/- 0.0037 GeV
    Double_t b       = 0.3517; // +/- 0.009 GeV^-1

    if(!(flag_input == 0 || flag_input == 1 || flag_input == 2))
    {
        cout << "WARNING: Get_muB function, flag_input is wrong. Either 0 or 1. Parametrization is used." << endl;
        flag_input = 0;
    }

    if(sqrt_sNN >= 0.0 && (flag_input == 0 || flag_input == 2))
    {
        muB = a/(1.0+b*sqrt_sNN);
    }
    if(flag_input == 2) // correct the muB values for 0-80%
    {
        Double_t scale_fac[2][2];  // [62.4,200][val,error]
        Double_t muB_mean[2][2]; // [62.4,200][val,error]
        for(Int_t ebeam = 0; ebeam < 2; ebeam++)
        {
            for(Int_t val_err = 0; val_err < 2; val_err++) // calculate mean muB + error
            {
                muB_mean[ebeam][val_err] = 0.0; // approximation for 0-80%
                Double_t total_weight = 0.0;
                for(Int_t cent = 0; cent < 9; cent++)
                {
                    Double_t weight = 1.0;
                    if(cent >= 7) weight = 0.5;
                    if(val_err == 0) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent];
                    if(val_err == 1) muB_mean[ebeam][val_err] += weight*muB_Energy[ebeam][val_err][cent]*weight*muB_Energy[ebeam][val_err][cent];
                    total_weight += weight;
                }
                if(total_weight > 0.0)
                {
                    if(val_err == 0) muB_mean[ebeam][val_err] /= total_weight;
                    if(val_err == 1)
                    {
                        muB_mean[ebeam][val_err] = TMath::Sqrt(muB_mean[ebeam][val_err]);
                        muB_mean[ebeam][val_err] /= total_weight;
                    }
                }
            }
            scale_fac[ebeam][0] = muB_mean[ebeam][0]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0); // <muB>/0-10% = 0-80%/0-10%
            Double_t term1 = muB_mean[ebeam][1]/((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7])/2.0);
            Double_t term2 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][8]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            Double_t term3 = muB_mean[ebeam][0]*muB_Energy[ebeam][1][7]*2.0/TMath::Power((muB_Energy[ebeam][0][8]+muB_Energy[ebeam][0][7]),2);
            scale_fac[ebeam][1] = TMath::Sqrt(term1*term1+term2*term2+term3*term3); // error on scaling factor
        }
        Double_t mean_scale_fac[2] = {(scale_fac[0][0]+scale_fac[1][0])/2.0,TMath::Sqrt(scale_fac[0][1]*scale_fac[0][1]+scale_fac[1][1]*scale_fac[1][1])/2.0}; // [value,error]
        Double_t original_muB = muB;
        muB     = muB * mean_scale_fac[0]; // 0-80% estimation
        muB_err = muB * mean_scale_fac[1]; // 0-80% error estimation
        //cout << "muB = " << original_muB << ", muB(0-80%) = " << muB << ", muB_err(0-80%) = " << muB_err << ", <scale> = " << mean_scale_fac[0] << ", <scale_err> = " << mean_scale_fac[1] << endl;
        //cout << "scale_fac(62 GeV) = " << scale_fac[0][0] << ", scale_fac(200 GeV) = " << scale_fac[1][0] << endl;
        //cout << "scale_fac_err(62 GeV) = " << scale_fac[0][1] << ", scale_fac_err(200 GeV) = " << scale_fac[1][1] << endl;
    }

    if(flag_input == 1)
    {
        // Preliminary muB parameters from Lokesh with strange particles for 0-80% Au+Au collisions
        // 7.7 :  3.61377e-01   1.51360e-02
        // 11.5:  2.59978e-01   1.06312e-02
        // 39  :  8.81846e-02   4.73062e-03
        // 200 :  22 -> Thats a guess

        // Fittet in TmuB_fit.cc macro with One_over_x_FitFunc function
        Double_t par0 = 2.13495e+00;
        Double_t par1 = 4.35916e-04;
        Double_t par2 = 8.67894e-01;
        muB = par0*TMath::Power(1.0/sqrt_sNN,par2) + par1;
    }

    return muB;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_new_SymbolExt(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size, Float_t line_width)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->SetLineWidth(line_width);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_hist_new_Symbol(TH1D* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TH1D* ge_clone_A = (TH1D*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TH1D* ge_clone_B = (TH1D*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TH1D* ge_clone_C = (TH1D*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err,
                                Double_t y_min_err, Double_t y_max_err,
                                Int_t style, Int_t color, Float_t size)
{
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);

    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.5*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_header(FILE* HTML_file)
{
    fprintf(HTML_file,"%s \n","<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">");
    fprintf(HTML_file,"%s \n","<html>");
    fprintf(HTML_file,"%s \n","<head>");


    fprintf(HTML_file,"%s \n","  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">");
    fprintf(HTML_file,"%s \n","  <title>STAR Physics Database</title>");
    fprintf(HTML_file,"%s \n","</head>");
    fprintf(HTML_file,"%s \n"," <body bgcolor=\"white\">");

    fprintf(HTML_file,"%s \n","<h2><span class=\"normaltext\"><font color=\"#ff8c00\" face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\">");
    fprintf(HTML_file,"%s \n","Observation of an energy-dependent difference in elliptic flow between particles and anti-particles in relativistic heavy ion collisions");
    fprintf(HTML_file,"%s \n","</font>");

    fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\"><br> <br> </font></span></h2>");

    //fprintf(HTML_file,"%s \n","<hr><h3><font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">Figure 1");
    //fprintf(HTML_file,"%s \n","</font></h3>");
    //fprintf(HTML_file,"%s \n","<font face=\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\" color=\"#0000cd\">");
    //fprintf(HTML_file,"%s \n","<A HREF=\"fig1.png\"> png </A> | <A HREF=\"fig1.eps\"> eps </A> </br>");
    fprintf(HTML_file,"%s \n","</font>");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_title(TString Label, FILE* HTML_file)
{
    TString label_html = "<hr><h3><font face=\\\"Arial,Helvetica,Geneva,Swiss,SunSans-Regular\\\" color=\\\"#0000cd\\\">";
    label_html += Label.Data();
    label_html += " </font> </h3> </hr>";
    //label_html += "\"";
    fprintf(HTML_file,"%s \n",label_html.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Write_HTML_tables(TString labelX, TString labelY, TGraphAsymmErrors* data_stat, TGraphAsymmErrors* syst_errors, TGraphAsymmErrors* glob_syst_errors, TString Label, FILE* HTML_file)
{
    fprintf(HTML_file,"%s %s %s \n","<tr><td>",Label.Data(),"</td> </tr>");
    TString label_html = "<table border=1> <tr><th> ";
    label_html += labelX.Data();
    label_html += " </th> <th> ";
    label_html += labelY.Data();
    label_html += " </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>";

    //fprintf(HTML_file,"%s \n","<table border=1> <tr><th> p<sub>T</sub> (GeV/c) </th> <th> v<sub>2</sub> </th> <th> stat. err. </th> <th> syst. low </th> <th> syst. high </th> <th> syst. glob. </th> </tr>");

    fprintf(HTML_file,"%s \n",label_html.Data());
    Int_t flag_glob_syst = 1;
    if(glob_syst_errors == syst_errors) flag_glob_syst = -1;

    for(Int_t epoint = 0; epoint < data_stat->GetN(); epoint++)
    {
        Double_t x_val, y_val, x_syst, stat_error, low_syst, high_syst, glob_syst;
        data_stat                     ->GetPoint(epoint,x_val,y_val);
        stat_error = data_stat        ->GetErrorYhigh(epoint);
        high_syst  = syst_errors      ->GetErrorYhigh(epoint);
        low_syst   = syst_errors      ->GetErrorYlow(epoint);
        glob_syst  = glob_syst_errors ->GetErrorYhigh(0);
        if(flag_glob_syst == -1) glob_syst = 0.0;

        fprintf(HTML_file,"%s %f %s %f %s %f %s %f %s %f %s %f %s \n","<tr> <td> ",x_val," </td> <td> ",y_val," </td> <td> ",stat_error," </td> <td> ",
                low_syst," </td> <td> ",high_syst," </td> <td> ",glob_syst,"</td></tr>");

    }
    fprintf(HTML_file,"%s \n","</table><hr><table border=1>");
}
//----------------------------------------------------------------------------------------


#if 0
//----------------------------------------------------------------------------------------
class Calc_Write_HEP_data
{
    TH2D* h2D_map;
    std::vector<Double_t> vec_height;
    std::vector<TPolyMarker3D*> vec_pm;
    Double_t x_val_min, x_val_max, y_val_min, y_val_max;
public:
    void label_dataset(TString lab_dataset);
    void label_dscomment(TString lab_dscomment);
    void label_reackey(TString lab_reackey);
    void clear();
};

void Calc_Map_height::add_map(TH2D* h2D_map_in)
{
    h2D_map = (TH2D*)h2D_map_in->Clone("h2D_map");
}
//----------------------------------------------------------------------------------------
# endif


//----------------------------------------------------------------------------------------
Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the minimum distance vector of a point to a straight given as parametric straight x = base + n * dir

    TVector3 diff = base-point;
    TVector3 dir_norm = dir;
    dir_norm *= (1.0/dir.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir - diff;

    return dist_vec;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_det_index(Int_t sensor_id)
{
    // Determines detector index (see definition below) based on the sonsor id
    // 0 = inner pixel left
    // 1 = outer pixel left
    // 2 = IST left
    // 3 = inner pixel right
    // 4 = outer pixel right
    // 5 = IST right

    Int_t HFT_det_index = -1;
    if(
       (sensor_id >= 1   && sensor_id <= 10)  ||
       (sensor_id >= 41  && sensor_id <= 50)  ||
       (sensor_id >= 81  && sensor_id <= 90)  ||
       (sensor_id >= 121 && sensor_id <= 130) ||
       (sensor_id >= 161 && sensor_id <= 170)
      )
    {
        HFT_det_index = 0;
    }

    if(
       (sensor_id >= 11   && sensor_id <= 40)  ||
       (sensor_id >= 51   && sensor_id <= 80)  ||
       (sensor_id >= 91   && sensor_id <= 120) ||
       (sensor_id >= 131  && sensor_id <= 160) ||
       (sensor_id >= 171  && sensor_id <= 200)
      )
    {
        HFT_det_index = 1;
    }
    if(
       sensor_id >= 1001 && sensor_id <= 1072
      )
    {
        HFT_det_index = 2;
    }
    if(
       (sensor_id >= 1+200   && sensor_id <= 10+200)  ||
       (sensor_id >= 41+200  && sensor_id <= 50+200)  ||
       (sensor_id >= 81+200  && sensor_id <= 90+200)  ||
       (sensor_id >= 121+200 && sensor_id <= 130+200) ||
       (sensor_id >= 161+200 && sensor_id <= 170+200)
      )
    {
        HFT_det_index = 3;
    }
    if(
       (sensor_id >= 11+200   && sensor_id <= 40+200)  ||
       (sensor_id >= 51+200   && sensor_id <= 80+200)  ||
       (sensor_id >= 91+200   && sensor_id <= 120+200) ||
       (sensor_id >= 131+200  && sensor_id <= 160+200) ||
       (sensor_id >= 171+200  && sensor_id <= 200+200)
      )
    {
        HFT_det_index = 4;
    }
    if(
       sensor_id >= 1073 && sensor_id <= 1145
      )
    {
        HFT_det_index = 5;
    }
    return HFT_det_index;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_histogram(TH1D* hist, TString plot_option, TString label_x, TString label_y, Double_t x_start, Double_t x_stop,
                    Double_t y_start, Double_t y_stop, Double_t title_offset_x, Double_t title_offset_y, Double_t label_offset_x,
                    Double_t label_offset_y, Double_t label_size, Double_t title_size, Int_t N_div)
{
    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(title_offset_x);
    hist->GetYaxis()->SetTitleOffset(title_offset_y);
    hist->GetXaxis()->SetLabelOffset(label_offset_x);
    hist->GetYaxis()->SetLabelOffset(label_offset_y);
    hist->GetXaxis()->SetLabelSize(label_size);
    hist->GetYaxis()->SetLabelSize(label_size);
    hist->GetXaxis()->SetTitleSize(title_size);
    hist->GetYaxis()->SetTitleSize(title_size);
    hist->GetXaxis()->SetNdivisions(N_div,'N');
    hist->GetYaxis()->SetNdivisions(N_div,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitle(label_x.Data());
    hist->GetYaxis()->SetTitle(label_y.Data());
    if(!(x_start == 0 && x_stop == 0)) hist->GetXaxis()->SetRangeUser(x_start,x_stop);
    if(!(y_start == 0 && y_stop == 0)) hist->GetYaxis()->SetRangeUser(y_start,y_stop);
    hist->DrawCopy(plot_option.Data());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t get_HFT_pixel_sector(Int_t sensor_id, Int_t &IST_ladder)
{
    // Determines pixel sector based on the sonsor id

    // PXL
    Int_t HFT_pxl_sector = -1;
    if(sensor_id >= 1 && sensor_id <= 400)
    {
        HFT_pxl_sector = (Int_t)((sensor_id-1)/40);
    }

    // IST -> 24 ladders, 6 sensors per ladder -> 144 sensors in total, starting from sensor_id = 1001
    if(sensor_id >= 1001 && sensor_id <= 1144)
    {
        IST_ladder = (Int_t)((sensor_id-1001)/6);
    }
    else IST_ladder = -1;

    return HFT_pxl_sector;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_STAR_3D()
{

    TGeoManager *geom = new TGeoManager("geom","My 3D Project");
    //------------------Creat materials------------------------------
    TGeoMaterial *vacuum = new TGeoMaterial("vacuum",0,0,0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.84,26.7,7.87);
    Fe->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_outer_tube = new TGeoMaterial("M_outer_tube",55.84,26.7,7.87);
    M_outer_tube->SetTransparency(93); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_IDS = new TGeoMaterial("M_IDS",55.84,26.7,7.87);
    M_IDS       ->SetTransparency(80); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_beampipe = new TGeoMaterial("M_beampipe",55.84,26.7,7.87);
    M_beampipe       ->SetTransparency(70); // higher value means more transparent, 100 is maximum

    TGeoMaterial *M_Pixel_support = new TGeoMaterial("M_Pixel_support",55.84,26.7,7.87);
    M_Pixel_support    ->SetTransparency(70); // higher value means more transparent, 100 is maximum


    //------------------Create media---------------------------------
    TGeoMedium *Air = new TGeoMedium("Air",0,vacuum);
    TGeoMedium *Iron = new TGeoMedium("Iron",1,Fe);
    TGeoMedium *Me_outer_tube = new TGeoMedium("Me_outer_tube",1,M_outer_tube);
    TGeoMedium *Me_IDS        = new TGeoMedium("Me_IDS",1,M_IDS);
    TGeoMedium *Me_beampipe   = new TGeoMedium("Me_beampipe",1,M_beampipe);
    TGeoMedium *Me_Pixel_support   = new TGeoMedium("Me_Pixel_support",1,M_Pixel_support);

    //------------------Create TOP volume----------------------------
    TGeoVolume *top = geom->MakeBox("top",Air,500,500,500);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);
    // If you want to see the boundary, please input the number, 1 instead of 0.
    // Like this, geom->SetTopVisible(1);


    TGeoVolume *inner_field_tube       = geom->MakeTube("inner_field_tube",Iron,49.5,50.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *outer_field_tube       = geom->MakeTube("outer_field_tube",Me_outer_tube,199.5,200.0,200.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_central_part       = geom->MakeTube("IDS_central_part",Me_IDS,42.8/2.0,43.0/2.0,56.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_side_parts         = geom->MakeTube("IDS_side_parts",Me_IDS,79.3/2.0,79.5/2.0,(222.7-64.0)/2.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *IDS_connection_parts_R = geom->MakeCone("IDS_connection_parts_R",Me_IDS,(64.0-56.0)/2.0,42.8/2.0,43.0/2.0,79.3/2.0,79.5/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *IDS_connection_parts_L = geom->MakeCone("IDS_connection_parts_L",Me_IDS,(64.0-56.0)/2.0,79.3/2.0,79.5/2.0,42.8/2.0,43.0/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *beampipe_central_part       = geom->MakeTube("beampipe_central_part",Me_beampipe,4.05/2.0,4.15/2.0,141.5);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_side_parts         = geom->MakeTube("beampipe_side_parts",Me_beampipe,9.52/2.0,9.62/2.0,100.0);  // r_min, r_max, dz (half of total length)
    TGeoVolume *beampipe_connection_parts_R = geom->MakeCone("beampipe_connection_parts_R",Me_beampipe,(191.5-141.5)/2.0,4.05/2.0,4.15/2.0,9.52/2.0,9.62/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max
    TGeoVolume *beampipe_connection_parts_L = geom->MakeCone("beampipe_connection_parts_L",Me_beampipe,(191.5-141.5)/2.0,9.52/2.0,9.62/2.0,4.05/2.0,4.15/2.0); // dz (half of total length), r1_min, r1_max, r2_min, r2_max

    TGeoVolume *Pixel_support       = geom->MakeTube("Pixel_support",Me_Pixel_support,21.8/2.0,22.0/2.0,56.0);  // r_min, r_max, dz (half of total length)

    inner_field_tube       ->SetLineColor(4);
    outer_field_tube       ->SetLineColor(kRed-8);
    IDS_central_part       ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_side_parts         ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_R ->SetLineColor(2);  // Inner Detector Support (IDS)
    IDS_connection_parts_L ->SetLineColor(2);  // Inner Detector Support (IDS)

    beampipe_central_part       ->SetLineColor(3);  // (beampipe)
    beampipe_side_parts         ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_R ->SetLineColor(3);  // (beampipe)
    beampipe_connection_parts_L ->SetLineColor(3);  // (beampipe)

    Pixel_support ->SetLineColor(kYellow-3);  // (pixel support)

    top->AddNodeOverlap(inner_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(outer_field_tube,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,64.0+(222.7-64.0)/2.0));
    top->AddNodeOverlap(IDS_side_parts,1,new TGeoTranslation(0,0,-(64.0+(222.7-64.0)/2.0)));
    top->AddNodeOverlap(IDS_connection_parts_R,1,new TGeoTranslation(0,0,56.0+(64.0-56.0)/2.0));
    top->AddNodeOverlap(IDS_connection_parts_L,1,new TGeoTranslation(0,0,-(56.0+(64.0-56.0)/2.0)));

    top->AddNodeOverlap(beampipe_central_part,1,new TGeoTranslation(0,0,0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,191.5+100.0));
    top->AddNodeOverlap(beampipe_side_parts,1,new TGeoTranslation(0,0,-(191.5+100.0)));
    top->AddNodeOverlap(beampipe_connection_parts_R,1,new TGeoTranslation(0,0,141.4+(191.5-141.5)/2.0));
    top->AddNodeOverlap(beampipe_connection_parts_L,1,new TGeoTranslation(0,0,-(141.4+(191.5-141.5)/2.0)));

    top->AddNodeOverlap(Pixel_support,1,new TGeoTranslation(0,0,0));

    top->DrawClone("ogl");



    const Int_t n_TPC_points = 50;
    TPolyLine3D   *TPC_endcaps[4];
    TPolyLine3D   *TPC_tube[4];
    TPolyLine3D   *TPC_tube_lines[n_TPC_points+1];

    Float_t radius_table[4] = {200,200,3.81,3.81};
    Float_t z_val_table[4]  = {200,-200,200,-200};

    Float_t radius_table_tube[4] = {50,50,50,50};
    Float_t z_val_table_tube[4]  = {200,-200,100,-100};

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_endcaps[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_endcaps[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
        }
        TPC_endcaps[r]->SetLineStyle(0);
        TPC_endcaps[r]->SetLineColor(28); // 28
        TPC_endcaps[r]->SetLineWidth(2);
        TPC_endcaps[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < 4; r++)
    {
        TPC_tube[r] = new TPolyLine3D();
        Float_t radius   = radius_table_tube[r];
        Float_t x_offset = 0.0;
        Float_t y_offset = 0.0;
        Float_t z_tpc_val   = z_val_table_tube[r];
        for(Int_t t = 0; t < n_TPC_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_TPC_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            TPC_tube[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
            if(r == 0 && (t%4 == 0))
            {
                TPC_tube_lines[t] = new TPolyLine3D();
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_tpc_val);
                TPC_tube_lines[t]->SetNextPoint(x_tpc_val,y_tpc_val,z_val_table_tube[r+1]);
                TPC_tube_lines[t]->SetLineStyle(0);
                TPC_tube_lines[t]->SetLineColor(28); // 28
                TPC_tube_lines[t]->SetLineWidth(1);
                //TPC_tube_lines[t]->DrawClone("ogl");
            }
        }
        TPC_tube[r]->SetLineStyle(0);
        TPC_tube[r]->SetLineColor(28); // 28
        TPC_tube[r]->SetLineWidth(2);
        TPC_tube[r]->DrawClone("ogl");
    }

    TPolyLine3D   *BeamLine;
    BeamLine       = new TPolyLine3D(2);
    BeamLine   ->SetPoint(0,0,0,-550);
    BeamLine   ->SetPoint(1,0,0,550);
    BeamLine   ->SetLineStyle(0);
    BeamLine   ->SetLineColor(4);
    BeamLine   ->SetLineWidth(2);
    BeamLine   ->DrawClone("ogl");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TH1F* Modify_hist(TH1F* hist_input, Double_t mean_shift, Double_t width_fac, Double_t width_shift, Long64_t sample_val)
{
    // Modifies an input histogram by scaling it with width_fac (with scaling center width_shift), shifting it with mean_shift
    // sample val is a threshold at which a pure sampling of bin contents is changed to a sampling with max of sample_val and weights
    // this is important if the histogram spans over many magnitudes -> pure sampling would take too much time

    TString histname = hist_input->GetName();
    histname += "_out";
    TH1F* hist_output = (TH1F*)hist_input->Clone(histname.Data());
    hist_output->Reset();

    TRandom3 r3b_hist;
    r3b_hist.SetSeed(0); // seed for random number generator changes every second

    // Get Integral
    Double_t integral_input = hist_input ->Integral(1,hist_input->GetNbinsX());

    // Get lowest bin content
    Double_t bin_cont_min = hist_input ->GetBinContent(hist_input->GetMaximumBin());
    for(Int_t ibin = 1; ibin <= hist_input->GetNbinsX(); ibin++)
    {
        Double_t bin_cont  = hist_input ->GetBinContent(ibin);
        if(bin_cont > 0.0 && bin_cont < bin_cont_min) bin_cont_min = bin_cont;
    }

    for(Int_t ibin = 1; ibin <= hist_input->GetNbinsX(); ibin++)
    {
        Double_t bin_cont  = hist_input ->GetBinContent(ibin);
        Double_t bin_width = hist_input ->GetBinWidth(ibin);
        Double_t bin_cent  = hist_input ->GetBinCenter(ibin);

        Double_t weight_sample = 1.0;
        if(bin_cont_min > 0.0) bin_cont /= bin_cont_min;
        if(bin_cont > sample_val && sample_val > 0)
        {
            weight_sample = ((Double_t)bin_cont)/((Double_t)sample_val);
            bin_cont      = sample_val;
        }

        for(Long64_t iloop = 0; iloop < (Long64_t)bin_cont; iloop++)
        {
            Double_t random = r3b_hist.Rndm();
            random *= bin_width;
            random += bin_cent - (bin_width/2.0) - width_shift;
            random *= width_fac;
            random += mean_shift + width_shift;
            hist_output ->Fill(random,weight_sample);
        }
    }

    if(bin_cont_min > 0.0) hist_output->Scale(bin_cont_min);

    return hist_output;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_Detector_2D(Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1,
                             const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                             Float_t x_offset = 0.0, Float_t y_offset = 0.0
                            )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine   *tp_Circles[n_radii];
    TPolyLine   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColor(color); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        tp_Circles[r]->DrawClone("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(line_style);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(line_width);
        tp_Radial[r]->DrawClone("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_2D(Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1,
                             const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                             Float_t x_offset = 0.0, Float_t y_offset = 0.0, Int_t fill_color = 1, Double_t transparency = 0.0
                            )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine   *tp_Circles[n_radii];
    TPolyLine   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            //cout << ", x: " << x_tpc_val << ", y: " << y_tpc_val << endl;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }

        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColorAlpha(color,transparency); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        tp_Circles[r]->SetFillColorAlpha(fill_color,transparency);
        //tp_Circles[r]->DrawClone("oglf");
        tp_Circles[r]->Draw("f");
        tp_Circles[r]->Draw("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->Draw("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_2D_new(Float_t radius_in = 1.0, Float_t radius_out = 2.0,const Int_t n_radii = 1,
                        const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                        Float_t x_offset = 0.0, Float_t y_offset = 0.0, Int_t fill_color = 1,
                        Double_t transparency_line = 0.0, Double_t transparency_fill = 0.0
                       )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine   *tp_Circles[n_radii];
    TPolyLine   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            //cout << ", x: " << x_tpc_val << ", y: " << y_tpc_val << endl;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }

        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColorAlpha(color,transparency_line); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        tp_Circles[r]->SetFillColorAlpha(fill_color,transparency_fill);
        //tp_Circles[r]->DrawClone("oglf");
        tp_Circles[r]->Draw("f");
        tp_Circles[r]->Draw("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->Draw("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_Ring(Float_t radius_in = 1.0, Float_t radius_out = 2.0, Float_t color = 2,
                      Int_t line_style = 1, Int_t line_width = 1,
                      Float_t x_offset = 0.0, Float_t y_offset = 0.0, Int_t fill_color = 1,
                      Double_t transparency_line = 0.0, Double_t transparency_fill = 0.0,
                      Double_t phi_start = 0.0, Double_t phi_stop = 2.0*TMath::Pi()
                     )
{
    const Int_t n_radii  = 2;
    const Int_t n_points = 50;
    TPolyLine* tp_Circles_fill;
    vector<TPolyLine*> tp_Circles_lines;
    tp_Circles_lines.resize(2); // inner, outer
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));
    Double_t delta_phi = fabs(phi_stop - phi_start);
    if(delta_phi <= 0.0) delta_phi = 2.0*TMath::Pi();

    tp_Circles_fill = new TPolyLine();
    Double_t x_first, y_first;
    Double_t phi_val_last = 0.0;
    for(Int_t r = 0; r < n_radii; r++)
    {
        tp_Circles_lines[r] = new TPolyLine();
        radius_table[r] = radius_in + r*delta_radius;
        Float_t radius   = radius_table[r];
        if(r == 1) phi_start = phi_val_last;
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t dir = 1.0;
            if(r == 1) dir = -1.0;
            Float_t phi_val = dir*((Float_t)t/(Float_t)n_points)*delta_phi + phi_start;
            phi_val_last = phi_val;
            if(t == n_points) phi_val -= dir*0.0001;
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            //cout << ", x: " << x_tpc_val << ", y: " << y_tpc_val << endl;
            tp_Circles_fill->SetNextPoint(x_tpc_val,y_tpc_val);
            tp_Circles_lines[r] ->SetNextPoint(x_tpc_val,y_tpc_val);
            if(r == 0 && t == 0)
            {
                x_first = x_tpc_val;
                y_first = y_tpc_val;
            }
            if(r == 1 && t == n_points) tp_Circles_fill->SetNextPoint(x_first,y_first);
        }

        tp_Circles_lines[r]->SetLineStyle(line_style);
        tp_Circles_lines[r]->SetLineColorAlpha(color,transparency_line); // 28
        tp_Circles_lines[r]->SetLineWidth(line_width);
        tp_Circles_lines[r]->Draw("ogl");
    }

    tp_Circles_fill->SetLineStyle(line_style);
    tp_Circles_fill->SetLineColorAlpha(color,transparency_line); // 28
    tp_Circles_fill->SetLineWidth(0);
    tp_Circles_fill->SetFillColorAlpha(fill_color,transparency_fill);
    tp_Circles_fill->Draw("f");

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_Circle_3D(Float_t radius_in = 1, Float_t radius_out = 2,const Int_t n_radii = 1,
                             const Int_t n_delta_phi = 2, Float_t color = 2, Int_t line_style = 1, Int_t line_width = 1,
                             Float_t x_offset = 0.0, Float_t y_offset = 0.0, Float_t z_offset = 0.0, Int_t fill_color = 1, Double_t transparency = 0.0
                            )
{
    Float_t z = 0.0;
    const Int_t n_points = 50;
    TPolyLine3D   *tp_Circles[n_radii];
    TPolyLine3D   *tp_Radial[n_delta_phi];
    Float_t radius_table[n_radii];
    Float_t delta_radius;
    if(n_radii > 1) {delta_radius = (radius_out-radius_in)/((Float_t)(n_radii-1));}
    else{delta_radius = 0.0;}
    Float_t delta_phi    = 2.0*TMath::Pi()/((Float_t)n_delta_phi);
    Float_t z_tpc_val    = z;

    for(Int_t r = 0; r < n_radii; r++)
    {
        radius_table[r] = radius_in + r*delta_radius;
        tp_Circles[r] = new TPolyLine3D();
        Float_t radius   = radius_table[r];
        for(Int_t t = 0; t < n_points+1; t++)
        {
            Float_t phi_val = ((Float_t)t/(Float_t)n_points)*(2.0*TMath::Pi());
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            //cout << ", x: " << x_tpc_val << ", y: " << y_tpc_val << endl;
            tp_Circles[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_offset);
        }

        tp_Circles[r]->SetLineStyle(line_style);
        tp_Circles[r]->SetLineColor(color); // 28
        tp_Circles[r]->SetLineWidth(line_width);
        //tp_Circles[r]->DrawClone("oglf");
        tp_Circles[r]->Draw("f");
        tp_Circles[r]->Draw("ogl");
    }

    for(Int_t r = 0; r < n_delta_phi; r++)
    {
        tp_Radial[r] = new TPolyLine3D();
        Float_t phi_val = r*delta_phi;
        for(Int_t t = 0; t < 2; t++)
        {
            Float_t radius;
            if(t == 0) {radius = radius_table[0];}
            else {radius = radius_table[n_radii-1];}
            Float_t x_tpc_val   = radius*TMath::Cos(phi_val)+x_offset;
            Float_t y_tpc_val   = radius*TMath::Sin(phi_val)+y_offset;
            tp_Radial[r]->SetNextPoint(x_tpc_val,y_tpc_val,z_offset);
        }
        tp_Radial[r]->SetLineStyle(0);
        tp_Radial[r]->SetLineColor(color); // 28
        tp_Radial[r]->SetLineWidth(1);
        tp_Radial[r]->Draw("ogl");
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Decay_mom(Double_t M_mother, Double_t M_daughterA, Double_t M_daughterB)
{
    Double_t mom = 0.0;

    Double_t E_term = M_mother*M_mother - M_daughterB*M_daughterB + M_daughterA*M_daughterA;
    Double_t E_daughterA = 0.0;
    if(E_term > 0.0)
    {
        E_daughterA = E_term/(2.0*M_mother);
    }
    else return -1.0;

    mom = TMath::Sqrt(E_daughterA*E_daughterA - M_daughterA*M_daughterA);

    return mom;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGraphAsymmErrors* Calc_feed_down_v2(TGraphAsymmErrors* tgae_input_v2, Int_t PID, Int_t Energy, Double_t T_BW, Double_t rho0_BW,
                        Double_t rhoa_BW, Double_t s2_BW)
{
    // PID:
    // 0 = pi
    // 1 = charged K
    // 2 = p
    // 3 = anti-p
    // 4 = Lambda
    // 5 = anti-Lambda
    // 6 = K0S

    // Energy:
    // 0 = 7.7 GeV
    // 1 = 11.5 GeV
    // 2 = 19.6 GeV
    // 3 = 27 GeV
    // 4 = 39 GeV
    // 5 = 62.4 GeV
    // 6 = 200 GeV
    // 7 = 2.76 TeV

    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);

    const Int_t N_Energies           = 8;
    const Int_t N_Particles_v2_input = 7;  // 0     1   2   3      4            5      6     7     8    9      10        11          12            13     14    15   16   17       18            19          20    21   22    23      24
    const Int_t N_Particles          = 25; // pi+, phi, K-, K+, anti-Lambda, anti-Xi, K0S, Lambda, Xi-, rho, eta(958), omega(782), Delta(1232)++, Sigma0, Xi0, Omega, p, anti-p, anti-Sigma0, anti-Delta++, gamma, phi, N*  eta(547)  N*
    const Int_t N_feed_down_channels = 7;  // max value for all particles
    const Int_t N_decay_products     = 3;  // max number of decay daughters

    Int_t PID_in_table[N_Particles_v2_input] = {0,2,16,17,7,4,6};
    Int_t PID_in = PID_in_table[PID];

    const Int_t feed_down_PID_table[N_Particles_v2_input][N_feed_down_channels] =
    {
        //{9,12,22,23,11,10,24}, // pi
        {24,10,23,12,22,11,-1}, // pi   // removed 22 (maybe not so bad), 11 is the bad guy
        {1,-1,-1,-1,-1,-1,-1}, // K
        {7,12,22,-1,-1,-1,-1}, // p
        {4,19,-1,-1,-1,-1,-1}, // anti-p
        {8,13,14,15,-1,-1,-1}, // Lambda
        {5,18,-1,-1,-1,-1,-1}, // anti-Lambda
        {21,-1,-1,-1,-1,-1,-1}, // K0S
    };

    const Int_t decay_products[N_Particles][N_decay_products] =
    {
        {-1,-1,-1}, // pi
        {2,3,-1},   // phi -> K+ + K-
        {-1,-1,-1}, // K-
        {-1,-1,-1}, // K+
        {0,17,-1}, // anti-Lambda -> pi + anti-p
        {0,4,-1}, // anti-Xi -> pi + anti-Lambda
        {0,0,-1}, // K0S -> pi + pi
        {0,16,-1}, // Lambda -> pi + proton
        {0,7,-1}, // Xi- -> pi + Lambda
        {0,0,-1}, // rho -> pi + pi
        {0,0,23}, // eta(958) -> pi + pi + eta
        {0,0,0}, // omega(782) -> pi+ + pi- + pi0
        {0,16,-1}, // Delta(1232)++ -> pi + proton
        {20,7,-1}, // Sigma0 -> gamma + Lambda
        {0,7,-1}, // Xi0 -> pi + Lambda
        {2,7,-1}, // Omega -> K + Lambda
        {-1,-1,-1}, // p
        {-1,-1,-1}, // anti-p
        {20,4,-1}, // anti-Sigma0 -> gamma + anti-Lambda
        {0,17,-1}, // anti-Delta++ -> pi + anti-proton
        {0,0,-1}, // gamma
        {6,6,-1}, // phi -> K0S + K0L
        {0,16,-1}, // N* -> pi + p
        {0,0,0}, // eta(547) -> pi+ + pi- + pi0
        {0,0,16}, // N* -> pi + pi + p
    };
    //                                                   0    1    2   3    4     5     6     7     8     9    10    11    12    13    14     15  16  17   18     19  20   21    22     23     24
    const Double_t decay_momentum_table[N_Particles] = {0.0,0.127,0.0,0.0,0.101,0.139,0.206,0.101,0.139,0.363,0.232,0.327,0.229,0.074,0.135,0.211,0.0,0.0,0.074,0.229,0.0,0.110,0.405, 0.174, 0.380};
    //                                          pi+     phi      K-      K+   anti-Lambda anti-Xi  K0S     Lambda    Xi-    rho    eta(958) omega Delta  Sigma0    Xi0    Omega     p      anti-p   aSigma0 aDelta gamma phi   N*   eta(547) N*
    const Double_t mass_table[N_Particles] = {0.13957,1.01946,0.493677,0.493677,1.115683,1.32131,0.497648,1.115683,1.32131,0.7755,0.95778,0.78265,1.232,1.192642,1.31483,1.67245,0.938272,0.938272,1.192642,1.232,0.0,1.01946,1.480,0.54751,1.480};

    const Long64_t N_sample = 200000;

    // Particle yields from SHM predictions -> THERMUS (Macro_THERMUS_predict.cc)
    Double_t Res_scale        = 1.0;  // rescale for N*
    Double_t Res_scale_NStar  = 0.4;  // rescale for N* 0.4
    Double_t Res_scale_Lambda = 0.2;  // rescale for Lambda   -> estimated dca efficiency
    Double_t Res_scale_Omega  = 0.4;  // rescale for Omega    -> estimated dca efficiency
    Double_t Res_scale_Xi     = 0.6;  // rescale for Xi       -> estimated dca efficiency
    Double_t Res_rho          = 0.5;  // rescale for rho      -> in medium decays, less flow
    Double_t Res_eta          = 0.23; // rescale for eta(547) -> branching ratio
    Double_t Res_omega        = 0.3; // rescale for eta(547) -> branching ratio 0.89 but some decay in medium due to short live time
    Double_t Res_eta_958      = 0.45; // rescale for eta(958) -> branching ratio
    Double_t particle_yields[N_Energies][N_Particles] =
    {
        // 7.7
        {18.842, 0.0809822, 0.984182, 2.03403, Res_scale_Lambda*0.012904, Res_scale_Xi*0.00180397, 1.49986, Res_scale_Lambda*1.97031, Res_scale_Xi*0.110732, Res_rho*1.08777, Res_eta_958*0.0766547,Res_omega*0.883634, 1.58338, 0.395304, 0.110232, Res_scale_Omega*0.00556682, 11.2899, 0.0307144,0.00256662,0.00469493,0.0,0.0809822,Res_scale*1.58338,Res_eta*0.529146,Res_scale_NStar*1.58338},
        // 11.5
        {24.7732,0.166343,1.92581,3.05456,Res_scale_Lambda*0.0665792,Res_scale_Xi*0.00848224,2.45721,Res_scale_Lambda*2.20594,Res_scale_Xi*0.154628,Res_rho*1.89612,Res_eta_958*0.152232,Res_omega*1.59568,1.47253,0.416078,0.155051,Res_scale_Omega*0.0105926,10.252,0.173699,0.0124554,0.0263733,0.0,0.166343,Res_scale*1.47253,Res_eta*0.843459,Res_scale_NStar*1.58338},
        // 19.6
        {31.4526,0.260379,3.02932,3.90747,Res_scale_Lambda*0.221877,Res_scale_Xi*0.0257517,3.41128,Res_scale_Lambda*1.98909,Res_scale_Xi*0.164424,Res_rho*2.69009,Res_eta_958*0.233442,Res_omega*2.30929,1.14568,0.358709,0.165846,Res_scale_Omega*0.0140665,7.91359,0.635199,0.0397831,0.0947728,0.0,0.260379,Res_scale*1.14568,Res_eta*1.13144,Res_scale_NStar*1.58338},
        // 27
        {33.7056,0.296082,3.42981,4.23668,Res_scale_Lambda*0.348153,Res_scale_Xi*0.0395389,3.76868,Res_scale_Lambda*1.7408,Res_scale_Xi*0.149879,Res_rho*2.97559,Res_eta_958*0.26391,Res_omega*2.56822,0.963581,0.309601,0.151293,Res_scale_Omega*0.0135663,6.66259,1.019,0.0616097,0.151298,0.0,0.296082,Res_scale*0.963581,Res_eta*1.23154,Res_scale_NStar*1.58338},
        // 39
        {35.2158,0.3212,3.70381,4.47285,Res_scale_Lambda*0.497978,Res_scale_Xi*0.0558188,4.01892,Res_scale_Lambda*1.48852,Res_scale_Xi*0.131371,Res_rho*3.17241,Res_eta_958*0.285248,Res_omega*2.7473,0.802959,0.262355,0.132621,Res_scale_Omega*0.0123063,5.56283,1.47752,0.0873729,0.218829,0.0,0.3212,Res_scale*0.802959,Res_eta*1.29967,Res_scale_NStar*1.58338},
        // 62.4
        {36.2148,0.337546,3.91726,4.58277,Res_scale_Lambda*0.659449,Res_scale_Xi*0.0724835,4.17698,Res_scale_Lambda*1.27079,Res_scale_Xi*0.115194,Res_rho*3.29889,Res_eta_958*0.299095,Res_omega*2.86261,0.667969,0.222697,0.116366,Res_scale_Omega*0.0111575,4.62849,1.99556,0.115131,0.294752,0.0,0.337546,Res_scale*0.667969,Res_eta*1.34311,Res_scale_NStar*1.58338},
        // 200
        {36.9634,0.348295,4.15469,4.54569,Res_scale_Lambda*0.877946,Res_scale_Xi*0.0925872,4.27385,Res_scale_Lambda*1.03121,Res_scale_Xi*0.0978851,Res_rho*3.38142,Res_eta_958*0.308185,Res_omega*2.93795,0.520315,0.179972,0.0990989,Res_scale_Omega*0.00999385,3.59353,2.76675,0.152901,0.406583,0.0,0.348295,Res_scale*0.520315,Res_eta*1.37131,Res_scale_NStar*1.58338},
        // 2760
        {36.9634,0.348295,4.15469,4.54569,Res_scale_Lambda*0.877946,Res_scale_Xi*0.0925872,4.27385,Res_scale_Lambda*1.03121,Res_scale_Xi*0.0978851,Res_rho*3.38142,0.308185,Res_omega*2.93795,0.520315,0.179972,0.0990989,Res_scale_Omega*0.00999385,3.59353,2.76675,0.152901,0.406583,0.0,0.348295,Res_scale*0.520315,Res_eta*1.37131,Res_scale_NStar*1.58338} // copy of 200 GeV
    };

    // Fits to preliminary BES data
    // T = 0.153 for K+ at 11.5
    // T = 0.207 for p  at 11.5

    // T = 0.2 for K+ at 62.4
    // T = 0.2 for p  at 62.4

    Double_t T_kin[N_Energies] = {0.18,0.2,0.2,0.2,0.2,0.2,0.25,0.3};

    TF1 *FlowFit         = new TF1("FlowFit",FlowFitFunc,0,3.15,5);
    TF1 *PtFit2_mod_x    = new TF1("PtFit2_mod_x",PtFitFunc2_mod_x,0.0,2.5,4);
    TF1* BlastWaveFit    = new TF1("BlastWaveFit",BlastWaveFitFunc_no_mass_array,0.0,2.5,6);
    for(Int_t i = 0; i < 6; i++)
    {
        BlastWaveFit ->SetParameter(i,0.0);
        BlastWaveFit ->SetParError(i,0.0);
    }
    for(Int_t i = 0; i < 4; i++)
    {
        PtFit2_mod_x ->SetParameter(i,0.0);
        PtFit2_mod_x ->SetParError(i,0.0);
    }

    PtFit2_mod_x ->SetParameter(2,100.0);         // amplitude
    PtFit2_mod_x ->SetParameter(3,0.0);           // shift

    for(Int_t x = 0; x < 5; x++)
    {
        FlowFit->ReleaseParameter(x);
        FlowFit->SetParError(x,0.0);
        FlowFit->SetParameter(x,0.0);
    }

    FlowFit->SetParameter(0,1.0);
    FlowFit->SetParameter(1,0.0);
    FlowFit->SetParameter(2,0.0);
    FlowFit->SetParameter(3,0.0);
    FlowFit->SetParameter(4,0.0);

    BlastWaveFit ->SetParameter(0,T_BW);
    BlastWaveFit ->SetParameter(1,rho0_BW);
    BlastWaveFit ->SetParameter(2,rhoa_BW);
    BlastWaveFit ->SetParameter(3,s2_BW);
    BlastWaveFit ->SetParameter(4,1.0);

    TGraphAsymmErrors* tgae_output_v2;
    tgae_output_v2 = (TGraphAsymmErrors*)tgae_input_v2->Clone("tgae_output_v2");
    TProfile* tp_v2_feed_down = new TProfile("tp_v2_feed_down","tp_v2_feed_down",20,0,2.5);

    // Loop over the feed down channels
    TLorentzVector pMother, pDaughterA, pDaughterB, pDaughterAB, pDaughterC;
    Double_t total_weigh_factor = 0.0;
    Double_t input_weigh_factor = particle_yields[Energy][PID_in];
    for(Int_t iFeed_down = 0; iFeed_down < N_feed_down_channels; iFeed_down++)
    {
        Int_t PID_mother    = feed_down_PID_table[PID][iFeed_down];
        Int_t PID_daughterA = decay_products[PID_mother][0];
        Int_t PID_daughterB = decay_products[PID_mother][1];
        Int_t PID_daughterC = decay_products[PID_mother][2];
        if(PID_mother == -1) break;
        Double_t Mass_mother    = mass_table[PID_mother];
        Double_t Mass_daughterA = mass_table[PID_daughterA];
        Double_t Mass_daughterB = mass_table[PID_daughterB];
        Double_t Mass_daughterC = -1.0;
        Int_t flag_three_body = 0; // flag for three body decay
        if(PID_daughterC >= 0)
        {
            Mass_daughterC = mass_table[PID_daughterC];
            flag_three_body = 1;
        }
        Double_t decay_momentum = decay_momentum_table[PID_mother];
        BlastWaveFit ->SetParameter(5,Mass_mother);
        PtFit2_mod_x ->SetParameter(1,T_kin[Energy]+Mass_mother/20.0); // Effective temperature -> small mass dependence
        PtFit2_mod_x ->SetParameter(0,Mass_mother);
        Double_t weigh_factor = particle_yields[Energy][PID_mother]; // yield from SHM

        total_weigh_factor += weigh_factor;

        //if(PID_daughterA == PID_in && PID_daughterB == PID_in) // double the weigh since both particles are used, BUT only one is analyzed...
        //{
        //     total_weigh_factor += weigh_factor;
        //}

#if 1
        cout << "iFeed_down = " << iFeed_down << ", PID_mother = " << PID_mother << ", PID_daughterA = " << PID_daughterA  << ", PID_daughterB = " << PID_daughterB
            << ", Mass_mother = " << Mass_mother << ", Mass_daughterA = " << Mass_daughterA << ", Mass_daughterB = " << Mass_daughterB << ", decay_momentum = " << decay_momentum
            << ", input_weigh_factor = " << input_weigh_factor << ", weigh_factor = " << weigh_factor
            << endl;
#endif



        //TH1F* h_Mother_pt = new TH1F("h_Mother_pt","h_Mother_pt",100,0,1.5);
        //for(Int_t iSample = 0; iSample < 1000000; iSample++)
        //{
        //    Double_t pt_val_mother = PtFit2_mod_x ->GetRandom();
        //    h_Mother_pt ->Fill(pt_val_mother);
        //}

        // Sampling
        for(Int_t iSample = 0; iSample < N_sample; iSample++)
        {
            if (iSample != 0  &&  iSample % 200 == 0)
                cout << "." << flush;
            if (iSample != 0  &&  iSample % 2000 == 0)
            {
                Float_t event_percent = 100.0*iSample/N_sample;
                cout << " " << iSample << " (" << event_percent << "%) " << "\n" << "==> Processing data, " << flush;
            }


            Double_t pt_val_mother = PtFit2_mod_x ->GetRandom();
            //Double_t pt_val_mother =  h_Mother_pt->GetRandom();
            Double_t v2_val_mother = BlastWaveFit ->Eval(pt_val_mother);
            FlowFit->SetParameter(2,v2_val_mother);
            Double_t Psi_phi_val_mother = FlowFit ->GetRandom();
            //cout << "Psi_phi_val_mother = " << Psi_phi_val_mother << endl;

            pMother.SetXYZM(pt_val_mother,0.0,0.0,Mass_mother); // Mother particle

            if(flag_three_body == 0) // two body decay
            {
                pDaughterA.SetXYZM(decay_momentum,0.0,0.0,Mass_daughterA);
                pDaughterB.SetXYZM(-decay_momentum,0.0,0.0,Mass_daughterB);

                Double_t anglex = r3b.Rndm()*TMath::Pi()*2;
                Double_t angley = r3b.Rndm()*TMath::Pi()*2;
                Double_t anglez = r3b.Rndm()*TMath::Pi()*2;

                pDaughterA.RotateZ(anglez);
                pDaughterA.RotateX(anglex);
                pDaughterA.RotateY(angley);

                pDaughterB.RotateZ(anglez);
                pDaughterB.RotateX(anglex);
                pDaughterB.RotateY(angley);

                pDaughterA.Boost(pMother.BoostVector());
                pDaughterB.Boost(pMother.BoostVector());

                Double_t AngleA = (pMother.BoostVector()).DeltaPhi(pDaughterA.BoostVector());
                Double_t AngleB = (pMother.BoostVector()).DeltaPhi(pDaughterB.BoostVector());

                Double_t PtA = pDaughterA.Pt();
                Double_t PtB = pDaughterB.Pt();

                //cout << "AngleA = " << AngleA << ", AngleB = " << AngleB << ", PtA = " << PtA << ", PtB = " << PtB << endl;

                AngleA = Psi_phi_val_mother - AngleA;
                AngleB = Psi_phi_val_mother - AngleB;

                if(AngleA < 0.0) AngleA += TMath::Pi();
                if(AngleA > TMath::Pi()) AngleA -= TMath::Pi();

                if(AngleB < 0.0) AngleB += TMath::Pi();
                if(AngleB > TMath::Pi()) AngleB -= TMath::Pi();

                Double_t cosA = TMath::Cos(2.0*AngleA);
                Double_t cosB = TMath::Cos(2.0*AngleB);

                //tp_v2_feed_down ->Fill(pt_val_mother,v2_val_mother); // test, working
                //tp_v2_feed_down ->Fill(pt_val_mother,TMath::Cos(2.0*Psi_phi_val_mother)); // test, working


                if(PID_daughterA == PID_in)
                {
                    tp_v2_feed_down ->Fill(PtA,cosA,weigh_factor);
                }
                if(PID_daughterB == PID_in)
                {
                    tp_v2_feed_down ->Fill(PtB,cosB,weigh_factor);
                }
            }

            if(flag_three_body == 1) // three body decay
            {
                // first decay
                Double_t low_limit_Mass_first = Mass_daughterA + Mass_daughterB;
                Double_t up_limit_Mass_first  = Mass_mother - Mass_daughterC;
                Double_t Mass_daughter_first  = r3b.Rndm()*(up_limit_Mass_first - low_limit_Mass_first) + low_limit_Mass_first;

                if(up_limit_Mass_first > low_limit_Mass_first)
                {
                    Double_t Decay_mom_first   = Decay_mom(Mass_mother,Mass_daughterC,Mass_daughter_first);

                    pDaughterAB.SetPxPyPzE(Decay_mom_first,0.0,0.0,TMath::Sqrt(Mass_daughterC*Mass_daughterC+Decay_mom_first*Decay_mom_first));
                    pDaughterC.SetPxPyPzE(-Decay_mom_first,0.0,0.0,TMath::Sqrt(Mass_daughter_first*Mass_daughter_first+Decay_mom_first*Decay_mom_first));

                    Double_t anglex = r3b.Rndm()*TMath::Pi()*2;
                    Double_t angley = r3b.Rndm()*TMath::Pi()*2;
                    Double_t anglez = r3b.Rndm()*TMath::Pi()*2;

                    pDaughterAB.RotateZ(anglez);
                    pDaughterAB.RotateX(anglex);
                    pDaughterAB.RotateY(angley);

                    pDaughterC.RotateZ(anglez);
                    pDaughterC.RotateX(anglex);
                    pDaughterC.RotateY(angley);

                    pDaughterAB.Boost(pMother.BoostVector());
                    pDaughterC.Boost(pMother.BoostVector());


                    // second decay
                    Double_t Decay_mom_second   = Decay_mom(Mass_daughter_first,Mass_daughterA,Mass_daughterB);


                    pDaughterA.SetXYZM(Decay_mom_second,0.0,0.0,Mass_daughterA);
                    pDaughterB.SetXYZM(-Decay_mom_second,0.0,0.0,Mass_daughterB);

                    anglex = r3b.Rndm()*TMath::Pi()*2;
                    angley = r3b.Rndm()*TMath::Pi()*2;
                    anglez = r3b.Rndm()*TMath::Pi()*2;

                    pDaughterA.RotateZ(anglez);
                    pDaughterA.RotateX(anglex);
                    pDaughterA.RotateY(angley);

                    pDaughterB.RotateZ(anglez);
                    pDaughterB.RotateX(anglex);
                    pDaughterB.RotateY(angley);

                    pDaughterA.Boost(pDaughterAB.BoostVector());
                    pDaughterB.Boost(pDaughterAB.BoostVector());

                    Double_t AngleA = (pMother.BoostVector()).DeltaPhi(pDaughterA.BoostVector());
                    Double_t AngleB = (pMother.BoostVector()).DeltaPhi(pDaughterB.BoostVector());
                    Double_t AngleC = (pMother.BoostVector()).DeltaPhi(pDaughterC.BoostVector());

                    Double_t PtA = pDaughterA.Pt();
                    Double_t PtB = pDaughterB.Pt();
                    Double_t PtC = pDaughterC.Pt();

                    //cout << "AngleA = " << AngleA << ", AngleB = " << AngleB << ", PtA = " << PtA << ", PtB = " << PtB << endl;

                    AngleA = Psi_phi_val_mother - AngleA;
                    AngleB = Psi_phi_val_mother - AngleB;
                    AngleC = Psi_phi_val_mother - AngleC;

                    if(AngleA < 0.0) AngleA += TMath::Pi();
                    if(AngleA > TMath::Pi()) AngleA -= TMath::Pi();

                    if(AngleB < 0.0) AngleB += TMath::Pi();
                    if(AngleB > TMath::Pi()) AngleB -= TMath::Pi();

                    if(AngleC < 0.0) AngleC += TMath::Pi();
                    if(AngleC > TMath::Pi()) AngleC -= TMath::Pi();

                    Double_t cosA = TMath::Cos(2.0*AngleA);
                    Double_t cosB = TMath::Cos(2.0*AngleB);
                    Double_t cosC = TMath::Cos(2.0*AngleC);

                    //tp_v2_feed_down ->Fill(pt_val_mother,v2_val_mother); // test, working
                    //tp_v2_feed_down ->Fill(pt_val_mother,TMath::Cos(2.0*Psi_phi_val_mother)); // test, working


                    if(PID_daughterA == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtA,cosA,weigh_factor);
                    }
                    if(PID_daughterB == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtB,cosB,weigh_factor);
                    }
                    if(PID_daughterC == PID_in)
                    {
                        tp_v2_feed_down ->Fill(PtC,cosC,weigh_factor);
                    }

                }
            }

        }
    }

    //Double_t diff_weigh_factor = input_weigh_factor - total_weigh_factor;
    Double_t diff_weigh_factor = input_weigh_factor;
    if(diff_weigh_factor <= 0.0) diff_weigh_factor = 1.0;

    if(PID == 0) // special treatment for pions, resonances have a 60% contribution, see http://arxiv.org/pdf/hep-ph/0407174.pdf, Fig. 1
    {
        total_weigh_factor = 0.6;
        input_weigh_factor = 0.4;
        diff_weigh_factor  = input_weigh_factor;
    }

    for(Int_t iPoint = 0; iPoint < tgae_output_v2->GetN(); iPoint++)
    {
        Double_t x_val, y_val;
        tgae_output_v2 ->GetPoint(iPoint,x_val,y_val);
        Double_t exl_ME = tgae_input_v2 ->GetErrorXlow(iPoint);
        Double_t exh_ME = tgae_input_v2 ->GetErrorXhigh(iPoint);
        Double_t eyl_ME = tgae_input_v2 ->GetErrorYlow(iPoint);
        Double_t eyh_ME = tgae_input_v2 ->GetErrorYhigh(iPoint);

        Double_t y_val_FD     = tp_v2_feed_down->GetBinContent(tp_v2_feed_down->FindBin(x_val));
        Double_t y_val_err_FD = tp_v2_feed_down->GetBinError(tp_v2_feed_down->FindBin(x_val));

        //Double_t y_val_out = (y_val*input_weigh_factor - y_val_FD*total_weigh_factor)/diff_weigh_factor;
        Double_t y_val_out = (y_val*(input_weigh_factor+total_weigh_factor) - y_val_FD*total_weigh_factor)/diff_weigh_factor;

        tgae_output_v2 ->SetPoint(iPoint,x_val,y_val_out); // original
        //tgae_output_v2 ->SetPoint(iPoint,x_val,y_val_FD);  // test
        //tgae_output_v2 ->SetPointError(iPoint,0.0,0.0,y_val_err_FD,y_val_err_FD); // test

        cout << "iPoint = " << iPoint << ", y_val_FD = " << y_val_FD << ", y_val_in = " << y_val << ", y_val_out = " << y_val_out
            << ", input_weigh_factor = " << input_weigh_factor << ", total_weigh_factor = " << total_weigh_factor << endl;
    }

    return tgae_output_v2;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// NEW
// Calculates the maximum/minimum envelope curve as systematic error for any number of added TGraphAsymmErrors/TH1F
// y-poins are the mean of the maximum/minimum
// Graphs need to have the same points in x
class Calc_syst_error_graph_hist
{
    TGraphAsymmErrors* tgae_stat;
    TGraphAsymmErrors* tgae_syst;
    vector<TGraphAsymmErrors*> vec_tgae;
    vector<TH1F*> vec_hist;
    Int_t flag_input, flag_add_stat_error;
public:
    void add_graph(TGraphAsymmErrors*);
    void add_hist(TH1F*);
    void add_stat_error(Int_t);
    void calculate_error();
    TGraphAsymmErrors* get_graph_stat();
    TGraphAsymmErrors* get_graph_syst();
    void clear();
};


void Calc_syst_error_graph_hist::add_graph(TGraphAsymmErrors* graph)
{
    vec_tgae.push_back(graph);
    flag_input = 0;
    flag_add_stat_error = 0;
}

void Calc_syst_error_graph_hist::add_hist(TH1F* hist)
{
    vec_hist.push_back(hist);
    flag_input = 1;
    flag_add_stat_error = 0;
}

void Calc_syst_error_graph_hist::add_stat_error(Int_t add_stat)
{
    if(add_stat == 0)
    {
        flag_add_stat_error = 0; // add statistical error to systematic value
    }
    if(add_stat == 1)
    {
        flag_add_stat_error = 1; // use maximum relative statistical error
    }
    if(add_stat >= 2)
    {
        flag_add_stat_error = 2; // use average of relative statistical error
    }
}

void Calc_syst_error_graph_hist::calculate_error()
{
    Int_t start_point, stop_point, size;
    if(flag_input == 0) // TGraphAsymmErrors
    {
        tgae_stat   = (TGraphAsymmErrors*)vec_tgae[0]->Clone("tgae_stat");
        tgae_syst   = (TGraphAsymmErrors*)vec_tgae[0]->Clone("tgae_syst");
        start_point = 0;
        stop_point  = vec_tgae[0]->GetN();
        size        = vec_tgae.size();
    }
    if(flag_input == 1) // TH1F
    {
        tgae_stat   = new TGraphAsymmErrors(vec_hist[0]);
        tgae_stat   ->SetName("tgae_stat");
        tgae_syst   = new TGraphAsymmErrors(vec_hist[0]);
        tgae_stat   ->SetName("tgae_syst");
        start_point = 1;
        stop_point  = vec_hist[0]->GetNbinsX()+1;
        size        = vec_hist.size();
    }

    //cout << "start_point: " << start_point << ", stop_point: " << stop_point << ", flag_input: " << flag_input << ", size: " << size << endl;

    for(Int_t i_point = start_point; i_point < stop_point; i_point++)
    {
        Double_t max_y_val, min_y_val, max_y_err, min_y_err, max_y_rel_err, min_y_rel_err, x_val, x_err_low, x_err_high;
        Int_t i_file_use_counter = 0;
        Double_t average_rel_high_error = 0.0;
        Double_t average_rel_low_error  = 0.0;
        for(Int_t i_file = 0; i_file < size; i_file++)
        {
            Double_t x,y;
            Double_t Xhigh,Xlow,Yhigh,Ylow;

            if(flag_input == 0) // TGraphAsymmErrors
            {
                vec_tgae[i_file]->GetPoint(i_point,x,y);
                Xhigh = vec_tgae[i_file]->GetErrorXhigh(i_point);
                Xlow  = vec_tgae[i_file]->GetErrorXlow(i_point);
                Yhigh = vec_tgae[i_file]->GetErrorYhigh(i_point);
                Ylow  = vec_tgae[i_file]->GetErrorYlow(i_point);
            }
            if(flag_input == 1) // TH1F
            {
                x     = vec_hist[i_file]->GetBinCenter(i_point);
                y     = vec_hist[i_file]->GetBinContent(i_point);
                Xlow  = vec_hist[i_file]->GetBinWidth(i_point)/2.0;
                Xhigh = Xlow;
                Yhigh = vec_hist[i_file]->GetBinError(i_point);
                Ylow  = Yhigh;
                //cout << "i_point: " << i_point << ", x: " << x << ", y: " << y << ", Xlow: " << Xlow << ", Ylow: " << Ylow << endl;
            }

            //if(y - Ylow < 1E-09 && x > 5.0) continue; // avoid huge negative error fluctuations
            i_file_use_counter++;

            if(i_file_use_counter == 1)
            {
                if(flag_add_stat_error == 0)
                {
                    max_y_val  = y + Yhigh;
                    min_y_val  = y - Ylow;
                }
                if(flag_add_stat_error >= 1)
                {
                    max_y_val  = y;
                    min_y_val  = y;
                }
                x_val          = x;
                x_err_low      = Xlow;
                x_err_high     = Xhigh;
                max_y_err      = Yhigh;
                min_y_err      = Ylow;
                max_y_rel_err  = Yhigh/y;
                min_y_rel_err  = Ylow/y;
            }
            else
            {
                if(flag_add_stat_error == 0)
                {
                    if( (y + Yhigh) > max_y_val)
                    {
                        max_y_val = y + Yhigh;
                        max_y_err = Yhigh;
                    }
                    if( (y - Ylow) < min_y_val)
                    {
                        min_y_val = y - Ylow;
                        min_y_err = Ylow;
                    }
                }
                if(flag_add_stat_error >= 1)
                {
                    if( y > max_y_val)
                    {
                        max_y_val = y;
                    }
                    if( y < min_y_val)
                    {
                        min_y_val = y;
                    }

                    if(y - Ylow < 1E-09) continue; // avoid large negative fluctuations in unfolding

                    // Absolute error
                    if(Yhigh > max_y_err)
                    {
                        max_y_err = Yhigh;
                    }
                    if(Ylow > min_y_err)
                    {
                        min_y_err = Ylow;
                    }

                    // Relative error
                    if((Yhigh/y) > max_y_rel_err)
                    {
                        max_y_rel_err = Yhigh/y;
                    }
                    if((Ylow/y) > min_y_rel_err)
                    {
                        min_y_rel_err = Ylow/y;
                    }

                    // Sum of relative errors
                    average_rel_high_error += Yhigh/y;
                    average_rel_low_error  += Ylow/y;
                }
            }

            //if(x > 10.0 && x < 11.2) cout << "i_file: " << i_file << ", Ylow: " << Ylow << ", Yhigh: " << Yhigh << endl;
        }
        Double_t average_y_val = (max_y_val+min_y_val)/2.0;
        if(i_file_use_counter > 0)
        {
            average_rel_low_error  /= (Double_t)i_file_use_counter;
            average_rel_high_error /= (Double_t)i_file_use_counter;
        }
        tgae_syst->SetPoint(i_point,x_val,average_y_val);
        tgae_stat->SetPoint(i_point,x_val,average_y_val);
        //cout << "i_point: " << i_point << ", x_val: " << x_val << ", average_y_val: " << average_y_val << ", min_y_val: " << min_y_val << ", max_y_val: " << max_y_val
        //    << ", min_y_err: " << min_y_err << ", max_y_err: " << max_y_err << endl;
        tgae_syst->SetPointError(i_point,x_err_low,x_err_high,average_y_val-min_y_val,max_y_val-average_y_val);
        if(flag_add_stat_error == 0)
        {
            tgae_stat->SetPointError(i_point,x_err_low,x_err_high,0.0,0.0);
        }
        if(flag_add_stat_error == 1)
        {
            //tgae_stat->SetPointError(i_point,x_err_low,x_err_high,min_y_err,max_y_err);
            tgae_stat->SetPointError(i_point,x_err_low,x_err_high,min_y_rel_err*average_y_val,max_y_rel_err*average_y_val);
            //cout << "i_point: " << i_point << ", x_val: " << x_val << ", average_y_val: " << average_y_val << ", min_y_err: " << min_y_err << endl;
        }
        if(flag_add_stat_error == 2)
        {
            //tgae_stat->SetPointError(i_point,x_err_low,x_err_high,min_y_err,max_y_err);
            tgae_stat->SetPointError(i_point,x_err_low,x_err_high,average_y_val*average_rel_low_error,average_y_val*average_rel_high_error);
            //cout << "i_point: " << i_point << ", x_val: " << x_val << ", average_y_val: " << average_y_val << ", min_y_err: " << min_y_err << endl;
        }
    }
}

TGraphAsymmErrors* Calc_syst_error_graph_hist::get_graph_stat()
{
    return tgae_stat;
}

TGraphAsymmErrors* Calc_syst_error_graph_hist::get_graph_syst()
{
    return tgae_syst;
}

void Calc_syst_error_graph_hist::clear()
{
    if(flag_input == 0) // TGraphAsymmErrors
    {
        vec_tgae.clear();
    }
    if(flag_input == 1) // TH1F
    {
        vec_hist.clear();
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calculates the errors of a ratio of numbers a and b based on Monte-Carlo
class Monte_Carlo_Error_Ratio
{
    Double_t nom, den, nom_err, den_err;
    Long64_t N_MC;
    TRandom  ran;
    Double_t RMS[2];
    Double_t ratio_nom_den;
public:
    void set_nominator_denominator(Double_t, Double_t);
    void set_nominator_denominator_error(Double_t, Double_t);
    void set_N_MC(Long64_t);
    void calc_errors();
    void truncate_RMS(Double_t);
    Double_t get_RMS_above();
    Double_t get_RMS_below();
    void clear();
    std::vector< std::vector<Double_t> > vec_RMS;
    std::vector< std::vector<Double_t> > vec_RMS_trunc;
};

void Monte_Carlo_Error_Ratio::set_nominator_denominator(Double_t a, Double_t b)
{
    nom = a;
    den = b;
}

void Monte_Carlo_Error_Ratio::set_nominator_denominator_error(Double_t a, Double_t b)
{
    nom_err = a;
    den_err = b;
}

void Monte_Carlo_Error_Ratio::set_N_MC(Long64_t i)
{
    N_MC = i;
}

void Monte_Carlo_Error_Ratio::calc_errors()
{
    vec_RMS.resize(2); // two RMS values are calculated, above and below mean value
    vec_RMS_trunc.resize(2);
    ran.SetSeed(0); // seed for random number generator changes every second
    Long64_t N_ratios = 0;
    RMS[0] = 0.0; // two RMS values are calculated, above and below mean value
    RMS[1] = 0.0;
    if(den > 0.0)
    {
        ratio_nom_den = nom/den;
    }
    else
    {
        RMS[0] = -1;
        RMS[1] = -1;
        cout << "nom: " << nom << ", den: " << den << endl;
    }
    if(RMS[0] >= 0.0)
    {
        for(Int_t i_MC = 0; i_MC < N_MC; i_MC++)
        {
            Double_t val_nom = ran.Gaus(nom,nom_err);
            Double_t val_den = ran.Gaus(den,den_err);

            if(val_den > 0.0 && val_nom > 0.0) // take only positive values (yields)
            {
                Double_t ratio_val = val_nom/val_den;
                Double_t RMS_val   = ratio_val - ratio_nom_den; // MC ratio - average
                Int_t RMS_above_below = 0; // above
                if(RMS_val < 0.0) RMS_above_below = 1; // below
                vec_RMS[RMS_above_below].push_back(RMS_val);
                RMS[RMS_above_below] += RMS_val*RMS_val;
#if 0
                cout << "i_MC: " << i_MC << ", RMS_above_below: " << RMS_above_below << ", val_nom: " << val_nom << ", nom_err: " << nom_err << ", val_den: " << val_den << ", den_err: " << den_err
                    << ", ratio_val: " << ratio_val << ", ratio_nom_den (mean): " << ratio_nom_den << ", N_ratios: " << N_ratios
                    << ", RMS_val: " << RMS[RMS_above_below]  << endl;
#endif
            }
        }
        for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
        {
            if(vec_RMS[i_above_below].size() > 0)
            {
                RMS[i_above_below] = TMath::Sqrt((1.0/((Double_t)vec_RMS[i_above_below].size()))*RMS[i_above_below]);
            }
            else
            {
                RMS[i_above_below] = -1;
                cout << "WARNING, RMS vector size below 1, old RMS used" << endl;
            }
        }
    }
    //cout << "RMS[0]: " << RMS[0] << endl;
}

void Monte_Carlo_Error_Ratio::truncate_RMS(Double_t nSigma)
{
    for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
    {
        //cout << "size RMS: " << vec_RMS[i_above_below].size() << endl;
        for(Int_t i_entry = 0; i_entry < vec_RMS[i_above_below].size(); i_entry++)
        {
            //cout << "i_entry: " << i_entry << ", nSigma: " << nSigma << ", RMS average: " << RMS[i_above_below]
            //    << ", diff: " << fabs(vec_RMS[i_above_below][i_entry] - ratio_nom_den) << ", RMS: "
            //    << vec_RMS[i_above_below][i_entry] << ", ratio_nom_den: " << ratio_nom_den << endl;
            //if( fabs(vec_RMS[i_above_below][i_entry] - ratio_nom_den) < nSigma*RMS[i_above_below])
            if( fabs(vec_RMS[i_above_below][i_entry]) < nSigma*RMS[i_above_below])
            {
                //cout << "ACCEPTED" << endl;
                vec_RMS_trunc[i_above_below].push_back(vec_RMS[i_above_below][i_entry]);
            }
            else
            {
                //cout << "i_entry: " << i_entry << ", REJECTED" << endl;
            }
        }
        //cout << "size RMS trunc: " << vec_RMS_trunc[i_above_below].size() << endl;
        if(vec_RMS_trunc[i_above_below].size() > 0) // make sure you always have entries
        {
            RMS[i_above_below] = 0.0;
            vec_RMS[i_above_below].clear();
            for(Int_t i_entry = 0; i_entry < vec_RMS_trunc[i_above_below].size(); i_entry++)
            {
                RMS[i_above_below] += vec_RMS_trunc[i_above_below][i_entry]*vec_RMS_trunc[i_above_below][i_entry];
                vec_RMS[i_above_below].push_back(vec_RMS_trunc[i_above_below][i_entry]);
            }
            RMS[i_above_below] = TMath::Sqrt((1.0/((Double_t)vec_RMS_trunc[i_above_below].size()))*RMS[i_above_below]);
        }
        vec_RMS_trunc[i_above_below].clear();
    }
}

Double_t Monte_Carlo_Error_Ratio::get_RMS_above()
{
    return RMS[0];
}
Double_t Monte_Carlo_Error_Ratio::get_RMS_below()
{
    return RMS[1];
}
void Monte_Carlo_Error_Ratio::clear()
{
    for(Int_t i_above_below = 0; i_above_below < 2; i_above_below++)
    {
        vec_RMS[i_above_below].clear();
        vec_RMS_trunc[i_above_below].clear();
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calculates the hit point of a ray with a TGeoVolume object
class C_Hit_on_Volume
{
    TGeoVolume* geo_volume;
    TVector3 TV_observer_point, TV_focal_point, TV_ray_vec, TV_hit;
    Int_t N_itt;
    Double_t f_scale;
    Double_t shift_z;

public:
    void set_TGeoVolume(TGeoVolume*); // object
    void set_start_point(Double_t, Double_t, Double_t); // observer point
    void set_stop_point(Double_t, Double_t, Double_t); // focal point
    void calc_hit_point();
    void set_N_itterations(Int_t);
    void set_scale_factor(Double_t);
    void set_shift_z(Double_t);
    TVector3 get_hit_point();
};

void C_Hit_on_Volume::set_TGeoVolume(TGeoVolume* i_geo_volume)
{
    // Default values
    N_itt   = 1000;
    f_scale = 5.0;
    shift_z = 0.0;
    geo_volume = i_geo_volume;
}

void C_Hit_on_Volume::set_start_point(Double_t x, Double_t y, Double_t z)
{
    TV_observer_point.SetXYZ(x,y,z);
}

void C_Hit_on_Volume::set_stop_point(Double_t x, Double_t y, Double_t z)
{
    TV_focal_point.SetXYZ(x,y,z);
}

void C_Hit_on_Volume::set_N_itterations(Int_t i_itt)
{
    N_itt = i_itt;
}

void C_Hit_on_Volume::set_scale_factor(Double_t i_scale)
{
    f_scale = i_scale;
}

void C_Hit_on_Volume::set_shift_z(Double_t i_shift_z)
{
    shift_z = i_shift_z;
}

void C_Hit_on_Volume::calc_hit_point()
{
    // Calculate normalized ray vector
    TV_ray_vec = TV_focal_point;
    TV_ray_vec -= TV_observer_point;
    TV_ray_vec *= 1.0/TV_ray_vec.Mag();

    Bool_t is_inside = kFALSE;
    for(Int_t i_sign = 0; i_sign < 2; i_sign++)
    {
        Double_t f_sign = 1.0;
        if(i_sign == 1) f_sign = -1.0;
        for(Int_t i_itt = 0; i_itt < N_itt; i_itt++)
        {
            TV_hit = TV_observer_point;
            TVector3 TV_ray_vec_scale = TV_ray_vec;
            TV_ray_vec_scale *= f_sign*f_scale*i_itt;
            TV_hit += TV_ray_vec_scale;

            Double_t vec_hit_point[3] = {TV_hit[0], TV_hit[1], TV_hit[2] + shift_z};
            is_inside = geo_volume->Contains(vec_hit_point);

            if(is_inside)
            {
                break;
            }
        }
        if(is_inside) break; // Don't proceed with negative direction if hit point was already found
    }
}

TVector3 C_Hit_on_Volume::get_hit_point()
{
    return TV_hit;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGraphAsymmErrors* ReflectTGAE(TGraphAsymmErrors* tgae_in, Double_t reflect_x_val)
{
    // Reflect the graph
    TGraphAsymmErrors* tgae_out = (TGraphAsymmErrors*)tgae_in->Clone("tgae_out");

    Double_t xA,xB,yA,yB;
    tgae_out->GetPoint(1,xA,yA);
    tgae_out->GetPoint(2,xB,yB);
    Double_t bin_width = fabs(xB-xA);
    //cout << "xA: " << xA << ", xB: " << xB << ", bin_width: " << bin_width << endl;
    for(Int_t i_point = 0; i_point < tgae_in->GetN(); i_point++)
    {
        Double_t x,y;
        tgae_in->GetPoint(i_point,x,y);
        Double_t Xhigh = tgae_in->GetErrorXhigh(i_point);
        Double_t Xlow  = tgae_in->GetErrorXlow(i_point);
        Double_t Yhigh = tgae_in->GetErrorYhigh(i_point);
        Double_t Ylow  = tgae_in->GetErrorYlow(i_point);

        Double_t x_new;
        if(x <= reflect_x_val) x_new = (reflect_x_val - x) + reflect_x_val;
        if(x > reflect_x_val)  x_new = reflect_x_val - (x - reflect_x_val);

        if(x <= reflect_x_val)
        {
            for(Int_t i_pointR = tgae_in->GetN()-1; i_pointR >= 0; i_pointR--)
            {
                Double_t xR,yR;
                tgae_out->GetPoint(i_pointR,xR,yR);
                Double_t XhighR = tgae_in->GetErrorXhigh(i_pointR);
                Double_t XlowR  = tgae_in->GetErrorXlow(i_pointR);
                Double_t YhighR = tgae_in->GetErrorYhigh(i_pointR);
                Double_t YlowR  = tgae_in->GetErrorYlow(i_pointR);
                if(xR > reflect_x_val)
                {
                    if(
                       ((xR-reflect_x_val - bin_width/2.0) <= (reflect_x_val-x)) &&
                       ((xR-reflect_x_val + bin_width/2.0) > (reflect_x_val-x))
                      )
                    {
                        //cout << "success, x: " << x << ", xR: " << xR << endl;
                        //cout << "val1: " << (xR-reflect_x_val - bin_width/2.0) << ", val2: " << (xR-reflect_x_val + bin_width/2.0) << ", valR: " << (reflect_x_val-x) << endl;
                        Double_t yR_new = (y+yR)/2.0;
                        Double_t Ylow_new  = TMath::Sqrt(TMath::Power(YlowR,2) + TMath::Power(Ylow,2))/2.0;
                        Double_t Yhigh_new = TMath::Sqrt(TMath::Power(YhighR,2) + TMath::Power(Yhigh,2))/2.0;
                        tgae_out->SetPoint(i_pointR,xR,yR_new);
                        tgae_out->SetPointError(i_pointR,XlowR,XhighR,Ylow_new,Yhigh_new);
                    }
                }
                else break;
            }
        }
        else break;
    }
    return tgae_out;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TGraphAsymmErrors* ScaleTGAE(TGraphAsymmErrors* tgae_in, Double_t scale_fac)
{
    // Scale the graph
    TGraphAsymmErrors* tgae_out = (TGraphAsymmErrors*)tgae_in->Clone("tgae_out");

    Double_t xA,xB,yA,yB;
    tgae_out->GetPoint(1,xA,yA);
    tgae_out->GetPoint(2,xB,yB);
    Double_t bin_width = fabs(xB-xA);
    //cout << "xA: " << xA << ", xB: " << xB << ", bin_width: " << bin_width << endl;
    for(Int_t i_point = 0; i_point < tgae_in->GetN(); i_point++)
    {
        Double_t x,y;
        tgae_in->GetPoint(i_point,x,y);
        Double_t Xhigh = tgae_in->GetErrorXhigh(i_point);
        Double_t Xlow  = tgae_in->GetErrorXlow(i_point);
        Double_t Yhigh = tgae_in->GetErrorYhigh(i_point);
        Double_t Ylow  = tgae_in->GetErrorYlow(i_point);

        tgae_out->SetPoint(i_point,x,y*scale_fac);
        tgae_out->SetPointError(i_point,Xlow,Xhigh,Ylow*scale_fac,Yhigh*scale_fac);
    }
    return tgae_out;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calculates average tgae and systematic error due to ME normalization
class Get_tgae_average_norm_syst
{
    vector<TGraphAsymmErrors*> vec_tgae; // input graphs
    TGraphAsymmErrors* tgae_central; // central graph + statistical errors
    TGraphAsymmErrors* tgae_syst; // systematic error due to normalization
    TGraphAsymmErrors* tgae_total_error; // systematic error due to normalization
public:
    void add_graph(TGraphAsymmErrors*);
    void calculate_average();
    TGraphAsymmErrors* get_syst_graph(); // systematic error due to normalization
    TGraphAsymmErrors* get_central_graph(); // includes statistical error
    TGraphAsymmErrors* get_total_err_graph(); // includes statistical and systematical errors
    void clear();
};

void Get_tgae_average_norm_syst::clear()
{
        vec_tgae.clear();
}

void Get_tgae_average_norm_syst::add_graph(TGraphAsymmErrors* graph)
{
    vec_tgae.push_back(graph);
}

TGraphAsymmErrors* Get_tgae_average_norm_syst::get_syst_graph()
{
    return tgae_syst;
}

TGraphAsymmErrors* Get_tgae_average_norm_syst::get_central_graph()
{
    return tgae_central;
}

TGraphAsymmErrors* Get_tgae_average_norm_syst::get_total_err_graph()
{
    return tgae_total_error;
}

void Get_tgae_average_norm_syst::calculate_average()
{
    Int_t size       = vec_tgae.size();
    tgae_central     = new TGraphAsymmErrors();
    tgae_syst        = new TGraphAsymmErrors();
    tgae_total_error = new TGraphAsymmErrors();

    for(Int_t i_point = 0; i_point < vec_tgae[0]->GetN(); i_point++)
    {
        Double_t x[size],y[size], XEhigh[size], XElow[size], YEhigh[size], YElow[size];
        Double_t y_min, y_max;
        for(Int_t i_range = 0; i_range < size; i_range++)
        {
            vec_tgae[i_range]->GetPoint(i_point,x[i_range],y[i_range]);
            XEhigh[i_range] = vec_tgae[i_range]->GetErrorXhigh(i_point);
            XElow[i_range]  = vec_tgae[i_range]->GetErrorXlow(i_point);
            YEhigh[i_range] = vec_tgae[i_range]->GetErrorYhigh(i_point);
            YElow[i_range]  = vec_tgae[i_range]->GetErrorYlow(i_point);
            if(i_range == 0)
            {
                y_min = y[i_range];
                y_max = y[i_range];
            }
            else
            {
                if(y[i_range] < y_min) y_min = y[i_range];
                if(y[i_range] > y_max) y_max = y[i_range];
            }
        }
        Double_t y_average = (y_max+y_min)/2.0;
        Double_t YEhigh_average = (YEhigh[0] + YEhigh[1] + YEhigh[2])/3.0;
        Double_t YElow_average  = (YElow[0]  + YElow[1]  + YElow[2])/3.0;

        tgae_central     ->SetPoint(i_point,x[0],y_average);
        tgae_central     ->SetPointError(i_point,XElow[0],XEhigh[0],YElow_average,YEhigh_average);

        tgae_syst        ->SetPoint(i_point,x[0],y_average);
        tgae_syst        ->SetPointError(i_point,XElow[0],XEhigh[0],fabs(y_average-y_min),fabs(y_average-y_max));

        Double_t YElow_total  = TMath::Sqrt(TMath::Power(YElow_average,2.0)  + TMath::Power(fabs(y_average-y_min),2.0));
        Double_t YEhigh_total = TMath::Sqrt(TMath::Power(YEhigh_average,2.0) + TMath::Power(fabs(y_average-y_max),2.0));
        tgae_total_error ->SetPoint(i_point,x[0],y_average);
        tgae_total_error ->SetPointError(i_point,XElow[0],XEhigh[0],YElow_total,YEhigh_total);
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Calculates the errors of a ratio of numbers a and b based on Monte-Carlo
class Calculate_TGAE_average_and_envelope
{
    vector<TGraphAsymmErrors*> vec_tgae; // input graphs
    TGraphAsymmErrors* tgae_central; // central graph + envelope as errors
    TGraphAsymmErrors* tgae_envelope; // includes statistical error
public:
    void add_graph(TGraphAsymmErrors*);
    void calculate_average();
    TGraphAsymmErrors* get_envelope_graph(); // central graph + envelope as errors
    TGraphAsymmErrors* get_central_graph(); // includes statistical error
    void clear();
};

void Calculate_TGAE_average_and_envelope::clear()
{
        vec_tgae.clear();
}

void Calculate_TGAE_average_and_envelope::add_graph(TGraphAsymmErrors* graph)
{
    vec_tgae.push_back(graph);
}

TGraphAsymmErrors* Calculate_TGAE_average_and_envelope::get_envelope_graph()
{
    return tgae_envelope;
}

TGraphAsymmErrors* Calculate_TGAE_average_and_envelope::get_central_graph()
{
    return tgae_central;
}

void Calculate_TGAE_average_and_envelope::calculate_average()
{
    Int_t size = vec_tgae.size();
    tgae_central  = new TGraphAsymmErrors();
    tgae_envelope = new TGraphAsymmErrors();

    for(Int_t i_point = 0; i_point < vec_tgae[0]->GetN(); i_point++)
    {
        Double_t average_y_val      = 0.0;
        Double_t weight_average     = 0.0;
        Double_t x_val              = 0.0;
        Double_t Xlow_val           = 0.0;
        Double_t Xhigh_val          = 0.0;
        Double_t Ylow_val_envelope  = 0.0;
        Double_t Yhigh_val_envelope = 0.0;
        Double_t Ylow_val_stat      = 0.0;
        Double_t Yhigh_val_stat     = 0.0;
        for(Int_t i_file = 0; i_file < size; i_file++)
        {
            Double_t x,y;
            Double_t Xhigh,Xlow,Yhigh,Ylow;

            vec_tgae[i_file]->GetPoint(i_point,x,y);
            Xhigh = vec_tgae[i_file]->GetErrorXhigh(i_point);
            Xlow  = vec_tgae[i_file]->GetErrorXlow(i_point);
            Yhigh = vec_tgae[i_file]->GetErrorYhigh(i_point);
            Ylow  = vec_tgae[i_file]->GetErrorYlow(i_point);

            if(i_file == 0)
            {
                Ylow_val_envelope  = y;
                Yhigh_val_envelope = y;
            }
            else
            {
                if(y > Yhigh_val_envelope) Yhigh_val_envelope = y;
                if(y < Ylow_val_envelope)  Ylow_val_envelope  = y;
            }

            if(Yhigh > 0.0)
            {
                Double_t weight = (Yhigh+Ylow)/2.0;
                weight_average += weight;
                average_y_val  += y*weight;

                Ylow_val_stat  += Ylow*weight;
                Yhigh_val_stat += Yhigh*weight;
            }
            x_val     = x;
            Xlow_val  = Xlow;
            Xhigh_val = Xhigh;
        }
        if(weight_average > 0.0)
        {
            average_y_val  /= weight_average;
            Ylow_val_stat  /= weight_average;
            Yhigh_val_stat /= weight_average;

            Yhigh_val_envelope -= average_y_val;
            Ylow_val_envelope  -= average_y_val;
        }

        //cout << "i_point: " << i_point << ", x_val: " << x_val << ", y_val: " << average_y_val << endl;

        tgae_central ->SetPoint(i_point,x_val,average_y_val);
        tgae_central ->SetPointError(i_point,Xlow_val,Xhigh_val,Ylow_val_stat,Yhigh_val_stat);

        tgae_envelope ->SetPoint(i_point,x_val,average_y_val);
        tgae_envelope ->SetPointError(i_point,Xlow_val,Xhigh_val,fabs(Ylow_val_envelope),fabs(Yhigh_val_envelope));
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void DrawBoxErrors(TGraphAsymmErrors* tgae_in, Int_t color, Double_t width, Int_t style,
                  Double_t x_offset)
{
    for(Int_t i_point = 0; i_point < tgae_in->GetN(); i_point++)
    {
        Double_t x,y;
        Double_t Xhigh,Xlow,Yhigh,Ylow;

        tgae_in->GetPoint(i_point,x,y);
        Xhigh = tgae_in->GetErrorXhigh(i_point);
        Xlow  = tgae_in->GetErrorXlow(i_point);
        Yhigh = tgae_in->GetErrorYhigh(i_point);
        Ylow  = tgae_in->GetErrorYlow(i_point);

        if(Yhigh != 0.0 && Ylow != 0)
        {
            TBox* box = new TBox();
            box->SetFillColor(color);
            box->SetFillStyle(style);
            box->SetX1(x_offset+x-width/2.0);
            box->SetX2(x_offset+x+width/2.0);
            box->SetY1(y+Yhigh);
            box->SetY2(y-Ylow);
            box->Draw("same");
        }
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t EPD_ADCfit(double *x, double *par)
{
return  par[1]*TMath::Exp(-x[0]/par[0])
        +par[5]*(
		par[7]*TMath::Poisson((x[0]-par[6])/par[2],par[3])
		/*
        +TMath::Poisson(2,par[4])*TMath::Poisson((x[0]-par[6])/par[2],2*par[3])
        +TMath::Poisson(3,par[4])*TMath::Poisson((x[0]-par[6])/par[2],3*par[3])
        +TMath::Poisson(4,par[4])*TMath::Poisson((x[0]-par[6])/par[2],4*par[3])
        +TMath::Poisson(5,par[4])*TMath::Poisson((x[0]-par[6])/par[2],5*par[3])
        +TMath::Poisson(6,par[4])*TMath::Poisson((x[0]-par[6])/par[2],6*par[3])
        +TMath::Poisson(7,par[4])*TMath::Poisson((x[0]-par[6])/par[2],7*par[3])
        +TMath::Poisson(8,par[4])*TMath::Poisson((x[0]-par[6])/par[2],8*par[3])
        +TMath::Poisson(9,par[4])*TMath::Poisson((x[0]-par[6])/par[2],9*par[3])
        +TMath::Poisson(10,par[4])*TMath::Poisson((x[0]-par[6])/par[2],10*par[3])
        +TMath::Poisson(11,par[4])*TMath::Poisson((x[0]-par[6])/par[2],11*par[3])
        +TMath::Poisson(12,par[4])*TMath::Poisson((x[0]-par[6])/par[2],12*par[3])
        +TMath::Poisson(13,par[4])*TMath::Poisson((x[0]-par[6])/par[2],13*par[3])
        +TMath::Poisson(14,par[4])*TMath::Poisson((x[0]-par[6])/par[2],14*par[3])
        +TMath::Poisson(15,par[4])*TMath::Poisson((x[0]-par[6])/par[2],15*par[3])
		*/
        );
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 BiotSavartFunc(Double_t current, TVector3 dl, TVector3 B_pos)
{
    // Calculated the B-field "B_field_vector" at position "B_pos" for a infinitesimal small conductor with length "dl" and current "current"
    TVector3 B_field_vector = dl.Cross(B_pos);
    Double_t scalar_constant = (mu_zero/(4.0*TMath::Pi()))*current/TMath::Power(B_pos.Mag(),3.0);
    B_field_vector *= scalar_constant;

    return B_field_vector;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 HelixPoint(Double_t t, Double_t radius, Double_t pitch)
{
    Double_t x = radius*TMath::Cos(t);
    Double_t y = radius*TMath::Sin(t);
    Double_t z = pitch*t;
    TVector3 helix_point(x,y,z);
    return helix_point;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// RMS of a Graph of x-axis with y-values as weights
Double_t GetRMS_Graph(TGraphAsymmErrors* tgraph, Double_t x_start, Double_t x_stop, Double_t &RMS_error)
{
    Int_t fNpoints = tgraph ->GetN();
    if (fNpoints <= 0) return 0;

    // Three loops are needed
    Double_t sumx = 0, sumx2 = 0, sum_weigh_x = 0, sum_weighs = 0;
    for (Int_t i = 0; i < fNpoints; i++)
    {
        Double_t x_val, y_val;
        tgraph ->GetPoint(i,x_val,y_val);
        Double_t exl = tgraph->GetErrorXlow(i);
        Double_t exh = tgraph->GetErrorXhigh(i);
        Double_t eyl = tgraph->GetErrorYlow(i);
        Double_t eyh = tgraph->GetErrorYhigh(i);
        Double_t ey_err = (eyl+eyh)/2.0;
        if(x_val >= x_start && x_val <= x_stop)
        {
            sumx            += y_val * x_val;
            sumx2           += y_val * x_val * x_val;
            sum_weigh_x     += y_val*x_val;
            sum_weighs      += y_val;
        }
    }
    if(sum_weighs <= 0) return 0;
    Double_t mean     = sumx / sum_weighs;
    Double_t rms2     = TMath::Abs(sumx2 / sum_weighs - mean * mean);
    Double_t rms      = TMath::Sqrt(rms2);

    // Second loop to calculate mean and rms2 error
    Double_t mean_err = 0, rms2_err = 0;
    for (Int_t i = 0; i < fNpoints; i++)
    {
        Double_t x_val, y_val;
        tgraph ->GetPoint(i,x_val,y_val);
        Double_t exl = tgraph->GetErrorXlow(i);
        Double_t exh = tgraph->GetErrorXhigh(i);
        Double_t eyl = tgraph->GetErrorYlow(i);
        Double_t eyh = tgraph->GetErrorYhigh(i);
        Double_t ey_err = (eyl+eyh)/2.0;
        if(x_val >= x_start && x_val <= x_stop)
        {
            mean_err += TMath::Power(((x_val*sum_weighs       - sum_weigh_x)/(sum_weighs*sum_weighs) * ey_err),2.0);
            rms2_err += TMath::Power(((x_val*x_val*sum_weighs - sumx2)/(sum_weighs*sum_weighs) * ey_err),2.0);
        }
    }
    Double_t rms_err = 0;
    rms2_err += mean_err; // mean_err still squared
    mean_err  = TMath::Sqrt(mean_err);
    rms_err   = rms2_err/(2.0*TMath::Sqrt(rms2)); // error calculation for sqrt(y) = dy*0.5/sqrt(y)
    RMS_error = rms2_err;

    return rms;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------
// Calculate integral for a TGraphAsymmErrors with constant bin widht
// Bin width is automatically estimated (!) and used for the integral
void Calc_Integral_tgae(TGraphAsymmErrors* tgae_in, Double_t x_start, Double_t x_stop, Double_t& integral, Double_t& integral_error)
{
    integral       = 0.0;
    integral_error = 0.0;
    Int_t N_points = tgae_in->GetN();
    if(N_points >= 2)
    {
        for(Int_t i_point = 0; i_point <  N_points; i_point++)
        {
            Double_t x_val,y_val;
            tgae_in ->GetPoint(i_point,x_val,y_val);
            if(x_val < x_start || x_val > x_stop) continue;
            Double_t eyl = tgae_in ->GetErrorYlow(i_point);
            Double_t eyh = tgae_in ->GetErrorYhigh(i_point);
            Double_t bin_width_tgae = 0;
            if(i_point == 0) // First point
            {
                Double_t x_val1,y_val1;
                tgae_in ->GetPoint(i_point+1,x_val1,y_val1);
                bin_width_tgae = fabs(x_val1 - x_val)/1.0;
            }
            if(i_point > 0 && i_point < (N_points-1)) // In between
            {
                Double_t x_val1,y_val1;
                tgae_in ->GetPoint(i_point-1,x_val1,y_val1);
                Double_t x_val2,y_val2;
                tgae_in ->GetPoint(i_point+1,x_val2,y_val2);
                bin_width_tgae = fabs(x_val2 - x_val)/2.0;
            }
            if(i_point == (N_points-1)) // Last point
            {
                Double_t x_val2,y_val2;
                tgae_in ->GetPoint(i_point-1,x_val2,y_val2);
                bin_width_tgae = fabs(x_val2 - x_val)/1.0;
            }
            //cout << "x_val: " << x_val << ", y_val: " << y_val << ", bin_width_tgae: " << bin_width_tgae << endl;
            integral       += y_val*bin_width_tgae;
            integral_error += TMath::Power(bin_width_tgae*(eyl+eyh/2.0),2.0);
        }
        integral_error = TMath::Sqrt(integral_error);
        //cout << "bin_width_tgae: " << bin_width_tgae << ", integral: " << integral << endl;
    }
}
//----------------------------------------------------------



//----------------------------------------------------------------------------------------
// Scale a TGraphAsymmErrors in width with scale_fac and center x_center
TGraphAsymmErrors* ScaleWidth_Graph(TGraphAsymmErrors* tgraph, Double_t x_center, Double_t scale_fac)
{
    TGraphAsymmErrors* tgraph_out = new TGraphAsymmErrors();
    Int_t fNpoints = tgraph ->GetN();

    for (Int_t i = 0; i < fNpoints; i++)
    {
        Double_t x_val, y_val;
        tgraph ->GetPoint(i,x_val,y_val);
        Double_t exl = tgraph->GetErrorXlow(i);
        Double_t exh = tgraph->GetErrorXhigh(i);
        Double_t eyl = tgraph->GetErrorYlow(i);
        Double_t eyh = tgraph->GetErrorYhigh(i);
        tgraph_out ->SetPoint(i,(x_val-x_center)*scale_fac + x_center,y_val);
        tgraph_out ->SetPointError(i,exl,exh,eyl,eyh);
    }
    return tgraph_out;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Projects a 2D histogram along a line
class Project_h2D
{
    TH2D* h2D_hist;
    TVector2 TV2_start, TV2_stop;
    Int_t N_bins;
    TH1D* h1D_proj;
public:
    void add_h2D(TH2D* h2D_hist_in);
    void define_line(TVector2 TV2_start_in, TVector2 TV2_stop_in);
    void set_N_bins(Int_t N_bins_in);
    TH1D* get_projection(TString name);
    TVector2 get_TV2_start();
    void clear();
    void clear_1D();
};

void Project_h2D::add_h2D(TH2D* h2D_hist_in)
{
    h2D_hist = (TH2D*)h2D_hist_in->Clone("h2D_hist_clone");
}

void Project_h2D::define_line(TVector2 TV2_start_in, TVector2 TV2_stop_in)
{
    TV2_start = TV2_start_in;
    TV2_stop  = TV2_stop_in;
}

void Project_h2D::set_N_bins(Int_t N_bins_in)
{
    N_bins = N_bins_in;
}

TH1D* Project_h2D::get_projection(TString name)
{
    Double_t x_val_start = TV2_start.Px();
    Double_t x_val_stop  = TV2_stop.Px();
    Double_t y_val_start = TV2_start.Py();
    Double_t y_val_stop  = TV2_stop.Py();

    // switch start and stop in case stop is smaller than start

    if(TV2_stop.Px() < x_val_start)
    {
        x_val_start = TV2_stop.Px();
        x_val_stop  = TV2_start.Px();

        y_val_start = TV2_stop.Py();
        y_val_stop  = TV2_start.Py();

        //TVector2 TV2_temp = TV2_stop;
        //TV2_stop  = TV2_start;
        //TV2_start = TV2_temp;
    }

    // Define vector from start to stop
    TVector2 TV2_length = TV2_stop - TV2_start;
    Double_t TV2_length_val = TMath::Sqrt(TV2_length.Px()*TV2_length.Px() + TV2_length.Py()*TV2_length.Py());
    Double_t phi = TV2_length.Phi(); // azimuthal angle
    //cout << "phi: " << phi << endl;

    // Flag which defines whether to do the projection along x or y axis
    Int_t flag_xy = 0;
    if(phi > TMath::Pi()/4.0 && phi <= 3.0*TMath::Pi()/4.0) flag_xy = 1;
    if(phi > 5.0*TMath::Pi()/4.0 && phi <= 7.0*TMath::Pi()/4.0) flag_xy = 1;

    Double_t length_proj = 0.0;
    if(flag_xy == 0) length_proj = fabs(TMath::Cos(phi))*TV2_length_val;
    if(flag_xy == 1) length_proj = fabs(TMath::Sin(phi))*TV2_length_val;

    Double_t h_start = 0.0;
    Double_t h_stop  = TV2_length.Mod();
    //Double_t h_start = -TV2_length.Mod()/2.0;
    //Double_t h_stop  = TV2_length.Mod()/2.0;

    Int_t N_2D_hist_bins_min = 0;
    if(flag_xy == 0)
    {
        Double_t x_val_start_bin = h2D_hist->GetXaxis()->GetBinCenter(1);
        for(Int_t i_bin_x = 1; i_bin_x <= h2D_hist->GetNbinsX(); i_bin_x++)
        {
            Double_t x_val = h2D_hist->GetXaxis()->GetBinCenter(i_bin_x);
            if(fabs(x_val-x_val_start_bin) < length_proj)
            {
                N_2D_hist_bins_min++;
            }
        }
    }

    if(flag_xy == 1)
    {
        Double_t y_val_start_bin = h2D_hist->GetYaxis()->GetBinCenter(1);
        for(Int_t i_bin_y = 1; i_bin_y <= h2D_hist->GetNbinsY(); i_bin_y++)
        {
            Double_t y_val = h2D_hist->GetYaxis()->GetBinCenter(i_bin_y);
            if(fabs(y_val-y_val_start_bin) < length_proj)
            {
                N_2D_hist_bins_min++;
            }
        }
    }

    N_2D_hist_bins_min--;
    //printf("N_bins: %d, N_2D_hist_bins_min: %d \n",N_bins,N_2D_hist_bins_min);

    if(N_bins > N_2D_hist_bins_min)
    {
        //cout << "WARNING in Project_h2D::get_projection, number of bins larger than 2D bins, set to: " << N_2D_hist_bins_min << endl;
        N_bins = N_2D_hist_bins_min;
    }

    h1D_proj = new TH1D(name.Data(),name.Data(),N_bins,h_start,h_stop);

    //cout << "h_start: " << h_start << ", h_stop: " << h_stop << endl;

    // Calculate slope and offset, check if the vector is going exactly into y-direction
    Int_t flag_perp = 0;
    Double_t m, b;
    if(TV2_stop.Px()-TV2_start.Px() == 0.0)
    {
        m = 1.0;
        b = 0.0;
        flag_perp = 1;
    }
    else
    {
        m = (TV2_stop.Py()-TV2_start.Py())/(TV2_stop.Px()-TV2_start.Px());
        b = TV2_start.Py() - m*TV2_start.Px();
    }

    if(flag_xy == 0) // go along x-axis to do projection
    {
        for(Int_t i_bin_x = 1; i_bin_x <= h2D_hist->GetNbinsX(); i_bin_x++)
        {
            Double_t x_val = h2D_hist->GetXaxis()->GetBinCenter(i_bin_x);
            if(x_val >= x_val_start && x_val < x_val_stop)
            {
                Double_t y_val = m*x_val + b;
                TVector2 TV2_pos(x_val,y_val);
                TVector2 TV2_pos_length = TV2_pos - TV2_start;
                Double_t h1D_axis_val = TV2_pos_length.Mod();

                Int_t i_bin_y = h2D_hist->GetYaxis()->FindBin(y_val);
                Double_t bin_cont = h2D_hist->GetBinContent(i_bin_x,i_bin_y);
                //Int_t h1D_bin = h1D_proj->FindBin(h1D_axis_val-TV2_length.Mod()/2.0);
                Int_t h1D_bin = h1D_proj->FindBin(h1D_axis_val);
                h1D_proj->SetBinContent(h1D_bin,bin_cont);
                //cout << "h1D_axis_val: " << h1D_axis_val << ", h1D_bin: " << h1D_bin << ", bin_cont: " << bin_cont << endl;
            }
        }
    }
    else  // go along y-axis to do projection
    {
        for(Int_t i_bin_y = 1; i_bin_y <= h2D_hist->GetNbinsY(); i_bin_y++)
        {
            Double_t y_val = h2D_hist->GetYaxis()->GetBinCenter(i_bin_y);
            if(y_val_start > y_val_stop)
            {
                Double_t y_val_temp  = y_val_start;
                y_val_start = y_val_stop;
                y_val_stop  = y_val_temp;
            }
            if(y_val >= y_val_start && y_val < y_val_stop)
            {
                Double_t x_val;
                if(!flag_perp)
                {
                    x_val = (y_val - b)/m;
                }
                else
                {
                    x_val = TV2_start.Px();
                }
                TVector2 TV2_pos(x_val,y_val);
                TVector2 TV2_pos_length = TV2_pos - TV2_start;
                Double_t h1D_axis_val = TV2_pos_length.Mod();

                Int_t i_bin_x = h2D_hist->GetXaxis()->FindBin(x_val);
                Double_t bin_cont = h2D_hist->GetBinContent(i_bin_x,i_bin_y);
                //Int_t h1D_bin = h1D_proj->FindBin(h1D_axis_val-TV2_length.Mod()/2.0);
                Int_t h1D_bin = h1D_proj->FindBin(h1D_axis_val);
                h1D_proj->SetBinContent(h1D_bin,bin_cont);
                //cout << "h1D_axis_val: " << h1D_axis_val << ", h1D_bin: " << h1D_bin << ", bin_cont: " << bin_cont << endl;
            }
        }
    }

    return h1D_proj;
}

TVector2 Project_h2D::get_TV2_start()
{
    return TV2_start;
}

void Project_h2D::clear()
{
    delete h2D_hist;
}

void Project_h2D::clear_1D()
{
    delete h1D_proj;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_2D_histo_and_canvas(TH2D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.22);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(max_val > min_val)
    {
        hist->GetZaxis()->SetRangeUser(min_val,max_val);
    }
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_histo_and_canvas(TH1D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_graph_and_canvas(TGraph* hist, TString name, Int_t x_size, Int_t y_size,
                                  Double_t min_val, Double_t max_val, TString option,
                                  Int_t style, Double_t size, Int_t color)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetTitle("");
    hist->SetMarkerSize(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(style);
    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    hist->Draw(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_profile_and_canvas(TProfile* hist, TString name, Int_t x_size, Int_t y_size,
                                  Double_t min_val, Double_t max_val, TString option,
                                  Int_t style, Double_t size, Int_t color)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->SetMarkerSize(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(style);
    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    hist->Draw(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TH2D* Subtract_plane_from_h2D(TH2D* h2D_in, TVector3 vec_x, TVector3 vec_y, Double_t min_elevation_ext)
{
    vec_x *= 1.0/vec_x.Mag();
    vec_y *= 1.0/vec_y.Mag();
    Double_t min_elevation = 100000.0;
    Double_t min_x_y_pos[2];
    for(Int_t i_bin_x = 1; i_bin_x <= h2D_in->GetNbinsX(); i_bin_x++)
    {
        for(Int_t i_bin_y = 1; i_bin_y <= h2D_in->GetNbinsY(); i_bin_y++)
        {
            Double_t bin_cont = h2D_in->GetBinContent(i_bin_x,i_bin_y);
            Double_t val_x = h2D_in->GetXaxis()->GetBinCenter(i_bin_x);
            Double_t val_y = h2D_in->GetYaxis()->GetBinCenter(i_bin_y);
            if(bin_cont > 0.0 && bin_cont < min_elevation)
            {
                min_elevation = bin_cont;
                min_x_y_pos[0] = val_x;
                min_x_y_pos[1] = val_y;
            }
        }
    }

    min_elevation = min_elevation_ext; // use the minimum value from extern
    TVector3 TV3_min;
    TV3_min.SetXYZ(min_x_y_pos[0],min_x_y_pos[1],min_elevation);

    TH2D* h2D_in_sub = (TH2D*)h2D_in->Clone("h2D_in_sub");

    for(Int_t i_bin_x = 1; i_bin_x <= h2D_in->GetNbinsX(); i_bin_x++)
    {
        for(Int_t i_bin_y = 1; i_bin_y <= h2D_in->GetNbinsY(); i_bin_y++)
        {
            Double_t bin_cont = h2D_in->GetBinContent(i_bin_x,i_bin_y);
            Double_t val_x = h2D_in->GetXaxis()->GetBinCenter(i_bin_x);
            Double_t val_y = h2D_in->GetYaxis()->GetBinCenter(i_bin_y);
            //if(bin_cont > 0.0)
            {
                //Double_t x_rel_val = val_x - min_x_y_pos[0];
                //Double_t y_rel_val = val_y - min_x_y_pos[1];
                Double_t x_rel_val = val_x;
                Double_t y_rel_val = val_y;

                //x_rel_val = TV3_min.X() + 10.0;
                //y_rel_val = TV3_min.Y();

                Double_t b_partA = (y_rel_val - TV3_min.Y())/vec_y.Y() + vec_x.Y()*(TV3_min.X() - x_rel_val)/(vec_x.X()*vec_y.Y());
                Double_t b_partB = 1.0 - (vec_y.X()*vec_x.Y())/(vec_x.X()*vec_y.Y());
                Double_t b = b_partA/b_partB;
                Double_t a = (x_rel_val-b*vec_y.X() - TV3_min.X())/vec_x.X();
                Double_t z = a*vec_x.Z() + b*vec_y.Z() + TV3_min.Z();
                //cout << "x_rel_val: " << x_rel_val << ", y_rel_val: " << y_rel_val << ", z: " << z << ", z_min: " << TV3_min.Z()
                //    << ", vec_x.Z(): " << vec_x.Z() << ", vec_y.Z(): " << vec_y.Z() << ", a: " << a << endl;

                Double_t new_bin_cont = bin_cont - z;
                //if(bin_cont == 0.0) new_bin_cont = -100000.0;
                h2D_in_sub->SetBinContent(i_bin_x,i_bin_y,new_bin_cont);
            }
        }
    }

    cout << "Plane subtracted" << endl;
    return h2D_in_sub;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Project_map_along_angles
{
    TH2D* h2D_map;
    std::vector<TH1D*> vec_h1D_proj;
    std::vector<Double_t> vec_angles;
    TVector2 TV2_center;
    std::vector< std::vector<TVector2> > TV2_proj_lines;
    Int_t N_bins;
    Double_t proj_length;
    TPolyMarker* PM_start;
    std::vector<TVector2> TV2_start_vectors;
public:
    void add_map(TH2D* h2D_map_in);
    void add_angle(Double_t angle_in); // in radiants
    void set_proj_length(Double_t proj_length_in);
    void set_center(TVector2 TV2_center_in);
    void set_N_bins(Int_t N_bins_in);
    std::vector<TH1D*> get_projections(TString name);
    void plot_lines(Int_t col, Double_t width, Int_t style);
    void plot_start_points(Int_t col, Double_t size, Int_t style);
    void plot_proj_index(Double_t size);
    void clear();
};

void Project_map_along_angles::add_map(TH2D* h2D_map_in)
{
    proj_length = 1.0;
    h2D_map = (TH2D*)h2D_map_in->Clone("h2D_map");
}

void Project_map_along_angles::add_angle(Double_t angle_in)
{
    vec_angles.push_back(angle_in);
}

void Project_map_along_angles::set_proj_length(Double_t proj_length_in)
{
    proj_length = proj_length_in;
}

void Project_map_along_angles::set_center(TVector2 TV2_center_in)
{
    TV2_center = TV2_center_in;
}

void Project_map_along_angles::set_N_bins(Int_t N_bins_in)
{
    N_bins = N_bins_in;
}

std::vector<TH1D*> Project_map_along_angles::get_projections(TString name)
{
    PM_start = new TPolyMarker();
    TV2_proj_lines.resize(2); // [start, stop]
    Project_h2D Proj_h2D;
    Proj_h2D.add_h2D(h2D_map);
    TVector2 TV2_start, TV2_stop, TV2_dir;
    for(Int_t i_angle = 0; i_angle < vec_angles.size(); i_angle++)
    {
        TV2_dir.Set(TMath::Cos(vec_angles[i_angle]),TMath::Sin(vec_angles[i_angle]));
        TV2_dir *= proj_length;
        TV2_start = TV2_dir;
        TV2_start += TV2_center;
        TV2_stop = TV2_center;
        TV2_stop -= TV2_dir;
        Proj_h2D.define_line(TV2_start,TV2_stop);
        Proj_h2D.set_N_bins(N_bins);
        HistName = "vec_h1D_proj_";
        HistName += name;
        HistName += "_";
        HistName += i_angle;
        HistNameB = HistName;
        HistNameB += "B";
        vec_h1D_proj.push_back((TH1D*)(Proj_h2D.get_projection(HistNameB.Data())->Clone(HistName.Data())));
        TV2_proj_lines[0].push_back(TV2_start);
        TV2_proj_lines[1].push_back(TV2_stop);
        TVector2 TV2_start_used = Proj_h2D.get_TV2_start();
        PM_start->SetNextPoint(TV2_start_used.Px(),TV2_start_used.Py());
        TV2_start_vectors.push_back(TV2_start_used);
    }
    Proj_h2D.clear();
    return vec_h1D_proj;
}

void Project_map_along_angles::plot_lines(Int_t col, Double_t width, Int_t style)
{
    for(Int_t i_proj = 0; i_proj < TV2_proj_lines[0].size(); i_proj++)
    {
        PlotLine(TV2_proj_lines[0][i_proj].Px(),TV2_proj_lines[1][i_proj].Px(),TV2_proj_lines[0][i_proj].Py(),TV2_proj_lines[1][i_proj].Py(),col,width,style); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    }
}

void Project_map_along_angles::plot_start_points(Int_t col, Double_t size, Int_t style)
{
    PM_start->SetMarkerStyle(style);
    PM_start->SetMarkerColor(col);
    PM_start->SetMarkerSize(size);
    PM_start->Draw();
}

void Project_map_along_angles::plot_proj_index(Double_t size)
{
    for(Int_t i_angle = 0; i_angle < TV2_start_vectors.size(); i_angle++)
    {
        HistName = "";
        HistName += i_angle;
        plotTopLegend((char*)HistName.Data(),TV2_start_vectors[i_angle].Px(),TV2_start_vectors[i_angle].Py(),size,1,0.0,42,0,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    }
}

void Project_map_along_angles::clear()
{
    delete PM_start;
    vec_angles.clear();
    vec_h1D_proj.clear();
    TV2_proj_lines.clear();
    TV2_start_vectors.clear();
}
//----------------------------------------------------------------------------------------






//----------------------------------------------------------------------------------------
class Calc_Map_height
{
    TH2D* h2D_map;
    std::vector<Double_t> vec_height;
    std::vector<TPolyMarker3D*> vec_pm;
    Double_t x_val_min, x_val_max, y_val_min, y_val_max;
public:
    void add_map(TH2D* h2D_map_in);
    void set_height_value(Double_t height_in);
    void set_xy_limits(Double_t x_val_min_in, Double_t x_val_max_in, Double_t y_val_min_in, Double_t y_val_max_in);
    void calc_iso_points();
    void draw_points(Int_t i_height, Double_t size, Int_t color, Int_t style);
    void draw_points_2D(Int_t i_height, Double_t size, Int_t color, Int_t style);
    void clear();
};

void Calc_Map_height::add_map(TH2D* h2D_map_in)
{
    h2D_map = (TH2D*)h2D_map_in->Clone("h2D_map");
}

void Calc_Map_height::set_height_value(Double_t height_in)
{
    vec_height.push_back(height_in);
}

void Calc_Map_height::set_xy_limits(Double_t x_val_min_in, Double_t x_val_max_in, Double_t y_val_min_in, Double_t y_val_max_in)
{
    x_val_min = x_val_min_in;
    x_val_max = x_val_max_in;
    y_val_min = y_val_min_in;
    y_val_max = y_val_max_in;
}

void Calc_Map_height::calc_iso_points()
{
    vec_pm.resize(vec_height.size());
    for(Int_t i_height = 0; i_height < vec_height.size(); i_height++)
    {
        TPolyMarker3D* pm_single = new TPolyMarker3D();
        vec_pm[i_height] = (TPolyMarker3D*)pm_single->Clone();
        Double_t bin_width_x = h2D_map->GetXaxis()->GetBinWidth(1);
        Double_t bin_width_y = h2D_map->GetYaxis()->GetBinWidth(1);
        for(Int_t i_bin_x  = 1; i_bin_x <= h2D_map->GetNbinsX(); i_bin_x++)
        {
            Double_t height_below = 0.0;
            Double_t height_above = 0.0;
            Int_t flag_below_above = 0;

            Double_t x_val = h2D_map->GetXaxis()->GetBinCenter(i_bin_x);
            for(Int_t i_bin_y  = 1; i_bin_y <= h2D_map->GetNbinsY(); i_bin_y++)
            {
                Double_t y_val = h2D_map->GetYaxis()->GetBinCenter(i_bin_y);
                Double_t bin_cont = h2D_map->GetBinContent(i_bin_x,i_bin_y);
                if(i_bin_y == 1)
                {
                    if(bin_cont <= vec_height[i_height])
                    {
                        height_below = bin_cont;
                        flag_below_above = 0;
                    }
                    else
                    {
                        height_above = bin_cont;
                        flag_below_above = 1;
                    }
                }
                else
                {
                    if(bin_cont <= vec_height[i_height] && flag_below_above == 1)
                    {
                        flag_below_above = 0;
                        if(x_val > x_val_min && x_val < x_val_max && y_val > y_val_min && y_val < y_val_max)
                        {
                            vec_pm[i_height]->SetNextPoint(x_val,y_val-bin_width_y/2.0,vec_height[i_height]);
                        }
                    }
                    if(bin_cont > vec_height[i_height] && flag_below_above == 0)
                    {
                        flag_below_above = 1;
                        if(x_val > x_val_min && x_val < x_val_max && y_val > y_val_min && y_val < y_val_max)
                        {
                            vec_pm[i_height]->SetNextPoint(x_val,y_val-bin_width_y/2.0,vec_height[i_height]);
                        }
                    }
                }
            }
        }
    }
}

void Calc_Map_height::draw_points(Int_t i_height, Double_t size, Int_t color, Int_t style)
{
    std::vector<TPolyMarker3D*> vec_pm_clone;
    vec_pm_clone.resize(vec_pm.size());
    for(Int_t i_height = 0; i_height < vec_pm.size(); i_height++)
    {
        vec_pm_clone[i_height] = (TPolyMarker3D*)vec_pm[i_height]->Clone();
    }

    vec_pm_clone[i_height]->SetMarkerStyle(style);
    vec_pm_clone[i_height]->SetMarkerColor(color);
    vec_pm_clone[i_height]->SetMarkerSize(size);
    vec_pm_clone[i_height]->Draw();
}

void Calc_Map_height::draw_points_2D(Int_t i_height, Double_t size, Int_t color, Int_t style)
{
    std::vector<TPolyMarker*> vec_pm_2D_clone;
    vec_pm_2D_clone.resize(vec_pm.size());
    for(Int_t i_height = 0; i_height < vec_pm.size(); i_height++)
    {
        vec_pm_2D_clone[i_height] = new TPolyMarker();

        for(Int_t i_point = 0; i_point < vec_pm[i_height]->GetN(); i_point++)
        {
            Double_t x,y,z;
            vec_pm[i_height]->GetPoint(i_point,x,y,z);
            vec_pm_2D_clone[i_height] ->SetNextPoint(x,y);
        }
    }

    vec_pm_2D_clone[i_height]->SetMarkerStyle(style);
    vec_pm_2D_clone[i_height]->SetMarkerColor(color);
    vec_pm_2D_clone[i_height]->SetMarkerSize(size);
    vec_pm_2D_clone[i_height]->Draw();
}

void Calc_Map_height::clear()
{
    vec_pm.clear();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Calc_rebin_TPM3D
{
    TPolyMarker3D* TPM3D_in;
    TPolyMarker3D* TPM3D_out;
    Int_t Grid_sizeX, Grid_sizeY;
    Double_t min_max_xy[2][2]; // [min,max][x,y]
public:
    void add_TPM3D(TPolyMarker3D* TPM3D_in_a);
    void set_Grid_size_XY(Int_t Grid_sizeX_a, Int_t Grid_sizeY_a);
    TPolyMarker3D* get_TPM3D_Grid_size();
};

void Calc_rebin_TPM3D::add_TPM3D(TPolyMarker3D* TPM3D_in_a)
{
    TPM3D_in = (TPolyMarker3D*)TPM3D_in_a->Clone();
}

void Calc_rebin_TPM3D::set_Grid_size_XY(Int_t Grid_sizeX_a, Int_t Grid_sizeY_a)
{
    Grid_sizeX = Grid_sizeX_a;
    Grid_sizeY = Grid_sizeY_a;
}

TPolyMarker3D* Calc_rebin_TPM3D::get_TPM3D_Grid_size()
{
    TPM3D_out = new TPolyMarker3D();

    // determine x and y range of data
    for(Int_t i_point = 0; i_point < TPM3D_in->GetN(); i_point++)
    {
        Double_t x,y,z;
        TPM3D_in->GetPoint(i_point,x,y,z);
        if(i_point == 0)
        {
            min_max_xy[0][0] = x;
            min_max_xy[1][0] = x;
            min_max_xy[0][1] = y;
            min_max_xy[1][1] = y;
        }
        else
        {
            if(x < min_max_xy[0][0]) min_max_xy[0][0] = x;
            if(x > min_max_xy[1][0]) min_max_xy[1][0] = x;
            if(y < min_max_xy[0][1]) min_max_xy[0][1] = y;
            if(y > min_max_xy[1][1]) min_max_xy[1][1] = y;
        }
    }

    Int_t N_bins_X = (Int_t)((min_max_xy[1][0] - min_max_xy[0][0])/Grid_sizeX);
    Int_t N_bins_Y = (Int_t)((min_max_xy[1][1] - min_max_xy[0][1])/Grid_sizeY);

    for(Int_t i_X = 0; i_X < N_bins_X; i_X++)
    {
        Double_t start_x = min_max_xy[0][0] + (Double_t)i_X*Grid_sizeX;
        Double_t stop_x  = min_max_xy[0][0] + (Double_t)i_X*Grid_sizeX + Grid_sizeX;

        for(Int_t i_Y = 0; i_Y < N_bins_Y; i_Y++)
        {
            Double_t start_y = min_max_xy[0][1] + (Double_t)i_Y*Grid_sizeY;
            Double_t stop_y  = min_max_xy[0][1] + (Double_t)i_Y*Grid_sizeY + Grid_sizeY;

            Double_t mean_x   = 0.0;
            Double_t mean_y   = 0.0;
            Int_t    n_points = 0;

            for(Int_t i_point = 0; i_point < TPM3D_in->GetN(); i_point++)
            {
                Double_t x,y,z;
                TPM3D_in->GetPoint(i_point,x,y,z);

                if(x > start_x && x <= stop_x && y > start_y && y <= stop_y)
                {
                    n_points++;
                    mean_x += x;
                    mean_y += y;
                }
            }
            if(n_points > 0)
            {
                mean_x /= (Double_t)n_points;
                mean_y /= (Double_t)n_points;
                TPM3D_out ->SetNextPoint(mean_x,mean_y,0.0);
            }
        }
    }

    return TPM3D_out;
}
//----------------------------------------------------------------------------------------



#if 1
//----------------------------------------------------------------------------------------
class Calc_radial_distribution
{
    TH2D* TH2D_in;
    std::vector< std::vector<Double_t> > vec_data;
    TPolyMarker3D* TPM3D_data;
    TVector2 TV2_center, TV2_dir;
    TH1D* TH1D_radial_projection;
    TH1D* TH1D_radial_projection_entries;
    Int_t flag_input_data; // 0 = TH2D, 1 = vec_data, 2 = TPolyMarker3D
public:
    void  add_TH2D(TH2D* TH2D_in_a);
    void  add_vector_data(std::vector< std::vector<Double_t> > vec_data_a);
    void  add_TPM3D_data(TPolyMarker3D* TPM3D_data_a);
    void  set_center_and_dir(TVector2 TV2_center_a, TVector2 TV2_dir_a);
    void  set_hist_dimensions(Int_t n_bins_a, Double_t start_a, Double_t stop_a);
    void  calc_radial_distribution();
    void  set_invert();
    TH1D* get_radial_distribution();
    void  set_Errors_to_zero();
    void  set_Zero_level();
    void  clear();
};

void Calc_radial_distribution::add_TH2D(TH2D* TH2D_in_a)
{
    TH2D_in = (TH2D*)TH2D_in_a->Clone("TH2D_in");
    flag_input_data = 0;
}

void Calc_radial_distribution::add_vector_data(std::vector< std::vector<Double_t> > vec_data_a)
{
    vec_data = vec_data_a;
    flag_input_data = 1;
}

void Calc_radial_distribution::add_TPM3D_data(TPolyMarker3D* TPM3D_data_a)
{
    TPM3D_data = (TPolyMarker3D*)TPM3D_data_a->Clone();
    flag_input_data = 2;
}

void Calc_radial_distribution::set_center_and_dir(TVector2 TV2_center_a, TVector2 TV2_dir_a)
{
    TV2_center = TV2_center_a;
    TV2_dir    = TV2_dir_a;
}

void Calc_radial_distribution::set_hist_dimensions(Int_t n_bins_a, Double_t start_a, Double_t stop_a)
{
    TH1D_radial_projection         = new TH1D("TH1D_radial_projection","TH1D_radial_projection",n_bins_a,start_a,stop_a);
    TH1D_radial_projection_entries = new TH1D("TH1D_radial_projection_entries","TH1D_radial_projection_entries",n_bins_a,start_a,stop_a);
}

void Calc_radial_distribution::calc_radial_distribution()
{
    Double_t slope = 0.0;
    Double_t Delta_x = TV2_dir.Px();
    Double_t Delta_y = TV2_dir.Py();
    if(Delta_x == 0.0) slope = 10000000.0;
    else
    {
        slope = Delta_y/Delta_x;
    }
    //y = mx + b
    Double_t b_val = TV2_center.Py() - slope*TV2_center.Px();

    //cout << "flag_input_data: " << flag_input_data << endl;
    //-------------------------------------------------------
    if(flag_input_data == 0) // TH2D input
    {
        Double_t min_val = TH2D_in->GetBinContent(TH2D_in->GetMinimumBin());
        for(Int_t bin_x = 1; bin_x <= TH2D_in->GetNbinsX(); bin_x++)
        {
            for(Int_t bin_y = 1; bin_y <= TH2D_in->GetNbinsY(); bin_y++)
            {
                Double_t x_val    = TH2D_in->GetXaxis()->GetBinCenter(bin_x);
                Double_t y_val    = TH2D_in->GetYaxis()->GetBinCenter(bin_y);
                Double_t bin_cont = TH2D_in->GetBinContent(bin_x,bin_y) ;
                Double_t distance_to_center = TMath::Sqrt(TMath::Power(x_val-TV2_center.Px(),2) + TMath::Power(y_val-TV2_center.Py(),2));
                Double_t bin_width = TH1D_radial_projection ->GetBinWidth(TH1D_radial_projection->FindBin(distance_to_center));
                //Double_t area = TMath::Pi()*(TMath::Power(distance_to_center+bin_width,2) - TMath::Power(distance_to_center,2));
                //Double_t area = distance_to_center;
                Double_t left_right_fact = 1.0;
                if(y_val < slope*x_val + b_val) left_right_fact = -1.0;
                Int_t bin = TH1D_radial_projection->FindBin(left_right_fact*distance_to_center);
                Double_t old_bin_cont = TH1D_radial_projection ->GetBinContent(bin);
                //if(area == 0.0)
                //{
                //    area = 1.0;
                //    cout << "WARNING! area = 0.0" << endl;
                //}
                //area = 1.0;
                //Double_t new_bin_cont = old_bin_cont + bin_cont/area;
                Double_t new_bin_cont = old_bin_cont + bin_cont;
                TH1D_radial_projection         ->SetBinContent(bin,new_bin_cont);

                // Analytic area calculation not precise engoug (bin center). Normalize directly by number of bins used.
                TH1D_radial_projection_entries ->SetBinContent(bin,TH1D_radial_projection_entries->GetBinContent(bin)+1.0);
                //if(bin == 100) cout << "new_bin_cont: " << new_bin_cont << ", area: " << area << endl;
                //if(bin == 105) cout << "105 new_bin_cont: " << new_bin_cont << ", area: " << area << endl;
                //TH1D_radial_projection ->Fill(left_right_fact*distance_to_center,bin_cont/area);
            }
        }
        TH1D_radial_projection->Divide(TH1D_radial_projection_entries);
    }
    //-------------------------------------------------------


    //-------------------------------------------------------
    if(flag_input_data == 1) // vector data input
    {
        cout << "data size: " << vec_data.size() << endl;
        for(Int_t i_data = 0; i_data < vec_data.size(); i_data++)
        {
            Double_t x_val    = vec_data[i_data][1];
            Double_t y_val    = vec_data[i_data][0];
            Double_t distance_to_center = TMath::Sqrt(TMath::Power(x_val-TV2_center.Px(),2) + TMath::Power(y_val-TV2_center.Py(),2));
            Double_t bin_width = TH1D_radial_projection ->GetBinWidth(TH1D_radial_projection->FindBin(distance_to_center));
            Double_t area = TMath::Pi()*(TMath::Power(distance_to_center+bin_width,2) - TMath::Power(distance_to_center,2));
            Double_t left_right_fact = 1.0;
            if(y_val < slope*x_val + b_val) left_right_fact = -1.0;
            TH1D_radial_projection ->Fill(left_right_fact*distance_to_center,1.0/area);
            //cout << "i_data: " << i_data << ", distance_to_center: " << distance_to_center << endl;
  
        }
    }
    //-------------------------------------------------------


    //-------------------------------------------------------
    if(flag_input_data == 2) // TPolyMarker3D data input
    {
        cout << "data size: " << TPM3D_data->GetN() << endl;
        for(Int_t i_data = 0; i_data < TPM3D_data->GetN(); i_data++)
        {

            Double_t x_val,y_val,z_val;
            TPM3D_data->GetPoint(i_data,x_val,y_val,z_val);
            Double_t distance_to_center = TMath::Sqrt(TMath::Power(x_val-TV2_center.Px(),2) + TMath::Power(y_val-TV2_center.Py(),2));
            Double_t bin_width = TH1D_radial_projection ->GetBinWidth(TH1D_radial_projection->FindBin(distance_to_center));
            Double_t area = TMath::Pi()*(TMath::Power(distance_to_center+bin_width,2) - TMath::Power(distance_to_center,2));
            Double_t left_right_fact = 1.0;
            if(y_val < slope*x_val + b_val) left_right_fact = -1.0;
            TH1D_radial_projection ->Fill(left_right_fact*distance_to_center,1.0/area);
            //cout << "i_data: " << i_data << ", distance_to_center: " << distance_to_center << endl;
  
        }
    }
    //-------------------------------------------------------
}

void Calc_radial_distribution::set_invert()
{
    for(Int_t i_bin = 1; i_bin <= TH1D_radial_projection->GetNbinsX(); i_bin++)
    {
        Double_t bin_cont = TH1D_radial_projection->GetBinContent(i_bin);
        TH1D_radial_projection->SetBinContent(i_bin,-1.0*bin_cont);
    }
}

TH1D* Calc_radial_distribution::get_radial_distribution()
{
    return TH1D_radial_projection;
}

void Calc_radial_distribution::set_Errors_to_zero()
{
    for(Int_t bin_x = 1; bin_x <= TH1D_radial_projection->GetNbinsX(); bin_x++)
    {
        TH1D_radial_projection->SetBinError(bin_x,0.0);
    }
}

void Calc_radial_distribution::set_Zero_level()
{
    Double_t min_bin_cont = TH1D_radial_projection->GetBinContent(TH1D_radial_projection->GetMinimumBin());
    for(Int_t bin_x = 1; bin_x <= TH1D_radial_projection->GetNbinsX(); bin_x++)
    {
        Double_t bin_cont = TH1D_radial_projection->GetBinContent(bin_x);
        TH1D_radial_projection ->SetBinContent(bin_x,bin_cont-min_bin_cont);
    }
}

void Calc_radial_distribution::clear()
{
    TH1D_radial_projection         ->Delete();
    TH1D_radial_projection_entries ->Delete();
    //TH2D_in                ->Delete();
}
//----------------------------------------------------------------------------------------
#endif



//----------------------------------------------------------------------------------------
TH1D* Reflect_TH1D(TH1D* TH1D_in, Double_t center)
{

    TH1D* TH1D_reflect = (TH1D*)TH1D_in->Clone("TH1D_reflect");
    TH1D_reflect->Reset();

    for(Int_t bin_x = 1; bin_x <= TH1D_in->GetNbinsX(); bin_x++)
    {
        Double_t bin_cont = TH1D_in->GetBinContent(bin_x);
        Double_t bin_cent = TH1D_in->GetBinCenter(bin_x);

        Double_t Delta_x   = center - bin_cent;
        Double_t new_x_val = center + Delta_x;
        Int_t    new_bin   = TH1D_reflect->FindBin(new_x_val);
        TH1D_reflect       ->SetBinContent(new_bin,bin_cont);
    }
    return TH1D_reflect;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t Calc_hist_diff(TH1D* TH1D_A, TH1D* TH1D_B, Double_t start, Double_t stop)
{
    // Calculates sum of bin-by-bin difference between two histograms of equal range and binning

    Double_t hist_diff = 0.0;
    for(Int_t bin_x = 1; bin_x <= TH1D_A->GetNbinsX(); bin_x++)
    {
        Double_t bin_cent   = TH1D_A->GetBinCenter(bin_x);
        Double_t bin_cont_A = TH1D_A->GetBinContent(bin_x);
        Double_t bin_cont_B = TH1D_B->GetBinContent(bin_x);

        if(bin_cent > start && bin_cent < stop)
        {
            hist_diff += fabs(bin_cont_A - bin_cont_B);
        }
    }
    return hist_diff;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TPolyMarker* convert_TPM3D_to_TPM(TPolyMarker3D* TPM3D_in)
{
    TPolyMarker* TPM_out = new TPolyMarker();

    for(Int_t i_point = 0; i_point < TPM3D_in->GetN(); i_point++)
    {
        Double_t x,y,z;
        TPM3D_in->GetPoint(i_point,x,y,z);

        TPM_out ->SetNextPoint(x,y);
    }

    return TPM_out;
}
//----------------------------------------------------------------------------------------


#if 1
//----------------------------------------------------------------------------------------
Double_t dist_to_2D_line(Double_t x_line_point, Double_t y_line_point, Double_t alpha,
                        Double_t x_point, Double_t y_point)
{
    TVector2 TV2_line, TV2_line_dir, TV2_point;
    TV2_line.Set(x_line_point,y_line_point);
    TV2_line_dir.Set(TMath::Cos(alpha),TMath::Sin(alpha));
    TV2_line_dir *= 1.0/TV2_line_dir.Mod();
    TV2_point.Set(x_point,y_point);

    TVector2 TV2_diff = TV2_point;
    TV2_diff -= TV2_line;
    Double_t proj = TV2_diff *= TV2_line_dir;
    TVector2 TV2_dir_proj = TV2_line_dir;
    TV2_dir_proj *= proj;
    TVector2 TV2_dist = TV2_diff;
    TV2_dist -= TV2_dir_proj;

    Double_t dist = TV2_dist.Mod();

    return dist;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
static Double_t x_values_minimize[7], y_values_minimize[7];
void fcn(int& nDim, Double_t* gout, Double_t& result, Double_t par[], int flg)
{
    result = 0.0;
    for(Int_t i_point = 0; i_point < 6; i_point++)
    {
        result += TMath::Power(dist_to_2D_line(par[0], par[1], par[2], x_values_minimize[i_point], y_values_minimize[i_point]),2);
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fit_line_to_points()
{
    TFitter* minimizer = new TFitter(3);
    minimizer->SetFCN(fcn);
    // arg1 – parameter numbe
    // arg2 – parameter name
    // arg3 – first guess at parameter value
    // arg4 – estimated distance to minimum
    // arg5, arg6 – ignore for now
    minimizer->SetParameter(0,"x line",2,1,0,0);
    minimizer->SetParameter(1,"y line",2,1,0,0);
    minimizer->SetParameter(2,"alpha dir",2,1,0,0);

    minimizer->ExecuteCommand("SIMPLEX",0,0);
    Double_t bestX     = minimizer->GetParameter(0);
    Double_t bestY     = minimizer->GetParameter(1);
    Double_t bestalpha = minimizer->GetParameter(2);

    cout << "bestX: " << bestX << ", bestY: " << bestY << ", bestalpha: " << bestalpha << endl;

    TPolyLine* TPL_line_fit = new TPolyLine();
    Double_t y_plot_start = 0.0;
    Double_t bestX_dir = TMath::Cos(bestalpha);
    Double_t bestY_dir = TMath::Sin(bestalpha);

    TPL_line_fit ->SetNextPoint(bestX-1200.0*bestX_dir,bestY-1200.0*bestY_dir);
    TPL_line_fit ->SetNextPoint(bestX+1200.0*bestX_dir,bestY+1200.0*bestY_dir);
    TPL_line_fit ->SetLineStyle(9);
    TPL_line_fit ->SetLineColorAlpha(kBlack,0.5);
    TPL_line_fit ->SetLineWidth(2);
    TPL_line_fit ->Draw("ogl");

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t dist_to_circle(Double_t x_circle_center, Double_t y_circle_center, Double_t radius,
                        Double_t x_point, Double_t y_point)
{
    Double_t dist_to_center = TMath::Sqrt(TMath::Power(x_point-x_circle_center,2) + TMath::Power(y_point-y_circle_center,2));
    Double_t dist = dist_to_center - radius;

    return dist;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fcn_circle(int& nDim, Double_t* gout, Double_t& result, Double_t par[], int flg)
{
    result = 0.0;
    for(Int_t i_point = 0; i_point < 7; i_point++)
    {
        result += TMath::Power(dist_to_circle(par[0], par[1], par[2], x_values_minimize[i_point], y_values_minimize[i_point]),2);
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void fit_circle_to_points()
{
    TFitter* minimizer = new TFitter(3);
    minimizer->SetFCN(fcn_circle);
    // arg1 – parameter numbe
    // arg2 – parameter name
    // arg3 – first guess at parameter value
    // arg4 – estimated distance to minimum
    // arg5, arg6 – ignore for no
    minimizer->SetParameter(0,"x circle center",2,1,0,0);
    minimizer->SetParameter(1,"y circle center",2,1,0,0);
    minimizer->SetParameter(2,"radius",2,1,0,0);

    minimizer->ExecuteCommand("SIMPLEX",0,0);
    Double_t bestX     = minimizer->GetParameter(0);
    Double_t bestY     = minimizer->GetParameter(1);
    Double_t radius    = minimizer->GetParameter(2);

    cout << "bestX: " << bestX << ", bestY: " << bestY << ", bestradius: " << radius << endl;

    Draw_Circle_2D_new(radius,radius,1,1,kBlack,9,3,bestX,bestY,0,0.25,0.0);

    /*
    TPolyLine* TPL_line_fit = new TPolyLine();
    Double_t y_plot_start = 0.0;
    Double_t bestX_dir = TMath::Cos(bestalpha);
    Double_t bestY_dir = TMath::Sin(bestalpha);

    TPL_line_fit ->SetNextPoint(bestX-1200.0*bestX_dir,bestY-1200.0*bestY_dir);
    TPL_line_fit ->SetNextPoint(bestX+1200.0*bestX_dir,bestY+1200.0*bestY_dir);
    TPL_line_fit ->SetLineStyle(9);
    TPL_line_fit ->SetLineColorAlpha(kBlack,0.5);
    TPL_line_fit ->SetLineWidth(2);
    TPL_line_fit ->Draw("ogl");
    */
}
//----------------------------------------------------------------------------------------
#endif

#if 1
//----------------------------------------------------------------------------------------
Double_t dist_to_3D_circle_earth(Double_t radius_to_center, Double_t longitude, Double_t latitude,
                                 Double_t point_longitude, Double_t point_latitude)
{
    // longitude = "x"
    // latitude  = "y"  Noth-South

    // Circle center transformation to cartesian
    Double_t x_c = radius_to_center * TMath::Cos(latitude) * TMath::Cos(longitude);
    Double_t y_c = radius_to_center * TMath::Cos(latitude) * TMath::Sin(longitude);
    Double_t z_c = radius_to_center * TMath::Sin(latitude);
    TVector3 tv3_circle, tv3_circle_norm;
    tv3_circle.SetXYZ(x_c,y_c,z_c);
    tv3_circle_norm = tv3_circle;
    tv3_circle_norm *= 1.0/tv3_circle_norm.Mag(); // normalized

    Double_t Radius_earth_AA = 6371.0;
    if(radius_to_center < 0.0 || radius_to_center > Radius_earth_AA) return -999999999.0;
    Double_t circle_radius = TMath::Sqrt(Radius_earth_AA*Radius_earth_AA - radius_to_center*radius_to_center);


    // Point transformation to cartesian
    point_latitude  *= TMath::DegToRad();
    point_longitude *= TMath::DegToRad();
    Double_t x_p = Radius_earth_AA * TMath::Cos(point_latitude) * TMath::Cos(point_longitude);
    Double_t y_p = Radius_earth_AA * TMath::Cos(point_latitude) * TMath::Sin(point_longitude);
    Double_t z_p = Radius_earth_AA * TMath::Sin(point_latitude);
    TVector3 tv3_point;
    tv3_point.SetXYZ(x_p,y_p,z_p);

    //cout << "point = {" << x_p << ", " << y_p << ", " << z_p << "}" << endl;
    //tv3_circle_norm.Print();

    Double_t proj_length = tv3_circle_norm.Dot(tv3_point);
    TVector3 tv3_proj = tv3_circle_norm;
    tv3_proj *= proj_length;
    TVector3 tv3_circle_radius = tv3_circle_norm;
    tv3_circle_radius *= radius_to_center;
    TVector3 tv3_diff_length = tv3_circle_radius;
    tv3_diff_length -= tv3_proj;
    Double_t diff_axis_length = tv3_diff_length.Mag();

    Double_t point_vec_length = tv3_point.Mag();
    Double_t radial_distance = TMath::Sqrt(point_vec_length*point_vec_length - proj_length*proj_length);

    Double_t dist = circle_radius - radial_distance;
    dist = TMath::Sqrt(dist*dist + diff_axis_length*diff_axis_length);
    //dist = diff_axis_length;

    //cout << "proj_length: " << proj_length << ", diff_axis_length: " << diff_axis_length << endl;

    return dist;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
static Double_t x_values_WGS_minimize[7], y_values_WGS_minimize[7];
void fcn_3D_circle_earth(int& nDim, Double_t* gout, Double_t& result, Double_t par[], int flg)
{
    result = 0.0;
    for(Int_t i_point = 0; i_point < 7; i_point++)
    {
        result += TMath::Power(dist_to_3D_circle_earth(par[0], par[1], par[2], x_values_WGS_minimize[i_point], y_values_WGS_minimize[i_point]),2);
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
static Double_t best_radius_3D, best_latitude_3D, best_longitude_3D;
void fit_3D_circle_earth_to_points()
{
    TFitter* minimizer = new TFitter(3);
    minimizer->SetFCN(fcn_3D_circle_earth);
    // arg1 – parameter numbe
    // arg2 – parameter name
    // arg3 – first guess at parameter value
    // arg4 – estimated distance to minimum
    // arg5, arg6 – ignore for now
    minimizer->SetParameter(0,"radius",3000.0,1,0,0);
    minimizer->SetParameter(1,"circle longitude",0.0,1,0,0);
    minimizer->SetParameter(2,"circle latitude",0.0,1,0,0);


    //minimizer->ExecuteCommand("SIMPLEX",0,0);
    //minimizer->ExecuteCommand("migrad",0,0);
    //minimizer->ExecuteCommand("minos",0,0);
    minimizer->ExecuteCommand("seek",0,0);

    best_radius_3D       = minimizer->GetParameter(0);
    best_longitude_3D    = minimizer->GetParameter(1);
    best_latitude_3D     = minimizer->GetParameter(2);

    cout << "best radius_3D: " << best_radius_3D << ", best_latitude_3D: " << best_latitude_3D << ", best_longitude_3D: " << best_longitude_3D << endl;

    //Draw_Circle_2D_new(radius,radius,1,1,kBlack,9,3,bestX,bestY,0,0.25,0.0);

    /*
    TPolyLine* TPL_line_fit = new TPolyLine();
    Double_t y_plot_start = 0.0;
    Double_t bestX_dir = TMath::Cos(bestalpha);
    Double_t bestY_dir = TMath::Sin(bestalpha);

    TPL_line_fit ->SetNextPoint(bestX-1200.0*bestX_dir,bestY-1200.0*bestY_dir);
    TPL_line_fit ->SetNextPoint(bestX+1200.0*bestX_dir,bestY+1200.0*bestY_dir);
    TPL_line_fit ->SetLineStyle(9);
    TPL_line_fit ->SetLineColorAlpha(kBlack,0.5);
    TPL_line_fit ->SetLineWidth(2);
    TPL_line_fit ->Draw("ogl");
    */
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_3D_circle_planet(Double_t radius_to_center, Double_t longitude, Double_t latitude,
                           Double_t line_width, Int_t line_style, Int_t line_color)
{
    cout << " " << endl;
    cout << "Draw_3D_circle_planet" << endl;
    Double_t Radius_earth_AA = 6371.0+2.0;
    Double_t circle_radius = TMath::Sqrt(Radius_earth_AA*Radius_earth_AA - radius_to_center*radius_to_center);

    // longitude = "x"
    // latitude  = "y"  Noth-South

    Double_t x_c = radius_to_center * TMath::Cos(latitude) * TMath::Cos(longitude);
    Double_t y_c = radius_to_center * TMath::Cos(latitude) * TMath::Sin(longitude);
    Double_t z_c = radius_to_center * TMath::Sin(latitude);

    TPolyLine3D* tpl3D_center_line = new TPolyLine3D();
    tpl3D_center_line->SetNextPoint(0.0,0.0,0.0);
    tpl3D_center_line->SetNextPoint(x_c,y_c,z_c);
    tpl3D_center_line->SetLineStyle(1);
    tpl3D_center_line->SetLineColor(kGreen); // 28
    tpl3D_center_line->SetLineWidth(4);
    //tpl3D_center_line->DrawClone("ogl");

    TVector3 tv3_circle, tv3_circle_norm;
    tv3_circle.SetXYZ(x_c,y_c,z_c);
    //tv3_circle.Print();
    tv3_circle_norm = tv3_circle;
    tv3_circle_norm *= 1.0/tv3_circle_norm.Mag();

    TVector3 tv3_circle_x = tv3_circle_norm.Orthogonal();
    TVector3 tv3_circle_y = tv3_circle_x.Cross(tv3_circle_norm);

    tv3_circle_x *= 1.0/tv3_circle_x.Mag();
    tv3_circle_y *= 1.0/tv3_circle_y.Mag();

    //cout << "radius_to_center: " << radius_to_center << ", latitude: " << latitude << ", longitude: " << longitude << endl;
    //cout << "circle_vec = {" << tv3_circle.X() << ", " << tv3_circle.Y() << ", " << tv3_circle.Z() << "}" << endl;
    //cout << "circle_vec_norm = {" << tv3_circle_norm.X() << ", " << tv3_circle_norm.Y() << ", " << tv3_circle_norm.Z() << "}" << endl;
    //cout << "x_vec = {" << tv3_circle_x.X() << ", " << tv3_circle_x.Y() << ", " << tv3_circle_x.Z() << "}" << endl;

    TPolyLine3D* tpl3D_circle = new TPolyLine3D();
    const Int_t N_points = 100;
    TVector3 tv3_circle_point;
    for(Int_t i_point = 0; i_point < N_points+1; i_point++)
    {
        Double_t alpha = 2.0*TMath::Pi()/((Double_t)N_points);
        Double_t x_val = TMath::Cos(alpha*(Double_t)i_point);
        Double_t y_val = TMath::Sin(alpha*(Double_t)i_point);

        TVector3 x_vec = tv3_circle_x;
        x_vec *= x_val;
        TVector3 y_vec = tv3_circle_y;
        y_vec *= y_val;

        tv3_circle_point = x_vec;
        tv3_circle_point += y_vec;
        tv3_circle_point *= circle_radius;
        tv3_circle_point += tv3_circle;

        tpl3D_circle ->SetNextPoint(tv3_circle_point.X(),tv3_circle_point.Y(),tv3_circle_point.Z());
        Double_t r_val = TMath::Sqrt(tv3_circle_point.X()*tv3_circle_point.X() + tv3_circle_point.Y()*tv3_circle_point.Y() + tv3_circle_point.Z()*tv3_circle_point.Z());
        //cout << "i_point: " << i_point << ", r_val: " << r_val << ", alpha: " << alpha << ", x_val: " << x_val << ", {" << tv3_circle_point.X() << ", " << tv3_circle_point.Y()
        //    << ", " << tv3_circle_point.Z() << "}" << endl;
    }

    tpl3D_circle->SetLineStyle(line_style);
    tpl3D_circle->SetLineColor(line_color); // 28
    tpl3D_circle->SetLineWidth(line_width);
    tpl3D_circle->DrawClone("ogl");

#if 0
    //----------------------------------
    // Point transformation to cartesian
    TPolyLine3D* tpl3D_point_line[7];

    Int_t color_table[7] = {kRed,kBlue,kGreen,kYellow,kOrange,kMagenta,kBlack};
    for(Int_t i_crater = 0; i_crater < 7; i_crater++)
    {
        Double_t point_longitude = x_values_WGS_minimize[i_crater];
        Double_t point_latitude  = y_values_WGS_minimize[i_crater];
        point_latitude  *= TMath::DegToRad();
        point_longitude *= TMath::DegToRad();
        Double_t x_p = Radius_earth_AA * TMath::Cos(point_latitude) * TMath::Cos(point_longitude);
        Double_t y_p = Radius_earth_AA * TMath::Cos(point_latitude) * TMath::Sin(point_longitude);
        Double_t z_p = Radius_earth_AA * TMath::Sin(point_latitude);
        TVector3 tv3_point;
        tv3_point.SetXYZ(x_p,y_p,z_p);

        tpl3D_point_line[i_crater] = new TPolyLine3D();
        tpl3D_point_line[i_crater]->SetNextPoint(0.0,0.0,0.0);
        tpl3D_point_line[i_crater]->SetNextPoint(x_p,y_p,z_p);
        tpl3D_point_line[i_crater]->SetLineStyle(line_style);
        tpl3D_point_line[i_crater]->SetLineColor(line_color); // 28
        tpl3D_point_line[i_crater]->SetLineWidth(line_width);
        tpl3D_point_line[i_crater]->DrawClone("ogl");
    }
    //----------------------------------
#endif
}
//----------------------------------------------------------------------------------------



#endif



//----------------------------------------------------------------------------------------
class Class_peak_finder
{
    TH1D* TH1D_in;
    Double_t search_width;
    Double_t search_min_diff; // minimum difference between signal and background
    Double_t exclude_width; // range to be excluded between two peaks
    Double_t bin_width;
    Int_t    N_maxima;
    Double_t min_radius;
    Double_t max_width_scan;
    Int_t    flag_peak_rim; // 0 = peak finding, 1 = rim finding
    Int_t    flag_debug = 0; // 1 = return debug output
    TPolyMarker* pm_peaks;
    std::vector< std::vector<Double_t> > vec_peak_positions;
public:
    void  add_TH1D(TH1D* TH1D_in_a);
    void  set_finding_parameters(Int_t flag_peak_rim_in, Double_t search_width_in,
                                 Double_t search_min_diff_in, Double_t exclude_width_in,
                                 Int_t N_maxima_in, Double_t min_radius_in, Double_t max_width_scan_in);
    void  find_peaks();
    TPolyMarker* get_polymarker();
    std::vector< std::vector<Double_t> > get_peak_positions();
    void set_debug(Int_t flag_debug_in);
    void  clear();
};

void Class_peak_finder::set_finding_parameters(Int_t flag_peak_rim_in, Double_t search_width_in,
                                               Double_t search_min_diff_in, Double_t exclude_width_in,
                                               Int_t N_maxima_in, Double_t min_radius_in, Double_t max_width_scan_in)
{
    flag_peak_rim   = flag_peak_rim_in;
    search_width    = search_width_in;
    search_min_diff = search_min_diff_in;
    exclude_width   = exclude_width_in;
    N_maxima        = N_maxima_in;
    min_radius      = min_radius_in;
    max_width_scan  = max_width_scan_in;
}

void Class_peak_finder::add_TH1D(TH1D* TH1D_in_a)
{
    TH1D_in = (TH1D*)TH1D_in_a->Clone("TH1D_in");
    bin_width = TH1D_in->GetBinWidth(1);
    pm_peaks = new TPolyMarker();
}

void Class_peak_finder::find_peaks()
{
    Int_t    max_bin   = TH1D_in->GetNbinsX();
    Double_t min_val_x = TH1D_in->GetBinCenter(1);
    Double_t bin_width = TH1D_in->GetBinWidth(1);
    Double_t max_val_x = TH1D_in->GetBinCenter(max_bin);

    std::vector< std::vector<Double_t> > vec_diff;
    vec_diff.resize(8);

    vec_peak_positions.resize(4); // radial distance, height, density, width

    // Loop over all bins
    for(Int_t i_bin = 1; i_bin <= max_bin; i_bin++)
    {
        Double_t center_val    = TH1D_in->GetBinCenter(i_bin);
        Double_t signal_height = TH1D_in->GetBinContent(i_bin);
        Double_t hist_start    = TH1D_in->GetBinCenter(1);

        if(center_val < min_radius) continue;

        // Determine search width with maximum density = number of points / area
        Double_t max_signal         = 0.0;
        Double_t max_width          = 0.0;
        Double_t max_signal_density = 0.0;
        Double_t max_veto           = 1.0;
        Double_t max_bkgr_left      = 0.0;
        Double_t max_bkgr_right     = 0.0;
        for(Int_t i_search = 0; i_search < 20; i_search++)
        {
            Double_t search_width_scan          = search_width + search_width*i_search;
            Int_t    search_width_scan_bin_half = TH1D_in->FindBin(hist_start + search_width_scan/2.0) - 1;
            if(search_width_scan_bin_half == 0) continue;
            search_width_scan = 2.0*search_width_scan_bin_half*bin_width + bin_width; // real bin width, factor 2 for left and right, last bin_width for center bin


            if(search_width_scan > max_width_scan)
            {
                //printf("i_bin: %d, i_search: %d, search_width_scan: %f, max_width_scan: %f \n",i_bin,i_search,search_width_scan,max_width_scan);
                break;
            }

            Int_t bin_left_bkgr    = i_bin - 2*search_width_scan_bin_half - 1;
            Int_t bin_left_signal  = i_bin - search_width_scan_bin_half;
            Int_t bin_right_signal = i_bin + search_width_scan_bin_half;
            Int_t bin_right_bkgr   = i_bin + 2*search_width_scan_bin_half + 1;

            Int_t N_bins_signal = fabs(bin_right_signal - bin_left_signal) + 1;
            Int_t N_bins_bkgr   = fabs(bin_left_signal - bin_left_bkgr);
            N_bins_bkgr += fabs(bin_right_bkgr - bin_right_signal);
            Double_t bkgr_scaling_factor = 1.0;
            if(N_bins_bkgr > 0.0) bkgr_scaling_factor = ((Double_t)N_bins_signal)/((Double_t)N_bins_bkgr);

            if(bin_left_bkgr < 1) break;
            if(bin_right_bkgr > max_bin) break;

            // Find minimum value in search range
            Double_t min_val_scan_range = 0.0;
            for(Int_t i_bin_search_min = bin_left_bkgr; i_bin_search_min <= bin_right_bkgr; i_bin_search_min++)
            {
                Double_t bin_cont = TH1D_in->GetBinContent(i_bin_search_min);
                if(bin_cont < min_val_scan_range) min_val_scan_range = bin_cont;
            }

            Double_t sum_signal     = 0.0;
            Double_t sum_bkgr       = 0.0;
            Double_t sum_bkgr_left  = 0.0;
            Double_t sum_bkgr_right = 0.0;
            for(Int_t i_bin_search_min = bin_left_bkgr; i_bin_search_min <= bin_right_bkgr; i_bin_search_min++)
            {
                Double_t bin_cont   = TH1D_in->GetBinContent(i_bin_search_min) - min_val_scan_range;
                if(i_bin_search_min >= bin_left_signal && i_bin_search_min <= bin_right_signal) sum_signal     += bin_cont;
                if(i_bin_search_min >= bin_left_bkgr && i_bin_search_min < bin_left_signal)     sum_bkgr_left  += bin_cont;
                if(i_bin_search_min > bin_right_signal && i_bin_search_min <= bin_right_bkgr)   sum_bkgr_right += bin_cont;
            }
            sum_bkgr_left  *= bkgr_scaling_factor;
            sum_bkgr_right *= bkgr_scaling_factor;
            sum_bkgr = sum_bkgr_left + sum_bkgr_right;


            if(flag_debug) printf("i_bin: %d, i_search: %d, search_width: %f, bkgr_scaling_factor: %f, bin_width: %d, bins: {%d,%d,%d,%d} \n",i_bin,i_search,search_width_scan,bkgr_scaling_factor,search_width_scan_bin_half,bin_left_bkgr,bin_left_signal,bin_right_signal,bin_right_bkgr);


            if(sum_bkgr > sum_signal) continue;
            if(sum_bkgr_left > sum_signal/2.0) continue;  // background has only half the width of signal -> factor 1/2.0
            if(sum_bkgr_right > sum_signal/2.0) continue;
            if(sum_bkgr == 0.0) continue;

            Double_t signal = sum_signal - sum_bkgr;
            Double_t signal_density = signal/search_width_scan;

            //printf("i_bin: %d, i_search: %d, center_val: %f, signal: %f, signal_density: %f,  width: %f \n",i_bin,i_search,center_val,signal,signal_density,search_width_scan);
     

            if(signal_density > max_signal_density)
            {
                max_signal         = signal;
                max_signal_density = signal_density;
                max_width          = search_width_scan;
                max_veto           = 0.0;
                max_bkgr_left      = sum_bkgr_left;
                max_bkgr_right     = sum_bkgr_right;
                //printf("i_bin: %d, i_search: %d, center_val: %f, max_signal: %f, max_signal_density: %f,  max_width: %f \n",i_bin,i_search,center_val,max_signal,max_signal_density,max_width);
            }
        }

        if(max_width <= 0.0) continue;

        vec_diff[0].push_back(center_val);
        vec_diff[1].push_back(signal_height);
        vec_diff[2].push_back(max_signal_density);
        vec_diff[3].push_back(max_signal);
        vec_diff[4].push_back(max_veto);
        vec_diff[5].push_back(max_bkgr_left);
        vec_diff[6].push_back(max_bkgr_right);
        vec_diff[7].push_back(max_width);

        if(flag_debug) printf("i_bin: %d, center_val: %f, signal height: %f, density: %f, signal: %f, bkgr_left: %f, bkgr_right: %f, width: %f \n",i_bin,center_val,signal_height,max_signal_density,max_signal,max_bkgr_left,max_bkgr_right,max_width);
    } // end of histogram bin loop

    // Find all maxima, exclude ranges within exclude_width after each loop
    Int_t max_i_diff = 0;
    while(max_i_diff >= 0 && vec_diff[0].size() > 0)
    {
        Double_t best_center_val     = 0.0;
        Double_t best_signal_height  = 0.0;
        Double_t best_signal_density = 0.0;
        Double_t best_signal         = 0.0;
        Double_t best_bkgr_left      = 0.0;
        Double_t best_bkgr_right     = 0.0;
        Double_t best_width          = 0.0;

        max_i_diff = -1;
        if(flag_peak_rim == 0) // peak finding
        {
            //cout << " " << endl;
            //cout << "-----------------------" << endl;
            // Find maximum signal to background above minimum defined
            for(Int_t i_diff = 0; i_diff < vec_diff[0].size(); i_diff++)
            {
                //printf("A density: %f, pos: %f, flag: %f, vec_diff size: %d \n",vec_diff[2][i_diff],vec_diff[0][i_diff],vec_diff[4][i_diff],(Int_t)vec_diff[0].size());
                //printf("   i_diff: %d, S/B: %f, veto: %f \n",i_diff,vec_diff[1][i_diff],vec_diff[3][i_diff]);
                if(vec_diff[2][i_diff] > best_signal_density
                   //&& vec_diff[1][i_diff] >= search_min_diff
                   && vec_diff[4][i_diff] < 1.0)
                {
                    // Check if this point with its width is already in range of a previously found point
                    // Avoid peaks on top of peaks
                    Int_t flag_in_range_of_accepted_point = 0;
                    for(Int_t i_already_accepted_point = 0; i_already_accepted_point < vec_peak_positions[0].size(); i_already_accepted_point++)
                    {
                        Double_t acc_center        = vec_peak_positions[0][i_already_accepted_point];
                        Double_t acc_width         = vec_peak_positions[3][i_already_accepted_point];
                        Double_t this_point_center = vec_diff[0][i_diff];
                        Double_t this_point_width  = vec_diff[7][i_diff];

                        // Make sure that the current peak doesn't share bins with already found peaks
                        if(fabs(this_point_center - acc_center) < (0.5*this_point_width + 0.5*acc_width + bin_width/2.0)) // sharing of half a bin is OK
                        {
                            flag_in_range_of_accepted_point = 1;
                            vec_diff[4][i_diff] = 1.0; // veto this point
                        }
                    }

                    if(!flag_in_range_of_accepted_point)
                    {
                        best_center_val     = vec_diff[0][i_diff];
                        best_signal_height  = vec_diff[1][i_diff];
                        best_signal_density = vec_diff[2][i_diff];
                        best_signal         = vec_diff[3][i_diff];
                        Double_t flag_range = vec_diff[4][i_diff];
                        best_bkgr_left      = vec_diff[5][i_diff];
                        best_bkgr_right     = vec_diff[6][i_diff];
                        best_width          = vec_diff[7][i_diff];
                        max_i_diff          = i_diff;
                        if(flag_debug) printf("    -> peak found at x_pos: %f, density: %f, flag_range: %f \n",best_center_val,best_signal_density,flag_range);
                    }
                }
            }
        }
        else // rim finding
        {
            /*
            // Find maximum signal to background above minimum defined
            for(Int_t i_diff = 0; i_diff < vec_diff[0].size(); i_diff++)
            {
                Double_t slope_left  = vec_diff[2][i_diff] - vec_diff[4][i_diff];
                Double_t slope_right = vec_diff[5][i_diff] - vec_diff[2][i_diff];

                Double_t slope_change = fabs(slope_right - slope_left);
                if(1)
                {
                    //printf("max_i_diff: %d, heights: {%f,%f,%f}, slopes: {%f,%f}, slope_change: %f \n",max_i_diff,vec_diff[4][i_diff],vec_diff[2][i_diff],vec_diff[5][i_diff],slope_left,slope_right,slope_change);
                    if(slope_change > max_diff &&  slope_change > search_min_diff)
                    {
                        max_diff_x = vec_diff[0][i_diff];
                        max_diff   = slope_change;
                        max_height = vec_diff[2][i_diff];
                        max_i_diff = i_diff;
                    }
                }
            }
            */
        }

        if(best_width > exclude_width) exclude_width = best_width;

        //cout << "start erase: "<<  i_diff_left_erase << ", stop erase: " << i_diff_right_erase << endl;

        if(max_i_diff >= 0)
        {
            //printf("polymarker set at: %f, height: %f \n",best_center_val,best_signal_height);
            pm_peaks ->SetNextPoint(best_center_val,best_signal_height);
            vec_peak_positions[0].push_back(best_center_val);
            vec_peak_positions[1].push_back(best_signal_height);
            vec_peak_positions[2].push_back(best_signal_density);
            vec_peak_positions[3].push_back(best_width);

            if(flag_debug) printf("max_i_diff: %d, x,height: {%f,%f}, best_density: %f, best_width: %f, S,Lbkgr,Rbkgr: {%f,%f,%f} \n",max_i_diff,best_center_val,best_signal_height,best_signal_density,best_width,best_signal,best_bkgr_left,best_bkgr_right);

            //cout << "max_diff_x: " << max_diff_x << ", max_diff: " << max_diff << ", max_i_diff: " << max_i_diff << ", size: " << vec_diff[0].size()
            //    << ", peak_total: " << vec_diff[2][max_i_diff] << ", bkgr_left: " << vec_diff[4][max_i_diff] << ", bkgr_right: " << vec_diff[5][max_i_diff] << endl;

        }

        if(vec_peak_positions[0].size() >= N_maxima) break;
    }
}

TPolyMarker* Class_peak_finder::get_polymarker()
{
    return pm_peaks;
}

std::vector< std::vector<Double_t> > Class_peak_finder::get_peak_positions()
{
    return vec_peak_positions;
}

void Class_peak_finder::set_debug(Int_t flag_debug_in)
{
    flag_debug = flag_debug_in;
}

void Class_peak_finder::clear()
{
    delete pm_peaks;
    vec_peak_positions.clear();
    delete TH1D_in;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Class_peak_bkgr_cleaner
{
    std::vector< std::vector<Double_t> > vec_2D_peak_points; // [point number][x,y,height]
    Int_t flag_debug;
    Int_t flag_plot;
    Int_t plot_bin;
    Int_t cat_AB_plot;
    Int_t x_bin_plot;
    Int_t y_bin_plot;
    Int_t cluster_plot;
    Int_t itt_cluster_separation_plot;
    Int_t vec_delta_bin_x[8];
    Int_t vec_delta_bin_y[8];
    Double_t corr_length; // distance between points to be merged for clusters
    Double_t min_walk_to_bird_ratio;
    Int_t N_itt_cluster_separation; // how many times the cluster separation loop is running itteratively
    TH2D* h2D_2D_peak_points = NULL;
    TH2D* h2D_2D_peak_points_clusters = NULL;
    TH2D* h2D_2D_bkgr_clusters = NULL;
    TH2D* h2D_2D_bkgr_clusters_itt = NULL;
    TH2D* h2D_cluster = NULL;
    TH2D* h2D_sub_cluster_map = NULL;
    TH2D* h2D_sub_cluster_map_out = NULL;
    Int_t cluster_index;
    Int_t min_cls_size;
    Int_t min_cluster_size_to_split;
    Int_t max_bin_x_init_full_scan;
    Int_t max_bin_y_init_full_scan;
    Int_t max_bin_x_target_full_scan;
    Int_t max_bin_y_target_full_scan;
    Int_t max_i_bin_total_init_full_scan;
    Int_t max_i_bin_total_target_full_scan;
    Double_t walk_to_bird_power;
    Int_t N_bins_total; // total number of cluster bins n_bins_x*n_bins_y
    Int_t N_new_split_cluster;
    Int_t max_cls_index;
    vector< vector< vector< vector< vector<Int_t> > > > > vec_walking_bins; // stored the full walking path per cluster
    vector< vector< vector< vector< vector<Int_t> > > > > vec_walking_bins_plot; // stored the full walking path per cluster
    vector< vector< vector<Int_t> > >  vec_init_target_xy_used; // stores the init and target bins per cluster
    vector< vector< vector<Int_t> > >  vec_min_max_x_y_cluster;
    Double_t max_ratio_full_scan;
    Double_t min_max_dist;
    vector<Int_t> vec_cluster_size;
    vector<Int_t> vec_bins;
    vector<Int_t> vec_bins_mod;
    Bool_t flag_cluster_plot_available;
public:
    void  set_debug(Int_t flag_debug_in);
    void  set_plot(Int_t flag_plot_in, Int_t plot_bin_in, Int_t cluster_plot_in,
                   Int_t cat_AB_plot_in, Int_t x_bin_plot_in, Int_t y_bin_plot_in,
                   Int_t itt_cluster_separation_plot_in);
    void  add_data(vector< vector<Double_t> > vec_2D_peak_points_in);
    void  add_data_TH2D(TH2D* h2D_2D_peak_points_in);
    void  set_parameters(Double_t corr_length_in, Int_t min_cls_size_in,
                         Double_t min_max_dist_in, Int_t min_cluster_size_to_split_in,
                         Double_t min_walk_to_bird_ratio_in,
                         Int_t N_itt_cluster_separation_in);
    void  calc_vec_bins(Int_t init_bin_x, Int_t init_bin_y, Int_t target_bin_x, Int_t target_bin_y, Int_t flip);
    void  fill_TH2();
    void  select_final_bkgr_clusters();
    TH2D* get_TH2D();
    TH2D* get_bkgr_TH2D();
    void  set_global_cluster_histogram(TH2D* h_global_cluster_histogram_in);
    void  identify_clusters();
    void  select_bkgr_clusters();
    Double_t calc_full_walk_way(vector< vector<Int_t> > vec_walk_path);
    void  do_full_walk(TH2D* h2D_map,Int_t i_cls_in);
    Int_t categorize_sub_clusters(TH2D* h2D_map,Int_t i_cls_in);
    TH2D* get_cluster_plot_map() const {return h2D_sub_cluster_map_out;}
    void  create_heat_maps();
    Bool_t IsClusterPlot() const {return flag_cluster_plot_available;}
    vector< vector< vector< vector< vector<Int_t> > > > > get_vec_walking_bins_plot();
    vector< vector<Double_t> > get_bkgr_cleaned_points();
    vector< vector<Double_t> > get_bkgr_points();
    void  clear();
};


vector< vector<Double_t> > Class_peak_bkgr_cleaner::get_bkgr_cleaned_points()
{
    vector< vector<Double_t> > vec_2D_peak_points_out;
    vector<Double_t> vec_peak_pos;
    vec_peak_pos.resize(5);
    for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
    {
        Double_t x_val   = vec_2D_peak_points[i_point][0];
        Double_t y_val   = vec_2D_peak_points[i_point][1];
        Double_t height  = vec_2D_peak_points[i_point][2];
        Double_t density = vec_2D_peak_points[i_point][3];
        Double_t width   = vec_2D_peak_points[i_point][4];

        Int_t bin = h2D_2D_bkgr_clusters_itt ->FindBin(x_val,y_val);
        if(h2D_2D_bkgr_clusters_itt ->GetBinContent(bin) <= -1.0)
        {
            Int_t point_bkgr = vec_2D_peak_points_out.size();
            //printf("signal point found, point: %d, x/y: {%f,%f} \n",point_bkgr,x_val,y_val);
            for(Int_t i_par = 0; i_par < 5; i_par++)
            {
                vec_peak_pos[i_par] = vec_2D_peak_points[i_point][i_par];
            }
            vec_2D_peak_points_out.push_back(vec_peak_pos);
        }
    }

    return vec_2D_peak_points_out;
}


vector< vector<Double_t> > Class_peak_bkgr_cleaner::get_bkgr_points()
{
    vector< vector<Double_t> > vec_2D_peak_points_out;
    vector<Double_t> vec_peak_pos;
    vec_peak_pos.resize(5);
    for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
    {
        Double_t x_val   = vec_2D_peak_points[i_point][0];
        Double_t y_val   = vec_2D_peak_points[i_point][1];
        Double_t height  = vec_2D_peak_points[i_point][2];
        Double_t density = vec_2D_peak_points[i_point][3];
        Double_t width   = vec_2D_peak_points[i_point][4];

        Int_t bin = h2D_2D_bkgr_clusters_itt ->FindBin(x_val,y_val);
        if(h2D_2D_bkgr_clusters_itt ->GetBinContent(bin) > -1.0)
        {
            Int_t point_bkgr = vec_2D_peak_points_out.size();
            //printf("background point found, point: %d, x/y: {%f,%f} \n",point_bkgr,x_val,y_val);
            for(Int_t i_par = 0; i_par < 5; i_par++)
            {
                vec_peak_pos[i_par] = vec_2D_peak_points[i_point][i_par];
            }
            vec_2D_peak_points_out.push_back(vec_peak_pos);
        }
    }

    return vec_2D_peak_points_out;
}

void Class_peak_bkgr_cleaner::set_debug(Int_t flag_debug_in)
{
    flag_debug = flag_debug_in;
}

void Class_peak_bkgr_cleaner::set_plot(Int_t flag_plot_in, Int_t plot_bin_in, Int_t cluster_plot_in,
                                       Int_t cat_AB_plot_in, Int_t x_bin_plot_in, Int_t y_bin_plot_in,
                                       Int_t itt_cluster_separation_plot_in)
{
    flag_plot                   = flag_plot_in;
    plot_bin                    = plot_bin_in;
    cluster_plot                = cluster_plot_in;
    cat_AB_plot                 = cat_AB_plot_in;
    x_bin_plot                  = x_bin_plot_in;
    y_bin_plot                  = y_bin_plot_in;
    itt_cluster_separation_plot = itt_cluster_separation_plot_in;
}

void Class_peak_bkgr_cleaner::add_data(vector< vector<Double_t> > vec_2D_peak_points_in)
{
    //vector< vector<Double_t> > vec_2D_peak_points(vec_2D_peak_points_in); // deep copy(?) -> clone all array members
    vec_2D_peak_points = vec_2D_peak_points_in; // [x,y,height,density]
}

void Class_peak_bkgr_cleaner::add_data_TH2D(TH2D* h2D_2D_peak_points_in)
{
    h2D_2D_peak_points = h2D_2D_peak_points_in;
    
    h2D_2D_peak_points ->GetXaxis()->SetRangeUser(h2D_2D_peak_points->GetXaxis()->GetBinCenter(1),h2D_2D_peak_points->GetXaxis()->GetBinCenter(h2D_2D_peak_points ->GetNbinsX()));
    h2D_2D_peak_points ->GetYaxis()->SetRangeUser(h2D_2D_peak_points->GetYaxis()->GetBinCenter(1),h2D_2D_peak_points->GetYaxis()->GetBinCenter(h2D_2D_peak_points ->GetNbinsX()));

    h2D_sub_cluster_map_out     = (TH2D*)h2D_2D_peak_points->Clone("h2D_sub_cluster_map_out");
}

void Class_peak_bkgr_cleaner::set_parameters(Double_t corr_length_in, Int_t min_cls_size_in,
                                             Double_t min_max_dist_in, Int_t min_cluster_size_to_split_in,
                                             Double_t min_walk_to_bird_ratio_in,
                                             Int_t N_itt_cluster_separation_in
                                            )
{
    // For constuctor
    vec_bins.resize(8);
    vec_bins_mod.resize(8);
    flag_cluster_plot_available = kFALSE;
    walk_to_bird_power = 1.0;

    corr_length               = corr_length_in;
    min_cls_size              = min_cls_size_in;
    min_max_dist              = min_max_dist_in;
    min_cluster_size_to_split = min_cluster_size_to_split_in;
    min_walk_to_bird_ratio    = min_walk_to_bird_ratio_in;
    N_itt_cluster_separation  = N_itt_cluster_separation_in;
    if(corr_length <= 0.0)
    {
        cout << "WARNING in Class_peak_bkgr_cleaner::set_parameters! corr_length <= 0.0, set to 1.0" << endl;
        corr_length = 1.0;
    }

    vec_delta_bin_x[0] = 1;
    vec_delta_bin_x[1] = 1;
    vec_delta_bin_x[2] = 0;
    vec_delta_bin_x[3] = -1;
    vec_delta_bin_x[4] = -1;
    vec_delta_bin_x[5] = -1;
    vec_delta_bin_x[6] = 0;
    vec_delta_bin_x[7] = 1;

    vec_delta_bin_y[0] = 0;
    vec_delta_bin_y[1] = 1;
    vec_delta_bin_y[2] = 1;
    vec_delta_bin_y[3] = 1;
    vec_delta_bin_y[4] = 0;
    vec_delta_bin_y[5] = -1;
    vec_delta_bin_y[6] = -1;
    vec_delta_bin_y[7] = -1;
}

void Class_peak_bkgr_cleaner::fill_TH2()
{
    // Determine range of points in x and y directions
    Double_t min_max_xy[2][2] = {0.0};
    for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
    {
        Double_t x_val  = vec_2D_peak_points[i_point][0];
        Double_t y_val  = vec_2D_peak_points[i_point][1];

        if(i_point == 0)
        {
            min_max_xy[0][0] = x_val; // min x
            min_max_xy[1][0] = x_val; // max x
            min_max_xy[0][1] = y_val; // min y
            min_max_xy[1][1] = y_val; // max y
        }
        else
        {
            if(x_val < min_max_xy[0][0]) min_max_xy[0][0] = x_val;
            if(x_val > min_max_xy[1][0]) min_max_xy[1][0] = x_val;
            if(y_val < min_max_xy[0][1]) min_max_xy[0][1] = y_val;
            if(y_val > min_max_xy[1][1]) min_max_xy[1][1] = y_val;
        }
    }

    // Determine number of bins needed ot match corrleation length
    Double_t length_x = fabs(min_max_xy[1][0] - min_max_xy[0][0]);
    Double_t length_y = fabs(min_max_xy[1][1] - min_max_xy[0][1]);
    Int_t    bins_x   = (Int_t)(length_x/(corr_length*2.0));
    Int_t    bins_y   = (Int_t)(length_y/(corr_length*2.0));


    // Fill 2D histogram
    h2D_2D_peak_points          = new TH2D("h2D_2D_peak_points","h2D_2D_peak_points",bins_x,min_max_xy[0][0],min_max_xy[1][0],bins_y,min_max_xy[0][1],min_max_xy[1][1]);

    for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
    {
        Double_t x_val  = vec_2D_peak_points[i_point][0];
        Double_t y_val  = vec_2D_peak_points[i_point][1];
        h2D_2D_peak_points ->Fill(x_val,y_val);
    }
    h2D_sub_cluster_map_out     = (TH2D*)h2D_2D_peak_points->Clone("h2D_sub_cluster_map_out");
}


void Class_peak_bkgr_cleaner::identify_clusters()
{
    h2D_2D_peak_points_clusters = (TH2D*)h2D_2D_peak_points->Clone("h2D_2D_peak_points_clusters");
    h2D_2D_peak_points_clusters ->Reset();
    h2D_2D_bkgr_clusters = (TH2D*)h2D_2D_peak_points->Clone("h2D_2D_bkgr_clusters");
    h2D_2D_bkgr_clusters ->Reset();

    for(Int_t bin_x = 1; bin_x <= h2D_2D_peak_points_clusters->GetNbinsX(); bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= h2D_2D_peak_points_clusters->GetNbinsY(); bin_y++)
        {
            h2D_2D_peak_points_clusters ->SetBinContent(bin_x,bin_y,-1.0);
            h2D_2D_bkgr_clusters        ->SetBinContent(bin_x,bin_y,-1.0);
        }
    }

    cluster_index = 0;

    Int_t bin_x_nb_shift[8] = {-1,-1,-1,0,0,+1,+1,+1};
    Int_t bin_y_nb_shift[8] = {+1,0,-1,+1,-1,+1,0,-1};

    Int_t max_x_bin = h2D_2D_peak_points->GetNbinsX();
    Int_t max_y_bin = h2D_2D_peak_points->GetNbinsY();

    printf("max_x/y_bin: {%d,%d} \n",max_x_bin,max_y_bin);

#if 1
    for(Int_t bin_x = 2; bin_x < max_x_bin; bin_x++) // exclude edge bins
    {
        for(Int_t bin_y = 2; bin_y < max_y_bin; bin_y++) // exclude edge bins
        {
            Double_t bin_cont = h2D_2D_peak_points                 ->GetBinContent(bin_x,bin_y);
            Int_t    cls_indx = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y);
            Int_t    flag_first_bin = 1;
            if(bin_cont > 0.0 && cls_indx == -1)
            {
                h2D_2D_peak_points_clusters ->SetBinContent(bin_x,bin_y,cluster_index); // flag this bin as used

                vector< vector<Int_t> > vec_bin_index_nb; // stores the bin indices of the neighbouring bins with positive bin content
                vec_bin_index_nb.resize(2); // x_bin, y_bin
                Int_t cluster_size = 0;
                while((Int_t)vec_bin_index_nb[0].size() > 0 || flag_first_bin) // start of cluster search
                {
                    Int_t bin_x_use = bin_x;
                    Int_t bin_y_use = bin_y;
                    if(!flag_first_bin)
                    {
                        bin_x_use = vec_bin_index_nb[0][0]; // always take the first element in the vector, will be deleted at the end
                        bin_y_use = vec_bin_index_nb[1][0];
                    }

                    //printf("x,y: {%d,%d}, flag_first_bin: %d, cluster_index: %d, bin_cont: %f \n",bin_x_use,bin_y_use,flag_first_bin,cluster_index,bin_cont);
                    // Check all 8 neighbouring bins for bin content > 0.0
                    for(Int_t i_nb_bin = 0; i_nb_bin < 8; i_nb_bin++)
                    {
                        Int_t bin_x_nb = bin_x_use + bin_x_nb_shift[i_nb_bin];
                        Int_t bin_y_nb = bin_y_use + bin_y_nb_shift[i_nb_bin];

                        if(
                           bin_x_nb    < 1
                           || bin_y_nb < 1
                           || bin_x_nb > max_x_bin+1
                           || bin_y_nb > max_y_bin+1
                          )
                        {
                            continue;
                        }


                        Double_t bin_cont_nb = h2D_2D_peak_points->GetBinContent(bin_x_nb,bin_y_nb);
                        Int_t    cls_indx_nb = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x_nb,bin_y_nb);

                        //printf("-----------------------> bin_cont_nb: %f, cls_indx_nb: %d \n",bin_cont_nb,cls_indx_nb);
                        if(bin_cont_nb > 0.0 && cls_indx_nb == -1)
                        {
                            //printf("--------> i_nb_bin: %d, bin_cont_nb: %f \n",i_nb_bin,bin_cont_nb);
                            h2D_2D_peak_points_clusters ->SetBinContent(bin_x_nb,bin_y_nb,cluster_index);
                            vec_bin_index_nb[0].push_back(bin_x_nb);
                            vec_bin_index_nb[1].push_back(bin_y_nb);
                            cluster_size++;
                        }
                    }

                    if(!flag_first_bin)
                    {
                        vec_bin_index_nb[0].erase(vec_bin_index_nb[0].begin(),vec_bin_index_nb[0].begin()+1);
                        vec_bin_index_nb[1].erase(vec_bin_index_nb[1].begin(),vec_bin_index_nb[1].begin()+1);
                    }

                    flag_first_bin = 0;
                }
                //printf("cluster_index: %d, cluster_size: %d \n",cluster_index,cluster_size);
                vec_cluster_size.push_back(cluster_size);
                cluster_index++;
            }
        }
    }

#endif
}


void Class_peak_bkgr_cleaner::set_global_cluster_histogram(TH2D* h_global_cluster_histogram_in)
{
    printf("Class_peak_bkgr_cleaner::set_global_cluster_histogram \n");
    if(h2D_2D_peak_points_clusters) h2D_2D_peak_points_clusters ->Delete();
    h2D_2D_peak_points_clusters = (TH2D*)h_global_cluster_histogram_in->Clone("h2D_2D_peak_points_clusters");
}



void Class_peak_bkgr_cleaner::select_bkgr_clusters()
{
    cout << "" << endl;
    printf("Class_peak_bkgr_cleaner::select_bkgr_clusters() \n");
    Int_t max_x_bin = h2D_2D_peak_points_clusters->GetNbinsX();
    Int_t max_y_bin = h2D_2D_peak_points_clusters->GetNbinsY();

    if(!h2D_2D_bkgr_clusters) h2D_2D_bkgr_clusters = (TH2D*)h2D_2D_peak_points_clusters->Clone("h2D_2D_bkgr_clusters");
    h2D_2D_bkgr_clusters ->Reset();

    // Determine maximum cluster index on map
    Int_t max_cluster_index = 0;
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Int_t cls_index = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y);
            if(cls_index > max_cluster_index) max_cluster_index = cls_index;
        }
    }
    printf("Number of clusters: %d \n",max_cluster_index+1);


    vec_cluster_size.clear();
    vec_cluster_size.resize(max_cluster_index+1);
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        vec_cluster_size[i_cls] = 0;
    }

    // Determine cluster size
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Int_t cls_index = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y);
            //printf("cls_index: %f, size: %d \n",cls_index,(Int_t)vec_cluster_size.size());
            if(cls_index >= 0) vec_cluster_size[cls_index]++;
        }
    }
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        //printf("i_cls: %d, cluster_size: %d \n",i_cls,vec_cluster_size[i_cls]);
    }


    vector< vector<Double_t> > vec_min_max_dist;
    vec_min_max_dist.resize((Int_t)vec_cluster_size.size());
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        vec_min_max_dist[i_cls].resize(4); // min, max distance to center, average distance, fraction of points outside
        vec_min_max_dist[i_cls][0] = 9999999.0;
        vec_min_max_dist[i_cls][1] = 0.0;
        vec_min_max_dist[i_cls][2] = 0.0;
        vec_min_max_dist[i_cls][3] = 0.0;
    }


    // Determine average, minimium and maximum distance to center for each cluster
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Double_t pos_x = h2D_2D_peak_points_clusters ->GetXaxis()->GetBinCenter(bin_x);
            Double_t pos_y = h2D_2D_peak_points_clusters ->GetYaxis()->GetBinCenter(bin_y);
            Double_t radius = TMath::Sqrt(pos_x*pos_x + pos_y*pos_y);

            Int_t cls_indx = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y); // start from 0 for real clusters
            if(cls_indx < 0) continue;
            Int_t cls_size = vec_cluster_size[cls_indx];

            if(cls_size < min_cls_size) continue; // select only clusters with a minimum amount of points
            if(radius < vec_min_max_dist[cls_indx][0]) vec_min_max_dist[cls_indx][0] = radius;
            if(radius > vec_min_max_dist[cls_indx][1]) vec_min_max_dist[cls_indx][1] = radius;
            vec_min_max_dist[cls_indx][2] += radius;
        }
    }

    // Calcualte average distance of cluster to center
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        Double_t cls_size = (Double_t)vec_cluster_size[i_cls];
        if(cls_size > 0.0) vec_min_max_dist[i_cls][2] /= cls_size;
    }

    // Calculate fraction of points outside allowed distance limit
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Double_t pos_x = h2D_2D_peak_points_clusters ->GetXaxis()->GetBinCenter(bin_x);
            Double_t pos_y = h2D_2D_peak_points_clusters ->GetYaxis()->GetBinCenter(bin_y);
            Double_t radius = TMath::Sqrt(pos_x*pos_x + pos_y*pos_y);

            Int_t cls_indx = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y); // start from 0 for real clusters
            if(cls_indx < 0) continue;
            Int_t cls_size = vec_cluster_size[cls_indx];

            if(cls_size < min_cls_size) continue; // select only clusters with a minimum amount of points
            Double_t average_radius = vec_min_max_dist[cls_indx][2];
            if(fabs(radius - average_radius) > 0.5*min_max_dist)
            {
                vec_min_max_dist[cls_indx][3] += 1.0;
            }
        }
    }

    // Calcualte fraction of points outside allowed limit
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        Double_t cls_size = (Double_t)vec_cluster_size[i_cls];
        if(cls_size > 0.0) vec_min_max_dist[i_cls][3] /= cls_size;
    }


    // Check if minimum to maximum distance exceeds allowed limit for a good (ring like) cluster
    vector<Int_t> vec_flag_bkgr_cluster;
    vec_flag_bkgr_cluster.resize((Int_t)vec_cluster_size.size());
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        Int_t cls_size = vec_cluster_size[i_cls]; // starts from 0 for real clusters
        if(cls_size < min_cls_size) continue; // select only clusters with a minimum amount of points

        Double_t delta_dist = fabs(vec_min_max_dist[i_cls][1] - vec_min_max_dist[i_cls][0]);
        Double_t fraction_outside = vec_min_max_dist[i_cls][3];
        //if(delta_dist > min_max_dist)
        if(fraction_outside > 0.1)
        {
            vec_flag_bkgr_cluster[i_cls] = 1;
        }
        else vec_flag_bkgr_cluster[i_cls] = 0;
        //printf("i_cls: %d , cls_size: %d, delta_dist: %f, flag: %d \n",i_cls,cls_size,delta_dist,vec_flag_bkgr_cluster[i_cls]);
    }


    // Fill background cluster histogram
    for(Int_t bin_x = 1; bin_x <= h2D_2D_bkgr_clusters->GetNbinsX(); bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= h2D_2D_bkgr_clusters->GetNbinsY(); bin_y++)
        {
            h2D_2D_bkgr_clusters        ->SetBinContent(bin_x,bin_y,-1.0);
        }
    }
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Int_t cls_indx = (Int_t)h2D_2D_peak_points_clusters ->GetBinContent(bin_x,bin_y);
            if(cls_indx < 0) continue;
            Int_t flag_bkgr_cluster = vec_flag_bkgr_cluster[cls_indx];
            if(flag_bkgr_cluster)
            {
                h2D_2D_bkgr_clusters ->SetBinContent(bin_x,bin_y,cls_indx);
                //printf("bin x,y: {%d,%d}, cls_indx: %d, flag_bgkr_cluster: %d \n",bin_x,bin_y,cls_indx,flag_bkgr_cluster);
            }

        }
    }


    h2D_2D_bkgr_clusters_itt = (TH2D*)h2D_2D_bkgr_clusters->Clone("h2D_2D_bkgr_clusters_itt");

}

void Class_peak_bkgr_cleaner::calc_vec_bins(Int_t init_bin_x, Int_t init_bin_y, Int_t target_bin_x, Int_t target_bin_y, Int_t flip)
{
    // Calculates the order of the neighbouring bins where to go next for a pathfinding algorithm
    // 3  2  1
    // 4  X  0
    // 5  6  7

    Int_t delta_bin_x = target_bin_x - init_bin_x;
    Int_t delta_bin_y = target_bin_y - init_bin_y;

    // Calculate order in which bins should be searched
    // 8 neighbouring bins, 16 direction orderings -> angle 22.5 degrees
    Double_t angle_dir = TMath::ATan2(delta_bin_y,delta_bin_x)*TMath::RadToDeg(); // [-180.0,180.0]
    if(angle_dir < 0.0) angle_dir += 360.0;
    Double_t delta_angle  = 45.0;

    Double_t angle_dir_mod = angle_dir + 22.5;
    if(angle_dir_mod > 360.0) angle_dir_mod -= 360.0;

    Int_t angle_bin     = (Int_t)(angle_dir_mod/delta_angle);
    Int_t angle_bin_dir = (Int_t)(angle_dir_mod/(0.5*delta_angle));
    Int_t angle_sign = angle_bin_dir%2;
    //printf("angle_dir: %f, angle_dir_mod: %f, angle_bin_dir: %d, angle_sign: %d \n",angle_dir,angle_dir_mod,angle_bin_dir,angle_sign);
    Int_t angle_dir_sign     = -1;
    Int_t angle_dir_sign_mod = -1;
    if(angle_sign == 1) angle_dir_sign     = 1;
    if(angle_sign == 1) angle_dir_sign_mod = 1;
    //printf("angle_bin: %d, angle_dir_sign: %d \n",angle_bin,angle_dir_sign);
    vec_bins[0]     = angle_bin;
    vec_bins_mod[0] = angle_bin;
    //printf("   ---> angle_bin: %d, angle_dir_mod: %f \n",angle_bin,angle_dir_mod);
    for(Int_t i_bin_dir = 1; i_bin_dir < 8; i_bin_dir++)
    {
        Int_t delta_bin = (i_bin_dir-1)/2 + 1;
        //printf("i_bin_dir: %d, delta_bin: %d \n",i_bin_dir,delta_bin);
        Int_t bin_dir = angle_bin + angle_dir_sign*delta_bin;
        angle_dir_sign *= -1;
        if(bin_dir < 0) bin_dir += 8;
        if(bin_dir > 7) bin_dir -= 8;
        vec_bins[i_bin_dir] = bin_dir;

        Int_t bin_dir_mod = angle_bin + (-1)*angle_dir_sign_mod*delta_bin;
        angle_dir_sign_mod *= -1;
        if(bin_dir_mod < 0) bin_dir_mod += 8;
        if(bin_dir_mod > 7) bin_dir_mod -= 8;
        vec_bins_mod[i_bin_dir] = bin_dir_mod;
        //printf("   ---> angle_dir: %f, angle_dir_sign: %d, angle_bin: %d, i_bin_dir: %d, bin_dir: %d, bin_dir_mod: %d \n",angle_dir,angle_dir_sign,angle_bin,i_bin_dir,bin_dir,bin_dir_mod);
    }

    if(flip)
    {
        for(Int_t i_dir = 0; i_dir < 8; i_dir++)
        {
            vec_bins[i_dir] = vec_bins_mod[i_dir];
        }
    }
}
//---------------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------------------
Double_t Class_peak_bkgr_cleaner::calc_full_walk_way(vector< vector<Int_t> > vec_walk_path)
{
    Double_t diagonal_move = TMath::Sqrt(2.0); // TMath::Sqrt(1*1 + 1*1)

    // Calculates the true walking path based on the walking bins
    Int_t N_bins_walk = (Int_t)vec_walk_path[0].size();
    Double_t walk_path = 0.0;
    for(Int_t i_walk = 1; i_walk < N_bins_walk; i_walk++)
    {
        Int_t bin_x_A = vec_walk_path[0][i_walk-1];
        Int_t bin_y_A = vec_walk_path[1][i_walk-1];
        Int_t bin_x_B = vec_walk_path[0][i_walk];
        Int_t bin_y_B = vec_walk_path[1][i_walk];

        Int_t delta_x = abs(bin_x_B - bin_x_A);
        Int_t delta_y = abs(bin_y_B - bin_y_A);
        Int_t sum_delta_xy = delta_x + delta_y;

        if(sum_delta_xy == 0) walk_path += 0.0; // shouldn't happen
        if(sum_delta_xy == 1) walk_path += 1.0; // one bin to left/right/up/down
        if(sum_delta_xy > 1) // diagonal move
        {
            walk_path += diagonal_move;
        }

        //printf("i_walk: %d, N_bins_walk: %d, bin_x/y_A -> B: {%d,%d} -> {%d,%d}, walk_path: %f \n",i_walk,N_bins_walk,bin_x_A,bin_y_A,bin_x_B,bin_y_B,walk_path);
    }

    return walk_path;
}
//---------------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------------------
void Class_peak_bkgr_cleaner::do_full_walk(TH2D* h2D_map, Int_t i_cls)
{
    cout << "Class_peak_bkgr_cleaner::do_full_walk" << endl;

    // Create linear 1D list of global bins with bin content > 0
    vector< vector<Int_t> > vec_list_bin_xy_global; // stores the bin_x and bin_y information of bins with bin_cont > 0
    vector<Int_t> vec_list_global_bins; // stores the corresponding global bin
    vec_list_bin_xy_global.resize(2);

    Int_t n_bins_x     = h2D_map->GetNbinsX();
    Int_t n_bins_y     = h2D_map->GetNbinsY();
    N_bins_total = n_bins_x * n_bins_y;

    for(Int_t bin_x = 1; bin_x <= n_bins_x; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= n_bins_y; bin_y++)
        {
            Double_t bin_cont   = h2D_map ->GetBinContent(bin_x,bin_y);
            if(bin_cont > 0)
            {
                Int_t global_bin = h2D_map ->GetBin(bin_x,bin_y);
                vec_list_global_bins.push_back(global_bin);
                vec_list_bin_xy_global[0].push_back(bin_x);
                vec_list_bin_xy_global[1].push_back(bin_y);

                //printf("Linear array, bin_xy: {%d,%d} \n",bin_x,bin_y);
            }
        }
    }

    max_bin_x_init_full_scan   = 0;
    max_bin_y_init_full_scan   = 0;
    max_bin_x_target_full_scan = 0;
    max_bin_y_target_full_scan = 0;
    max_ratio_full_scan        = 0.0;

    vec_walking_bins.clear();
    vec_walking_bins.resize(N_bins_total); // init bins
    for(Int_t i_bin_total_init = 0; i_bin_total_init < N_bins_total; i_bin_total_init++)
    {
        vec_walking_bins[i_bin_total_init].resize(N_bins_total); // target bins
        for(Int_t i_bin_total_target = 0; i_bin_total_target < N_bins_total; i_bin_total_target++)
        {
            vec_walking_bins[i_bin_total_init][i_bin_total_target].resize(2);
            for(Int_t i_flip_walk_vector = 0; i_flip_walk_vector < 2; i_flip_walk_vector++)
            {
                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector].resize(2); // x bin, y bin
            }
        }
    }

    vec_init_target_xy_used.clear();
    vec_init_target_xy_used.resize(2);
    for(Int_t i_flip_walk_vector = 0; i_flip_walk_vector < 2; i_flip_walk_vector++)
    {
        vec_init_target_xy_used[i_flip_walk_vector].resize(2); // init, target
    }

    for(Int_t i_flip_walk_vector = 0; i_flip_walk_vector < 2; i_flip_walk_vector++)
    {
        // Loop over all bin combinations
        for(Int_t i_global_bin_A = 0; i_global_bin_A < (Int_t)vec_list_bin_xy_global[0].size(); i_global_bin_A++)
        {
            Int_t init_bin_x = vec_list_bin_xy_global[0][i_global_bin_A];
            Int_t init_bin_y = vec_list_bin_xy_global[1][i_global_bin_A];
            //for(Int_t i_global_bin_B = i_global_bin_A+1; i_global_bin_B < (Int_t)vec_list_bin_xy_global[0].size(); i_global_bin_B++)
            for(Int_t i_global_bin_B = 0; i_global_bin_B < (Int_t)vec_list_bin_xy_global[0].size(); i_global_bin_B++) // A->B and B->A might not be the same, do full scan!
            {
                if(i_global_bin_A == i_global_bin_B) continue;
                Int_t target_bin_x = vec_list_bin_xy_global[0][i_global_bin_B];
                Int_t target_bin_y = vec_list_bin_xy_global[1][i_global_bin_B];

                //---------------------------------------
                Int_t flag_target_found   = 0;
                Int_t walking_step        = 0;
                Int_t previous_bin_x      = init_bin_x;
                Int_t previous_bin_y      = init_bin_y;
                Int_t flag_Z_bin_previous = 0;
                Int_t flag_Z_bin          = 0;

                vector< vector<Int_t> > vec_flag_bin_used_as_close_bin;
                vec_flag_bin_used_as_close_bin.resize(2);

                vector< vector<Int_t> > vec_walking_steps;
                vec_walking_steps.resize(n_bins_x);
                for(Int_t bin_x = 0; bin_x < n_bins_x; bin_x++)
                {
                    vec_walking_steps[bin_x].resize(n_bins_y);
                    for(Int_t bin_y = 0; bin_y < n_bins_y; bin_y++)
                    {
                        vec_walking_steps[bin_x][bin_y] = 0;
                    }
                }
                vec_walking_steps[previous_bin_x-1][previous_bin_y-1] = 1; // histogram bins start from 1, flag init bin as used

                Int_t i_bin_total_init   = (init_bin_x-1)   + n_bins_x*(init_bin_y-1);
                Int_t i_bin_total_target = (target_bin_x-1) + n_bins_x*(target_bin_y-1);

                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].push_back(init_bin_x);
                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].push_back(init_bin_y);


                //--------------------------------------------------------
                // Store global init and target bins
                vec_init_target_xy_used[i_flip_walk_vector][0].push_back(i_bin_total_init);
                vec_init_target_xy_used[i_flip_walk_vector][1].push_back(i_bin_total_target);

                //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                //{
                //    printf("i_bin_total_init: %d, i_bin_total_target: %d \n",i_bin_total_init,i_bin_total_target);
                //}
                //--------------------------------------------------------

                // Single init target finding
                while(!flag_target_found)
                {
                    calc_vec_bins(previous_bin_x,previous_bin_y,target_bin_x,target_bin_y,i_flip_walk_vector); // gives "best" next bin direction(s), sets vec_bins

                    if(flag_Z_bin) flag_Z_bin_previous = 1;
                    else flag_Z_bin_previous = 0;
                    flag_Z_bin = 0;

                    //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                    //{
                    //    printf("Start searching for bin direction \n");
                    //}

                    for(Int_t i_bin_dir = 0; i_bin_dir < 8; i_bin_dir++) // ordered 8 next possible directions
                    {
                        Int_t flag_dead_end = 0;
                        Int_t flag_X_bin    = 0;

                        Int_t bin_dir = vec_bins[i_bin_dir];
                        Int_t bin_x_search = previous_bin_x + vec_delta_bin_x[bin_dir];
                        Int_t bin_y_search = previous_bin_y + vec_delta_bin_y[bin_dir];

                        //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                        //{
                        //    printf("i_bin_dir: %d, bin_dir: %d, delta_x/y: {%d,%d}, bin_x/y_search: {%d,%d} \n",i_bin_dir,bin_dir,vec_delta_bin_x[bin_dir],vec_delta_bin_y[bin_dir],bin_x_search,bin_y_search);
                        //}


                        Int_t flag_out_of_bounce = 0;
                        if(bin_x_search < 1 || bin_y_search < 1 || bin_x_search > n_bins_x || bin_y_search > n_bins_y)
                        {
                            flag_out_of_bounce = 1;

                            //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                            //{
                            //    printf("Out of bounce \n");
                            //}

                            if(i_bin_dir != 7) continue;
                            else flag_dead_end = 1;
                        }
                        Double_t bin_cont = h2D_cluster ->GetBinContent(bin_x_search,bin_y_search);
                        if(bin_cont <= 0.0)
                        {
                            //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                            //{
                            //    printf("Bin entry is ZERO \n");
                            //}
                            if(i_bin_dir != 7) continue;
                            else flag_dead_end = 1;
                        }

                        if(!flag_out_of_bounce && vec_walking_steps[bin_x_search-1][bin_y_search-1]) // histogram bins start from 1
                        {
                            //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                            //{
                            //    printf("Bin was already found \n");
                            //}
                            if(i_bin_dir != 7) continue; // bin was already used
                            else flag_dead_end = 1;
                        }




                        if(!flag_out_of_bounce && !flag_X_bin && !flag_Z_bin) vec_walking_steps[bin_x_search-1][bin_y_search-1] = 1; // histogram bins start from 1

                        if(!flag_dead_end) // good bin so far
                        {
                            //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                            //{
                            //    printf("GOOD \n");
                            //}
                            vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].push_back(bin_x_search);
                            vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].push_back(bin_y_search);
                        }


                        if(flag_dead_end)
                        {
                            //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                            //{
                            //    printf("Dead end \n");
                            //}

                            Int_t walking_step_previous = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size() - 2;
                            if(walking_step_previous < 0) //break; // CHECK
                            {
                                walking_step_previous = 0;
                            }
                            else
                            {
                                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].erase(vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].end()-1,vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].end());
                                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].erase(vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].end()-1,vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].end());
                            }

                            bin_x_search = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0][walking_step_previous];
                            bin_y_search = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1][walking_step_previous];

                        }


                        previous_bin_x = bin_x_search;
                        previous_bin_y = bin_y_search;

                        Int_t walking_step_used = (Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size() - 1; // init is already first element
                        walking_step++;


                        //if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                        //{
                        //    printf("last i_bin_dir: %d, bin_dir: %d, delta_x/y: {%d,%d}, bin_x/y_search: {%d,%d}, walking_step_used: %d \n",i_bin_dir,bin_dir,vec_delta_bin_x[bin_dir],vec_delta_bin_y[bin_dir],bin_x_search,bin_y_search,walking_step_used);
                        //}


                        if((bin_x_search == target_bin_x && bin_y_search == target_bin_y)
                           || walking_step > 250
                          )
                        {
                            flag_target_found = 1;
                        }
                        break;
                    }  // end of ordered 8 directions
                } // end while(!flag_target_found)
                //---------------------------------------


                //---------------------------------------
                // Optimize found path
                // Check if a neighbouring bin (not the last one) was found -> shorter path
                Int_t flag_end_optimize = 0;

                while(!flag_end_optimize)
                {
                    Int_t flag_erase = 0;

                    // Loop over all steps
                    for(Int_t i_walk = 0; i_walk  < (Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size(); i_walk++)
                    {
                        Int_t bin_x = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0][i_walk];
                        Int_t bin_y = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1][i_walk];

                        // Check all following steps, starting from the overnext one
                        for(Int_t i_walk_next = i_walk+2; i_walk_next  < (Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size(); i_walk_next++)
                        {
                            Int_t bin_x_next = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0][i_walk_next];
                            Int_t bin_y_next = vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1][i_walk_next];

                            if(fabs(bin_x_next - bin_x) <= 1 && fabs(bin_y_next - bin_y) <= 1) // close bin found
                            {
                                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].erase(vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].begin()+i_walk+1,vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].begin()+i_walk_next);
                                vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].erase(vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].begin()+i_walk+1,vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].begin()+i_walk_next);

                                flag_erase = 1;
                            }
                            if(flag_erase) break;
                        }
                        if(flag_erase) break;

                        if(i_walk == ((Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size() - 1)) flag_end_optimize = 1;
                    }
                }
                // End of optimization

                if(i_cls == 3 && init_bin_x == 35 && init_bin_y == 33 && target_bin_x == 34 && target_bin_y == 33)
                {
                    printf("size after optimization: %d \n",(Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].size() - 1);
                }
                //---------------------------------------


            } // end global_bin_B
        } // end global_bin_A
    } // end of flip direction vector

    //---------------------------------------
    printf("-------------------------- Check for biggest difference -------------------------- \n");

    Int_t max_i_bin_total_init   = 0;
    Int_t max_i_bin_total_target = 0;
    Int_t max_path_total         = 0;
    Int_t max_bin_x_init         = 0; // max_bin_x/y_init/target define the two bins which have the largest separation -> used as reference to define sub clusters
    Int_t max_bin_y_init         = 0;
    Int_t max_bin_x_target       = 0;
    Int_t max_bin_y_target       = 0;
    Int_t max_idx_init           = 0;
    Int_t max_idx_target         = 0;

    Double_t max_ratio = 0.0;
    for(Int_t i_path = 0; i_path < (Int_t)vec_init_target_xy_used[0][0].size(); i_path++)
    {
        Int_t i_bin_total_init   = vec_init_target_xy_used[0][0][i_path];
        Int_t i_bin_total_target = vec_init_target_xy_used[0][1][i_path];

        // Four different walk ways per singl init-target bin combination
        // A->B and B->A bins can be different
        // flip is changing the direction vectors of the walking algorithm
        Double_t total_walk_way_A_to_B      = calc_full_walk_way(vec_walking_bins[i_bin_total_init][i_bin_total_target][0]);
        Double_t total_walk_way_B_to_A      = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init][0]);
        Double_t total_walk_way_A_to_B_flip = calc_full_walk_way(vec_walking_bins[i_bin_total_init][i_bin_total_target][1]);
        Double_t total_walk_way_B_to_A_flip = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init][1]);

        //printf("A->B: %f, B->A: %f, A->B flip: %f, B->A flip: %f \n",total_walk_way_A_to_B,total_walk_way_B_to_A,total_walk_way_A_to_B_flip,total_walk_way_B_to_A_flip);

        //vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][0].push_back(bin_x_search);
        //vec_walking_bins[i_bin_total_init][i_bin_total_target][i_flip_walk_vector][1].push_back(bin_y_search);

        Int_t    bin_x_init   = vec_walking_bins[i_bin_total_init][i_bin_total_target][0][0][0]; // first entry in walking path
        Int_t    bin_y_init   = vec_walking_bins[i_bin_total_init][i_bin_total_target][0][1][0];
        Int_t    bin_x_target = vec_walking_bins[i_bin_total_init][i_bin_total_target][0][0][(Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][0][0].size()-1]; // last entry in walking path
        Int_t    bin_y_target = vec_walking_bins[i_bin_total_init][i_bin_total_target][0][1][(Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][0][1].size()-1];
        Double_t delta_x      = (Double_t)(bin_x_target - bin_x_init);
        Double_t delta_y      = (Double_t)(bin_y_target - bin_y_init);
        Double_t bird_way     = (Double_t)TMath::Sqrt(delta_x*delta_x + delta_y*delta_y); // direct path, like the bird would fly :)


        //if(i_cls == 3 && bin_x_init == 35 && bin_y_init == 33 && bin_x_target == 34 && bin_y_target == 33)
        //{
        //    printf("i_path: %d, size at path: %d, i_bin_total_init: %d, i_bin_total_target: %d  \n",i_path,(Int_t)vec_walking_bins[i_bin_total_init][i_bin_total_target][0][0].size() - 1,i_bin_total_init,i_bin_total_target);
        //}

        if(bird_way <= 0) bird_way = 1.0; // Shouldn't happen...
        //Double_t ratio_walk_to_bird = total_walk_way/bird_way;
        Double_t ratio_walk_to_bird_A_to_B      = TMath::Power(total_walk_way_A_to_B,walk_to_bird_power)/bird_way;
        Double_t ratio_walk_to_bird_B_to_A      = TMath::Power(total_walk_way_B_to_A,walk_to_bird_power)/bird_way;
        Double_t ratio_walk_to_bird_A_to_B_flip = TMath::Power(total_walk_way_A_to_B_flip,walk_to_bird_power)/bird_way;
        Double_t ratio_walk_to_bird_B_to_A_flip = TMath::Power(total_walk_way_B_to_A_flip,walk_to_bird_power)/bird_way;

        // Take the minimum for all four combinations, A->B, B->A and with the flipped direction vectors for the searching algorithm
        // They have all the same bird way, so the minimum of the ratio is also the optimal (shortest) path
        Double_t min_A = TMath::Min(ratio_walk_to_bird_A_to_B,ratio_walk_to_bird_B_to_A);
        Double_t min_B = TMath::Min(min_A,ratio_walk_to_bird_A_to_B_flip);
        Double_t min_C = TMath::Min(min_B,ratio_walk_to_bird_B_to_A_flip);

        Double_t ratio_walk_to_bird = min_C;

        //printf("A->B: %f, B->A: %f, A->B flip: %f, B->A flip: %f, ratio_walk_to_bird: %f \n",ratio_walk_to_bird_A_to_B,ratio_walk_to_bird_B_to_A,ratio_walk_to_bird_A_to_B_flip,ratio_walk_to_bird_B_to_A_flip,ratio_walk_to_bird);

        if(ratio_walk_to_bird > max_ratio)
        {
            max_ratio              = ratio_walk_to_bird;
            max_i_bin_total_init   = i_bin_total_init;
            max_i_bin_total_target = i_bin_total_target;
            max_path_total         = i_path;
            max_bin_x_init         = bin_x_init;
            max_bin_y_init         = bin_y_init;
            max_bin_x_target       = bin_x_target;
            max_bin_y_target       = bin_y_target;

        }

        //if(i_cls == 3) printf("full loop, i_cls: %d, i_path: %d, total_walk_way_A_to_B: %f, bird_way: %f, ratio: %f, init: {%d,%d}, target: {%d,%d} \n",i_cls,i_path,total_walk_way_A_to_B,bird_way,ratio_walk_to_bird,bin_x_init,bin_y_init,bin_x_target,bin_y_target);
    }

    max_bin_x_init_full_scan         = max_bin_x_init;
    max_bin_y_init_full_scan         = max_bin_y_init;
    max_bin_x_target_full_scan       = max_bin_x_target;
    max_bin_y_target_full_scan       = max_bin_y_target;
    max_ratio_full_scan              = max_ratio;
    max_i_bin_total_init_full_scan   = max_i_bin_total_init;
    max_i_bin_total_target_full_scan = max_i_bin_total_target;


    printf("========> full loop: max_ratio: %f, max_i_bin_total_init: %d, max_i_bin_total_target: %d, init bins x,y: {%d,%d}, target bins x,y: {%d,%d}, max_idx_init \n",max_ratio,max_i_bin_total_init,max_i_bin_total_target,max_bin_x_init,max_bin_y_init,max_bin_x_target,max_bin_y_target);

    //---------------------------------------

}
//---------------------------------------------------------------------------------------------------------------------



#if 1
//---------------------------------------------------------------------------------------------------------------------
Int_t Class_peak_bkgr_cleaner::categorize_sub_clusters(TH2D* h2D_map, Int_t i_cls)
{
    // Categorize two sub clusters
    cout << "" << endl;
    printf("Class_peak_bkgr_cleaner::categorize_sub_clusters: start categorize two sub clusters \n");


    Int_t n_bins_x     = h2D_map->GetNbinsX();
    Int_t n_bins_y     = h2D_map->GetNbinsY();
    h2D_sub_cluster_map = new TH2D("h2D_sub_cluster_map","h2D_sub_cluster_map",n_bins_x,0,n_bins_x,n_bins_y,0,n_bins_y);

    // Set init and target bins which define the categories
    h2D_sub_cluster_map ->SetBinContent(max_bin_x_init_full_scan,max_bin_y_init_full_scan,1);
    h2D_sub_cluster_map ->SetBinContent(max_bin_x_target_full_scan,max_bin_y_target_full_scan,2);


    for(Int_t i_bin_total_target = 0; i_bin_total_target < N_bins_total; i_bin_total_target++)
    {
        // Those are the two bins determined earlier in the do_full_walk function. They are the reference bins which define what is in the  catergory X and category Y sub clusters
        Int_t i_bin_total_init_cat_X   = max_i_bin_total_init_full_scan;
        Int_t i_bin_total_init_cat_Y   = max_i_bin_total_target_full_scan;

        if(i_bin_total_init_cat_X == i_bin_total_target) continue;
        if(i_bin_total_init_cat_Y == i_bin_total_target) continue;

        // Four different walk ways per singl init-target bin combination
        // A->B and B->A bins can be different
        // flip is changing the direction vectors of the walking algorithm

        // Keep always the same target bin.
        Double_t total_walk_way_A_to_B_cat_X      = calc_full_walk_way(vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0]);
        //printf("i_bin_total_target: %d, total_walk_way_A_to_B_cat_X: %f \n",i_bin_total_target,total_walk_way_A_to_B_cat_X);
        if(total_walk_way_A_to_B_cat_X <= 0.0) continue;

        Double_t total_walk_way_B_to_A_cat_X      = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init_cat_X][0]);
        if(total_walk_way_B_to_A_cat_X <= 0.0) continue;

        Double_t total_walk_way_A_to_B_flip_cat_X = calc_full_walk_way(vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][1]);
        if(total_walk_way_A_to_B_flip_cat_X <= 0.0) continue;

        Double_t total_walk_way_B_to_A_flip_cat_X = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init_cat_X][1]);
        if(total_walk_way_B_to_A_flip_cat_X <= 0.0) continue;

        //printf("walk ways cat X: {%f,%f,%f,%f} \n",total_walk_way_A_to_B_cat_X,total_walk_way_B_to_A_cat_X,total_walk_way_A_to_B_flip_cat_X,total_walk_way_B_to_A_flip_cat_X);

        // Take the minimum for all four combinations, A->B, B->A and with the flipped direction vectors for the searching algorithm
        // They have all the same bird way, so the minimum of the ratio is also the optimal (shortest) path
        Double_t min_A_cat_X              = TMath::Min(total_walk_way_A_to_B_cat_X,total_walk_way_B_to_A_cat_X);
        Double_t min_B_cat_X              = TMath::Min(min_A_cat_X,total_walk_way_A_to_B_flip_cat_X);
        Double_t min_total_walk_way_cat_X = TMath::Min(min_B_cat_X,total_walk_way_B_to_A_flip_cat_X);

        Double_t total_walk_way_A_to_B_cat_Y = calc_full_walk_way(vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0]);
        if(total_walk_way_A_to_B_cat_Y <= 0.0) continue;

        Double_t total_walk_way_B_to_A_cat_Y = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init_cat_Y][0]);
        if(total_walk_way_B_to_A_cat_Y <= 0.0) continue;

        Double_t total_walk_way_A_to_B_flip_cat_Y = calc_full_walk_way(vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][1]);
        if(total_walk_way_A_to_B_flip_cat_Y <= 0.0) continue;

        Double_t total_walk_way_B_to_A_flip_cat_Y = calc_full_walk_way(vec_walking_bins[i_bin_total_target][i_bin_total_init_cat_Y][1]);
        if(total_walk_way_B_to_A_flip_cat_Y <= 0.0) continue;

        // Take the minimum for all four combinations, A->B, B->A and with the flipped direction vectors for the searching algorithm
        // They have all the same bird way, so the minimum of the ratio is also the optimal (shortest) path
        Double_t min_A_cat_Y              = TMath::Min(total_walk_way_A_to_B_cat_Y,total_walk_way_B_to_A_cat_Y);
        Double_t min_B_cat_Y              = TMath::Min(min_A_cat_Y,total_walk_way_A_to_B_flip_cat_Y);
        Double_t min_total_walk_way_cat_Y = TMath::Min(min_B_cat_Y,total_walk_way_B_to_A_flip_cat_Y);

        // Calculate bird way category X
        //printf("N_bins_total: %d, i_bin_total_init_cat_X: %d, i_bin_total_target: %d \n",N_bins_total,i_bin_total_init_cat_X,i_bin_total_target);
        Int_t    bin_x_init_cat_X   = vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][0][0];
        Int_t    bin_y_init_cat_X   = vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][1][0];
        Int_t    bin_x_target_cat_X = vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][0][(Int_t)vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][0].size()-1];
        Int_t    bin_y_target_cat_X = vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][1][(Int_t)vec_walking_bins[i_bin_total_init_cat_X][i_bin_total_target][0][0].size()-1];
        Double_t delta_x_cat_X      = (Double_t)(bin_x_target_cat_X - bin_x_init_cat_X);
        Double_t delta_y_cat_X      = (Double_t)(bin_y_target_cat_X - bin_y_init_cat_X);
        Double_t bird_way_cat_X     = (Double_t)TMath::Sqrt(delta_x_cat_X*delta_x_cat_X + delta_y_cat_X*delta_y_cat_X); // direct path, like the bird would fly :)

        // Calculate bird way catergory Y
        Int_t    bin_x_init_cat_Y   = vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][0][0];
        Int_t    bin_y_init_cat_Y   = vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][1][0];
        Int_t    bin_x_target_cat_Y = vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][0][(Int_t)vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][0].size()-1];
        Int_t    bin_y_target_cat_Y = vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][1][(Int_t)vec_walking_bins[i_bin_total_init_cat_Y][i_bin_total_target][0][0].size()-1];
        Double_t delta_x_cat_Y      = (Double_t)(bin_x_target_cat_Y - bin_x_init_cat_Y);
        Double_t delta_y_cat_Y      = (Double_t)(bin_y_target_cat_Y - bin_y_init_cat_Y);
        Double_t bird_way_cat_Y     = (Double_t)TMath::Sqrt(delta_x_cat_Y*delta_x_cat_Y + delta_y_cat_Y*delta_y_cat_Y); // direct path, like the bird would fly :)

        if(bird_way_cat_X <= 0) bird_way_cat_X = 1.0; // Shouldn't happen...
        Double_t ratio_walk_to_bird_cat_X = TMath::Power(min_total_walk_way_cat_X,walk_to_bird_power)/bird_way_cat_X;

        if(bird_way_cat_Y <= 0) bird_way_cat_Y = 1.0; // Shouldn't happen...
        Double_t ratio_walk_to_bird_cat_Y = TMath::Power(min_total_walk_way_cat_Y,walk_to_bird_power)/bird_way_cat_Y;

        Int_t cat_XY = 2;

        if(ratio_walk_to_bird_cat_X < ratio_walk_to_bird_cat_Y) cat_XY = 1;
        if((min_total_walk_way_cat_X/bird_way_cat_X) < 2.0 && (min_total_walk_way_cat_Y/bird_way_cat_Y) < 2.0)
        {
            if(min_total_walk_way_cat_X < min_total_walk_way_cat_Y) cat_XY = 1;
            else cat_XY = 2;
        }

        //if(i_cls == 14)
        //{
        //    printf("cat X, x/y init->target: {%d,%d}->{%d,%d}, cat Y, x/y init->target: {%d,%d}->{%d,%d} \n",bin_x_init_cat_X,bin_y_init_cat_X,bin_x_target_cat_X,bin_y_target_cat_X,bin_x_init_cat_Y,bin_y_init_cat_Y,bin_x_target_cat_Y,bin_y_target_cat_Y);
        //    printf("cat X/Y bird: {%f,%f}, cat X/Y walk: {%f,%f} \n",bird_way_cat_X,bird_way_cat_Y,min_total_walk_way_cat_X,min_total_walk_way_cat_Y);
        //}

        h2D_sub_cluster_map ->SetBinContent(bin_x_target_cat_X,bin_y_target_cat_X,cat_XY);
    }

    // Update global cluster map with sub clusters
    // printf("--------------> i_cls: %d \n",i_cls);
    N_new_split_cluster++;
    for(Int_t bin_x_target_use = 1; bin_x_target_use <= h2D_sub_cluster_map->GetNbinsX(); bin_x_target_use++)
    {
        for(Int_t bin_y_target_use = 1; bin_y_target_use <= h2D_sub_cluster_map->GetNbinsY(); bin_y_target_use++)
        {
            Int_t sub_cluster_idx = (Int_t)h2D_sub_cluster_map->GetBinContent(bin_x_target_use,bin_y_target_use);
            if(sub_cluster_idx == 0) continue;

            //printf("bin_xy: {%d,%d} sub_cluster_idx: %d \n",bin_x_target_use,bin_y_target_use,sub_cluster_idx);

            Int_t x_bin_full = bin_x_target_use + vec_min_max_x_y_cluster[i_cls][0][0] - 1;
            Int_t y_bin_full = bin_y_target_use + vec_min_max_x_y_cluster[i_cls][1][0] - 1;

            Int_t global_cluster_idx = (Int_t)h2D_2D_bkgr_clusters ->GetBinContent(x_bin_full,y_bin_full);

            Int_t new_idx = global_cluster_idx;
            //if(sub_cluster_idx == 2) new_idx = max_cls_index + 1 + global_cluster_idx;
            if(sub_cluster_idx == 2) new_idx = max_cls_index + N_new_split_cluster;

            h2D_2D_bkgr_clusters_itt ->SetBinContent(x_bin_full,y_bin_full,new_idx);
        }
    }

    return 1;
}
//---------------------------------------------------------------------------------------------------------------------
#endif


//---------------------------------------------------------------------------------------------------------------------
void Class_peak_bkgr_cleaner::create_heat_maps()
{
    //---------------------------
    // Identify ranges in x and y for each cluster and caluclate how large the clusters are

    for(Int_t itt_cluster_separation = 0; itt_cluster_separation < N_itt_cluster_separation; itt_cluster_separation++)
    {
        N_new_split_cluster = 0;
        printf("Class_peak_bkgr_cleaner::create_heat_maps: itt_cluster_separation: %d, itt_cluster_separation_plot: %d \n",itt_cluster_separation,itt_cluster_separation_plot);
        if(itt_cluster_separation > 0)
        {
            if(h2D_2D_bkgr_clusters) h2D_2D_bkgr_clusters ->Delete();
            h2D_2D_bkgr_clusters = (TH2D*)h2D_2D_bkgr_clusters_itt ->Clone("h2D_2D_bkgr_clusters");
        }
        Int_t max_x_bin = h2D_2D_bkgr_clusters ->GetNbinsX();
        Int_t max_y_bin = h2D_2D_bkgr_clusters ->GetNbinsY();

        max_cls_index = 0;

        vector<Int_t> vec_N_cluster;
        for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
        {
            for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
            {
                Double_t pos_x     =  h2D_2D_bkgr_clusters ->GetXaxis()->GetBinCenter(bin_x);
                Double_t pos_y     =  h2D_2D_bkgr_clusters ->GetYaxis()->GetBinCenter(bin_y);
                Int_t    cls_index =  (Int_t)h2D_2D_bkgr_clusters ->GetBinContent(bin_x,bin_y);
                if(cls_index < 0) continue;

                //printf("bin_x/y: {%d,%d}, cls_index: %d \n",bin_x,bin_y,cls_index);

                if(cls_index > max_cls_index) max_cls_index = cls_index;

                if((Int_t)vec_N_cluster.size() <= cls_index)
                {
                    vec_N_cluster.resize(cls_index+1);
                    vec_min_max_x_y_cluster.resize(cls_index+1);
                    for(Int_t i_x_y = 0; i_x_y < 2; i_x_y++)
                    {
                        vec_min_max_x_y_cluster[cls_index].resize(2);
                        for(Int_t i_min_max = 0; i_min_max < 2; i_min_max++)
                        {
                            vec_min_max_x_y_cluster[cls_index][i_x_y].resize(2);
                        }
                    }
                    vec_min_max_x_y_cluster[cls_index][0][0] = bin_x; // bin_x min
                    vec_min_max_x_y_cluster[cls_index][0][1] = bin_x; // bin_x max
                    vec_min_max_x_y_cluster[cls_index][1][0] = bin_y; // bin_y min
                    vec_min_max_x_y_cluster[cls_index][1][1] = bin_y; // bin_y max
                }
                if((Int_t)vec_min_max_x_y_cluster[cls_index].size() == 0)
                {
                    for(Int_t i_x_y = 0; i_x_y < 2; i_x_y++)
                    {
                        vec_min_max_x_y_cluster[cls_index].resize(2);
                        for(Int_t i_min_max = 0; i_min_max < 2; i_min_max++)
                        {
                            vec_min_max_x_y_cluster[cls_index][i_x_y].resize(2);
                        }
                    }
                    vec_min_max_x_y_cluster[cls_index][0][0] = bin_x; // bin_x min
                    vec_min_max_x_y_cluster[cls_index][0][1] = bin_x; // bin_x max
                    vec_min_max_x_y_cluster[cls_index][1][0] = bin_y; // bin_y min
                    vec_min_max_x_y_cluster[cls_index][1][1] = bin_y; // bin_y max
                }



                if(bin_x < vec_min_max_x_y_cluster[cls_index][0][0]) vec_min_max_x_y_cluster[cls_index][0][0] = bin_x;
                if(bin_x > vec_min_max_x_y_cluster[cls_index][0][1]) vec_min_max_x_y_cluster[cls_index][0][1] = bin_x;
                if(bin_y < vec_min_max_x_y_cluster[cls_index][1][0]) vec_min_max_x_y_cluster[cls_index][1][0] = bin_y;
                if(bin_y > vec_min_max_x_y_cluster[cls_index][1][1]) vec_min_max_x_y_cluster[cls_index][1][1] = bin_y;
                vec_N_cluster[cls_index]++;
            }
        }
        //---------------------------


        //---------------------------
        for(Int_t i_cls = 0; i_cls < (Int_t)vec_N_cluster.size(); i_cls++)
        {
            Int_t flag_walk_to_bird_ratio = 1;
            Int_t cls_size = vec_N_cluster[i_cls];
            //printf("i_cls: %d, out of %d, cls_size: %d \n",i_cls,(Int_t)vec_N_cluster.size(),cls_size);

            if(cls_size < min_cluster_size_to_split) continue;
            cout << "" << endl;
            cout << "===========================================================" << endl;
            //printf("i_cls: %d, size: %d \n",i_cls,vec_N_cluster[i_cls]);

            //printf("bin x range: {%d,%d} \n",vec_min_max_x_y_cluster[i_cls][0][0],vec_min_max_x_y_cluster[i_cls][0][1]);
            //printf("bin y range: {%d,%d} \n",vec_min_max_x_y_cluster[i_cls][1][0],vec_min_max_x_y_cluster[i_cls][1][1]);

            // Create a 2D size optimized histogram only for current cluster i_cls
            Int_t n_bins_x = vec_min_max_x_y_cluster[i_cls][0][1] - vec_min_max_x_y_cluster[i_cls][0][0] + 1;
            Int_t n_bins_y = vec_min_max_x_y_cluster[i_cls][1][1] - vec_min_max_x_y_cluster[i_cls][1][0] + 1;
            h2D_cluster = new TH2D("h2D_cluster","h2D_cluster",n_bins_x,0,n_bins_x,n_bins_y,0,n_bins_y);

            Int_t N_h2D_cluster_size = 0;
            for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
            {
                for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
                {
                    Int_t cls_index =  (Int_t)h2D_2D_bkgr_clusters ->GetBinContent(bin_x,bin_y);
                    //if(cls_index == 0) printf("cls_index: %d, bin_x/y: {%d,%d} \n",cls_index,bin_x,bin_y);
                    if(cls_index != i_cls) continue;
                    Double_t pos_x     =  h2D_2D_bkgr_clusters ->GetXaxis()->GetBinCenter(bin_x);
                    Double_t pos_y     =  h2D_2D_bkgr_clusters ->GetYaxis()->GetBinCenter(bin_y);

                    h2D_cluster ->SetBinContent(bin_x - vec_min_max_x_y_cluster[i_cls][0][0] + 1,bin_y - vec_min_max_x_y_cluster[i_cls][1][0] + 1,1);
                    N_h2D_cluster_size++;
                }
            }

            do_full_walk(h2D_cluster,i_cls);
            if(max_ratio_full_scan < min_walk_to_bird_ratio)
            {
                if(h2D_cluster) h2D_cluster ->Delete();
                continue;
            }
            Int_t cat_return = categorize_sub_clusters(h2D_cluster,i_cls);

            if(i_cls == cluster_plot && itt_cluster_separation == itt_cluster_separation_plot)
            {
                flag_cluster_plot_available = kTRUE;
                h2D_sub_cluster_map_out         = (TH2D*)h2D_sub_cluster_map->Clone("h2D_sub_cluster_map_out");
                vec_walking_bins_plot = vec_walking_bins;
            }
            delete h2D_sub_cluster_map;

            h2D_cluster         ->Delete();
        } // end of cluster loop

        vec_min_max_x_y_cluster.clear();
        //---------------------------

        cout << "end of cluster separation itteration loop" << endl;
    } // end of cluster separation itteration loop

}

vector< vector< vector< vector< vector<Int_t> > > > >  Class_peak_bkgr_cleaner::get_vec_walking_bins_plot()
{
    return vec_walking_bins_plot;
}

void Class_peak_bkgr_cleaner::select_final_bkgr_clusters()
{
    cout << "" << endl;
    printf("Class_peak_bkgr_cleaner::select_final_bkgr_clusters() \n");

    //h2D_2D_bkgr_clusters_itt
    Int_t max_x_bin =  h2D_2D_bkgr_clusters_itt->GetNbinsX();
    Int_t max_y_bin =  h2D_2D_bkgr_clusters_itt->GetNbinsY();

    vec_cluster_size.clear();

    // Determine number of clusters and cluster sizes
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Int_t cls_indx = (Int_t)h2D_2D_bkgr_clusters_itt ->GetBinContent(bin_x,bin_y); // start from 0 for real clusters
            if(cls_indx < 0) continue;
            if(cls_indx >= vec_cluster_size.size())
            {
                vec_cluster_size.resize(cls_indx+1);
                vec_cluster_size[cls_indx] = 1;
            }
            vec_cluster_size[cls_indx]++;
        }
    }


    vector< vector<Double_t> > vec_min_max_dist;
    vec_min_max_dist.resize((Int_t)vec_cluster_size.size());
    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        vec_min_max_dist[i_cls].resize(2); // min, max distance to center
        vec_min_max_dist[i_cls][0] = 9999999.0;
        vec_min_max_dist[i_cls][1] = 0.0;
    }



    // Determine minimium and maximum distance to center for each cluster
    for(Int_t bin_x = 1; bin_x <= max_x_bin; bin_x++)
    {
        for(Int_t bin_y = 1; bin_y <= max_y_bin; bin_y++)
        {
            Double_t pos_x = h2D_2D_bkgr_clusters_itt ->GetXaxis()->GetBinCenter(bin_x);
            Double_t pos_y = h2D_2D_bkgr_clusters_itt ->GetYaxis()->GetBinCenter(bin_y);
            Double_t radius = TMath::Sqrt(pos_x*pos_x + pos_y*pos_y);

            Int_t cls_indx = (Int_t)h2D_2D_bkgr_clusters_itt ->GetBinContent(bin_x,bin_y); // start from 0 for real clusters
            if(cls_indx < 0) continue;
            Int_t cls_size = vec_cluster_size[cls_indx];

            if(cls_size < min_cls_size) continue; // select only clusters with a minimum amount of points
            if(radius < vec_min_max_dist[cls_indx][0]) vec_min_max_dist[cls_indx][0] = radius;
            if(radius > vec_min_max_dist[cls_indx][1]) vec_min_max_dist[cls_indx][1] = radius;
        }
    }

    for(Int_t i_cls = 0; i_cls < (Int_t)vec_cluster_size.size(); i_cls++)
    {
        //printf("i_cls: %d, cluster size: %d, min/max radius: {%f,%f} \n",i_cls,vec_cluster_size[i_cls],vec_min_max_dist[i_cls][0],vec_min_max_dist[i_cls][1]);
    }

}


TH2D* Class_peak_bkgr_cleaner::get_TH2D()
{
    return h2D_2D_peak_points;
}

TH2D* Class_peak_bkgr_cleaner::get_bkgr_TH2D()
{
    return h2D_2D_bkgr_clusters_itt;
}

void Class_peak_bkgr_cleaner::clear()
{
    if(h2D_2D_peak_points) h2D_2D_peak_points          ->Delete();
    h2D_2D_peak_points = NULL;
    if(h2D_2D_peak_points_clusters) h2D_2D_peak_points_clusters ->Delete();
    h2D_2D_peak_points_clusters = NULL;
    if(h2D_2D_bkgr_clusters) h2D_2D_bkgr_clusters        ->Delete();
    h2D_2D_bkgr_clusters = NULL;
    //vec_2D_peak_points.clear();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Class_2D_crater_finder
{
    std::vector< std::vector<Double_t> > vec_2D_peak_points; // [point number][x,y,height]
    Double_t search_width;
    Double_t min_radius;
    Double_t max_radius;
    Double_t search_step_size;
    Double_t x_center, y_center;
    Double_t exclude_width;
    Int_t    N_maxima;
    Double_t min_points_per_area; // minimum number of points per area for a ring to be considered
    vector< vector<Double_t> > vec_radii_peak_positions_quality;
    vector< vector<Double_t> > vec_radii_peak_positions_radius;
    vector< vector<Double_t> > vec_radii_peak_positions_radius_merge;
    vector< vector< vector<Double_t> > > vec_radii_partial_ring_angles;
    vector< vector< vector<Double_t> > > vec_radii_partial_ring_angles_radius;
    vector< vector< vector<Double_t> > > vec_radii_partial_ring_angles_merge;
    vector< vector<Double_t> > vec_radii_peak_positions_radius_merge_copy;
    vector< vector< vector<Double_t> > > vec_radii_partial_ring_angles_merge_copy;
    vector<Int_t> index_quality_to_radius;
    vector<Int_t> index_radius_to_quality;
    Int_t flag_debug;
    Int_t param_radii;

    Double_t ring_quality;
    Double_t ring_points;
public:
    void set_debug(Int_t flag_debug_in);
    void add_data(std::vector< std::vector<Double_t> > vec_2D_peak_points_in);
    void set_finding_parameters(Double_t search_width_in,Double_t min_radius_in,Double_t max_radius_in,
                                Double_t search_step_size_in,
                                Double_t min_points_per_area_in, Double_t exclude_width_in,
                                Double_t x_center_in,Double_t y_center_in,Int_t N_maxima_in);
    void order_radii();
    void find_craters();
    void Merge_close_rings(Double_t ring_dist_merge);
    vector< vector<Double_t> > get_radii_peak_positions();
    vector< vector<Double_t> > get_radii_peak_positions_merged();
    vector< vector<Double_t> > get_radii_peak_positions_ordered_by_radius();
    vector< vector< vector<Double_t> > > get_radii_partial_ring_angles();
    vector< vector< vector<Double_t> > > get_radii_partial_ring_angles_merged();
    Double_t get_ring_quality();
};

void Class_2D_crater_finder::set_debug(Int_t flag_debug_in)
{
    flag_debug = flag_debug_in;
}

void Class_2D_crater_finder::add_data(std::vector< std::vector<Double_t> > vec_2D_peak_points_in)
{
    //cout << "Class_2D_crater_finder::add_data" << endl;
    vec_2D_peak_points = vec_2D_peak_points_in; // [x,y,height,density]
}

void Class_2D_crater_finder::set_finding_parameters(Double_t search_width_in,Double_t min_radius_in,Double_t max_radius_in,
                                                    Double_t search_step_size_in,
                                                    Double_t min_points_per_area_in, Double_t exclude_width_in,
                                                    Double_t x_center_in,Double_t y_center_in,Int_t N_maxima_in)
{
    //cout << "Class_2D_crater_finder::set_finding_parameters" << endl;
    search_width        = search_width_in;
    min_radius          = min_radius_in;
    max_radius          = max_radius_in;
    search_step_size    = search_step_size_in;
    x_center            = x_center_in;
    y_center            = y_center_in;
    min_points_per_area = min_points_per_area_in;
    exclude_width       = exclude_width_in;
    N_maxima            = N_maxima_in;

    ring_quality        = 0.0;
    ring_points         = 0.0;
    param_radii         = 12;
}


void Class_2D_crater_finder::order_radii()
{
    //------------------
    // Order radii from inner to outer, in vec_radii_peak_positions_quality they are ordered by quality
    // This is important to calculate the background in between the rings
    vec_radii_peak_positions_radius.clear();
    vec_radii_peak_positions_radius.resize(param_radii);
    vec_radii_partial_ring_angles_radius.clear();
    index_quality_to_radius.clear();
    index_radius_to_quality.clear();

    vector<Int_t> flag_used_smallest_radius;
    flag_used_smallest_radius.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    index_quality_to_radius.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    index_radius_to_quality.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)flag_used_smallest_radius.size(); i_2D_radius++)
    {
        flag_used_smallest_radius[i_2D_radius] = 0;
    }
    Int_t N_used_radii = 0;
    while(N_used_radii < (Int_t)vec_radii_peak_positions_quality[0].size())
    {
        Double_t smallest_radius       = 900000000.0;
        Int_t    index_smallest_radius = 0;
        for(Int_t i_2D_radius = 0; i_2D_radius < vec_radii_peak_positions_quality[0].size(); i_2D_radius++)
        {
            Double_t radius = vec_radii_peak_positions_quality[0][i_2D_radius];
            if(radius < smallest_radius && !flag_used_smallest_radius[i_2D_radius])
            {
                smallest_radius       = radius;
                index_smallest_radius = i_2D_radius;
            }
        }
        flag_used_smallest_radius[index_smallest_radius] = 1;
        for(Int_t i_param = 0; i_param < param_radii; i_param++)
        {
            vec_radii_peak_positions_radius[i_param].push_back(vec_radii_peak_positions_quality[i_param][index_smallest_radius]);
        }
        Int_t index_radius = vec_radii_peak_positions_radius[0].size() - 1;
        index_quality_to_radius[index_smallest_radius] = index_radius; // input is index of best quality ordered array, output is index of smallest radius ordered array
        index_radius_to_quality[index_radius]          = index_smallest_radius;
        vec_radii_partial_ring_angles_radius.push_back(vec_radii_partial_ring_angles[index_smallest_radius]);
        N_used_radii++;
    }
    //------------------
}


void Class_2D_crater_finder::find_craters()
{
    //cout << "Class_2D_crater_finder::find_craters" << endl;
    // Find maximum radius
    Double_t max_radius_points = 0.0;
    for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
    {
        Double_t x_val  = vec_2D_peak_points[i_point][0] - x_center;
        Double_t y_val  = vec_2D_peak_points[i_point][1] - y_center;
        Double_t height = vec_2D_peak_points[i_point][2];

        Double_t point_radius = TMath::Sqrt(x_val*x_val + y_val*y_val);
        if(point_radius > max_radius_points) max_radius_points = point_radius;
    }

    vector<Double_t> vec_angle_partial_ring_start_stop;
    vec_angle_partial_ring_start_stop.resize(2); // start, stop

    //----------------------------------------------------------------------
    // Search in rings within inner radius and outer radius
    vector<Double_t> vec_quality;
    vector< vector< vector<Double_t> > > vec_partial_rings_angle;
    vector<Double_t> vec_points_partial_ring;
    vector<Double_t> vec_density;
    vector<Double_t> vec_radius;
    vector<Double_t> vec_half_width;
    vector<Int_t>    vec_N_points;
    vector<Double_t> vec_RMS_height;
    vector<Double_t> vec_RMS_width;
    vector<Double_t> vec_mean_height;
    vector<Double_t> vec_mean_width;
    vector<Double_t> vec_max_height;
    vector<Double_t> vec_min_height;
    vector<Double_t> vec_ext_density;
    Int_t i_radius = 0;
    Double_t radius_low  = i_radius*search_step_size;
    while(radius_low < max_radius_points && radius_low < max_radius) // scan radius
    {
        radius_low  = i_radius*search_step_size + min_radius;
        Double_t radius_high = radius_low + search_width/2.0;
        Double_t radius_mean = (radius_high + radius_low)/2.0;
        //Double_t area        = TMath::Pi()*(radius_high*radius_high - radius_low*radius_low);

        //printf("i_radius: %d, radius_mean: %f, min_radius: %f \n",i_radius,radius_mean,min_radius);

        // Determine search width with maximum density
        Double_t max_quality     = 0.0;
        Double_t max_density     = 0.0;
        Double_t max_half_width  = 0.0;
        Double_t max_RMS_height  = 0.0;
        Double_t max_RMS_width   = 0.0;
        Int_t max_N_points       = 0;
        Double_t max_mean_height = 0.0;
        Double_t max_mean_width  = 0.0;
        Double_t max_max_height  = 0.0;
        Double_t max_min_height  = 0.0;
        for(Int_t i_search = 0; i_search < 20; i_search++) // scan width
        {
            Double_t search_width_half_scan = search_width/2.0 + (search_width/2.0)*i_search;

            Double_t radius_low_scan  = radius_mean - search_width_half_scan;
            Double_t radius_high_scan = radius_mean + search_width_half_scan;
            Double_t area             = 2.0*search_width_half_scan;

            if(area <= 0.0) continue;
            if(radius_low_scan <= 0.0) continue;

            Double_t N_points_in_range  = 0.0;
            Double_t RMS_height         = 10000000.0;
            Double_t RMS_width          = 10000000.0;
            Double_t mean_height        = 0.0;
            Double_t mean_width         = 0.0;
            Double_t max_height         = 0.0;
            Double_t min_height         = 0.0;
            Double_t full_average_width_from_peaks = 0.0;

            Double_t sum_weights_itt2 = 0.0;
            for(Int_t i_itt = 0; i_itt < 3; i_itt++) // one itteration truncated mean, second itteration is with new width
            {
                // i_itt == 0 -> get mean and sigma from peak width
                // i_itt == 1 -> apply mean and sigma from peak width from i_itt0 and get truncated mean and sigma
                // i_itt == 2 -> apply truncated mean and sigma from i_itt1, modify search range, cut peaks within new search range
                if(i_itt == 2) // get new width from point width RMS
                {
                    full_average_width_from_peaks = 0.7*mean_width + 0.0*RMS_width; // 0.5
                    radius_low_scan  = radius_mean - full_average_width_from_peaks/2.0;
                    radius_high_scan = radius_mean + full_average_width_from_peaks/2.0;
                    area             = full_average_width_from_peaks;

                    if(area <= 0.0) area = 10000000.0;
                    if(radius_low_scan <= 0.0) radius_low_scan = 0.0;
                }

                N_points_in_range           = 0.0;
                Double_t sum_weights        = 0.0;
                Double_t sum_height         = 0.0;
                Double_t sum_heightsq       = 0.0;
                Double_t sum_width          = 0.0;
                Double_t sum_widthsq        = 0.0;
                Double_t sum_weight_rho_sq  = 0.0;
                
                for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
                {
                    Double_t x_val        = vec_2D_peak_points[i_point][0] - x_center;
                    Double_t y_val        = vec_2D_peak_points[i_point][1] - y_center;
                    Double_t height       = vec_2D_peak_points[i_point][2];
                    Double_t density_peak = vec_2D_peak_points[i_point][3]; // from 1D peak finder
                    Double_t width_peak   = vec_2D_peak_points[i_point][4]; // from 1D peak finder, full width

                    Double_t point_radius = TMath::Sqrt(x_val*x_val + y_val*y_val);
                    Double_t rho_sq       = TMath::Power(point_radius - radius_mean,2);

                    Double_t weight       = 1.0;
                    //Double_t weight       = density_peak;

                    sum_weight_rho_sq  += weight*rho_sq;

                    if(point_radius >= radius_low_scan
                       && point_radius <= radius_high_scan
                       //&& fabs(height - mean_height) < 2.0*RMS_height
                       && (fabs(width_peak - mean_width) < 2.0*RMS_width)
                      )
                    {
                        if(N_points_in_range <= 0.0)
                        {
                            max_height = height;
                            min_height = height;
                        }
                        else
                        {
                            if(height > max_height) max_height = height;
                            if(height < min_height) min_height = height;
                        }
                        N_points_in_range += 1.0;
                        sum_weights       += weight;
                        sum_height        += weight*height;
                        sum_heightsq      += weight*height*height;
                        sum_width         += weight*width_peak;
                        sum_widthsq       += weight*width_peak*width_peak;
                        Double_t average_local_width = -1.0;
                        if(N_points_in_range > 0.0) average_local_width = sum_width/sum_weights;
                        //if(radius_mean == 58.25 && i_search == 0 && i_itt == 2) printf("i_search: %d, i_itt: %d, width_peak: %f, pos: {%f,%f}, average_local_width: %f, points_in_range: %f \n",i_search,i_itt,width_peak,x_val,y_val,average_local_width,sum_weights);
                    }
                } // end of loop over points

                if(i_itt == 2) sum_weights_itt2 = sum_weights;
                if(sum_weights > 0.0)
                {
                    mean_height = sum_height/sum_weights;
                    mean_width  = sum_width/sum_weights;
                    RMS_height  = TMath::Sqrt((1.0/sum_weights)*(sum_heightsq - 2.0*mean_height*sum_height + sum_weights*mean_height*mean_height));
                    RMS_width   = TMath::Sqrt((1.0/sum_weights)*(sum_widthsq - 2.0*mean_width*sum_width + sum_weights*mean_width*mean_width));
                    //if(radius_mean == 58.25 && i_search == 0 && i_itt == 2) printf("  i_search: %d, i_itt: %d, mean_width: %f, RMS_width: %f \n",i_search,i_itt,mean_width,RMS_width);
                }
                //printf("radius_mean: %f, i_search: %d, search_width_half_scan: %f i_itt: %d, mean_peak_width: %f, RMS_peak_width: %f \n",radius_mean,i_search,search_width_half_scan,i_itt,mean_width,RMS_width);
            } // end of truncation

            if(sum_weights_itt2 <= 0.0) continue;
            Double_t density = sum_weights_itt2;
            Double_t quality = -1.0;
            quality = density;
            //printf("   i_search: %d, density: %f, area: %f, sum_weights_itt2: %f, search_width_half_scan: %f, RMS_height: %f \n",i_search,density,area,sum_weights_itt2,search_width_half_scan,RMS_height);
            if(quality > max_quality)
            {
                max_quality     = quality;
                max_density     = density;
                //max_half_width  = search_width_half_scan;
                max_half_width  = full_average_width_from_peaks/2.0;
                max_N_points    = N_points_in_range;
                max_RMS_height  = RMS_height;
                max_RMS_width   = RMS_width;
                max_mean_height = mean_height;
                max_mean_width  = mean_width;
                max_max_height  = max_height;
                max_min_height  = min_height;
                //printf("      -> radius_mean: %f, i_search: %d, max_density: %f, max_half_width: %f, RMS_width: %f \n",radius_mean,i_search,max_density,max_half_width,RMS_width);
            }
        } // end of scan width


        //-----------------------------------
        // Partial ring detection
        //printf("radius: %f, half_width: %f \n",radius_mean,max_half_width);
        // Make new point vector for points within detected ring, add angle information
        vector< vector<Double_t> > vec_points_ring;
        vec_points_ring.resize(4); // x,y,radius,angle
        Double_t radius_low_scan  = radius_mean - max_half_width;
        Double_t radius_high_scan = radius_mean + max_half_width;
        for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
        {
            Double_t x_val        = vec_2D_peak_points[i_point][0] - x_center;
            Double_t y_val        = vec_2D_peak_points[i_point][1] - y_center;
            Double_t point_radius = TMath::Sqrt(x_val*x_val + y_val*y_val);
            Double_t angle        = TMath::ATan2(y_val,x_val)*TMath::RadToDeg(); // angle in degree
            if(angle < 0.0) angle += 360.0;

            if(point_radius >= radius_low_scan && point_radius <= radius_high_scan)
            {
                vec_points_ring[0].push_back(x_val);
                vec_points_ring[1].push_back(y_val);
                vec_points_ring[2].push_back(point_radius);
                vec_points_ring[3].push_back(angle);
            }

        }

        //------------------
        // Order vec_points_ring from lowest to highest angle
        vector< vector<Double_t> > vec_points_ring_ordered_angle;
        vec_points_ring_ordered_angle.resize(4); // x,y,radius,angle
        vector<Int_t> flag_used_smallest_angle;
        flag_used_smallest_angle.resize((Int_t)vec_points_ring[0].size());
        for(Int_t i_point = 0; i_point < vec_points_ring[0].size(); i_point++)
        {
            flag_used_smallest_angle[i_point] = 0;
        }
        Int_t N_used_angle = 0;
        while(N_used_angle < (Int_t)vec_points_ring[0].size())
        {
            Double_t smallest_angle       = 900000000.0;
            Int_t    index_smallest_angle = 0;
            for(Int_t i_point = 0; i_point < vec_points_ring[0].size(); i_point++)
            {
                Double_t angle = vec_points_ring[3][i_point]; // angle in degree
                if(angle < smallest_angle && !flag_used_smallest_angle[i_point])
                {
                    smallest_angle       = angle;
                    index_smallest_angle = i_point;
                }
            }
            flag_used_smallest_angle[index_smallest_angle] = 1;
            for(Int_t i_param = 0; i_param < 4; i_param++)
            {
                vec_points_ring_ordered_angle[i_param].push_back(vec_points_ring[i_param][index_smallest_angle]);
            }
            N_used_angle++;
        }
        vec_points_ring.clear();
        //------------------


        //------------------
        // Loop over all points which are ordered from lowest to highest angle
        // Take always the next Delta_points_angle points and check the angular density: points/Delta_angle
        // If the angular density is above the threshold then take the next Delta_points_angle points
        Double_t N_points_per_angle_threshold = 0.7; // 0.7
        Int_t Delta_points_angle = 10;
        Int_t N_points_angle = vec_points_ring_ordered_angle[0].size();
        Int_t i_point_start = 0;
        Double_t Sum_angle_partial_rings  = 0.0;
        Int_t    Sum_points_partial_rings = 0;
        vector< vector<Double_t> > vec_array_angle_start_stop;
        while(i_point_start < N_points_angle)
        {
            Double_t angle_start = vec_points_ring_ordered_angle[3][i_point_start]; // angle in degree
            Double_t angle_start_partial_ring = angle_start;
            Int_t i_point_stop = i_point_start + Delta_points_angle;
            Int_t flag_break = 0;
            Int_t N_total_points_in_partial_ring = 0;
            Double_t angle_stop = angle_start;
            Double_t previous_angle_stop = angle_stop;
            Double_t N_points_per_angle = 0.0;
            while(i_point_stop < N_points_angle)
            {
                if(i_point_stop >= N_points_angle)
                {
                    i_point_stop = (N_points_angle - 1);
                    flag_break = 1;
                }
                angle_stop  = vec_points_ring_ordered_angle[3][i_point_stop]; // angle in degree
                Double_t Delta_angle = fabs(angle_stop - angle_start);
                if(Delta_angle > 360.0) Delta_angle -= 360.0; // angle integration range can be up to 360 degree
                N_points_per_angle = 0.0;
                if(Delta_angle > 0.0)
                {
                    N_points_per_angle = ((Double_t)Delta_points_angle)/Delta_angle;
                }
                N_total_points_in_partial_ring = (i_point_stop - i_point_start) + 1;
                angle_start = angle_stop;
                if(N_points_per_angle < N_points_per_angle_threshold) flag_break = 1;
                if(flag_break) break;
                i_point_stop += Delta_points_angle;
                previous_angle_stop = angle_stop;
            }
            i_point_start = i_point_stop + 1;
            Double_t Delta_angle_partial_ring = fabs(previous_angle_stop - angle_start_partial_ring);
            Double_t Total_points_per_angle = 0.0;
            if(Delta_angle_partial_ring > 0.0)
            {
                Total_points_per_angle = ((Double_t)N_total_points_in_partial_ring)/Delta_angle_partial_ring;
            }
            if(Delta_angle_partial_ring > 20.0 && Total_points_per_angle > N_points_per_angle_threshold) // 20.0, 1.0
            {
                Sum_angle_partial_rings  += Delta_angle_partial_ring;
                Sum_points_partial_rings += N_total_points_in_partial_ring;
                vec_angle_partial_ring_start_stop[0] = angle_start_partial_ring;
                vec_angle_partial_ring_start_stop[1] = previous_angle_stop;
                vec_array_angle_start_stop.push_back(vec_angle_partial_ring_start_stop);
                //printf("   ->partial ring: {%f,%f}, points: %d, Total_points_per_angle: %f \n",angle_start_partial_ring,previous_angle_stop,N_total_points_in_partial_ring,Total_points_per_angle);
            }
        }
        Double_t Extrapolated_points_full_ring = 0.0;
        if(Sum_angle_partial_rings > 90.0) // 60.0
        {
            Extrapolated_points_full_ring = 360.0*((Double_t)Sum_points_partial_rings)/Sum_angle_partial_rings;
        }
        //Extrapolated_points_full_ring /= (2.0*max_half_width);
        //printf("Sum_angle_partial_rings: %f, Sum_points_partial_rings: %d, Extrapolated_points_full_ring: %f \n",Sum_angle_partial_rings,Sum_points_partial_rings,Extrapolated_points_full_ring);
        //cout << "" << endl;
        //------------------
        // End of partical ring detection
        //-----------------------------------


        //-----------------------------------
        // Store optimal values for each ring
        //vec_quality.push_back(max_quality);
        //Double_t used_quality = Extrapolated_points_full_ring*(Sum_angle_partial_rings/360.0)/(2.0*max_half_width);
        Double_t used_quality = (Double_t)Sum_points_partial_rings/(2.0*max_half_width);
        //Double_t used_quality = (Double_t)Sum_points_partial_rings;
        //Double_t used_quality = Extrapolated_points_full_ring/(2.0*max_half_width);
        //Double_t used_quality = Extrapolated_points_full_ring;  // NORMX
        vec_quality.push_back(used_quality);
        vec_partial_rings_angle.push_back(vec_array_angle_start_stop);
        vec_points_partial_ring.push_back(Sum_points_partial_rings);
        vec_density.push_back(max_density);
        vec_radius.push_back(radius_mean);
        vec_half_width.push_back(max_half_width);
        vec_N_points.push_back(max_N_points);
        vec_RMS_height.push_back(max_RMS_height);
        vec_RMS_width.push_back(max_RMS_width);
        vec_mean_height.push_back(max_mean_height);
        vec_mean_width.push_back(max_mean_width);
        vec_max_height.push_back(max_max_height);
        vec_min_height.push_back(max_min_height);
        vec_ext_density.push_back(Extrapolated_points_full_ring);

        vec_array_angle_start_stop.clear();
        //-----------------------------------

        i_radius++;
    } // end of radius scan
    //----------------------------------------------------------------------


    //----------------------------------------------------------------------
    // determine min and max radii
    Double_t min_max_radii[2];
    for(Int_t i_point = 0; i_point < vec_radius.size(); i_point++)
    {
        if(i_point == 0)
        {
            min_max_radii[0] = vec_radius[i_point];
            min_max_radii[1] = vec_radius[i_point];
        }
        else
        {
            if(vec_radius[i_point] < min_max_radii[0]) min_max_radii[0] = vec_radius[i_point];
            if(vec_radius[i_point] > min_max_radii[1]) min_max_radii[1] = vec_radius[i_point];
        }
    }
    //----------------------------------------------------------------------


    //----------------------------------------------------------------------
    // Identify highest densities
    Int_t N_max_densities = N_maxima;
    vec_radii_peak_positions_quality.clear();
    vec_radii_peak_positions_quality.resize(param_radii);
    vec_radii_peak_positions_radius.clear();
    vec_radii_peak_positions_radius.resize(param_radii);
    vec_radii_partial_ring_angles.clear();
    vec_radii_partial_ring_angles_radius.clear();
    Int_t N_max_densities_found = 0;
    vector<Double_t> vec_exclude_low;
    vector<Double_t> vec_exclude_high;
    vec_exclude_low.push_back(-1.0);
    vec_exclude_high.push_back(-1.0);
    while(N_max_densities_found < N_max_densities)
    {
        if(N_max_densities_found > 10) break;
        Double_t best_radius      = 0.0;
        Double_t best_quality     = 0.0;
        Double_t best_density     = 0.0;
        Double_t best_half_width  = 0.0;
        Int_t    best_point       = -1;
        Int_t    best_N_points    = 0;
        Int_t    best_N_points_PR = 0;
        Double_t best_RMS_height  = 0.0;
        Double_t best_RMS_width   = 0.0;
        Double_t best_mean_height = 0.0;
        Double_t best_mean_width  = 0.0;
        Double_t best_max_height  = 0.0;
        Double_t best_min_height  = 0.0;
        Double_t best_ext_density = 0.0;
        for(Int_t i_point = 0; i_point < vec_radius.size(); i_point++)
        {
            Double_t radius          = vec_radius[i_point];
            Double_t half_width      = vec_half_width[i_point];
            Double_t quality         = vec_quality[i_point];
            Double_t density         = vec_density[i_point];
            Int_t    N_points_PR     = vec_points_partial_ring[i_point];
            Int_t    N_points        = vec_N_points[i_point];
            Double_t RMS_height      = vec_RMS_height[i_point];
            Double_t RMS_width       = vec_RMS_width[i_point];
            Double_t mean_height     = vec_mean_height[i_point];
            Double_t mean_width      = vec_mean_width[i_point];
            Double_t max_height      = vec_max_height[i_point];
            Double_t min_height      = vec_min_height[i_point];
            Double_t ext_density     = vec_ext_density[i_point];

            //if(N_max_densities_found == 0) printf("i_point: %d, radius: %f, half_widht: %f, density: %f, N_points: %d, RMS_height: %f \n",i_point,radius,half_width,density,N_points,RMS_height);

            Int_t flag_exclude = 0;
            for(Int_t i_exclude = 0; i_exclude < vec_exclude_high.size(); i_exclude++)
            {
                if((radius + half_width) >= vec_exclude_low[i_exclude] && (radius-half_width) <= vec_exclude_high[i_exclude])
                {
                    flag_exclude = 1;
                    break;
                }
            }

            //if(density > best_density && !flag_exclude)
            if(quality > best_quality && !flag_exclude)
            {
                best_radius      = radius;
                best_quality     = quality;
                best_density     = density;
                best_half_width  = half_width;
                best_point       = i_point;
                best_N_points    = N_points;
                best_N_points_PR = N_points_PR;
                best_RMS_height  = RMS_height;
                best_RMS_width   = RMS_width;
                best_mean_height = mean_height;
                best_mean_width  = mean_width;
                best_max_height  = max_height;
                best_min_height  = min_height;
                //printf("best_density: %f \n",best_density);
            }
        }

        //printf("best_radius: %f, best_quality: %f, best_density: %f, best_half_width: %f, best_N_points: %d, best_RMS_height: %f, best_RMS_width: %f \n",best_radius,best_quality,best_density,best_half_width,best_N_points,best_RMS_height,best_RMS_width);
        if(best_half_width < exclude_width) best_half_width = exclude_width;
        if(best_point < 0) break;
        vec_radii_peak_positions_quality[0].push_back(vec_radius[best_point]);
        vec_radii_peak_positions_quality[1].push_back(vec_density[best_point]);
        vec_radii_peak_positions_quality[2].push_back(vec_RMS_height[best_point]);
        vec_radii_peak_positions_quality[3].push_back(vec_half_width[best_point]);
        vec_radii_peak_positions_quality[4].push_back(vec_RMS_width[best_point]);
        vec_radii_peak_positions_quality[5].push_back(vec_quality[best_point]);
        vec_radii_peak_positions_quality[6].push_back(vec_mean_height[best_point]);
        vec_radii_peak_positions_quality[7].push_back(vec_mean_width[best_point]);
        vec_radii_peak_positions_quality[8].push_back(vec_points_partial_ring[best_point]);
        vec_radii_peak_positions_quality[9].push_back(vec_max_height[best_point]);
        vec_radii_peak_positions_quality[10].push_back(vec_min_height[best_point]);
        vec_radii_peak_positions_quality[11].push_back(vec_ext_density[best_point]);
        vec_radii_partial_ring_angles.push_back(vec_partial_rings_angle[best_point]);
        vec_exclude_low.push_back(best_radius-1.0*best_half_width);
        vec_exclude_high.push_back(best_radius+1.0*best_half_width);
        N_max_densities_found++;
    }
    //----------------------------------------------------------------------

    order_radii();

    //cout << "" << endl;
    //cout << "-------------------------------" << endl;
    //for(Int_t i_2D_radius = 0; i_2D_radius < vec_radii_peak_positions_radius[0].size(); i_2D_radius++)
    //{
    //    printf("index radius, ordered by radius: %d, radius: %f, N_points: %f \n",i_2D_radius,vec_radii_peak_positions_radius[0][i_2D_radius],vec_radii_peak_positions_radius[1][i_2D_radius]);
    //}
    //cout << "-------------------------------" << endl;
    //cout << "" << endl;
    //------------------


    //------------------
    // Calculate number of peaks in between rings -> should be almost empty for good craters
    // Here the full space between the rings is used
    vector< vector<Double_t> > vec_points_in_between;
    vec_points_in_between.resize((Int_t)vec_radii_peak_positions_radius[0].size());
    vector< vector<Double_t> > vec_quality_in_between;
    vec_quality_in_between.resize((Int_t)vec_radii_peak_positions_radius[0].size());


    // example
    // [0,1,2,3,4,5,6] // radius
    // [4,1,6,0,2,3,5] // quality

    // [0,1,2,3,4,5,6] // quality
    // [3,1,4,5,0,6,2] // radius

    if(flag_debug == 2)
    {
        cout << " " << endl;
        for(Int_t i_index = 0; i_index < index_radius_to_quality.size(); i_index++)
        {
            printf("i_index: %d, index_radius_to_quality: %d \n",i_index,index_radius_to_quality[i_index]);
        }
        cout << " " << endl;
        for(Int_t i_index = 0; i_index < index_quality_to_radius.size(); i_index++)
        {
            printf("i_index: %d, index_quality_to_radius: %d \n",i_index,index_quality_to_radius[i_index]);
        }
        cout << " " << endl;
    }



    // Consecutive sum over radii. Use the first (quality) N radii and get the corresponding background.
    for(Int_t i_2D_radius_qual_cons = 0; i_2D_radius_qual_cons < (Int_t)vec_radii_peak_positions_quality[0].size(); i_2D_radius_qual_cons++)
    {
        if(flag_debug == 2) cout << "" << endl;
        if(flag_debug == 2) printf("i_2D_radius_qual_cons: %d \n",i_2D_radius_qual_cons);
        Int_t N_radii_used = 0;
        for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_radius[0].size()-1; i_2D_radius++)
        {
            // Get first index according to radius
            Int_t index_quality = index_radius_to_quality[i_2D_radius];
            if(index_quality > i_2D_radius_qual_cons) continue; // Use only radii for the first i_2D_radius_qual_cons quality radii
            Int_t i_2D_radiusA = i_2D_radius; // radii ordered
            if(flag_debug == 2) printf(" i_2D_radius: %d, index_quality: %d, i_2D_radius_qual_cons: %d \n",i_2D_radius,index_quality,i_2D_radius_qual_cons);

            // Get second index according to radius
            Int_t i_2D_radiusB    = -1;
            Int_t index_quality_X = -1;
            for(Int_t i_2D_radius_X = i_2D_radius+1; i_2D_radius_X < (Int_t)vec_radii_peak_positions_radius[0].size(); i_2D_radius_X++)
            {
                index_quality_X = index_radius_to_quality[i_2D_radius_X];
                if(index_quality_X > i_2D_radius_qual_cons)
                {
                    index_quality_X = -1;
                    continue; // Use only radii for the first i_2D_radius_qual_cons quality radii
                }
                else
                {
                    i_2D_radiusB = i_2D_radius_X; // radii ordered
                    break;
                }
            }

            if(flag_debug == 2) printf("  i_2D_radiusA/B (radii ordered): {%d,%d}, index_quality: {%d,%d} \n",i_2D_radiusA,i_2D_radiusB,index_quality,index_quality_X);

            //-------------------------------------
            // Get points from center to first ring
            Double_t halfwidthA     = -1.0;
            Double_t halfwidthB     = -1.0;
            Double_t radiusA        = -1.0;
            Double_t radiusB        = -1.0;
            Double_t N_points_in_between = 0.0;

            if(N_radii_used == 0)
            {
                halfwidthA  = vec_radii_peak_positions_radius[3][i_2D_radiusA];
                radiusA     = vec_radii_peak_positions_radius[0][i_2D_radiusA] - halfwidthA;
                for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
                {
                    Double_t x_val        = vec_2D_peak_points[i_point][0] - x_center;
                    Double_t y_val        = vec_2D_peak_points[i_point][1] - y_center;
                    Double_t point_radius = TMath::Sqrt(x_val*x_val + y_val*y_val);

                    if(point_radius <= radiusA)
                    {
                        N_points_in_between += 1.0;
                    }
                }
                Double_t total_width_in_between = radiusA;
                if(total_width_in_between <= 1.0) total_width_in_between = 100000.0; // Ignore this background ring, its too small

                if(flag_debug == 2) printf("    First ring: radii+/-width: {%f,%f}, widths: {%f,%f} \n",0.0,radiusA,0.0,halfwidthA);
                vec_points_in_between[i_2D_radius_qual_cons].push_back(N_points_in_between);
                vec_quality_in_between[i_2D_radius_qual_cons].push_back(N_points_in_between/total_width_in_between); // NORMX
            }
            //-------------------------------------



            //-------------------------------------
            // Get points in between rings
            if(i_2D_radiusB <= 0) continue;
            halfwidthA  = vec_radii_peak_positions_radius[3][i_2D_radiusA];
            radiusA     = vec_radii_peak_positions_radius[0][i_2D_radiusA] + halfwidthA;
            halfwidthB  = vec_radii_peak_positions_radius[3][i_2D_radiusB];
            radiusB     = vec_radii_peak_positions_radius[0][i_2D_radiusB] - halfwidthB;

            if(flag_debug == 2) printf("    radii: {%f,%f}, radii+/-width: {%f,%f}, widths: {%f,%f} \n",vec_radii_peak_positions_radius[0][i_2D_radiusA],vec_radii_peak_positions_radius[0][i_2D_radiusB],radiusA,radiusB,halfwidthA,halfwidthB);
            N_points_in_between = 0.0;
            for(Int_t i_point = 0; i_point < vec_2D_peak_points.size(); i_point++)
            {
                Double_t x_val        = vec_2D_peak_points[i_point][0] - x_center;
                Double_t y_val        = vec_2D_peak_points[i_point][1] - y_center;
                Double_t point_radius = TMath::Sqrt(x_val*x_val + y_val*y_val);

                if(point_radius >= radiusA && point_radius <= radiusB)
                {
                    N_points_in_between += 1.0;
                }
            }
            Double_t total_width_in_between = radiusB - radiusA;
            if(total_width_in_between <= 1.0) total_width_in_between = 100000.0; // Ignore this background ring, its too small

            vec_points_in_between[i_2D_radius_qual_cons].push_back(N_points_in_between);
            vec_quality_in_between[i_2D_radius_qual_cons].push_back(N_points_in_between/total_width_in_between); // NORMX
            //-------------------------------------

            Int_t index_vec_quality = (Int_t)vec_quality_in_between[i_2D_radius_qual_cons].size() - 1;
            if(flag_debug == 2) printf("--> radii: {%f,%f}, N_points: %f, quality: %f, width: %f \n",radiusA,radiusB,N_points_in_between,vec_quality_in_between[i_2D_radius_qual_cons][index_vec_quality],total_width_in_between);
            N_radii_used++;
        } // End of ring loop
    } // End of consecutive points in between loop

    //cout << "" << endl;
    //cout << "-------------------------------" << endl;
    //for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_points_in_between.size(); i_2D_radius++)
    //{
    //    printf("N_points_in_between: %f, quality in between: %f \n",vec_points_in_between[i_2D_radius],vec_quality_in_between[i_2D_radius]);
    //}
    //cout << "-------------------------------" << endl;
    //cout << "" << endl;
    //------------------



    //------------------
    cout << "" << endl;
    cout << "-------------------------------" << endl;
    vector<Double_t> sum_ring_quality;
    vector<Double_t> sum_ring_points;
    sum_ring_quality.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    sum_ring_points.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    for(Int_t i_2D_radius = 0; i_2D_radius < vec_radii_peak_positions_quality[0].size(); i_2D_radius++)
    {
        cout << "i_2D_radius: " << i_2D_radius << ", radius: " << vec_radii_peak_positions_quality[0][i_2D_radius]
            << ", density: : " << vec_radii_peak_positions_quality[1][i_2D_radius]
            << ", quality: : " << vec_radii_peak_positions_quality[5][i_2D_radius] << ", RMS_height: " << vec_radii_peak_positions_quality[2][i_2D_radius]
            << ", RMS_width: " << vec_radii_peak_positions_quality[4][i_2D_radius]
            << ", <height>: " << vec_radii_peak_positions_quality[6][i_2D_radius]
            << ", half_width: " << vec_radii_peak_positions_quality[3][i_2D_radius]
            << ", N_points_PR: " << vec_radii_peak_positions_quality[8][i_2D_radius]
            << ", max_height: " << vec_radii_peak_positions_quality[9][i_2D_radius]
            << ", min_height: " << vec_radii_peak_positions_quality[10][i_2D_radius]
            << ", N_points_ext: " << vec_radii_peak_positions_quality[11][i_2D_radius] << endl;

        if(vec_radii_peak_positions_quality[2][i_2D_radius] > 0.0)
        {
            //ring_quality += vec_radii_peak_positions_quality[1][i_2D_radius]/vec_radii_peak_positions_quality[2][i_2D_radius]; // density/RMS_height
            //ring_quality += vec_radii_peak_positions_quality[1][i_2D_radius];
            if(i_2D_radius < 100)
            {
                ring_quality += vec_radii_peak_positions_quality[5][i_2D_radius];
                sum_ring_quality[i_2D_radius] = ring_quality; // consecutive sum of ring qualities
                ring_points += vec_radii_peak_positions_quality[1][i_2D_radius];
                sum_ring_points[i_2D_radius] = ring_points; // consecutive sum of ring points
            }
        }
    }
    //for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_points_in_between.size(); i_2D_radius++)
    //{
    //    ring_quality -= vec_points_in_between[i_2D_radius];
    //}
    for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_quality_in_between.size(); i_2D_radius++)
    {
        //ring_quality -= vec_quality_in_between[i_2D_radius];
    }
    cout << "-------------------------------" << endl;
    cout << "" << endl;

    // Loop over consecutive sums of backgroud and signal and determine best number of rings
    Int_t    N_radii_max_quality = -1;
    Double_t max_quality_radii   = -1000000.0;
    Double_t max_significance    = -1000000.0;
    for(Int_t i_2D_radius = 0; i_2D_radius < vec_radii_peak_positions_quality[0].size(); i_2D_radius++) // Number of radii used -> consecutive sums
    {
        // Calculate sum of width for all rings used in this loop
        Double_t sum_ring_widths     = 0.0;
        Double_t max_radius_used     = 0.0;
        Double_t max_half_width_used = 0.0;
        for(Int_t i_width_sum = 0; i_width_sum < (i_2D_radius+1); i_width_sum++)
        {
            Double_t radius_used     = vec_radii_peak_positions_quality[0][i_width_sum];
            Double_t half_width_used = vec_radii_peak_positions_quality[3][i_width_sum];
            sum_ring_widths += 2.0*half_width_used;
            if(radius_used > max_radius_used)
            {
                max_radius_used      = radius_used;
                max_half_width_used  = half_width_used;
            }
        }
        Double_t sum_bkgr_widths = max_radius_used + max_half_width_used - sum_ring_widths;

        vec_radii_peak_positions_quality[0][i_2D_radius];

        Double_t bkgr_sum        = 0.0;
        Double_t bkgr_sum_points = 0.0;
        for(Int_t i_bkgr_sum = 0; i_bkgr_sum < vec_quality_in_between[i_2D_radius].size(); i_bkgr_sum++)
        {
            bkgr_sum        += vec_quality_in_between[i_2D_radius][i_bkgr_sum];
            bkgr_sum_points += vec_points_in_between[i_2D_radius][i_bkgr_sum];
        }
        Double_t this_ring_quality = sum_ring_quality[i_2D_radius];
        Double_t this_ring_points  = sum_ring_points[i_2D_radius];
        if(i_2D_radius > 0)
        {
            this_ring_quality = sum_ring_quality[i_2D_radius] - sum_ring_quality[i_2D_radius-1];
            this_ring_points  = sum_ring_points[i_2D_radius] - sum_ring_points[i_2D_radius-1];
        }
        Double_t ring_quality_consecutive = sum_ring_quality[i_2D_radius] - bkgr_sum;
        Double_t ring_significance = -1.0;
        if(sum_ring_widths < 0.0 || sum_bkgr_widths < 0.0) continue;
        //Double_t rel_ring_points = sum_ring_points[i_2D_radius]/sum_ring_widths;
        //Double_t rel_bkgr_points = bkgr_sum_points/sum_bkgr_widths;
        Double_t rel_ring_points = sum_ring_points[i_2D_radius];
        Double_t rel_bkgr_points = bkgr_sum_points;
        if((rel_ring_points + rel_bkgr_points) > 0.0)
        {
            ring_significance = rel_ring_points/TMath::Sqrt(rel_ring_points + rel_bkgr_points);
        }
        Double_t average_bkgr = bkgr_sum/((Double_t)(i_2D_radius + 1));
        //ring_quality_consecutive *= 1.0/(i_2D_radius + 10.0);
        Double_t ratio_average_signal_bkgr = 1000.0;
        if(average_bkgr > 0.0)
        {
            ratio_average_signal_bkgr = this_ring_quality/average_bkgr;
        }
        printf("sum_ring_widths: %f, qual.: %f, <S/B>: %f, sum_bkgr_widths: %f, ring points: %f, bkgr points: %f, rel_ring_points: %f, rel_bkgr_points: %f, ring sign.: %f \n",sum_ring_widths,ring_quality_consecutive,ratio_average_signal_bkgr,sum_bkgr_widths,sum_ring_points[i_2D_radius],bkgr_sum_points,rel_ring_points,rel_bkgr_points,ring_significance);
        //printf("N_radii used for quality: %d, ring sum: %f, background sum: %f, quality: %f, average_bkgr: %f, this_ring_quality: %f, ratio_average_signal_bkgr: %f, ring_significance: %f \n",i_2D_radius,sum_ring_quality[i_2D_radius],bkgr_sum,ring_quality_consecutive,average_bkgr,this_ring_quality,ratio_average_signal_bkgr,ring_significance);
        if(ring_significance > max_significance && ratio_average_signal_bkgr > 1.25)
        //if(ring_quality_consecutive > max_quality_radii && ratio_average_signal_bkgr > 1.25) // 2.5
        {
            max_quality_radii   = ring_quality_consecutive;
            max_significance    = ring_significance;
            N_radii_max_quality = i_2D_radius;
        }
    } // end of loop over consecutive sums


    N_radii_max_quality += 1; // i_2D_radius starts from 0, one good ring: 0+1 = 1

    cout << "" << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    printf("=========> N_radii_max_quality: %d \n",N_radii_max_quality);
    cout << "---------------------------------------------------------------------------" << endl;
    cout << "" << endl;

    Double_t cons_sum_points_ext = 0.0;
    vector<Int_t> vec_index_best_points_ext;
    vector<Int_t> vec_flag_best_points_ext;
    vec_index_best_points_ext.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    vec_flag_best_points_ext.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_quality[0].size(); i_2D_radius++) // Number of radii used -> consecutive sums
    {
        vec_flag_best_points_ext[i_2D_radius] = 0;
    }
    Int_t Index_sum = 0;
    while(Index_sum  < (Int_t)vec_radii_peak_positions_quality[0].size())
    {
        Double_t max_N_points_ext = 0.0;
        Int_t    max_index        = 0;
        for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_quality[0].size(); i_2D_radius++) // Number of radii used -> consecutive sums
        {
            Double_t radius     = vec_radii_peak_positions_quality[0][i_2D_radius];
            Double_t half_width = vec_radii_peak_positions_quality[3][i_2D_radius];
            if(!vec_flag_best_points_ext[i_2D_radius] && half_width > 0.0)
            {
                Double_t N_points_ext = vec_radii_peak_positions_quality[11][i_2D_radius]*TMath::Power(radius,0.3)/(2.0*half_width);
                if(N_points_ext >= max_N_points_ext)
                {
                    max_N_points_ext = N_points_ext;
                    max_index        = i_2D_radius;
                }
            }
        }
        vec_flag_best_points_ext[max_index]  = 1;
        vec_index_best_points_ext[Index_sum] = max_index;
        Index_sum++;
    }
    for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_quality[0].size(); i_2D_radius++) // Number of radii used -> consecutive sums
    {
        printf("index radius (quality): %d, radius: %f, N_points_ext: %f \n",i_2D_radius,vec_radii_peak_positions_quality[0][i_2D_radius],vec_radii_peak_positions_quality[11][i_2D_radius]);
    }
    cout << "" << endl;
    vector<Double_t> vec_cons_sum_N_points;
    //vec_cons_sum_N_points.resize((Int_t)vec_radii_peak_positions_quality[0].size());
    vec_cons_sum_N_points.resize(N_radii_max_quality);
    //for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_quality[0].size(); i_2D_radius++) // Number of radii used -> consecutive sums
    for(Int_t i_2D_radius = 0; i_2D_radius < N_radii_max_quality; i_2D_radius++) // Number of radii used -> consecutive sums
    {
        Int_t index_cons_points = vec_index_best_points_ext[i_2D_radius];
        Double_t radius     = vec_radii_peak_positions_quality[0][index_cons_points];
        Double_t half_width = vec_radii_peak_positions_quality[3][index_cons_points];
        printf("i_2D_radius: %d, radius: %f, half_width: %f, index_cons_points: %d \n",i_2D_radius,radius,half_width,index_cons_points);
        if(half_width > 0.0)
        {
            Double_t ext_quality = vec_radii_peak_positions_quality[11][index_cons_points]*TMath::Power(radius,0.3)/(2.0*half_width);
            cons_sum_points_ext += ext_quality;
            vec_cons_sum_N_points[i_2D_radius] = cons_sum_points_ext;
            printf("index radius (ext points): %d, index_cons_points: %d, radius: %f, ext_quality: %f, N_points_ext: %f, cons_sum_points_ext: %f \n",i_2D_radius,index_cons_points,vec_radii_peak_positions_quality[0][index_cons_points],ext_quality,vec_radii_peak_positions_quality[11][index_cons_points],cons_sum_points_ext);
        }
    }
    cout << "---------------------------------------------------------------------------" << endl;
    cout << "" << endl;

    //ring_quality = max_quality_radii;
    //ring_quality = max_significance;
    if((Int_t)vec_cons_sum_N_points.size() >= 4) ring_quality = vec_cons_sum_N_points[3];
    else
    {
        if((Int_t)vec_cons_sum_N_points.size() <= 0) ring_quality = 0.0;
        else ring_quality = vec_cons_sum_N_points[(Int_t)vec_cons_sum_N_points.size()-1];
    }
    printf("N_radii_max_quality: %d, max_quality_radii: %f, max_significance: %f, ring_quality: %f \n",N_radii_max_quality,max_quality_radii,max_significance,ring_quality);

    for(Int_t i_param = 0; i_param < param_radii; i_param++)
    {
        vec_radii_peak_positions_quality[i_param].resize(N_radii_max_quality);
        vec_radii_partial_ring_angles.resize(N_radii_max_quality);
    }
    order_radii();


    vec_radii_peak_positions_radius_merge      = vec_radii_peak_positions_radius;
    vec_radii_partial_ring_angles_merge        = vec_radii_partial_ring_angles_radius;
    vec_radii_peak_positions_radius_merge_copy = vec_radii_peak_positions_radius;
    vec_radii_partial_ring_angles_merge_copy   = vec_radii_partial_ring_angles_radius;
    //------------------

}

void Class_2D_crater_finder::Merge_close_rings(Double_t ring_dist_merge)
{
    // Merge rings which are seperated (radius + width) by less than ring_dist_merge
    //printf("Merge rings called \n");

    vec_radii_partial_ring_angles_merge.clear();
    vec_radii_peak_positions_radius_merge.clear();
    vec_radii_peak_positions_radius_merge.resize(param_radii);
    vector<Double_t> vec_angle_partial_ring_start_stop;
    vec_angle_partial_ring_start_stop.resize(2); // start, stop

    vector<Int_t> vec_flag_ring_merged;
    vec_flag_ring_merged.resize((Int_t)vec_radii_peak_positions_radius_merge_copy[0].size());
    for(Int_t i_2D_radiusA = 0; i_2D_radiusA < (Int_t)vec_radii_peak_positions_radius_merge_copy[0].size(); i_2D_radiusA++)
    {
        vec_flag_ring_merged[i_2D_radiusA] = 0;
    }
    for(Int_t i_2D_radiusA = 0; i_2D_radiusA < (Int_t)vec_radii_peak_positions_radius_merge_copy[0].size(); i_2D_radiusA++)
    {
        if(vec_flag_ring_merged[i_2D_radiusA]) continue;
        Double_t radiusA      = vec_radii_peak_positions_radius_merge_copy[0][i_2D_radiusA];
        Double_t densityA     = vec_radii_peak_positions_radius_merge_copy[1][i_2D_radiusA];
        Double_t RMS_heightA  = vec_radii_peak_positions_radius_merge_copy[2][i_2D_radiusA];
        Double_t halfwidthA   = vec_radii_peak_positions_radius_merge_copy[3][i_2D_radiusA];
        Double_t RMS_widthA   = vec_radii_peak_positions_radius_merge_copy[4][i_2D_radiusA];
        Double_t qualityA     = vec_radii_peak_positions_radius_merge_copy[5][i_2D_radiusA];
        Double_t mean_heightA = vec_radii_peak_positions_radius_merge_copy[6][i_2D_radiusA];
        Double_t mean_widthA  = vec_radii_peak_positions_radius_merge_copy[7][i_2D_radiusA];
        Double_t N_points_PRA = vec_radii_peak_positions_radius_merge_copy[8][i_2D_radiusA];
        Double_t max_heightA  = vec_radii_peak_positions_radius_merge_copy[9][i_2D_radiusA];
        Double_t min_heightA  = vec_radii_peak_positions_radius_merge_copy[10][i_2D_radiusA];
        Double_t ext_densityA = vec_radii_peak_positions_radius_merge_copy[11][i_2D_radiusA];


        Int_t flag_mergeB = 0;
        for(Int_t i_2D_radiusB = i_2D_radiusA+1; i_2D_radiusB < vec_radii_peak_positions_radius_merge_copy[0].size(); i_2D_radiusB++)
        {
            if(i_2D_radiusB > (Int_t)vec_radii_peak_positions_radius_merge_copy[0].size()) continue;
            Double_t radiusB      = vec_radii_peak_positions_radius_merge_copy[0][i_2D_radiusB];
            Double_t densityB     = vec_radii_peak_positions_radius_merge_copy[1][i_2D_radiusB];
            Double_t RMS_heightB  = vec_radii_peak_positions_radius_merge_copy[2][i_2D_radiusB];
            Double_t halfwidthB   = vec_radii_peak_positions_radius_merge_copy[3][i_2D_radiusB];
            Double_t RMS_widthB   = vec_radii_peak_positions_radius_merge_copy[4][i_2D_radiusB];
            Double_t qualityB     = vec_radii_peak_positions_radius_merge_copy[5][i_2D_radiusB];
            Double_t mean_heightB = vec_radii_peak_positions_radius_merge_copy[6][i_2D_radiusB];
            Double_t mean_widthB  = vec_radii_peak_positions_radius_merge_copy[7][i_2D_radiusB];
            Double_t N_points_PRB = vec_radii_peak_positions_radius_merge_copy[8][i_2D_radiusB];
            Double_t max_heightB  = vec_radii_peak_positions_radius_merge_copy[9][i_2D_radiusB];
            Double_t min_heightB  = vec_radii_peak_positions_radius_merge_copy[10][i_2D_radiusB];
            Double_t ext_densityB = vec_radii_peak_positions_radius_merge_copy[11][i_2D_radiusB];

            if(fabs((radiusB-halfwidthB) - (radiusA+halfwidthA)) < ring_dist_merge) // Merge rings
            {
                cout << "" << endl;
                printf("Merge rings, radii: {%f,%f}, widths: {%f,%f} \n",radiusA,radiusB,halfwidthA,halfwidthB);
                Double_t radius_merge      = ((radiusA - halfwidthA) + (radiusB + halfwidthB))/2.0;
                Double_t density_merge     = (densityA + densityB);
                Double_t RMS_height_merge  = (RMS_heightA + RMS_heightB)/2.0;
                Double_t halfwidth_merge   = (radiusB + halfwidthB) - radius_merge;
                Double_t RMS_width_merge   = (RMS_widthA + RMS_widthB)/2.0;
                Double_t quality_merge     = (qualityA + qualityB);
                Double_t mean_height_merge = (mean_heightA + mean_heightB)/2.0;
                Double_t mean_width_merge  = (mean_widthA + mean_widthB)/2.0;
                Double_t N_points_PR_merge = (N_points_PRA + N_points_PRB);
                Double_t max_height_merge  = TMath::Max(max_heightA,max_heightB);
                Double_t min_height_merge  = TMath::Min(min_heightA,min_heightB);
                Double_t ext_density_merge = (ext_densityA + ext_densityB)/2.0;


                vec_radii_peak_positions_radius_merge[0].push_back(radius_merge);
                vec_radii_peak_positions_radius_merge[1].push_back(density_merge);
                vec_radii_peak_positions_radius_merge[2].push_back(RMS_height_merge);
                vec_radii_peak_positions_radius_merge[3].push_back(halfwidth_merge);
                vec_radii_peak_positions_radius_merge[4].push_back(RMS_width_merge);
                vec_radii_peak_positions_radius_merge[5].push_back(quality_merge);
                vec_radii_peak_positions_radius_merge[6].push_back(mean_height_merge);
                vec_radii_peak_positions_radius_merge[7].push_back(mean_width_merge);
                vec_radii_peak_positions_radius_merge[8].push_back(N_points_PR_merge);
                vec_radii_peak_positions_radius_merge[9].push_back(max_height_merge);
                vec_radii_peak_positions_radius_merge[10].push_back(min_height_merge);
                vec_radii_peak_positions_radius_merge[11].push_back(ext_density_merge);


                vector< vector<Double_t> > vec_array_angle_start_stop_B_copy = vec_radii_partial_ring_angles_merge_copy[i_2D_radiusB];
                for(Int_t i_angle_partial_ringA = 0; i_angle_partial_ringA < (Int_t)vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA].size(); i_angle_partial_ringA++)
                {
                    Double_t angle_startA        = vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA][i_angle_partial_ringA][0];
                    Double_t angle_stopA         = vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA][i_angle_partial_ringA][1];
                    Double_t angle_start_merge   = angle_startA;
                    Double_t angle_stop_merge    = angle_stopA;
                    Int_t    flag_merge_A        = 0;
                    Int_t    first_index_merge_B = 0;
                    for(Int_t i_angle_partial_ringB = 0; i_angle_partial_ringB < (Int_t)vec_array_angle_start_stop_B_copy.size(); i_angle_partial_ringB++)
                    {
                        Double_t angle_startB = vec_array_angle_start_stop_B_copy[i_angle_partial_ringB][0];
                        Double_t angle_stopB  = vec_array_angle_start_stop_B_copy[i_angle_partial_ringB][1];

                        //printf("Angles A: {%f,%f}, angles B: {%f,%f}, flag_merge_A: %d \n",angle_start_merge,angle_stop_merge,angle_startB,angle_stopB,flag_merge_A);
                        if((angle_startB <= angle_stop_merge && angle_startB >= angle_start_merge)
                           || (angle_stopB <= angle_stop_merge && angle_stopB >= angle_start_merge)
                           || (angle_start_merge <= angle_stopB && angle_start_merge >= angle_startB)
                           || (angle_stop_merge <= angle_stopB && angle_stop_merge >= angle_startB)
                          )
                        {
                            angle_start_merge = TMath::Min(angle_start_merge,angle_startB);
                            angle_stop_merge  = TMath::Max(angle_stop_merge,angle_stopB);
                            //printf("  -->Merge angles: {%f,%f} \n",angle_start_merge,angle_stop_merge);
                            vec_angle_partial_ring_start_stop[0] = angle_start_merge;
                            vec_angle_partial_ring_start_stop[1] = angle_stop_merge;
                            if(flag_merge_A == 0) first_index_merge_B = i_angle_partial_ringB;
                            flag_merge_A++;
                        }
                    }
                    if(flag_merge_A == 0)
                    {
                        vec_angle_partial_ring_start_stop[0] = angle_start_merge;
                        vec_angle_partial_ring_start_stop[1] = angle_stop_merge;
                        vec_array_angle_start_stop_B_copy.push_back(vec_angle_partial_ring_start_stop);
                        //printf("    New entry created: {%f,%f} \n",angle_start_merge,angle_stop_merge);
                    }
                    else
                    {
                        //printf("    -->Erase indices: [%d,%d) \n",first_index_merge_B,first_index_merge_B+flag_merge_A);
                        vec_array_angle_start_stop_B_copy[first_index_merge_B] = vec_angle_partial_ring_start_stop;
                        vec_array_angle_start_stop_B_copy.erase(vec_array_angle_start_stop_B_copy.begin()+first_index_merge_B+1, vec_array_angle_start_stop_B_copy.begin()+first_index_merge_B+flag_merge_A);
                    }
                }

                vec_radii_partial_ring_angles_merge.push_back(vec_array_angle_start_stop_B_copy);
                //vec_array_angle_start_stop.clear();
                vec_flag_ring_merged[i_2D_radiusA] = 1;
                vec_flag_ring_merged[i_2D_radiusB] = 1;
                flag_mergeB = 1;
            }
            if(flag_mergeB) break;
        }
        if(!vec_flag_ring_merged[i_2D_radiusA]) // "A" ring was not merged with one of the "B" rings -> simply copy "A" ring
        {
            vec_radii_peak_positions_radius_merge[0].push_back(radiusA);
            vec_radii_peak_positions_radius_merge[1].push_back(densityA);
            vec_radii_peak_positions_radius_merge[2].push_back(RMS_heightA);
            vec_radii_peak_positions_radius_merge[3].push_back(halfwidthA);
            vec_radii_peak_positions_radius_merge[4].push_back(RMS_widthA);
            vec_radii_peak_positions_radius_merge[5].push_back(qualityA);
            vec_radii_peak_positions_radius_merge[6].push_back(mean_heightA);
            vec_radii_peak_positions_radius_merge[7].push_back(mean_widthA);
            vec_radii_peak_positions_radius_merge[8].push_back(N_points_PRA);
            vec_radii_peak_positions_radius_merge[9].push_back(max_heightA);
            vec_radii_peak_positions_radius_merge[10].push_back(min_heightA);
            vec_radii_peak_positions_radius_merge[11].push_back(ext_densityA);
            vec_radii_partial_ring_angles_merge.push_back(vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA]);
            //printf("---> i_2D_radiusA: %d, size of vec_radii_partial_ring_angles_merge: %d, size of vec_radii_peak_positions_radius_merge: %d, size of vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA]: %d \n",i_2D_radiusA,(Int_t)vec_radii_partial_ring_angles_merge.size(),(Int_t)vec_radii_peak_positions_radius_merge[0].size(),(Int_t)vec_radii_partial_ring_angles_merge_copy[i_2D_radiusA].size());
        }
    }

    vec_radii_partial_ring_angles_merge_copy   = vec_radii_partial_ring_angles_merge;
    vec_radii_peak_positions_radius_merge_copy = vec_radii_peak_positions_radius_merge;
}

vector< vector<Double_t> > Class_2D_crater_finder::get_radii_peak_positions() // ordered by quality
{
    return vec_radii_peak_positions_quality;
}

vector< vector<Double_t> > Class_2D_crater_finder::get_radii_peak_positions_merged() // ordered by quality
{
    cout << "" << endl;
    cout << "Merged ring properties" << endl;
    for(Int_t i_2D_radius = 0; i_2D_radius < (Int_t)vec_radii_peak_positions_radius_merge[0].size(); i_2D_radius++)
    {
        printf("Merged radius: %f, N_points: %f, quality: %f \n",vec_radii_peak_positions_radius_merge[0][i_2D_radius],vec_radii_peak_positions_radius_merge[1][i_2D_radius],vec_radii_peak_positions_radius_merge[5][i_2D_radius]);
    }
    cout << "" << endl;

    return vec_radii_peak_positions_radius_merge;
}

vector< vector<Double_t> > Class_2D_crater_finder::get_radii_peak_positions_ordered_by_radius()
{
    return vec_radii_peak_positions_radius;
}

vector< vector< vector<Double_t> > > Class_2D_crater_finder::get_radii_partial_ring_angles() // ordered by quality
{
    return vec_radii_partial_ring_angles;
}

vector< vector< vector<Double_t> > > Class_2D_crater_finder::get_radii_partial_ring_angles_merged() // ordered by quality
{
    //printf("size of vec_radii_partial_ring_angles_merge: %d \n",(Int_t)vec_radii_partial_ring_angles_merge.size());
    //for(Int_t i_r = 0; i_r < (Int_t)vec_radii_partial_ring_angles_merge.size(); i_r++)
    //{
    //    printf("i_r: %d, size: %d \n",i_r,(Int_t)vec_radii_partial_ring_angles_merge[i_r].size());
    //}
    return vec_radii_partial_ring_angles_merge;
}

Double_t Class_2D_crater_finder::get_ring_quality()
{
    return ring_quality;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Class_height_scan
{
    TH1D* TH1D_in;
    Double_t search_width, search_height;
    Double_t bin_width;
    TPolyMarker* pm_peaks;
    std::vector< std::vector<Double_t> > vec_peak_positions;
public:
    void  add_TH1D(TH1D* TH1D_in_a);
    void  set_finding_parameters(Double_t search_width_in, Double_t search_height_in);
    void  find_peaks();
    TPolyMarker* get_polymarker();
    std::vector< std::vector<Double_t> > get_peak_positions(); // [point][radial position,"height"]
    void  clear();
};

void Class_height_scan::set_finding_parameters(Double_t search_width_in, Double_t search_height_in)
{
    search_width    = search_width_in;
    search_height   = search_height_in;
}

void Class_height_scan::add_TH1D(TH1D* TH1D_in_a)
{
    TH1D_in = (TH1D*)TH1D_in_a->Clone("TH1D_in");
    bin_width = TH1D_in->GetBinWidth(1);
    pm_peaks = new TPolyMarker();
}

void Class_height_scan::find_peaks()
{
    Int_t max_bin = TH1D_in->GetNbinsX();
    Double_t min_val_x = TH1D_in->GetBinCenter(1);
    Double_t max_val_x = TH1D_in->GetBinCenter(max_bin);

    std::vector< std::vector<Double_t> > vec_diff;
    vec_diff.resize(7);

    vec_peak_positions.resize(2); // radial distance, height

    for(Int_t i_bin = 1; i_bin <= max_bin; i_bin++)
    {
        Double_t center_val = TH1D_in->GetBinCenter(i_bin);
        Double_t signal_val = TH1D_in->GetBinContent(i_bin);

        Double_t start_value_background_left  = center_val - search_width/2.0 - search_width/2.0; // half for the peak, half for the background
        Double_t stop_value_background_left   = center_val - search_width/2.0; // half for the peak, left side
        Double_t stop_value_background_right  = center_val + search_width/2.0 + search_width/2.0; // half for the peak, half for the background
        Double_t start_value_background_right = center_val + search_width/2.0; // half for the peak, rigth side

        Int_t bin_bgkr_left_start  = TH1D_in->FindBin(start_value_background_left);
        Int_t bin_bgkr_left_stop   = TH1D_in->FindBin(stop_value_background_left);
        Int_t bin_bgkr_right_start = TH1D_in->FindBin(start_value_background_right);
        Int_t bin_bgkr_right_stop  = TH1D_in->FindBin(stop_value_background_right);

        // Calculate background on the left side
        if(bin_bgkr_left_start == 0)     bin_bgkr_left_start = 1;
        if(bin_bgkr_left_stop > max_bin) bin_bgkr_left_stop  = max_bin;
        Int_t bins_bkgr_left_used = 0;
        Double_t bkgr_left = 0.0;
        for(Int_t i_bin_bkgr = bin_bgkr_left_start; i_bin_bkgr <= bin_bgkr_left_stop; i_bin_bkgr++)
        {
            bkgr_left += TH1D_in->GetBinContent(i_bin_bkgr);
            bins_bkgr_left_used++;
        }

        // Calculate background on the right side
        if(bin_bgkr_right_start == 0)     bin_bgkr_right_start = 1;
        if(bin_bgkr_right_stop > max_bin) bin_bgkr_right_stop  = max_bin;
        Int_t bins_bkgr_right_used = 0;
        Double_t bkgr_right = 0.0;
        for(Int_t i_bin_bkgr = bin_bgkr_right_start; i_bin_bkgr <= bin_bgkr_right_stop; i_bin_bkgr++)
        {
            bkgr_right += TH1D_in->GetBinContent(i_bin_bkgr);
            bins_bkgr_right_used++;
        }

        // Calculate total average background
        Double_t bkgr_total = bkgr_left + bkgr_right;
        Int_t bins_bkgr_total_used = bins_bkgr_left_used + bins_bkgr_right_used;
        if(bins_bkgr_total_used > 0)
        {
            bkgr_total /= (Double_t)(bins_bkgr_total_used);
        }

        // Calculate total signal
        Int_t bin_peak_left  = bin_bgkr_left_stop   + 1;
        Int_t bin_peak_right = bin_bgkr_right_start - 1;
        if(bin_peak_left == 0)  bin_peak_left = 1;
        if(bin_peak_right > max_bin)  bin_peak_right = max_bin;
        Int_t bins_peak_used = 0;
        Double_t peak_total = 0.0;
        for(Int_t i_bin_peak = bin_peak_left; i_bin_peak <= bin_peak_right; i_bin_peak++)
        {
            peak_total += TH1D_in->GetBinContent(i_bin_peak);
            bins_peak_used++;
        }
        if(bins_peak_used > 0)
        {
            peak_total /= (Double_t)bins_peak_used;
        }

        // Veto peak if left or right side background is larger than signal
        Int_t veto_bkgr = 0;
        if(bins_bkgr_left_used > 0 && bins_bkgr_right_used > 0)
        {
            bkgr_left /= (Double_t)bins_bkgr_left_used;
            bkgr_right /= (Double_t)bins_bkgr_right_used;
            if(bkgr_left > peak_total && bkgr_right > peak_total) veto_bkgr = 1.0;
            if(bkgr_left < peak_total && bkgr_right < peak_total) veto_bkgr = 1.0;
            if(
               !(
                (bkgr_left > search_height && bkgr_right < search_height) ||
                (bkgr_left < search_height && bkgr_right > search_height)
               )
              )
            {
                veto_bkgr = 1.0;
            }
        }

        Double_t diff_peak_bkgr = peak_total - search_height;
        vec_diff[0].push_back(center_val);
        vec_diff[1].push_back(diff_peak_bkgr); // S/B
        vec_diff[2].push_back(peak_total);
        vec_diff[3].push_back(veto_bkgr);
        vec_diff[4].push_back(bkgr_left);
        vec_diff[5].push_back(bkgr_right);
        vec_diff[6].push_back(signal_val);
    }

    // Find all maxima, exclude ranges within search_width after earch loop
    Int_t    max_i_diff = 0;
    while(max_i_diff >= 0 && vec_diff[0].size() > 0)
    {
        Double_t max_diff   = 1000000.0;
        Double_t max_diff_x = 0.0;
        Double_t max_height = 0.0;
        Double_t max_signal = 0.0;

        // Find maximum signal to background above minimum defined
        max_i_diff = -1;
        for(Int_t i_diff = 0; i_diff < vec_diff[0].size(); i_diff++)
        {
            //if(vec_diff[1][i_diff] > max_diff && vec_diff[1][i_diff] >= search_height && vec_diff[3][i_diff] < 1.0)
            if(fabs(vec_diff[1][i_diff]) < max_diff && vec_diff[3][i_diff] < 1.0)
            {
                max_diff_x = vec_diff[0][i_diff];
                max_diff   = vec_diff[1][i_diff];
                max_height = vec_diff[2][i_diff];
                max_i_diff = i_diff;
                max_signal = vec_diff[6][i_diff];
            }
        }

        // Find range within search_width to be erased
        Int_t i_diff_left_erase = 0;
        for(Int_t i_diff = max_i_diff; i_diff >= 0; i_diff--)
        {
            if(fabs(max_diff_x - vec_diff[0][i_diff]) > search_width/2.0)
            {
                i_diff_left_erase = i_diff;
                break;
            }
        }

        Int_t i_diff_right_erase = vec_diff[0].size();
        for(Int_t i_diff = max_i_diff; i_diff <= vec_diff[0].size(); i_diff++)
        {
            if(fabs(max_diff_x - vec_diff[0][i_diff]) > search_width/2.0)
            {
                i_diff_right_erase = i_diff;
                break;
            }
        }


        //cout << "start erase: "<<  i_diff_left_erase << ", stop erase: " << i_diff_right_erase << endl;
      
        if(max_i_diff >= 0)
        {
            pm_peaks ->SetNextPoint(max_diff_x,max_height);
            vec_peak_positions[0].push_back(max_diff_x);
            vec_peak_positions[1].push_back(max_signal);

            //cout << "max_diff_x: " << max_diff_x << ", max_diff: " << max_diff << ", max_i_diff: " << max_i_diff << ", size: " << vec_diff[0].size()
            //    << ", peak_total: " << vec_diff[2][max_i_diff] << ", bkgr_left: " << vec_diff[4][max_i_diff] << ", bkgr_right: " << vec_diff[5][max_i_diff] << endl;

        }

        vec_diff[0].erase(vec_diff[0].begin() + i_diff_left_erase, vec_diff[0].begin() + i_diff_right_erase);
        vec_diff[1].erase(vec_diff[1].begin() + i_diff_left_erase, vec_diff[1].begin() + i_diff_right_erase);
        vec_diff[2].erase(vec_diff[2].begin() + i_diff_left_erase, vec_diff[2].begin() + i_diff_right_erase);
        vec_diff[3].erase(vec_diff[3].begin() + i_diff_left_erase, vec_diff[3].begin() + i_diff_right_erase);
        vec_diff[4].erase(vec_diff[4].begin() + i_diff_left_erase, vec_diff[4].begin() + i_diff_right_erase);
        vec_diff[5].erase(vec_diff[5].begin() + i_diff_left_erase, vec_diff[5].begin() + i_diff_right_erase);
        vec_diff[6].erase(vec_diff[6].begin() + i_diff_left_erase, vec_diff[6].begin() + i_diff_right_erase);
    }
}

TPolyMarker* Class_height_scan::get_polymarker()
{
    return pm_peaks;
}

std::vector< std::vector<Double_t> > Class_height_scan::get_peak_positions()
{
    return vec_peak_positions;
}

void Class_height_scan::clear()
{
    delete pm_peaks;
    vec_peak_positions.clear();
    delete TH1D_in;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Draw_polygon(vector< vector<Double_t> > vec_points, Int_t color_fill, Int_t color_line, Double_t alpha_fill,
              Double_t alpha_line, Int_t style_line, Int_t line_width)
{
    TPolyLine* pl = new TPolyLine();

    for(Int_t i_point = 0; i_point < vec_points.size(); i_point++)
    {
        pl ->SetNextPoint(vec_points[i_point][0],vec_points[i_point][1]);
    }
    pl ->SetFillColorAlpha(color_fill,alpha_fill);
    pl ->SetLineStyle(style_line);
    pl ->SetLineColorAlpha(color_line,alpha_line);
    pl ->SetLineWidth(line_width);
    pl ->Draw("f");
    pl ->Draw("ogl");
}
//----------------------------------------------------------------------------------------


