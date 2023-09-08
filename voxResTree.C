#define voxResTree_cxx
#include "voxResTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>




void voxResTree::Loop(Int_t N_events_loop, Int_t file_selected)
{
    printf("voxResTree::Loop, file_selected: %d, global_position: %d \n",file_selected,global_position);

    if(CheckBox_sectoraverage[0] ->GetState() == kButtonDown) printf("sector average for file: %d \n",0);
    if(CheckBox_sectoraverage[1] ->GetState() == kButtonDown) printf("sector average for file: %d \n",1);
    if(CheckBox_sectoraverage[2] ->GetState() == kButtonDown) printf("sector recovery for file: %d \n",0);
    if(CheckBox_sectoraverage[3] ->GetState() == kButtonDown) printf("sector recovery for file: %d \n",1);

    Double_t sign_invert[3]       = {1.0,1.0,1.0};
    Int_t    flag_gaussfilter     = 0;
    Int_t    flag_sectoraverage   = 0;
    Int_t    flag_sector_recovery = 0;
    Int_t flag_low_radii_extrapolation = 0;
    if(CheckBox_invert_X[file_selected]        ->GetState() == kButtonDown) sign_invert[0]       = -1.0;
    if(CheckBox_invert_Y[file_selected]        ->GetState() == kButtonDown) sign_invert[1]       = -1.0;
    if(CheckBox_invert_Z[file_selected]        ->GetState() == kButtonDown) sign_invert[2]       = -1.0;
    if(CheckBox_gaussfilter[file_selected]     ->GetState() == kButtonDown) flag_gaussfilter     = 1;
    if(CheckBox_sectoraverage[file_selected]   ->GetState() == kButtonDown) flag_sectoraverage   = 1;
    if(CheckBox_sectoraverage[file_selected+2] ->GetState() == kButtonDown) flag_sector_recovery = 1;
    if(CheckBox_low_radii_extrapolation        ->GetState() == kButtonDown) flag_low_radii_extrapolation = 1;
    if (fChain == 0) return;

    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_h_Distortions[i_xyz]           ->Reset();
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->Reset();
    }
    h2D_DZ_vs_Z       ->Reset();
    TP_DZ_vs_Z        ->Reset();
    TP_DX_vs_R        ->Reset();
    TP_DY_vs_R        ->Reset();
    h2D_DZ_vs_Z_trunc ->Reset();
    h2D_DX_vs_sector  ->Reset();
    TP_DX_vs_sector[file_selected]->Reset();
    TP_DY_vs_sector[file_selected]->Reset();
    TP_DX_vs_sector[file_selected+2]->Reset();
    TP_DY_vs_sector[file_selected+2]->Reset();
    for(Int_t i_hist = 0; i_hist < 20; i_hist++)
    {
        TP_DX_vs_sector_AC[file_selected][i_hist] ->Reset();
        TP_DY_vs_sector_AC[file_selected][i_hist] ->Reset();
    }
    h2D_DY_vs_sector  ->Reset();
    h2D_DZ_vs_sector  ->Reset();
    h2D_DY_vs_DX      ->Reset();
    h2D_DZ_vs_DX      ->Reset();
    h2D_DX_vs_stat    ->Reset();
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        h2D_DXYZS_vs_radius[i_xyz] ->Reset();
        TP_DXYZS_vs_radius[i_xyz]  ->Reset();
    }

    for(Int_t i_yz = 0; i_yz < 15; i_yz++)  // y over z bin
    {
        for(Int_t i_AC = 0; i_AC < 2; i_AC++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                TP_DXYZS_vs_radius_Y_AC[i_yz][i_AC][i_xyz] ->Reset();
                for(Int_t i_sector = 0; i_sector < 18; i_sector++)
                {
                    TP_DXYZS_vs_radius_Y_AC_sec[i_sector][i_yz][i_AC][i_xyz] ->Reset();

                    for(Int_t i_z = 0; i_z < 5; i_z++)
                    {
                        TP_DXYZS_vs_radius_Y_AC_sec_z[i_z][i_sector][i_yz][i_AC][i_xyz] ->Reset();
                    }
                }
            }
        }
    }


    h2D_DX_vs_radius  ->Reset();        //NEW
    h2D_DY_vs_stat    ->Reset();
    //TP_DZ_vs_Z_trunc  ->Reset();

    for(Int_t i_trunc = 0; i_trunc < 2; i_trunc++)
    {
        for(Int_t i_radius = 0; i_radius < 4; i_radius++)
        {
            vec_h2D_DZ_vs_Z[i_trunc][i_radius] ->Reset();
            if(i_trunc == 0) vec_TP_DZ_vs_Z[i_trunc][i_radius] ->Reset();
        }
    }

    for(Int_t i_z_bin = 0; i_z_bin < 10; i_z_bin++)
    {
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_h2D_DY_Y_vs_X[file_selected][i_z_bin][i_xyz]   ->Reset();
            vec_h2D_DSY_Y_vs_X[file_selected][i_z_bin][i_xyz]  ->Reset();
            vec_h2D_stat_Y_vs_X[file_selected][i_z_bin][i_xyz] ->Reset();
        }
    }

    for(Int_t i_sector = 0; i_sector < 9; i_sector++)
    {
        for(Int_t i_phi = 0; i_phi < 15; i_phi++)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                vec_h2D_DY_X_vs_Z[file_selected][i_sector][i_phi][i_xyz]   ->Reset();
                vec_h2D_DSY_X_vs_Z[file_selected][i_sector][i_phi][i_xyz]  ->Reset();
                vec_h2D_stat_X_vs_Z[file_selected][i_sector][i_phi][i_xyz] ->Reset();
            }
        }
    }


    for(Int_t i_z = 0; i_z < 2; i_z++)
    {
        for(Int_t i_zx = 0; i_zx < 5; i_zx++)
        {
            vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->Reset();
            vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->Reset();
        }
    }




    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t N_entries_loop = N_events_loop;
    if(N_events_loop == -1 || N_events_loop > nentries) N_entries_loop = nentries;

    //Float_t vec_DXYZ_vox[3][36][152][15][5]        = {0.0};
    //Float_t vec_DXYZ_vox_sec[3][36][152][15][5]    = {0.0};
    //Float_t vec_DXYZ_vox_sec_GF[3][36][152][15][5] = {0.0}; // after Gaussian filter
    vector<vector<vector<vector<vector<Float_t>>>>> vec_DXYZ_vox;
    vector<vector<vector<vector<vector<Float_t>>>>> vec_DXYZ_vox_sec;
    vector<vector<vector<vector<vector<Float_t>>>>> vec_DXYZ_vox_GF;
    vector<vector<vector<vector<vector<Float_t>>>>> vec_DXYZ_vox_sec_GF;
    vec_DXYZ_vox.resize(3);
    vec_DXYZ_vox_sec.resize(3);
    vec_DXYZ_vox_GF.resize(3);
    vec_DXYZ_vox_sec_GF.resize(3);
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_DXYZ_vox[i_xyz].resize(36);
        vec_DXYZ_vox_sec[i_xyz].resize(36);
        vec_DXYZ_vox_GF[i_xyz].resize(36);
        vec_DXYZ_vox_sec_GF[i_xyz].resize(36);
        for(Int_t i_sector = 0; i_sector < 36; i_sector++)
        {
            vec_DXYZ_vox[i_xyz][i_sector].resize(152);
            vec_DXYZ_vox_sec[i_xyz][i_sector].resize(152);
            vec_DXYZ_vox_GF[i_xyz][i_sector].resize(152);
            vec_DXYZ_vox_sec_GF[i_xyz][i_sector].resize(152);
            for(Int_t voxX = 0; voxX < 152; voxX++)
            {
                vec_DXYZ_vox[i_xyz][i_sector][voxX].resize(15);
                vec_DXYZ_vox_sec[i_xyz][i_sector][voxX].resize(15);
                vec_DXYZ_vox_GF[i_xyz][i_sector][voxX].resize(15);
                vec_DXYZ_vox_sec_GF[i_xyz][i_sector][voxX].resize(15);
                for(Int_t voxY = 0; voxY < 15; voxY++)
                {
                    vec_DXYZ_vox[i_xyz][i_sector][voxX][voxY].resize(5);
                    vec_DXYZ_vox_sec[i_xyz][i_sector][voxX][voxY].resize(5);
                    vec_DXYZ_vox_GF[i_xyz][i_sector][voxX][voxY].resize(5);
                    vec_DXYZ_vox_sec_GF[i_xyz][i_sector][voxX][voxY].resize(5);
                    for(Int_t voxZ = 0; voxZ < 5; voxZ++)
                    {
                        vec_DXYZ_vox[i_xyz][i_sector][voxX][voxY][voxZ] = 0.0;
                        vec_DXYZ_vox_sec[i_xyz][i_sector][voxX][voxY][voxZ] = 0.0;
                        vec_DXYZ_vox_GF[i_xyz][i_sector][voxX][voxY][voxZ] = 0.0;
                        vec_DXYZ_vox_sec_GF[i_xyz][i_sector][voxX][voxY][voxZ] = 0.0;
                    }
                }
            }
        }
    }
    //printf("value: %4.3f \n",vec_DXYZ_vox_GF[0][0][0][0][0]);


    N_bins_X_GF = TGNum_DeltaX_GF->GetNumberEntry()->GetNumber();
    N_bins_Y_GF = TGNum_DeltaY_GF->GetNumberEntry()->GetNumber();
    N_bins_Z_GF = TGNum_DeltaZ_GF->GetNumberEntry()->GetNumber();
    sigma_GF    = TGNum_sigma_GF ->GetNumberEntry()->GetNumber();

    scale_XYZ[0] = TGNum_scale_X[file_selected] ->GetNumberEntry()->GetNumber();
    scale_XYZ[1] = TGNum_scale_Y[file_selected] ->GetNumberEntry()->GetNumber();
    scale_XYZ[2] = TGNum_scale_Z[file_selected] ->GetNumberEntry()->GetNumber();

    offset_XYZ[0] = TGNum_scale_X[file_selected+2] ->GetNumberEntry()->GetNumber();
    offset_XYZ[1] = TGNum_scale_Y[file_selected+2] ->GetNumberEntry()->GetNumber();
    offset_XYZ[2] = TGNum_scale_Z[file_selected+2] ->GetNumberEntry()->GetNumber();

    vector<vector<vector<Double_t>>> vec_GKernel;
    vec_GKernel.resize(N_bins_X_GF*2+1);
    for(Int_t i_X = 0; i_X < (Int_t)vec_GKernel.size(); i_X++)
    {
        vec_GKernel[i_X].resize(N_bins_Y_GF*2+1);
        for(Int_t i_Y = 0; i_Y < (Int_t)vec_GKernel[i_X].size(); i_Y++)
        {
            vec_GKernel[i_X][i_Y].resize(N_bins_Z_GF*2+1);
        }
    }
    //double GKernel[5][5];
    //FilterCreation(GKernel);
    vec_FilterCreation(vec_GKernel,N_bins_X_GF,N_bins_Y_GF,N_bins_Z_GF,sigma_GF);


#if 0
    for(Int_t i_X = 0; i_X < (Int_t)vec_GKernel.size(); i_X++)
    {
        for(Int_t i_Y = 0; i_Y < (Int_t)vec_GKernel[i_X].size(); i_Y++)
        {
            for(Int_t i_Z = 0; i_Z < (Int_t)vec_GKernel[i_X][i_Y].size(); i_z++)
            {
                printf("X/Y: {%d, %d, %d}, value: %4.3f \n",i_X,i_Y,i_Z,vec_GKernel[i_X][i_Y][i_Z]);
            }
        }
    }
#endif


    Long64_t nbytes = 0, nb = 0;

    //------------------------------------------------------------------
    // Gaussian filtering
    if(flag_gaussfilter || flag_sectoraverage || flag_sector_recovery)
    {
        printf("\n-- Gaus or Sectoraverage/recovery --\n");
        for(Long64_t jentry = 0; jentry < N_entries_loop; jentry++)
        {
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;

            Int_t bvox_Z = (Int_t)bvox[0];  // Z/X index 0..4
            Int_t bvox_F = (Int_t)bvox[1];  // Y/X index 0..14
            Int_t bvox_X = (Int_t)bvox[2];  // X index   0..151

            Float_t DX = D[0];
            Float_t DY = D[1];
            Float_t DZ = D[2];

            Float_t DSX = DS[0];
            Float_t DSY = DS[1];
            Float_t DSZ = DS[2];

            Int_t sector = (Int_t)bsec;
            Int_t i_z = 1;
            Float_t sign_z = 1.0;
            if(sector >= 18)
            {
                sign_z = -1.0;
                i_z    = 0;
            }

            vec_DXYZ_vox_sec[0][sector][bvox_X][bvox_F][bvox_Z] = DX;
            vec_DXYZ_vox_sec[1][sector][bvox_X][bvox_F][bvox_Z] = DY;
            vec_DXYZ_vox_sec[2][sector][bvox_X][bvox_F][bvox_Z] = DZ;

            vec_DXYZ_vox[0][sector][bvox_X][bvox_F][bvox_Z] = DX;
            vec_DXYZ_vox[1][sector][bvox_X][bvox_F][bvox_Z] = DY;
            vec_DXYZ_vox[2][sector][bvox_X][bvox_F][bvox_Z] = DZ;
        }
    }


    if(flag_sectoraverage || flag_sector_recovery)
    {
        printf("\n-- Sectoraverage -- \n");
        for(Int_t voxX = 0; voxX < 152; voxX++)
        {
            for(Int_t voxY = 0; voxY < 15; voxY++)
            {
                for(Int_t voxZ = 0; voxZ < 5; voxZ++)
                {
                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        for(Int_t i_AC = 0; i_AC < 2; i_AC++)
                        {
                            Double_t average_value = 0.0;
                            Double_t weight_sum    = 0.0;
                            for(Int_t i_sector = (0+i_AC*18); i_sector < (18+i_AC*18); i_sector++)
                            {
                                if(CheckBox_sectors_used[i_sector] ->GetState() == kButtonDown)
                                {
                                    weight_sum    += 1.0;
                                    average_value += vec_DXYZ_vox_sec[i_xyz][i_sector][voxX][voxY][voxZ];
                                }
                            }
                            for(Int_t i_sector = (0+i_AC*18); i_sector < (18+i_AC*18); i_sector++)
                            {
                                if(weight_sum > 0.0) vec_DXYZ_vox_sec[i_xyz][i_sector][voxX][voxY][voxZ] = average_value/weight_sum;
                            }
                        }
                    }
                }
            }
        }
    }



    printf("flag_gaussfilter: %d, flag_sector_recovery: %d \n",flag_gaussfilter,flag_sector_recovery);


    if(flag_gaussfilter)
    {
        printf("\n-- Gaussfilter -- \n");
        // XALEX
        for(Int_t i_sector = 0; i_sector < 36; i_sector++)
        {
            //Double_t arr_values[3][5][5]      = {0.0};
            //Double_t arr_values_used[3][5][5] = {0.0};

            vector<vector<vector<vector<Float_t>>>> arr_values;
            vector<vector<vector<vector<Float_t>>>> arr_values_sec;
            vector<vector<vector<vector<Float_t>>>> arr_values_used;
            arr_values.resize(3);
            arr_values_sec.resize(3);
            arr_values_used.resize(3);


            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                arr_values[i_xyz].resize(N_bins_X_GF*2+1);
                arr_values_sec[i_xyz].resize(N_bins_X_GF*2+1);
                arr_values_used[i_xyz].resize(N_bins_X_GF*2+1);
                for(Int_t i_X = 0; i_X < (Int_t)arr_values[i_xyz].size(); i_X++)
                {
                    arr_values[i_xyz][i_X].resize(N_bins_Y_GF*2+1);
                    arr_values_sec[i_xyz][i_X].resize(N_bins_Y_GF*2+1);
                    arr_values_used[i_xyz][i_X].resize(N_bins_Y_GF*2+1);
                    for(Int_t i_Y = 0; i_Y < (Int_t)arr_values[i_xyz][i_X].size(); i_Y++)
                    {
                        arr_values[i_xyz][i_X][i_Y].resize(N_bins_Z_GF*2+1);
                        arr_values_sec[i_xyz][i_X][i_Y].resize(N_bins_Z_GF*2+1);
                        arr_values_used[i_xyz][i_X][i_Y].resize(N_bins_Z_GF*2+1);
                        for(Int_t i_Z = 0; i_Z < (Int_t)arr_values[i_xyz][i_X][i_Y].size(); i_Z++)
                        {
                            arr_values[i_xyz][i_X][i_Y][i_Z]      = 0.0;
                            arr_values_sec[i_xyz][i_X][i_Y][i_Z]  = 0.0;
                            arr_values_used[i_xyz][i_X][i_Y][i_Z] = 0.0;
                        }
                    }
                }
            }



            for(Int_t voxX = 0; voxX < 152; voxX++)
            {
                for(Int_t voxY = 0; voxY < 15; voxY++)
                {
                    for(Int_t voxZ = 0; voxZ < 5; voxZ++)
                    {
                        Float_t sum_weight[3]     = {0.0};
                        Float_t sum_values[3]     = {0.0};
                        Float_t sum_values_sec[3] = {0.0};
                        for(Int_t index_voxXB = -N_bins_X_GF; index_voxXB <= N_bins_X_GF; index_voxXB++)
                        {
                            Int_t voxXB = voxX + index_voxXB;
                            if(voxXB < 0) continue;
                            if(voxXB > 151) continue;
                            for(Int_t index_voxYB = -N_bins_Y_GF; index_voxYB <= N_bins_Y_GF; index_voxYB++)
                            {
                                Int_t voxYB = voxY + index_voxYB;
                                if(voxYB < 0) continue;
                                if(voxYB > 14) continue;
                                for(Int_t index_voxZB = -N_bins_Z_GF; index_voxZB <= N_bins_Z_GF; index_voxZB++)
                                {
                                    Int_t voxZB = voxZ + index_voxZB;
                                    if(voxZB < 0) continue;
                                    if(voxZB > 4) continue;
                                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                    {
                                        arr_values_used[i_xyz][index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF] = 1.0;
                                        arr_values_sec[i_xyz][index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF]  = vec_GKernel[index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF]*vec_DXYZ_vox_sec[i_xyz][i_sector][voxXB][voxYB][voxZB];
                                        arr_values[i_xyz][index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF]      = vec_GKernel[index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF]*vec_DXYZ_vox[i_xyz][i_sector][voxXB][voxYB][voxZB];
                                        sum_weight[i_xyz]     += vec_GKernel[index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF];
                                        sum_values_sec[i_xyz] += arr_values_sec[i_xyz][index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF];
                                        sum_values[i_xyz]     += arr_values[i_xyz][index_voxXB+N_bins_X_GF][index_voxYB+N_bins_Y_GF][index_voxZB+N_bins_Z_GF];
                                    }
                                }
                            }
                        }

                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            if(sum_weight[i_xyz] > 0.0)
                            {
                                sum_values[i_xyz]     /= sum_weight[i_xyz];
                                sum_values_sec[i_xyz] /= sum_weight[i_xyz];
                                vec_DXYZ_vox_sec_GF[i_xyz][i_sector][voxX][voxY][voxZ] = sum_values_sec[i_xyz];
                                vec_DXYZ_vox_GF[i_xyz][i_sector][voxX][voxY][voxZ]     = sum_values[i_xyz];
                            }
                        }
                    }
                }
            }


        }
    }
    //------------------------------------------------------------------




    //------------------------------------------------------------------
    // Loop over voxels which have no entry
    if(flag_sector_recovery || flag_sectoraverage)
    {
        for(Int_t sector = 0; sector < 36; sector++)
        {
            for(Int_t bvox_X = 0; bvox_X < 152; bvox_X++)
            {
                stat[2] = RowX[bvox_X];
                for(Int_t bvox_F = 0; bvox_F < 15; bvox_F++)
                {
                    for(Int_t bvox_Z = 0; bvox_Z < 5; bvox_Z++)
                    {
                        Int_t index_average_map  = bvox_X+152*bvox_F+152*15*bvox_Z+152*15*5*sector;

                        Int_t i_z = 1;
                        Int_t AC_side = 0;
                        Float_t sign_z = 1.0;
                        if(sector >= 18)
                        {
                            sign_z = -1.0;
                            i_z    = 0;
                            AC_side = 1;

                        }
                        if(index_average_map > 410400) printf("WARNING: index_average_map out of range! \n");

                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            if(vec_DXYZ_vox[i_xyz][sector][bvox_X][bvox_F][bvox_Z] == 0.0)
                            {
                                if(flag_gaussfilter)
                                {
                                    D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec_GF[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                                }
                                else
                                {
                                    D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                                }
                            }
                            //if(sector == 30) printf("sector: %d, flag_sector_recovery: %d index_average_map: %d, value: %4.3f \n",sector,flag_sector_recovery,index_average_map,D[i_xyz]);
                        }



                        //--------------------------------------------------
                        //h2D_DXS_vs_radius ->Fill(sign_z*stat[2],DSX);
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            h2D_DXYZS_vs_radius[i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
                            TP_DXYZS_vs_radius[i_xyz]  ->Fill(sign_z*stat[2],D[i_xyz]);
                            TP_DXYZS_vs_radius_Y_AC[bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
                            TP_DXYZS_vs_radius_Y_AC_sec[sector%18][bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
                            TP_DXYZS_vs_radius_Y_AC_sec_z[bvox_Z][sector%18][bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
                        }
                        //--------------------------------------------------





                        //--------------------------------------------------
                        if(flag_low_radii_extrapolation)
                        {
                            //printf("flag_low_radii_extrapolation: %d \n",flag_low_radii_extrapolation);
                            Int_t flag_do_extrapolation = 0;
                            Float_t x_val_ext = sign_z*stat[2];
                            if(!AC_side) // A side
                            {
                                if(x_val_ext < TGNum_fit_min_A->GetNumberEntry()->GetNumber())
                                {
                                    flag_do_extrapolation = 1;
                                }
                            }
                            else // C side
                            {
                                if(x_val_ext > TGNum_fit_max_C->GetNumberEntry()->GetNumber())
                                {
                                    flag_do_extrapolation = 1;
                                }
                            }

                            if(flag_do_extrapolation)
                            {
                                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                {
                                    //D[i_xyz] = func_PolyFitFunc_xyz_AC[i_xyz][AC_side] ->Eval(x_val_ext);
                                    if(TGR_select_low_radius_extr[0] ->GetState() == kButtonDown)
                                    {
                                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC[bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                                    }
                                    if(TGR_select_low_radius_extr[1] ->GetState() == kButtonDown)
                                    {
                                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC_sec[sector%18][bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                                    }
                                    if(TGR_select_low_radius_extr[2] ->GetState() == kButtonDown)
                                    {
                                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC_sec_z[bvox_Z][sector%18][bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                                    }
                                }
                            }
                        }
                        //--------------------------------------------------




                        //--------------------------------------------------
                        //printf("index_average_map: %d, DXorig: %4.3f \n",index_average_map,DXorig);
                        vec_VoxRes[file_selected][index_average_map].D[0]      = (float)D[0];
                        vec_VoxRes[file_selected][index_average_map].D[1]      = (float)D[1];
                        vec_VoxRes[file_selected][index_average_map].D[2]      = (float)D[2];

                        vec_VoxRes[file_selected][index_average_map].DS[0]     = (float)D[0];
                        vec_VoxRes[file_selected][index_average_map].DS[1]     = (float)D[1];
                        vec_VoxRes[file_selected][index_average_map].DS[2]     = (float)D[2];

                        vec_VoxRes[file_selected][index_average_map].DC[0]     = (float)D[0];
                        vec_VoxRes[file_selected][index_average_map].DC[1]     = (float)D[1];
                        vec_VoxRes[file_selected][index_average_map].DC[2]     = (float)D[2];

                        vec_VoxRes[file_selected][index_average_map].E[0]      = (float)0.0;
                        vec_VoxRes[file_selected][index_average_map].E[1]      = (float)0.0;
                        vec_VoxRes[file_selected][index_average_map].E[2]      = (float)0.0;

                        vec_VoxRes[file_selected][index_average_map].stat[0]   = (float)1000;
                        vec_VoxRes[file_selected][index_average_map].stat[1]   = (float)1000;
                        vec_VoxRes[file_selected][index_average_map].stat[2]   = (float)1000;
                        vec_VoxRes[file_selected][index_average_map].stat[3]   = (float)1000; // number of entries used


                        vec_VoxRes[file_selected][index_average_map].EXYCorr   = 1.0;
                        vec_VoxRes[file_selected][index_average_map].dYSigMAD  = 1.0;
                        vec_VoxRes[file_selected][index_average_map].dZSigLTM  = 1.0;

                        vec_VoxRes[file_selected][index_average_map].bvox[0]   = bvox_Z;
                        vec_VoxRes[file_selected][index_average_map].bvox[1]   = bvox_F;
                        vec_VoxRes[file_selected][index_average_map].bvox[2]   = bvox_X;
                        vec_VoxRes[file_selected][index_average_map].bsec      = sector;
                        vec_VoxRes[file_selected][index_average_map].flags     = 7;
                        //--------------------------------------------------


                        Int_t bin_y_min, bin_y_max, bin_x_min, bin_x_max;

                        //--------------------------------------------------
                        // R vs. Z plots
                        Float_t X_pad_position     = RowX[bvox_X];
                        Float_t radius = X_pad_position;
                        Float_t X_pad_position_min = RowX[bvox_X];
                        Float_t X_pad_position_max = RowX[bvox_X+1];
                        Float_t Z_pad_position_min = sign_z*(z_over_x_voxel_binning[bvox_Z]*X_pad_position);
                        Float_t Z_pad_position_max = sign_z*(z_over_x_voxel_binning[bvox_Z+1]*X_pad_position);




                        //--------------------------------------------------
                        // Y vs. X plots

                        //Int_t bvox_Z = (Int_t)bvox[0];  // Z/X index 0..4
                        //Int_t bvox_F = (Int_t)bvox[1];  // Y/X index 0..14
                        //Int_t bvox_X = (Int_t)bvox[2];  // X index   0..151


                        Int_t bvox_z_offset = 5;
                        if(sign_z < 0.0) bvox_z_offset = 4; // sign is negative: bvox_Z = 4 -> z_bin = 0, sign is positive: bvox_z = 0 -> z_bin = 5
                        Int_t z_bin = bvox_z_offset + (Int_t)(sign_z*bvox_Z);

                        // local TPC coordinate system
                        X_pad_position_min = RowX[bvox_X];
                        X_pad_position_max = RowX[bvox_X+1];

                        // range: -0.17..0.17 in 15 bins
                        Double_t y_over_x_value_min = (0.34/15.0)*bvox_F - 0.17;
                        Double_t y_over_x_value_max = (0.34/15.0)*(bvox_F+1) - 0.17;
                        Double_t Y_pad_position_min = y_over_x_value_min*X_pad_position_min;
                        Double_t Y_pad_position_max = y_over_x_value_max*X_pad_position_min;

                        bin_y_min = h2D_Y_vs_X_TPC_sector ->GetYaxis()->FindBin(X_pad_position_min);
                        bin_y_max = h2D_Y_vs_X_TPC_sector ->GetYaxis()->FindBin(X_pad_position_max);
                        bin_x_min = h2D_Y_vs_X_TPC_sector ->GetXaxis()->FindBin(Y_pad_position_min);
                        bin_x_max = h2D_Y_vs_X_TPC_sector ->GetXaxis()->FindBin(Y_pad_position_max);

                        TVector2 vec_2D_local, vec_2D_global;
                        for(Int_t i_bin_x = bin_x_min; i_bin_x < bin_x_max; i_bin_x++)
                        {
                            Double_t x_val = h2D_Y_vs_X_TPC_sector->GetXaxis()->GetBinCenter(i_bin_x);
                            for(Int_t i_bin_y = bin_y_min; i_bin_y < bin_y_max; i_bin_y++)
                            {
                                Double_t y_val = h2D_Y_vs_X_TPC_sector->GetYaxis()->GetBinCenter(i_bin_y);
                                vec_2D_local.SetX(-x_val);
                                vec_2D_local.SetY(y_val);
                                Double_t phi_val = (-90.0 + (sector%18)*20.0 + 10.0)*TMath::DegToRad(); // rotation from local to global TPC coordinate system
                                vec_2D_global = vec_2D_local.Rotate(phi_val);

                                Int_t i_bin_x_rot = vec_h2D_DY_Y_vs_X[file_selected][z_bin][0]  ->GetXaxis()->FindBin(vec_2D_global.X());
                                Int_t i_bin_y_rot = vec_h2D_DY_Y_vs_X[file_selected][z_bin][0]  ->GetYaxis()->FindBin(vec_2D_global.Y());

                                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                                {
                                    Double_t bin_cont = vec_h2D_DY_Y_vs_X[file_selected][z_bin][i_xyz]   ->GetBinContent(i_bin_x_rot,i_bin_y_rot);
                                    //if(bin_cont > 0.0) printf("WARNING! Bin already filled! \n");
                                    vec_h2D_DY_Y_vs_X[file_selected][z_bin][i_xyz]   ->SetBinContent(i_bin_x_rot,i_bin_y_rot,D[i_xyz]);
                                    vec_h2D_DSY_Y_vs_X[file_selected][z_bin][i_xyz]  ->SetBinContent(i_bin_x_rot,i_bin_y_rot,DS[i_xyz]);
                                }
                                vec_h2D_stat_Y_vs_X[file_selected][z_bin][0] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,stat[3]);
                                vec_h2D_stat_Y_vs_X[file_selected][z_bin][1] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,flags);
                                vec_h2D_stat_Y_vs_X[file_selected][z_bin][2] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,dYSigMAD);
                            }
                        }
                        //--------------------------------------------------




                    }
                }
            }
        }
    }
    //------------------------------------------------------------------




#if 1
    //------------------------------------------------------------------------
    Float_t maxDX = 0.0;
    Float_t maxDY = 0.0;
    Float_t maxDZ = 0.0;
    nbytes = 0;
    nb = 0;
    for(Long64_t jentry = 0; jentry < N_entries_loop; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        Int_t bvox_Z = (Int_t)bvox[0];  // Z/X index 0..4
        Int_t bvox_F = (Int_t)bvox[1];  // Y/X index 0..14
        Int_t bvox_X = (Int_t)bvox[2];  // X index   0..151


        Float_t DSX = DS[0];
        Float_t DSY = DS[1];
        Float_t DSZ = DS[2];


        Int_t sector = (Int_t)bsec;
        Int_t i_z = 1;
        Int_t AC_side = 0;
        Float_t sign_z = 1.0;
        if(sector >= 18)
        {
            sign_z = -1.0;
            i_z    = 0;
            AC_side = 1;

        }

        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            D[i_xyz] *= scale_XYZ[i_xyz]*sign_invert[i_xyz] + offset_XYZ[i_xyz];
        }

        Float_t DXorig = D[0];
        Float_t DYorig = D[1];
        Float_t DZorig = D[2];

        if(flag_sectoraverage && !flag_sector_recovery)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
            }
        }

        if(flag_sector_recovery)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                if(CheckBox_sectors_used[sector] ->GetState() == kButtonDown)
                {
                    D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                }
                else
                {
                    D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                    //printf("recovery of sector: %d, xyz: %d, voxel: {%d, %d, %d}, value: %4.3f \n",sector,i_xyz,bvox_X,bvox_F,bvox_Z,D[i_xyz]);
                }
            }
        }

        if(flag_gaussfilter)
        {
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_GF[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                if(flag_sectoraverage && !flag_sector_recovery) D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec_GF[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                if(flag_sector_recovery)
                {
                    if(CheckBox_sectors_used[sector] ->GetState() == kButtonDown)
                    {
                        D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_GF[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                    }
                    else
                    {
                        D[i_xyz] = scale_XYZ[i_xyz]*sign_invert[i_xyz]*vec_DXYZ_vox_sec_GF[i_xyz][sector][bvox_X][bvox_F][bvox_Z] + offset_XYZ[i_xyz];
                    }
                }
            }
        }



        //--------------------------------------------------
        //h2D_DXS_vs_radius ->Fill(sign_z*stat[2],DSX);
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            h2D_DXYZS_vs_radius[i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
            TP_DXYZS_vs_radius[i_xyz]  ->Fill(sign_z*stat[2],D[i_xyz]);
            TP_DXYZS_vs_radius_Y_AC[bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
            TP_DXYZS_vs_radius_Y_AC_sec[sector%18][bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
            TP_DXYZS_vs_radius_Y_AC_sec_z[bvox_Z][sector%18][bvox_F][AC_side][i_xyz] ->Fill(sign_z*stat[2],D[i_xyz]);
        }
        //--------------------------------------------------



        //--------------------------------------------------
        if(flag_low_radii_extrapolation)
        {
            //printf("flag_low_radii_extrapolation: %d \n",flag_low_radii_extrapolation);
            Int_t flag_do_extrapolation = 0;
            Float_t x_val_ext = sign_z*stat[2];
            if(!AC_side) // A side
            {
                if(x_val_ext < TGNum_fit_min_A->GetNumberEntry()->GetNumber())
                {
                    flag_do_extrapolation = 1;
                }
            }
            else // C side
            {
                if(x_val_ext > TGNum_fit_max_C->GetNumberEntry()->GetNumber())
                {
                    flag_do_extrapolation = 1;
                }
            }

            if(flag_do_extrapolation)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    //D[i_xyz] = func_PolyFitFunc_xyz_AC[i_xyz][AC_side] ->Eval(x_val_ext);
                    if(TGR_select_low_radius_extr[0] ->GetState() == kButtonDown)
                    {
                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC[bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                    }
                    if(TGR_select_low_radius_extr[1] ->GetState() == kButtonDown)
                    {
                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC_sec[sector%18][bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                    }
                    if(TGR_select_low_radius_extr[2] ->GetState() == kButtonDown)
                    {
                        D[i_xyz] = func_PolyFitFunc_xyz_Y_AC_sec_z[bvox_Z][sector%18][bvox_F][AC_side][i_xyz] ->Eval(x_val_ext);
                    }
                }
            }
        }
        //--------------------------------------------------




        Float_t DX = D[0];
        Float_t DY = D[1];
        Float_t DZ = D[2];

        //if(bvox_Z == 0 && bvox_X == 120 && sector == 0)
        //{
        //    printf("bvox_F: %d, D: {%4.3f, %4.3f, %4.3f} \n",bvox_F,DX,DY,DZ);
        //}


        //--------------------------------------------------
        Int_t index_average_map  = bvox_X+152*bvox_F+152*15*bvox_Z+152*15*5*sector;
        if(index_average_map > 410400) printf("WARNING: index_average_map out of range! \n");

        //printf("index_average_map: %d, DXorig: %4.3f \n",index_average_map,DXorig);
        vec_VoxRes[file_selected][index_average_map].D[0]      = (float)DXorig;
        vec_VoxRes[file_selected][index_average_map].D[1]      = (float)DYorig;
        vec_VoxRes[file_selected][index_average_map].D[2]      = (float)DZorig;

        vec_VoxRes[file_selected][index_average_map].DS[0]     = (float)DX;
        vec_VoxRes[file_selected][index_average_map].DS[1]     = (float)DY;
        vec_VoxRes[file_selected][index_average_map].DS[2]     = (float)DZ;

        vec_VoxRes[file_selected][index_average_map].DC[0]     = (float)DX;
        vec_VoxRes[file_selected][index_average_map].DC[1]     = (float)DY;
        vec_VoxRes[file_selected][index_average_map].DC[2]     = (float)DZ;

        vec_VoxRes[file_selected][index_average_map].E[0]      = (float)E[0];
        vec_VoxRes[file_selected][index_average_map].E[1]      = (float)E[1];
        vec_VoxRes[file_selected][index_average_map].E[2]      = (float)E[2];

        vec_VoxRes[file_selected][index_average_map].stat[0]   = (float)stat[0];
        vec_VoxRes[file_selected][index_average_map].stat[1]   = (float)stat[1];
        vec_VoxRes[file_selected][index_average_map].stat[2]   = (float)stat[2];
        vec_VoxRes[file_selected][index_average_map].stat[3]   = (float)stat[3]; // number of entries used


        vec_VoxRes[file_selected][index_average_map].EXYCorr   = EXYCorr;
        vec_VoxRes[file_selected][index_average_map].dYSigMAD  = dYSigMAD;
        vec_VoxRes[file_selected][index_average_map].dZSigLTM  = dZSigLTM;

        vec_VoxRes[file_selected][index_average_map].bvox[0]   = bvox[0];
        vec_VoxRes[file_selected][index_average_map].bvox[1]   = bvox[1];
        vec_VoxRes[file_selected][index_average_map].bvox[2]   = bvox[2];
        vec_VoxRes[file_selected][index_average_map].bsec      = bsec;
        vec_VoxRes[file_selected][index_average_map].flags     = 7;
        //--------------------------------------------------


        Int_t bin_y_min, bin_y_max, bin_x_min, bin_x_max;

        //--------------------------------------------------
        // R vs. Z plots
        Float_t X_pad_position     = RowX[bvox_X];
        Float_t radius = X_pad_position;
        Float_t X_pad_position_min = RowX[bvox_X];
        Float_t X_pad_position_max = RowX[bvox_X+1];
        Float_t Z_pad_position_min = sign_z*(z_over_x_voxel_binning[bvox_Z]*X_pad_position);
        Float_t Z_pad_position_max = sign_z*(z_over_x_voxel_binning[bvox_Z+1]*X_pad_position);


        vec_TP_DZ_vs_DX_tanTheta[i_z][bvox_Z] ->Fill(DSX,DSZ);
        vec_TP_DY_vs_DX_tanTheta[i_z][bvox_Z] ->Fill(DSX,DSY);

        //if(sector == 4) printf("X_pad: {%4.3f, %4.3f}, Z_pad: {%4.3f, %4.3f}, vox: {%4.3f, %4.3f, %4.3f}, entries: %4.3f \n",X_pad_position_min,X_pad_position_max,Z_pad_position_min,Z_pad_position_max,stat[0],stat[1],stat[2],stat[3]);

        //Double_t Z_global = stat[0] + 0.5*(Z_pad_position_max + Z_pad_position_min);
        Double_t Z_global = stat[0]*stat[2];
        //printf("stat[0]: %4.3f, stat[1]: %4.3f, stat[2]: %4.3f \n",stat[0],stat[1],stat[2]);
        if(stat[3] > 100 && flags < 20)
        {
            h2D_DZ_vs_Z ->Fill(Z_global,DZ);
            TP_DZ_vs_Z  ->Fill(Z_global,DZ);

            vec_h2D_DZ_vs_Z[0][0] ->Fill(Z_global,DZ);
            vec_TP_DZ_vs_Z[0][0]  ->Fill(Z_global,DZ);
            if(radius < 110.0)
            {
                vec_h2D_DZ_vs_Z[0][1] ->Fill(Z_global,DZ);
                vec_TP_DZ_vs_Z[0][1]  ->Fill(Z_global,DZ);
            }
            if(radius >= 110.0 && radius < 200.0)
            {
                vec_h2D_DZ_vs_Z[0][2] ->Fill(Z_global,DZ);
                vec_TP_DZ_vs_Z[0][2]  ->Fill(Z_global,DZ);
            }
            if(radius >= 200.0)
            {
                vec_h2D_DZ_vs_Z[0][3] ->Fill(Z_global,DZ);
                vec_TP_DZ_vs_Z[0][3]  ->Fill(Z_global,DZ);
            }
        }


        Int_t plot_sector = sector%9;
        Int_t plot_sector_opposite = plot_sector + 9;

        //if(bvox_F == 2)
        {
            if(sector == plot_sector || sector == plot_sector_opposite || sector == (plot_sector + 18) || sector == (plot_sector_opposite + 18))
            {
                if(sector == plot_sector_opposite || sector == (plot_sector_opposite + 18))
                {
                    X_pad_position_min *= -1.0;
                    X_pad_position_max *= -1.0;
                }

                bin_y_min = h2D_DY_X_vs_Z ->GetYaxis()->FindBin(X_pad_position_min);
                bin_y_max = h2D_DY_X_vs_Z ->GetYaxis()->FindBin(X_pad_position_max);
                bin_x_min = h2D_DY_X_vs_Z ->GetXaxis()->FindBin(Z_pad_position_min);
                bin_x_max = h2D_DY_X_vs_Z ->GetXaxis()->FindBin(Z_pad_position_max);


                Int_t bin_dummy = 0;
                if(bin_x_min > bin_x_max)
                {
                    bin_dummy = bin_x_max;
                    bin_x_max = bin_x_min;
                    bin_x_min = bin_dummy;

                }
                if(bin_y_min > bin_y_max)
                {
                    bin_dummy = bin_y_max;
                    bin_y_max = bin_y_min;
                    bin_y_min = bin_dummy;

                }

                for(Int_t i_bin_x = bin_x_min; i_bin_x < bin_x_max; i_bin_x++)
                {
                    for(Int_t i_bin_y = bin_y_min; i_bin_y < bin_y_max; i_bin_y++)
                    {
                        h2D_DY_X_vs_Z                          ->SetBinContent(i_bin_x,i_bin_y,DY);
                        //if(plot_sector == 0 && bvox_F == 0) printf("plot_sector: %d, bvox_F: %d \n",plot_sector,bvox_F);
                        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                        {
                            vec_h2D_DY_X_vs_Z[file_selected][plot_sector][bvox_F][i_xyz]  ->SetBinContent(i_bin_x,i_bin_y,D[i_xyz]);
                            vec_h2D_DSY_X_vs_Z[file_selected][plot_sector][bvox_F][i_xyz] ->SetBinContent(i_bin_x,i_bin_y,DS[i_xyz]);
                        }
                        vec_h2D_stat_X_vs_Z[file_selected][plot_sector][bvox_F][0] ->SetBinContent(i_bin_x,i_bin_y,stat[3]);
                        vec_h2D_stat_X_vs_Z[file_selected][plot_sector][bvox_F][1] ->SetBinContent(i_bin_x,i_bin_y,flags);
                        vec_h2D_stat_X_vs_Z[file_selected][plot_sector][bvox_F][2] ->SetBinContent(i_bin_x,i_bin_y,dYSigMAD);
                    }
                }
            }
        }
        //--------------------------------------------------




        //--------------------------------------------------
        // Y vs. X plots

        //Int_t bvox_Z = (Int_t)bvox[0];  // Z/X index 0..4
        //Int_t bvox_F = (Int_t)bvox[1];  // Y/X index 0..14
        //Int_t bvox_X = (Int_t)bvox[2];  // X index   0..151


        Int_t bvox_z_offset = 5;
        if(sign_z < 0.0) bvox_z_offset = 4; // sign is negative: bvox_Z = 4 -> z_bin = 0, sign is positive: bvox_z = 0 -> z_bin = 5
        Int_t z_bin = bvox_z_offset + (Int_t)(sign_z*bvox_Z);

        // local TPC coordinate system
        X_pad_position_min = RowX[bvox_X];
        X_pad_position_max = RowX[bvox_X+1];

        // range: -0.17..0.17 in 15 bins
        Double_t y_over_x_value_min = (0.34/15.0)*bvox_F - 0.17;
        Double_t y_over_x_value_max = (0.34/15.0)*(bvox_F+1) - 0.17;
        Double_t Y_pad_position_min = y_over_x_value_min*X_pad_position_min;
        Double_t Y_pad_position_max = y_over_x_value_max*X_pad_position_min;

        bin_y_min = h2D_Y_vs_X_TPC_sector ->GetYaxis()->FindBin(X_pad_position_min);
        bin_y_max = h2D_Y_vs_X_TPC_sector ->GetYaxis()->FindBin(X_pad_position_max);
        bin_x_min = h2D_Y_vs_X_TPC_sector ->GetXaxis()->FindBin(Y_pad_position_min);
        bin_x_max = h2D_Y_vs_X_TPC_sector ->GetXaxis()->FindBin(Y_pad_position_max);

        TVector2 vec_2D_local, vec_2D_global;
        for(Int_t i_bin_x = bin_x_min; i_bin_x < bin_x_max; i_bin_x++)
        {
            Double_t x_val = h2D_Y_vs_X_TPC_sector->GetXaxis()->GetBinCenter(i_bin_x);
            for(Int_t i_bin_y = bin_y_min; i_bin_y < bin_y_max; i_bin_y++)
            {
                Double_t y_val = h2D_Y_vs_X_TPC_sector->GetYaxis()->GetBinCenter(i_bin_y);
                vec_2D_local.SetX(-x_val);
                vec_2D_local.SetY(y_val);
                Double_t phi_val = (-90.0 + (sector%18)*20.0 + 10.0)*TMath::DegToRad(); // rotation from local to global TPC coordinate system
                vec_2D_global = vec_2D_local.Rotate(phi_val);

                Int_t i_bin_x_rot = vec_h2D_DY_Y_vs_X[file_selected][z_bin][0]  ->GetXaxis()->FindBin(vec_2D_global.X());
                Int_t i_bin_y_rot = vec_h2D_DY_Y_vs_X[file_selected][z_bin][0]  ->GetYaxis()->FindBin(vec_2D_global.Y());

                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    Double_t bin_cont = vec_h2D_DY_Y_vs_X[file_selected][z_bin][i_xyz]   ->GetBinContent(i_bin_x_rot,i_bin_y_rot);
                    //if(bin_cont > 0.0) printf("WARNING! Bin already filled! \n");
                    vec_h2D_DY_Y_vs_X[file_selected][z_bin][i_xyz]   ->SetBinContent(i_bin_x_rot,i_bin_y_rot,D[i_xyz]);
                    vec_h2D_DSY_Y_vs_X[file_selected][z_bin][i_xyz]  ->SetBinContent(i_bin_x_rot,i_bin_y_rot,DS[i_xyz]);
                }
                vec_h2D_stat_Y_vs_X[file_selected][z_bin][0] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,stat[3]);
                vec_h2D_stat_Y_vs_X[file_selected][z_bin][1] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,flags);
                vec_h2D_stat_Y_vs_X[file_selected][z_bin][2] ->SetBinContent(i_bin_x_rot,i_bin_y_rot,dYSigMAD);
            }
        }
        //--------------------------------------------------


        if(DX > maxDX) maxDX = DX;
        if(DY > maxDY) maxDY = DY;
        if(DZ > maxDZ) maxDZ = DZ;

        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_h_Distortions[i_xyz]           ->Fill(D[i_xyz]);
            vec_h2D_Distortions_vs_voxZ[i_xyz] ->Fill(sign_z*(bvox_Z+0.5),D[i_xyz]);
        }


        //if(sector == 0)
        {
            h2D_DX_vs_stat ->Fill(stat[3],D[0]);
            h2D_DY_vs_stat ->Fill(stat[3],D[1]);
        }

        h2D_DX_vs_radius  ->Fill(sign_z*stat[2],DX);


        if(bvox[0] == 0 && stat[3] > 0)
        {
            Int_t x_bin = h2D_DX_vs_sector->GetXaxis()->FindBin(bsec+((Double_t)bvox[1])/15.0);
            Int_t y_bin = h2D_DX_vs_sector->GetYaxis()->FindBin(D[0]);
            //h2D_DX_vs_sector ->SetBinContent(x_bin,y_bin,bvox[2]);
            h2D_DX_vs_sector ->SetBinContent(x_bin,y_bin,stat[2]); // radial component (stat[2]) as color

            x_bin = h2D_DY_vs_sector->GetXaxis()->FindBin(bsec+((Double_t)bvox[1])/15.0);
            y_bin = h2D_DY_vs_sector->GetYaxis()->FindBin(D[1]);
            h2D_DY_vs_sector ->SetBinContent(x_bin,y_bin,stat[2]);

            x_bin = h2D_DZ_vs_sector->GetXaxis()->FindBin(bsec+((Double_t)bvox[1])/15.0);
            y_bin = h2D_DZ_vs_sector->GetYaxis()->FindBin(D[2]);
            h2D_DZ_vs_sector ->SetBinContent(x_bin,y_bin,stat[2]);
        }

        if(bvox[0] == 0 && stat[3] > 0)
        {
            if(bvox[2] >= 130 && bvox[2] <= 134)
            {
                TP_DX_vs_sector[file_selected] ->Fill(bsec+((Double_t)bvox[1])/15.0,D[0]);
                TP_DY_vs_sector[file_selected] ->Fill(bsec+((Double_t)bvox[1])/15.0,D[1]);
            }
            if(bvox[2] >= 10 && bvox[2] <= 20)
            {
                TP_DX_vs_sector[file_selected+2] ->Fill(bsec+((Double_t)bvox[1])/15.0,D[0]);
                TP_DY_vs_sector[file_selected+2] ->Fill(bsec+((Double_t)bvox[1])/15.0,D[1]);
            }
        }

        //if(sector == 21 && stat[2] > 200.0)
        {
            h2D_DY_vs_DX ->Fill(D[0],D[1]);
            h2D_DZ_vs_DX ->Fill(D[0],D[2]);
        }

        TP_DX_vs_R        ->Fill(sign_z*stat[2],D[0]);
        TP_DY_vs_R        ->Fill(sign_z*stat[2],D[1]);

        //if(sign_z > 0) printf("jentry: %lld, bvox_Z: %d, bvox_F: %d, bvox_X: %d, dist: {%4.3f, %4.3f, %4.3f}, sign_z*bvox_Z+1: %4.2f \n",jentry,bvox_Z,bvox_F,bvox_X,DX,DY,DZ,sign_z*(bvox_Z+1));
        // if (Cut(ientry) < 0) continue;
    }
    printf("max dist: {%4.3f, %4.3f, %4.3f} \n",maxDX,maxDY,maxDZ);
    //------------------------------------------------------------------------






    //------------------------------------------------------------------------
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetFillColor(10);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetTopMargin(0.1);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetBottomMargin(0.2);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetRightMargin(0.05);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetLeftMargin(0.2);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetTicks(1,1);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetGrid(0,0);
        can_vec_h_Distortions ->cd(i_xyz + 1);
        can_vec_h_Distortions ->cd(i_xyz + 1)->SetLogy(1);
        vec_h_Distortions[i_xyz] ->SetFillColor(kBlack);
        vec_h_Distortions[i_xyz] ->SetFillStyle(3002);
        vec_h_Distortions[i_xyz] ->SetLineColor(kBlack);
        vec_h_Distortions[i_xyz] ->GetXaxis()->CenterTitle();
        vec_h_Distortions[i_xyz] ->GetYaxis()->CenterTitle();
        vec_h_Distortions[i_xyz] ->SetStats(0);
        vec_h_Distortions[i_xyz] ->SetTitle("");
        vec_h_Distortions[i_xyz] ->GetXaxis()->SetTitleOffset(1.2);
        vec_h_Distortions[i_xyz] ->GetYaxis()->SetTitleOffset(1.3);
        vec_h_Distortions[i_xyz] ->GetXaxis()->SetLabelSize(0.06);
        vec_h_Distortions[i_xyz] ->GetYaxis()->SetLabelSize(0.06);
        vec_h_Distortions[i_xyz] ->GetXaxis()->SetTitleSize(0.06);
        vec_h_Distortions[i_xyz] ->GetYaxis()->SetTitleSize(0.06);
        vec_h_Distortions[i_xyz] ->GetXaxis()->SetNdivisions(505,'N');
        vec_h_Distortions[i_xyz] ->GetYaxis()->SetNdivisions(505,'N');
        vec_h_Distortions[i_xyz] ->GetXaxis()->SetTitle(arr_label_xyz[i_xyz].Data());
        vec_h_Distortions[i_xyz] ->GetYaxis()->SetTitle("entries");
        vec_h_Distortions[i_xyz] ->DrawCopy("hist");
        can_vec_h_Distortions ->cd(i_xyz + 1)->Update();
    }
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetFillColor(10);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetTopMargin(0.1);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetBottomMargin(0.2);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetRightMargin(0.05);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetLeftMargin(0.2);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetTicks(1,1);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetGrid(0,0);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1);
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->SetLogz(1);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->CenterTitle();
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->CenterTitle();
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->SetStats(0);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->SetTitle("");
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->SetTitleOffset(1.2);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->SetTitleOffset(1.3);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->SetLabelSize(0.06);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->SetLabelSize(0.06);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->SetTitleSize(0.06);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->SetTitleSize(0.06);
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->SetNdivisions(505,'N');
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->SetNdivisions(505,'N');
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetXaxis()->SetTitle("vox Z");
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetYaxis()->SetTitle(arr_label_xyz[i_xyz].Data());
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->GetZaxis()->SetTitle("entries");
        vec_h2D_Distortions_vs_voxZ[i_xyz] ->DrawCopy("colz");
        can_vec_h2D_Distortions_vs_voxZ ->cd(i_xyz + 1)->Update();
    }
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    can_h2D_DX_vs_sector ->cd();
    can_h2D_DX_vs_sector ->cd()->SetFillColor(10);
    can_h2D_DX_vs_sector ->cd()->SetTopMargin(0.1);
    can_h2D_DX_vs_sector ->cd()->SetBottomMargin(0.2);
    can_h2D_DX_vs_sector ->cd()->SetRightMargin(0.15);
    can_h2D_DX_vs_sector ->cd()->SetLeftMargin(0.2);
    can_h2D_DX_vs_sector ->cd()->SetTicks(1,1);
    can_h2D_DX_vs_sector ->cd()->SetGrid(0,0);
    can_h2D_DX_vs_sector ->cd();
    can_h2D_DX_vs_sector ->cd()->SetLogz(0);
    h2D_DX_vs_sector ->GetXaxis()->CenterTitle();
    h2D_DX_vs_sector ->GetYaxis()->CenterTitle();
    h2D_DX_vs_sector ->SetStats(0);
    h2D_DX_vs_sector ->SetTitle("");
    h2D_DX_vs_sector ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DX_vs_sector ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DX_vs_sector ->GetXaxis()->SetLabelSize(0.06);
    h2D_DX_vs_sector ->GetYaxis()->SetLabelSize(0.06);
    h2D_DX_vs_sector ->GetXaxis()->SetTitleSize(0.06);
    h2D_DX_vs_sector ->GetYaxis()->SetTitleSize(0.06);
    h2D_DX_vs_sector ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_sector ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_sector ->GetXaxis()->SetTitle("sector");
    h2D_DX_vs_sector ->GetYaxis()->SetTitle("#DeltaX (cm)");
    h2D_DX_vs_sector ->GetZaxis()->SetTitle("radius (cm)");
    h2D_DX_vs_sector ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DX_vs_sector ->DrawCopy("colz");

    func_modulation ->SetParameter(0,1.0);
    func_modulation ->SetParameter(1,TMath::Pi()*2.0/9.0);
    func_modulation ->SetParameter(2,-(TMath::Pi()/2.0) - 0.9);
    func_modulation ->SetParameter(3,-0.35);
    func_modulation ->SetLineColor(kAzure-2);
    func_modulation ->SetLineWidth(4);
    func_modulation ->SetLineStyle(1);
    //func_modulation ->DrawCopy("same");

    can_h2D_DX_vs_sector ->Update();
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    can_TP_DX_vs_sector->cd(1);
    TP_DX_vs_sector[0] ->GetXaxis()->CenterTitle();
    TP_DX_vs_sector[0] ->GetYaxis()->CenterTitle();
    TP_DX_vs_sector[0] ->SetStats(0);
    TP_DX_vs_sector[0] ->SetTitle("");
    TP_DX_vs_sector[0] ->GetXaxis()->SetTitleOffset(1.3);
    TP_DX_vs_sector[0] ->GetYaxis()->SetTitleOffset(0.3);
    TP_DX_vs_sector[0] ->GetXaxis()->SetLabelSize(0.06);
    TP_DX_vs_sector[0] ->GetYaxis()->SetLabelSize(0.06);
    TP_DX_vs_sector[0] ->GetXaxis()->SetTitleSize(0.06);
    TP_DX_vs_sector[0] ->GetYaxis()->SetTitleSize(0.06);
    TP_DX_vs_sector[0] ->GetXaxis()->SetNdivisions(505,'N');
    TP_DX_vs_sector[0] ->GetYaxis()->SetNdivisions(505,'N');
    TP_DX_vs_sector[0] ->GetXaxis()->SetTitle("sector");
    TP_DX_vs_sector[0] ->GetYaxis()->SetTitle("#DeltaX (cm)");
    TP_DX_vs_sector[0] ->GetYaxis()->SetRangeUser(-1.8,1.8);
    TP_DX_vs_sector[0] ->SetLineColor(kBlack);
    TP_DX_vs_sector[0] ->DrawCopy("hist");
    TP_DX_vs_sector[1] ->SetLineColor(kRed);
    TP_DX_vs_sector[1] ->DrawCopy("same hist");

    TP_DX_vs_sector[2] ->SetLineColor(kGray+1);
    TP_DX_vs_sector[2] ->DrawCopy("same hist");
    TP_DX_vs_sector[3] ->SetLineColor(kRed+2);
    TP_DX_vs_sector[3] ->DrawCopy("same hist");

    for(Int_t i_sector = 0; i_sector < 36; i_sector++)
    {
        PlotLine(0.0 + i_sector,0.0 + i_sector,-1.5,1.5,kGray+2,1,3); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    }


    can_TP_DX_vs_sector->cd(2);
    TP_DY_vs_sector[0] ->GetXaxis()->CenterTitle();
    TP_DY_vs_sector[0] ->GetYaxis()->CenterTitle();
    TP_DY_vs_sector[0] ->SetStats(0);
    TP_DY_vs_sector[0] ->SetTitle("");
    TP_DY_vs_sector[0] ->GetXaxis()->SetTitleOffset(1.3);
    TP_DY_vs_sector[0] ->GetYaxis()->SetTitleOffset(0.3);
    TP_DY_vs_sector[0] ->GetXaxis()->SetLabelSize(0.06);
    TP_DY_vs_sector[0] ->GetYaxis()->SetLabelSize(0.06);
    TP_DY_vs_sector[0] ->GetXaxis()->SetTitleSize(0.06);
    TP_DY_vs_sector[0] ->GetYaxis()->SetTitleSize(0.06);
    TP_DY_vs_sector[0] ->GetXaxis()->SetNdivisions(505,'N');
    TP_DY_vs_sector[0] ->GetYaxis()->SetNdivisions(505,'N');
    TP_DY_vs_sector[0] ->GetXaxis()->SetTitle("sector");
    TP_DY_vs_sector[0] ->GetYaxis()->SetTitle("#DeltaY (cm)");
    TP_DY_vs_sector[0] ->GetYaxis()->SetRangeUser(-1.8,1.8);
    TP_DY_vs_sector[0] ->SetLineColor(kBlack);
    TP_DY_vs_sector[0] ->DrawCopy("hist");
    TP_DY_vs_sector[1] ->SetLineColor(kRed);
    TP_DY_vs_sector[1] ->DrawCopy("same hist");

    // small radii
    TP_DY_vs_sector[2] ->SetLineColor(kGray+1);
    TP_DY_vs_sector[2] ->DrawCopy("same hist");
    TP_DY_vs_sector[3] ->SetLineColor(kRed+2);
    TP_DY_vs_sector[3] ->DrawCopy("same hist");

    for(Int_t i_sector = 0; i_sector < 36; i_sector++)
    {
        PlotLine(0.0 + i_sector,0.0 + i_sector,-1.5,1.5,kGray+2,1,3); // Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle
    }

    can_TP_DX_vs_sector->cd(1) ->Update();
    can_TP_DX_vs_sector->cd(2) ->Update();
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    can_h2D_DY_vs_sector ->cd();
    can_h2D_DY_vs_sector ->cd()->SetFillColor(10);
    can_h2D_DY_vs_sector ->cd()->SetTopMargin(0.1);
    can_h2D_DY_vs_sector ->cd()->SetBottomMargin(0.2);
    can_h2D_DY_vs_sector ->cd()->SetRightMargin(0.15);
    can_h2D_DY_vs_sector ->cd()->SetLeftMargin(0.2);
    can_h2D_DY_vs_sector ->cd()->SetTicks(1,1);
    can_h2D_DY_vs_sector ->cd()->SetGrid(0,0);
    can_h2D_DY_vs_sector ->cd();
    can_h2D_DY_vs_sector ->cd()->SetLogz(0);
    h2D_DY_vs_sector ->GetXaxis()->CenterTitle();
    h2D_DY_vs_sector ->GetYaxis()->CenterTitle();
    h2D_DY_vs_sector ->SetStats(0);
    h2D_DY_vs_sector ->SetTitle("");
    h2D_DY_vs_sector ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DY_vs_sector ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DY_vs_sector ->GetXaxis()->SetLabelSize(0.06);
    h2D_DY_vs_sector ->GetYaxis()->SetLabelSize(0.06);
    h2D_DY_vs_sector ->GetXaxis()->SetTitleSize(0.06);
    h2D_DY_vs_sector ->GetYaxis()->SetTitleSize(0.06);
    h2D_DY_vs_sector ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_sector ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_sector ->GetXaxis()->SetTitle("sector");
    h2D_DY_vs_sector ->GetYaxis()->SetTitle("#DeltaY (cm)");
    h2D_DY_vs_sector ->GetZaxis()->SetTitle("radius (cm)");
    h2D_DY_vs_sector ->GetYaxis()->SetRangeUser(-2.5,2.5);
    h2D_DY_vs_sector ->DrawCopy("colz");

    func_modulation ->SetParameter(0,0.5);
    func_modulation ->SetParameter(1,TMath::Pi()*2.0/9.0);
    func_modulation ->SetParameter(2,-(TMath::Pi()/2.0) - 0.9);
    func_modulation ->SetParameter(3,-0.72);
    func_modulation ->SetLineColor(kAzure-2);
    func_modulation ->SetLineWidth(4);
    func_modulation ->SetLineStyle(1);
    func_modulation ->DrawCopy("same");

    can_h2D_DY_vs_sector ->Update();
    //------------------------------------------------------------------------


    //------------------------------------------------------------------------
    can_h2D_DZ_vs_sector ->cd();
    can_h2D_DZ_vs_sector ->cd()->SetFillColor(10);
    can_h2D_DZ_vs_sector ->cd()->SetTopMargin(0.1);
    can_h2D_DZ_vs_sector ->cd()->SetBottomMargin(0.2);
    can_h2D_DZ_vs_sector ->cd()->SetRightMargin(0.15);
    can_h2D_DZ_vs_sector ->cd()->SetLeftMargin(0.2);
    can_h2D_DZ_vs_sector ->cd()->SetTicks(1,1);
    can_h2D_DZ_vs_sector ->cd()->SetGrid(0,0);
    can_h2D_DZ_vs_sector ->cd();
    can_h2D_DZ_vs_sector ->cd()->SetLogz(0);
    h2D_DZ_vs_sector ->GetXaxis()->CenterTitle();
    h2D_DZ_vs_sector ->GetYaxis()->CenterTitle();
    h2D_DZ_vs_sector ->SetStats(0);
    h2D_DZ_vs_sector ->SetTitle("");
    h2D_DZ_vs_sector ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DZ_vs_sector ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DZ_vs_sector ->GetXaxis()->SetLabelSize(0.06);
    h2D_DZ_vs_sector ->GetYaxis()->SetLabelSize(0.06);
    h2D_DZ_vs_sector ->GetXaxis()->SetTitleSize(0.06);
    h2D_DZ_vs_sector ->GetYaxis()->SetTitleSize(0.06);
    h2D_DZ_vs_sector ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DZ_vs_sector ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DZ_vs_sector ->GetXaxis()->SetTitle("sector");
    h2D_DZ_vs_sector ->GetYaxis()->SetTitle("#DeltaZ (cm)");
    h2D_DZ_vs_sector ->GetZaxis()->SetTitle("radius (cm)");
    h2D_DZ_vs_sector ->GetYaxis()->SetRangeUser(-2.5,2.5);
    h2D_DZ_vs_sector ->DrawCopy("colz");
    can_h2D_DZ_vs_sector ->Update();
    //------------------------------------------------------------------------




    //------------------------------------------------------------------------
    can_h2D_DY_vs_DX ->cd();
    can_h2D_DY_vs_DX ->cd()->SetFillColor(10);
    can_h2D_DY_vs_DX ->cd()->SetTopMargin(0.1);
    can_h2D_DY_vs_DX ->cd()->SetBottomMargin(0.2);
    can_h2D_DY_vs_DX ->cd()->SetRightMargin(0.15);
    can_h2D_DY_vs_DX ->cd()->SetLeftMargin(0.2);
    can_h2D_DY_vs_DX ->cd()->SetTicks(1,1);
    can_h2D_DY_vs_DX ->cd()->SetGrid(0,0);
    can_h2D_DY_vs_DX ->cd();
    can_h2D_DY_vs_DX ->cd()->SetLogz(0);
    h2D_DY_vs_DX ->GetXaxis()->CenterTitle();
    h2D_DY_vs_DX ->GetYaxis()->CenterTitle();
    h2D_DY_vs_DX ->SetStats(0);
    h2D_DY_vs_DX ->SetTitle("");
    h2D_DY_vs_DX ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DY_vs_DX ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DY_vs_DX ->GetXaxis()->SetLabelSize(0.06);
    h2D_DY_vs_DX ->GetYaxis()->SetLabelSize(0.06);
    h2D_DY_vs_DX ->GetXaxis()->SetTitleSize(0.06);
    h2D_DY_vs_DX ->GetYaxis()->SetTitleSize(0.06);
    h2D_DY_vs_DX ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_DX ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_DX ->GetXaxis()->SetTitle("#DeltaX (cm)");
    h2D_DY_vs_DX ->GetYaxis()->SetTitle("#DeltaY (cm)");
    h2D_DY_vs_DX ->GetZaxis()->SetTitle("entries");
    //h2D_DY_vs_DX ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DY_vs_DX ->DrawCopy("colz");
    can_h2D_DY_vs_DX ->Update();
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    can_h2D_DZ_vs_DX ->cd();
    can_h2D_DZ_vs_DX ->cd()->SetFillColor(10);
    can_h2D_DZ_vs_DX ->cd()->SetTopMargin(0.1);
    can_h2D_DZ_vs_DX ->cd()->SetBottomMargin(0.2);
    can_h2D_DZ_vs_DX ->cd()->SetRightMargin(0.15);
    can_h2D_DZ_vs_DX ->cd()->SetLeftMargin(0.2);
    can_h2D_DZ_vs_DX ->cd()->SetTicks(1,1);
    can_h2D_DZ_vs_DX ->cd()->SetGrid(0,0);
    can_h2D_DZ_vs_DX ->cd();
    can_h2D_DZ_vs_DX ->cd()->SetLogz(0);
    h2D_DZ_vs_DX ->GetXaxis()->CenterTitle();
    h2D_DZ_vs_DX ->GetYaxis()->CenterTitle();
    h2D_DZ_vs_DX ->SetStats(0);
    h2D_DZ_vs_DX ->SetTitle("");
    h2D_DZ_vs_DX ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DZ_vs_DX ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DZ_vs_DX ->GetXaxis()->SetLabelSize(0.06);
    h2D_DZ_vs_DX ->GetYaxis()->SetLabelSize(0.06);
    h2D_DZ_vs_DX ->GetXaxis()->SetTitleSize(0.06);
    h2D_DZ_vs_DX ->GetYaxis()->SetTitleSize(0.06);
    h2D_DZ_vs_DX ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DZ_vs_DX ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DZ_vs_DX ->GetXaxis()->SetTitle("#DeltaX (cm)");
    h2D_DZ_vs_DX ->GetYaxis()->SetTitle("#DeltaZ (cm)");
    h2D_DZ_vs_DX ->GetZaxis()->SetTitle("entries");
    //h2D_DZ_vs_DX ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DZ_vs_DX ->DrawCopy("colz");
    can_h2D_DZ_vs_DX ->Update();
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    can_h2D_DX_vs_stat ->cd();
    can_h2D_DX_vs_stat ->cd()->SetFillColor(10);
    can_h2D_DX_vs_stat ->cd()->SetTopMargin(0.1);
    can_h2D_DX_vs_stat ->cd()->SetBottomMargin(0.2);
    can_h2D_DX_vs_stat ->cd()->SetRightMargin(0.15);
    can_h2D_DX_vs_stat ->cd()->SetLeftMargin(0.2);
    can_h2D_DX_vs_stat ->cd()->SetTicks(1,1);
    can_h2D_DX_vs_stat ->cd()->SetGrid(0,0);
    can_h2D_DX_vs_stat ->cd();
    can_h2D_DX_vs_stat ->cd()->SetLogz(1);
    h2D_DX_vs_stat ->GetXaxis()->CenterTitle();
    h2D_DX_vs_stat ->GetYaxis()->CenterTitle();
    h2D_DX_vs_stat ->SetStats(0);
    h2D_DX_vs_stat ->SetTitle("");
    h2D_DX_vs_stat ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DX_vs_stat ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DX_vs_stat ->GetXaxis()->SetLabelSize(0.06);
    h2D_DX_vs_stat ->GetYaxis()->SetLabelSize(0.06);
    h2D_DX_vs_stat ->GetXaxis()->SetTitleSize(0.06);
    h2D_DX_vs_stat ->GetYaxis()->SetTitleSize(0.06);
    h2D_DX_vs_stat ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_stat ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_stat ->GetXaxis()->SetTitle("statistics");
    h2D_DX_vs_stat ->GetYaxis()->SetTitle("#DeltaX (cm)");
    h2D_DX_vs_stat ->GetZaxis()->SetTitle("entries");
    //h2D_DX_vs_stat ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DX_vs_stat ->DrawCopy("colz");
    can_h2D_DX_vs_stat ->Update();
    //------------------------------------------------------------------------




    //------------------------------------------------------------------------
    can_h2D_DY_vs_stat ->cd();
    can_h2D_DY_vs_stat ->cd()->SetFillColor(10);
    can_h2D_DY_vs_stat ->cd()->SetTopMargin(0.1);
    can_h2D_DY_vs_stat ->cd()->SetBottomMargin(0.2);
    can_h2D_DY_vs_stat ->cd()->SetRightMargin(0.15);
    can_h2D_DY_vs_stat ->cd()->SetLeftMargin(0.2);
    can_h2D_DY_vs_stat ->cd()->SetTicks(1,1);
    can_h2D_DY_vs_stat ->cd()->SetGrid(0,0);
    can_h2D_DY_vs_stat ->cd();
    can_h2D_DY_vs_stat ->cd()->SetLogz(1);
    h2D_DY_vs_stat ->GetXaxis()->CenterTitle();
    h2D_DY_vs_stat ->GetYaxis()->CenterTitle();
    h2D_DY_vs_stat ->SetStats(0);
    h2D_DY_vs_stat ->SetTitle("");
    h2D_DY_vs_stat ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DY_vs_stat ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DY_vs_stat ->GetXaxis()->SetLabelSize(0.06);
    h2D_DY_vs_stat ->GetYaxis()->SetLabelSize(0.06);
    h2D_DY_vs_stat ->GetXaxis()->SetTitleSize(0.06);
    h2D_DY_vs_stat ->GetYaxis()->SetTitleSize(0.06);
    h2D_DY_vs_stat ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_stat ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DY_vs_stat ->GetXaxis()->SetTitle("statistics");
    h2D_DY_vs_stat ->GetYaxis()->SetTitle("#DeltaY (cm)");
    h2D_DY_vs_stat ->GetZaxis()->SetTitle("entries");
    //h2D_DY_vs_stat ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DY_vs_stat ->DrawCopy("colz");
    can_h2D_DY_vs_stat ->Update();
    //------------------------------------------------------------------------


#if 0
    // -----------------------------------------------------------------------
    //                    X-Residuals (DX) vs Radius (R)
    //------------------------------------------------------------------------
    can_TP_DX_vs_R ->cd();
    can_TP_DX_vs_R ->cd()->SetFillColor(10);
    can_TP_DX_vs_R ->cd()->SetTopMargin(0.1);
    can_TP_DX_vs_R ->cd()->SetBottomMargin(0.2);
    can_TP_DX_vs_R ->cd()->SetRightMargin(0.15);
    can_TP_DX_vs_R ->cd()->SetLeftMargin(0.2);
    can_TP_DX_vs_R ->cd()->SetTicks(1,1);
    can_TP_DX_vs_R ->cd()->SetGrid(0,0);
    can_TP_DX_vs_R ->cd();
    can_TP_DX_vs_R ->cd()->SetLogz(0);
    TP_DX_vs_R ->GetXaxis()->CenterTitle();
    TP_DX_vs_R ->GetYaxis()->CenterTitle();
    TP_DX_vs_R ->SetStats(0);
    TP_DX_vs_R ->SetTitle("");
    TP_DX_vs_R ->GetXaxis()->SetTitleOffset(1.2);
    TP_DX_vs_R ->GetYaxis()->SetTitleOffset(1.3);
    TP_DX_vs_R ->GetXaxis()->SetLabelSize(0.06);
    TP_DX_vs_R ->GetYaxis()->SetLabelSize(0.06);
    TP_DX_vs_R ->GetXaxis()->SetTitleSize(0.06);
    TP_DX_vs_R ->GetYaxis()->SetTitleSize(0.06);
    TP_DX_vs_R ->GetXaxis()->SetNdivisions(505,'N');
    TP_DX_vs_R ->GetYaxis()->SetNdivisions(505,'N');
    TP_DX_vs_R ->GetXaxis()->SetTitle("R (cm)");
    TP_DX_vs_R ->GetYaxis()->SetTitle("#DeltaX (cm)");
    //TP_DX_vs_R ->GetYaxis()->SetRangeUser(-4.5,4.5);
    TP_DX_vs_R ->DrawCopy("hist");
    can_TP_DX_vs_R ->Update();
    //------------------------------------------------------------------------
#endif


#if 0
    //------------------------------------------------------------------------
    can_TP_DY_vs_R ->cd();
    can_TP_DY_vs_R ->cd()->SetFillColor(10);
    can_TP_DY_vs_R ->cd()->SetTopMargin(0.1);
    can_TP_DY_vs_R ->cd()->SetBottomMargin(0.2);
    can_TP_DY_vs_R ->cd()->SetRightMargin(0.15);
    can_TP_DY_vs_R ->cd()->SetLeftMargin(0.2);
    can_TP_DY_vs_R ->cd()->SetTicks(1,1);
    can_TP_DY_vs_R ->cd()->SetGrid(0,0);
    can_TP_DY_vs_R ->cd();
    can_TP_DY_vs_R ->cd()->SetLogz(0);
    TP_DY_vs_R ->GetXaxis()->CenterTitle();
    TP_DY_vs_R ->GetYaxis()->CenterTitle();
    TP_DY_vs_R ->SetStats(0);
    TP_DY_vs_R ->SetTitle("");
    TP_DY_vs_R ->GetXaxis()->SetTitleOffset(1.2);
    TP_DY_vs_R ->GetYaxis()->SetTitleOffset(1.3);
    TP_DY_vs_R ->GetXaxis()->SetLabelSize(0.06);
    TP_DY_vs_R ->GetYaxis()->SetLabelSize(0.06);
    TP_DY_vs_R ->GetXaxis()->SetTitleSize(0.06);
    TP_DY_vs_R ->GetYaxis()->SetTitleSize(0.06);
    TP_DY_vs_R ->GetXaxis()->SetNdivisions(505,'N');
    TP_DY_vs_R ->GetYaxis()->SetNdivisions(505,'N');
    TP_DY_vs_R ->GetXaxis()->SetTitle("R (cm)");
    TP_DY_vs_R ->GetYaxis()->SetTitle("#DeltaY (cm)");
    //TP_DY_vs_R ->GetYaxis()->SetRangeUser(-4.5,4.5);
    TP_DY_vs_R ->DrawCopy("hist");
    can_TP_DY_vs_R ->Update();
    //------------------------------------------------------------------------
#endif


#if 0
    // -----------------------------------------------------------------------
    can_h2D_and_h1D_DX_vs_radius ->Divide(1,2);
    can_h2D_and_h1D_DX_vs_radius ->cd(1);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetFillColor(10);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetTopMargin(0.1);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetBottomMargin(0.2);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetRightMargin(0.15);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetLeftMargin(0.2);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetTicks(1,1);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetGrid(0,0);
    can_h2D_and_h1D_DX_vs_radius ->cd(1);
    can_h2D_and_h1D_DX_vs_radius ->cd(1)->SetLogz(1);
    h2D_DX_vs_radius ->GetXaxis()->CenterTitle();
    h2D_DX_vs_radius ->GetYaxis()->CenterTitle();
    h2D_DX_vs_radius ->SetStats(0);
    h2D_DX_vs_radius ->SetTitle("");
    h2D_DX_vs_radius ->GetXaxis()->SetTitleOffset(1.2);
    h2D_DX_vs_radius ->GetYaxis()->SetTitleOffset(1.3);
    h2D_DX_vs_radius ->GetXaxis()->SetLabelSize(0.06);
    h2D_DX_vs_radius ->GetYaxis()->SetLabelSize(0.06);
    h2D_DX_vs_radius ->GetXaxis()->SetTitleSize(0.06);
    h2D_DX_vs_radius ->GetYaxis()->SetTitleSize(0.06);
    h2D_DX_vs_radius ->GetXaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_radius ->GetYaxis()->SetNdivisions(505,'N');
    h2D_DX_vs_radius ->GetXaxis()->SetTitle("radius (cm)");
    h2D_DX_vs_radius ->GetYaxis()->SetTitle("#DeltaX (cm)");
    h2D_DX_vs_radius ->GetZaxis()->SetTitle("entries");
    //h2D_DX_vs_radius ->GetYaxis()->SetRangeUser(-4.5,4.5);
    h2D_DX_vs_radius ->DrawCopy("colz");

    can_h2D_and_h1D_DX_vs_radius ->cd(2);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetFillColor(10);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetTopMargin(0.1);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetBottomMargin(0.2);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetRightMargin(0.15);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetLeftMargin(0.2);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetTicks(1,1);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetGrid(0,0);
    can_h2D_and_h1D_DX_vs_radius ->cd(2);
    can_h2D_and_h1D_DX_vs_radius ->cd(2)->SetLogz(0);
    TP_DX_vs_R ->DrawCopy("hist");
    can_h2D_and_h1D_DX_vs_radius ->Update();
    //------------------------------------------------------------------------
#endif


    //------------------------------------------------------------------------
    Int_t arr_color[2][5] =
    {
        {kRed,kRed+1,kRed+2,kRed+3,kRed+4},
        {kAzure,kAzure-2,kAzure-3,kAzure-3,kAzure-4}
    };

    can_vec_TP_DZ_vs_DX_tanTheta ->cd();
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetFillColor(10);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetTopMargin(0.1);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetBottomMargin(0.2);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetRightMargin(0.15);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetLeftMargin(0.2);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetTicks(1,1);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetGrid(0,0);
    can_vec_TP_DZ_vs_DX_tanTheta ->cd();
    can_vec_TP_DZ_vs_DX_tanTheta ->cd()->SetLogz(0);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->CenterTitle();
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->CenterTitle();
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->SetStats(0);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->SetTitle("");
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitleOffset(1.2);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitleOffset(1.3);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->SetLabelSize(0.06);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetLabelSize(0.06);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitleSize(0.06);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitleSize(0.06);
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->SetNdivisions(505,'N');
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetNdivisions(505,'N');
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitle("#DeltaX (cm)");
    vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitle("#DeltaZ (cm)");
    //vec_TP_DZ_vs_DX_tanTheta[0][0] ->GetYaxis()->SetRangeUser(-4.5,4.5);
    for(Int_t i_z = 0; i_z < 2; i_z++)
    {
        for(Int_t i_zx = 0; i_zx < 5; i_zx++)
        {
            vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerStyle(20);
            vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerSize(0.8);
            vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerColor(arr_color[i_z][i_zx]);

            if(i_z == 0 && i_zx == 0)
            {
                vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->DrawCopy("P");
            }
            else
            {
              vec_TP_DZ_vs_DX_tanTheta[i_z][i_zx] ->DrawCopy("same P");
            }
        }
    }
    can_vec_TP_DZ_vs_DX_tanTheta ->Update();


    can_vec_TP_DY_vs_DX_tanTheta ->cd();
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetFillColor(10);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetTopMargin(0.1);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetBottomMargin(0.2);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetRightMargin(0.15);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetLeftMargin(0.2);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetTicks(1,1);
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetGrid(0,0);
    can_vec_TP_DY_vs_DX_tanTheta ->cd();
    can_vec_TP_DY_vs_DX_tanTheta ->cd()->SetLogz(0);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->CenterTitle();
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->CenterTitle();
    vec_TP_DY_vs_DX_tanTheta[0][0] ->SetStats(0);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->SetTitle("");
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitleOffset(1.2);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitleOffset(1.3);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->SetLabelSize(0.06);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetLabelSize(0.06);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitleSize(0.06);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitleSize(0.06);
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->SetNdivisions(505,'N');
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetNdivisions(505,'N');
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetXaxis()->SetTitle("#DeltaX (cm)");
    vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetTitle("#DeltaY (cm)");
    //vec_TP_DY_vs_DX_tanTheta[0][0] ->GetYaxis()->SetRangeUser(-4.5,4.5);
    for(Int_t i_z = 0; i_z < 2; i_z++)
    {
        for(Int_t i_zx = 0; i_zx < 5; i_zx++)
        {
            vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerStyle(20);
            vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerSize(0.8);
            vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->SetMarkerColor(arr_color[i_z][i_zx]);

            if(i_z == 0 && i_zx == 0)
            {
                vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->DrawCopy("P");
            }
            else
            {
              vec_TP_DY_vs_DX_tanTheta[i_z][i_zx] ->DrawCopy("same P");
            }
        }
    }
    can_vec_TP_DY_vs_DX_tanTheta ->Update();

    //------------------------------------------------------------------------


    Update_DY_X_vs_Z();


#endif
}
