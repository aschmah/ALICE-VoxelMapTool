

void Ana_voxels()
{
    gROOT->ProcessLine(".L voxResTree.C++");
    gROOT->ProcessLine("voxResTree t");
    TString proc_init = "t.Init(";
    //proc_init += input_file;
    proc_init += ")";
    gROOT->ProcessLine(proc_init.Data());
}