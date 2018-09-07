
#include <iostream>
#include <fstream>

#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatElasticity2D.h"
#include "TPorousElasticity.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzanalysis.h"
//#include "pznonlinanalysis.h"



// Geometry utilities
TPZGeoMesh * ReadGeometry();
void PrintGeometry(TPZGeoMesh * gmesh);

// Computational utilities
TPZCompMesh * DeformationMesh(TPZGeoMesh * gmesh, int p_order);
TPZAnalysis * Analysis(TPZCompMesh * cmesh);

// Post-process utilities
void PostProcess(TPZAnalysis *an);


int main() 
{
    int p_order = 1;
    TPZGeoMesh * gmesh = ReadGeometry();
    PrintGeometry(gmesh);
    
    TPZCompMesh * cmesh = DeformationMesh(gmesh, p_order);
    TPZAnalysis * analysis = Analysis(cmesh);
    
    TPZFMatrix<STATE> x(analysis->Solution()), dx;
    x.Zero();
    REAL tol = 1.0e-6;
    int n_it = 20;
    bool stop_criterion_Q = false;
    REAL norm_res;
    for (int i = 0; i < n_it; i++) {
        analysis->Assemble();
        analysis->Rhs() *= -1.0;
        analysis->Solve();
        dx = analysis->Solution();
        x += dx;
        analysis->LoadSolution(x);
        analysis->AssembleResidual();
        norm_res = Norm(analysis->Rhs());
        stop_criterion_Q = norm_res < tol;
        if (stop_criterion_Q) {
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i << std::endl;
            break;
        }
    }
    
    if (stop_criterion_Q == false) {
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }

    PostProcess(analysis);
    std::cout << "Execution complete." << std::endl;
    return 0;
}

TPZGeoMesh * ReadGeometry(){

    TPZGeoMesh * gmesh = new TPZGeoMesh;
    std::string grid("Wellbore.msh");
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Wellbore section");
    gmesh->SetName(name);
    
    return gmesh;
}

void PrintGeometry(TPZGeoMesh * gmesh){
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << "geometry" << ".txt";
    vtk_name    << "geometry" << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

TPZCompMesh * DeformationMesh(TPZGeoMesh * gmesh, int p_order){
    
    unsigned int dim  = gmesh->Dimension();
    const std::string name("Porous Elasticity on wellbore ");
    REAL to_MPa = 1.0;
    
    // Setting up attributes
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName(name);
    cmesh->SetDefaultOrder(p_order);
    cmesh->SetDimModel(dim);
    
    // Frist running with linear elasticity.
    unsigned int rock_id = 1;
//    TPZMatElasticity2D * rock = new TPZMatElasticity2D(rock_id);
//    REAL Ey = 29269.0*to_MPa;
//    REAL nu = 0.20300;
//    rock->SetPlaneStrain();
//    rock->SetElasticity(Ey, nu);

    TPorousElasticity * rock = new TPorousElasticity(rock_id);
    
    STATE nu = 0.203;
    STATE kappa = 0.024;
    STATE pt_el = 5.835;
    STATE e_0 = 0.34;
    STATE p_0 = 0.0;
    TPZManVector<STATE,6> s_0(6);
    
    rock->SetPlaneStrain();
    rock->SetPorousElasticity(kappa, pt_el, e_0, p_0);
    rock->SetPoissonRatioConstant(nu);
    rock->SetInitialStress(s_0);
    
//    TPZFNMatrix<36,STATE> De(6,6);
//    TPZFNMatrix<6,STATE> epsilon_vec(6,1),sigma_vec(6,1);
//    De.Zero();
//    TPZTensor<STATE> epsilon, sigma;
//    epsilon.XX() = -0.00005;
//    epsilon.XY() = -0.000003;
//    epsilon.XZ() = -0.000004;
//    epsilon.YY() = -0.000005;
//    epsilon.YZ() = -0.000006;
//    epsilon.ZZ() = -0.00007;
//
//    rock->Sigma(epsilon, sigma);
//    rock->De(epsilon, De);
//    
//    epsilon.CopyTo(epsilon_vec);
//    De.Multiply(epsilon_vec, sigma_vec);
//    
//    
//    sigma.Print(std::cout);
//    sigma_vec.Print(std::cout);
//    De.Print(std::cout);
    
    unsigned int bc_i_id, bc_e_id, bc_index, bc_id;
    
    bc_i_id = 2;
    bc_index = 5;
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 1.0*to_MPa;
    TPZBndCond * bc_i = rock->CreateBC(rock, bc_i_id, bc_index, val1, val2);
    
    bc_e_id = 3;
    bc_index = 5;
    val2(0,0) = 0.0*to_MPa;
    TPZBndCond * bc_e = rock->CreateBC(rock, bc_e_id, bc_index, val1, val2);
    
    bc_id = 4;
    bc_index = 7;
    val2(0,0) = 0.0;
    TPZBndCond * bc_ux_fixed = rock->CreateBC(rock, bc_id, bc_index, val1, val2);
    
    bc_id = 5;
    bc_index = 8;
    val2(0,0) = 0.0;
    TPZBndCond * bc_uy_fixed = rock->CreateBC(rock, bc_id, bc_index, val1, val2);
    
    cmesh->InsertMaterialObject(rock);
    cmesh->InsertMaterialObject(bc_i);
    cmesh->InsertMaterialObject(bc_e);
    cmesh->InsertMaterialObject(bc_ux_fixed);
    cmesh->InsertMaterialObject(bc_uy_fixed);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZAnalysis * Analysis(TPZCompMesh * cmesh){
    
    int numofThreads = 0;
    TPZAnalysis * analysis = new TPZAnalysis(cmesh,true);
    TPZSkylineStructMatrix matrix(cmesh);
//    TPZSymetricSpStructMatrix matrix(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    matrix.SetNumThreads(numofThreads);
    analysis->SetStructuralMatrix(matrix);
    analysis->SetSolver(step);
    return analysis;
}

void PostProcess(TPZAnalysis *an)
{
    const int dim = an->Mesh()->Dimension();
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "Wellbore.vtk";
    
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("SigmaZ");
    vecnames.Push("Displacement");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(div);
}

