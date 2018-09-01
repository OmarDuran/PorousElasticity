
#include <iostream>

#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

// Geometry utilities
TPZGeoMesh * ReadGeometry();
void PrintGeometry(TPZGeoMesh * gmesh);

int main() 
{
    TPZGeoMesh * gmesh = ReadGeometry();
    PrintGeometry(gmesh);

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
