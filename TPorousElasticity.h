//
//  TPorousElasticity.h
//  PorousElasticity
//
//  Created by Omar Durán on 9/5/18.
//

#ifndef TPorousElasticity_h
#define TPorousElasticity_h

#include <stdio.h>
#include <stdio.h>
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZTensor.h"

class TPorousElasticity : public TPZMaterial {
    
private:
    
    /// Logarithmic bulk modulus
    STATE m_kappa;
    
    /// Elastic tensile strengh
    STATE m_pt_el;
    
    /// Initial void ratio
    STATE m_e_0;
    
    /// Initial equivalent pressure stress
    STATE m_p_0;
    
    /// Poisson ratio
    STATE m_nu;
    
    /// Second lamé parameter
    STATE m_mu;
    
    /// Directive for define constant shear modulus calculations
    bool m_is_G_constant_Q;
    
    /// Plain stress directive
    int m_plane_stress;
    
    ///  Intial stress in Voigt form https://en.wikipedia.org/wiki/Voigt_notation
    TPZManVector<STATE,6>  m_sigma_0;
    
public:
    
    virtual int ClassId() const;
    
    TPorousElasticity();
    
    TPorousElasticity(int matid);
    
    TPorousElasticity &operator=(const TPorousElasticity &other);
    
    virtual ~TPorousElasticity();
    
    TPorousElasticity(const TPorousElasticity &other);
    
    void Print(std::ostream & out);
    
    std::string Name() { return "TPorousElasticity"; }
    
    int Dimension() const {return 2;}
    
    int NStateVariables();
    
    void SetPorousElasticity(STATE kappa, STATE pt_el, STATE e_0, STATE p_0)
    {
        m_kappa = kappa;
        m_pt_el = pt_el;
        m_e_0 = e_0;
        m_p_0 = p_0;
    }
    
    void SetShearModulusConstant(STATE G){
        m_is_G_constant_Q = true;
        m_mu = G;
    }

    void SetPoissonRatioConstant(STATE nu){
        m_is_G_constant_Q = false;
        m_nu = nu;
    }
    
    void SetPlaneStrain(){
        m_plane_stress = 0;
    }
    
    void SetPlaneStress(){
        m_plane_stress = 1;
    }
    
    void SetInitialStress(TPZManVector<STATE,6> & sigma_0){
        m_sigma_0 = sigma_0;
    }
    
    void FillDataRequirements(TPZMaterialData &data);
    
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) {
        DebugStop();
    }
    
    void Write(TPZStream &buf, int withclassid) const{
        DebugStop();
    }
    
    void Read(TPZStream &buf, void *context){
        DebugStop();
    }
    
    void De(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De);
    
    void De_Shear_constant(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De);
    
    void De_Poisson_constant(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De);
    
    void G(TPZTensor<STATE> &epsilon, STATE & G, STATE & dG_desp_vol);
    
    void Poisson(TPZTensor<STATE> &epsilon, STATE & nu, STATE & dnu_desp_vol);
    
    void Sigma(TPZTensor<STATE> &epsilon, TPZTensor<STATE> &sigma);
    
};


#endif /* TPorousElasticity_h */
