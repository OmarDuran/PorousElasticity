//
//  TPorousElasticity.cpp
//  PorousElasticity
//
//  Created by Omar Durán on 9/5/18.
//

#include "TPorousElasticity.h"

int TPorousElasticity::ClassId() const {
    return Hash("TPorousElasticity");
}

TPorousElasticity::TPorousElasticity(){

    m_kappa = 0.0;
    m_pt_el = 0.0;
    m_e_0 = 0.0;
    m_p_0 = 0.0;
    m_nu = 0.0;
    m_mu = 0.0;
    m_is_G_constant_Q = false;
    m_sigma_0.Resize(6);
}

TPorousElasticity::TPorousElasticity(int matid) : TPZMaterial(matid){

    m_kappa = 0.0;
    m_pt_el = 0.0;
    m_e_0 = 0.0;
    m_p_0 = 0.0;
    m_nu = 0.0;
    m_mu = 0.0;
    m_is_G_constant_Q = false;
    m_plane_stress = 0;
    m_sigma_0.Resize(6);
    
}

TPorousElasticity & TPorousElasticity::operator=(const TPorousElasticity &other){
   
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_kappa = other.m_kappa;
    m_pt_el = other.m_pt_el;
    m_e_0 = other.m_e_0;
    m_p_0 = other.m_p_0;
    m_nu = other.m_nu;
    m_mu = other.m_mu;
    m_is_G_constant_Q = other.m_is_G_constant_Q;
    m_plane_stress = other.m_plane_stress;
    m_sigma_0 = other.m_sigma_0;
    
    return *this;
}

/** @brief Copy constructor */
TPorousElasticity::TPorousElasticity(const TPorousElasticity &other){
    
    m_kappa = other.m_kappa;
    m_pt_el = other.m_pt_el;
    m_e_0 = other.m_e_0;
    m_p_0 = other.m_p_0;
    m_nu = other.m_nu;
    m_mu = other.m_mu;
    m_is_G_constant_Q = other.m_is_G_constant_Q;
    m_plane_stress = other.m_plane_stress;
    m_sigma_0 = other.m_sigma_0;
    
}

TPorousElasticity::~TPorousElasticity(){

}

int TPorousElasticity::NStateVariables(){
    return 2;
}

void TPorousElasticity::Print(std::ostream & out){
    
    out << "Material Name : " << Name() << "\n";
    out << "Properties for elasticity: \n";
    out << "\t Logarithmic bulk modulus   = "                                 << m_kappa        << std::endl;
    out << "\t Elastic tensile strengh   = "                                  << m_pt_el        << std::endl;
    out << "\t Initial void ratio   = "                                       << m_e_0        << std::endl;
    out << "\t Initial equivalent pressure stress   = "                       << m_p_0        << std::endl;
    out << "\t Poisson ratio   = "                                            << m_nu        << std::endl;
    out << "\t Second lamé parameter   = "                                    << m_mu        << std::endl;
    out << "\t G modulus is constant Q = "                                    << m_is_G_constant_Q << std::endl;
    out << "Plane Problem (m_plane_stress = 0, for Plane Strain conditions) " << m_plane_stress << std::endl;
    out << "\t Initial stress  = "  << std::endl;
    out << "\t m_sigma_0_xx   = "   << m_sigma_0[0] << std::endl;
    out << "\t m_sigma_0_yy   = "   << m_sigma_0[1] << std::endl;
    out << "\t m_sigma_0_zz   = "   << m_sigma_0[2] << std::endl;
    out << "\t m_sigma_0_yz   = "   << m_sigma_0[3] << std::endl;
    out << "\t m_sigma_0_xz   = "   << m_sigma_0[4] << std::endl;
    out << "\t m_sigma_0_xy   = "   << m_sigma_0[5] << std::endl;
    out << "Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}

void TPorousElasticity::FillDataRequirements(TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNeighborSol = true;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
}

void TPorousElasticity::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}

void TPorousElasticity::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> sol_u = data.sol[0];
    TPZFMatrix<STATE> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    TPZFNMatrix<4,STATE> val1loc(bc.Val1()),val2loc(bc.Val2());

    int phru = phiu.Rows();
    short in,jn;
    STATE v2[3];
    TPZFMatrix<STATE> &v1 = val1loc;
    v2[0] = val2loc(0,0);    //    Ux displacement or Tnx
    v2[1] = val2loc(1,0);    //    Uy displacement or Tny
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    switch (bc.Type())
    {
        case 0 :
        {
            //    Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)    += BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            break;
        }
        case 1 :
        {
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //    Tnx
                ef(2*in+1,0)    += -1.0*v2[1]*phiu(in,0)*weight;        //    Tny
            }
            break;
        }
        case 2 :
        {
            //    Mixed condition for each state variable no used here
            //    Elasticity Equation
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += v1(i,j)*data.sol[0][j];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += weight * (v2[0]-res(0,0)) * phiu(in,0);
                ef(2*in+1,0) += weight * (v2[1]-res(1,0)) * phiu(in,0);
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf< this->Dimension(); jdf++)
                    {
                        ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*phiu(in,0)*phiu(jn,0)*weight;
                        //      Not Complete with val2? HERE! PHIL!!!!
                        //      DebugStop();
                    }
                }
            }
            
            break;
        }
        case 3 :
        {
            //    Null Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*( v2[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)    += BIGNUMBER*( v2[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)       += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            break;
        }
        case 4 :
        {
            //    Stress Field as Neumann condition for each state variable
            //    Elasticity Equation
            
            for(in = 0; in < this->Dimension(); in ++){
                v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);
            }
            
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //    Tnx
                ef(2*in+1,0)    += -1.0*v2[1]*phiu(in,0)*weight;        //    Tny
            }
            
            break;
        }
        case 5 :
            //    Normal Pressure condition Pressure value Should be inserted in v2[0]
            //    Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
                for(int jn=0; jn< phru; jn++)
                {
                    for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf < this->Dimension(); jdf++)
                    {
                        ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
                        //      Not Complete with val2? HERE! PHIL!!!!
                        //      DebugStop();
                    }
                }
            }
        }
            break;
        case 6 :
            //    Normal Pressure condition Pressure value Should be inserted in v2[0]
            //    Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
                for(int jn=0; jn< phru; jn++)
                {
                    for(int idf=0; idf < this->Dimension(); idf++) for(int jdf=0; jdf < this->Dimension(); jdf++)
                    {
                        ek(2*in+idf,2*jn+jdf) += v1(idf,jdf)*data.normal[idf]*data.normal[jdf]*phiu(in,0)*phiu(jn,0)*weight;
                        //      Not Complete
                        //      DebugStop();
                    }
                }
            }
        }
            break;
        case 7 :
        {
            //    Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;    // X displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in,2*jn)        += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // X displacement
                }
            }
            
            break;
        }
        case 8 :
        {
            //    Dirichlet condition for uy
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)    += BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;    // y displacement Value
                
                for (jn = 0 ; jn < phru; jn++)
                {
                    //    Contribution for Stiffness Matrix
                    ek(2*in+1,2*jn+1)    += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;    // Y displacement
                }
            }
            
            break;
        }
        default:
        {
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }
}

void TPorousElasticity::De(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De){

    STATE lambda, G, dGdespv;
    this->G(epsilon,G, dGdespv);
//    G = (3*(((1 + m_e_0)*(m_p_0 + m_pt_el))/m_kappa)*(1 - 2*m_nu))/(2*(1 + m_nu));
    lambda = (2.0*G*m_nu)/(1.0-2.0*m_nu);
    
    // Line 0
    De.PutVal(_XX_,_XX_, lambda + 2. * G);
    De.PutVal(_XX_,_YY_, lambda);
    De.PutVal(_XX_,_ZZ_, lambda);
    
    // Line 1
    De.PutVal(_XY_,_XY_, 2. * G);
    
    // Line 2
    De.PutVal(_XZ_,_XZ_, 2. * G);
    
    // Line 3
    De.PutVal(_YY_,_XX_, lambda);
    De.PutVal(_YY_,_YY_, lambda + 2. * G);
    De.PutVal(_YY_,_ZZ_, lambda);
    
    // Line 4
    De.PutVal(_YZ_,_YZ_, 2. * G);
    
    // Line 5
    De.PutVal(_ZZ_,_XX_, lambda);
    De.PutVal(_ZZ_,_YY_, lambda);
    De.PutVal(_ZZ_,_ZZ_, lambda + 2. * G);
    
    /// Nonlinear correction
    TPZManVector<STATE> nl(6);
    TPZFMatrix<STATE> De_nl(6,6,0.0);
    nl[0] = (2*dGdespv*epsilon.YY()*m_nu)/(1 - 2*m_nu) + (2*dGdespv*epsilon.ZZ()*m_nu)/(1 - 2*m_nu) +
    epsilon.XX()*(2*dGdespv + (2*dGdespv*m_nu)/(1 - 2*m_nu));
    nl[1] = 2*dGdespv*epsilon.XY();
    nl[2] = 2*dGdespv*epsilon.XZ();
    nl[3] = (2*dGdespv*epsilon.XX()*m_nu)/(1 - 2*m_nu) +
    (2*dGdespv*epsilon.ZZ()*m_nu)/(1 - 2*m_nu) +
    epsilon.YY()*(2*dGdespv + (2*dGdespv*m_nu)/(1 - 2*m_nu));
    nl[4] =  2*dGdespv*epsilon.YZ();
    nl[5] = (2*dGdespv*epsilon.XX()*m_nu)/(1 - 2*m_nu) + (2*dGdespv*epsilon.YY()*m_nu)/(1 - 2*m_nu) + epsilon.ZZ()*(2*dGdespv + (2*dGdespv*m_nu)/(1 - 2*m_nu));
    
    for (int i =0 ; i < 6; i++) {
        De_nl(i,0) = nl[i];
        De_nl(i,3) = nl[i];
        De_nl(i,5) = nl[i];
    }
    De+=De_nl;
}

void TPorousElasticity::G(TPZTensor<STATE> &epsilon, STATE & G, STATE & dGdesp_vol){
    
    STATE epsv = epsilon.I1();
    G = (3*(1 + m_e_0)*(1 + epsv)*(1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*m_kappa*(1 + m_nu));
    
    dGdesp_vol = (-3*pow(1 + m_e_0,2)*(1 + epsv)*
                  (1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*
     pow(m_kappa,2)*(1 + m_nu)) +
    (3*(1 + m_e_0)*(1 - 2*m_nu)*(m_p_0 + m_pt_el))/
    (2.*exp(((1 + m_e_0)*epsv)/m_kappa)*
     m_kappa*(1 + m_nu));
}

void TPorousElasticity::Sigma(TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma){
 
    STATE trace = epsilon.I1();
    STATE lambda, G, dGdesp_vol;
    this->G(epsilon,G, dGdesp_vol);
//    G = (3*(((1 + m_e_0)*(m_p_0 + m_pt_el))/m_kappa)*(1 - 2*m_nu))/(2*(1 + m_nu));
    lambda = (2.0*G*m_nu)/(1.0-2.0*m_nu);
    
    sigma.Identity();
    sigma.Multiply(trace, lambda);
    sigma.Add(epsilon, 2. * G);
    
}

void TPorousElasticity::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // Getting weight functions
    TPZFMatrix<REAL>  & phi_u     =  data.phi;
    TPZFMatrix<REAL>  & dphi_u    =  data.dphix;
    int n_phi_u = phi_u.Rows();
    int first_u  = 0;
    
    TPZFNMatrix<40,REAL> grad_phi_u(3,n_phi_u);
    TPZAxesTools<REAL>::Axes2XYZ(data.dphix, grad_phi_u, data.axes);
    TPZFMatrix<STATE> dsol_u    = data.dsol[0];
    
    REAL dvdx,dvdy,dudx,dudy;
    REAL duxdx,duxdy,duydx,duydy;
    
    //  Gradient for ux
    duxdx = dsol_u(0,0)*data.axes(0,0)+dsol_u(1,0)*data.axes(1,0); // dux/dx
    duxdy = dsol_u(0,0)*data.axes(0,1)+dsol_u(1,0)*data.axes(1,1); // dux/dy
    
    //  Gradient for uy
    duydx = dsol_u(0,1)*data.axes(0,0)+dsol_u(1,1)*data.axes(1,0); // duy/dx
    duydy = dsol_u(0,1)*data.axes(0,1)+dsol_u(1,1)*data.axes(1,1); // duy/dy
    
    STATE s0_xx = m_sigma_0[0];
    STATE s0_yy = m_sigma_0[1];
    STATE s0_xy = m_sigma_0[5];
    
    TPZTensor<STATE> epsilon;
    TPZTensor<STATE> sigma;
    
    TPZFNMatrix<9,STATE> eps, grad_u(3,3,0.0),grad_u_t;
    grad_u(0,0) = duxdx;
    grad_u(0,1) = duxdy;
    grad_u(1,0) = duydx;
    grad_u(1,1) = duydy;
    
    grad_u.Transpose(&grad_u_t);
    eps = 0.5*(grad_u + grad_u_t);
    
    epsilon.XX() = eps(0,0);
    epsilon.XY() = eps(0,1);
    epsilon.YY() = eps(1,1);
    Sigma(epsilon, sigma);
    
    
    TPZFNMatrix<36,STATE> De(6,6,0.0);
    De.Zero();
    this->De(epsilon,De);
    
    TPZFNMatrix<4,STATE> Deriv(2, 2);
    STATE val;
    for(int iu = 0; iu < n_phi_u; iu++ )
    {
        dvdx = grad_phi_u(0,iu);
        dvdy = grad_phi_u(1,iu);
        
        ef(2*iu + first_u)     +=    -1.0 * weight * (s0_xx*dvdx + s0_xy*dvdy);    // x direction
        ef(2*iu+1 + first_u)   +=    -1.0 * weight * (s0_xy*dvdx + s0_yy*dvdy);    // y direction
        
        if (m_plane_stress == 1)
        {
            /* Plain stress state */
            DebugStop();
        }
        else
        {
            /* Plain Strain State */
            
            ef(2*iu + first_u)     +=    weight * (sigma.XX()*dvdx + sigma.XY()*dvdy);    // x direction
            ef(2*iu+1 + first_u)   +=    weight * (sigma.XY()*dvdx + sigma.YY()*dvdy);    // y direction
            
        }
        
        for(int ju = 0; ju < n_phi_u; ju++)
        {
            
//            dudx = grad_phi_u(0,ju);
//            dudy = grad_phi_u(1,ju);
        
            
            if (this->m_plane_stress == 1)
            {
                DebugStop();
            }
            else
            {
                
                /* Plain Strain State */
//                ek(2*iu + first_u,2*ju + first_u)         += weight*    (De(0,0)*dudx*dvdx    + De(1,1)*(0.5*dudy)*dvdy);
//
//                ek(2*iu + first_u,2*ju+1 + first_u)       += weight*    (De(0,3)*dudy*dvdx    + De(1,1)*(0.5*dudx)*dvdy);
//
//                ek(2*iu+1 + first_u,2*ju + first_u)       += weight*    (De(3,0)*dvdy*dudx    + De(1,1)*(0.5*dudy)*dvdx);
//
//                ek(2*iu+1 + first_u,2*ju+1 + first_u)     += weight*    (De(3,3)*dvdy*dudy    + De(1,1)*(0.5*dudx)*dvdx);
                
                for (int ud = 0; ud < 2; ud++) {
                    for (int vd = 0; vd < 2; vd++) {
                        Deriv(vd, ud) = grad_phi_u(vd, iu) * grad_phi_u(ud, ju);
                    }
                }
                
                val = 2. * De(_XX_, _XX_) * Deriv(0, 0);
                val += De(_XX_, _XY_) * Deriv(0, 1);
                val += 2. * De(_XY_, _XX_) * Deriv(1, 0);
                val += De(_XY_, _XY_) * Deriv(1, 1);
                val *= 0.5;
                ek(2*iu + first_u,2*ju + first_u) += weight * val;
                
                val = De(_XX_, _XY_) * Deriv(0, 0);
                val += 2. * De(_XX_, _YY_) * Deriv(0, 1);
                val += De(_XY_, _XY_) * Deriv(1, 0);
                val += 2. * De(_XY_, _YY_) * Deriv(1, 1);
                val *= 0.5;
                ek(2*iu + first_u,2*ju+1 + first_u) += weight * val;
                
                val = 2. * De(_XY_, _XX_) * Deriv(0, 0);
                val += De(_XY_, _XY_) * Deriv(0, 1);
                val += 2. * De(_YY_, _XX_) * Deriv(1, 0);
                val += De(_YY_, _XY_) * Deriv(1, 1);
                val *= 0.5;
                ek(2*iu+1 + first_u,2*ju + first_u) += weight * val;
                
                val = De(_XY_, _XY_) * Deriv(0, 0);
                val += 2. * De(_XY_, _YY_) * Deriv(0, 1);
                val += De(_YY_, _XY_) * Deriv(1, 0);
                val += 2. * De(_YY_, _YY_) * Deriv(1, 1);
                val *= 0.5;
                ek(2*iu+1 + first_u,2*ju+1 + first_u) += weight * val;
                
            }
        }
    }
    
}

void TPorousElasticity::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake;
    TPZMaterial::Contribute(data, weight, ef);
}

void TPorousElasticity::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL>  &phiu = data.phi;
    TPZManVector<STATE,3> sol_u = data.sol[0];
    TPZFMatrix<STATE> dsol_u = data.dsol[0];
    
    REAL ux = sol_u[0];
    REAL uy = sol_u[1];
    
    TPZFNMatrix<4,STATE> val1loc(bc.Val1()),val2loc(bc.Val2());
    
    int phru = phiu.Rows();
    short in,jn;
    STATE v2[3];
    TPZFMatrix<STATE> &v1 = val1loc;
    v2[0] = val2loc(0,0);    //    Ux displacement or Tnx
    v2[1] = val2loc(1,0);    //    Uy displacement or Tny
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    switch (bc.Type())
    {
        case 0 :
        {
            //    Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)    += BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;    // y displacement Value
                
            }
            
            break;
        }
        case 1 :
        {
            //    Neumann condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //    Tnx
                ef(2*in+1,0)    += -1.0*v2[1]*phiu(in,0)*weight;        //    Tny
            }
            break;
        }
        case 2 :
        {
            //    Mixed condition for each state variable no used here
            //    Elasticity Equation
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += v1(i,j)*data.sol[0][j];
            }
            
            for(in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += weight * (v2[0]-res(0,0)) * phiu(in,0);
                ef(2*in+1,0) += weight * (v2[1]-res(1,0)) * phiu(in,0);
                
            }
            
            break;
        }
        case 3 :
        {
            //    Null Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)      += BIGNUMBER*( v2[0])*phiu(in,0)*weight;    // X displacement Value
                ef(2*in+1,0)    += BIGNUMBER*( v2[1])*phiu(in,0)*weight;    // y displacement Value
                
            }
            
            break;
        }
        case 4 :
        {
            //    Stress Field as Neumann condition for each state variable
            //    Elasticity Equation
            
            for(in = 0; in < this->Dimension(); in ++){
                v2[in] = ( v1(in,0) * data.normal[0] + v1(in,1) * data.normal[1]);
            }
            
            for(in = 0 ; in <phru; in++)
            {
                //    Normal Tension Components on neumann boundary
                ef(2*in,0)      += -1.0*v2[0]*phiu(in,0)*weight;        //    Tnx
                ef(2*in+1,0)    += -1.0*v2[1]*phiu(in,0)*weight;        //    Tny
            }
            
            break;
        }
        case 5 :
            //    Normal Pressure condition Pressure value Should be inserted in v2[0]
            //    Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
            }
        }
            break;
        case 6 :
            //    Normal Pressure condition Pressure value Should be inserted in v2[0]
            //    Elasticity Equation
        {
            TPZFNMatrix<2,STATE> res(2,1,0.);
            for(int i=0; i<2; i++) for(int j=0; j<2; j++)
            {
                res(i,0) += data.normal[i]*bc.Val1()(i,j)*data.sol[0][j]*data.normal[j];
            }
            for(int in = 0 ; in < phru; in++)
            {
                ef(2*in+0,0) += (v2[0]*data.normal[0]-res(0,0)) * phiu(in,0) * weight ;
                ef(2*in+1,0) += (v2[0]*data.normal[1]-res(1,0)) * phiu(in,0) * weight ;
            }
        }
            break;
        case 7 :
        {
            //    Dirichlet condition for each state variable
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in,0)        += BIGNUMBER*(ux - v2[0])*phiu(in,0)*weight;    // X displacement Value
                
            }
            
            break;
        }
        case 8 :
        {
            //    Dirichlet condition for uy
            //    Elasticity Equation
            for(in = 0 ; in < phru; in++)
            {
                //    Contribution for load Vector
                ef(2*in+1,0)    += BIGNUMBER*(uy - v2[1])*phiu(in,0)*weight;    // y displacement Value
                
            }
            
            break;
        }
        default:
        {
            PZError << "TPZMatElasticity2D::ContributeBC error - Wrong boundary condition type" << std::endl;
            DebugStop();
        }
            break;
    }
    
}

int TPorousElasticity::VariableIndex(const std::string &name){
    //    Elasticity Variables
    if(!strcmp("Displacement",name.c_str()))                return    1;
    if(!strcmp("SolidPressure",name.c_str()))                return    2;
    if(!strcmp("SigmaX",name.c_str()))                        return    3;
    if(!strcmp("SigmaY",name.c_str()))                        return    4;
    if(!strcmp("SigmaZ",name.c_str()))                        return    5;
    if(!strcmp("TauXY",name.c_str()))                        return    6;
    if(!strcmp("EpsX",name.c_str()))                        return    7;
    if(!strcmp("EpsY",name.c_str()))                        return    8;
    if(!strcmp("EpsZ",name.c_str()))                        return    9;
    if(!strcmp("EpsXY",name.c_str()))                        return    10;
    //    PZError << "TPZMatElasticity2D::VariableIndex Error\n";
    
    return TPZMaterial::VariableIndex(name);
}

int TPorousElasticity::NSolutionVariables(int var){
    if(var == 1)    return 3;
    if(var == 2)    return 1;
    if(var == 3)    return 1;
    if(var == 4)    return 1;
    if(var == 5)    return 1;
    if(var == 6)    return 1;
    if(var == 7)    return 1;
    if(var == 8)    return 1;
    if(var == 9)    return 1;
    if(var == 10)    return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPorousElasticity::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    Solout.Resize(this->NSolutionVariables(var));
    
    TPZManVector<STATE,3> SolU, SolP;
    TPZFNMatrix <6,STATE> DSolU, DSolP;
    TPZFNMatrix <9> axesU, axesP;
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(3);
    TPZFMatrix<STATE> flux(5,1);
    
    if (data.sol.size() != 1) {
        DebugStop();
    }
    
    SolU    =    data.sol[0];
    DSolU    =    data.dsol[0];
    axesU    =    data.axes;
    
    
    //    Displacements
    if(var == 1 || var == 0){
        Solout[0] = SolU[0];
        Solout[1] = SolU[1];
        if(var==1) Solout[2] = 0.0;
        return;
    }
    
    m_mu = (3*(((1 + m_e_0)*(m_p_0 + m_pt_el))/m_kappa)*(1 - 2*m_nu))/(2*(1 + m_nu));
    STATE m_lambda = (2.0*m_mu*m_nu)/(1.0-2.0*m_nu);
    
    STATE s0_xx = m_sigma_0[0];
    STATE s0_yy = m_sigma_0[1];
    STATE s0_zz = m_sigma_0[2];
    STATE s0_xy = m_sigma_0[5];
    
    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL Tau, DSolxy[2][2];
    REAL divu;
    
    DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0); // dUx/dx
    DSolxy[1][0] = DSolU(0,0)*axesU(0,1)+DSolU(1,0)*axesU(1,1); // dUx/dy
    
    DSolxy[0][1] = DSolU(0,1)*axesU(0,0)+DSolU(1,1)*axesU(1,0); // dUy/dx
    DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1); // dUy/dy
    
    divu = DSolxy[0][0]+DSolxy[1][1]+0.0;
    
    epsx = DSolxy[0][0];// du/dx
    epsy = DSolxy[1][1];// dv/dy
    epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
    REAL C11 = 4*(m_mu)*(m_lambda+m_mu)/(m_lambda+2*m_mu);
    REAL C22 = 2*(m_mu)*(m_lambda)/(m_lambda+2*m_mu);
    
    if (this->m_plane_stress)
    {
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*m_mu*epsxy;
    }
    else
    {
        SigX = ((m_lambda + 2*m_mu)*(epsx) + (m_lambda)*epsy);
        SigY = ((m_lambda + 2*m_mu)*(epsy) + (m_lambda)*epsx);
        SigZ = m_lambda*divu;
        Tau = 2.0*m_mu*epsxy;
    }
    
    
    //    Hydrostatic stress
    if(var == 2)
    {
        Solout[0] = SigX+SigY+SigZ;
        return;
    }
    
    //    Effective Stress x-direction
    if(var == 3) {
        Solout[0] = SigX + s0_xx;
        return;
    }
    
    //    Effective Stress y-direction
    if(var == 4) {
        Solout[0] = SigY + s0_yy;
        return;
    }
    
    //    Effective Stress y-direction
    if(var == 5) {
        Solout[0] = SigZ + s0_zz;
        return;
    }
    
    //    Shear Stress
    if(var == 6) {
        Solout[0] = Tau + s0_xy;
        return;
    }
    
    // epsx
    if (var == 7) {
        Solout[0] = epsx;
    }
    
    // epsy
    if (var == 8) {
        Solout[0] = epsy;
    }
    
    // epsz
    if (var == 9) {
        if (m_plane_stress) {
            Solout[0] = -m_nu*(epsx+epsy);
        }
        else
        {
            Solout[0] = 0.;
        }
    }
    
    // epsxy
    if (var == 10) {
        Solout[0] = epsxy;
    }
}


