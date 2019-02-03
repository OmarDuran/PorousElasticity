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

void TPorousElasticity::De_Shear_constant(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De){
    
    STATE lambda, nu, dnu_desp_vol;
    this->Poisson(epsilon, nu, dnu_desp_vol);
    lambda = (2.0*m_mu*nu)/(1.0-2.0*nu);
    
    // Line 0
    De.PutVal(_XX_,_XX_, lambda + 2. * m_mu);
    De.PutVal(_XX_,_YY_, lambda);
    De.PutVal(_XX_,_ZZ_, lambda);
    
    // Line 1
    De.PutVal(_XY_,_XY_, 2. * m_mu);
    
    // Line 2
    De.PutVal(_XZ_,_XZ_, 2. * m_mu);
    
    // Line 3
    De.PutVal(_YY_,_XX_, lambda);
    De.PutVal(_YY_,_YY_, lambda + 2. * m_mu);
    De.PutVal(_YY_,_ZZ_, lambda);
    
    // Line 4
    De.PutVal(_YZ_,_YZ_, 2. * m_mu);
    
    // Line 5
    De.PutVal(_ZZ_,_XX_, lambda);
    De.PutVal(_ZZ_,_YY_, lambda);
    De.PutVal(_ZZ_,_ZZ_, lambda + 2. * m_mu);
    
    /// Nonlinear correction
    TPZFMatrix<STATE> De_nl(6,6,0.0);
    REAL denominator = 2.0 * m_mu * ( epsilon.XX() + epsilon.YY() + epsilon.ZZ() ) * dnu_desp_vol;
    REAL constant = ((denominator)/(1.0-2.0*nu))*(1.0/(1.0-2.0*nu));
    
    // Line 0
    De_nl.PutVal(_XX_, _XX_, 1.0);
    De_nl.PutVal(_XX_, _YY_, 1.0);
    De_nl.PutVal(_XX_, _ZZ_, 1.0);
    
    // Line 3
    De_nl.PutVal(_YY_, _XX_, 1.0);
    De_nl.PutVal(_YY_, _YY_, 1.0);
    De_nl.PutVal(_YY_, _ZZ_, 1.0);
    
    // Line 5
    De_nl.PutVal(_ZZ_, _XX_, 1.0);
    De_nl.PutVal(_ZZ_, _YY_, 1.0);
    De_nl.PutVal(_ZZ_, _ZZ_, 1.0);
    
    De+=constant*De_nl;
    
}

void TPorousElasticity::De_Poisson_constant(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De){
    
    STATE lambda, G, dGdespv;
    this->G(epsilon,G, dGdespv);
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
    TPZFMatrix<STATE> De_nl(6,6,0.0);
    REAL constant = (2.0/(-1.0+2.0*m_nu));
    
    // Line 0
    REAL l0_val = constant * ( epsilon.XX() * (m_nu - 1.0) - m_nu * (epsilon.YY() + epsilon.ZZ())) * dGdespv;
    De_nl.PutVal(_XX_, _XX_, l0_val);
    De_nl.PutVal(_XX_, _YY_, l0_val);
    De_nl.PutVal(_XX_, _ZZ_, l0_val);
    
    // Line 1
    REAL l1_val = 2.0 * dGdespv * epsilon.XY();
    De_nl.PutVal(_XY_, _XX_, l1_val);
    De_nl.PutVal(_XY_, _YY_, l1_val);
    De_nl.PutVal(_XY_, _ZZ_, l1_val);
    
    // Line 2
    REAL l2_val = 2.0 * dGdespv * epsilon.XZ();
    De_nl.PutVal(_XZ_, _XX_, l2_val);
    De_nl.PutVal(_XZ_, _YY_, l2_val);
    De_nl.PutVal(_XZ_, _ZZ_, l2_val);
    
    // Line 3
    REAL l3_val = constant * ( epsilon.YY() * (m_nu - 1.0) - m_nu * (epsilon.XX() + epsilon.ZZ())) * dGdespv;
    De_nl.PutVal(_YY_, _XX_, l3_val);
    De_nl.PutVal(_YY_, _YY_, l3_val);
    De_nl.PutVal(_YY_, _ZZ_, l3_val);
    
    // Line 4
    REAL l4_val = 2.0 * dGdespv * epsilon.YZ();
    De_nl.PutVal(_YZ_, _XX_, l4_val);
    De_nl.PutVal(_YZ_, _YY_, l4_val);
    De_nl.PutVal(_YZ_, _ZZ_, l4_val);
    
    // Line 5
    REAL l5_val = constant * ( epsilon.ZZ() * (m_nu - 1.0) - m_nu * (epsilon.XX() + epsilon.YY())) * dGdespv;
    De_nl.PutVal(_ZZ_, _XX_, l5_val);
    De_nl.PutVal(_ZZ_, _YY_, l5_val);
    De_nl.PutVal(_ZZ_, _ZZ_, l5_val);
    
    De+=De_nl;
    
}

void TPorousElasticity::De(TPZTensor<STATE> &epsilon, TPZFMatrix<STATE> & De){

    if (m_is_G_constant_Q) {
        De_Shear_constant(epsilon, De);
    }else{
        De_Poisson_constant(epsilon, De);
    }
    
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

void TPorousElasticity::Poisson(TPZTensor<STATE> &epsilon, STATE & nu, STATE & dnu_desp_vol){
    
    STATE epsv = epsilon.I1();
    nu = -1 + (9*(1 + epsv)*(1 + m_e_0)*(m_pt_el + m_p_0))/(6*(1 + epsv)*(1 + m_e_0)*(m_pt_el + m_p_0) +
                                                            2*exp((epsv*(1 + m_e_0))/m_kappa)*m_kappa*m_mu);
    
    dnu_desp_vol = (-9*exp((epsv*(1 + m_e_0))/ m_kappa)*(1 + m_e_0)*(1 + epsv + m_e_0 + epsv*m_e_0 - m_kappa)*
                    (m_pt_el + m_p_0)* m_mu)/(2.*pow(3*(1 + epsv)*(1 + m_e_0)*(m_pt_el + m_p_0) + exp((epsv*(1 + m_e_0))/                m_kappa)*m_kappa*m_mu,2));
}

void TPorousElasticity::Sigma(TPZTensor<STATE> & epsilon, TPZTensor<STATE> & sigma){
 
    STATE trace = epsilon.I1();
    
    if (m_is_G_constant_Q) {
        STATE lambda, nu, dnu_desp_vol;
        this->Poisson(epsilon,nu, dnu_desp_vol);
        lambda = (2.0*m_mu*nu)/(1.0-2.0*nu);
        sigma.Identity();
        sigma.Multiply(trace, lambda);
        sigma.Add(epsilon, 2. * m_mu);
    }else{
        STATE lambda, G, dGdesp_vol;
        this->G(epsilon,G, dGdesp_vol);
        lambda = (2.0*G*m_nu)/(1.0-2.0*m_nu);
        sigma.Identity();
        sigma.Multiply(trace, lambda);
        sigma.Add(epsilon, 2. * G);
    }
    
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
    if(!strcmp("SigmaX",name.c_str()))                        return    2;
    if(!strcmp("SigmaY",name.c_str()))                        return    3;
    if(!strcmp("SigmaZ",name.c_str()))                        return    4;
    return TPZMaterial::VariableIndex(name);
}

int TPorousElasticity::NSolutionVariables(int var){
    if(var == 1)    return 3;
    if(var == 2)    return 1;
    if(var == 3)    return 1;
    if(var == 4)    return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPorousElasticity::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    Solout.Resize(this->NSolutionVariables(var));
    
    
    // Getting weight functions
    TPZManVector<REAL,3> u = data.sol[0];
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
    
    if (data.sol.size() != 1) {
        DebugStop();
    }
    
    //    Displacements
    if(var == 1 || var == 0){
        Solout[0] = u[0];
        Solout[1] = u[1];
        if(var==1) Solout[2] = 0.0;
        return;
    }
    
    if(var == 2)
    {
        Solout[0] = sigma.XX();
        return;
    }
    
    if(var == 3) {
        Solout[0] = sigma.YY();
        return;
    }
    
    if(var == 4) {
        Solout[0] = sigma.ZZ();
        return;
    }

}


