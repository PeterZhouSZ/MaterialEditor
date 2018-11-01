#include "fitmat.h"
//#include "ceres/ceres.h"

#include "unsupported/Eigen/ArpackSupport"
#include "optimization.h"

#include <fstream>

using namespace std;
using namespace Eigen;

FitMat::FitMat()
{

}

void FitMat::LoadMesh(std::string filename)
{
    if(filename.substr(filename.length() - 4) == ".abq")
        m_mesh.load(filename);
    else
    {
        tetrahedronize(filename);
        m_mesh.loadTetgenFiles(filename);
    }

    // Calculate Laplacian matrix
    std::vector<std::set<int>> vToE(m_mesh.nv); // vertex to element connnection

    for(int el = 0; el < m_mesh.nt; el++)
        for(int i = 0; i < 4; i++)
        {
            int vid = m_mesh.m_F(el,i);
            vToE[vid].insert(el);
        }

    std::vector<std::set<int>> EtoE(m_mesh.nt); // element to element connection
    for(auto& vE : vToE)
        for(int el : vE)
            EtoE[el].insert(vE.begin(),vE.end());

    // Initialize the slot in Laplacian
    std::vector<Eigen::Triplet<double>> triplets; 
    for(auto& EE : EtoE)
        for(int elA : EE)
            for(int elB : EE)
                triplets.push_back(Eigen::Triplet<double>(elA,elB));

    DSparse L(m_mesh.nt,m_mesh.nt);
    L.setFromTriplets(triplets.begin(),triplets.end());

    for(int el = 0; el < m_mesh.nt; el++)
    {
        std::set<int>& EE = EtoE[el];
        for(int elB : EE)
            if(elB != el)
                L.coeffRef(el,elB) = -1; // the adjacent term
        L.coeffRef(el,el) = EE.size() - 1; // the degree term
    }

    m_Lap = L;

}

void FitMat::Init(std::vector<int>& fixPoints,
                  std::vector<int>& visPoints, std::vector<vec3>& visDisp,
                  std::vector<vec3>& appliedForce)
{
    int nt = m_mesh.nt;
    int nv = m_mesh.nv;

    // Initialize mapping between global index and free vertex index
    m_mesh.m_dofID.resize(3 * nv);

    m_mesh.m_dofID.assign(3 * nv, 0);
    for(int vid : fixPoints)
        for(int i = 0; i < 3; i++)
            m_mesh.m_dofID[3 * vid + i] = -1; // set to -1 if the vertex is fixed

    int neq = 0;
    for(int i = 0; i < 3 * nv; i++)
        if(m_mesh.m_dofID[i] >= 0)
            m_mesh.m_dofID[i] = neq++; // set the global index (ignoring the fixed vertex)

    m_neq = neq;

    m_K = DSparse(neq,neq);
    std::vector<Eigen::Triplet<double>> triplets;

    // initialize m_K's slot
    for(int el = 0; el < nt; el++)
    {
        for(int ii = 0; ii < 4; ii++)
            for(int jj = 0; jj < 3; jj++)
                for(int mm = 0; mm < 4; mm++)
                    for(int nn = 0; nn < 3; nn++)
                    {
                        int vidA = m_mesh.m_F(el,ii);
                        int vidB = m_mesh.m_F(el,mm);
                        int gidA = m_mesh.m_dofID[3 * vidA + jj];
                        int gidB = m_mesh.m_dofID[3 * vidB + nn];

                        if(gidA >= 0 && gidB >= 0)
                        {
                            triplets.push_back(Eigen::Triplet<double>(gidA,gidB));
                        }
                    }
    }

    m_K.setFromTriplets(triplets.begin(),triplets.end());

    // Assemble 
    m_visPoints = visPoints;
    m_visDisp.resize(m_visPoints.size() * 3);
    m_visDofs.resize(m_visPoints.size() * 3);
    for(int i = 0; i < m_visPoints.size(); i++)
        for(int j = 0; j < 3; j++)
        {
            int id = m_mesh.m_dofID[3 * m_visPoints[i] + j];
            m_visDofs[3 * i + j] = id; // mapping from global to local for visible points
            m_visDisp[3 * i + j] = visDisp[i][j]; // displacement of visible points
        }

    // Initialize the user-defined force at non-fixed vertex
    // For fixed vertex, we can simply removed it to make the optimization problem unconstrainted
    m_fext.resize(neq);

    for(int vid = 0; vid < nv; vid++)
        for(int i = 0; i < 3; i++)
        {
            int gid = m_mesh.m_dofID[3 * vid + i];
            if(gid >= 0)
                m_fext[gid] = appliedForce[vid][i];
        }
}

// Compute energy and gradient using model reduction
void FitMat::ComputeGradientSub(double* Esub,double& energy,double* grad,double* Jac)
{
    DSparse K = m_K;
    DSparse DK = m_K;
    DVector ubat(m_neq);
    DVector fbat = m_fext;
    DVector dubat(m_visDofs.size());
    DVector LtU(m_neq);
    LtU.setZero();
    int nv = m_mesh.nv;
    int nt = m_mesh.nt;
    int ne = m_Phi.cols();
    DVector Efull = m_Phi * Eigen::Map<DVector>(Esub,ne);

    Eigen::SparseLU<DSparse> solver;
    {
        UpdateStiffness(K, Efull.data());
        solver.compute(K);
        ubat =  solver.solve(fbat); // ubat = u* = K^(-1) * f

        for(int i = 0; i < m_visDofs.size(); i++)
            dubat[i] = ubat[m_visDofs[i]] - m_visDisp[i]; // du = u* - u

        energy = 0.5 * dubat.squaredNorm(); // || K^(-1) * f - u ||^2


        if(grad == nullptr) return;

        for(int i = 0; i < ne; i++)
            grad[i] = 0.0;

        for(int ir = 0; ir < ne; ir++)
        {
            UpdateStiffness(DK, m_Phi.data() + ir * ne); // compute each p(K)/p(q_i)
                                                        // p denote partial derivative
            DVector Df = DK * ubat; // DK: p(K)/p(q_i)   ubat: K^(-1) * f
            Df = solver.solve(Df); // Df = K^(-1) * p(K)/p(q_i) * K^(-1) * f

            for(int i = 0; i < m_visDofs.size(); i++)
                grad[ir] += -Df[m_visDofs[i]] * dubat[i]; // chain rule
        }
    }
}

// Compute the energy and gradient without model reduction
void FitMat::ComputeGradient(double* E,double& energy,double* grad,double* Jac)
{
    DSparse K = m_K;
    DVector ubat(m_neq); // u
    DVector fbat = m_fext; // f
    DVector dubat(m_neq); // u of free vertices
    DVector LtU(m_visDofs.size()); // L is the selected observing vertices, here L contains all the vertices
    LtU.setZero();
    int nv = m_mesh.nv;
    int nt = m_mesh.nt;

    DVector EV(nt);
    for(int el = 0; el < nt; el++)
        EV[el] = E[el];

    Eigen::SparseLU<DSparse> solver;
    {
        UpdateStiffness(K,E);
        solver.compute(K);
        ubat = solver.solve(fbat); // Ku = f -> u* = K^-1 * f

        // Calculate energy of F
        for(int i = 0; i < m_visDofs.size(); i++)
            dubat[i] = ubat[m_visDofs[i]] - m_visDisp[i]; // delta(u) = u* - u

        energy = 0.5 * dubat.squaredNorm();
        DVector LapE = m_Lap * EV;
        double energyLap = 0.5 * m_lambda * LapE.squaredNorm(); // lambda * E^T * L * E

        energy += energyLap;

        if(grad == nullptr) return;

        // Calculate gradient of f
        Eigen::Map<DVector>(grad,nt) = m_lambda * (m_Lap * LapE); // gradient of lambda * E^T * E

        for(int i = 0; i < m_visDofs.size(); i++)
            LtU[m_visDofs[i]] = dubat[i];
        LtU = solver.solve(LtU); // L^t * u = K^-1 * (K^-1 * f - u) 

        for(int el = 0; el < nt; el++)
        {
            DMatrix& Ke = m_mesh.m_KE[el];
            DVector ltue(12); // element's (L^T * u)_e
            ltue.setZero();
            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 3; j++)
                {
                    int vid = m_mesh.m_F(el, i);
                    int gid = 3 * vid + j; // global id
                    int lid = m_mesh.m_dofID[gid]; // local id
                    if(lid >= 0) // not fixed
                        ltue[3 * i + j] = LtU[lid];
                }

            ltue = Ke * ltue; // grad(K)_E * K^(-1) * (K^-1 * f - u) 

            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 3; j++)
                {
                    int vid = m_mesh.m_F(el,i);
                    int gid = 3 * vid + j;
                    int lid = m_mesh.m_dofID[gid];
                    if(lid >= 0)
                        grad[el] += -1.0 * ltue[3 * i + j] * ubat[lid];
                        // -1.0 * K^(-1)grad(K)_E * K^(-1)f * (K^-1 * f - u) 
                }
        }
    }

}

void function1_grad(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
{
    FitMat* fm = (FitMat*)(ptr);
    double* E = const_cast<double*>(&x[0]);
    fm->ComputeGradient(E,func,&grad[0],nullptr);

    std::cout<<"energy: "<<func<<std::endl;
}

void function1_gradSub(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *ptr)
{
    FitMat* fm = (FitMat*)(ptr);
    double* E = const_cast<double*>(&x[0]);
    DMatrix& Phi = fm->m_Phi;
    DMatrix& PhiT = fm->m_PhiT;
    DVector Efull = Phi * Eigen::Map<DVector>(E,Phi.cols());
    DVector gfull(Phi.rows());
    // fm->ComputeGradient(Efull.data(),func,gfull.data(),nullptr);
    // DVector gsub = PhiT * gfull;
    //for(int i = 0; i < x.length(); i++)
        //grad[i] = gsub[i];

    fm->ComputeGradientSub(E,func,&grad[0],nullptr);

    std::cout<<"energy: "<<func<<std::endl;
}

void FitMat::DoOptimizeSubLap(int nr)
{
    int nt = m_mesh.nt;
    int nv = m_mesh.nv;

    ComputeEigenMode(nr);
    int ne = m_Phi.cols();

    DVector E(ne);
    DVector dE(ne);

    double Y = 1.0e6;
    double nu = 0.33;
    double Lambda = nu * Y / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double Mu = Y / (2.0 * (1.0 + nu));
    m_mesh.m_KE.resize(nt);
    for(int el = 0; el < nt; el++)
        m_mesh.m_KE[el] = Lambda * m_mesh.m_KLamda[el] + Mu * m_mesh.m_KMu[el];

    E.setZero();
    E[0] = 0.8e6 / m_Phi(0,0) / Y;

    alglib::real_1d_array x;
    x.setlength(ne);
    x.setcontent(ne,E.data());

    double epsg = 0.0000000001;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0.1;
    alglib::ae_int_t maxits = 30;
    alglib::mincgstate state;
    alglib::mincgreport rep;

    // first run
    alglib::mincgcreate(x, state);
    alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
    alglib::mincgsetstpmax(state, stpmax);
    alglib::mincgoptimize(state, function1_gradSub,nullptr,this);
    alglib::mincgresults(state, x, rep);

    cout << "Optimized Ke:"  << endl;
    for(int i = 0; i < ne; i++)
        cout << x[i] << endl;

    // Map Ee to color
    double *Ee = const_cast<double*>(&x[0]);
    DVector Ev = ComputeEvSub(Ee, ne);
    m_mesh.visualizeEv(Ev);
    DVector EFull = m_Phi * Eigen::Map<DVector>(Ee, m_Phi.cols());
    m_mesh.visualizeEe(EFull);
}

void FitMat::DoOptimizeFull()
{
    int nt = m_mesh.nt;
    int nv = m_mesh.nv;

    DVector E(nt);
    DVector dE(nt);

    double Y = 1.0e6; // Young's modulus
    double nu = 0.33; // Poission's ratio
    double Lambda = nu * Y / ((1.0 + nu) * (1.0 - 2.0 * nu)); // lambda coefficient
    double Mu = Y / (2.0 * (1.0 + nu)); // lambda coefficient
    m_mesh.m_KE.resize(nt);
    for(int el = 0; el < nt; el++)
        m_mesh.m_KE[el] = Lambda * m_mesh.m_KLamda[el] + Mu * m_mesh.m_KMu[el];

    E.setConstant(0.8e6 / Y);

    alglib::real_1d_array x;
    x.setlength(nt);
    x.setcontent(nt, E.data());

    double epsg = 0.0000000001;
    double epsf = 0;
    double epsx = 0;
    double stpmax = 0.1;
    alglib::ae_int_t maxits = 30;
    alglib::mincgstate state;
    alglib::mincgreport rep;

    // first run
    alglib::mincgcreate(x, state);
    alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
    alglib::mincgsetstpmax(state, stpmax);
    alglib::mincgoptimize(state, function1_grad,nullptr,this);
    alglib::mincgresults(state, x, rep);

    cout << "Optimized Ke:"  << endl;
    for(int i = 0; i < nt; i++)
        cout << x[i] << endl;

    // Map Et to color
    double *Et = const_cast<double*>(&x[0]);
    DVector Ev = ComputeEvFull(Et, nt);
    m_mesh.visualizeEv(Ev);
    DVector EFull = Eigen::Map<DVector>(Et, nt);
    m_mesh.visualizeEe(EFull);
}

void FitMat::UpdateStiffness(DSparse &K, double* E)
{
    K *= 0.0;
    for(int el = 0; el < m_mesh.nt; el++)
    {
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 3; j++)
            {
                int vidA = m_mesh.m_F(el,i);
                int gidA = 3 * vidA + j;
                int lidA = m_mesh.m_dofID[gidA];
                if(lidA < 0) continue;

                for(int m = 0; m < 4; m++)
                    for(int n = 0; n < 3; n++)
                    {
                        int vidB = m_mesh.m_F(el,m);
                        int gidB = 3 * vidB + n;
                        int lidB = m_mesh.m_dofID[gidB];
                        if(lidB < 0) continue;
                        K.coeffRef(lidA,lidB) += m_mesh.m_KE[el](3 * i + j, 3 * m + n) * E[el];
                    }
            }
    }
}

void FitMat::ComputeEigenMode(int nr)
{
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<DSparse,Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>, true> aparckSolver;

    aparckSolver.compute(m_Lap, nr, "SM");

    DMatrix Phi = aparckSolver.eigenvectors();
    DVector eigVal = aparckSolver.eigenvalues();

    std::cout<<"eigVal:"<<std::endl<<eigVal<<std::endl;

    m_Phi = Phi;
    m_PhiT = Phi.transpose();

    m_eigVal = eigVal;
}

DVector FitMat::ComputeEvSub(double* E, int ne)
{
    DVector EFull = m_Phi * Eigen::Map<DVector>(E, m_Phi.cols());
    return ComputeEvFull(EFull.data(), m_Phi.cols());
}

DVector FitMat::ComputeEvFull(double* E, int nt)
{
    int nv = m_mesh.nv;

    DVector Ev = DVector::Zero(nv);
    IVector vCount = IVector::Zero(nv);
    for(int i = 0; i < nt; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            Ev(m_mesh.m_F(i, j)) += E[i];
            vCount(m_mesh.m_F(i, j)) += 1;
        }
    }

    for(int i = 0; i < nv; i++)
    {
        if(vCount(i) != 0)
            Ev(i) /= vCount(i);
    }

    cout << "Ev: " << endl << Ev << endl;

    return Ev;
}


