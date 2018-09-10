#ifndef FITMAT_H
#define FITMAT_H

#include "settings.h"

struct TetMesh
{

    void load(std::string filename);

    void PreCompute();

    double getVolume(int el);

    void visualizeEv(DVector &Ev);

    void visualizeEe(DVector &Ee);

    DMatrix& getElementColor();

    void draw();

    void draw(float alpha);

    void calcTetCenter();

    DMatrix& getTetCenter();

    DMatrix m_V0;
    DMatrix m_V;
    IMatrix m_F;
    DMatrix m_N;
    DMatrix m_C;

    DMatrix m_FCenter;
    DMatrix m_Ce; // color per tetrahedrons

    std::vector<DMatrix> m_KLamda,m_KMu;
    std::vector<DMatrix> m_KE;

    int nv = 0,nt = 0;
    std::vector<int>  m_dofID;
};

class FitMat
{
public:
    FitMat();

    void LoadMesh(std::string filename);

    void Init(std::vector<int>& fixPoints,
              std::vector<int>& visPoints, std::vector<vec3>& visDisp,
              std::vector<vec3>& appliedForce);

    void DoOptimizeFull();
    void DoOptimizeSubLap(int nr = 10);

    void ComputeEigenMode(int nr);

    void ComputeGradient(double* E,double& energy,double* grad,double* Jac);
    void ComputeGradientSub(double* E,double& energy,double* grad,double* Jac);

    // Compute average color of each vertex
    DVector computeEvSub(double* E, int ne);
    DVector computeEvFull(double* E, int nt);
protected:
    void UpdateStiffness(DSparse& K,double* E);

public:

    TetMesh          m_mesh;

    DVector          m_fext;

    std::vector<int> m_visPoints;
    DVector          m_visDisp;
    std::vector<int> m_visDofs;

    std::vector<int> m_fixPoints;

    DSparse          m_K;

    int m_neq = 0;
    DSparse          m_Lap;
    double           m_lambda = 0e-5;

    DMatrix          m_Phi,m_PhiT;
    DVector          m_eigVal;



};

#endif // FITMAT_H
