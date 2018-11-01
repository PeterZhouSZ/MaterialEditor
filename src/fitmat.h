#ifndef FITMAT_H
#define FITMAT_H

#include "settings.h"
#include "tetmesh.h"

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
    DVector ComputeEvSub(double* E, int ne);
    DVector ComputeEvFull(double* E, int nt);

    TetMesh          m_mesh;

    DVector          m_fext;

    std::vector<int> m_visPoints; // visible points: points that are not fixed
    DVector          m_visDisp; // displacement of visible points
    std::vector<int> m_visDofs; // mapping from global to local of visible points

    std::vector<int> m_fixPoints;

    DSparse          m_K;

    int m_neq = 0;
    DSparse          m_Lap;
    double           m_lambda = 0e-5;

    DMatrix          m_Phi, m_PhiT;
    DVector          m_eigVal;

protected:
    void UpdateStiffness(DSparse& K,double* E);

};

#endif // FITMAT_H
