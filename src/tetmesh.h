#ifndef TETMESH_H
#define TETMESH_H
#include "settings.h"

struct TetMesh
{

    void load(std::string filename);

    void loadTetgenFiles(std::string modelName);

    void PreCompute();

    double getVolume(int el);

    void visualizeEv(DVector &Ev);

    void visualizeEe(DVector &Ee);

    DMatrix& getElementColor();

    void draw();

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

#endif // TETMESH_H