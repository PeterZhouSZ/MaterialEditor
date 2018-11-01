#include "tetmesh.h"

#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
#include <igl/opengl2/draw_mesh.h>
//#include "ceres/ceres.h"

#include "abqparser.h"
#include <fstream>

using namespace std;
using namespace Eigen;
using namespace igl;

void TetMesh::load(std::string filename)
{
    //const char* etype="C3D4";
    //char ** elementType = &etype;
    printf("Parsing %s...\n",filename);fflush(NULL);

    std::vector<int> v(4);

    // parse the abaqus file
    ABQParser abqParser;

    if (abqParser.open(filename.c_str()) != 0)
    {
        printf("Error: could not open file %s.\n",filename);
        return;
    }

    int mode = 0;

    // pass 1: determine number of vertices and number of elements

    int numVertices = 0;
    int numElements = 0;


    char s[4096];

    while (abqParser.getNextLine(s) != 0)
    {
        if ((mode == 0) && (strncmp(s,"*NODE",5)==0))
        {
            mode = 1;
            continue;
        }

        // seek for *ELEMENT
        if ((mode == 1) && (strncmp(s,"*ELEMENT",8)==0))
        {
            mode = 2;
            continue;
        }

        if ((mode == 2) && (s[0] == '*'))
        {
            mode = 3; // reached *ELSET
        }

        // in mode 1, increase numVertices counter
        if (mode == 1)
            numVertices++;

        // in mode 2, read element info
        if (mode == 2)
            numElements++;


    }

    abqParser.rewindToStart();


    m_V.resize(numVertices,3);
    m_F.resize(numElements,4);

    // pass 3: read vertices and elements

    abqParser.rewindToStart();

    numElements = 0;
    numVertices = 0;
    mode = 0;


    while (abqParser.getNextLine(s) != 0)
    {
        //printf("%s\n", s);
        // now, next input line is in s

        // seek for *NODE
        if ((mode == 0) && (strncmp(s,"*NODE",5)==0))
        {
            mode = 1;
            continue;
        }

        // seek for *ELEMENT
        if ((mode == 1) && (strncmp(s,"*ELEMENT",8)==0))
        {
            mode = 2;
            continue;
        }

        if ((mode == 2) && (s[0] == '*'))
        {
            mode = 3; // reached *ELSET
        }

        // in mode 1, read the vertex info
        if (mode == 1)
        {
            int index;
            double x,y,z;
            sscanf(s,"%d,%lf,%lf,%lf",&index,&x,&y,&z);

            m_V.row(numVertices) = vec3(x,y,z);

            numVertices++;
        }

        // in mode 2, read element info
        if (mode == 2)
        {
            int index;
            sscanf(s,"%d",&index);
            char * ch = s;
            for(int i=0; i<4; i++)
            {
                // seek for next comma
                while((*ch != ',') && (*ch != 0))
                    ch++;

                if (*ch == 0)
                {
                    printf("Error parsing line %s.\n", s);
                    //free(v);
                    throw 12;
                }

                ch++;
                sscanf(ch,"%d",&v[i]);
            }

            for (int k=0; k<4; k++)
            {
                v[k]--; // compensate for fact that vertices are 1-numbered in abq files
            }



            for(int j=0; j<4; j++) // vertices are 1-indexed in .ele files
            {

                m_F(numElements,j) = v[j];
            }

            numElements++;
        }
    }
    abqParser.close();

    // Compute the normals
    per_vertex_normals(m_V, m_F, m_N);

    nv = m_V.rows();
    nt = m_F.rows();

    m_V0 = m_V;
    PreCompute();
    calcTetCenter();

    std::cout<<"number of verties:  "<<nv<<std::endl;
    std::cout<<"number of elements:  "<<nt<<std::endl;
}

void TetMesh::loadTetgenFiles(std::string modelName)
{
    // Read vertices
    ifstream finVtx(_MODEL_PREFIX_ + modelName + "/" + modelName + ".1.node");
    if(!finVtx.is_open())
        cout << "Unable to open vertex file" << endl;
    int row, col, tmp, row_index;
    finVtx >> row >> col >> tmp >> tmp;
    m_V.resize(row, col);
    for(int i = 0; i < row; i++)
    {
        finVtx >> row_index;
        for(int j = 0; j < col; j++)
            finVtx >> m_V(i, j);
    }

    // read elements/tetrahedrons
    ifstream finTet(_MODEL_PREFIX_ + modelName + "/" + modelName + ".1.ele");
    if(!finTet.is_open())
        cout << "Unable to open element file" << endl;
    finTet >> row >> col >> tmp;
    m_F.resize(row, col);
    for(int i = 0; i < row; i++)
    {
        finTet >> row_index;
        for(int j = 0; j < col; j++)
            finTet >> m_F(i, j);
    }

    // Compute the normals
    per_vertex_normals(m_V, m_F, m_N);

    nv = m_V.rows();
    nt = m_F.rows();

    m_V0 = m_V;
    PreCompute();
    calcTetCenter();

    std::cout<<"number of verties:  "<<nv<<std::endl;
    std::cout<<"number of elements:  "<<nt<<std::endl;
}

// Adapted from Vega FEM - corotationalLinearFEM
// Symbol representation reference to the original paper
void TetMesh::PreCompute()
{
    m_KLamda.resize(nt);
    m_KMu.resize(nt);

    // Compute linear tensor of elasticity's gradient at lambda
    double lambda = 1.0;
    double mu = 0;
    double E_lambda[36] = { lambda + 2 * mu, lambda, lambda, 0, 0, 0,
                     lambda, lambda + 2 * mu, lambda, 0, 0, 0,
                     lambda, lambda, lambda + 2 * mu, 0, 0, 0,
                     0, 0, 0, mu, 0, 0,
                     0, 0, 0, 0, mu, 0,
                     0, 0, 0, 0, 0, mu };

    // Compute linear tensor of elasticity's gradient at mu
    lambda = 0.0;
    mu = 1.0;
    double E_mu[36] = { lambda + 2 * mu, lambda, lambda, 0, 0, 0,
                     lambda, lambda + 2 * mu, lambda, 0, 0, 0,
                     lambda, lambda, lambda + 2 * mu, 0, 0, 0,
                     0, 0, 0, mu, 0, 0,
                     0, 0, 0, 0, mu, 0,
                     0, 0, 0, 0, 0, mu };

    for(int el = 0; el < nt; el++)
    {
        double MInv[16];

        {
            DMatrix M(4,4);
            for(int ii = 0; ii < 4; ii++)
                for(int jj = 0; jj < 3; jj++)
                {
                    int vid = m_F(el,ii);
                    M(jj,ii) = m_V(vid,jj);
                    M(3,ii) = 1.0;
                }

            M = M.inverse();
            for(int ii = 0; ii < 4; ii++)
                for(int jj = 0; jj < 4; jj++)
                    MInv[4 * ii + jj] = M(ii,jj);
        }

        // Compute Be, constant matrix that only depends on the tetrahedron’s rest shape
        double B[72] =
        { MInv[0], 0, 0, MInv[4], 0, 0, MInv[8], 0, 0, MInv[12], 0, 0,
          0, MInv[1], 0, 0, MInv[5], 0, 0, MInv[9], 0, 0, MInv[13], 0,
          0, 0, MInv[2], 0, 0, MInv[6], 0, 0, MInv[10], 0, 0, MInv[14],
          MInv[1], MInv[0], 0, MInv[5], MInv[4], 0, MInv[9], MInv[8],
          0, MInv[13], MInv[12], 0, 0, MInv[2], MInv[1], 0, MInv[6], MInv[5],
          0, MInv[10], MInv[9], 0, MInv[14], MInv[13], MInv[2], 0, MInv[0],
          MInv[6], 0, MInv[4], MInv[10], 0, MInv[8], MInv[14], 0, MInv[12] };


        // EB = E * B
        double EB[72];

        // Precompute Ke's gradient at Young's modulus
        {
            for(int ii = 0; ii < 72; ii++)
            EB[ii] = 0.0;

            for (int i=0; i<6; i++)
                for (int j=0; j<12; j++)
                    for (int k=0; k<6; k++)
                        EB[12 * i + j] += E_lambda[6 * i + k] * B[12 * k + j];

            DMatrix& KE_lambda = m_KLamda[el];
            KE_lambda.resize(12,12);
            KE_lambda.setZero();
            for (int i=0; i<12; i++)
                for (int j=0; j<12; j++)
                    for (int k=0; k<6; k++)
                        KE_lambda(i,j) += B[12 * k + i] * EB[12 * k + j];

            KE_lambda *= getVolume(el);
        }

        // Precompute Ke's gradient at Poission's ratio
        {
            for(int ii = 0; ii < 72; ii++)
            EB[ii] = 0.0;

            for (int i=0; i<6; i++)
                for (int j=0; j<12; j++)
                    for (int k=0; k<6; k++)
                        EB[12 * i + j] += E_mu[6 * i + k] * B[12 * k + j];

            DMatrix& KE_mu = m_KMu[el];
            KE_mu.resize(12,12);
            KE_mu.setZero();
            for (int i=0; i<12; i++)
                for (int j=0; j<12; j++)
                    for (int k=0; k<6; k++)
                        KE_mu(i,j) += B[12 * k + i] * EB[12 * k + j];

            KE_mu *= getVolume(el);
        }
    }
}

double TetMesh::getVolume(int el)
{
    //return (1.0 / 6 * fabs( dot(*a - *d, cross(*b - *d, *c - *d)) ));
    vec3 a = m_V.row(m_F(el,0));
    vec3 b = m_V.row(m_F(el,1));
    vec3 c = m_V.row(m_F(el,2));
    vec3 d = m_V.row(m_F(el,3));
    return 1.0/6 * std::abs((a - d).dot((b - d).cross(c - d)));
}

void TetMesh::visualizeEv(DVector &Ev)
{
    jet(Ev, false, m_C);
}

void TetMesh::visualizeEe(DVector &Ee)
{
    jet(Ee, false, m_Ce);
}

DMatrix& TetMesh::getElementColor()
{
    return m_Ce;
}

void TetMesh::draw()
{
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    opengl2::draw_mesh(m_V, m_F, m_N, m_C);
}

void TetMesh::calcTetCenter()
{
    m_FCenter = DMatrix::Zero(m_F.rows(), 3);
    for(int i = 0; i < m_F.rows(); i++)
    {
        for(int j = 0; j < m_F.cols(); j++)
            m_FCenter.row(i) += m_V.row(m_F(i, j));
        m_FCenter.row(i) /= 4.0;
    }
}

DMatrix& TetMesh::getTetCenter()
{
    return m_FCenter;
}
