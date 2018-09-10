#include "viewer.h"
using namespace Eigen;

Viewer::Viewer()
{

}

void Viewer::load(std::string filename)
{
    m_fitmat.LoadMesh(filename);
    DMatrix& V = m_fitmat.m_mesh.m_V;
    vec3 minB = V.colwise().minCoeff();
    vec3 maxB = V.colwise().maxCoeff();

    camera()->setSceneBoundingBox(qglviewer::Vec(minB[0],minB[1],minB[2]),qglviewer::Vec(maxB[0],maxB[1],maxB[2]));
    camera()->showEntireScene();
    update();

    // Load shader
    // loadShader("phong");
}

void Viewer::init()
{
    setBackgroundColor(QColor(255,255,255));
    glDisable(GL_LIGHTING);

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // initShader();
}

//void Viewer::initShader()
//{
//    shader = Shader("shaders/phong/phong.vert", "shaders/phong/phong.frag");

//    // transparent VAO
//    unsigned int VAO, VBO;
//    glGenVertexArrays(1, &VAO);
//    glGenBuffers(1, &VBO);
//    glBindVertexArray(VAO);
//    glBindBuffer(GL_ARRAY_BUFFER, VBO);
//    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
//    glEnableVertexAttribArray(0);
//    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
//    glEnableVertexAttribArray(1);
//    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
//    glBindVertexArray(0);
//}

void Viewer::mousePressEvent(QMouseEvent* e)
{
    if((e->modifiers() & Qt::ControlModifier) && e->button() == Qt::RightButton)
    {
        raise();
        makeCurrent();

        bool found;
        qglviewer::Vec p = camera()->pointUnderPixel(e->pos(),found);
        if(found)
        {
            double minDist = std::numeric_limits<double>::max();
            DMatrix& V = m_fitmat.m_mesh.m_V;

            vec3 p3(p[0],p[1],p[2]);

            int minID = 0;
            for(int i = 0; i < V.rows(); i++)
            {
                vec3 v = V.row(i);
                vec3 dp = p3 - v;
                if(dp.norm() < minDist)
                {
                    minDist = dp.norm();
                    minID = i;
                }

            }

            m_selectedVerts.push_back(minID);
        }
    }

    QGLViewer::mousePressEvent(e);
}

void Viewer::draw()
{
    if(drawMode == FRAME)
        this->drawFrame();
    else if(drawMode == SURFACE)
        this->drawSurface();
    else if(drawMode == INTERNAL)
        this->drawInternal();
}

void Viewer::drawInternal()
{
    // Draw the opaque color
    // configure global opengl state
    // -----------------------------
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    TetMesh& m = m_fitmat.m_mesh;

    std::map<double, int> sorted;
    DMatrix& tetCenter = m_fitmat.m_mesh.getTetCenter();
    qglviewer::Vec cameraPosVec = camera()->position();
    Vector3d cameraPos = Vector3d(cameraPosVec[0], cameraPosVec[1], cameraPosVec[2]);
    for (unsigned int i = 0; i < m_fitmat.m_mesh.nt; i++)
    {
        Vector3d center = tetCenter.row(i);
        double distance = (cameraPos - center).norm();
        sorted[distance] = i;
    }

    DMatrix color = m_fitmat.m_mesh.getElementColor();

    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0,1.0);
    glBegin(GL_TRIANGLES);
    for(auto it = sorted.rbegin(); it != sorted.rend(); ++it)
    {
        int el = it->second;
        glColor4f(color(el, 0), color(el, 1), color(el, 2), 0.08);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                int lid = (i + j) % 4;
                int gid = m.m_F(el,lid);
                vec3 v = m.m_V.row(gid);
                glVertex3d(v[0],v[1],v[2]);
            }  
        }
            
    }
    glEnd();

    // Draw the outline
    glDisable(GL_POLYGON_OFFSET_FILL);

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(0.5,0.5);

    glBegin(GL_TRIANGLES);
    for(auto it = sorted.rbegin(); it != sorted.rend(); ++it)
    {
        int el = it->second;
        glColor4f(color(el, 0), color(el, 1), color(el, 2), 0.5);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                int lid = (i + j) % 4;
                int gid = m.m_F(el,lid);
                vec3 v = m.m_V.row(gid);
                glVertex3d(v[0],v[1],v[2]);
            }  
        }
            
    }
    glEnd();
}

void Viewer::drawSurface()
{
    glDisable(GL_BLEND);

    m_fitmat.m_mesh.draw();

    // Draw the line
//    TetMesh& m = m_fitmat.m_mesh;

//    glDisable(GL_POLYGON_OFFSET_FILL);

//    glColor3f(0.0,0.0,0.0);
//    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
//    glEnable(GL_POLYGON_OFFSET_LINE);
//    glPolygonOffset(0.5,0.5);

//    glBegin(GL_TRIANGLES);
//    for(int el = 0; el < m.nt; el++)
//    {
//        for(int i = 0; i < 4; i++)
//            for(int j = 0; j < 3; j++)
//            {
//                int lid = (i + j) % 4;
//                int gid = m.m_F(el,lid);
//                vec3 v = m.m_V.row(gid);
//                glVertex3d(v[0],v[1],v[2]);
//            }
//    }
//    glEnd();
}

void Viewer::drawFrame()
{
    glDisable(GL_BLEND);

    TetMesh& m = m_fitmat.m_mesh;

    glColor3f(0.8,0.8,0.8);

    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0,1.0);
    glBegin(GL_TRIANGLES);
    for(int el = 0; el < m.nt; el++)
    {
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 3; j++)
            {
                int lid = (i + j) % 4;
                int gid = m.m_F(el,lid);
                vec3 v = m.m_V.row(gid);
                glVertex3d(v[0],v[1],v[2]);
            }
    }
    glEnd();

    glDisable(GL_POLYGON_OFFSET_FILL);

    glColor3f(0.0,0.0,0.0);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(0.5,0.5);

    glBegin(GL_TRIANGLES);
    for(int el = 0; el < m.nt; el++)
    {
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 3; j++)
            {
                int lid = (i + j) % 4;
                int gid = m.m_F(el,lid);
                vec3 v = m.m_V.row(gid);
                glVertex3d(v[0],v[1],v[2]);
            }
    }
    glEnd();

    glColor3f(0.0,1.0,0.0);

    glDisable(GL_POLYGON_OFFSET_LINE);
    glPointSize(8);
    glBegin(GL_POINTS);

    for(int vid : m_selectedVerts)
    {
        vec3 v = m.m_V.row(vid);
        glVertex3d(v[0],v[1],v[2]);
    }

    glEnd();

    glColor3f(1.0,0.0,0.0);

    glDisable(GL_POLYGON_OFFSET_LINE);
    glPointSize(8);
    glBegin(GL_POINTS);

    for(int vid : m_fixedVerts)
    {
        vec3 v = m.m_V.row(vid);
        glVertex3d(v[0],v[1],v[2]);
    }

    glEnd();

    glColor3f(0.0,0.0,1.0);

    glDisable(GL_POLYGON_OFFSET_LINE);
    glPointSize(8);
    glBegin(GL_POINTS);

    for(int vid : m_actedVerts)
    {
        vec3 v = m.m_V.row(vid);
        glVertex3d(v[0],v[1],v[2]);
    }

    glEnd();

}

// void Viewer::loadShader(std::string shaderName)
// {
//     const string shaderPrefix = "~/graphics/InvMat/src/shaders/";
//     opengl::create_shader_program(shaderPrefix + shaderName + "/" + shaderName + ".vert",
//                                   shaderPrefix + shaderName + "/" + shaderName + ".frag",
//                                   this->shaders);
// }

// void Viewer::useShader(std::string shaderName)
// {
//     GLuint id = this->shaders[shaderName];
//     glUseProgram(id);
// }


void Viewer::keyPressEvent(QKeyEvent* e)
{
    if(e->key() == Qt::Key_F)
    {
        m_fixedVerts = m_selectedVerts;
        m_selectedVerts.clear();
        update();
        return;
    }

    if(e->key() == Qt::Key_L)
    {

        m_actedVerts = m_selectedVerts;
        m_selectedVerts.clear();
        update();
        return;
    }

    QGLViewer::keyPressEvent(e);
}
