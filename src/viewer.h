#ifndef VIEWER_H
#define VIEWER_H

#include "shader.h"
#include "fitmat.h"
#include "QGLViewer/qglviewer.h"

// #include <igl/opengl/create_shader_program.h>

class Viewer : public QGLViewer
{
public:
    Viewer();

public:
    void init();
    void draw();

    void load(std::string filename);

private:
    void drawFrame();
    void drawSurface();
    void drawInternal();

    void initShader();

    Shader shader;

    // void loadShader(std::string shaderName);
    // void useShader(std::string shaderName);

    // std::map<std::string,GLuint> shaders;

public:
    virtual void mousePressEvent(QMouseEvent* e);
    virtual void keyPressEvent(QKeyEvent* e);

public:
    FitMat m_fitmat;

    std::vector<int> m_selectedVerts;
    std::vector<int> m_fixedVerts;
    std::vector<int> m_actedVerts;

    enum DrawMode {FRAME, SURFACE, INTERNAL};
    enum DrawMode drawMode = FRAME;
};

#endif // VIEWER_H
