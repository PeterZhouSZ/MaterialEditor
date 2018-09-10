#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::pressAction(QAction * act)
{
    if(act == ui->actionLoad)
    {
        QString filename = QFileDialog::getOpenFileName(this,"Open","/home/wenjing/Dropbox/varquadic_mesh/","*.abq");
        if(filename != "")
        {
            m_viewer->load(filename.toStdString());
            std::cout<<"input file: "<<filename.toStdString()<<std::endl;
        }

    }
}

void MainWindow::doOptimize()
{
    double fx = ui->xLineEdit->text().toDouble();
    double fy = ui->yLineEdit->text().toDouble();
    double fz = ui->zLineEdit->text().toDouble();
    double mag = ui->magLineEdit->text().toDouble();



    vec3 fext(fx,fy,fz);
    if(fext.norm() > 1e-10)
        fext.normalize();
    fext *= mag;

    TetMesh& m = m_viewer->m_fitmat.m_mesh;
    int nv = m.nv;
    int nt = m.nt;

    DMatrix& V  = m.m_V;
    IMatrix& F = m.m_F;
    DMatrix& V0 = m.m_V0;

    std::vector<int>& fixedPoints = m_viewer->m_fixedVerts;
    std::vector<int>& actPoints = m_viewer->m_actedVerts;

    std::vector<int> visPoints;
    std::vector<vec3> visDisp;

    std::vector<vec3> appliedForce(nv,vec3::Zero());
    for(int vid : actPoints)
        appliedForce[vid] = fext;
    std::vector<bool> flag(nv,false);
    for(int vid: fixedPoints)
        flag[vid] = true;
    for(int i = 0; i < nv; i++)
    {
        if(flag[i]) continue;

        visPoints.push_back(i);
        vec3 vd = V.row(i) - V0.row(i);
        visDisp.push_back(vd);

    }
    m_viewer->m_fitmat.Init(fixedPoints,visPoints,visDisp,appliedForce);

    int nr = ui->nRLineEdit->text().toInt();
    double lambda = ui->regLineEdit->text().toDouble();
    m_viewer->m_fitmat.m_lambda = lambda;

    if(ui->useSubCheckBox->isChecked())
        m_viewer->m_fitmat.DoOptimizeSubLap(nr);
    else
        m_viewer->m_fitmat.DoOptimizeFull();

}

void MainWindow::applyForce()
{
    m_viewer->drawMode = Viewer::DrawMode::FRAME;

    double fx = ui->xLineEdit->text().toDouble();
    double fy = ui->yLineEdit->text().toDouble();
    double fz = ui->zLineEdit->text().toDouble();
    double mag = ui->magLineEdit->text().toDouble();

    double Y = ui->yLineEdit_2->text().toDouble();
    double rho = ui->rhoLineEdit->text().toDouble();

    vec3 fext(fx,fy,fz);
    if(fext.norm() > 1e-10)
        fext.normalize();
    fext *= mag;

    TetMesh& m = m_viewer->m_fitmat.m_mesh;
    int nv = m.nv;
    int nt = m.nt;

    DMatrix& V  = m.m_V;
    IMatrix& F = m.m_F;

    std::vector<int>& fixedPoints = m_viewer->m_fixedVerts;
    std::vector<int>& actPoints = m_viewer->m_actedVerts;

    m.m_dofID.resize(3 * nv);
    m.m_dofID.assign(3 * nv, 0);
    for(int vid : fixedPoints)
        for(int i = 0; i < 3; i++)
            m.m_dofID[3 * vid + i] = -1;

    int neq = 0;
    for(int i = 0; i < 3 * nv; i++)
        if(m.m_dofID[i] >= 0)
            m.m_dofID[i] = neq++;

    DSparse K(neq,neq);
    DVector f(neq);
    DVector u(neq);

    f.setZero();
    u.setZero();

    //double Lambda = 0.0;
    //double Mu = 0.0;

    double E = Y;
    double nu = rho;
    double Lambda = nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double Mu = E / (2.0 * (1.0 + nu));

    for(int el = 0; el < nt; el++)
    {
        DMatrix Ke = m.m_KLamda[el] * Lambda + m.m_KMu[el] * Mu;

        for(int ii = 0; ii < 4; ii++)
            for(int jj = 0; jj < 3; jj++)
                for(int mm = 0; mm < 4; mm++)
                    for(int nn = 0; nn < 3; nn++)
                    {
                        int vidA = F(el,ii);
                        int vidB = F(el,mm);
                        int gidA = m.m_dofID[3 * vidA + jj];
                        int gidB = m.m_dofID[3 * vidB + nn];

                        if(gidA >= 0 && gidB >= 0)
                            K.coeffRef(gidA, gidB) += Ke(3 * ii + jj, 3 * mm + nn);
                    }

    }

    for(int vid : actPoints)
        for(int i = 0; i < 3; i++)
        {
            int gid = m.m_dofID[3 * vid + i];
            if(gid >= 0)
                f[gid] = fext[i];
        }

    Eigen::SparseLU<DSparse> solver;
    solver.compute(K);
    u = solver.solve(f);

    for(int i = 0; i < nv; i++)
        for(int j = 0; j < 3; j++)
        {
            int gid = m.m_dofID[3 * i + j];
            if(gid >= 0)
                m.m_V(i,j) = m.m_V0(i,j) + u[gid];
            else
                m.m_V(i,j) = m.m_V0(i,j);
        }

    //std::cout<<"unorm: "<<u.norm()<<std::endl;
    this->statusBar()->showMessage("unorm: " + QString::number(u.norm()));

    m_viewer->update();
}

void MainWindow::setDrawSurface()
{
    m_viewer->drawMode = Viewer::DrawMode::SURFACE;
    m_viewer->update();
}

void MainWindow::setDrawInternal()
{
    m_viewer->drawMode = Viewer::DrawMode::INTERNAL;
    m_viewer->update();
}
