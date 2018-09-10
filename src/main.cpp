#include "mainwindow.h"
#include <QApplication>
#include "viewer.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    Viewer viewer;
    w.setCentralWidget(&viewer);
    w.m_viewer = &viewer;

    std::string filename = "/home/serena/graphics/InvMat/models/bunny.abq";
    viewer.load(filename);

    w.show();

    return a.exec();
}
