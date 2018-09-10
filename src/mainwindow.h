#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "viewer.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void pressAction(QAction*);
    void applyForce();
    void doOptimize();
    void setDrawInternal();
    void setDrawSurface();
private:
    Ui::MainWindow *ui;

public:
    Viewer* m_viewer;
};

#endif // MAINWINDOW_H
