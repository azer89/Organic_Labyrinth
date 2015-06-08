#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <iostream>

#include "SystemParams.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    animTimer = new QTimer(this); connect(animTimer, SIGNAL(timeout()), this, SLOT(AnimationThread()));
    connect(ui->actionStart,	 SIGNAL(triggered()), this, SLOT(SimulationStart()));
}

MainWindow::~MainWindow()
{
    delete ui;
    delete animTimer;
}

void MainWindow::AnimationThread()
{
    this->ui->widget->GetGLWidget()->repaint();
    if(this->ui->widget->GetGLWidget()->IsCalculationDone())
    {
        animTimer->stop();
    }
    //else
    //{
    //
    //}
}

void MainWindow::SimulationStart()
{
    std::cout << "simulation start\n";
    this->ui->widget->GetGLWidget()->StartEvolution();
    this->ui->widget->GetGLWidget()->repaint();
    animTimer->start(0);
}
