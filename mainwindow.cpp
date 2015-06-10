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
    float duration = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;

    ui->timeLabel->setText("Time: " + QString::number(duration));
    ui->pointsLabel->setText("# Points: " + QString::number(this->ui->widget->GetGLWidget()->GetNPoints()));

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

    startTime = std::clock();

    animTimer->start(0);
}
