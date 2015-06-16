#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

#include <ctime>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QTimer* animTimer;

    std::clock_t startTime;

private slots:
    // thread to open image
    void AnimationThread();
    void SimulationStart();

    void FileOpen();
};

#endif // MAINWINDOW_H
