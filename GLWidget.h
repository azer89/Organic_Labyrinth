#ifndef GLWIDGET_H
#define GLWIDGET_H

//#include "stdafx.h"
//#include "MyPoint.h"
//#include "MyLine.h"
//#include "MyIndexedLine.h"

#include <QMatrix4x4>
#include <QtOpenGL/QGLWidget>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

class GLWidget : public QGLWidget
{
    Q_OBJECT

private:
    bool    _isMouseDown;
    float   _zoomFactor;
    QPoint  _scrollOffset;

    // image size
    int _img_width;
    int _img_height;

    // shader
    QOpenGLShaderProgram* _shaderProgram;

    // centroids
    QOpenGLBuffer _pointsVbo;
    QOpenGLVertexArrayObject _pointsVao;

    // for rendering
   int _mvpMatrixLocation;
   int _colorLocation;
   QMatrix4x4 _perspMatrix;
   QMatrix4x4 _transformMatrix;

protected:
    // qt event
    bool event( QEvent * event );
    // init opengl
    void initializeGL();
    // draw
    void paintGL();

    void resizeGL(int width, int height);


public:

    // constructor
    GLWidget( QGLFormat format, QWidget *parent = 0);
    // destructor
    ~GLWidget();

    QSize GetCanvasSize() { return QSize(_img_width, _img_height); }

    // zoom in handle
    void ZoomIn();
    // zoom out handle
    void ZoomOut();
    // set zoom value
    void SetZoom(int val){this->_zoomFactor = val;}
    // get zoom value
    float GetZoomFactor() { return this->_zoomFactor; }

    // set horizontal scroll position
    void HorizontalScroll(int val);
    // set vertical scroll position
    void VerticalScroll(int val);
    // get scroll position (horizontal and vertical)
    QPoint GetScrollOffset() {return this->_scrollOffset;}

    // mouse press
    void mousePressEvent(int x, int y);
    // mouse move
    void mouseMoveEvent(int x, int y);
    // mouse release
    void mouseReleaseEvent(int x, int y);
    // mouse double click
    void mouseDoubleClick(int x, int y);
};

#endif // GLWIDGET_H
