#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QMatrix4x4>
#include <QtOpenGL/QGLWidget>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>

#include "MyPoint.h"
#include "MyLine.h"



class GLWidget : public QGLWidget
{
    Q_OBJECT

private:
    bool    _isMouseDown;
    float   _zoomFactor;
    QPoint  _scrollOffset;

    int _currentIter;
    int _iterStatus; // -1 do nothing or stop
                     //  0 init
                     //  1 running

    // random
    std::random_device _rd;
    std::mt19937 _e2;
    std::normal_distribution<> _dist;

    // image size
    int _img_width;
    int _img_height;

    // shader
    QOpenGLShaderProgram* _shaderProgram;

    // points
    std::vector<MyPoint>        _points;
    QOpenGLBuffer               _pointsVbo;
    QOpenGLVertexArrayObject    _pointsVao;

    // lines
    QOpenGLBuffer               _linesVbo;
    QOpenGLVertexArrayObject    _linesVao;

    // the selected point
    //bool                        _drawSelPoint;
    //MyPoint                     _selPoint;
    //QOpenGLBuffer               _selPointsVbo;
    //QOpenGLVertexArrayObject    _selPointsVao;

    // image
    QOpenGLTexture* _texture;
    QOpenGLBuffer _imageVbo;
    QOpenGLVertexArrayObject _imageVao;
    QImage _imgOriginal;
    QImage _imgColor;
    //QSize _img_size;

    // for rendering
    int         _mvpMatrixLocation;
    int         _colorLocation;
    QMatrix4x4  _perspMatrix;
    QMatrix4x4  _transformMatrix;

    float _runningTime;

    // right points (debug)
    //std::vector<MyPoint> _rPoints;
    //QOpenGLBuffer _rPointsVbo;
    //QOpenGLVertexArrayObject _rPointsVao;

    // left points (debug)
    //std::vector<MyPoint> _lPoints;
    //QOpenGLBuffer _lPointsVbo;
    //QOpenGLVertexArrayObject _lPointsVao;

    // right lines (debug)
    //std::vector<MyLine> _rLines;
    //QOpenGLBuffer _rLinesVbo;
    //QOpenGLVertexArrayObject _rLinesVao;

    // left lines (debug)
    //std::vector<MyLine> _lLines;
    //QOpenGLBuffer _lLinesVbo;
    //QOpenGLVertexArrayObject _lLinesVao;

    std::vector<std::vector<MyPoint>> _allCPoints;

private:
   void InitCurve();
   void EvolveCurve();
   void PaintCurve();
   void CreateCurveVAO();

   void ResampleCurve();
   void UniformResample(std::vector<MyPoint>& oriCurve, std::vector<MyPoint>& resampleCurve, double maxDist);
   float GetRandomNumber();

   //void GetClosestSegments(int ptIndex, std::vector<MyLine>& cLines);
   void GetClosestPoints();
   void GetClosestPoints(MyPoint curPt,
                         std::vector<MyPoint>& cPoints);

   //void GetClosestSegments(int ptIndex, std::vector<MyLine>& rLines, std::vector<MyLine>& lLines);
   //void GetClosestPoints(MyPoint curPt,
   //                      std::vector<MyLine> rLines,
   //                      std::vector<MyLine> lLines,
   //                      std::vector<MyPoint>& rPoints,
   //                      std::vector<MyPoint>& lPoints);

   QImage LoadImageAsGrayscale(QString img);

   MyPoint GetClosestPointToALine(MyPoint v, MyPoint w, MyPoint p);

   //MyPoint GetAttractionRepulsion1(int ptIdx);
   //MyPoint GetAttractionRepulsion2(int ptIdx);
   MyPoint GetAttractionRepulsion(int ptIdx);
   float GetLennardJones(float r);

   void PreparePointsVAO(std::vector<MyPoint> points, QOpenGLBuffer* ptsVbo, QOpenGLVertexArrayObject* ptsVao, QVector3D vecCol);
   void PrepareLinesVAO(std::vector<MyLine> lines, QOpenGLBuffer* linesVbo, QOpenGLVertexArrayObject* linesVao, QVector3D vecCol);

   void SaveToSvg();

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

    void SetImage(QString img);

    void SetRunningTime(float val) { this->_runningTime = val; }
    int GetCurrentIter() { return _currentIter; }
    int GetNPoints(){return _points.size();}
    bool IsCalculationDone(){ return _iterStatus == -1; }

    void StartEvolution(){ _iterStatus = 0; }

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
