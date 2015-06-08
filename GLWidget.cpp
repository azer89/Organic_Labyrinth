
#include "GLWidget.h"

#include <iostream>
#include <random>

#include <QGLFormat>

#include "SystemParams.h"
#include "VertexData.h"

GLWidget::GLWidget(QGLFormat format, QWidget *parent) :
    QGLWidget(format, parent),
    _isMouseDown(false),
    _zoomFactor(25.0),
    _shaderProgram(0),
    _img_width(10),
    _img_height(10),
    _iterStatus(-1),
    _drawSelPoint(false)

{
#if (QT_VERSION >= QT_VERSION_CHECK(5, 1, 0))
    std::cout << "Qt version >= 5.1.0\n";
#endif
}


GLWidget::~GLWidget()
{
    if(_shaderProgram) delete _shaderProgram;
}

void GLWidget::initializeGL()
{
    QGLFormat glFormat = QGLWidget::format();
    if (!glFormat.sampleBuffers())
        { std::cerr << "Could not enable sample buffers." << std::endl; return; }

    glShadeModel(GL_SMOOTH);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor( 1.0, 1.0, 1.0, 1.0 );
    glEnable(GL_DEPTH_TEST);

    _shaderProgram = new QOpenGLShaderProgram();
    if (!_shaderProgram->addShaderFromSourceFile(QOpenGLShader::Vertex, "../OrganicLabyrinth/shader.vert"))
        { std::cerr << "Cannot load vertex shader." << std::endl; return; }

    if (!_shaderProgram->addShaderFromSourceFile(QOpenGLShader::Fragment, "../OrganicLabyrinth/shader.frag"))
        { std::cerr << "Cannot load fragment shader." << std::endl; return; }

    if ( !_shaderProgram->link() )
        { std::cerr << "Cannot link shaders." << std::endl; return; }


    _shaderProgram->bind();

    _mvpMatrixLocation = _shaderProgram->uniformLocation("mvpMatrix");
    _colorLocation = _shaderProgram->attributeLocation("vertexColor");

    //InitCurve();
}


bool GLWidget::event( QEvent * event )
{
    return QGLWidget::event(event);
}


// This is an override function from Qt but I can't find its purpose
void GLWidget::resizeGL(int width, int height)
{
}

void GLWidget::paintGL()
{
    if(_iterStatus == 0)
    {
        _currentIter = 0;
        _drawSelPoint = false;
        InitCurve();
        std::cout << "init curve done\n";
    }

    if(_iterStatus == 1)
    {
        EvolveCurve();
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, this->width(),  this->height());

    int current_width = width();
    int current_height = height();

    // Set orthographic Matrix
    QMatrix4x4 orthoMatrix;

    orthoMatrix.ortho(0.0 +  _scrollOffset.x(),
                      (float)current_width +  _scrollOffset.x(),
                      (float)current_height + _scrollOffset.y(),
                      0.0 + _scrollOffset.y(),
                      -100, 100);

    // Translate the view to the middle
    QMatrix4x4 transformMatrix;
    transformMatrix.setToIdentity();
    transformMatrix.scale(_zoomFactor);

    _shaderProgram->setUniformValue(_mvpMatrixLocation, orthoMatrix * transformMatrix);

    PaintCurve();
}

// Mouse is pressed
void GLWidget::mousePressEvent(int x, int y)
{
    _isMouseDown = true;

    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;

    // your stuff
    //std::cout << dx << " " << dy << "\n";

    /*
    MyPoint _selPoint;
    QOpenGLBuffer _sPointsVbo;
    QOpenGLVertexArrayObject _sPointsVao;
    */

    // fix me: should use kd-tree

    //if(_iterStatus != -1) return;

    //std::cout << "mouse pressed\n";


    _drawSelPoint = true;
    int idx = 0;
    MyPoint mousePt(dx, dy);
    float minDist = std::numeric_limits<float>::max();
    for(size_t a = 0; a < _points.size(); a++)
    {
        float dist = mousePt.Distance(_points[a]);
        if(dist < minDist)
        {
            minDist = dist;
            _selPoint = _points[a];
            idx = a;
        }
    }

    GetClosestSegments(idx, _rLines, _lLines);

    //std::cout << minDist << " " << _selPoint.x << " " << _selPoint.y << "\n";

    QVector3D vecCol = QVector3D(0.0, 0.0, 1.0);
    std::vector<MyPoint> points;
    points.push_back(_selPoint);
    PreparePointsVAO(points, &_selPointsVbo, &_selPointsVao, vecCol);

    /*
    if(_selPointsVao.isCreated()) { _selPointsVao.destroy(); }

    _selPointsVao.create();
    _selPointsVao.bind();

    QVector<VertexData> pointVertices;
    pointVertices.append(VertexData(QVector3D(_selPoint.x, _selPoint.y,  0), QVector2D(), vecCol));

    _selPointsVbo.create();
    _selPointsVbo.bind();
    _selPointsVbo.allocate(pointVertices.data(), pointVertices.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    _selPointsVao.release();
    */


    // L LINES
    vecCol = QVector3D(0.0, 1.0, 1.0);
    PrepareLinesVAO(_lLines, &_lLinesVbo, &_lLinesVao, vecCol);
    /*
    if(_lLinesVao.isCreated())
    {
        _lLinesVao.destroy();
    }

    _lLinesVao.create();
    _lLinesVao.bind();

    QVector<VertexData> linesVertices;
    for(size_t a = 0; a < _lLines.size(); a++)
    {
        linesVertices.append(VertexData(QVector3D(_lLines[a].XA, _lLines[a].YA,  0), QVector2D(), vecCol));
        linesVertices.append(VertexData(QVector3D(_lLines[a].XB, _lLines[a].YB,  0), QVector2D(), vecCol));

        std:: cout << _lLines[a].XA << _lLines[a].YA << "\n";
        std:: cout << _lLines[a].XB << _lLines[a].YB << "\n";
    }

    _lLinesVbo.create();
    _lLinesVbo.bind();
    _lLinesVbo.allocate(linesVertices.data(), linesVertices.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    _lLinesVao.release();
    */

    // R LINES
    vecCol = QVector3D(1.0, 0.0, 1.0);
    PrepareLinesVAO(_rLines, &_rLinesVbo, &_rLinesVao, vecCol);
    /*
    if(_rLinesVao.isCreated())
    {
        _rLinesVao.destroy();
    }

    _rLinesVao.create();
    _rLinesVao.bind();

    linesVertices.clear();
    for(size_t a = 0; a < _rLines.size(); a++)
    {
        linesVertices.append(VertexData(QVector3D(_rLines[a].XA, _rLines[a].YA,  0), QVector2D(), vecCol));
        linesVertices.append(VertexData(QVector3D(_rLines[a].XB, _rLines[a].YB,  0), QVector2D(), vecCol));
    }

    _rLinesVbo.create();
    _rLinesVbo.bind();
    _rLinesVbo.allocate(linesVertices.data(), linesVertices.size() * sizeof(VertexData));

    offset = 0;

    vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    _rLinesVao.release();
    */
    this->repaint();
}


// Mouse is moved
void GLWidget::mouseMoveEvent(int x, int y)
{
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;

    // your stuff

    this->repaint();
}


// Mouse is released
void GLWidget::mouseReleaseEvent(int x, int y)
{
    _isMouseDown = false;
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;

    // your stuff

    this->repaint();
}

void GLWidget::mouseDoubleClick(int x, int y)
{
    double dx = x + _scrollOffset.x();
    dx /= _zoomFactor;

    double dy = y + _scrollOffset.y();
    dy /= _zoomFactor;
}

void GLWidget::HorizontalScroll(int val) { _scrollOffset.setX(val); }
void GLWidget::VerticalScroll(int val) { _scrollOffset.setY(val); }
void GLWidget::ZoomIn() { this->_zoomFactor += 0.05f; }
void GLWidget::ZoomOut() { this->_zoomFactor -= 0.05f; if(this->_zoomFactor < 0.1f) _zoomFactor = 0.1f; }


void GLWidget::InitCurve()
{
    _points.clear();
    //QVector<VertexData> vertices;

    MyPoint centerPt(this->_img_width / 2, this->_img_height / 2);


    float addValue = (M_PI * 2.0 / (float)SystemParams::init_slices);
    for(float a = 0.0; a < M_PI * 2.0; a += addValue)
    {
        float xPt = centerPt.x + SystemParams::radius * sin(a);
        float yPt = centerPt.y + SystemParams::radius * cos(a);
        _points.push_back(MyPoint(xPt, yPt));
        //vertices.append(VertexData(QVector3D(xPt, yPt,  -10), QVector2D(), vecCol));
    }


    // resampling
    ResampleCurve();

    _iterStatus = 1;

}

void GLWidget::CreatePointVAO()
{
    QVector3D vecCol = QVector3D(1.0, 0.0, 0.0);
    PreparePointsVAO(_points, &_pointsVbo, &_pointsVao, vecCol);
    // POINTS VAO
    /*
    if(_pointsVao.isCreated())
    {
        _pointsVao.destroy();
    }

    _pointsVao.create();
    _pointsVao.bind();

    QVector<VertexData> pointVertices;
    for(size_t a = 0; a < _points.size(); a++)
    {
        pointVertices.append(VertexData(QVector3D(_points[a].x, _points[a].y,  0), QVector2D(), vecCol));
    }

    _pointsVbo.create();
    _pointsVbo.bind();
    _pointsVbo.allocate(pointVertices.data(), pointVertices.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    _pointsVao.release();
    */

    // LINE VAO
    if(_linesVao.isCreated())
    {
        _linesVao.destroy();
    }

    _linesVao.create();
    _linesVao.bind();

    vecCol = QVector3D(0.0, 1.0, 0.0);
    QVector<VertexData> linesVertices;
    for(size_t a = 0; a < _points.size(); a++)
    {
        linesVertices.append(VertexData(QVector3D(_points[a].x, _points[a].y,  0), QVector2D(), vecCol));
        if(a < _points.size() - 1)
        {
            linesVertices.append(VertexData(QVector3D(_points[a+1].x, _points[a+1].y,  0), QVector2D(), vecCol));
        }
    }
    linesVertices.append(VertexData(QVector3D(_points[0].x, _points[0].y,  0), QVector2D(), vecCol));

    _linesVbo.create();
    _linesVbo.bind();
    _linesVbo.allocate(linesVertices.data(), linesVertices.size() * sizeof(VertexData));

    quintptr offset = 0;
    //offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    _linesVao.release();
}

void GLWidget::ResampleCurve()
{
    std::vector<MyPoint> tempPoints;
    UniformResample(_points, tempPoints, SystemParams::D);
    _points.clear();
    for(size_t a = 0; a < tempPoints.size(); a++)
    {
        _points.push_back(tempPoints[a]);
    }
}

void GLWidget::EvolveCurve()
{
    std::cout << "iteration " << _currentIter << "\n";
    if(_currentIter == SystemParams::max_iter - 1)
    {
        _iterStatus = -1;
        _drawSelPoint = false;
    }

    // 1 - brownian
    // fix me, delta is static
    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint pt = _points[a];
        //double fBrownian = SystemParams::f_b * RandomNumber() * 1.0 * SystemParams::D;
        pt.x += SystemParams::f_b * (float)RandomNumber() * 0.05f * SystemParams::D;
        pt.y += SystemParams::f_b * (float)RandomNumber() * 0.05f * SystemParams::D;

        _points[a] = pt;

    }

    ResampleCurve();


    // 2 - fairing
    // fix me, delta is static
    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint curPt = _points[a];
        MyPoint prevPt;
        MyPoint nextPt;

        if(a == 0) { prevPt = _points[_points.size() - 1]; }
        else { prevPt = _points[a - 1]; }

        if(a < _points.size() - 1) { nextPt = _points[a + 1]; }
        else { nextPt = _points[0]; }

        float xFactor = (prevPt.x + nextPt.x) / 2.0f;
        float yFactor = (prevPt.y + nextPt.y) / 2.0f;

        xFactor -= curPt.x;
        yFactor -= curPt.y;

        curPt.x += SystemParams::f_f * xFactor;
        curPt.y += SystemParams::f_f * yFactor;

        _points[a] = curPt;
    }

    ResampleCurve();

    CreatePointVAO();
    //this->repaint();

    _currentIter++;
}


void GLWidget::PaintCurve()
{
    if(_points.size() == 0)
    {
        return;
    }

    int use_color_location = _shaderProgram->uniformLocation("use_color");
    _shaderProgram->setUniformValue(use_color_location, (GLfloat)1.0);


    if(_drawSelPoint)
    {
        glPointSize(15.0f);
        _selPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, 1);
        _selPointsVao.release();

        _lLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _lLines.size() * 2);
        _lLinesVao.release();

        _rLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _rLines.size() * 2);
        _rLinesVao.release();
    }

    glPointSize(5.0f);
    _pointsVao.bind();
    glDrawArrays(GL_POINTS, 0, _points.size());
    _pointsVao.release();

//    _linesVao.bind();
//    glDrawArrays(GL_LINES, 0, _points.size() * 2);
//    _linesVao.release();





}

void GLWidget::UniformResample(std::vector<MyPoint>& oriCurve, std::vector<MyPoint>& resampleCurve, double maxDist)
{

    // fix me
    MyPoint lastPt = oriCurve[0];
    oriCurve.push_back(lastPt);

    resampleCurve.clear();
    //for(int a = 0; a <= N; a++) { resampleCurve.push_back(MyPoint(0,0)); }
    //resampleCurve[0].x = oriCurve[0].x;
    //resampleCurve[0].y = oriCurve[0].y;
    //resampleCurve[N].x = oriCurve[oriCurve.size()-1].x;
    //resampleCurve[N].y = oriCurve[oriCurve.size()-1].y;

    //double pl_length = CurveLength(oriCurve);

    //double resample_size = pl_length / (double) N;
    double resample_size = maxDist;

    int curr = 0;
    double dist = 0.0;
    //int i = 1;
    //for (int i = 1; i < N; )
    while(curr < oriCurve.size() - 1)
    {
        double last_dist = oriCurve[curr].Distance(oriCurve[curr+1]);

        dist += last_dist;

        if (dist >= resample_size)
        {
            //put a point on line
            double _d = last_dist - (dist-resample_size);
            MyPoint cp(oriCurve[curr].x, oriCurve[curr].y);

            //MyPoint cp1(oriCurve[curr+1].x, oriCurve[curr+1].y);

            //MyPoint cp1(oriCurve[0].x, oriCurve[0].y);
            //if(curr <  oriCurve.size() - 1)
            //{
                MyPoint  cp1 = MyPoint(oriCurve[curr+1].x, oriCurve[curr+1].y);
            //}

            MyPoint dirv = cp1-cp;
            dirv = dirv * (1.0 / dirv.Length());

            //resampleCurve[i] = cp + dirv * _d;
            //i++;
            MyPoint insertPt1 = cp + dirv * _d;
            resampleCurve.push_back(insertPt1);

            dist = last_dist - _d; //remaining dist

            //if remaining dist to next point needs more sampling... (within some epsilon)
            while (dist - resample_size > 1e-8 )
            {
                //resampleCurve[i] = resampleCurve[i-1] + dirv * resample_size;
                MyPoint insertPt2 = resampleCurve[resampleCurve.size() - 1] + dirv * resample_size;
                resampleCurve.push_back(insertPt2);

                dist -= resample_size;
                //i++;
            }
        }
        curr++;
    }

    //oriCurve.pop_back();
}

int GLWidget::RandomNumber()
{
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> dist(0, SystemParams::dist_std_dev);

    return std::round(dist(e2));
}

void GLWidget::GetClosestSegments(int ptIndex, std::vector<MyLine>& rLines, std::vector<MyLine>& lLines)
{
    // clear
    rLines.clear();
    lLines.clear();

    int ptSize = _points.size();
    MyPoint oriPt = _points[ptIndex];

    // find left (backward)
    int curIdx = ptIndex;
    bool shouldStop = false;
    for(size_t a = 0; a < SystemParams::search_a_r && !shouldStop; a++)
    {
        MyPoint pt1 = _points[curIdx];
        curIdx--;
        if(curIdx < 0){ curIdx = ptSize - 1; }
        MyPoint pt2 = _points[curIdx];

        if(pt2.Distance(oriPt) > SystemParams::radius_a_r) {shouldStop = true;}
        else { lLines.push_back(MyLine(pt1.x, pt1.y, pt2.x, pt2.y)); }
    }


    // find right (upward)
    curIdx = ptIndex;
    shouldStop = false;
    for(size_t a = 0; a < SystemParams::search_a_r && !shouldStop; a++)
    {
        MyPoint pt1 = _points[curIdx];
        curIdx++;
        if(curIdx == ptSize){ curIdx = 0; }
        MyPoint pt2 = _points[curIdx];

        if(pt2.Distance(oriPt) > SystemParams::radius_a_r) {shouldStop = true;}
        else { rLines.push_back(MyLine(pt1.x, pt1.y, pt2.x, pt2.y)); }
    }

    std::cout << "r " << rLines.size() << " l " << lLines.size() << "\n";
}

MyPoint GLWidget::GetClosestPointToALine(MyPoint v, MyPoint w, MyPoint p)
{
    // Return minimum distance between line segment vw and point p
    double l2 = v.DistanceSquared(w);					   // i.e. |w-v|^2 -  avoid a sqrt
    if (l2 > -1e-8 && l2 < 1e-8) return v;   // v == w case

    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line.
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    double t = (p - v).Dot(w - v) / l2;

    if (t < 0.0)	  { return  v; }       // Beyond the 'v' end of the segment
    else if (t > 1.0) { return  w; }        // Beyond the 'w' end of the segment
    MyPoint projection = v + (w - v) * t;     // Projection falls on the segment
    return projection;
}

void GLWidget::PreparePointsVAO(std::vector<MyPoint> points, QOpenGLBuffer* ptsVbo, QOpenGLVertexArrayObject* ptsVao, QVector3D vecCol)
{
    if(ptsVao->isCreated())
    {
        ptsVao->destroy();
    }

    ptsVao->create();
    ptsVao->bind();

    QVector<VertexData> data;
    for(size_t a = 0; a < points.size(); a++)
    {
        data.append(VertexData(QVector3D(points[a].x, points[a].y,  0), QVector2D(), vecCol));
    }

    ptsVbo->create();
    ptsVbo->bind();
    ptsVbo->allocate(data.data(), data.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    ptsVao->release();
}

//void GLWidget::PrepareLinesVAO(std::vector<MyPoint> points, QOpenGLBuffer* linesVbo, QOpenGLVertexArrayObject* linesVao, QVector3D vecCol)
//{
//}

void GLWidget::PrepareLinesVAO(std::vector<MyLine> lines, QOpenGLBuffer* linesVbo, QOpenGLVertexArrayObject* linesVao, QVector3D vecCol)
{
    if(linesVao->isCreated())
    {
        linesVao->destroy();
    }

    linesVao->create();
    linesVao->bind();

    QVector<VertexData> data;
    for(size_t a = 0; a < lines.size(); a++)
    {
        data.append(VertexData(QVector3D(lines[a].XA, lines[a].YA,  0), QVector2D(), vecCol));
        data.append(VertexData(QVector3D(lines[a].XB, lines[a].YB,  0), QVector2D(), vecCol));
    }

    linesVbo->create();
    linesVbo->bind();
    linesVbo->allocate(data.data(), data.size() * sizeof(VertexData));

    quintptr offset = 0;

    int vertexLocation = _shaderProgram->attributeLocation("vert");
    _shaderProgram->enableAttributeArray(vertexLocation);
    _shaderProgram->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(VertexData));

    offset += sizeof(QVector3D);
    offset += sizeof(QVector2D);

    _shaderProgram->enableAttributeArray(_colorLocation);
    _shaderProgram->setAttributeBuffer(_colorLocation, GL_FLOAT, offset, 3, sizeof(VertexData));

    linesVao->release();
}
