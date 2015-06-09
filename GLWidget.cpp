
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
//#if (QT_VERSION >= QT_VERSION_CHECK(5, 1, 0))
//    std::cout << "Qt version >= 5.1.0\n";
//#endif
}


GLWidget::~GLWidget()
{
    if(_shaderProgram) delete _shaderProgram;
}

void GLWidget::initializeGL()
{
    QGLFormat glFormat = QGLWidget::format();
    if (!glFormat.sampleBuffers()) { std::cerr << "Could not enable sample buffers." << std::endl; return; }

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
    // fix me: should use kd-tree
    //if(_iterStatus != -1) return;

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
    GetClosestPoints(_selPoint, _rLines, _lLines, _rPoints, _lPoints);

    // SELECTED POINTS
    QVector3D vecCol = QVector3D(0.0, 0.0, 1.0);
    std::vector<MyPoint> points;
    points.push_back(_selPoint);
    PreparePointsVAO(points, &_selPointsVbo, &_selPointsVao, vecCol);

    // R LINES
    vecCol = QVector3D(1.0, 0.0, 1.0);
    PrepareLinesVAO(_rLines, &_rLinesVbo, &_rLinesVao, vecCol);
    this->repaint();

    // L LINES
    vecCol = QVector3D(0.0, 1.0, 1.0);
    PrepareLinesVAO(_lLines, &_lLinesVbo, &_lLinesVao, vecCol);

    // R POINTS
    vecCol = QVector3D(1.0, 0.0, 1.0);
    PreparePointsVAO(_rPoints, &_rPointsVbo, &_rPointsVao, vecCol);

    // L POINTS
    vecCol = QVector3D(0.0, 1.0, 1.0);
    PreparePointsVAO(_lPoints, &_lPointsVbo, &_lPointsVao, vecCol);
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
void GLWidget::ZoomIn() { this->_zoomFactor += 0.5f; }
void GLWidget::ZoomOut() { this->_zoomFactor -= 0.5f; if(this->_zoomFactor < 0.1f) _zoomFactor = 0.1f; }


void GLWidget::InitCurve()
{
    _points.clear();

    MyPoint centerPt(this->_img_width / 2, this->_img_height / 2);

    float addValue = (M_PI * 2.0 / (float)SystemParams::init_slices);
    for(float a = 0.0; a < M_PI * 2.0; a += addValue)
    {
        float xPt = centerPt.x + SystemParams::radius * sin(a);
        float yPt = centerPt.y + SystemParams::radius * cos(a);
        _points.push_back(MyPoint(xPt, yPt));
    }

    // resampling
    ResampleCurve();

    _iterStatus = 1;

}

void GLWidget::CreateCurveVAO()
{
    // POINTS VAO
    QVector3D vecCol = QVector3D(0.0, 0.0, 0.0);
    PreparePointsVAO(_points, &_pointsVbo, &_pointsVao, vecCol);

    // LINES VAO
    vecCol = QVector3D(0.0, 1.0, 0.0);
    std::vector<MyLine> lines;
    for(size_t a = 0; a < _points.size(); a++)
    {
        if(a < _points.size() - 1) { lines.push_back(MyLine(_points[a].x, _points[a].y, _points[a+1].x, _points[a+1].y)); }
        else { lines.push_back(MyLine(_points[a].x, _points[a].y, _points[0].x, _points[0].y)); }
    }
    PrepareLinesVAO(lines, &_linesVbo, &_linesVao, vecCol);
}

void GLWidget::ResampleCurve()
{
    std::vector<MyPoint> tempPoints;

    // fix me
    float dFactor = 0.9f;

    // stupid hack
    MyPoint firstPt = _points[0];
    _points.push_back(firstPt);

    // a really simple resampling (a version from the paper)
    for(size_t a = 0; a < _points.size() - 1;)
    {
        MyPoint pt1 = _points[a];
        MyPoint pt2 = _points[a+1];

        tempPoints.push_back(pt1);

        if(pt1.Distance(pt2) > SystemParams::D * dFactor)   // split
        {
            MyPoint newPt = pt1 + (pt2 - pt1) * 0.5;
            tempPoints.push_back(newPt);
            a++;
        }
        else    // delete / skip next
        {
            a += 2;
        }
    }

    // rebuild curve
    _points.clear();
    for(size_t a = 0; a < tempPoints.size(); a++)
        { _points.push_back(tempPoints[a]); }

    /*
    // uniform resampling (more precise)
    std::vector<MyPoint> tempPoints;
    UniformResample(_points, tempPoints, SystemParams::D);
    _points.clear();
    for(size_t a = 0; a < tempPoints.size(); a++)
    {
        _points.push_back(tempPoints[a]);
    }*/
}

void GLWidget::EvolveCurve()
{
    //std::cout << "iteration " << _currentIter << "\n";
    if(_currentIter == SystemParams::max_iter - 1)
    {
        _iterStatus = -1;
        _drawSelPoint = false;
    }

    std::vector<MyPoint> tempPoints;
    for(size_t a = 0; a < _points.size(); a++)
    {
        tempPoints.push_back(MyPoint(_points[a].x, _points[a].y));
    }

    // 1 --- BROWNIAN MOTION
    // fix me, delta is static
    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint pt = _points[a];

        //float xB = SystemParams::f_b * (float)RandomNumber() * 1.0f * SystemParams::D;
        //float yB = SystemParams::f_b * (float)RandomNumber() * 1.0f * SystemParams::D;

        float randX = RandomNumber() * cos(RandomNumber());
        float randY = RandomNumber() * sin(RandomNumber());

        //std::cout << randX << " " << randY << "\n";

        float xB = SystemParams::f_b * randX * 1.0f * SystemParams::D;
        float yB = SystemParams::f_b * randY * 1.0f * SystemParams::D;

        tempPoints[a] += MyPoint(xB, yB);

    }

    // 2 --- FAIRING
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

        tempPoints[a] += MyPoint(SystemParams::f_f * xFactor, SystemParams::f_f * yFactor);
    }

    // 3 --- ATTRACTION - REPULSION
    for(size_t a = 0; a < _points.size(); a++)
    {
        MyPoint pt = GetAttractionRepulsion2(a);
        std::cout << "AR " << pt.x << " - " << pt.y << "\n";
        tempPoints[a] += pt;
    }

    // update
    for(size_t a = 0; a < _points.size(); a++)
    {
        _points[a] = tempPoints[a];
    }

    ResampleCurve();

    CreateCurveVAO();

    _currentIter++;
}


void GLWidget::PaintCurve()
{
    if(_points.size() == 0) { return; }

    int use_color_location = _shaderProgram->uniformLocation("use_color");
    _shaderProgram->setUniformValue(use_color_location, (GLfloat)1.0);

    glPointSize(5.0f);
    _pointsVao.bind();
    glDrawArrays(GL_POINTS, 0, _points.size());
    _pointsVao.release();

    if(_drawSelPoint)
    {
        glPointSize(15.0f);
        _selPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, 1);
        _selPointsVao.release();

        _rPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, _rPoints.size());
        _rPointsVao.release();

        _lPointsVao.bind();
        glDrawArrays(GL_POINTS, 0, _lPoints.size());
        _lPointsVao.release();

        glLineWidth(2.0f);
        _lLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _lLines.size() * 2);
        _lLinesVao.release();

        glLineWidth(2.0f);
        _rLinesVao.bind();
        glDrawArrays(GL_LINES, 0, _rLines.size() * 2);
        _rLinesVao.release();
    }

    _linesVao.bind();
    glDrawArrays(GL_LINES, 0, _points.size() * 2);
    _linesVao.release();
}

void GLWidget::UniformResample(std::vector<MyPoint>& oriCurve, std::vector<MyPoint>& resampleCurve, double maxDist)
{
    // fix me
    MyPoint lastPt = oriCurve[0];
    oriCurve.push_back(lastPt);

    resampleCurve.clear();

    double resample_size = maxDist;

    int iter = 0;
    double dist = 0.0;

    while(iter < oriCurve.size() - 1)
    {
        double last_dist = oriCurve[iter].Distance(oriCurve[iter + 1]);

        dist += last_dist;

        if (dist >= resample_size)
        {
            //put a point on line
            double _d = last_dist - (dist - resample_size);
            MyPoint cp(oriCurve[iter].x, oriCurve[iter].y);

            MyPoint  cp1 = MyPoint(oriCurve[iter+1].x, oriCurve[iter+1].y);

            MyPoint dirv = cp1 - cp;
            dirv = dirv * (1.0 / dirv.Length());

            MyPoint insertPt1 = cp + dirv * _d;
            resampleCurve.push_back(insertPt1);

            dist = last_dist - _d; //remaining dist

            //if remaining dist to next point needs more sampling... (within some epsilon)
            while (dist - resample_size > 1e-8 )
            {
                MyPoint insertPt2 = resampleCurve[resampleCurve.size() - 1] + dirv * resample_size;
                resampleCurve.push_back(insertPt2);
                dist -= resample_size;
            }
        }
        iter++;
    }

    //oriCurve.pop_back();
}

float GLWidget::RandomNumber()
{
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<> dist(0, SystemParams::dist_std_dev);

    return ((float)std::round(dist(e2))) / ((float)SystemParams::dist_std_dev);
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

    //std::cout << "r " << rLines.size() << " l " << lLines.size() << "\n";
}

MyPoint GLWidget::GetAttractionRepulsion2(int ptIdx)
{
    MyPoint curPt = _points[ptIdx];

    std::vector<MyLine> rLines;
    std::vector<MyLine> lLines;
    std::vector<MyPoint> rPoints;
    std::vector<MyPoint> lPoints;

    GetClosestSegments(ptIdx, rLines, lLines);
    GetClosestPoints(curPt, rLines, lLines, rPoints, lPoints);

    float aX = 0.0f;
    float aY = 0.0f;

    // right
    for(size_t a = 0; a < rPoints.size(); a++)
    {
        MyPoint minVec = curPt - rPoints[a];
        float lVec = minVec.Length();

        // fix me
        if(lVec >= 0.1f * SystemParams::radius_a_r)
        {
            float lj = GetLennardJones(lVec / (SystemParams::D * 1.0f));
            MyPoint fij = (minVec / lVec) * lj;
            aX += fij.x;
            aY += fij.y;
        }

    }

    // left
    for(size_t a = 0; a < lPoints.size(); a++)
    {
        MyPoint minVec = curPt - lPoints[a];
        float lVec = minVec.Length();

        // fix me
        if(lVec >= 0.1f * SystemParams::radius_a_r)
        {
            float lj = GetLennardJones(lVec / (SystemParams::D * 1.0f));
            MyPoint fij = (minVec / lVec) * lj;
            aX += fij.x;
            aY += fij.y;
        }
    }

    return MyPoint(aX, aY) * SystemParams::f_a;
}

// this doesn't work
MyPoint GLWidget::GetAttractionRepulsion1(int ptIdx)
{
    MyPoint curPt = _points[ptIdx];

    std::vector<MyLine> rLines;
    std::vector<MyLine> lLines;
    std::vector<MyPoint> rPoints;
    std::vector<MyPoint> lPoints;

    GetClosestSegments(ptIdx, rLines, lLines);
    GetClosestPoints(curPt, rLines, lLines, rPoints, lPoints);

    float aX = 0.0f;
    float aY = 0.0f;

    // right f_i_j
    for(size_t a = 0; a < rPoints.size(); a++)
    {
        MyPoint minVec = curPt - rPoints[a];
        float lVec = minVec.Length();

        float lj = GetLennardJones(lVec / (SystemParams::D * 1.0f));
        //float lj = 0.1f;

        if(lj <= 1e-8) continue;

        MyPoint fij = minVec / lVec * lj;

        //std::cout << "lVec lj fij " << lVec << " " << lj << " (" << fij.x << ", " << fij.y << ")\n";

        //if(isnan(fij.x)) fij.x = 0;
        //if(isnan(fij.y)) fij.y = 0;

        aX += fij.x;
        aY += fij.y;
    }

    // left
    for(size_t a = 0; a < lPoints.size(); a++)
    {
        MyPoint minVec = curPt - lPoints[a];
        float lVec = minVec.Length();

        float lj = GetLennardJones(lVec / (SystemParams::D * 1.0f));
        //float lj = 0.1f;

        if(lj <= 1e-8) continue;

        MyPoint fij = minVec / lVec * lj;

        //std::cout << "lVec lj fij " << lVec << " " << lj << " (" << fij.x << ", " << fij.y << ")\n";

        //if(isnan(fij.x)) fij.x = 0;
        //if(isnan(fij.y)) fij.y = 0;

        aX += fij.x;
        aY += fij.y;
    }

    return MyPoint(aX, aY) * SystemParams::f_a;
}

// this might work
float GLWidget::GetLennardJones(float r)
{
    float dljr = SystemParams::delta_l_j / r;
    float pow_6_dljr = pow(dljr, 6.0);
    return (pow_6_dljr * pow_6_dljr) - pow_6_dljr;
}

void GLWidget::GetClosestPoints(MyPoint curPt,
                                std::vector<MyLine> rLines,
                                std::vector<MyLine> lLines,
                                std::vector<MyPoint>& rPoints,
                                std::vector<MyPoint>& lPoints)
{
    rPoints.clear();
    lPoints.clear();

    // r points
    for(size_t a = 0; a < rLines.size(); a++)
    {
        MyPoint closestPt = GetClosestPointToALine(MyPoint(rLines[a].XA, rLines[a].YA), MyPoint(rLines[a].XB, rLines[a].YB), curPt);
        rPoints.push_back(closestPt);
    }

    // l points
    for(size_t a = 0; a < lLines.size(); a++)
    {
        MyPoint closestPt = GetClosestPointToALine(MyPoint(lLines[a].XA, lLines[a].YA), MyPoint(lLines[a].XB, lLines[a].YB), curPt);
        lPoints.push_back(closestPt);
    }
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
